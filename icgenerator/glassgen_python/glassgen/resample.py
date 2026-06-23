"""Native particle split / merge for multi-resolution glass relaxation.

Ports the geometry of SPH_methods.split_particless / merge_my_neighbor
(setupfiles/SPH_methods.py) but operates directly on glassgen's (pos, mass,
aux-dict) arrays and uses glassgen's own NeighborList for the
nearest-neighbour directions - no dependency on the SPH 21-column format.

split() doubles the particle count per pass (binary 1->2, nn_ort placement:
children sit at parent +/- ortho*dr*dist where ortho is perpendicular to the
nearest-neighbour direction). A target multiplier of 2^p is reached by
repeating the binary pass p times, rebuilding neighbours each pass so each
split follows the *current* local geometry. An adaptive separation floor
(SPH_methods.py:871-888) skips neighbours closer than the relaxed-glass
minimum spacing so children are never placed on top of a just-split sibling.

merge() is the inverse (mutual-nearest pairing). Both are seed operators:
rho and h are placeholders, re-derived by the next relax().

Per-field rules (rule table; extend by adding field names):
  mass        : split -> halve     ; merge -> sum
  pos         : split -> +/- offset ; merge -> mass-weighted centroid
  HALVE_FIELDS: split -> halve     ; merge -> sum      (extensive: u, B, momenti)
  everything else (vel, rho, h, soft, metals, spin, tform, phi, ...):
                split -> inherit   ; merge -> mass-weighted mean
  (vel as a mass-weighted mean is exactly momentum conservation.)
"""
import numpy as np
from numba import njit, prange

from .neighbors import NeighborList

# Fields that are halved on split / summed on merge (extensive per particle),
# matching SPH_methods.split_particless (u/2, B/2, momenti/2 at lines 948-967).
HALVE_FIELDS = frozenset({'u', 'B', 'momenti'})


@njit(cache=True, parallel=True)
def _nn_dirs(pos, idx, h, box, periodic, nsmooth, dist_floor):
    """For each particle find, among its kNN candidates, the primary split
    neighbour (smallest min-image distance whose d/h >= dist_floor, i.e. at
    least the relaxed-glass spacing) and a secondary neighbour (next smallest)
    used to orient the orthogonal split plane.

    Returns vnn (N,3) min-image vector to the primary, dr (N,) its length, and
    vort (N,3) vector to the secondary (the ort_neigh_vector of SPH_methods).
    Falls back to the farthest candidate when none clears the floor."""
    n = pos.shape[0]
    k = idx.shape[1]
    bx, by, bz = box[0], box[1], box[2]
    vnn = np.zeros((n, 3))
    vort = np.zeros((n, 3))
    dr = np.zeros(n)
    for i in prange(n):
        floor = dist_floor * h[i]
        # primary: smallest d with d >= floor; track farthest as fallback
        bp = -1
        bpd = 1e30
        bfar = -1
        bfard = -1.0
        # secondary: smallest d > 0 (any), distinct from primary slot
        for c in range(k):
            j = idx[i, c]
            if j < 0:
                continue
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            if periodic:
                if dx > 0.5 * bx:
                    dx -= bx
                elif dx < -0.5 * bx:
                    dx += bx
                if dy > 0.5 * by:
                    dy -= by
                elif dy < -0.5 * by:
                    dy += by
                if dz > 0.5 * bz:
                    dz -= bz
                elif dz < -0.5 * bz:
                    dz += bz
            d = np.sqrt(dx * dx + dy * dy + dz * dz)
            if d <= 0.0:
                continue
            if d > bfard:
                bfard = d
                bfar = c
            if d >= floor and d < bpd:
                bpd = d
                bp = c
        if bp < 0:
            bp = bfar  # nothing cleared the floor: take the farthest
        # secondary = smallest d distinct from the chosen primary
        bs = -1
        bsd = 1e30
        for c in range(k):
            if c == bp:
                continue
            j = idx[i, c]
            if j < 0:
                continue
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            if periodic:
                if dx > 0.5 * bx:
                    dx -= bx
                elif dx < -0.5 * bx:
                    dx += bx
                if dy > 0.5 * by:
                    dy -= by
                elif dy < -0.5 * by:
                    dy += by
                if dz > 0.5 * bz:
                    dz -= bz
                elif dz < -0.5 * bz:
                    dz += bz
            d = np.sqrt(dx * dx + dy * dy + dz * dz)
            if d <= 0.0:
                continue
            if d < bsd:
                bsd = d
                bs = c
        if bp >= 0:
            j = idx[i, bp]
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            if periodic:
                if dx > 0.5 * bx:
                    dx -= bx
                elif dx < -0.5 * bx:
                    dx += bx
                if dy > 0.5 * by:
                    dy -= by
                elif dy < -0.5 * by:
                    dy += by
                if dz > 0.5 * bz:
                    dz -= bz
                elif dz < -0.5 * bz:
                    dz += bz
            vnn[i, 0] = dx
            vnn[i, 1] = dy
            vnn[i, 2] = dz
            dr[i] = np.sqrt(dx * dx + dy * dy + dz * dz)
        if bs >= 0:
            j = idx[i, bs]
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            if periodic:
                if dx > 0.5 * bx:
                    dx -= bx
                elif dx < -0.5 * bx:
                    dx += bx
                if dy > 0.5 * by:
                    dy -= by
                elif dy < -0.5 * by:
                    dy += by
                if dz > 0.5 * bz:
                    dz -= bz
                elif dz < -0.5 * bz:
                    dz += bz
            vort[i, 0] = dx
            vort[i, 1] = dy
            vort[i, 2] = dz
    return vnn, dr, vort


@njit(cache=True, parallel=True)
def _place_children(pos, vnn, dr, vort, box, periodic, dist, out):
    """Fill out (2N,3): out[2i]/out[2i+1] = pos[i] +/- ortho*dr*dist, ortho a
    unit vector perpendicular to the nearest-neighbour direction (cross with
    the secondary direction; arbitrary-axis fallback when parallel). Children
    wrapped into the box (single image: offset = dist*dr <= ~0.4 h << L/2)."""
    n = pos.shape[0]
    bx, by, bz = box[0], box[1], box[2]
    hbx, hby, hbz = 0.5 * bx, 0.5 * by, 0.5 * bz
    for i in prange(n):
        d = dr[i]
        if d <= 0.0:
            # no usable neighbour: tiny deterministic split along x
            ox, oy, oz = 1.0, 0.0, 0.0
            disp = 1e-6
        else:
            ux, uy, uz = vnn[i, 0] / d, vnn[i, 1] / d, vnn[i, 2] / d
            # ortho = u x vort
            ox = uy * vort[i, 2] - uz * vort[i, 1]
            oy = uz * vort[i, 0] - ux * vort[i, 2]
            oz = ux * vort[i, 1] - uy * vort[i, 0]
            nrm = np.sqrt(ox * ox + oy * oy + oz * oz)
            if nrm < 1e-8:
                # parallel: cross u with the smallest-magnitude axis
                ax, ay, az = abs(ux), abs(uy), abs(uz)
                if ax < ay and ax < az:
                    tx, ty, tz = 1.0, 0.0, 0.0
                elif ay < az:
                    tx, ty, tz = 0.0, 1.0, 0.0
                else:
                    tx, ty, tz = 0.0, 0.0, 1.0
                ox = uy * tz - uz * ty
                oy = uz * tx - ux * tz
                oz = ux * ty - uy * tx
                nrm = np.sqrt(ox * ox + oy * oy + oz * oz)
            ox, oy, oz = ox / nrm, oy / nrm, oz / nrm
            disp = d * dist
        x1 = pos[i, 0] + ox * disp
        y1 = pos[i, 1] + oy * disp
        z1 = pos[i, 2] + oz * disp
        x2 = pos[i, 0] - ox * disp
        y2 = pos[i, 1] - oy * disp
        z2 = pos[i, 2] - oz * disp
        if periodic:
            if x1 > hbx:
                x1 -= bx
            elif x1 < -hbx:
                x1 += bx
            if y1 > hby:
                y1 -= by
            elif y1 < -hby:
                y1 += by
            if z1 > hbz:
                z1 -= bz
            elif z1 < -hbz:
                z1 += bz
            if x2 > hbx:
                x2 -= bx
            elif x2 < -hbx:
                x2 += bx
            if y2 > hby:
                y2 -= by
            elif y2 < -hby:
                y2 += by
            if z2 > hbz:
                z2 -= bz
            elif z2 < -hbz:
                z2 += bz
        out[2 * i, 0] = x1
        out[2 * i, 1] = y1
        out[2 * i, 2] = z1
        out[2 * i + 1, 0] = x2
        out[2 * i + 1, 1] = y2
        out[2 * i + 1, 2] = z2


def _interleave(arr, factor=2):
    """Repeat each row `factor` times contiguously: rows 0,0,1,1,... for the
    children of each parent (matches _place_children's 2i/2i+1 layout)."""
    return np.repeat(arr, factor, axis=0)


def _split_once(pos, mass, aux, box, nsmooth, margin, dist, halve):
    """One binary split pass (N -> 2N)."""
    n = len(pos)
    periodic = box is not None
    boxarr = np.ones(3) if box is None else np.asarray(box, dtype=np.float64)
    nl = NeighborList(boxarr if periodic else None, nsmooth, margin=margin,
                      build_transpose=False)
    idx, h, _, _ = nl.update(np.ascontiguousarray(pos, dtype=np.float64))
    # relaxed-glass minimum spacing / h (SPH_methods nearest_distanceoverhmin)
    dist_floor = 1.0 / (0.5 * nsmooth) ** (1.0 / 3.0)
    vnn, dr, vort = _nn_dirs(pos, idx, h, boxarr, periodic, nsmooth,
                             dist_floor)
    pos2 = np.empty((2 * n, 3))
    _place_children(pos, vnn, dr, vort, boxarr, periodic, dist, pos2)
    mass2 = _interleave(mass) * 0.5
    aux2 = {}
    for name, a in aux.items():
        a = np.asarray(a)
        rep = _interleave(a)
        if name in halve:
            rep = rep * 0.5
        aux2[name] = rep
    return pos2, mass2, aux2


def split(pos, mass, aux=None, box=None, factor=8, nsmooth=64, margin=8,
          dist=0.4, halve_fields=None):
    """Split every particle, multiplying the count by `factor` (a power of 2).

    pos : (N,3). mass : (N,). aux : dict {name: (N,) or (N,3)} of extra fields
    (vel, rho, u, soft, metals, B, spin, momenti, tform, h, phi, ...); names in
    halve_fields are halved (extensive), the rest inherited. halve_fields :
    override the default HALVE_FIELDS set (e.g. add B-field aux keys). box :
    (3,) periodic lengths or None (open). Returns (pos2, mass2, aux2) with
    N*factor rows. Total mass is conserved (each child carries mass/2^p)."""
    if aux is None:
        aux = {}
    halve = HALVE_FIELDS if halve_fields is None else frozenset(halve_fields)
    factor = int(factor)
    if factor < 1 or (factor & (factor - 1)) != 0:
        raise ValueError(f"factor must be a power of 2, got {factor}")
    passes = factor.bit_length() - 1
    pos = np.array(pos, dtype=np.float64)
    mass = np.array(mass, dtype=np.float64)
    aux = {k: np.array(v) for k, v in aux.items()}
    for _ in range(passes):
        pos, mass, aux = _split_once(pos, mass, aux, box, nsmooth, margin,
                                     dist, halve)
    return pos, mass, aux


@njit(cache=True)
def _merge_pairs(pos, mass, idx, h, box, periodic):
    """Mutual-nearest pairing within h_i+h_j. Returns a (N,) partner array:
    partner[i] = j (>i) if (i,j) merge, -1 otherwise; each particle in at most
    one pair. Serial (greedy claim) for determinism."""
    n = pos.shape[0]
    k = idx.shape[1]
    bx, by, bz = box[0], box[1], box[2]
    nn = np.full(n, -1, np.int64)
    nnd = np.full(n, 1e30)
    for i in range(n):
        for c in range(k):
            j = idx[i, c]
            if j < 0:
                continue
            dx = pos[j, 0] - pos[i, 0]
            dy = pos[j, 1] - pos[i, 1]
            dz = pos[j, 2] - pos[i, 2]
            if periodic:
                if dx > 0.5 * bx:
                    dx -= bx
                elif dx < -0.5 * bx:
                    dx += bx
                if dy > 0.5 * by:
                    dy -= by
                elif dy < -0.5 * by:
                    dy += by
                if dz > 0.5 * bz:
                    dz -= bz
                elif dz < -0.5 * bz:
                    dz += bz
            d = np.sqrt(dx * dx + dy * dy + dz * dz)
            if d > 0.0 and d < nnd[i]:
                nnd[i] = d
                nn[i] = j
    partner = np.full(n, -1, np.int64)
    claimed = np.zeros(n, np.bool_)
    for i in range(n):
        if claimed[i]:
            continue
        j = nn[i]
        if j < 0 or claimed[j]:
            continue
        if nn[j] == i and nnd[i] <= h[i] + h[j]:
            a, b = (i, j) if i < j else (j, i)
            partner[a] = b
            claimed[i] = True
            claimed[j] = True
    return partner


@njit(cache=True, parallel=True)
def _nearest(pos, idx, box, periodic):
    """Index of each particle's single nearest neighbour (excluding self),
    scanning its small candidate row. -1 if the row holds no other particle."""
    n, k = idx.shape
    bx, by, bz = box[0], box[1], box[2]
    nn = np.full(n, -1, np.int64)
    for i in prange(n):
        best = -1
        bestd = 1e30
        for c in range(k):
            j = idx[i, c]
            if j < 0 or j == i:
                continue
            dx = pos[i, 0] - pos[j, 0]
            dy = pos[i, 1] - pos[j, 1]
            dz = pos[i, 2] - pos[j, 2]
            if periodic:
                if dx > 0.5 * bx:
                    dx -= bx
                elif dx < -0.5 * bx:
                    dx += bx
                if dy > 0.5 * by:
                    dy -= by
                elif dy < -0.5 * by:
                    dy += by
                if dz > 0.5 * bz:
                    dz -= bz
                elif dz < -0.5 * bz:
                    dz += bz
            d2 = dx * dx + dy * dy + dz * dz
            if d2 < bestd:
                bestd = d2
                best = j
        nn[i] = best
    return nn


def _combine(pos, mass, aux, pa, pb, keep, halve, boxarr, periodic):
    """Combine paired particles (leaders `pa` with partners `pb`) and carry the
    unpaired `keep` through. Mass summed, pos = mass-weighted centroid (min-
    image), halve-fields summed, the rest mass-weighted mean. Shared by merge()
    and coarsen()."""
    m_a, m_b = mass[pa], mass[pb]
    m_sum = m_a + m_b
    off = pos[pb] - pos[pa]
    if periodic:
        off -= boxarr * np.rint(off / boxarr)
    cen = pos[pa] + (m_b / m_sum)[:, None] * off
    if periodic:
        cen -= boxarr * np.rint(cen / boxarr)
    new_pos = np.vstack([pos[keep], cen]) if len(keep) else cen
    new_mass = np.concatenate([mass[keep], m_sum])
    new_aux = {}
    for name, a in aux.items():
        if name in halve:
            merged = a[pa] + a[pb]
        else:
            w = (m_a / m_sum) if a.ndim == 1 else (m_a / m_sum)[:, None]
            wb = (m_b / m_sum) if a.ndim == 1 else (m_b / m_sum)[:, None]
            merged = w * a[pa] + wb * a[pb]
        new_aux[name] = np.concatenate([a[keep], merged]) if len(keep) else merged
    return new_pos, new_mass, new_aux


def merge(pos, mass, aux=None, box=None, target_factor=2, nsmooth=64,
          margin=8, halve_fields=None):
    """Merge particles toward count/target_factor (power of 2) by mutual-
    nearest pairing. Inverse of split: mass summed, pos mass-weighted centroid,
    halve_fields summed, the rest mass-weighted mean (vel -> momentum). Particles
    without a mutual partner are carried through unchanged."""
    if aux is None:
        aux = {}
    halve = HALVE_FIELDS if halve_fields is None else frozenset(halve_fields)
    target_factor = int(target_factor)
    if target_factor < 1 or (target_factor & (target_factor - 1)) != 0:
        raise ValueError(f"target_factor must be a power of 2, got "
                         f"{target_factor}")
    passes = target_factor.bit_length() - 1
    periodic = box is not None
    boxarr = np.ones(3) if box is None else np.asarray(box, dtype=np.float64)
    pos = np.array(pos, dtype=np.float64)
    mass = np.array(mass, dtype=np.float64)
    aux = {k: np.array(v, dtype=np.float64) for k, v in aux.items()}
    for _ in range(passes):
        n = len(pos)
        nl = NeighborList(boxarr if periodic else None, nsmooth,
                          margin=margin, build_transpose=False)
        idx, h, _, _ = nl.update(pos)
        partner = _merge_pairs(pos, mass, idx, h, boxarr, periodic)
        pa = np.where(partner >= 0)[0]      # pair leaders (a < b)
        pb = partner[pa]                    # their partners
        involved = np.zeros(n, dtype=bool)
        involved[pa] = True
        involved[pb] = True
        keep = np.where(~involved)[0]       # unpaired: carry through
        pos, mass, aux = _combine(pos, mass, aux, pa, pb, keep, halve,
                                  boxarr, periodic)
    return pos, mass, aux


def coarsen(pos, mass, aux=None, box=None, factor=8, k=2, margin=4,
            halve_fields=None):
    """Coarsen the particle set to EXACTLY count//factor by repeated mutual-
    nearest pairing, rebuilding the (k-nearest) neighbours each pass until the
    target is reached - the inverse-by-count of split(): ``split(coarsen(x, F),
    F)`` restores the original count (when count is divisible by F).

    This is the iterate-to-target form of SPH_methods.merge_particles (k=2,
    parallel mutual-nearest test, rebuild each round) but WITHOUT its h_i+h_j
    distance cap - the coarse set is only a relaxation seed, so merging the
    nearest pair regardless of separation is fine and converges faster. The
    final pass merges only the `need` CLOSEST mutual pairs so the target count
    is hit exactly. Mass is conserved; rho/h are placeholders, re-derived by the
    next relax().

    Intended for the progressive coarse seed: coarsen the full-res random pre-IC
    down (factor = the cascade's total split factor), relax, then split back to
    the original count - no separate low-res file, count tied 1:1 to the input."""
    if aux is None:
        aux = {}
    halve = HALVE_FIELDS if halve_fields is None else frozenset(halve_fields)
    factor = int(factor)
    if factor < 1 or (factor & (factor - 1)) != 0:
        raise ValueError(f"factor must be a power of 2, got {factor}")
    periodic = box is not None
    boxarr = np.ones(3) if box is None else np.asarray(box, dtype=np.float64)
    pos = np.array(pos, dtype=np.float64)
    mass = np.array(mass, dtype=np.float64)
    aux = {kk: np.array(v, dtype=np.float64) for kk, v in aux.items()}
    target = max(1, len(pos) // factor)
    while len(pos) > target:
        n = len(pos)
        nl = NeighborList(boxarr if periodic else None, k, margin=margin,
                          build_transpose=False)
        idx, _, _, _ = nl.update(pos)
        nn = _nearest(pos, idx, boxarr, periodic)
        # mutual-nearest leaders: nn[i] valid, i < nn[i], and nn[nn[i]] == i
        ii = np.arange(n)
        ok = nn >= 0
        back = np.where(ok, nn[np.where(ok, nn, 0)], -1)   # nn[nn[i]] (guarded)
        leader = ok & (ii < nn) & (back == ii)
        pa = np.where(leader)[0]
        pb = nn[pa]
        if len(pa) == 0:
            break                            # no mutual pairs left (degenerate)
        need = n - target                    # pairs to reach target exactly
        if len(pa) > need:
            # final pass: keep only the `need` closest pairs (most redundant)
            off = pos[pb] - pos[pa]
            if periodic:
                off -= boxarr * np.rint(off / boxarr)
            d2 = (off * off).sum(1)
            sel = np.argsort(d2)[:need]
            pa, pb = pa[sel], pb[sel]
        involved = np.zeros(n, dtype=bool)
        involved[pa] = True
        involved[pb] = True
        keep = np.where(~involved)[0]
        pos, mass, aux = _combine(pos, mass, aux, pa, pb, keep, halve,
                                  boxarr, periodic)
    return pos, mass, aux
