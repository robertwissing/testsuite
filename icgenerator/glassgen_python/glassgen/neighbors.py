"""Neighbor search: cKDTree kNN lists with Verlet-skin reuse.

The candidate list per particle is the k = nsmooth + margin nearest neighbors
(self excluded). h_i is half the distance to the nsmooth-th nearest, matching
gasoline where fBall encloses nSmooth particles (including self).

Between tree rebuilds the candidate set is reused: pair distances are always
recomputed from current positions, and the kernel cutoff (W = 0 at q >= 2)
makes stale candidates beyond 2h harmless. The list is rebuilt when particles
have moved far enough that a non-candidate could have entered some ball:
2 * max accumulated displacement > min over i of (d_far_i - 2h_i at build).
"""
import numpy as np
from numba import get_num_threads, get_thread_id, njit, prange
from scipy.spatial import cKDTree


@njit(cache=True)
def _select_d(d, n, kk):
    """Partial quickselect: afterwards the kk smallest of d[:n] occupy
    d[:kk] (unordered)."""
    lo = 0
    hi = n - 1
    while lo < hi:
        pivot = d[(lo + hi) >> 1]
        i1 = lo
        i2 = hi
        while i1 <= i2:
            while d[i1] < pivot:
                i1 += 1
            while d[i2] > pivot:
                i2 -= 1
            if i1 <= i2:
                d[i1], d[i2] = d[i2], d[i1]
                i1 += 1
                i2 -= 1
        if kk - 1 <= i2:
            hi = i2
        elif kk - 1 >= i1:
            lo = i1
        else:
            break


@njit(cache=True, parallel=True)
def _update_h(pos, idx, box, periodic, nsmooth, h, ws):
    # positions must be wrapped to [-L/2, L/2] when periodic (see
    # density_loop); ws is a per-thread distance workspace, (nthreads, k)
    n, k = idx.shape
    bx, by, bz = box[0], box[1], box[2]
    hbx, hby, hbz = 0.5 * bx, 0.5 * by, 0.5 * bz
    for i in prange(n):
        d2 = ws[get_thread_id()]
        for jj in range(k):
            j = idx[i, jj]
            dx = pos[i, 0] - pos[j, 0]
            dy = pos[i, 1] - pos[j, 1]
            dz = pos[i, 2] - pos[j, 2]
            if periodic:
                if dx > hbx:
                    dx -= bx
                elif dx < -hbx:
                    dx += bx
                if dy > hby:
                    dy -= by
                elif dy < -hby:
                    dy += by
                if dz > hbz:
                    dz -= bz
                elif dz < -hbz:
                    dz += bz
            d2[jj] = dx * dx + dy * dy + dz * dz
        # ball of nsmooth particles counting self -> radius is the
        # (nsmooth-1)-th non-self candidate: the largest of the nsmooth-1
        # smallest distances
        _select_d(d2, k, nsmooth - 1)
        far = d2[0]
        for u in range(1, nsmooth - 1):
            if d2[u] > far:
                far = d2[u]
        h[i] = 0.5 * np.sqrt(far)


@njit(cache=True, parallel=True)
def _transpose_fill_par(flat, k, indptr, indices, hist):
    """Parallel counting sort of the kNN edge list by target. Each thread
    histograms its chunk of edges (pass 1), per-target exclusive prefix
    sums over threads turn the histograms into disjoint write cursors, and
    pass 2 scatters race-free. Output ordering is identical to the serial
    fill (edges ascending within each target)."""
    m = flat.shape[0]
    n = indptr.shape[0] - 1
    nt = hist.shape[0]
    chunk = (m + nt - 1) // nt
    for t in prange(nt):
        h = hist[t]
        for j in range(n):
            h[j] = 0
        lo = t * chunk
        hi = min(lo + chunk, m)
        for e in range(lo, hi):
            h[flat[e]] += 1
    for j in prange(n):
        acc = 0
        for t in range(nt):
            acc += hist[t, j]
        indptr[j + 1] = acc
    indptr[0] = 0
    for j in range(n):
        indptr[j + 1] += indptr[j]
    for j in prange(n):
        acc = indptr[j]
        for t in range(nt):
            tmp = hist[t, j]
            hist[t, j] = acc
            acc += tmp
    for t in prange(nt):
        h = hist[t]
        lo = t * chunk
        hi = min(lo + chunk, m)
        for e in range(lo, hi):
            jj = flat[e]
            indices[h[jj]] = e // k
            h[jj] += 1


def _transpose_graph(idx):
    """CSR of the transposed kNN graph: for each i, the j with i in idx[j].
    Parallel counting-sort fill; a serial fill (or an argsort) here
    dominates the whole iteration cost at large N."""
    n, k = idx.shape
    flat = idx.ravel()
    # cap the histogram workspace (nt x n int64); 64 threads saturate the
    # memory bandwidth this pass is bound by anyway
    nt = min(get_num_threads(), 64)
    hist = np.empty((nt, n), dtype=np.int64)
    indptr = np.zeros(n + 1, dtype=np.int64)
    indices = np.empty(flat.shape[0], dtype=np.int32)
    _transpose_fill_par(flat, k, indptr, indices, hist)
    return indptr, indices


@njit(cache=True, parallel=True)
def _count_stale(disp, idx, skin):
    """Per-particle Verlet staleness: particle i's candidate list is stale
    when its own accumulated displacement plus the largest displacement
    among its candidates exceeds its own skin (d_far_i - 2h_i at build).
    Unlike a global max-step vs min-skin test, this keeps regions of very
    different h (e.g. a 360:1 collapse sphere vs its ambient) from
    invalidating each other's lists."""
    n, k = idx.shape
    nstale = 0
    for i in prange(n):
        di = disp[i]
        si = skin[i]
        m = 0.0
        for jj in range(k):
            dj = disp[idx[i, jj]]
            if dj > m:
                m = dj
                if di + m > si:
                    break
        if di + m > si:
            nstale += 1
    return nstale


@njit(cache=True, parallel=True)
def _cell_ids(pos, box, ncell, cid):
    n = pos.shape[0]
    nx, ny, nz = ncell[0], ncell[1], ncell[2]
    sx = box[0] / nx
    sy = box[1] / ny
    sz = box[2] / nz
    for i in prange(n):
        cx = int((pos[i, 0] + 0.5 * box[0]) / sx)
        cy = int((pos[i, 1] + 0.5 * box[1]) / sy)
        cz = int((pos[i, 2] + 0.5 * box[2]) / sz)
        cx = min(max(cx, 0), nx - 1)
        cy = min(max(cy, 0), ny - 1)
        cz = min(max(cz, 0), nz - 1)
        cid[i] = (cx * ny + cy) * nz + cz


@njit(cache=True, parallel=True)
def _gather3(src, perm, dst):
    for i in prange(perm.shape[0]):
        j = perm[i]
        dst[i, 0] = src[j, 0]
        dst[i, 1] = src[j, 1]
        dst[i, 2] = src[j, 2]


@njit(cache=True)
def _fill_cells(cid, indptr, indices):
    cursor = indptr[:-1].copy()
    for m in range(cid.shape[0]):
        c = cid[m]
        indices[cursor[c]] = m
        cursor[c] += 1


@njit(cache=True)
def _select_k(d, j, n, kk):
    """Partial quickselect on d[:n] (j permuted alongside): afterwards the
    kk smallest distances occupy d[:kk] (unordered)."""
    lo = 0
    hi = n - 1
    while lo < hi:
        pivot = d[(lo + hi) >> 1]
        i1 = lo
        i2 = hi
        while i1 <= i2:
            while d[i1] < pivot:
                i1 += 1
            while d[i2] > pivot:
                i2 -= 1
            if i1 <= i2:
                d[i1], d[i2] = d[i2], d[i1]
                j[i1], j[i2] = j[i2], j[i1]
                i1 += 1
                i2 -= 1
        if kk - 1 <= i2:
            hi = i2
        elif kk - 1 >= i1:
            lo = i1
        else:
            break


@njit(cache=True, parallel=True)
def _grid_knn(posg, orig, box, ncell, cellptr, todo, rsearch,
              nsmooth, k, bufd, bufj, idx, h, dfar2, fail):
    """Exact kNN via periodic cell lists. posg holds positions gathered into
    cell order (posg[p] is particle orig[p]); candidates within a cell are
    then contiguous in memory, and the periodic image shift is resolved per
    cell offset instead of per pair. todo/fail/rsearch are indexed by the
    cell-order slot p, outputs (idx, h, dfar2) by the original particle id.
    Particles whose search sphere holds < k neighbors are flagged in fail
    (caller retries them with a larger rsearch)."""
    nx, ny, nz = ncell[0], ncell[1], ncell[2]
    sx = box[0] / nx
    sy = box[1] / ny
    sz = box[2] / nz
    B = bufd.shape[1]
    for t in prange(todo.shape[0]):
        p = todo[t]
        i = orig[p]
        tid = get_thread_id()
        bd = bufd[tid]
        bj = bufj[tid]
        xi = posg[p, 0]
        yi = posg[p, 1]
        zi = posg[p, 2]
        r = rsearch[p]
        cut2 = r * r
        # home cell (positions live in [-L/2, L/2], cells in [0, L))
        cx = min(max(int((xi + 0.5 * box[0]) / sx), 0), nx - 1)
        cy = min(max(int((yi + 0.5 * box[1]) / sy), 0), ny - 1)
        cz = min(max(int((zi + 0.5 * box[2]) / sz), 0), nz - 1)
        # offsets per axis; if the search range spans the axis, enumerate
        # every cell exactly once instead (avoids periodic double-visits)
        # and fall back to per-pair min-image on that axis
        nrx = int(np.ceil(r / sx))
        nry = int(np.ceil(r / sy))
        nrz = int(np.ceil(r / sz))
        wrapx = 2 * nrx + 1 >= nx
        wrapy = 2 * nry + 1 >= ny
        wrapz = 2 * nrz + 1 >= nz
        if wrapx:
            cx = 0
            xlo, xhi = 0, nx - 1
        else:
            xlo, xhi = -nrx, nrx
        if wrapy:
            cy = 0
            ylo, yhi = 0, ny - 1
        else:
            ylo, yhi = -nry, nry
        if wrapz:
            cz = 0
            zlo, zhi = 0, nz - 1
        else:
            zlo, zhi = -nrz, nrz
        cnt = 0
        for ox in range(xlo, xhi + 1):
            rcx = cx + ox
            ccx = rcx % nx
            shx = 0.0
            if rcx >= nx:
                shx = box[0]
            elif rcx < 0:
                shx = -box[0]
            for oy in range(ylo, yhi + 1):
                rcy = cy + oy
                ccy = rcy % ny
                shy = 0.0
                if rcy >= ny:
                    shy = box[1]
                elif rcy < 0:
                    shy = -box[1]
                crow = (ccx * ny + ccy) * nz
                for oz in range(zlo, zhi + 1):
                    rcz = cz + oz
                    ccz = rcz % nz
                    shz = 0.0
                    if rcz >= nz:
                        shz = box[2]
                    elif rcz < 0:
                        shz = -box[2]
                    c = crow + ccz
                    for ptr in range(cellptr[c], cellptr[c + 1]):
                        if ptr == p:
                            continue
                        dx = xi - posg[ptr, 0] - shx
                        dy = yi - posg[ptr, 1] - shy
                        dz = zi - posg[ptr, 2] - shz
                        if wrapx:
                            dx -= box[0] * np.rint(dx / box[0])
                        if wrapy:
                            dy -= box[1] * np.rint(dy / box[1])
                        if wrapz:
                            dz -= box[2] * np.rint(dz / box[2])
                        d2 = dx * dx + dy * dy + dz * dz
                        if d2 >= cut2:
                            continue
                        if cnt == B:
                            # buffer full: keep the k nearest so far and
                            # tighten the cutoff to their max
                            _select_k(bd, bj, B, k)
                            cnt = k
                            far = 0.0
                            for u in range(k):
                                if bd[u] > far:
                                    far = bd[u]
                            cut2 = far
                            if d2 >= cut2:
                                continue
                        bd[cnt] = d2
                        bj[cnt] = ptr
                        cnt += 1
        if cnt < k:
            fail[p] = 1
            continue
        _select_k(bd, bj, cnt, k)
        far = 0.0
        for u in range(k):
            if bd[u] > far:
                far = bd[u]
        dfar2[i] = far
        for u in range(k):
            idx[i, u] = orig[bj[u]]
        # ball of nsmooth particles counting self -> radius is the
        # (nsmooth-1)-th non-self neighbor: largest of the nsmooth-1
        # smallest (allocation-free, destroys bd/bj pairing - idx is
        # already written)
        _select_d(bd, k, nsmooth - 1)
        hf = bd[0]
        for u in range(1, nsmooth - 1):
            if bd[u] > hf:
                hf = bd[u]
        h[i] = 0.5 * np.sqrt(hf)


def cell_order(pos, box):
    """Permutation sorting particles by cell id on the same grid the
    neighbor build uses. Reordering all particle arrays with it before a
    rebuild makes every gather in the SPH loops cache-local. Returns None
    when the grid would be degenerate (tiny problems)."""
    n = len(pos)
    box = np.asarray(box, dtype=np.float64)
    ncell = np.maximum(
        1, (box / (5.0 * box.prod() / n) ** (1.0 / 3.0)).astype(np.int64))
    if (ncell < 3).any():
        return None
    cid = np.empty(n, dtype=np.int64)
    _cell_ids(pos, box, ncell, cid)
    counts = np.bincount(cid, minlength=int(ncell.prod()))
    indptr = np.zeros(len(counts) + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    perm = np.empty(n, dtype=np.int32)
    _fill_cells(cid, indptr, perm)
    return perm


class NeighborList:
    """kNN candidate lists with Verlet-skin reuse across iterations."""

    def __init__(self, box, nsmooth, margin=8, build_transpose=True):
        self.box = None if box is None else np.asarray(box, dtype=np.float64)
        self.periodic = box is not None
        self.nsmooth = nsmooth
        self.k = nsmooth + margin
        self.build_transpose = build_transpose
        self.idx = None
        self.h = None
        self.rev_indptr = None
        self.rev_indices = None
        self._disp = 0.0  # max accumulated displacement since build
        self._skin = 0.0  # min_i (d_far_i - 2h_i) at build
        self.nbuilds = 0
        self.last_build = None  # 'grid' or 'tree'
        self._hws = None  # per-thread distance workspace for _update_h
        self.skin_pp = None  # per-particle skin (d_far - 2h at build)
        self.pdisp = None  # caller-attached per-particle displacements
        self._stale = None  # memoized needs_rebuild verdict

    def _boxarr(self):
        if self.periodic:
            return self.box
        return np.ones(3)  # unused when periodic is False

    def track_displacements(self, disp):
        """Attach a per-particle displacement accumulator (the caller adds
        each step's per-particle magnitude, e.g. via stepper.apply_step).
        Enables the per-particle staleness test; reset at every rebuild."""
        self.pdisp = disp

    def needs_rebuild(self):
        """True when the next update() will rebuild the candidate lists
        (lets callers reorder their arrays for locality beforehand)."""
        if self.idx is None:
            return True
        if self.pdisp is not None and self.skin_pp is not None:
            if self._stale is None:
                self._stale = bool(
                    _count_stale(self.pdisp, self.idx, self.skin_pp))
            return self._stale
        return 2.0 * self._disp > self._skin

    def update(self, pos):
        """Refresh candidate lists (if needed) and h for current positions."""
        if self.needs_rebuild():
            self._build(pos)
        else:
            if self._hws is None:
                self._hws = np.empty((get_num_threads(), self.k))
            _update_h(pos, self.idx, self._boxarr(), self.periodic,
                      self.nsmooth, self.h, self._hws)
        return self.idx, self.h, self.rev_indptr, self.rev_indices

    def _build(self, pos):
        if self.periodic and self._build_grid(pos):
            self.last_build = 'grid'
        else:
            self._build_tree(pos)
            self.last_build = 'tree'
        self._disp = 0.0
        if self.pdisp is not None:
            self.pdisp[:] = 0.0
        self._stale = False
        if self.build_transpose:
            self.rev_indptr, self.rev_indices = _transpose_graph(self.idx)
        else:
            self.rev_indptr = np.zeros(len(self.idx) + 1, dtype=np.int64)
            self.rev_indices = np.empty(0, dtype=np.int32)
        self.nbuilds += 1

    def _build_tree(self, pos):
        if self.periodic:
            # cKDTree wants coordinates in [0, L); the box is centered on 0
            shifted = np.mod(pos + 0.5 * self.box, self.box)
            tree = cKDTree(shifted, boxsize=self.box)
            dist, idx = tree.query(shifted, k=self.k + 1, workers=-1)
        else:
            tree = cKDTree(pos)
            dist, idx = tree.query(pos, k=self.k + 1, workers=-1)
        # drop self (column 0; query returns sorted distances, self first);
        # the ball holds nsmooth particles counting self, so its radius is
        # the (nsmooth-1)-th non-self distance
        self.idx = np.ascontiguousarray(idx[:, 1:], dtype=np.int32)
        if self.h is None:
            self.h = np.empty(len(pos))
        self.h[:] = 0.5 * dist[:, self.nsmooth - 1]
        self.skin_pp = dist[:, -1] - 2.0 * self.h
        self._skin = float(np.min(self.skin_pp))

    def _build_grid(self, pos):
        """Numba cell-list kNN build; returns False when the problem is too
        small for a sensible grid (caller falls back to the kd-tree)."""
        n = len(pos)
        box = self.box
        if n <= self.k + 1:
            return False
        # target ~5 particles per cell
        ncell = np.maximum(
            1, (box / (5.0 * box.prod() / n) ** (1.0 / 3.0)).astype(np.int64))
        if (ncell < 3).any():
            return False
        cid = np.empty(n, dtype=np.int64)
        _cell_ids(pos, box, ncell, cid)
        counts = np.bincount(cid, minlength=int(ncell.prod()))
        cellptr = np.zeros(len(counts) + 1, dtype=np.int64)
        np.cumsum(counts, out=cellptr[1:])
        orig = np.empty(n, dtype=np.int32)
        _fill_cells(cid, cellptr, orig)
        # gather positions into cell order: the candidate scans in
        # _grid_knn then stream through memory instead of hopping randomly
        posg = np.empty_like(pos)
        _gather3(pos, orig, posg)

        # search radius: expected k-th neighbor distance plus slack; the
        # fail/retry loop below backstops any underestimate (e.g. stale h
        # after the large early-sweep moves). Capped at the min-image
        # envelope: beyond it the scan would degenerate to all-pairs.
        scale = 1.25 * (self.k / self.nsmooth) ** (1.0 / 3.0)
        rcap = 0.5 * float(np.sqrt((box * box).sum()))
        if self.h is not None:
            rsearch = np.minimum((2.0 * self.h[orig]) * scale, rcap)
        else:
            ndens = n / box.prod()
            rk = (3.0 * self.k / (4.0 * np.pi * ndens)) ** (1.0 / 3.0)
            rsearch = np.full(n, min(scale * rk, rcap))

        nt = get_num_threads()
        B = max(4 * self.k, 256)
        bufd = np.empty((nt, B))
        bufj = np.empty((nt, B), dtype=np.int64)
        idx = np.empty((n, self.k), dtype=np.int32)
        hbuf = np.empty(n)
        dfar2 = np.empty(n)
        fail = np.zeros(n, dtype=np.uint8)
        todo = np.arange(n, dtype=np.int64)
        for _ in range(60):
            _grid_knn(posg, orig, box, ncell, cellptr, todo, rsearch,
                      self.nsmooth, self.k, bufd, bufj, idx, hbuf, dfar2,
                      fail)
            todo = todo[fail[todo] == 1]
            if todo.size == 0:
                break
            fail[todo] = 0
            rsearch[todo] = np.minimum(rsearch[todo] * 1.6, rcap)
        else:
            raise RuntimeError('grid kNN did not converge')
        self.idx = idx
        if self.h is None:
            self.h = np.empty(n)
        self.h[:] = hbuf
        self.skin_pp = np.sqrt(dfar2) - 2.0 * hbuf
        self._skin = float(np.min(self.skin_pp))
        return True

    def record_move(self, max_step):
        """Accumulate the largest displacement applied this iteration."""
        self._disp += max_step
        self._stale = None  # displacements changed; re-evaluate staleness
