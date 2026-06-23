"""Fixed-radius ball-gather neighbour search (CSR), for the density-based h
modes (measured/target). Unlike the kNN backend it does no k-th-element
selection: each particle gathers every neighbour inside 2*h_i, with h_i the
density-based smoothing length (known a priori from rho0 for 'target', from
the previous iteration's rho for 'measured'). This is what the C force pass
(smReSmooth) actually does - there is no reason to run a kNN for these modes.

Neighbour counts vary per particle, so the lists are CSR (indptr/indices)
rather than the dense (N, k) kNN array. The nsmoothmin floor is enforced by
expanding a particle's radius until its ball holds >= nsmoothmin neighbours
(fail/retry, as in the kNN build) - this replaces the kNN-h-based floor of the
dense path, since no kNN h exists here.
"""
import numpy as np
from numba import get_num_threads, get_thread_id, njit, prange

from .neighbors import _cell_ids, _fill_cells, _gather3


@njit(cache=True, parallel=True)
def _ball_scan(posg, orig, box, ncell, cellptr, todo, rsearch, nsmoothmin,
               count_only, indptr, cnt, hout, fail, indices):
    """One cell-list pass. count_only: just count neighbours of each todo
    particle within rsearch and flag those with < nsmoothmin (fail). Else
    (fill): write each particle's neighbour ids into indices[indptr[i]:].
    Positions are in cell order (posg[p] is particle orig[p]); coordinates
    live in [-L/2, L/2]."""
    nx, ny, nz = ncell[0], ncell[1], ncell[2]
    sx = box[0] / nx
    sy = box[1] / ny
    sz = box[2] / nz
    for t in prange(todo.shape[0]):
        p = todo[t]
        i = orig[p]
        xi = posg[p, 0]
        yi = posg[p, 1]
        zi = posg[p, 2]
        r = rsearch[p]
        cut2 = r * r
        cx = min(max(int((xi + 0.5 * box[0]) / sx), 0), nx - 1)
        cy = min(max(int((yi + 0.5 * box[1]) / sy), 0), ny - 1)
        cz = min(max(int((zi + 0.5 * box[2]) / sz), 0), nz - 1)
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
        c = 0
        w = 0 if count_only else indptr[i]
        for ox in range(xlo, xhi + 1):
            rcx = cx + ox
            ccx = rcx % nx
            shx = box[0] if rcx >= nx else (-box[0] if rcx < 0 else 0.0)
            for oy in range(ylo, yhi + 1):
                rcy = cy + oy
                ccy = rcy % ny
                shy = box[1] if rcy >= ny else (-box[1] if rcy < 0 else 0.0)
                crow = (ccx * ny + ccy) * nz
                for oz in range(zlo, zhi + 1):
                    rcz = cz + oz
                    ccz = rcz % nz
                    shz = box[2] if rcz >= nz else (-box[2] if rcz < 0
                                                    else 0.0)
                    cc = crow + ccz
                    for ptr in range(cellptr[cc], cellptr[cc + 1]):
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
                        if dx * dx + dy * dy + dz * dz >= cut2:
                            continue
                        if count_only:
                            c += 1
                        else:
                            indices[w] = orig[ptr]
                            w += 1
        if count_only:
            if c < nsmoothmin:
                fail[p] = 1
            else:
                cnt[i] = c
                hout[i] = 0.5 * r


@njit(cache=True, parallel=True)
def _ball_store(posg, orig, box, ncell, cellptr, rsearch, nt, chunk,
                tbuf, tcap, loc_start, cnt, overflow):
    """Fused phase-2 walk: each thread walks the cells of its contiguous slot
    chunk *once*, appending neighbour ids into a thread-private buffer
    tbuf[t*tcap:] and recording per-particle (loc_start, cnt). A later scatter
    copies the buffer into the global CSR once indptr is known. This replaces
    the separate count+fill cell-walks (the distance work is done once). If a
    thread exceeds its cap (the a-priori estimate was too small for a sharp
    region) it sets overflow[t]; the caller then falls back to the 2-pass
    build for that rebuild. Coordinates are in cell order, [-L/2, L/2]."""
    nx, ny, nz = ncell[0], ncell[1], ncell[2]
    sx = box[0] / nx
    sy = box[1] / ny
    sz = box[2] / nz
    n = orig.shape[0]
    for t in prange(nt):
        lo = t * chunk
        hi = min(lo + chunk, n)
        base = t * tcap
        w = 0
        for p in range(lo, hi):
            i = orig[p]
            xi = posg[p, 0]
            yi = posg[p, 1]
            zi = posg[p, 2]
            r = rsearch[p]
            cut2 = r * r
            cx = min(max(int((xi + 0.5 * box[0]) / sx), 0), nx - 1)
            cy = min(max(int((yi + 0.5 * box[1]) / sy), 0), ny - 1)
            cz = min(max(int((zi + 0.5 * box[2]) / sz), 0), nz - 1)
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
            start = w
            ovf = False
            for ox in range(xlo, xhi + 1):
                rcx = cx + ox
                ccx = rcx % nx
                shx = box[0] if rcx >= nx else (-box[0] if rcx < 0 else 0.0)
                for oy in range(ylo, yhi + 1):
                    rcy = cy + oy
                    ccy = rcy % ny
                    shy = box[1] if rcy >= ny else (-box[1] if rcy < 0
                                                    else 0.0)
                    crow = (ccx * ny + ccy) * nz
                    for oz in range(zlo, zhi + 1):
                        rcz = cz + oz
                        ccz = rcz % nz
                        shz = box[2] if rcz >= nz else (-box[2] if rcz < 0
                                                        else 0.0)
                        cc = crow + ccz
                        for ptr in range(cellptr[cc], cellptr[cc + 1]):
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
                            if dx * dx + dy * dy + dz * dz >= cut2:
                                continue
                            if w >= tcap:
                                ovf = True
                                break
                            tbuf[base + w] = orig[ptr]
                            w += 1
                        if ovf:
                            break
                    if ovf:
                        break
                if ovf:
                    break
            if ovf:
                overflow[t] = 1
                break
            loc_start[p] = start
            cnt[i] = w - start


@njit(cache=True, parallel=True)
def _ball_scatter(orig, nt, chunk, tbuf, tcap, loc_start, indptr, indices):
    """Copy each thread-private buffer segment into the global CSR row for its
    particle (indptr already prefix-summed from the stored counts)."""
    n = orig.shape[0]
    for t in prange(nt):
        lo = t * chunk
        hi = min(lo + chunk, n)
        base = t * tcap
        for p in range(lo, hi):
            i = orig[p]
            src = base + loc_start[p]
            dst = indptr[i]
            c = indptr[i + 1] - dst
            for e in range(c):
                indices[dst + e] = tbuf[src + e]


@njit(cache=True, parallel=True)
def _ball_esrc(indptr, esrc):
    # per-edge source row: esrc[e] = i for e in indptr[i]:indptr[i+1].
    # Each row owns a disjoint contiguous range -> race-free over rows.
    n = indptr.shape[0] - 1
    for i in prange(n):
        for e in range(indptr[i], indptr[i + 1]):
            esrc[e] = i


@njit(cache=True, parallel=True)
def _ball_transpose_par(indices, esrc, tindptr, tindices, hist):
    """Parallel counting sort of the CSR edge list by target j, mirroring the
    kNN _transpose_fill_par. esrc[e] supplies the source row i (the kNN path
    recovers it as e // k from its dense array; the ball graph is ragged, so
    the row map is materialised). Each thread histograms its chunk of edges,
    per-target prefix sums over threads form disjoint write cursors, and the
    scatter is race-free. Output ordering matches the serial fill (sources
    ascending within each target)."""
    m = indices.shape[0]
    n = tindptr.shape[0] - 1
    nt = hist.shape[0]
    chunk = (m + nt - 1) // nt
    for t in prange(nt):
        h = hist[t]
        for j in range(n):
            h[j] = 0
        lo = t * chunk
        hi = min(lo + chunk, m)
        for e in range(lo, hi):
            h[indices[e]] += 1
    for j in prange(n):
        acc = 0
        for t in range(nt):
            acc += hist[t, j]
        tindptr[j + 1] = acc
    tindptr[0] = 0
    for j in range(n):
        tindptr[j + 1] += tindptr[j]
    for j in prange(n):
        acc = tindptr[j]
        for t in range(nt):
            tmp = hist[t, j]
            hist[t, j] = acc
            acc += tmp
    for t in prange(nt):
        h = hist[t]
        lo = t * chunk
        hi = min(lo + chunk, m)
        for e in range(lo, hi):
            jj = indices[e]
            tindices[h[jj]] = esrc[e]
            h[jj] += 1


def _transpose_csr(indptr, indices, n):
    """Transpose a CSR neighbour graph: for each i, the j whose ball contains
    i. Used to supply the h_j side of the symmetric force walk. Parallel
    counting-sort fill - a serial transpose dominated the rebuild at large N
    (~20M edges at 311k)."""
    m = indices.shape[0]
    esrc = np.empty(m, dtype=np.int32)
    _ball_esrc(indptr, esrc)
    nt = min(get_num_threads(), 64)
    hist = np.empty((nt, n), dtype=np.int64)
    tindptr = np.zeros(n + 1, dtype=np.int64)
    tindices = np.empty(m, dtype=np.int32)
    _ball_transpose_par(indices, esrc, tindptr, tindices, hist)
    return tindptr, tindices


class BallNeighborList:
    """Fixed-radius CSR ball gather for density-based h, with Verlet-skin
    reuse across iterations and cell-order reorder support.

    Each rebuild gathers a *padded* ball (radius 2*h*(1+skin)) so the cached
    lists stay valid for a few iterations of motion + h drift; the density/
    force loops mask to the exact 2*h with the kernel cutoff. The list is
    reused until any particle's slack R_cached - 2*h_now is eaten by twice the
    accumulated displacement (covers both motion and - for structured targets
    - the per-particle h growth as particles cross the density gradient).
    h is floored to keep >= nsmoothmin members (retry-expand at build).
    """

    def __init__(self, box, nsmoothmin=32, skin=0.3, build_transpose=True):
        self.box = np.asarray(box, dtype=np.float64)
        self.nsmoothmin = nsmoothmin
        self.skin = skin
        self.build_transpose = build_transpose
        self.nbuilds = 0
        self.indptr = self.indices = self.tip = self.tii = None
        self.hout = None          # current (possibly reused) kernel h
        self.hfloor = None        # per-particle floor level from last build
        self.rcached = None       # per-particle gather radius from last build
        self._disp = 0.0          # max accumulated displacement since build
        self.pdisp = None

    def track_displacements(self, disp):
        self.pdisp = disp

    def _slack_min(self, h):
        # smallest per-particle headroom R_cached - 2*h_now (negative once a
        # particle's kernel has grown past its cached gather radius)
        return float(np.min(self.rcached - 2.0 * h))

    def needs_rebuild(self, h_in):
        """Per-particle Verlet-skin staleness test. Particle i's cached ball
        (radius R_cached_i) stays valid while every neighbour stays inside it;
        the relative approach of i and any neighbour k is bounded by
        pdisp[i] + pdisp[k] <= pdisp[i] + dmax. So a rebuild is needed iff some
        i has  pdisp[i] + dmax > slack_i = R_cached_i - 2*h_now_i.

        This is strictly tighter than the old global test (2*dmax >
        min_i slack_i): there, one fast interface mover (large dmax) coupled to
        an unrelated small-slack particle forced a full rebuild. Pairing each
        particle's own slack with its own displacement cuts the rebuild rate
        in sharp-gradient regions. Falls back to the global bound when no
        per-particle displacement array is wired (standalone use)."""
        if self.indptr is None:
            return True
        hout = np.maximum(h_in, self.hfloor)
        if self.pdisp is None:
            return 2.0 * self._disp > self._slack_min(hout)
        dmax = float(self.pdisp.max())
        # rebuild iff max_i (pdisp[i] + 2*h_now_i - R_cached_i) + dmax > 0
        worst = float(np.max(self.pdisp + 2.0 * hout - self.rcached))
        return worst + dmax > 0.0

    def record_move(self, max_step):
        self._disp += max_step

    def update(self, pos, h_in, rebuild=None):
        """Return (indptr, indices, hout, t_indptr, t_indices). Rebuilds the
        ball when the cached lists are stale, else reuses them with h updated
        to the current density-based value (floored). `rebuild` lets the caller
        pass a precomputed needs_rebuild() decision: the driver evaluates it
        once *before* the reorder (which permutes pos/h_in but not the cached
        rcached/pdisp), so update() must not re-test against the now-misaligned
        state. None -> decide internally (standalone use, no reorder)."""
        if rebuild is None:
            rebuild = self.needs_rebuild(h_in)
        if rebuild:
            self._build(pos, h_in)
        else:
            self.hout = np.maximum(h_in, self.hfloor)
        return self.indptr, self.indices, self.hout, self.tip, self.tii

    def _fused_gather(self, posg, orig, box, ncell, cellptr, rfull, cnt_unpad,
                      indptr):
        """Single-walk padded gather: store neighbour ids into per-thread
        buffers, prefix-sum the counts into indptr, then scatter into a tight
        CSR. Returns the indices array, or None if a thread overflowed its
        a-priori cap (caller falls back to the 2-pass build). cnt_unpad holds
        the phase-1 unpadded counts (orig order); the padded ball is bounded
        by the volume ratio (1+skin)^3 with headroom."""
        n = orig.shape[0]
        nt = min(get_num_threads(), 64)
        chunk = (n + nt - 1) // nt
        # per-slot padded estimate from the unpadded count; size each thread's
        # cap to its chunk's estimated edges with generous headroom
        est_slot = cnt_unpad[orig] * ((1.0 + self.skin) ** 3) * 1.5 + 32.0
        thread_est = np.zeros(nt, dtype=np.int64)
        for t in range(nt):
            lo = t * chunk
            hi = min(lo + chunk, n)
            if lo < hi:
                thread_est[t] = est_slot[lo:hi].sum()
        tcap = int(thread_est.max()) + 64 if n else 1
        try:
            tbuf = np.empty(nt * tcap, dtype=np.int32)
        except (MemoryError, ValueError):
            return None
        loc_start = np.zeros(n, dtype=np.int64)
        cntf = np.zeros(n, dtype=np.int64)
        overflow = np.zeros(nt, dtype=np.uint8)
        _ball_store(posg, orig, box, ncell, cellptr, rfull, nt, chunk,
                    tbuf, tcap, loc_start, cntf, overflow)
        if overflow.any():
            return None
        np.cumsum(cntf, out=indptr[1:])
        indices = np.empty(int(indptr[-1]), dtype=np.int32)
        _ball_scatter(orig, nt, chunk, tbuf, tcap, loc_start, indptr, indices)
        return indices

    def _build(self, pos, h_in):
        n = len(pos)
        box = self.box
        ncell = np.maximum(
            1, (box / (5.0 * box.prod() / n) ** (1.0 / 3.0)).astype(np.int64))
        cid = np.empty(n, dtype=np.int64)
        _cell_ids(pos, box, ncell, cid)
        counts = np.bincount(cid, minlength=int(ncell.prod()))
        cellptr = np.zeros(len(counts) + 1, dtype=np.int64)
        np.cumsum(counts, out=cellptr[1:])
        orig = np.empty(n, dtype=np.int32)
        _fill_cells(cid, cellptr, orig)
        posg = np.empty_like(pos)
        _gather3(pos, orig, posg)

        rcap = 0.5 * float(np.sqrt((box * box).sum()))
        # phase 1: find the floored kernel h - gather within 2*h, expanding
        # where < nsmoothmin members (rsearch indexed by cell-order slot p)
        rsearch = np.minimum(2.0 * h_in[orig], rcap)
        cnt = np.zeros(n, dtype=np.int64)
        hout = np.empty(n)
        fail = np.zeros(n, dtype=np.uint8)
        todo = np.arange(n, dtype=np.int64)
        _di = np.empty(0, dtype=np.int32)
        _dip = np.zeros(1, dtype=np.int64)
        for _ in range(60):
            _ball_scan(posg, orig, box, ncell, cellptr, todo, rsearch,
                       self.nsmoothmin, True, _dip, cnt, hout, fail, _di)
            todo = todo[fail[todo] == 1]
            if todo.size == 0:
                break
            fail[todo] = 0
            rsearch[todo] = np.minimum(rsearch[todo] * 1.6, rcap)
        else:
            raise RuntimeError('ball gather did not reach nsmoothmin')
        # phase 2: gather the *padded* ball (2*h_out*(1+skin)) for reuse;
        # density/force mask back to 2*h_out via the kernel cutoff. Fused
        # store+scatter (one cell-walk) with a 2-pass fallback on overflow.
        rfull = np.minimum(rsearch * (1.0 + self.skin), rcap)
        indptr = np.zeros(n + 1, dtype=np.int64)
        indices = self._fused_gather(posg, orig, box, ncell, cellptr,
                                     rfull, cnt, indptr)
        if indices is None:                # overflow -> robust 2-pass build
            cntf = np.zeros(n, dtype=np.int64)
            _ball_scan(posg, orig, box, ncell, cellptr,
                       np.arange(n, dtype=np.int64), rfull, 0, True, _dip,
                       cntf, np.empty(n), np.zeros(n, dtype=np.uint8), _di)
            indptr = np.zeros(n + 1, dtype=np.int64)
            np.cumsum(cntf, out=indptr[1:])
            indices = np.empty(int(indptr[-1]), dtype=np.int32)
            _ball_scan(posg, orig, box, ncell, cellptr,
                       np.arange(n, dtype=np.int64), rfull, 0, False, indptr,
                       cntf, np.empty(n), np.zeros(n, dtype=np.uint8), indices)

        if self.build_transpose:
            tip, tii = _transpose_csr(indptr, indices, n)
        else:
            tip = np.zeros(n + 1, dtype=np.int64)
            tii = np.empty(0, dtype=np.int32)
        self.indptr, self.indices, self.tip, self.tii = indptr, indices, \
            tip, tii
        self.hout = hout
        # floor level: where the retry expanded h above the density value
        self.hfloor = np.where(hout > h_in * 1.001, hout, 0.0)
        # cached gather radius per particle (orig order): rfull/2 in h-units
        rc = np.empty(n)
        rc[orig] = rfull
        self.rcached = rc
        self._disp = 0.0
        if self.pdisp is not None:
            self.pdisp[:] = 0.0
        self.nbuilds += 1
