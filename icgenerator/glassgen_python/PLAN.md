# Python glass-relaxation IC generator (`glassgen`)

## Status (2026-06-22) — Wendland C2 kernel + nsmoothmin pairing-instability fix

**Kernel switched M4 cubic spline → Wendland C2 (3D); M4 removed.**
`kernels.py` now provides only `w_wc2`/`dw_wc2` (+ `W0 = 21/16`, the central
shape value). The Wendland normalization `21/16` is baked into the kernel
*shape*, so every density/force loop keeps the gasoline `fNorm = 1/(pi h^3)`
machinery byte-unchanged — only the shape and the explicit `W(0)` self-terms
(`s = mass[i] * W0`, Q0 self `mass[i]/rho[i] * W0`) differ. Call sites renamed
in `sph_core.py`/`diagnostics.py`. 60/60 tests pass; the four e2e quality tests
moved from `nsmooth=32` → `64` (Wendland's proper operating point, already the
`relax`/benchmark default) with tolerances matching its honest accuracy.

**Why: the `nsmoothmin`-near-`nSmooth` density runaway (root-caused this
session).** On the **ball (measured/target) path**, density is computed with the
density-based h and the next h is recomputed from rho (`h=(m·hfact3/rho)^⅓`),
closing a **rho→h→rho feedback loop** (the knn path is immune — its density uses
the geometric kNN h, no feedback). At the consistent fixed point each particle's
2h ball holds ~`nSmooth` neighbours *distributed around* nSmooth, so a floor near
nSmooth (e.g. 31 vs nSmooth=32) fell in the middle of that distribution and fired
the ball-expansion (`ball.py` ×1.6 radius → ~4× volume → **~131 neighbours**) on
~40% of particles. With the **M4 cubic spline (pairing instability onset ~60
neighbours)** that collective over-expansion clumped particles → rho spike → h↓
→ rho_sph↑ runaway (traced rho 4 → 15,849 → **55,439**, then polish froze the
wreckage). Confirmed by sweep: frac with >60 nbrs went 0%/5%/40% for
nsmoothmin 8/24/31, err 0.004/0.004/**29.6**. **Wendland C2 has no pairing
instability**, so the expansion can't clump: nsmoothmin=31 → err 0.10 (no
runaway), and at nSmooth=64 nsmoothmin 8 vs 32 are identical (0.0099/0.0100).
**`nsmoothmin` is no longer a stability-critical parameter.**

**Tradeoff:** Wendland needs more neighbours for equal density accuracy —
uniform-glass err ~0.010 at nSmooth=64 (vs M4 ~0.002), ~0.004 at nSmooth=100
(≈ old M4@32). nSmooth=64 is the default/middle ground; bump to ~100 for
M4-parity.

**Other changes this session (all in `benchmark.py`/`relax.py`/`progressive.py`):**
- `--no-stop` now removes **all** early-stop criteria (zeroes `stop_frac`,
  `icr0_stop`, `active_stop`, `polish_iter`) so a run goes the full `--max-iter`
  budget. Previously it only zeroed `stop_frac`, so the `icr0_stop=3e-3` finish
  still cut runs short.
- Fixed the **coarse-skip bug** in `default_levels`/`cascade_levels`
  (`nx_full ≤ nx_coarse`): the single `'full'` level now honours the passed
  `icr0_stop`/`polish_iter` instead of silently inheriting `relax()`'s 3e-3
  default (why `--icr0-stop 0` was being ignored).
- **Dual diagnostic**: `RelaxResult` snapshots the glass at polish entry
  (`pos/h/rho_prepolish`, `polish_iter0`); benchmark always prints both
  `non-polish @it…` and `polish @it…` blocks, exposing polish-phase degradation
  (e.g. the measured icr0-collapse freeze: force jumped 1.2→10 across polish on
  the M4 360:1 case — a separate finding; the non-knn variants enter polish
  *before* converging and the asymmetric polish rates collapse icr0 to ~0,
  freezing an unconverged glass, whereas knn enters polish already converged).
- `nsmoothmin` default 8 → **32**; guard added: `if nsmoothmin >= nsmooth:
  nsmoothmin = nsmooth-1` (NB this guard only prevents the literal `≥`
  degeneracy — with M4 a near-nSmooth floor was still unsafe; with Wendland it's
  moot since there's no pairing instability).

## Status (2026-06-12)

Steps 1–8 implemented and locally verified: all 27 tests in `tests/` pass
(`python -m pytest tests/`), and the CLI round-tripped a real pre-IC
(`test_cases/shocktube/datafiles/shocktube64_latticeRyu_1b.00000`, 157k gas
particles, anisotropic periodic box from the .param): non-position columns
bit-identical, velocities zeroed, rho/hsmooth updated, positions wrapped,
`-oi` snapshots written. Note: open (non-periodic) domains have no
confinement, so the open path is smoke-tested for finiteness only
(`test_open_boundary_smoke`); meaningful glass generation needs a periodic
box.

**gdicgen comparison DONE (2026-06-13, `comparison/COMPARISON_REPORT.md`)**:
both codes driven with the exact runtest ICCONFIG, symmetric gdforce, quality
recomputed with glassgen on both glasses. Two cases — currentsheet64_rand
(37k, uniform profile 1) and mhdcollapse50_randrd360 (311k, tanh sphere 360:1,
density table matching the C analytic profile to 3e-5). Result: functionally
equivalent — Q1→1 on both, radial ρ(r) profiles lie on top of each other and
on the target, density-error/Q1/NN histograms overlap in shape; Python 6–10×
faster wall clock (caveat: C is a verbose debug build). One systematic gap:
the C glass reaches a ~1.5–3× tighter density-error floor. Root cause pinned
in source — C `DenDVDX` (smoothfcn.c:6080-6082) overrides h to density-based
`(m·hfact3/ρ)^⅓` before the force resmooth; **Python v1 uses kNN h**, whose
geometric shot noise sets a residual floor (Python denserr plateaus, so it's
an operator floor not under-convergence). Python's zero-order force error E0
is ~2× *lower* than C's in both cases. This answers the kNN-vs-density-h
question PLAN step 5 deferred. **Recommended follow-up**: add density-based h
(after density_loop set h=(m·hfact3/ρ)^⅓, use for force resmooth) — expected
to close the floor gap while keeping the speed + lower E0.

**Density-based h option + benchmark (2026-06-13)**: `relax(hmode=...)` /
CLI `--h-mode` adds `knn` (default v1), `measured` (h=(m·hfact3/ρ_sph)^⅓, the
C DenDVDX choice) and `target` (h=(m·hfact3/ρ0)^⅓, analytic noise-free) — only
the **force kernel + move step** use it; density and pP stay on kNN h, exactly
matching the *inconsistent* hacked-gasoline path (verified: ICGenerator,
smoothfcn.c:5804, gathers a fixed BALL2=4·hdens² ball — a `smReSmooth`, not a
kNN — and reuses the DenDVDX kNN-h density for pP). A `nsmoothmin` floor
(default 8) caps how small h may get (`h ≥ h_knn·(nsmoothmin/nSmooth)^⅓`) so an
analytic-h particle sitting in a sparse region whose ρ0 says "dense" still
gathers real neighbours instead of an empty ball. Findings: on **uniform**
targets density-based h clearly tightens the floor (currentsheet density CoV
0.0077 knn → 0.0050 measured / 0.0057 target, toward C's 0.0022) at the cost
of higher E0 (the density-vs-order tradeoff the C source itself comments on,
smoothfcn.c:5826). On a **sharp 360:1 interface** (mhd50) it improves the bulk
median/p95 (target best) but *worsens* the interface max (target err_max 0.45
vs knn 0.25 — the large h-jump across the discontinuity); `measured` is the
safe pick for discontinuous targets, `target` for smooth ones. Runtime: as
implemented density-based h has **no per-iteration win** (still does the kNN
build + reuse; slightly slower from the extra cbrt and, for target-h, larger
outer-region balls + more rebuilds); the net effect is via convergence
iteration count (measured converged faster on both mhd12 and mhd50). The real
neighbour-search win is unrealized — see Future tasks (target-h ball gather).

`benchmark.py` (package root) is the easy perf+accuracy harness: cases
mhdcollapse12 (synth low-res structured), sedov64 (synth uniform),
mhdcollapse50, shocktube256 (real pre-ICs; targets ported from each .param's
dICdens* profiles 1/2/3); variants knn / knn-1side / measured / target /
target-1side; reports ms/it, Mparticle·it/s throughput, density err + Q1, and
optional `--with-c` to time gasoline.gdicgen. `spatial.py` renders the radial
ρ(r) profile + a particle-distribution slice for any glass.

**Performance rework (2026-06-12, after 256³ = 16.7M-particle trial)**: the
periodic neighbor build no longer uses cKDTree but a numba cell-grid kNN
(`_grid_knn`): positions gathered into cell order (streaming candidate
scans), periodic image shift resolved per cell offset, per-particle
fail/retry for sparse regions; kd-tree kept for open domains and degenerate
grids (< 3 cells/axis, exercised as exact-match fallback). The transposed
kNN graph is filled by a parallel counting sort (`_transpose_fill_par`) —
the old serial fill was 42 s at 16.7M, silently dominating every rebuild.
The force loop is now a top-level @njit(cache=True) with integer flavour
dispatch (flavours.ig_pair) instead of a per-flavour closure factory, so it
JIT-compiles once ever instead of ~30-60 s on every process start.
Measured at 16.7M, 128 cores: full rebuild 108 s -> ~11 s; density loop
7 s; force loop 9 s. A `glassrelax` wrapper script (package root) runs the
CLI with --param/-o derived from the pre-IC filename.

**Round 2 (same day, target: SPH-EXA per-iteration parity, ~6 s at
23.9M)**: branch-based min-image replaced rint(dx/box) in all pair loops
(positions kept wrapped; relax wraps its input), per-thread workspaces +
in-place quickselect removed all per-particle allocations, and relax()
gained reorder=True: particle arrays are permuted into cell order at each
rebuild (cell_order + numba gathers), making every SPH gather cache-local;
outputs are mapped back to input order. The numpy glue (step+wrap+max
reduction, pP+density error) moved into fused numba kernels
(stepper.apply_step, stepper.pseudo_pressure) - the rint wrap alone was
1.4 s/iter. Search radii are capped at the min-image envelope (an
uncapped stale-h radius once degenerated to an O(N^2) scan). Measured at
16.7M: density 0.14 s, force 0.46 s, h-update 0.47 s, full rebuild ~4 s;
end-to-end 6.1 s/iter with a rebuild every iteration, 4.2 s/iter with
--margin 16 (CLI flag; wider Verlet skin -> rebuild ~every other
iteration), vs SPH-EXA's 5.3-7.3 s at 23.9M on the same node.

**--one-sided** (relax(one_sided=True)): SPH-EXA-style transpose-free
force walk - both kernel sides evaluated from i's own candidate list, no
transposed graph built. Pair forces are then not exactly antisymmetric;
with margin m the deviation needs h_j/h_i > ((nsmooth+m)/nsmooth)^(1/3)
across one kernel radius and sits at the kernel edge where dW -> 0.
Measured force-field difference vs the symmetric walk: ~6e-4 relative L2
on a random pre-IC, ~1e-4 on a relaxed glass; converged glass quality is
statistically identical (test_one_sided_force_glass_quality), and the
256-cubed denserr/q1 traces agree to ~0.4%. 2.7 s/iter at 16.7M with
--margin 16 (vs 4.2 symmetric) - faster per particle than SPH-EXA. Keep
the default (symmetric) for the gdicgen validation; --one-sided is the
production-speed mode. Sharp-discontinuity check (28-cubed, 30:1 tanh step
over ~1 spacing, mhdcollapse-like): both modes converge with equivalent
interface statistics (band err 0.26 sym / 0.21 one-sided, bulk 0.02) -
the interface error is kernel-smoothing bias common to both operators
(and to the C code), not asymmetry; one-sided is safe for discontinuous
targets too. Remaining follow-up ideas: partial rebuilds
(per-particle skin/displacement tracking), numba.cuda.

## Context

The glass-relaxation IC generator lives in a hacked Gasoline build inside the testsuite repo (github.com/robertwissing/testsuite) at `icgenerator/gasoline/` (`-DICGEN`), invoked as `gasoline.gdicgen` from `runtest.sh:1007`. It relaxes a particle pre-IC into a glass matching a target density profile via damped pseudo-pressure-gradient descent. Density profiles are hardcoded in C (`pkdUpdateDensity`), and experimenting with SPH gradient operators means recompiling.

Goal: a numba-accelerated Python reimplementation where:
- target density comes from a **density table** (1D/2D/3D, the existing `densitytable_xdr` format from `setupfiles/IC_denstable.py`) or any Python **callable**, plus a coordinate mapping (x, y, z, signed-x, spherical r, cylindrical r) — no analytic profiles built in; IC_setup*.py scripts generate tables instead
- the SPH force part has **pluggable flavours** (gradient operators / symmetrizations); initially the standard kernel gradient with GDFORCE `1/(ρᵢρⱼ)` and standard `1/ρ²` variants; ISPH matrix-corrected gradients addable later without restructuring
- delivered as an importable **library** (numpy in/out, callable from `setupfiles/IC_setup*.py`) plus a thin **CLI** (tipsy pre-IC + density table in → tipsy glass out) that can functionally replace `gasoline.gdicgen` in runtest.sh

**Working location**: `/mn/stornext/d17/extragalactic/personal/robertwi/Projects/ICGEN/sims/OPTKERNELSIMS/a/testsuite/`. All paths are repo-relative. Key anchors in this clone: ICdir/analysisdir at runtest.sh:9–11, ICCONFIG at runtest.sh:44, gdicgen invocation at runtest.sh:1007; `ICGenerator` at icgenerator/gasoline/smoothfcn.c:5794; density-based-h line (`fBall2 = 4·hdens²`) active at smoothfcn.c:6082.

Confirmed available: numba 0.60.0, scipy 1.13.1 (cKDTree with `boxsize` + `workers`).

## Reference algorithm (verified in C source)

Per iteration (effectively damped steepest descent):

1. **Neighbors/h**: h_i = ½ · distance enclosing nSmooth neighbors; recomputed every iteration; periodic minimum image.
   ⚠️ Copy difference: in this repo's copy, `DenDVDX` actively sets `fBall2 = 4·hdens²` with `hdens = (m·hfact3/ρ)^⅓`, `hfact3 = 3·nSmooth/(32π)` — density-based h used by the subsequent force resmooth; in the older u3 copy this line is commented (pure kNN h). v1 implements **kNN h** (simple, regular arrays); density-based h noted as follow-up option (needs radius queries/CSR lists). Final default pinned in step 5 verification against C output.
2. **SPH density**: ρ_i = Σ m_j W(r_ij, h_i), M4 cubic spline, gasoline convention (BALL2 = (2h)², fNorm = 1/(2π)·ih³).
3. **Target density** ρ⁰_i from table/callable at mapped coordinate (table interp = port of `pkdGetDensityFromTable`).
4. **Pseudo-pressure**: pP_i = (ρ_i/ρ⁰_i)^rhopow (= ICrhopow, default 2.0).
5. **Pair force** (`ICGenerator`): ICvel_p −= (ig_p·dx + ig_q·dx)·(dW/dr-term)·m_q, symmetric on q; ig from flavour (GDFORCE: pP/(ρ_pρ_q); standard: pP/ρ²). ISPH (future): dx → −(C·dx) with kernel value instead of derivative.
6. **Displacement normalization** (pkdICGetFmax/NormalizeFmax/NormalizeVel): fmaxpart_i = max |component|; ICfmax_i = log(fmaxpart_i·T+1)/log(fmax·T+1), T=0.1; speednorm = (1−1/1000)·ICfmax + 1/1000; ICvel_i ← unit(ICvel_i)·speednorm.
7. **Move**: x_i += ICvel_i · ICR0 · h_i; periodic wrap.
8. **Adaptive step control** (msrICGetFmax + restart logic in msrCalcEandL): Q1 = rms |ICvel|; step-zero rescale ICR0 /= h_avg, ICRLOOP0 *= ICR0; bad move (Q1 grew) → ICR0 /= (1+0.8·rate), good → ICR0 *= (1+0.16·rate); polish phase uses (1+0.1·rate)/(1+0.002·rate); nfac = 0.025·(1/nSmooth)^⅓; when ICR0 < nfac: re-anneal to ICRLOOP0 if converging, else halve ICRLOOP0 and restart; finish when ICRLOOP0 < 2·nfac after polish. runtest.sh:44 defaults: ICR0=1.0, ICRLOOP0=0.5, ICR0Rate=1.5, ICrhopow=2.0, n=4800 iters, oi=4800.
9. **Diagnostics** (optional, same pair loop): Q1_i (partition of unity → 1), E0_i (zero-order force error), Q4 tensor.

## Layout

```
icgenerator/glassgen_python/
  PLAN.md           # this plan document
  glassgen/
    __init__.py     # exports relax, RelaxResult, TableDensity, CallableDensity
    kernels.py      # @njit M4 cubic spline w(q), dw(q), gasoline normalization
    neighbors.py    # cKDTree kNN -> (idx, dist, h); Verlet-skin reuse; transpose graph
    density.py      # DensityField: TableDensity (XDR-compatible), CallableDensity, mappings
    flavours.py     # registry: gdforce, rho2 (ig/dterm/weight @njit pieces per flavour)
    sph_core.py     # @njit prange loops: density, force-loop factory, diagnostics
    stepper.py      # @njit displacement normalization + StepController (pure Python)
    relax.py        # relax() driver, RelaxResult, IterStats
    cli.py          # python -m glassgen.cli : tipsy in/out
  tests/            # pytest: test_kernels, test_density, test_neighbors, test_stepper, test_relax_e2e
```

Import path for consumers: `setupfiles/IC_setup*.py` scripts add the package root once (`sys.path.insert(0, .../icgenerator/glassgen_python)` relative to the script location) and `import glassgen`; the CLI runs as `PYTHONPATH=icgenerator/glassgen_python python -m glassgen.cli ...` from runtest.sh.

## Key design decisions

**Neighbor search with Verlet-skin reuse** (key to scaling — tree queries, not the pair loop, dominate at large N): scipy `cKDTree` (`boxsize` periodic, `workers=-1`), `tree.query(k=nsmooth+1+margin)` → dense `(N, k)` int32 array; h_i = ½·dist over the nsmooth nearest. The list is **reused across iterations**: track accumulated per-particle displacement since last build, rebuild only when `max(disp) > skin·h` (margin neighbors + distance re-sort keep h and the gather set correct between rebuilds). Early sweep rebuilds often; the long polish phase rarely — ~5–10× cut in neighbor cost. Target: 128³ (2.1M) full relaxation ~1–2 h on one fat node; single-node ceiling ~5–10M particles, beyond which a future `numba.cuda` force loop (factory design keeps it swappable) or the MPI C version takes over.

**Race-free parallel force loop**: gather-only reformulation. Particle i receives (a) terms from j ∈ idx[i] (h_i kernel) and (b) terms from the **transposed kNN graph** (j's whose ball contains i, h_j kernel), built as CSR in O(N·k) numpy. `prange` over i writes only `icvel[i]` — no atomics. Costs 2 kernel evals per unordered pair; buys perfect parallelism.

**Flavour dispatch**: a single `flavour` parameter names the full operator combination — `"gdforce"` (default), `"rho2"`, later `"isph-gdforce"` etc. Each registry entry in `flavours.py` bundles its `ig(...)`, `dterm(...)`, `weight(...)` @njit pieces; `sph_core.get_force_loop(flavour)` builds and caches the @njit(parallel=True) loop with those pieces inlined as closures. Adding ISPH later = one `dterm_isph` + kernel-value `weight` + registry entry + a C-matrix accumulation loop filling a `cmat` array already in the force-loop signature (dummy `(1,1)` when unused). No restructuring.

**Density field**: evaluated in plain Python/numpy once per particle per iteration (O(N), cheap) — callables have zero numba constraints. `TableDensity` does @njit lin/bi/trilinear interp, clamped at edges; `from_xdr()` reads the existing `densitytable_xdr` format with a **vendored struct-based reader** (`IC_denstable.py` uses `xdrlib`, removed in Python 3.13). Mappings `"x","y","z","r_sph","r_cyl","signed_x"` = C direction codes 1–6; mapping applies to 1D tables, 2D/3D tables index (x,y)/(x,y,z) directly as in C.

**Step controller**: faithful port of the C state machine (constants identical: 0.8/0.16, polish 0.1/0.002, T=0.1, maxspeed=1000, nfac), isolated as a pure-Python `StepController` class with state {icr0, icrloop0, q1_prev, phase ∈ SWEEP/POLISH} — unit-testable and swappable. No redesign: the C schedule is what's validated.

**API**:
```python
relax(pos, mass, density, box=None, nsmooth=64,
      flavour="gdforce", rhopow=2.0,
      icr0=1.0, icrloop0=0.5, icr0rate=1.5, max_iter=4800,
      diagnostics=False, callback=None, callback_every=100) -> RelaxResult
# RelaxResult: .pos, .h, .rho, .niter, .converged, .history
```
float64 internally, tipsy f32 at edges, `cache=True` on all @njit.

**CLI** (`python -m glassgen.cli`): `-o out -I pre.00000 --table densitytable_xdr --direction r_sph -s 64 -n 4800 -oi 4800 --icr0/--icrloop0/--rhopow/--icr0rate [--param pre.param | --box LX LY LZ | --open]`. Reuses `setupfiles/readtipsy.py`; box parsed from .param (à la `SPH_methods.read_box_size_from_param`); writes `${out}_IC` + optional `${out}.NNNNN` snapshots via callback — functional replacement for the gdicgen line at runtest.sh:1007 with ICCONFIG defaults (runtest.sh:44).

## Implementation order (each step verified before the next)

1. **kernels.py** — verify: ∫W dV = 1 (quadrature), dW vs finite diff, W(r≥2h)=0.
2. **neighbors.py** — verify: brute-force O(N²) match at N=500 (periodic+open); transpose-graph correctness; lattice h sanity; skin-reuse equivalence to fresh rebuild.
3. **density.py** — verify: XDR round-trip vs `IC_denstable.write_density_table_xdr` output; linear-table interp exact; callable vs table agreement.
4. **sph_core density loop** — verify: uniform lattice/glass ρ = n·m to ~0.1%; Q1→1; cross-check `SPH_methods.calculate_density` at small N.
5. **flavours + force loop** — verify: against a one-iteration O(N²) pure-numpy serial reference (written in tests) to fp tolerance, both flavours (gdforce, rho2); sign test (force points overdense→underdense). Pin the `pair_dwdr` h-convention (owner-h vs symmetrized) and the kNN- vs density-based-h question here, against C behavior.
6. **stepper.py** — verify: hand-traced Q1 sequence exercising step-zero rescale, shrink, grow, re-anneal, halve-restart, polish, finish.
7. **relax.py** — verify E2E (local, small N): (a) 32³ random pre-IC, periodic box, uniform ρ⁰ → σ(ρ)/ρ < ~1%, Q1→1; (b) 1D sigmoid/step table at 32³ → binned ρ(x) matches table; (c) 64³ timing to extrapolate 128³.
8. **cli.py** — verify: round-trip a real `datafiles/$prefile.00000` from a test case; non-position columns unchanged; then user runs a full-size runtest.sh comparison vs `gasoline.gdicgen` (statistical equivalence of density-error histograms / radial distributions — particle-for-particle identity not expected).

## Files referenced (repo-relative)
- C reference: `icgenerator/gasoline/` — `ICGenerator` (smoothfcn.c:5794), normalization (`pkdICGetFmax`/`pkdICNormalizeFmax`/`pkdICNormalizeVel` in pkd.c), table interp `pkdGetDensityFromTable` + `pkdUpdateDensity` (pkd.c), step control `msrICGetFmax` + restart logic (master.c)
- Reused Python: `setupfiles/IC_denstable.py` (table format; note duplicate `read_density_table_xdr` defs — second wins), `setupfiles/readtipsy.py` (tipsy I/O), `setupfiles/SPH_methods.py` (kernel/param-parsing reference), `setupfiles/SPH_evolve.py` (prior SPH force-loop reference)
- Pipeline target: `runtest.sh` lines 9–11 (ICdir/analysisdir), 44 (ICCONFIG), 1007 (gdicgen invocation)
- Consumer pattern: `setupfiles/IC_setup_evrard.py` (already builds tables via IC_denstable)

## Future tasks

1. **accdisk benchmark case + density-table setup** (hardest test). accdisk
   uses `dICdensprofile=7` (the disk profile, `dir=5` r_cyl) which
   `benchmark.target_from_param` does not port (it raises for profile≥4). The
   clean fix per the design (no analytic profiles built in): **add a density-
   table generation step to `setupfiles/IC_setup_accdisk.py`** (à la
   IC_setup_mhdcollapse / IC_setup_evrard, using `IC_denstable.generate_density_table`
   + `write_density_table_xdr`) that tabulates the profile-7 sigmadisk field
   (pkd.c:7443-7489), then point a `benchmark.py` `accdisk64` case at that
   table via `TableDensity.from_xdr`. Until then accdisk is the deferred case
   in the CASES dict.

2. **Ball-gather backend DONE (2026-06-13, `ball.py`)**. `gather='auto'` (CLI
   `--gather`) routes measured/target to a fixed-radius CSR **ball gather**
   (no kNN, no k-selection) — the CORRPARTITION gasoline path: h sets a *mean*
   nSmooth via the h-ρ relation `h=(m·hfact3/ρ)^⅓` (mean ~64 verified), and
   density+pP+force all use that same h (fully consistent, unlike the
   inconsistent dense kNN-h path). target h is a priori from ρ0; measured
   bootstraps from ρ0 then uses the previous ρ. Floor = expand the ball until
   ≥nsmoothmin members (retry). Verlet reuse (padded radius, per-particle
   slack R_cached−2h covering motion *and* h growth) + cell-order reorder:
   ~8–28% rebuild rate. Result: **measured-ball is faster than kNN** on the
   hard structured case (mhd12: 14.9 vs 15.5 ms/it, 24 vs 30 s); target wins
   on total time (far fewer iters) but is slower per-iter — NOT because of
   ρ0 (measured evaluates ρ0 per-iter for pP too), but because target h is set
   by the *sharp* ρ0: outer particles at the 360:1 interface get a large h
   whose 2h ball reaches into the dense core (mhd12 nnb tail p99 264 / max 387
   vs measured's 67/68), inflating every neighbour loop, and ρ0(pos) changing
   fast across the gradient doubles the rebuild rate (57% vs 28%). measured h
   tracks the smooth actual density → clean ~64 neighbours, fewer rebuilds.
   Per-phase timing now printed (RelaxResult.timings): **nbr dominates
   ~55–80%**, frc ~10–20%, dens ~10%.

   **Rebuild optimisations DONE (2026-06-13, ball.py)**: two of the flagged
   items landed, both exact (neighbour sets + transpose verified identical to
   the old path, convergence trajectory bit-for-bit unchanged).
   (a) **Parallel CSR transpose** (`_transpose_csr` now uses `_ball_esrc` +
   `_ball_transpose_par`): a per-edge source-row map `esrc` (built parallel
   over rows, each row owns a disjoint range) lets the counting-sort reuse the
   exact kNN `_transpose_fill_par` structure — the kNN path recovered i as
   `e // k` from its dense array; the ball graph is ragged so the row map is
   materialised. This was *the* dominant serial cost at scale (~50M edges in
   the mhd50 target-h build): transpose ~150 ms parallel vs the old serial
   fill. (b) **Fused single-pass build** (`_fused_gather` → `_ball_store` +
   `_ball_scatter`): phase-2 now walks the cell neighbourhood *once*, storing
   ids into per-thread buffers sized from the phase-1 unpadded counts ×
   (1+skin)³ headroom, then prefix-sums and scatters into a tight CSR — vs the
   old count-walk + fill-walk (two walks). An overflow flag (a thread exceeds
   its a-priori cap in a sharp region) falls back to the robust 2-pass build,
   so it is always correct. Net: nbr 18.2→17.2 s on the realistic mhd50
   measured run (~5%, the saved count-walk minus the scatter memcpy);
   transpose change is the larger win. **Full mhd50 311k convergence
   head-to-head (the run that motivated this — was 745 ball vs 469 kNN, ball
   59% slower): now measured-ball 538.8 s (3134 it, 171.9 ms/it, err_mean
   0.0108) vs kNN 512.2 s (3389 it, 151.2 ms/it, err_mean 0.0117).** So ball
   does *not* yet beat kNN on wall — but the gap closed from 59% to ~5%.
   Notably measured converges in *fewer* iterations (3134 vs 3389) and to
   *better* density error (0.0108 vs 0.0117, the density-h quality edge); kNN
   still wins on wall purely via per-iter speed (151 vs 172 ms/it), i.e. the
   ball rebuild. nbr is still ~85% at ~58% rebuild rate, so the per-iter gap
   is now entirely rebuild *frequency*. mhd12 measured-ball 9.9 ms/it (was
   14.9) and converges in fewer iters than kNN (1621 vs 1962).

   **Per-particle stale test DONE (2026-06-13, ball.py needs_rebuild)**: the
   rebuild trigger was global — `2*dmax > min_i slack_i` (dmax = max mover,
   slack_i = R_cached_i - 2*h_now_i) — so one fast interface mover coupled to
   an unrelated small-slack particle forced a full rebuild. Replaced with the
   per-particle Verlet condition `exists i: pdisp[i] + dmax > slack_i` (the
   per-particle displacement `pdisp[i]` was already accumulated by
   stepper.apply_step). Provably correct (a neighbour k within 2h_i now was
   within 2h_i + pdisp[i] + pdisp[k] <= 2h_i + pdisp[i] + dmax at build, so if
   the test says reuse then k is in the cached ball — locked in by
   test_ball.py::test_per_particle_stale_never_drops_a_neighbour, a brute-force
   superset invariant on every reuse step) and provably <= the old trigger rate
   (new-true => old-true). One wiring fix: needs_rebuild() was evaluated twice
   per iter (gate the reorder, then again inside update()); the reorder
   permutes pos/h_in but not the cached rcached/pdisp, so the second, tighter
   test could read misaligned arrays. The driver now evaluates it **once** on
   the pre-reorder state and passes the decision into `update(pos, h_in,
   rebuild=...)`. Effect on mhd50: rebuild rate 58% -> 48% on a short probe,
   but over a full run the polish-phase cut is much larger and per-iter drops
   a lot: **measured 171.9 -> 113.7 ms/it, now FASTER per-iter than kNN
   (152.0)** — per iteration the ball finally beats the kNN backend.

   **BUT end-to-end is a wash, and why is the real open item.** Full mhd50
   back-to-back (same node load; kNN reproduced at 3389 it / 515 s, confirming
   the earlier 469-vs-512 gap was pure load): per-particle **measured ran the
   full 4800-iter cap [NOT CONVERGED], 545.6 s**, vs the pre-per-particle
   measured that *converged* at 3134 it / 538.8 s. Quality at the 4800 cap is
   nonetheless **equivalent to converged kNN** (err_mean 0.0117 = kNN 0.0117,
   err_p95 0.0221 < kNN 0.0259, Q1 0.9947) — the glass is effectively done; the
   StepController's formal stop (icr0 < stop_frac*nfac in polish) just never
   trips on this trajectory within budget. mhd12 shows the same direction:
   niter 1621 -> 1772 (+9%). The per-particle test perturbs the convergence
   *trajectory* (density/force are provably identical per-iter — lists stay
   complete + kernel-masked — so this is FP-summation-order via the changed
   rebuild->reorder schedule, plus staler hfloor between rarer rebuilds feeding
   h_avg into the controller), and the formal stop trips later.

   **Root cause pinned (2026-06-13, diag instrumentation) — it is the stop
   THRESHOLD, not the controller.** The state machine works correctly: on the
   per-particle measured mhd50 run icrloop0 halves the full 12 times
   (40.07 -> 0.0098), re-anneals 21 times, and polish IS entered (~iter 1864).
   So the earlier "re-anneals forever / never halves" guess was WRONG. The
   benchmark just uses the strict `stop_frac=0.05`, i.e. stop only once polish
   AND icr0 < 0.05*nfac ~ 3.1e-4. Polish starts with icr0 ~ 0.0098, so the run
   must then grind icr0 down ~30x via the *gentle* polish rates (up 1.003 /
   down 0.87 per move) driven by the noisy q1 sequence — a biased random-walk
   threshold crossing. That crossing is what the per-particle trajectory
   stretches past the 4800 cap; **with stop_frac=2.0 (stop at polish entry) the
   exact same per-particle run converges at 1864 it.** Crucially the polish
   grind is past useful quality: denserr is already ~0.011 and FLAT/noisy
   through it (0.0108 @3134, 0.0117 @4800 — non-monotonic), so those thousands
   of deep-polish iters do not improve the glass. So the per-particle test did
   not "break" convergence; it exposed that stop_frac=0.05 is a noisy
   random-walk crossing well past quality convergence. Net: the large per-iter
   win (113.7 vs kNN 152 ms/it) is real but currently spent on quality-useless
   deep-polish iterations. **Fix = stop on a quality plateau (the glass is done
   at polish entry), not on icr0 < 0.05*nfac.**

   **Plateau stop DONE (2026-06-17, `stepper.PlateauStop`).** Implemented as a
   plateau on the **residual-force RMS `frms`** (the C ICvelAvg, the quantity
   each iteration actually descends) rather than denserr/Q1 as first proposed:
   frms is already computed every iteration at zero extra cost, IS the
   force-equilibrium objective, and reaches its noisy floor at the same polish
   entry where denserr goes flat (root-cause data above: denserr 0.0108@3134 ->
   0.0117@4800, non-monotonic — so frms-plateau == quality-plateau here, and
   frms avoids an extra diagnostic pass). `PlateauStop(tol=1e-3, window=50,
   patience=300)`: a trailing-median smooth over `window` kills per-iteration
   noise, then a best-so-far + `patience` non-improving-samples test fires.
   Wired into `relax()` as the PRIMARY finish (`force_stop=True` default,
   `force_tol/force_window/force_patience` knobs); the old icr0 < stop_frac*nfac
   crossing is kept as a backstop and max_iter as a hard floor.
   `RelaxResult.stop_reason` reports which fired ('force_plateau' / 'icr0' /
   'max_iter'). pstop is fed ONLY in the polish phase, so the window fills from
   polish entry. Flows through every path by default (cli.py/runtest, progressive
   coarse+fine, benchmark — benchmark adds `--no-force-stop`). Verified: 55/55
   tests, incl. `test_force_plateau_stop_banks_deep_polish_on_ball_path` (ball
   measured-h, gradient target, icr0 backstop disabled so only the plateau can
   stop early) — fires at iter 1403 vs the 2500 cap at byte-equal denserr
   (mean 0.0631, p95 0.2030 both). On uniform targets polish is short and the
   icr0 backstop still trips first (force_plateau is specifically for the long
   gradient-target polish grind). With this stop measured banks its per-iter
   edge over kNN. Remaining ball optimisations: skin tuning, numba.cuda.

3. **Setup + analysis framework end-to-end DONE (2026-06-14, benchmark.py)**.
   `benchmark.py` now drives the REAL test-suite pipeline instead of synth boxes
   / committed pre-ICs. The `synth`/`file` CASES are replaced by a single
   `kind='setup'` schema (`dict(test=, nx=, extra=[])`); `gen_pre_ic()` shells
   out to the runtest.sh contract `IC_createsetup.py <test> <nx> <vm> 1 <out>
   [extra]` into a `_pipeline/` work dir (datafiles/+initruns/), **reusing an
   existing `<prefile>.00000` rather than regenerating** (the runtest "Skipping
   creation" cache; `--force-gen` overrides). Box + target come from the
   generated `.param` via the existing `box_from_param`/`target_from_param` (no
   analytic profiles hardcoded). New flags: `--gen-nx` (resolution sweep),
   `--gen-vm`, `--work-dir`, `--force-gen`, `--analyze`, `--analyze-res/-backend`.
   Cases: mhdcollapse12/sedov64/mhdcollapse50/shocktube256/currentsheet64
   (accdisk still gated on its density table, Future task #1).

   **Analysis: in-process, glass-appropriate (NOT IC_analysis_<test>.py).** The
   per-test scripts assume an evolved t>0 run — sedov's analytic divides by t
   (`D=0.4*R/t` -> ZeroDivisionError on a t=0 glass) and the metric series are
   meaningless on a static IC. So `--analyze` writes the relaxed glass as a tipsy
   snapshot (`write_glass_snapshot`, reusing cli.py's round-trip + copying the
   pre-IC aux) and calls `analyze_glass()`, which uses the framework's render +
   profile BUILDING BLOCKS directly (`render_panels` grid SPH density map +
   `binned_profile` of SPH rho vs the TARGET rho0 along the dICdensdir
   coordinate + a density-error histogram) — robust for every test (verified
   sedov uniform + mhdcollapse tanh-sphere r_sph). Figures land in
   `_pipeline/analysis/<run>_{render,profile,denserr}.png`. `_pipeline/` is
   gitignored. pytest still 35/35.

4. **Density-floor gap follow-up**: even with measured-h, C stays ~2× tighter
   on density (currentsheet 0.0050 vs 0.0022). Candidate causes to check:
   C polish depth / iteration schedule, effective nSmooth, or a second
   density consistency pass. Quantify before deciding it matters.

5. **Multi-resolution / progressive relaxation (IMPLEMENTED 2026-06-14)**.
   Idea: relax a COARSE particle set first to capture the
   large-scale density field cheaply, then `split` particles up to full
   resolution and run a FINE, LOCAL relaxation (step ≤ h, no cross-box jumps)
   that only fixes small-scale structure. The coarse glass already encodes the
   large-scale field, so the fine phase likely needs no `ICRLOOP0` re-anneal
   loop (the loop exists precisely to remove large-scale error during a single
   full-res run — see Future task #2's deep-polish waste). Faster + cleaner than
   one monolithic relax, esp. for density-gradient ICs (mhdcollapse 360:1).

   **Decisions locked with user:**
   - **Split/merge** = new glassgen-native module `glassgen/resample.py` (pure
     arrays, no I/O), algorithm ported from `SPH_methods.split_particless`
     (SPH_methods.py:754-983, binary 1→2, `nn_ort` placement `parent ±
     ortho·dr·0.4`, adaptive-k separation floor `≥ 1/(0.5·nsmooth)^(1/3)·h`) and
     `merge_my_neighbor` (:1216-1435), but rewritten to round-trip ALL fields
     (mass/pos/vel/rho/u/soft/metals/B/spin/momenti/tform/h) via a declarative
     per-field rule table (HALVE: mass,u,B,momenti; INHERIT: vel,rho,h,soft,
     metals,spin,tform; OFFSET/centroid: pos). rho/h re-derived by next relax(),
     so split only needs sane seed pos+mass. Use glassgen's OWN
     `NeighborList.update` for the nn direction (a new `@njit _nn_vectors`
     scanning each kNN row for the two smallest min-image d²; idx rows aren't
     distance-sorted on the grid path) — do NOT import SPH_methods. Reach a
     target multiplier by repeating binary split log2(factor) times (assert
     power of 2), rebuilding neighbors each pass; mass→parent/2^p auto-matches
     full-res per-particle mass when factor=(nx_full/nx_coarse)^3 (assert).
   - **Standalone CLI** `glassgen/resample_cli.py` + `glassresample` wrapper
     (callable from main folder alongside runtest.sh): TIPSY in/out reusing
     cli.py round-trip; `-I -o --op{split,merge} --factor 8 -s 64 --dist 0.4
     (--param|--box|--open)`. RISK: N changes on split, so aux sibling files
     (B/spin/...) must be RE-EMITTED with new count header + correct per-label
     dtype (readtipsy `intauxfiles`/`writetipsyaux`/`writetipsyVecaux`) — cannot
     byte-copy like the relax path (benchmark.py:207-219).
   - **Progressive driver** `glassgen/progressive.py`:
     `relax_progressive(pos, mass, density, box, *, nx_full, nx_coarse=64,
     levels=None, nsmooth=64, coarse_seed=None, **relax_kw) -> RelaxResult`.
     Per-level dict `{name, split_factor, move_cap, reanneal, max_iter, icr0,
     icrloop0, stop_frac}`; default 2-level (coarse: split_factor=1,
     move_cap=False, reanneal=True, big icr0/cross-box; fine: split_factor=
     (nx_full/nx_coarse)^3→pow2, move_cap=True, reanneal=False, small icr0).
     **Coarse-skip: if nx_full ≤ nx_coarse → single full-res level** (so nx=32
     behaves like today; nx=128 & nx=256 share the SAME coarse run at nx_coarse,
     differing only in fine split_factor 8 vs 64). N-level generic over `levels`
     so a MEDIUM stage = insert a third dict. Coarse seed: core API accepts
     `coarse_seed` (keeps driver I/O-free + unit-testable), uniform-random
     fallback; pipeline supplies a density-consistent seed by shelling to
     IC_createsetup at nx_coarse (option A, à la benchmark.py:gen_pre_ic).
     Export from `glassgen/__init__.py`.
   - **Movement cap**: add trailing scalar `max_step_frac` to `apply_step`
     (stepper.py:67; default 0.0 = current behaviour). When >0: clamp each
     component to `max_step_frac·h[i]` and use a single rint wrap instead of the
     multi-box while-loop (step ≤ h ≪ L/2 ⇒ one image). `relax()` gains
     `confine=False` → passes `1.0 if confine else 0.0` at relax.py:243 (covers
     knn+ball). Keep new arg LAST (numba sig stability); assert/clamp if h→L/2.
   - **Optional re-anneal**: `relax(reanneal=True)` → `StepController`. In
     `update_restart` keep the `icrloop0<2·nfac→icfmax=-1` polish entry, then
     `if not self.reanneal: return` before the ladder (polish + stop_frac exit
     still reachable; just no icr0 re-inflation/halving). Default True =
     unchanged.
   - **Pipeline exposure** at the relaxation-DRIVER layer, NOT a new entry code
     (don't touch IC_createsetup `0/1/glass` contract): cli.py gets
     `--progressive --nx-coarse 64 --coarse-iter --[no-]fine-confine
     --[no-]fine-reanneal`; runtest.sh adds a `GLASSMODE` string at the relax
     invocation (:1007), empty = single-stage. The second IC_createsetup
     (:1031) consuming `${glassfile}_IC` is unchanged.

   **Order:** resample.py + (apply_step cap & reanneal flag) in parallel →
   progressive.py → resample_cli.py → cli/runtest wiring → benchmark variant +
   this PLAN.md status. **Verify:** tests/test_resample.py (mass conservation,
   count=N·factor, child locality incl. box-face wrap, field rules,
   merge(split)≈conserved), tests/test_stepper.py (max_step_frac=1.0 clamps ≤h &
   never multi-box-wraps; reanneal=False still reaches polish/stop_frac),
   tests/test_progressive.py (small uniform converges; coarse-skip == single
   relax). Benchmark `--progressive` on sedov64 (uniform) + mhdcollapse50
   (gradient, the key case): progressive quality (mean/p95 denserr, Q1) within
   ~1.1× of single full-res and faster total wall (coarse+split+fine). Keep
   `python -m pytest tests/` green. **Full design also at**
   `~/.claude/plans/idea-for-icgenerator-evolve-keen-wall.md`.

   **Built (all shipped, 50/50 tests pass):** `glassgen/resample.py` (native
   split/merge, per-field rules, glassgen NeighborList for nn_ort directions,
   factor=2^p; `halve_fields` override), `glassgen/progressive.py`
   (`relax_progressive` + `default_levels`, fixed nx_coarse, coarse-skip when
   nx_full<=nx_coarse), `glassgen/resample_cli.py` + `glassresample` wrapper
   (tipsy in/out, aux files re-emitted with new N — verified on the real 157k
   shocktube64 pre-IC: x8 -> 1.26M, mass + B-flux conserved, all 4 aux files).
   stepper.apply_step gained trailing `max_step_frac` (clamp each component to
   frac*h + single-image wrap; relax `confine=`), StepController gained
   `reanneal` (relax `reanneal=`); confine passes box=None to the controller to
   drop the periodic box-amplification so capped steps start at ~icr0*h.
   cli.py `--progressive/--nx-coarse/--nx-full/--no-fine-confine/
   --fine-reanneal`; runtest.sh opt-in `GLASSMODE` (default "" = C gdicgen
   unchanged) + `ICCONFIGPY`; benchmark.py `--progressive/--nx-coarse` variant
   (coarse seed via IC_createsetup at nx_coarse, option A). Tests:
   test_resample.py, test_progressive.py, extended test_stepper.py.

   **First findings (smoke, nx_full=16):** the fine confined+no-reanneal phase
   is CORRECT — on the uniform sedov case progressive matches single-relax
   quality exactly (err_mean 0.0066 vs 0.0065; Q1/E0/E1 identical), confirming
   the cap + box=None controller + polish-on-icr0<nfac path produce an
   equivalent glass. On uniform it is slightly SLOWER (6.4 vs 3.7 s) — no
   large-scale structure to save, just coarse-stage overhead, as expected.
   On the steep mhdcollapse 360:1 gradient at this tiny resolution progressive
   was faster (4.9 vs 11.1 s) but WORSE quality (err_mean 0.124 vs 0.052):
   nx_coarse=8 (~1k particles) cannot capture a 360:1 sphere, and the fine
   phase is local-only by design (moves <= h) so it can't repair large-scale
   error the coarse missed. Design premise made concrete: **the win needs (a)
   a real resolution gap (nx_full >> nx_coarse) and (b) nx_coarse adequate for
   the density dynamic range.**

   **Production-scale sweep (mhdcollapse50, 360:1 gradient, knn, nsmooth=64),
   nx_coarse {12,25,32} -> nx=50, vs single-relax baseline. denserr/force as
   mean/p95/max:**

   | config        |    N    | niter | wall  | denserr mn/p95/mx    | force mn/p95/mx     |
   |---------------|---------|-------|-------|----------------------|---------------------|
   | baseline      | 312,210 | 2765  | 478.8 | 0.0114/0.0252/0.248  | 0.47/1.68/19.8      |
   | c=25 (x8)     | 310,944 | 1501  | 177.6 | 0.0111/0.0239/0.228  | 0.63/1.86/15.5      |
   | c=32 (x4)     | 325,544 | 2393  | 290.2 | 0.0103/0.0229/0.226  | 0.49/1.40/16.5      |
   | c=12 (x64)    | 270,080 | 2142  | 174.4 | 0.0267/0.0572/0.263  | 0.81/2.57/28.5      |

   **c=25 (x8, N matched to baseline) is the sweet spot: 2.6x faster
   (177.6 vs 478.8 s) at equal-or-better quality** (denserr lower across
   mean/p95/max, incl. a *shorter* worst-case interface tail 0.228 vs 0.248;
   radial rho(r) and a slice render are visually identical to baseline -
   `_pipeline/analysis/mhd50_compare_{density,render}.png`). c=32 (x4) is 1.65x
   faster, similar quality. **c=12 (x64) confirms the too-coarse failure
   mode**: a 12^3 (~4.3k) seed can't resolve 360:1, so the local fine phase
   can't repair it -> denserr AND residual force both clearly worse (the new
   `force` diagnostic tracks quality: 0.81 vs ~0.5 for adequate coarse).
   nx_coarse must scale with the density dynamic range; 25-32 was enough for
   360:1, 12 was not.

   **Coarse-vs-fine breakdown (c=25, from RelaxResult.levels):** coarse stage
   38,868 particles -> denserr_mean 0.0295 (p95 0.086) in 1931 it / 60 s; fine
   stage 310,944 -> 0.0111 (p95 0.024) in 1501 it / 120 s. So the cheap coarse
   stage captures the large-scale field to ~3% in 1/3 of the total time, and
   the fine stage halves the bulk error while sharpening the interface. (Note:
   the per-particle force MAGNITUDE is not comparable across stages - it scales
   with the differing per-particle mass/h - so compare denserr across stages,
   force within a stage/resolution.)

   **Diagnostics added (this session):** `diagnostics.residual_force` (actual
   descent force |icvel| with real pP, the convergence/stop basis, distinct
   from E0's uniform-pP operator error); benchmark/cli/RelaxResult now report
   denserr + force as mean/p95/max with an N column; `relax_progressive`
   records per-stage end-state error in `RelaxResult.levels`
   (name/n/niter/converged/wall + denserr & force mn/p95/mx).
   `_scratch/compare_density.py` builds the baseline-vs-progressive density
   comparison figure (radial profile + denserr distribution overlay +
   side-by-side render).

   **nx=100 (~2.5M) result (mhdcollapse, 360:1, knn), progressive nx_coarse=25
   (x64) vs baseline:** progressive **2.6x faster AND better quality** - wall
   2092.8 s vs 5468.4 s, denserr mean/p95 0.0080/0.0186 vs 0.0091/0.0217, Q1
   4.56e-3 vs 4.64e-3, E1 2.90e-2 vs 3.08e-2, max ~equal (0.080). Stages:
   coarse 38,868 -> denserr 0.0289 in 1975 it / 66 s; fine 2.49M -> 0.0080 in
   4767 it / 2019 s. So the proven-adequate nx=25 seed (39k) captures 360:1 to
   ~3% in ~1 min even before a x64 jump, and the fine stage refines to better
   than baseline. **Key mechanism uncovered: the fine stage runs at 423 ms/it
   vs the baseline's 1360 ms/it - ~3x CHEAPER PER ITERATION at the same 2.5M
   size.** The confinement (capped sub-h moves) keeps the Verlet neighbour
   lists valid far longer, so neighbour rebuilds (the dominant cost) fire much
   less often. The fine even ran MORE iters than baseline (4767 vs 4021) yet
   finished in 2.6x less wall - at scale the win is mostly *cheaper iterations*,
   not just fewer. (force_mean higher for progressive, 0.63 vs 0.44 - the
   shallower force-polish seen at nx=50 - but denserr is better; force_max lower
   30.5 vs 38.3.) The ~2.6x ratio held from nx=50 to nx=100.

   Remaining TODO: sweep nx_coarse vs gradient steepness; try nx_coarse=50 (x8)
   at nx=100 to compare x8-vs-x64 split depth; 128/256; decide fine reanneal
   on/off (off was fine through nx=100). A confirming partial-rebuild / skin
   tuning of the *baseline* could shrink its rebuild cost too, but the confined
   fine phase gets that win for free.

6. **Cascade prolongation + stop-criterion deep-dive (IN PROGRESS, 2026-06-18).**
   Two strands worked this session; cascade is shipped+tested, the stop tuning
   has a clear next step (resume here tomorrow).

   **(a) Multigrid-style cascade split — DONE, 58/58 tests.**
   `progressive.cascade_levels(nx_full, nx_coarse, inter_iter=80)` +
   benchmark `--cascade [--inter-iter N]` (implies `--progressive`). Instead of
   one big x(nx_full/nx_coarse)^3 jump, climb to full res one binary (x2) split
   at a time, running `inter_iter` confined iters after each split, then a final
   full confined relax. Schedule: coarse(sf=1,full) -> k-1 x (sf=2, inter_iter
   confined smooth) -> fine(sf=2, full), k=log2(total factor). Exported from
   `glassgen/__init__`. Tests: test_cascade_levels_one_split_per_stage,
   test_cascade_progressive_converges.
   **Post-split diagnostic — DONE.** `relax_progressive` now records a
   `'<name>.split'` entry in `RelaxResult.levels` (via `_measure`: a zero-step
   `relax(max_iter=1, icr0=0)` to get rho/h on the static split set with the
   fine stage's own conventions, then `_level_diag`). Shows the quality the fine
   stage INHERITS from the split. test_progressive_records_post_split_diagnostic.
   **Trajectory dump — DONE.** benchmark `--history-every N` prints the fine
   stage's per-iteration (it, denserr, frms/force, icr0, polish) from
   `RelaxResult.history` (which relax already records) — lets us read where
   metrics actually plateau vs where the stop fires.

   **Cascade findings (mhdcollapse50 360:1, knn, nx_coarse=25; sedov64 uniform):**
   - The split damage IS much gentler one-x2-at-a-time: each post-split denserr
     ~0.04 / force ~10 vs the single x8 jump's **0.21 / 53-569** (~5x gentler).
   - But total iterations barely move because the fine-stage niter is dominated
     by the NOISY icr0 stop, not the starting quality (see (b)). mhd50: jump
     fine 1501it/119s (denserr 0.0111); cascade-80 fine 1308it/100s (0.0113);
     cascade-10 fine **3122it**/200s (0.0108) — identical fine.split start AND
     end quality as cascade-80 (0.0379 vs 0.0378; 0.0108 vs 0.0113) yet 2.4x the
     iters. sedov: cascade-10 fastest (fine 2180, 90.7s) vs jump 2699/109.6s.
     So inter_iter is case-dependent and the niter comparison is confounded by
     the stop noise.

   **(b) Stop criterion is too conservative — root-caused, fix half-built.**
   Full-budget trajectory (mhd50 jump, `--no-stop --no-force-stop`, 2500 it):
   - **frms does NOT plateau** — declines smoothly+monotonically to 2500
     (9.1e-3 -> 1.33e-3, still dropping ~25% over the last 1200 it). So the
     `PlateauStop` (patience=300) never fires before the icr0 backstop: on a
     monotone signal with tol=1e-3 the smoothed best is beaten every ~7 it so
     `wait` never accumulates. **The plateau is the wrong detector shape here.**
   - **denserr saturates early** (~it 1250 within 5% of the 0.0108 asymptote;
     p95 0.0237 even earlier); the deep polish past that is quality-useless.
   - **icr0 is noisy, everywhere in polish** (confirmed by probe): each polish
     iter icr0 takes one of two MULTIPLICATIVE moves — grow x1.003 (frms fell)
     or shrink x0.870 (frms rose), by the SIGN of that iter's frms change. frms
     is noisy -> sawtooth the WHOLE polish, constant *relative* amplitude
     independent of icr0 level (|dlog icr0| ~0.009-0.011 in both the high- and
     low-icr0 terciles). Self-regulating (95.8% grow/4.2% shrink at quasi-equil,
     slight net-down drift that grows as icr0 shrinks). So the icr0<stop_frac*nfac
     stop time is a noisy random-walk threshold crossing (1308 vs 3122 for the
     SAME quality) — NOT fixable by raising stop_frac (jitter isn't a low-icr0
     artifact).
   - **Key physics (user):** frms is the COMPLETE convergence measure — via
     pP=(rho/rho0)^rhopow it already embeds denserr (a density mismatch IS the
     residual pressure-gradient force) AND carries the SPH operator errors
     (gradient/partition-of-unity) that bare denserr can't see. denserr is also
     interface-biased (the 360:1 jump dominates its tail). So stop on FORCE, not
     denserr.

   **`force_target` (ForceReductionStop) — BUILT, needs ref-anchor fix.**
   `stepper.ForceReductionStop(target, window=50)`: capture a median-smoothed
   reference `ref` a window into polish, finish once smoothed frms < target*ref.
   Larger target = earlier = less conservative. Wired into `relax(force_target=
   0.0)` (0 = keep PlateauStop; >0 selects it; stop_reason='force_target') +
   benchmark `--force-target F`. Unit-tested (3 tests). Reproducible (smooth
   signal) where icr0 isn't.
   **SWEEP RESULT (mhd50 jump) — all OVER-FREEZE (the very thing to avoid):**
   f=0.30 -> fine 302it, denserr 0.0131, force **1.59**; f=0.25 -> 388it,
   0.0128, force 1.36; f=0.20 -> 532it, 0.0124, force **1.11**. Asymptote is
   denserr 0.0108 / force ~0.5, so ALL three stop ~3-8x too early with force
   2-3x the converged value (monotone trend: smaller f -> more iters/better but
   even 0.20 is still clearly under-converged). **ROOT CAUSE: `ref` is anchored too early.** The fine stage
   (reanneal=False) enters polish at **iter 23**, so `ref` is the median of
   iters ~23-73 — captured during the STEEP initial frms drop, when the glass is
   barely relaxed. Reducing from that near-unrelaxed baseline triggers while
   still under-converged.

   **NEXT (tomorrow):** fix the ref anchor so force_target measures reduction
   from a STABLE baseline, not the polish-entry transient. Options to try:
   (i) anchor ref to the PEAK/early-max smoothed frms rather than the first
   window; (ii) delay ref capture N iters after polish entry (or after frms's
   initial fast drop flattens); (iii) anchor to an absolute frms scale
   (e.g. a target relative to E0/operator floor). Then re-sweep f on mhd50 jump
   + cascade to find the least-conservative f that holds denserr/force within
   ~4% of asymptote without freezing (force should land ~0.5-0.6, denserr
   ~0.011). Also consider the COMPLEMENTARY rate lever (user): the 46:1 polish
   asymmetry (0.1 shrink vs 0.002 grow, stepper.py:181-182) drives the freeze;
   gentler shrink/stronger grow keeps icr0 mobile -> frms drops faster -> good
   quality in fewer iters. icr0rate (default 1.5) scales both but keeps the
   ratio; to change freeze, expose the two constants separately. Data shows we
   are NOT currently freezing prematurely (frms live at 2500), so the stop is
   the primary lever, rates secondary. Keep `python -m pytest tests/` green
   (58 now). Uncommitted: glassgen_python is still untracked in this clone.

   **(c) Active-fraction stop — BUILT + CALIBRATED (2026-06-18/19), 64 tests.**
   `stepper.ActiveFractionStop(eps, window)`: stop when the smoothed fraction of
   still-MOVING particles falls below `eps`. A particle is "active" while its
   NORMALISED step `|dx_i|/h_i = speednorm_i * icr0` exceeds `active_move_frac`
   (the log-normalisation already soft-freezes the quiet bulk ~1000x). This is
   the global/resolution-invariant form of the per-particle freeze: it counts
   who is still doing non-negligible work, not an absolute force level, and a
   median smooth over `window` kills the icr0 sawtooth so the crossing is
   REPRODUCIBLE (unlike the raw icr0 backstop's noisy 1308-vs-3122 crossing).
   Wired into `relax(active_stop=False, active_move_frac=3e-4, active_eps=0.08)`
   (stop_reason='active_set'; requires force_stop=True) + benchmark
   `--active-stop/--active-move-frac/--active-eps`. `relax.ACTIVE_THRESHOLDS =
   (1e-2,3e-3,1e-3,3e-4,1e-4)` are recorded each polish iter into
   `RelaxResult.active_history` for calibration (benchmark `--history-every`
   prints the table).
   **Calibration (mhd50 jump fine stage, 2000-it no-stop run):** denserr is a
   slow monotone grind to asymptote ~0.0108 (NO knee; frms never plateaus —
   confirms PlateauStop is the wrong shape and explains the force_target
   ref-anchor pain). The active-fraction THRESHOLD choice matters: `1e-2/3e-3`
   saturate to ~0 by it ~250 (stop way too early); the old default `1e-3`
   crosses eps=0.01 at it ~550 -> denserr +13% (the under-converged zone the
   force_target sweep was rejected for); **`3e-4` stays resolvable (0.08-0.11)
   through the quality plateau (it 1200-1600)** -> real control. So defaults
   CHANGED to move_frac=3e-4 / eps=0.08 (relax + benchmark).
   **Verified:** `--active-move-frac 3e-4 --active-eps 0.08` fires at fine-stage
   **it 1208**, stop_reason=active_set, denserr 0.0114 (+5.6%, p95 0.0241 ~=
   asymptote, max 0.226 < asymptote), force 0.70, Q1 0.0059, fine wall 118 s vs
   235 s for the 2000-it budget (2x). Strictly better than force_target (denserr
   0.0114<0.0124, force 0.70<1.1) AND reproducible where the icr0 backstop
   wasn't. **This is the recommended finish for the long gradient-target polish.**

   **(d) icr0-RATE asymmetry is INERT in the confined fine region — settled
   (2026-06-19), the user's "rate lever" hypothesis is a DEAD END here.**
   Exposed the two polish coefficients separately: `StepController(...,
   polish_down=0.1, polish_up=0.002)` -> `relax(polish_down=, polish_up=)` ->
   forwarded by relax_progressive; defaults reproduce the C ×0.870 shrink /
   ×1.003 grow (unit-verified: `(0.8,0.16)` correctly applies ×0.4545/×1.240).
   `_scratch/sweep_polish_rates.py` (modes fine/single, records icr0+frms+denserr
   p95/p99 per polish iter from a SINGLE shared coarse+split start; both stops
   disabled, fixed budget). **7-pair sweep (down 0.02-0.1, up 0.002-0.02) on
   mhd50 fine: frms(it) is RATE-INVARIANT** (≤4% spread past it 250; denserr p99
   identical to 3 digits; plot `_pipeline/analysis/polish_rate_sweep_fine.png`).
   icr0(it) is ALSO nearly invariant — all pairs self-regulate onto the same
   curve, NOT fanning out as first predicted.
   **Mechanism = negative-feedback loop (the real answer).** The grow/shrink
   decision is `sign(frms_prev - frms)`; icr0 above the critical step ->
   overshoot -> frms rises -> shrink; below -> progress -> frms falls -> grow.
   So icr0 self-regulates to a SETPOINT = the critical step size (~2.8-3e-3·h
   at equilibrium, a slowly-DECLINING target as the glass stiffens, starting
   ~1e-2·h at polish entry). The setpoint is set by the force-landscape
   stiffness (density-gradient curvature / kernel / nsmooth), NOT the rate
   constants. The rate only sets HOW VIOLENTLY icr0 jitters around the setpoint
   and the shrink-FRACTION needed for zero net drift.
   **Unrefined-rate confirmation (0.8/0.16 = the sweep-phase rates, applied in
   polish):** icr0 per-iter |Δln| jumps 0.006 -> 0.338 (**55x bigger jitter**),
   range widens [2.3,4.2]e-3 -> [1.25,7.2]e-3 (~3x), shrink-fraction self-adjusts
   2.2% -> 21.4% (drift stays ~0: 0.786·ln1.24 + 0.214·ln0.4545 ≈ 0), yet the
   SETPOINT (median 2.8-2.9e-3·h) and frms(it) are UNCHANGED (frms even ~5%
   lower, not worse). Plot `_pipeline/analysis/polish_rate_unrefined.png`.
   **Conclusion: the polish rate asymmetry changes neither convergence (frms)
   nor termination (icr0 stays parked far above the stop threshold and does NOT
   drift down in the reanneal=False fine stage) — it is NOT a quality/speed
   lever.** Termination must read the shared frms/quality curve directly = the
   active-fraction stop (c). Keep polish_down/polish_up at C defaults; knobs
   exposed if ever needed. (Single-relax/reanneal=True regime — where icr0 spans
   a wide range via box-amp + the re-anneal ladder — is the only place the rate
   could still matter; NOT yet tested, needs a full multi-k-iter convergence run
   per pair. `MODE=single` in the sweep script is ready for it.)

   **State:** 64/64 tests green. relax + benchmark defaults updated
   (active_move_frac 3e-4, active_eps 0.08). New sweep harness + plots under
   `_scratch/`/`_pipeline/analysis/` (both gitignored). glassgen_python still
   untracked in this clone. **NEXT options:** (1) make active_stop the pipeline
   default finish (cli.py/runtest GLASSMODE) now that it's calibrated; (2) test
   the rate lever in the single-relax regime (MODE=single, ~30 min); (3) the
   original force_target ref-anchor fix is now MOOT — active-fraction supersedes
   it; (4) commit glassgen_python.

7. **Momentum / Nesterov accelerated move step — THE next experiment
   (2026-06-19, scoped not built).** Task #6(d) proved the icr0 rate is inert:
   the current move is *first-order damped steepest descent* (`x_i += icvel_i ·
   icr0 · h_i`, no memory) and we are already at its speed limit — the setpoint
   is the max stable step, and SD's convergence is intrinsically capped by the
   problem condition number κ (error contracts by at most ~(κ−1)/(κ+1) per step;
   the 360:1 gradient makes κ large → the slow frms tail). **Tuning cannot beat
   SD's intrinsic rate; a better METHOD can.** The cheapest upgrade is momentum:
   give the descent memory so the cross-valley zig-zag cancels and the
   along-valley push accumulates (overdamped marble-in-honey → ball-with-inertia
   coasting the valley floor). Improves the ceiling ~κ → ~√κ, attacking exactly
   the slow tail. Complementary to progressive/multigrid (which already removes
   the large-scale/low-frequency slow modes — the 2.6× win; momentum attacks the
   remaining LOCAL slow modes in the fine stage).

   **Build:** add a trailing `momentum=0.0` knob to `stepper.apply_step`
   (default 0 = today's behaviour, numba sig stability — keep arg LAST). Keep a
   per-particle `disp_prev` workspace (input-order, permuted alongside pos on
   reorder). Heavy-ball: `disp_i = β·disp_prev_i + icvel_i ; x_i += disp_i·icr0·
   h_i ; disp_prev_i = disp_i`. Note `icvel` is ALREADY the normalised descent
   direction (post normalize_displacements), so β blends normalised directions —
   fine for heavy-ball. Nesterov look-ahead (`grad` at `x+β·v`) needs the force
   recomputed at the shifted positions = a second SPH pass per iter (2× force
   cost); START with heavy-ball (free) and only try Nesterov if heavy-ball
   overshoots. `relax(momentum=0.0)` passes it through (covers knn+ball; forward
   via relax_progressive **relax_kw). Interaction to watch: momentum + the
   confine cap (max_step_frac) — clamp the FINAL disp, not icvel; and momentum +
   the icr0 feedback controller (frms-sign decisions now see momentum-smoothed
   moves — may shift the setpoint, that's OK as long as frms(it) drops faster).
   Reset/decay disp_prev on a neighbour rebuild? Probably not needed (positions
   continuous); test both.

   **Verify/experiment (reuse the EXACT harness from #6d):**
   `_scratch/sweep_polish_rates.py` already records frms(it)+denserr p95/p99 from
   a shared coarse+split start with stops disabled — clone it to sweep `β ∈
   {0, 0.5, 0.7, 0.9, 0.95}` on the mhd50 fine stage and overlay frms(it) vs the
   β=0 baseline (the `polish_rate_sweep_fine.png` curve). WIN = β>0 frms(it)
   visibly BELOW baseline at matched iters (fewer iters to the same frms), with
   denserr p95/p99 no worse and no divergence (β too high → overshoot, frms
   oscillates/blows up — that's the β ceiling). Then check it composes with the
   active-fraction stop (should fire earlier on the faster curve) and re-run the
   mhd50/sedov progressive benchmark for end-to-end wall. Add a
   tests/test_stepper.py case (momentum=0 byte-identical to current;
   momentum>0 conserves nothing but must stay finite + reduce frms on a simple
   quadratic bowl). Keep `python -m pytest tests/` green.

   **Rationale recap for cold-start:** rate/stop levers are EXHAUSTED (#6). frms
   is still in genuine slow descent at the 2500 cap (not at the E0 operator
   floor — ~25% drop over the last 1200 iters), so there IS headroom an
   accelerated optimiser can capture. Momentum is the highest-payoff, lowest-cost
   remaining lever.

   **HEAVY-BALL BUILT + SWEPT (2026-06-19), 66/66 tests, the hypothesis HOLDS.**
   `stepper.apply_step` gained trailing `disp_prev, momentum` (numba sig stable;
   `momentum=0` takes the old branch, byte-identical — `test_apply_step_momentum_
   zero_byte_identical`). Heavy-ball `d=β·d_prev+icvel ; x+=d·icr0·h ; d_prev=d`
   on the ALREADY-normalised icvel (so β blends normalised descent directions);
   d_prev stored UNCLAMPED (bounded 1/(1−β) since |icvel|≤1 — `test_apply_step_
   momentum_accumulates_and_stays_finite` checks the 1/(1−β) fixed point), the
   confine cap clamps only the applied move. `relax(momentum=0.0)` allocates the
   per-particle d_prev workspace, **permutes it alongside pos on every reorder**
   (both knn + ball paths; NOT reset on rebuild — positions continuous), forwards
   through relax_progressive via **relax_kw. Sweep harness `_scratch/sweep_
   momentum.py` (clone of sweep_polish_rates.py: ONE shared coarse+split start,
   confined fine stage reanneal=False, both stops disabled, fixed budget) +
   overlay plot `_pipeline/analysis/momentum_sweep_fine.png`.

   **mhd50 fine sweep (β∈{0,.5,.7,.9,.95}, 1500-it budget, 310,944 particles)
   — frms(it), the descent objective:**

   | β    | frms@100 | frms@250 | frms@500 | frms@750 | frms@1500 | denserr p99 |
   |------|----------|----------|----------|----------|-----------|-------------|
   | 0.00 | 9.26e-3  | 4.89e-3  | 3.03e-3  | 2.42e-3  | 1.63e-3   | 0.127       |
   | 0.50 | 6.22e-3  | 2.85e-3  | 1.88e-3  | 1.52e-3  | **1.26e-3** (−23%) | 0.130 |
   | 0.70 | 5.92e-3  | 2.75e-3  | 2.06e-3  | 1.71e-3  | 1.46e-3 (−10%) | 0.128  |
   | 0.90 | 1.33e-2  | 5.40e-3  | 3.45e-3  | 2.80e-3  | 2.11e-3 (+29%) | 0.129  |
   | 0.95 | 2.22e-2  | 1.01e-2  | 5.79e-3  | 4.23e-3  | 2.89e-3 (+77%) | 0.129  |

   **Verdict: momentum works, β=0.5 is the robust pick.** It is the textbook
   heavy-ball picture (~κ → ~√κ): β=0.5 sits strictly BELOW the β=0 baseline at
   EVERY checkpoint (no crossover) — it reaches baseline's it=750 frms (2.42e-3)
   by ~it 300, i.e. **~2–2.5× fewer iters to a matched frms**, and ends 23% lower
   at the cap. β=0.7 is marginally better very early but its icr0 random-walks to
   ~1e-12 by the end (the frms-sign controller over-shrinks once the momentum-
   smoothed move stops changing sign — the controller-interaction flagged in the
   scope; harmless to quality but kills late progress), so 0.5 wins on the tail.
   **β≥0.9 is past the ceiling**: a visible frms overshoot at polish entry (peak
   ~0.4–0.5 in the plot) that never recovers — exactly the predicted divergence
   onset (still finite, just slower). denserr p99 is FLAT at 0.127–0.130 for all
   β (the interface tail is kernel-bias-limited, untouched by the optimiser) — so
   **the acceleration is free of quality cost.** The win is the descent RATE, not
   a lower floor (baseline reaches the same floor given enough iters); it is
   exactly complementary to progressive/multigrid (which removes the low-freq
   modes) — momentum attacks the remaining local slow modes.

   **β=0.3 added (the gap between 0 and 0.5) — the conservative sweet spot.**
   Same shared mhd50 fine start (β=0 reproduced 1.632e-3 bit-for-bit, confirming
   determinism). frms@{100,250,500,750,1500}:

   | β    | @100    | @250    | @500    | @750    | @1500   | icr0_final | denserr p99 |
   |------|---------|---------|---------|---------|---------|------------|-------------|
   | 0.00 | 9.26e-3 | 4.89e-3 | 3.03e-3 | 2.42e-3 | 1.63e-3 | 4.1e-4     | 0.1270      |
   | 0.30 | 7.53e-3 | 3.46e-3 | 2.22e-3 | 1.78e-3 | 1.40e-3 (−14%) | 7.1e-4 | 0.1282 |
   | 0.50 | 6.22e-3 | 2.85e-3 | 1.88e-3 | 1.52e-3 | 1.26e-3 (−23%) | 4.6e-4 | 0.1295 |
   | 0.70 | 5.92e-3 | 2.75e-3 | 2.06e-3 | 1.71e-3 | 1.46e-3 (−10%) | 1.0e-12| 0.1280 |

   β=0.3 sits cleanly BETWEEN baseline and 0.5 at EVERY checkpoint, monotone, NO
   overshoot, and has the BEST quality/stability profile of the three: denserr
   p99 tracks baseline (0.1282 vs 0.1270; it even dips LOWEST, ~0.1265, around it
   300-500 in the plot) where β=0.5 drifts up to 0.1295 in the tail, and its icr0
   stays the most MOBILE (7.1e-4 final, highest of all — no collapse, vs β=0.7's
   1e-12 freeze). So the spectrum is: **β=0.3 = conservative (safe, least quality
   drift, icr0 stays alive, still reaches baseline's it=750 frms ~2x faster);
   β=0.5 = most aggressive frms reduction (−23%) at a hair more denserr drift;
   β≥0.7 risks the icr0 collapse / β≥0.9 overshoots.** Recommend β=0.3 as the
   default and 0.5 as the speed setting. Plot `momentum_sweep_fine.png` now shows
   the 0/0.3/0.5/0.7 useful range.

   **Per-level wiring DONE (67/67 tests).** momentum is applied to the POST-SPLIT
   confined (reanneal=False) stages only — the coarse box-wide reanneal=True stage
   stays β=0 (its wide-icr0 + re-anneal-ladder regime is untested for momentum and
   it already moves fast via box amplification; user call). `default_levels(...,
   momentum=)` injects into the `fine` defaults; `cascade_levels(..., momentum=)`
   into the `fine` AND every `cascadeN` intermediate smooth (so each x2 split's
   local damage heals faster and the final stage inherits a lower-frms glass);
   `relax_progressive(..., momentum=)` threads it to default_levels when building
   its own schedule (kept OUT of **relax_kw so it can't leak onto the coarse
   stage); `benchmark.py --momentum F`. The level-dict pass-through (lvl_kw ->
   relax) carries it with no extra plumbing. Test
   test_momentum_lands_on_post_split_levels_only locks the invariant (coarse never
   gets it; fine_kw without a momentum key doesn't clobber the injected value).

   **Absolute icr0 stop `Icr0Stop` BUILT + ANALYSED (2026-06-19), 68/68 tests.**
   User-driven decision: do NOT stop on absolute frms (it scales with resolution/
   mass/h and the density gradient, so a fixed value doesn't travel). Stop on
   **icr0**, which IS the maximum step as a FRACTION of h (move = icvel·icr0·h,
   icvel normalised so max component ≤ 1), hence resolution/gradient/N-invariant —
   "stop once even the most-active particle moves < threshold·h". `stepper.Icr0Stop
   (threshold=3e-3, window=50)`: a trailing-MEDIAN smooth over the icr0 sawtooth
   (the controller grows/shrinks icr0 every iter on the sign of the frms change, so
   a BARE crossing is a noisy one-sample event — e.g. β=0.3 dips below 3e-3 at it
   500 then bounces back to 3.1e-3 at 750). Wired `relax(icr0_stop=0.0)` as the
   TOP-priority polish finish (stop_reason='icr0_abs'); per-level in default_levels/
   cascade_levels/relax_progressive (POST-SPLIT fine stage only — the coarse
   box-amplified icr0 is in box units, an absolute h-fraction threshold doesn't
   apply there) + `benchmark.py --icr0-stop F`. Ties to #6d: in the refine stage
   icr0 self-regulates to the critical-step SETPOINT, so icr0<3e-3 == "the critical
   step has shrunk to a negligible h-fraction == the glass has reached its
   force-equilibrium stiffness".

   **THE UNIFYING RESULT — "how good do you want your glass" is ONE knob system.**
   Read off the mhd50 fine sweep (iters to a smoothed frms target, w=25):

   | frms target |  β=0 | β=0.3 | β=0.5 | β=0.7 | best (speedup) |
   |-------------|------|-------|-------|-------|----------------|
   | 5.0e-3 loose|  250 |  169  |  138  | **131** | β=0.7 (1.9×) |
   | 3.4e-3      |  431 |  267  |  208  | **188** | β=0.7 (2.3×) |
   | 2.5e-3      |  715 |  429  | **315**| 315   | 0.5/0.7 (2.3×) |
   | 2.0e-3      | 1032 |  623  | **458**| 548   | β=0.5 (2.25×) |
   | 1.6e-3 tight| >1500|  887  | **671**| 938   | β=0.5 |

   The three things we kept circling — momentum β, the icr0 stop threshold, and
   "glass quality" — are ONE system: **the icr0 stop threshold IS the quality dial**
   (icr0<3e-3 ≈ a frms ~2-3.4e-3 glass; tighten to icr0<1e-3 for ~1.6e-3), and
   **momentum buys ~2× fewer iters to whatever level you pick.** CROSSOVER: for a
   LOOSE glass (frms ≳ 3e-3) higher β wins (β=0.7 2.3×, and icr0<3e-3 catches it
   exactly there); for a TIGHT glass (frms ≤ 2e-3) **β=0.5 wins and β=0.7 LOSES**
   (1.6e-3 needs β=0.7 938 it vs β=0.5 671 — the overshoot/icr0-collapse that was
   harmless for a loose glass now stalls the deep tail). β=0.3 = monotone safe
   middle (~1.5-1.67× everywhere, no collapse risk).

   **The β=0.7 icr0-collapse, correctly diagnosed (NOT "premature").** At its
   icr0<3e-3 fire (it 175) β=0.7 is actually the BEST glass of all betas AT THAT
   ITER (frms 3.40e-3 / denserr 0.0128, lowest of all) — it descended fastest. The
   stop fired right when it reached frms 3.4e-3 (it hits 3.4e-3 at it 188 ≈ the fire
   point), so for a 3.4e-3-quality target it is NOT premature, just 2.7× faster than
   baseline's same-frms stop (it 479). The real artefact: at it 175 β=0.7 has icr0
   2.44e-3 while β=0.5 at the SAME convergence (frms 3.67e-3) still has icr0 4.71e-3
   — momentum shrank icr0 ahead of the actual stiffening, so icr0 UNDER-reports the
   remaining headroom and caps β=0.7's reachable depth via this gate at baseline
   level. β=0.3/0.5 fire DEEPER (frms 2.0-2.2e-3) precisely because their icr0 stays
   a faithful convergence proxy (no collapse). Per-β denserr at the icr0<3e-3 fire:
   β=0 it479 +12% above its asymptote, β=0.3 it501 +7%, β=0.5 it429 +7%, β=0.7 it175
   +14%; p99 (interface tail, kernel-limited) is at asymptote for all.

   **icr0<3e-3 vs the active-set stop (same physics, different distribution point).**
   Both ask "are steps now a negligible fraction of h" (all in h-units). icr0<3e-3
   watches the single GLOBAL-MAX step (≈icr0), smoothed over time → aggressive,
   one knob, fires ~it 480 at denserr +7-12%. The active-set (move_frac 3e-4 / eps
   0.08) watches the FRACTION of particles still moving > move_frac·h, smoothed over
   BOTH population and time → conservative, two knobs, fires ~it 1208 at denserr
   +5.6% (≈2.5× more iters). The log-normalisation freezes the quiet bulk ~1000×
   before the interface, so active_frac ≈ "what fraction of the interface shell is
   still working" and decays smoothly to 0. Two points on one speed↔thoroughness
   axis. icr0<3e-3 is trustworthy despite watching only the loudest particle because
   of the setpoint physics (#6d): the loudest particle's step falling below 3e-3·h
   means the whole landscape stiffened, not one quiet particle fooling you.

   **Recommendation: β=0.5 default** (best-or-tied at every target ≤2.5e-3, no
   collapse downside; fires EARLIER than baseline at the icr0 stop AND at better
   frms), quality set by the icr0 threshold (3e-3 fast / 1e-3 tight). β=0.7 only for
   a deliberately quick loose glass; β=0.3 if a monotone no-surprises run is wanted.
   Sweep harness `_scratch/sweep_momentum.py` (GLASS_BETAS env override) + overlay
   plot `_pipeline/analysis/momentum_sweep_fine.png`.

   **NEXT:** (1) run the end-to-end progressive benchmark (mhd50 + sedov) with
   `--momentum 0.5 --icr0-stop 3e-3` to bank the fewer-iters win as WALL-CLOCK and
   confirm final-glass quality at full convergence (the sweep was fixed-budget,
   stops OFF); compare against `--active-stop` for the conservative finish; (2)
   confirm β + threshold travel to intermediate resolutions / sedov uniform; (3)
   consider auto-detecting the β≥0.7 icr0 collapse or just cap the default at 0.5;
   (4) Nesterov only if a higher β is ever wanted. Uncommitted: glassgen_python
   still untracked in this clone.
