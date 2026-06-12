# Analysis smoke suite

End-to-end regression tests for the `setupfiles/IC_analysis_*.py` scripts. Each
test stands up a tiny **synthetic** run on disk (no gasoline, no IC build), then
drives the real analysis CLI as a subprocess and checks it:

1. runs to completion (exit 0) and writes its plot(s);
2. blesses a reference baseline with `--save-reference`;
3. self-compares against that baseline + `--reg-tol` and **passes** (residual ~0);
4. trips `--reg-tol` (exit 2) once its data is perturbed away from the baseline.

This is the captured form of the ad-hoc synthetic runs the analyses were
validated against in development, and the safety net the framework
de-duplication refactors rely on.

## Running

From this directory, `setupfiles/analysis_smoketests/` (no third-party deps;
headless via `MPLBACKEND=Agg`):

```bash
python -m unittest test_analyses
python -m unittest test_analyses.TestSedov -v
```

`pytest setupfiles/analysis_smoketests/` also works if pytest is installed
(plain `unittest.TestCase`s).

**CI / one-shot runner:** `bash setupfiles/analysis_smoketests/run.sh` runs the
whole suite from any working directory (it cd's to the testsuite root, sets
`MPLBACKEND=Agg`, puts `setupfiles/` on `PYTHONPATH`, and exits non-zero on any
failure -- so a CI job is just that one line). Pass a test id to run a subset:
`bash setupfiles/analysis_smoketests/run.sh test_analyses.TestSedov`.

## Files

- `synthdata.py` — writes a synthetic run directory: `make_run(parent, runname,
  snaps, log_series=..., log_params=...)`. Snapshots are big-endian tipsy
  (`<runname>.NNNNN`), aux fields are component/scalar files
  (`.BFieldx/y/z`, `.DivB`, …), and the log is a 52-column `<runname>.log` with
  the `# UNITS/BOX/GAS PHYSICS` header lines (dTuFac, dxPeriod, …). All readable
  by `IC_analysis_general` (loaddata / read_vector_aux / read_log_series /
  read_log_param). Gas columns: `0 mass | 1-3 xyz | 4-6 v | 7 rho | 8 temp |
  9 h | 10 metals | 11 phi`.
- `test_analyses.py` — `AnalysisSmokeBase` does the 4-step bless/gate dance; one
  subclass per covered test sets `test = "<name>"` and `build_run(parent,
  diverge=False)` (the `diverge` data perturbation is what must trip the gate).

## Currently covered (ALL 25 analyses)

- **run_cli metric-vs-time** (framework `--reference`/reg-tol gate): **kh**, **rt**,
  **mti**, **square**, **taylorgreen**, **blob** (cloud-survival mass fractions,
  normalised to t=0).
- **standalone energy/divergence series** from the dense `.log`: **mhdloop**
  (Emag), **orzag** (Emag), **mhdrotor** (Emag+Ekin), **mhddivtest** (divBerr),
  **currentsheet** (Emag/Ekin/Eth + divBerr), **gresho** (Ekin),
  **mhdcollapse** (rho_peak).
- **standalone binned-profile reference** (`binned_profile`): **sedov** (+ a
  `--render` test of the ke_density/thermal_pressure factories), **noh**,
  **zeldovich**, **shocktube**, **evrard**, **polytrope** (the last two use the
  geometric-log-grid / median binning).
- **standalone per-snapshot / profile**: **alfven** (profile + energy-block ref,
  parse_resolution, the `_binned_profile` wrapper), **isentropicvortex** (E_vort),
  **coldkeplerian** (ring sigma_r + L_z), **wave** (sound eigenmode RMS),
  **mhdisowave** (longitudinal pulse vx/drho profile — NOTE: its reference keeps
  empty bins as NaN on a fixed `[-2,2]` grid, so the synthetic particles must
  DENSELY cover `[-2,2]` or `residuals` returns NaN and the gate can't trip).
- **standalone amplitude-vs-scale-factor**: **cosmowave** (case 3 Fourier
  projection of vx/rho/Bz onto the seeded mode; raw amplitudes vs the scale
  factor `a` = snapshot time).

The 6 analyses with a `--convergence` mode (**alfven, gresho, isentropicvortex,
wave, coldkeplerian, mhdisowave**) ALSO have a `test_convergence` that drives that
path on two synthetic resolutions (n_x = 32, 64; `build_run(runname=...)` gives
each a distinct leading `<test><nx>`) and asserts a `_convergence.png` is written
— via the shared `AnalysisSmokeBase.check_convergence` helper.

Every analysis with a CLI is now covered. The 7 runtest setups with NO analysis
script (accdisk, areablob, mhdbalsaravortex, mhdloopshear, mri, planet,
planetcollision) have nothing to test yet.

## Adding a test

Subclass `AnalysisSmokeBase`, set `test`, and implement `build_run` to write a
synthetic run whose analysis is meaningful, with a `diverge=True` branch that
shifts the blessed diagnostic enough to trip `--reg-tol`. Override
`expected_plots` if the script's `--save` filename suffix isn't `*.png`-generic.

**Gotcha — normalised metrics:** when the blessed series is normalised to its
t=0 value (e.g. `E/E0`, `Ekin/Ekin0`, `E_vort/E_vort0`), a uniform scale across
all snapshots leaves the ratio unchanged and will NOT trip the gate. The
`diverge=True` branch must change the series *shape* — perturb only the later
snapshot(s), or change the relative values. Non-normalised metrics (raw divBerr,
spurious Mach, mode amplitudes) trip on a plain scale. The shared `mhd_snaps` /
`lattice3d` / `_grid2d` helpers cover most geometries.
