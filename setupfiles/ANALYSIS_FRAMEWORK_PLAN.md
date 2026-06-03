# Analysis Framework — Plan & Progress

> **Restart instructions (if you hit the token limit):**
> Start a new Claude Code session in the `testsuite/` repo and say:
> *"Read `setupfiles/ANALYSIS_FRAMEWORK_PLAN.md` and continue implementing from
> the Progress checklist. We're building `IC_analysis_framework.py` and wiring
> **kh** in first."*
> That reads one file instead of re-exploring the whole codebase. The Progress
> checklist below is the source of truth for what's done vs. pending — keep it
> updated as you go.

## Progress checklist

- [ ] 1. This durable plan/restart doc committed in repo (`setupfiles/ANALYSIS_FRAMEWORK_PLAN.md`)
- [ ] 2. `setupfiles/IC_analysis_framework.py` — shared driver module
- [ ] 3. `setupfiles/IC_analysis_kh.py` — rewritten onto the framework
- [ ] 4. `runtest.sh` — `-doreference` (`-r`) flag + analysis hook
- [ ] 5. Verified end-to-end (parity, save-reference, auto-discover+residuals, no-reference, flag defaults)

**Follow-up (separate change, not this one):** convert `IC_analysis_blob.py` and
`IC_analysis_cosmowave.py` using kh as the template; add numeric pass/fail gate
(`--tol` + non-zero exit) when CI lands.

---

## Context

The repo runs ~39 simulation test cases via `runtest.sh`, each built by an
`IC_setup_<test>.py` script with parameters passed as `-a`…`-h` flags. Three
post-processing scripts exist — `IC_analysis_blob.py`, `IC_analysis_kh.py`,
`IC_analysis_cosmowave.py` — but they are ad-hoc: each duplicates argument
parsing, hardcodes the physical parameters of the run (e.g. blob hardcodes
`clouddens=100, mach=1.5`, kh hardcodes `mode/boundval`), and they do not know
what parameters generated the run nor compare against a saved baseline.

Goal: a **general analysis framework** so every test can be analyzed the same way.
An analysis run must:

1. Take the directory (or directories) to analyze as positional arguments.
2. Take the **same `-a`…`-h` parameters with the same defaults** as that test's
   `runtest.sh` setup function (no sidecar, no `.param` parsing — the operator
   declares what was run, analysis derives physics from them).
3. **Always** overlay/compare against the **analytic solution** for that test.
4. Support a **reference (regression baseline)** workflow:
   - `runtest.sh -doreference` (switch `-r`): after a sim completes, run the
     analysis and **save** the computed metric series as `<runname>_reference`.
   - A normal analysis run **auto-discovers** a matching `*_reference` (same
     test/res/params/code via the run-folder name), plots input + reference +
     analytic together, and prints **residuals of input vs. reference**.
   - CI/CD (later): bless baselines once with `-r`; after a code change, rerun
     and residuals are the regression signal.

## Decisions (confirmed with user)

- Param input: `-a`…`-h` flags mirroring `runtest.sh` defaults exactly.
- Analytic comparison: always on.
- Reference: saved by `runtest.sh -r`; auto-discovered + overlaid on normal runs;
  report residuals (input vs reference). **No pass/fail threshold yet.**
- First test wired in: **kh**.
- Output now: residuals + plot, no exit-code gate.

## Files

- **New:** `setupfiles/IC_analysis_framework.py` — shared driver, builds on
  existing helpers in `setupfiles/IC_analysis_general.py` (`loaddata`,
  `find_files`, `try_aux`, `read_bfield`, `sin_cos_amp`, `amp_mag`,
  `setup_rcparams`, `finish_figure`).
- **Rewrite:** `setupfiles/IC_analysis_kh.py`.
- **Edit:** `runtest.sh` — `-r` flag + analysis hook.
- Leave `IC_analysis_blob.py` / `IC_analysis_cosmowave.py` unchanged for now.

## 1. `IC_analysis_framework.py`

Test-agnostic driver. A test plugs in via: its `-a`…`-h` param spec
(letter → name + default, matching `runtest.sh`), an `analyze` callback, an
`analytic` callback, and metric metadata.

- `build_parser(param_spec)`: positional `inputs` (dir/prefix/glob → `find_files`);
  one option per spec entry as `-a/--<name>` … `-h/--<name>` with runtest default;
  plus `--save <prefix>`, `--save-reference`, `--reference <path>`,
  `--no-analytic`, `--labels`.
- `run_series(files, analyze, params)`: same `find_files`+`loaddata` loop as the
  current `plot_series_for_pattern`; call `analyze(tgdata, time, params, snap)`
  per snapshot; return ordered table `{"x": ..., "<metric>": ...}`.
- Reference I/O (JSON, numpy-friendly): `save_reference(table, path)`,
  `reference_name_for(input_dir)` (= basename + `_reference`),
  `find_reference(input_dir, explicit=None)`.
- `residuals(input_table, ref_table, metric)`: `np.interp` input onto reference
  x-grid; return `max|Δ|` and L2. Printed, not thresholded.
- `plot_compare(input_tables, ref_table, analytic_fn, metric_meta, save)`: per
  metric, plot input + reference + analytic; reuse `setup_rcparams`/`finish_figure`.

## 2. `IC_analysis_kh.py`

Mirror `runtest.sh` `setup_kh` (runtest.sh:344-367) flag-for-flag:

| flag | name    | default     | used? |
|------|---------|-------------|-------|
| -a   | rhodiff | 2.0         | yes (density contrast) |
| -b   | mach    | sqrt(6)/5   | yes (analytic timing) |
| -c   | smooth  | 1           | passthrough |
| -d   | B0      | 0.0         | passthrough |
| -e   | Bdir    | 1           | passthrough |
| -f   | avisc   | 0.01        | passthrough |
| -g   | bvisc   | 2           | passthrough |
| -h   | steps   | 200         | passthrough |

- Keep `analyze_kh(tgdata)` physics (mode amplitude `M` via sin/cos
  exponential-weighted projection; top-5% kinetic energy). Geometric constants
  `mode=4*pi`, `boundval=0.25` stay (documented). Wrap as `analyze` callback →
  `{"x": time, "M": M, "Ekinmax": Ekinmax}`.
- `analytic`: existing growth line `0.01*exp(pi*(t-0.35))` over the linear window
  — overlay on the `M` panel.
- Metric metadata: `M` (semilogy, "Linear mode amplitude"), `Ekinmax` (semilogy,
  "Max kinetic energy").
- `main()` ~15 lines using the framework.
- Parity: `python IC_analysis_kh.py <dir>` (default flags) reproduces today's plot.

## 3. `runtest.sh`

- `getopts` (runtest.sh:111) is `":ia:b:c:d:e:f:g:h:"`; add no-arg flag `r` →
  `doreference=1`. (bash getopts can't parse literal `-doreference`; switch is
  `-r`, documented in usage/README.)
- After sim run (after runtest.sh:1107): if `doreference=1` AND
  `setupfiles/IC_analysis_${testsim}.py` exists, run:
  `python ${analysisdir}/IC_analysis_${testsim}.py "$runname" <-a..-h from args[]> --save-reference`
  → writes `${runname}_reference` in `test_cases/$testsim/`.
  Forward `args[0..7]` as `-a..-h` so analysis params always match the IC params.
- Guard on file-exists so the other ~36 tests are unaffected until converted.

## Verification

1. **Parity:** kh run dir under `test_cases/kh/` (or `./runtest.sh kh 64 128 2 0`);
   `python setupfiles/IC_analysis_kh.py <dir>` matches pre-refactor M plot.
2. **Save reference:** `./runtest.sh -r kh 64 128 2 0` → `test_cases/kh/<runname>_reference` exists with the series.
3. **Auto-discover + residuals:** rerun analysis with reference present → overlays
   input+reference+analytic, prints `M: max|res|=…  L2=…`. Exit 0.
4. **No-reference path:** dir with no matching reference → input+analytic plot,
   notes "no reference found", no crash.
5. **Flag defaults:** no `-a..-h` vs explicit `-a 2.0 -b …` → defaults equal
   `setup_kh`, overrides apply.
