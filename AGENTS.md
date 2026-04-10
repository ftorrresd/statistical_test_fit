# AGENTS.md

## Environment
- Run from the repo root inside a Combine-enabled CMSSW area. The driver scripts call `gSystem.Load("libRooFit")` and `gSystem.Load("libHiggsAnalysisCombinedLimit")`, so `cmsenv` must already be active.
- No repo-local `pyproject.toml`, `requirements*.txt`, CI workflow, or pre-commit config exists. Do not guess a lint/test pipeline.

## Real Entrypoints
- Use `python3 pseudodata.py` for toy studies and `python3 realdata.py` for the real-data workflow.
- `pseudodata.py` defaults to `--fits-to-run all`, which runs 1D first and 2D second in the same process. Use `--fits-to-run 1d` or `--fits-to-run 2d` when debugging; a 1D failure prevents the 2D step from starting.
- The top-level scripts only parse args, clear outputs, load ROOT libs, and dispatch into `statistical_test_fit/`. Most real logic lives in `statistical_test_fit/fit_1d.py`, `fit_2d.py`, and `fit_2d_data.py`.
- `statistical_test_fit/__init__.py` uses lazy `__getattr__` imports on purpose so pseudodata runs do not pull in the heavier real-data stack at import time.

## Workflow Gotchas
- Both driver scripts delete all files under `plots/` except `.gitkeep` at startup. Do not keep manual artifacts there.
- `run_fit_2d_data()` also deletes repo-root `*.json` caches unless `--use-cache` is passed. Use `python3 realdata.py --use-cache` for iterative work if you want to preserve cached fit parameters.
- Relaxed background-family selection is on by default. Add `--strict-mode` to either driver if a family should abort the run instead of falling back to the best fit-quality-passing candidate.
- `--use-cache` only affects the JSON-backed parameter reuse. `run_fit_2d_data()` still rebuilds the resonant Higgs/Z workspaces and their plots on every run.
- `dimuon_non_correlated()` snapshots `inputs/selected_Run2_dimuon_non_correlated_renamed_branch.root` on every real-data run; `inputs/` is not effectively read-only.
- `fit_2d_data.py` builds the final `RooMultiPdf` in memory and prints the workspace, but does not write that final workspace / multipdf to disk.
- The 2D toy flow in `fit_2d.py` currently generates pure background: the generation model is called with `0.0` Z and Higgs signal fractions even though nominal `z_sigfrac` / `h_sigfrac` values are defined nearby.

## Environment-Specific Code Constraints
- `statistical_test_fit/normalization_fit.py` intentionally uses ROOT plotting plus NumPy; do not reintroduce `matplotlib` casually. This repo runs inside CMSSW Python environments where NumPy / matplotlib ABI mismatches are a real failure mode.
- `RooExpPoly` is imported behind a guard in `statistical_test_fit/bkg_model.py`. The exponential family is not safe to enable blindly on every ROOT build.

## Verification
- The smallest verified syntax check is `python3 -m compileall pseudodata.py realdata.py statistical_test_fit`.
- Focused smoke commands:
  - `python3 pseudodata.py --fits-to-run 1d --events 1000 --seed 1`
  - `python3 pseudodata.py --fits-to-run 2d --events 1000 --seed 1`
  - `python3 realdata.py --use-cache`
