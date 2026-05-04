# AGENTS.md

## Environment
- Run from the repo root inside a Combine-enabled CMSSW area. The scripts expect `cmsenv` to already be active before they load RooFit / Combine libraries.
- There is no repo-local `pyproject.toml`, `requirements*.txt`, CI workflow, or pre-commit config. Do not invent a lint, typecheck, or test pipeline.

## Real Entrypoints
- Use `python3 scripts/pseudodata.py` for toy studies, `python3 scripts/signal.py` for standalone signal fits, `python3 scripts/resonant_background.py` for standalone resonant-background fits, `python3 scripts/non_resonant_background.py` for the real-data dimuon-plus-sideband workflow, and `python3 scripts/build_bundled_workspace.py` for Combine packaging.
- The driver scripts mostly parse args, clear workflow-specific outputs, configure ROOT, and dispatch into `statistical_test_fit/`. Most real logic lives in `fit_1d.py`, `fit_2d.py`, `fit_2d_data.py`, `signal_modeling.py`, and `resonant_bkg_modeling.py`.
- `statistical_test_fit/__init__.py` uses lazy `__getattr__` imports on purpose so pseudodata imports do not pull in the heavier real-data stack.
- `scripts/signal.py` intentionally contains a stdlib `signal` shim because the filename shadows Python's `signal` module. Do not remove that shim unless the script is renamed.

## Workflow Order
- The single-card builder depends on persisted outputs from three earlier steps, in this order: `python3 scripts/signal.py`, `python3 scripts/resonant_background.py`, `python3 scripts/non_resonant_background.py`, then `python3 scripts/build_bundled_workspace.py`.
- `scripts/pseudodata.py` defaults to `--fits-to-run all`, which runs 1D first and 2D second in the same process. Use `--fits-to-run 1d` or `--fits-to-run 2d` for focused debugging; a 1D failure prevents the 2D step from starting.
- `scripts/resonant_background.py` uses cache by default. The recompute flag is `--skip-cache`; there is no `--use-cache` flag for this script.
- `scripts/non_resonant_background.py` deletes repo-root `*.json` caches unless `--use-cache` is passed. Use `--use-cache` for iterative work if you want to preserve `upsilon_model_params.json`.
- Relaxed family selection is on by default in the pseudodata and real-data background scans. Add `--strict-mode` when a family should fail the run instead of falling back to the best fit-quality-passing candidate.

## Output Gotchas
- `scripts/pseudodata.py` and `scripts/non_resonant_background.py` only clear `plots/fit_1d`, `plots/fit_2d`, and `plots/fit_2d_data`. `scripts/signal.py` only clears `plots/signal_fit/`, and `scripts/resonant_background.py` only clears `plots/resonant_background/`.
- `dimuon_non_correlated()` overwrites `inputs/selected_Run2_dimuon_non_correlated_renamed_branch.root`; `inputs/` is not effectively read-only.
- `fit_2d_data.py` now writes `non_resonant_background_workspace.root` plus the non-resonant summary JSON. The final simultaneous Combine card is still produced later by `scripts/build_bundled_workspace.py`.
- `scripts/build_bundled_workspace.py` writes directly into `./datacards` by default and does not clear that directory first. Stale legacy per-channel files can coexist with the current single-card outputs until you clean them.
- Use `python3 scripts/clean_outputs.py --dry-run` to inspect generated artifacts and `python3 scripts/clean_outputs.py` to wipe them. That script removes the whole `datacards/` and `plots/` trees plus repo-root generated ROOT/JSON files.
- Treat `datacards/` as generated output, not as the source of truth for the current card layout. The authoritative implementation is `scripts/build_bundled_workspace.py`.

## Combine Packaging
- `scripts/build_bundled_workspace.py` writes one simultaneous card with public processes `H_1S`, `H_2S`, `H_3S`, `Z_1S`, `Z_2S`, `Z_3S`, `non_resonant_bkg`, `resonant_H_bkg`, and `resonant_Z_bkg`.
- Builder validation is on by default. It runs `text2workspace.py` first and then a blind `combine -M AsymptoticLimits` smoke test.
- Unless `--skip-validation` is used, the builder also writes `validation.log`, `validation_workspace.root`, and `higgsCombine*.root` next to `workspace.root`, `datacard.txt`, and the generated `README.md`.
- The non-resonant side of the builder depends on `plots/fit_2d_data/non_resonant_fit_summary.json` containing the strict/relaxed family selections and candidate normalization estimates. If those keys are missing, rerun `python3 scripts/non_resonant_background.py` with current code.

## Code Constraints
- `statistical_test_fit/normalization_fit.py` intentionally uses ROOT plotting plus NumPy; do not reintroduce `matplotlib` casually in this CMSSW environment.
- `RooExpPoly` is guarded in `statistical_test_fit/bkg_model.py`; do not enable the exponential family blindly on every ROOT build.
- Multiprocessing uses `ProcessPoolExecutor` with `spawn` in `parallel_utils.py`. Keep ROOT imports and ROOT configuration inside worker bodies; moving them to module scope breaks spawned workers.
- The 2D toy flow in `fit_2d.py` currently generates pure background: nominal `z_sigfrac` / `h_sigfrac` values exist, but the generation call still passes `0.0` for both signal fractions.

## Verification
- Smallest syntax check: `python3 -m compileall scripts/pseudodata.py scripts/non_resonant_background.py scripts/signal.py scripts/resonant_background.py scripts/build_bundled_workspace.py scripts/clean_outputs.py statistical_test_fit`
- Focused smoke commands:
  - `python3 scripts/pseudodata.py --fits-to-run 1d --events 1000 --seed 1`
  - `python3 scripts/pseudodata.py --fits-to-run 2d --events 1000 --seed 1`
  - `python3 scripts/signal.py`
  - `python3 scripts/resonant_background.py`
  - `python3 scripts/build_bundled_workspace.py --output-dir /tmp/statistical_test_fit_single_card`
  - For a faster builder-only check that skips Combine validation: `python3 scripts/build_bundled_workspace.py --skip-validation --output-dir /tmp/statistical_test_fit_single_card`
