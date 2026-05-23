# AGENTS.md

## Environment
- Run commands from the repo root after `cmsenv` in a Combine-enabled CMSSW area; scripts load `libRooFit` and `libHiggsAnalysisCombinedLimit` at runtime.
- There is no repo-local `pyproject.toml`, `requirements*.txt`, CI workflow, pre-commit config, or task runner. Do not invent lint/type/test commands.

## Entrypoints
- `python3 scripts/prepare_input_files.py` creates `inputs_weighted/`, rewriting MC `weight` branches to the nominal product, removing data `weight*` branches, and copying other inputs unchanged; export `STATISTICAL_TEST_FIT_INPUT_DIR=inputs_weighted` before downstream real-data/MC workflows to use it.
- `python3 scripts/pseudodata.py` runs toy 1D/2D scans; default `--fits-to-run all` runs 1D first and 2D second, so a 1D failure blocks 2D.
- `python3 scripts/signal.py` runs the six signal fits; keep its stdlib `signal` shim because the filename shadows Python's `signal` module.
- `python3 scripts/resonant_background.py` runs resonant backgrounds; cache is on by default and the recompute flag is `--skip-cache`.
- `python3 scripts/non_resonant_background.py` runs the real-data dimuon plus sideband non-resonant workflow; without `--use-cache` it deletes repo-root `*.json` caches.
- `python3 scripts/build_bundled_workspace.py` packages persisted fit outputs into one workspace/card; run it only after signal, resonant, and non-resonant outputs exist.
- `python3 scripts/limits.py`, `python3 scripts/branching_fraction_table.py`, `python3 scripts/bias_study.py`, and `python3 scripts/correlation_study.py` are post-card or diagnostic workflows.
- Driver scripts mostly parse args, clear workflow-specific outputs, configure ROOT, and dispatch into `statistical_test_fit/`; keep `statistical_test_fit/__init__.py` lazy so light workflows do not import the heavier real-data stack.

## Workflow Order
- Full Combine chain with prepared inputs: `scripts/prepare_input_files.py`, export `STATISTICAL_TEST_FIT_INPUT_DIR=inputs_weighted`, then `scripts/signal.py`, then `scripts/resonant_background.py`, then `scripts/non_resonant_background.py`, then `scripts/build_bundled_workspace.py`, then optional `scripts/limits.py` and `scripts/branching_fraction_table.py` or `scripts/bias_study.py`.
- Relaxed family selection is default in pseudodata and real-data background scans; pass `--strict-mode` when a failed family should abort instead of falling back to the best fit-quality-passing candidate.
- `scripts/limits.py` defaults to `limits/`, not `datacards/limits/`; `--run-name` reuses and clears that run directory if it already exists.
- `scripts/branching_fraction_table.py` defaults to the newest `limits/*/blind_limits_summary.json` and compiles PDF through Singularity unless `--skip-compile-pdf` is used.
- `scripts/bias_study.py` defaults to expensive toy jobs (`--toys 1000`, all POI schemes); use `--dataset-strategy asimov`, `--quick`, `--skip-fits`, or `--plot-only` for focused checks.

## Outputs
- `scripts/pseudodata.py` and `scripts/non_resonant_background.py` clear only `plots/fit_1d`, `plots/fit_2d`, and `plots/fit_2d_data`; `scripts/signal.py` clears only `plots/signal_fit`; `scripts/resonant_background.py` clears only `plots/resonant_background`.
- `dimuon_non_correlated()` overwrites `inputs/selected_Run2_dimuon_non_correlated_renamed_branch.root`; `inputs/` is not effectively read-only.
- With `STATISTICAL_TEST_FIT_INPUT_DIR=inputs_weighted`, `dimuon_non_correlated()` writes the renamed dimuon file inside `inputs_weighted/` instead of touching `inputs/`.
- `fit_2d_data.py` writes `non_resonant_background_workspace.root` and `plots/fit_2d_data/non_resonant_fit_summary.json`; the final simultaneous Combine card is still produced by `scripts/build_bundled_workspace.py`.
- `scripts/build_bundled_workspace.py` clears generated `datacards/` contents before writing there by default; custom `--output-dir` locations are not cleared automatically.
- Treat `inputs_weighted/`, `datacards/`, `limits/`, `bias_study/`, `plots/`, repo-root generated ROOT files, and repo-root generated JSON files as outputs; inspect with `python3 scripts/clean_outputs.py --dry-run` before wiping with `python3 scripts/clean_outputs.py`.

## Combine Packaging
- The single-card builder writes public processes `H_1S`, `H_2S`, `H_3S`, `Z_1S`, `Z_2S`, `Z_3S`, `non_resonant_bkg`, `resonant_H_bkg`, and `resonant_Z_bkg`.
- The builder requires `inputs/yields_nevents.json` and derives asymmetric `lnN` rows for `pu_r`, `trg`, `muon_id`, `muon_iso`, `ph_id`, `ele_veto`, `pdf_alpha_s_weight`, and `l1_prefiring`; only `lumi` and `HZ_xs_sc` remain hardcoded.
- The builder requires `inputs/mass_systematics_summary.json` and turns signal `mean_boson` / DCB `sigma_boson` into functions of `CMS_sig_mmg_muon_cor`, `CMS_sig_mmg_photon_E_scale`, and `CMS_sig_mmg_photon_E_smearing` param nuisances.
- Builder validation is on by default and runs only `text2workspace.py datacard.txt -m 125 -o validation_workspace.root`; use `--skip-validation` for workspace/card generation only.
- Builder outputs are `workspace.root`, `datacard.txt`, `bundle_summary.html`, `bundle_summary.json`, plus `validation.log` and `validation_workspace.root` when validation is enabled.
- `bundle_summary.*` is a dark persisted-artifact validation report: it reopens the written workspace and parses the written datacard before summarizing objects, nuisances, and consistency checks.
- The builder requires `plots/fit_2d_data/non_resonant_fit_summary.json` to contain strict/relaxed selections and candidate normalization estimates; if those keys are missing, rerun `scripts/non_resonant_background.py`.
- `from_mauricio/.../combine_helpers.py` is not a runtime input to the builder; legacy `lnN` values are embedded in `scripts/build_bundled_workspace.py`.

## Code Constraints
- `statistical_test_fit/normalization_fit.py` intentionally uses ROOT plotting plus NumPy; do not add matplotlib there casually in this CMSSW environment.
- `RooExpPoly` is guarded in `statistical_test_fit/bkg_model.py`; do not enable the exponential family blindly on every ROOT build.
- Multiprocessing uses `ProcessPoolExecutor` with `spawn`; keep ROOT imports and ROOT configuration inside worker bodies unless you verify spawned workers still start cleanly.
- The 2D toy flow in `fit_2d.py` currently generates pure background: nominal `z_sigfrac` / `h_sigfrac` values exist, but generation still passes `0.0` for both.

## Verification
- Smallest syntax check: `python3 -m compileall scripts statistical_test_fit`
- Focused toy checks: `python3 scripts/pseudodata.py --fits-to-run 1d --events 1000 --seed 1` and `python3 scripts/pseudodata.py --fits-to-run 2d --events 1000 --seed 1`.
- Builder smoke without touching default `datacards/`: `python3 scripts/build_bundled_workspace.py --skip-validation --output-dir /tmp/statistical_test_fit_single_card`.
- Builder validation smoke: `python3 scripts/build_bundled_workspace.py --output-dir /tmp/statistical_test_fit_single_card`.
