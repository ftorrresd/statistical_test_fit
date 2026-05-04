# Statistical Test Fit for H/Z to Y(nS) + Photon

## Overview

This repository contains RooFit-based background studies for `H/Z -> Upsilon(nS) + gamma`.

There are four top-level workflows:

- `scripts/pseudodata.py`: toy studies for 1D and 2D background-model selection
- `scripts/signal.py`: standalone signal-model fits for the six Run2 `H/Z -> Upsilon(nS) + gamma` samples
- `scripts/resonant_background.py`: standalone resonant-background and control-region normalization workflow
- `scripts/non_resonant_background.py`: renamed from `realdata.py`; runs the dimuon fit plus the real-data non-resonant sideband workflow and final `RooMultiPdf` construction

There is also one packaging step for Combine inputs:

- `scripts/build_bundled_workspace.py`: writes one shared RooWorkspace plus one single simultaneous parametric datacard covering `H_1S`, `H_2S`, `H_3S`, `Z_1S`, `Z_2S`, `Z_3S`, `non_resonant_bkg`, `resonant_H_bkg`, and `resonant_Z_bkg`

There is no `main.py` in this repository.

## Setup

This code is meant to run inside a Combine-enabled CMSSW area.

Example setup:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git -c advice.detachedHead=false clone --depth 1 --branch v10.2.1 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
scramv1 b clean
scramv1 b
cd ../..
git clone git@github.com:ftorrresd/statistical_test_fit.git
cd statistical_test_fit
```

For later sessions:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
cd statistical_test_fit
```

The driver scripts load these libraries at runtime:

- `libRooFit`
- `libHiggsAnalysisCombinedLimit`

## How To Run

Run from the repository root.

### Pseudodata Workflow

Default run:

```bash
python3 scripts/pseudodata.py
```

Important behavior:

- Default `--fits-to-run` is `all`
- `all` runs 1D first and 2D second in the same process
- If the 1D step raises, the 2D step does not start

Useful examples:

```bash
python3 scripts/pseudodata.py --fits-to-run 1d
python3 scripts/pseudodata.py --fits-to-run 2d
python3 scripts/pseudodata.py --fits-to-run all --events 20000 --seed 42 --nbins 60
python3 scripts/pseudodata.py --fits-to-run 2d --events 1000 --seed 1
python3 scripts/pseudodata.py --fits-to-run 2d --events 1000 --seed 1 --strict-mode
```

CLI options:

- `--fits-to-run {1d,2d,all}`: choose which pseudodata workflow to run; default is `all`
- `--events INT`: number of generated events; default is `10000`
- `--cheb STR`: comma-separated Chebychev coefficients for the generated 1D background; default is `-0.2,0.2,-0.1`
- `-s, --seed INT`: random seed for RooFit generation
- `--nbins INT`: binning used in control and summary plots; default is `60`
- `--strict-mode`: abort if a family fails strict selection; the default is relaxed fallback to the best fit-quality-passing candidate

### Signal Workflow

Run the standalone signal-model fits:

```bash
python3 scripts/signal.py
```

CLI options:

- `--nbins INT`: binning used in the saved signal projection plots; default is `60`
- `--workers INT`: number of worker processes for the six independent signal fits; default is `min(os.cpu_count(), 6)`

Inputs used for the signal model:

- `inputs/mass_H_HToUps1SG_Run2.root`
- `inputs/mass_H_HToUps2SG_Run2.root`
- `inputs/mass_H_HToUps3SG_Run2.root`
- `inputs/mass_Z_ZToUpsilon1SGamma_Run2.root`
- `inputs/mass_Z_ZToUpsilon2SGamma_Run2.root`
- `inputs/mass_Z_ZToUpsilon3SGamma_Run2.root`

The signal RooFit model follows `from_mauricio/HZUpsilonPhotonRun2Statistics/signal_modeling.py` in shape:

- boson: `RooCBShape + Gaussian` combined with `RSUM`
- upsilon: `RooDoubleCB`
- full signal model: `PROD(signal_model_boson, signal_model_upsilon)`

The signal fit uses signal-specific boson windows while keeping the repo-wide Upsilon range:

- Z signal boson fit: `70` to `120` GeV
- Higgs signal boson fit: `100` to `150` GeV
- `upsilon_mass`: `8` to `12` GeV

Saved plots are zoomed to more useful windows around the Z/H peak and the selected Upsilon state.
`signal.py` lists the six jobs before starting and updates a progress bar as each sample finishes.

### Resonant-Background Workflow

Run the resonant-background workflow on its own:

```bash
python3 scripts/resonant_background.py
python3 scripts/resonant_background.py --skip-cache
```

CLI options:

- `--nbins INT`: binning used in the saved resonant-background plots; default is `60`
- `--skip-cache`: recompute `resonant_background_model_Z_params.json` and `NormParams_CR*.json` instead of using the default cache-backed behavior
- `--workers INT`: number of worker processes for the independent resonant-background stages

`resonant_background.py` runs in stages:

1. Higgs workspace, Z workspace, and Z boson fit in parallel
2. CR1-CR4 normalization fits in parallel
3. normalization extrapolation in the parent process

Each stage prints its job list first and updates a progress bar as jobs complete.

### Non-Resonant Background Workflow

Default run:

```bash
python3 scripts/non_resonant_background.py
```

Useful examples:

```bash
python3 scripts/non_resonant_background.py --use-cache
python3 scripts/non_resonant_background.py --nbins 60
python3 scripts/non_resonant_background.py --use-cache --strict-mode
```

CLI options:

- `--nbins INT`: binning used in control and summary plots; default is `60`
- `--workers INT`: number of worker processes for the parallel background-candidate scan
- `--use-cache`: reuse cached dimuon fit parameters (`upsilon_model_params.json`) instead of recomputing them
- `--strict-mode`: abort if a family fails strict selection; the default is relaxed fallback to the best fit-quality-passing candidate

`scripts/non_resonant_background.py` is the renamed `realdata.py` entrypoint and now only runs the real-data pieces:

1. dimuon non-correlated fit
2. real-data sideband fit and `RooMultiPdf` construction

`scripts/non_resonant_background.py` runs the dimuon fit serially, then fits the independent background candidates in parallel and finally rebuilds the winning PDFs in the parent process for plotting and `RooMultiPdf` assembly.
It prints the job list for each stage and shows a progress bar as jobs complete.

Signal and resonant-background preparation are standalone workflows and are no longer called by `non_resonant_background.py`:

- `python3 scripts/signal.py`
- `python3 scripts/resonant_background.py`

### Single Datacard Builder

After the signal, resonant-background, and non-resonant-background inputs exist, build the single simultaneous Combine workspace and card with:

```bash
python3 scripts/build_bundled_workspace.py
```

Useful example:

```bash
python3 scripts/build_bundled_workspace.py --strict-mode --output-dir datacards
```

CLI options:

- `--output-dir PATH`: directory where `workspace.root`, `datacard.txt`, and a generated `README.html` summary are written
- `--workspace-name STR`: RooWorkspace name stored in the ROOT file; default is `combined_workspace`
- `--workspace-file-name STR`: bundled ROOT filename; default is `workspace.root`
- `--datacard-file-name STR`: datacard filename; default is `datacard.txt`
- `--strict-mode`: require strict non-resonant family selection instead of the default relaxed selection
- `--skip-validation`: skip the validation step
- `--signal-mass-label STR`: mass label passed to `text2workspace.py` and `combine` during validation; default is `125`

By default the builder also validates the produced card by running:

- `text2workspace.py`
- `combine -M AsymptoticLimits --run blind`

Validation outputs are written next to the datacard as `validation.log`, `validation_workspace.root`, and the usual `higgsCombine*.root` file when the smoke run succeeds.

### Limit Driver

After `scripts/build_bundled_workspace.py` has produced `datacards/datacard.txt` and `datacards/workspace.root`, run the expected-limit automation with:

```bash
python3 scripts/limits.py
```

Default behavior:

- Builds the six independent process-POI scheme plus `z_grouped` and `h_grouped` by default
- In `z_grouped`, all Z signal processes share one POI while the H signal strengths are profiled individually; `h_grouped` does the reverse
- Runs both `AsymptoticLimits --run blind` and `HybridNew --LHCmode LHC-limits`
- Uses HybridNew expected quantiles `0.16,0.5,0.84` by default
- Does not pass `--dataset` or `--bypassFrequentistFit` to `HybridNew`
- By default, the POI range is always `0,1000000`, including in quick mode
- Optional non-default `--hybrid-range-from-asymptotic` runs asymptotic limits first and sets each HybridNew target/quantile range to `0,min(1000000,10*r_asymp)`
- In asymptotic-derived range mode, retryable HybridNew minimization warnings immediately trigger 10x r-range retries with `--cminDefaultMinimizerStrategy 2` when that job finishes, up to `--hybrid-range-max-retries` (default `3`); exhausted warning retries are reported prominently but return code `0` jobs remain successful
- Optional non-default `--quick` mode adds HybridNew options `--rRelAcc 0.10 --rAbsAcc 10 --clsAcc 0.02 -T 100`
- Runs all ready commands in parallel within each dependency wave; adaptive HybridNew retries are submitted per completed job; use `--workers N` to cap concurrency
- Prints the complete job manifest before each dependency wave starts
- Clears the run directory and each per-job working directory before staging inputs and running commands
- Creates one working directory per Combine call with staged inputs, ROOT outputs, `command.txt`, `stdout.txt`, `stderr.txt`, `result.json`, and a dark `summary.html`
- Reports each job duration as `DD:HH:MM:SS:MS` in the live log, per-job JSON/HTML, and aggregate summary
- Writes an aggregate `blind_limits_summary.json` and dark `README.html` under `datacards/limits/run_YYYYmmdd_HHMMSS/`

Useful examples:

```bash
python3 scripts/limits.py --workers 8
python3 scripts/limits.py --quick --workers 8
python3 scripts/limits.py --poi-scheme six --methods asymptotic
python3 scripts/limits.py --poi-scheme z_grouped --methods hybrid --hybrid-toys 1000 --cls-acc 0.005
python3 scripts/limits.py --methods both --hybrid-range-from-asymptotic
python3 scripts/limits.py --poi-scheme grouped --methods hybrid
python3 scripts/limits.py --run-name test_blind_limits
```

### Branching-Fraction Limit Table

After `scripts/limits.py` has written a `blind_limits_summary.json`, make a LaTeX table of theory branching fractions and expected/observed branching-fraction limits with:

```bash
python3 scripts/branching_fraction_table.py --summary datacards/limits/<run>/blind_limits_summary.json
```

Default behavior:

- Writes grouped H, grouped Z, and individual six-POI tables for both `asymptotic` and `hybrid_lhc` methods by default
- Uses one row for each grouped H/Z table, with the theory branching fraction equal to the sum of the three states in that boson group
- Includes raw observed/expected signal-strength columns plus branching-fraction limit columns
- Marks upper-limit cells with `<`
- Multiplies each signal-strength limit by the corresponding theory branching fraction
- Represents the one-standard-deviation expected band as deltas, e.g. `99.8^{+10.1}_{-0.9}`
- Writes `branching_fraction_limits.table.tex`, `branching_fraction_limits.tex`, and `branching_fraction_limits.json` into a `tables/` directory next to the summary
- Leaves observed cells as `--` when no observed limit is available in the JSON

To compile the standalone LaTeX file to PDF, use Singularity with the Docker LaTeX image:

```bash
python3 scripts/branching_fraction_table.py \
  --summary datacards/limits/<run>/blind_limits_summary.json \
  --compile-pdf
```

Useful examples:

```bash
python3 scripts/branching_fraction_table.py --method asymptotic
python3 scripts/branching_fraction_table.py --scheme z_grouped
python3 scripts/branching_fraction_table.py --scheme h_grouped
python3 scripts/branching_fraction_table.py --compile-pdf --latex-image docker://ghcr.io/xu-cheng/texlive-full:latest
```

## Workflow Notes

- `scripts/pseudodata.py` only clears `plots/fit_1d`, `plots/fit_2d`, and `plots/fit_2d_data` at startup.
- `scripts/non_resonant_background.py` only clears `plots/fit_1d`, `plots/fit_2d`, and `plots/fit_2d_data` at startup.
- `scripts/signal.py` only clears `plots/signal_fit/`.
- `scripts/resonant_background.py` only clears `plots/resonant_background/`.
- `scripts/non_resonant_background.py` deletes repo-root `*.json` cache files unless `--use-cache` is passed.
- `--use-cache` on `scripts/non_resonant_background.py` now only affects the dimuon JSON-backed parameter reuse (`upsilon_model_params.json`).
- `dimuon_non_correlated()` snapshots `inputs/selected_Run2_dimuon_non_correlated_renamed_branch.root` on every real-data run, so `inputs/` is not effectively read-only.
- Relaxed model-family selection is enabled by default, so a bad family falls back to the best fit-quality-passing candidate after printing a large warning block.
- `--strict-mode` restores abort-on-failure behavior and raises a `RuntimeError` when a family cannot satisfy strict selection.
- `fit_2d_data.py` builds the final `RooMultiPdf` in memory and prints the workspace, but does not write that final workspace / multipdf to disk.
- The 2D toy flow currently generates pure background: `fit_2d.py` passes `0.0` Z and Higgs signal fractions into the generation model even though nominal `z_sigfrac` / `h_sigfrac` values are defined nearby.

## Outputs

Typical outputs include:

- `plots/fit_1d/`
- `plots/fit_2d/`
- `plots/fit_2d_data/`
- `plots/signal_fit/`
- `plots/resonant_background/`
- `upsilon_model_params.json`
- `signal_workspace_H_HToUps1SG_Run2.root`
- `signal_workspace_H_HToUps2SG_Run2.root`
- `signal_workspace_H_HToUps3SG_Run2.root`
- `signal_workspace_Z_ZToUpsilon1SGamma_Run2.root`
- `signal_workspace_Z_ZToUpsilon2SGamma_Run2.root`
- `signal_workspace_Z_ZToUpsilon3SGamma_Run2.root`
- `resonant_background_fit_HiggsDalitz.root`
- `resonant_background_fit_ZGamma.root`

`non_resonant_background.py` only produces and reuses the dimuon cache in the repo root:

- `upsilon_model_params.json`

## Verification

Minimal syntax check:

```bash
python3 -m compileall scripts/pseudodata.py scripts/non_resonant_background.py scripts/signal.py scripts/resonant_background.py scripts/build_bundled_workspace.py statistical_test_fit
```

Focused smoke commands:

```bash
python3 scripts/pseudodata.py --fits-to-run 1d --events 1000 --seed 1
python3 scripts/pseudodata.py --fits-to-run 2d --events 1000 --seed 1
python3 scripts/pseudodata.py --events 5000 --seed 42 --nbins 60
python3 scripts/signal.py
python3 scripts/resonant_background.py --use-cache
python3 scripts/build_bundled_workspace.py --output-dir /tmp/statistical_test_fit_single_card
```

## Environment Notes

- `statistical_test_fit/normalization_fit.py` uses ROOT plotting plus NumPy. Avoid reintroducing `matplotlib` casually in this CMSSW environment; NumPy / matplotlib ABI mismatches are a real failure mode here.
- The exponential family depends on `RooExpPoly`, which is not available in every ROOT build used with this repo.
