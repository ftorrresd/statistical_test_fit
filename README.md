# Statistical Test Fit for H/Z to Y(nS) + Photon

## Overview

This repository contains RooFit-based background studies for `H/Z -> Upsilon(nS) + gamma`.

There are four top-level workflows:

- `pseudodata.py`: toy studies for 1D and 2D background-model selection
- `signal.py`: standalone signal-model fits for the six Run2 `H/Z -> Upsilon(nS) + gamma` samples
- `resonant_background.py`: standalone resonant-background and control-region normalization workflow
- `realdata.py`: real-data dimuon fit plus 2D sideband workflow with final `RooMultiPdf` construction

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
python3 pseudodata.py
```

Important behavior:

- Default `--fits-to-run` is `all`
- `all` runs 1D first and 2D second in the same process
- If the 1D step raises, the 2D step does not start

Useful examples:

```bash
python3 pseudodata.py --fits-to-run 1d
python3 pseudodata.py --fits-to-run 2d
python3 pseudodata.py --fits-to-run all --events 20000 --seed 42 --nbins 60
python3 pseudodata.py --fits-to-run 2d --events 1000 --seed 1
python3 pseudodata.py --fits-to-run 2d --events 1000 --seed 1 --strict-mode
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
python3 signal.py
```

CLI options:

- `--nbins INT`: binning used in the saved signal projection plots; default is `60`

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

This repo uses the existing project ranges for the observables:

- `boson_mass`: `57` to `200` GeV
- `upsilon_mass`: `8` to `12` GeV

Saved plots are zoomed to more useful windows around the Z/H peak and the selected Upsilon state.

### Resonant-Background Workflow

Run the resonant-background workflow on its own:

```bash
python3 resonant_background.py
python3 resonant_background.py --use-cache
```

CLI options:

- `--nbins INT`: binning used in the saved resonant-background plots; default is `60`
- `--use-cache`: reuse `resonant_background_model_Z_params.json` and `NormParams_CR*.json` instead of recomputing them

### Real-Data Workflow

Default run:

```bash
python3 realdata.py
```

Useful examples:

```bash
python3 realdata.py --use-cache
python3 realdata.py --nbins 60
python3 realdata.py --use-cache --strict-mode
```

CLI options:

- `--nbins INT`: binning used in control and summary plots; default is `60`
- `--use-cache`: reuse cached dimuon fit parameters (`upsilon_model_params.json`) instead of recomputing them
- `--strict-mode`: abort if a family fails strict selection; the default is relaxed fallback to the best fit-quality-passing candidate

`realdata.py` now only runs the real-data pieces:

1. dimuon non-correlated fit
2. real-data sideband fit and `RooMultiPdf` construction

Signal and resonant-background preparation are standalone workflows and are no longer called by `realdata.py`:

- `python3 signal.py`
- `python3 resonant_background.py --use-cache`

## Workflow Notes

- `pseudodata.py` and `realdata.py` delete files under `plots/` at startup, keeping only `.gitkeep`.
- `signal.py` only clears `plots/signal_fit/`.
- `resonant_background.py` only clears `plots/resonant_background/`.
- `realdata.py` deletes repo-root `*.json` cache files unless `--use-cache` is passed.
- `--use-cache` on `realdata.py` now only affects the dimuon JSON-backed parameter reuse (`upsilon_model_params.json`).
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

`realdata.py` only produces and reuses the dimuon cache in the repo root:

- `upsilon_model_params.json`

## Verification

Minimal syntax check:

```bash
python3 -m compileall pseudodata.py realdata.py signal.py resonant_background.py statistical_test_fit
```

Focused smoke commands:

```bash
python3 pseudodata.py --fits-to-run 1d --events 1000 --seed 1
python3 pseudodata.py --fits-to-run 2d --events 1000 --seed 1
python3 pseudodata.py --events 5000 --seed 42 --nbins 60
python3 signal.py
python3 resonant_background.py --use-cache
```

## Environment Notes

- `statistical_test_fit/normalization_fit.py` uses ROOT plotting plus NumPy. Avoid reintroducing `matplotlib` casually in this CMSSW environment; NumPy / matplotlib ABI mismatches are a real failure mode here.
- The exponential family depends on `RooExpPoly`, which is not available in every ROOT build used with this repo.
