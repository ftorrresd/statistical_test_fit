# Statistical Test Fit for H/Z to Y(nS) + Photon

## Overview

This repository contains RooFit-based background studies for `H/Z -> Upsilon(nS) + gamma`.

There are two top-level workflows:

- `pseudodata.py`: toy studies for 1D and 2D background-model selection
- `realdata.py`: real-data 2D workflow with auxiliary fits and `RooMultiPdf` construction

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
python3 pseudodata.py --fits-to-run 2d --events 1000 --seed 1 --relax-strict-mode
```

CLI options:

- `--fits-to-run {1d,2d,all}`: choose which pseudodata workflow to run; default is `all`
- `--events INT`: number of generated events; default is `10000`
- `--cheb STR`: comma-separated Chebychev coefficients for the generated 1D background; default is `-0.2,0.2,-0.1`
- `-s, --seed INT`: random seed for RooFit generation
- `--nbins INT`: binning used in control and summary plots; default is `60`
- `--relax-strict-mode`: if strict family selection fails, fall back to the best fit-quality-passing candidate instead of aborting

### Real-Data Workflow

Default run:

```bash
python3 realdata.py
```

Useful examples:

```bash
python3 realdata.py --use-cache
python3 realdata.py --nbins 60
python3 realdata.py --use-cache --relax-strict-mode
```

CLI options:

- `--nbins INT`: binning used in control and summary plots; default is `60`
- `--use-cache`: reuse cached fit parameters and normalization JSON files instead of recomputing them
- `--relax-strict-mode`: if strict family selection fails, fall back to the best fit-quality-passing candidate instead of aborting

## Workflow Notes

- Both driver scripts delete files under `plots/` at startup, keeping only `.gitkeep`.
- `realdata.py` deletes repo-root `*.json` cache files unless `--use-cache` is passed.
- Strict model-family selection is enabled by default. A bad family raises a `RuntimeError` after printing a large warning block.
- `--relax-strict-mode` keeps the workflow moving by selecting the best fit-quality-passing candidate even if no candidate passes the chi-square threshold.

## Outputs

Typical outputs include:

- `plots/fit_1d/`
- `plots/fit_2d/`
- `plots/fit_2d_data/`
- `upsilon_model_params.json`
- `resonant_background_fit_HiggsDalitz.root`
- `resonant_background_fit_ZGamma.root`

The real-data workflow also produces and reuses normalization and resonant-background cache files in the repo root.

## Verification

Minimal syntax check:

```bash
python3 -m compileall pseudodata.py realdata.py statistical_test_fit
```

Focused smoke commands:

```bash
python3 pseudodata.py --fits-to-run 1d --events 1000 --seed 1 --relax-strict-mode
python3 pseudodata.py --fits-to-run 2d --events 1000 --seed 1 --relax-strict-mode
python3 pseudodata.py --events 5000 --seed 42 --nbins 60 --relax-strict-mode
python3 realdata.py --use-cache --relax-strict-mode
```

## Environment Notes

- `statistical_test_fit/normalization_fit.py` uses ROOT plotting plus NumPy. Avoid reintroducing `matplotlib` casually in this CMSSW environment; NumPy / matplotlib ABI mismatches are a real failure mode here.
- The exponential family depends on `RooExpPoly`, which is not available in every ROOT build used with this repo.
