# Statistical Test Fit for H/Z to Y(nS) + Photon

## Setup

One should first setup [Combined](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/).

```
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

## Set env

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
cd statistical_test_fit
```

## Usage

```
python3 main.py"
```
