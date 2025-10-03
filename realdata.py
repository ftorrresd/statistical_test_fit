import argparse
import os
from enum import Enum

from ROOT import (  # type: ignore
    RooFit,  # type: ignore
    RooMsgService,  # type: ignore
    gROOT,  # type: ignore
    gSystem,  # type: ignore
)

from statistical_test_fit import run_fit_1d, run_fit_2d, run_fit_2d_data

gROOT.SetBatch(True)
RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)


def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(
        description="Background-only RooFit with a blinded signal region on real Data."
    )
    parser.add_argument("--nbins", type=int, default=60)

    args = parser.parse_args()

    # clear plots dir
    os.system("mkdir -p plots")
    os.system("mkdir -p plots/fit_1d")
    os.system("mkdir -p plots/fit_2d")
    os.system("mkdir -p plots/fit_2d_data")
    os.system(r'find plots -type f ! -name ".gitkeep" -delete')

    gSystem.Load("libRooFit")
    gSystem.Load("libHiggsAnalysisCombinedLimit")

    run_fit_2d_data(args)


if __name__ == "__main__":
    main()
