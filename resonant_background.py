import argparse
import os

from ROOT import (  # type: ignore
    RooFit,  # type: ignore
    RooMsgService,  # type: ignore
    gROOT,  # type: ignore
    gSystem,  # type: ignore
)

from statistical_test_fit import run_resonant_background

gROOT.SetBatch(True)
RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)


def main():
    parser = argparse.ArgumentParser(
        description="Standalone resonant-background modeling workflow."
    )
    parser.add_argument("--nbins", type=int, default=60)
    parser.add_argument("--use-cache", default=False, action="store_true")
    args = parser.parse_args()

    os.system("mkdir -p plots")
    os.system("mkdir -p plots/resonant_background")
    os.system(r'find plots/resonant_background -type f ! -name ".gitkeep" -delete')

    gSystem.Load("libRooFit")
    gSystem.Load("libHiggsAnalysisCombinedLimit")

    run_resonant_background(args)


if __name__ == "__main__":
    main()
