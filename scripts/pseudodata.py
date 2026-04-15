import argparse
import sys
from enum import Enum
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _clear_plot_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for child in path.rglob("*"):
        if child.is_file() and child.name != ".gitkeep":
            child.unlink()


from ROOT import (  # type: ignore
    RooFit,  # type: ignore
    RooMsgService,  # type: ignore
    gROOT,  # type: ignore
    gSystem,  # type: ignore
)

from statistical_test_fit import run_fit_1d, run_fit_2d

gROOT.SetBatch(True)
RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)


# Define the Enum for colors
class FitToRun(Enum):
    ONEDIM = "1d"
    TWODIM = "2d"
    ALL = "all"

    # Create a custom type for argparse that converts the string to an Enum member
    @staticmethod
    def arg_type(arg_value):
        try:
            # Get the Enum member that matches the argument value
            return FitToRun(arg_value.lower())
        except ValueError:
            # If the value doesn't match any Enum member, raise an error
            raise argparse.ArgumentTypeError(
                f"Invalid run config: '{arg_value}'. Choose from {list(c.value for c in FitToRun)}"
            )


def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(
        description="Background-only RooFit with a blinded signal region."
    )
    parser.add_argument(
        "--fits-to-run",
        type=FitToRun.arg_type,
        default=FitToRun.ALL,
        help="Choose a run configuration (pseudodata): 1d, 2d, all (default), data or none.",
    )
    parser.add_argument(
        "--events", type=int, default=10000, help="N events to generate (unextended)"
    )
    parser.add_argument(
        "--cheb", type=str, default="-0.2,0.2,-0.1", help="Chebychev coeffs c1,c2,..."
    )
    parser.add_argument("-s", "--seed", type=int)
    parser.add_argument("--nbins", type=int, default=60)
    parser.add_argument(
        "--strict-mode",
        action="store_true",
        default=False,
        help="Abort when a family fails strict selection instead of using the default relaxed fallback.",
    )

    args = parser.parse_args()

    # clear only the workflow-specific plot directories
    (REPO_ROOT / "plots").mkdir(parents=True, exist_ok=True)
    _clear_plot_dir(REPO_ROOT / "plots" / "fit_1d")
    _clear_plot_dir(REPO_ROOT / "plots" / "fit_2d")
    _clear_plot_dir(REPO_ROOT / "plots" / "fit_2d_data")

    gSystem.Load("libRooFit")
    gSystem.Load("libHiggsAnalysisCombinedLimit")

    if args.fits_to_run == FitToRun.ONEDIM or args.fits_to_run == FitToRun.ALL:
        run_fit_1d(args)

    if args.fits_to_run == FitToRun.TWODIM or args.fits_to_run == FitToRun.ALL:
        run_fit_2d(args)


if __name__ == "__main__":
    main()
