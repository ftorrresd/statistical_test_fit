import argparse
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _clear_plot_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for child in path.rglob("*"):
        if child.is_file() and child.name != ".gitkeep":
            child.unlink()


def main():
    from statistical_test_fit.root_runtime import configure_root

    # Set up the argument parser
    parser = argparse.ArgumentParser(
        description="Non-resonant background RooFit with a blinded signal region on real data."
    )
    parser.add_argument("--nbins", type=int, default=60, help="Number of bins for diagnostic plots.")
    parser.add_argument(
        "--chi2-nbins", type=int, default=60, help="Number of bins for chi-square goodness-of-fit computation."
    )
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--use-cache", default=False, action="store_true")
    parser.add_argument(
        "--strict-mode",
        action="store_true",
        default=False,
        help="Abort when a family fails strict selection instead of using the default relaxed fallback.",
    )
    parser.add_argument(
        "--chebychev-orders",
        type=int,
        nargs=2,
        default=[3, 9],
        metavar=("START", "STOP"),
        help="Inclusive order range for Chebychev polynomial candidates (default: 3 9).",
    )
    parser.add_argument(
        "--bernstein-orders",
        type=int,
        nargs=2,
        default=[3, 9],
        metavar=("START", "STOP"),
        help="Inclusive order range for Bernstein polynomial candidates (default: 3 9).",
    )

    args = parser.parse_args()

    # clear only the workflow-specific plot directories
    (REPO_ROOT / "plots").mkdir(parents=True, exist_ok=True)
    _clear_plot_dir(REPO_ROOT / "plots" / "fit_1d")
    _clear_plot_dir(REPO_ROOT / "plots" / "fit_2d")
    _clear_plot_dir(REPO_ROOT / "plots" / "fit_2d_data")

    configure_root()

    from statistical_test_fit import run_fit_2d_data

    run_fit_2d_data(args)


if __name__ == "__main__":
    main()
