import argparse
import os
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def main():
    from statistical_test_fit.mass_ranges import UPSILON_MASS_LOWER, UPSILON_MASS_UPPER
    from statistical_test_fit.root_runtime import configure_root
    from statistical_test_fit.correlation_study import (
        PLOT_DIR,
        RESULTS_JSON,
        SIGNAL_PROCESSES,
    )

    parser = argparse.ArgumentParser(
        description=(
            "Fit m_mumugamma signal shapes in m_mumu windows after merging "
            "the three Upsilon states per boson process."
        )
    )
    parser.add_argument("--nbins", type=int, default=60)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--windows", type=int, default=10)
    parser.add_argument("--m-mumu-min", type=float, default=UPSILON_MASS_LOWER)
    parser.add_argument("--m-mumu-max", type=float, default=UPSILON_MASS_UPPER)
    parser.add_argument(
        "--processes",
        nargs="+",
        choices=SIGNAL_PROCESSES,
        default=list(SIGNAL_PROCESSES),
        help="Signal boson processes to fit. Defaults to H and Z.",
    )
    parser.add_argument("--plot-dir", default=PLOT_DIR)
    parser.add_argument("--output-json", default=RESULTS_JSON)
    args = parser.parse_args()

    os.makedirs(args.plot_dir, exist_ok=True)
    for path in Path(args.plot_dir).iterdir():
        if path.is_file() and path.name != ".gitkeep":
            path.unlink()

    configure_root()

    from statistical_test_fit import run_correlation_study

    run_correlation_study(args)


if __name__ == "__main__":
    main()
