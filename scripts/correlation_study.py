import argparse
import os
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def main():
    from statistical_test_fit.correlation_study import (
        DEFAULT_M_MUMU_LOWER,
        DEFAULT_M_MUMU_UPPER,
        PLOT_DIR,
        RESULTS_JSON,
        SUMMARY_PLOT,
        SIGNAL_PROCESSES,
        make_correlation_summary_plot,
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
    parser.add_argument(
        "--m-mumu-hist-bins",
        type=int,
        default=None,
        help="Bins for the m_mumu background histogram in the summary plot.",
    )
    parser.add_argument("--m-mumu-min", type=float, default=DEFAULT_M_MUMU_LOWER)
    parser.add_argument("--m-mumu-max", type=float, default=DEFAULT_M_MUMU_UPPER)
    parser.add_argument(
        "--processes",
        nargs="+",
        choices=SIGNAL_PROCESSES,
        default=list(SIGNAL_PROCESSES),
        help="Signal boson processes to fit. Defaults to H and Z.",
    )
    parser.add_argument("--plot-dir", default=PLOT_DIR)
    parser.add_argument("--output-json", default=RESULTS_JSON)
    parser.add_argument("--summary-plot", default=SUMMARY_PLOT)
    parser.add_argument(
        "--exclude-summary-windows",
        nargs="*",
        type=int,
        default=(),
        help="Window indices to omit from the summary plot.",
    )
    parser.add_argument(
        "--only-summary-plot",
        default=False,
        action="store_true",
        help="Build only the summary plot from --output-json without rerunning fits.",
    )
    args = parser.parse_args()

    if args.only_summary_plot:
        make_correlation_summary_plot(
            args.output_json,
            args.summary_plot,
            tuple(args.processes),
            tuple(args.exclude_summary_windows),
        )
        return

    from statistical_test_fit.root_runtime import configure_root

    os.makedirs(args.plot_dir, exist_ok=True)
    for path in Path(args.plot_dir).iterdir():
        if path.is_file() and path.name != ".gitkeep":
            path.unlink()
    configure_root()

    from statistical_test_fit import run_correlation_study

    run_correlation_study(args)


if __name__ == "__main__":
    main()
