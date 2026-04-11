import argparse
import os


def main():
    from statistical_test_fit.root_runtime import configure_root

    # Set up the argument parser
    parser = argparse.ArgumentParser(
        description="Non-resonant background RooFit with a blinded signal region on real data."
    )
    parser.add_argument("--nbins", type=int, default=60)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--use-cache", default=False, action="store_true")
    parser.add_argument(
        "--strict-mode",
        action="store_true",
        default=False,
        help="Abort when a family fails strict selection instead of using the default relaxed fallback.",
    )

    args = parser.parse_args()

    # clear plots dir
    os.system("mkdir -p plots")
    os.system("mkdir -p plots/fit_1d")
    os.system("mkdir -p plots/fit_2d")
    os.system("mkdir -p plots/fit_2d_data")
    os.system(r'find plots -type f ! -name ".gitkeep" -delete')

    configure_root()

    from statistical_test_fit import run_fit_2d_data

    run_fit_2d_data(args)


if __name__ == "__main__":
    main()
