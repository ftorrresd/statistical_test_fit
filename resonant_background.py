import argparse
import os


def main():
    from statistical_test_fit.root_runtime import configure_root

    parser = argparse.ArgumentParser(
        description="Standalone resonant-background modeling workflow."
    )
    parser.add_argument("--nbins", type=int, default=60)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--use-cache", default=False, action="store_true")
    args = parser.parse_args()

    os.system("mkdir -p plots")
    os.system("mkdir -p plots/resonant_background")
    os.system(r'find plots/resonant_background -type f ! -name ".gitkeep" -delete')

    configure_root()

    from statistical_test_fit import run_resonant_background

    run_resonant_background(args)


if __name__ == "__main__":
    main()
