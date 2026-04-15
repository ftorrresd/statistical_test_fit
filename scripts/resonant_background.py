import argparse
import os
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def main():
    from statistical_test_fit.root_runtime import configure_root

    parser = argparse.ArgumentParser(
        description="Standalone resonant-background modeling workflow."
    )
    parser.add_argument("--nbins", type=int, default=60)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument(
        "--skip-cache",
        default=False,
        action="store_true",
        help="Recompute the resonant cached artifacts instead of reusing them.",
    )
    args = parser.parse_args()
    args.use_cache = not args.skip_cache

    os.system("mkdir -p plots")
    os.system("mkdir -p plots/resonant_background")
    os.system(r'find plots/resonant_background -type f ! -name ".gitkeep" -delete')

    configure_root()

    from statistical_test_fit import run_resonant_background

    run_resonant_background(args)


if __name__ == "__main__":
    main()
