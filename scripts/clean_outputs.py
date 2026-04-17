#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


TARGET_FILES = [
    REPO_ROOT / "combine_logger.out",
    REPO_ROOT / "non_resonant_background_workspace.root",
    REPO_ROOT / "resonant_background_fit_HiggsDalitz.root",
    REPO_ROOT / "resonant_background_fit_ZGamma.root",
    REPO_ROOT / "resonant_background_model_Z_params.json",
    REPO_ROOT / "upsilon_model_params.json",
]

TARGET_FILE_GLOBS = [
    "higgsCombine*.root",
    "NormParams_*.json",
    "signal_workspace_*.root",
    "validation.log",
    "validation_workspace.root",
]

TARGET_DIRS = [
    REPO_ROOT / "datacards",
    REPO_ROOT / "plots",
    REPO_ROOT / "__pycache__",
]


def log(message: str) -> None:
    print(f"[clean_outputs] {message}")


def relative_to_repo(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Remove files produced by the local workflows and bundling steps without "
            "touching inputs, source code, or external reference material."
        )
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be removed without deleting anything.",
    )
    return parser


def collect_targets() -> list[Path]:
    targets: list[Path] = []
    for path in TARGET_FILES:
        if path.exists():
            targets.append(path)
    for pattern in TARGET_FILE_GLOBS:
        targets.extend(sorted(REPO_ROOT.glob(pattern)))
    for path in TARGET_DIRS:
        if path.exists():
            targets.append(path)
    return sorted(set(targets))


def remove_path(path: Path, dry_run: bool) -> None:
    if dry_run:
        log(f"would remove {relative_to_repo(path)}")
        return
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()
    log(f"removed {relative_to_repo(path)}")


def main() -> int:
    args = build_parser().parse_args()
    targets = collect_targets()

    if not targets:
        log("no produced files found")
        return 0

    log("cleanup targets")
    for path in targets:
        log(f"  - {relative_to_repo(path)}")

    for path in targets:
        remove_path(path, dry_run=args.dry_run)

    if args.dry_run:
        log("dry run complete")
    else:
        log("cleanup complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
