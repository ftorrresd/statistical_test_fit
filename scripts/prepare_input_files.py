#!/usr/bin/env python3
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import sys
import time
import traceback


REQUIRED_WEIGHT_COLUMNS = [
    "weight_l1_prefiring",
    "weight_muon_id",
    "weight_muon_iso",
    "weight_photon_id",
    "weight_photon_electron_veto",
    "weight_pileup",
    "weight_trigger_sf",
    "weight_pdf_alpha_s_weight",
    "weight_generator",
]

MC_NAME_MARKERS = (
    "GluGluHTo",
    "ggH_HTo",
    "HToUps",
    "HToMuMuG",
    "ZGTo",
    "ZToUpsilon",
)

PRODUCT_EXPRESSION = " * ".join(
    f"double({column})" for column in REQUIRED_WEIGHT_COLUMNS
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a workflow-ready input directory where MC ROOT files have "
            "their weight branch redefined as the product of the nominal scale factors."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        default="inputs",
        type=Path,
        help="Directory containing the original input files.",
    )
    parser.add_argument(
        "--output-dir",
        default="inputs_weighted",
        type=Path,
        help="Directory where prepared copies are written.",
    )
    parser.add_argument(
        "--tree-name",
        default="Events",
        help="TTree name to rewrite in MC ROOT files.",
    )
    parser.add_argument(
        "--file-glob",
        default="*",
        help="Input file glob, relative to --input-dir. Useful for focused checks.",
    )
    parser.add_argument(
        "--n-workers",
        type=int,
        default=None,
        help="Number of parallel file-preparation workers. Default: min(number of jobs, nproc).",
    )
    parser.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Fail if an output file already exists. By default outputs are replaced.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned actions without writing files.",
    )
    parser.add_argument(
        "--no-copy-non-mc",
        action="store_true",
        help="Only write rewritten MC files instead of a complete prepared input directory.",
    )
    return parser.parse_args()


def is_mc_file(path: Path) -> bool:
    if path.suffix != ".root":
        return False
    name = path.name
    if name.startswith(("mass_H_", "mass_Z_")):
        return True
    return any(marker in name for marker in MC_NAME_MARKERS)


def is_weight_column(column: str) -> bool:
    return column == "weight" or column.startswith("weight_")


def open_root_file(path: Path):
    import ROOT  # type: ignore

    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {path}")
    return root_file


def branch_names(path: Path, tree_name: str) -> list[str]:
    root_file = open_root_file(path)
    try:
        tree = root_file.Get(tree_name)
        if not tree:
            raise RuntimeError(f"Missing tree {tree_name!r} in {path}")
        branches = tree.GetListOfBranches()
        return [branches.At(index).GetName() for index in range(branches.GetEntries())]
    finally:
        root_file.Close()


def has_tree(path: Path, tree_name: str) -> bool:
    root_file = open_root_file(path)
    try:
        return bool(root_file.Get(tree_name))
    finally:
        root_file.Close()


def validate_mc_file(path: Path, tree_name: str) -> list[str]:
    columns = branch_names(path, tree_name)
    missing = [column for column in REQUIRED_WEIGHT_COLUMNS if column not in columns]
    if "weight" not in columns:
        missing.append("weight")
    return missing


def string_vector(values: list[str]):
    import ROOT  # type: ignore

    vector = ROOT.std.vector("string")()
    for value in values:
        vector.push_back(value)
    return vector


def copy_other_root_keys(
    input_path: Path, output_path: Path, rewritten_tree_name: str
) -> None:
    import ROOT  # type: ignore

    source = ROOT.TFile.Open(str(input_path), "READ")
    target = ROOT.TFile.Open(str(output_path), "UPDATE")
    if not source or source.IsZombie():
        raise RuntimeError(f"Could not reopen source ROOT file: {input_path}")
    if not target or target.IsZombie():
        raise RuntimeError(f"Could not reopen target ROOT file: {output_path}")

    try:
        target.cd()
        keys = source.GetListOfKeys()
        for index in range(keys.GetEntries()):
            key = keys.At(index)
            obj = key.ReadObj()
            if obj.InheritsFrom(ROOT.TTree.Class()) and key.GetName() == rewritten_tree_name:
                continue
            target.cd()
            obj.Write(key.GetName(), ROOT.TObject.kOverwrite)
    finally:
        target.Close()
        source.Close()


def rewrite_mc_file(input_path: Path, output_path: Path, tree_name: str) -> None:
    import ROOT  # type: ignore

    columns = branch_names(input_path, tree_name)
    missing = [column for column in REQUIRED_WEIGHT_COLUMNS if column not in columns]
    if missing:
        raise RuntimeError(
            f"Cannot rewrite {input_path}: missing required columns {', '.join(missing)}"
        )
    if "weight" not in columns:
        raise RuntimeError(f"Cannot rewrite {input_path}: missing existing weight branch")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    options = ROOT.RDF.RSnapshotOptions()
    options.fMode = "RECREATE"
    dataframe = ROOT.RDataFrame(tree_name, str(input_path))
    rewritten = dataframe.Redefine("weight", PRODUCT_EXPRESSION)
    rewritten.Snapshot(tree_name, str(output_path), string_vector(columns), options)
    copy_other_root_keys(input_path, output_path, tree_name)


def strip_data_weight_columns(input_path: Path, output_path: Path, tree_name: str) -> None:
    import ROOT  # type: ignore

    columns = branch_names(input_path, tree_name)
    snapshot_columns = [column for column in columns if not is_weight_column(column)]
    removed_columns = [column for column in columns if is_weight_column(column)]
    if not removed_columns:
        copy_file(input_path, output_path)
        return
    if not snapshot_columns:
        raise RuntimeError(f"Cannot strip all columns from {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    options = ROOT.RDF.RSnapshotOptions()
    options.fMode = "RECREATE"
    dataframe = ROOT.RDataFrame(tree_name, str(input_path))
    dataframe.Snapshot(tree_name, str(output_path), string_vector(snapshot_columns), options)
    copy_other_root_keys(input_path, output_path, tree_name)


def copy_file(input_path: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_path, output_path)


def ensure_can_write(output_path: Path, overwrite: bool) -> None:
    if not output_path.exists():
        return
    if not overwrite:
        raise FileExistsError(
            f"Refusing to overwrite existing file {output_path}. "
            "Remove --no-overwrite to replace generated outputs."
        )
    output_path.unlink()


def resolve_worker_count(requested_workers: int | None, n_jobs: int) -> int:
    if n_jobs <= 0:
        return 0
    cpu_count = os.cpu_count() or 1
    if requested_workers is None:
        return max(1, min(cpu_count, n_jobs))
    return max(1, min(int(requested_workers), n_jobs))


def format_elapsed(seconds: float) -> str:
    return time.strftime("%H:%M:%S", time.gmtime(seconds))


def render_progress(
    completed: int,
    total: int,
    *,
    last_label: str,
    last_status: str,
    started_at: float,
) -> str:
    bar_width = 28
    fraction = 1.0 if total == 0 else completed / total
    filled = int(round(bar_width * fraction))
    bar = "#" * filled + "-" * (bar_width - filled)
    percent = 100.0 * fraction
    elapsed = format_elapsed(time.monotonic() - started_at)
    return (
        f"Prepare inputs [{bar}] {completed}/{total} {percent:5.1f}% "
        f"last={last_label} {last_status} elapsed={elapsed}"
    )


def action_label(action: str, input_path: Path) -> str:
    return f"{action}:{input_path.name}"


def run_action(payload: tuple[str, Path, Path, str]) -> str:
    action, input_path, output_path, tree_name = payload
    if action == "rewrite-mc-weight":
        rewrite_mc_file(input_path, output_path, tree_name)
    elif action == "strip-data-weights":
        strip_data_weight_columns(input_path, output_path, tree_name)
    elif action == "copy-unchanged":
        copy_file(input_path, output_path)
    else:
        raise ValueError(f"Unknown preparation action: {action}")
    return action


def run_actions(
    actions: list[tuple[str, Path, Path]], tree_name: str, requested_workers: int | None
) -> dict[str, int]:
    max_workers = resolve_worker_count(requested_workers, len(actions))
    print(f"\nUsing {max_workers} worker(s) for input preparation.")

    counts = {
        "rewrite-mc-weight": 0,
        "strip-data-weights": 0,
        "copy-unchanged": 0,
    }
    failures: list[tuple[str, str]] = []
    started_at = time.monotonic()

    if max_workers == 1:
        for completed, (action, input_path, output_path) in enumerate(actions, start=1):
            label = action_label(action, input_path)
            status = "OK"
            try:
                counts[run_action((action, input_path, output_path, tree_name))] += 1
            except Exception:
                status = "FAILED"
                failures.append((label, traceback.format_exc()))
            print(
                "\r"
                + render_progress(
                    completed,
                    len(actions),
                    last_label=label,
                    last_status=status,
                    started_at=started_at,
                ),
                end="",
                flush=True,
            )
        print()
    else:
        context = mp.get_context("spawn")
        with ProcessPoolExecutor(
            max_workers=max_workers, mp_context=context
        ) as executor:
            future_to_job = {
                executor.submit(run_action, (action, input_path, output_path, tree_name)): (
                    action,
                    input_path,
                    output_path,
                )
                for action, input_path, output_path in actions
            }
            for completed, future in enumerate(as_completed(future_to_job), start=1):
                action, input_path, _ = future_to_job[future]
                label = action_label(action, input_path)
                status = "OK"
                try:
                    counts[future.result()] += 1
                except Exception:
                    status = "FAILED"
                    failures.append((label, traceback.format_exc()))
                print(
                    "\r"
                    + render_progress(
                        completed,
                        len(actions),
                        last_label=label,
                        last_status=status,
                        started_at=started_at,
                    ),
                    end="",
                    flush=True,
                )
        print()

    if failures:
        print("\nInput preparation failures:")
        for label, failure in failures:
            print(f"--- {label} ---")
            print(failure)
        raise RuntimeError(f"Input preparation failed for {len(failures)} job(s).")

    return counts


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir

    if input_dir.resolve() == output_dir.resolve():
        raise ValueError("--output-dir must be different from --input-dir")
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    files = sorted(path for path in input_dir.glob(args.file_glob) if path.is_file())
    if not files:
        raise FileNotFoundError(f"No files matched {input_dir / args.file_glob}")

    actions: list[tuple[str, Path, Path]] = []
    validation_errors: list[str] = []
    for input_path in files:
        output_path = output_dir / input_path.relative_to(input_dir)
        if is_mc_file(input_path):
            missing = validate_mc_file(input_path, args.tree_name)
            if missing:
                validation_errors.append(
                    f"{input_path}: missing required columns {', '.join(missing)}"
                )
            actions.append(("rewrite-mc-weight", input_path, output_path))
        elif (
            input_path.suffix == ".root"
            and has_tree(input_path, args.tree_name)
            and any(
                is_weight_column(column)
                for column in branch_names(input_path, args.tree_name)
            )
            and not args.no_copy_non_mc
        ):
            actions.append(("strip-data-weights", input_path, output_path))
        elif not args.no_copy_non_mc:
            actions.append(("copy-unchanged", input_path, output_path))

    if validation_errors:
        raise RuntimeError("\n".join(validation_errors))

    if not actions:
        raise RuntimeError("No files selected for preparation")

    print("MC weight expression:")
    print(f"  weight = {PRODUCT_EXPRESSION}")
    print("\nPlanned actions:")
    for action, input_path, output_path in actions:
        print(f"  {action}: {input_path} -> {output_path}")

    if args.dry_run:
        print("\nDry run only; no files written.")
        return 0

    overwrite = not args.no_overwrite
    for _, _, output_path in actions:
        ensure_can_write(output_path, overwrite)

    counts = run_actions(actions, args.tree_name, args.n_workers)

    print("\nPrepared input directory complete.")
    print(f"  Rewritten MC ROOT files: {counts['rewrite-mc-weight']}")
    print(
        "  Data ROOT files with weight columns removed: "
        f"{counts['strip-data-weights']}"
    )
    print(f"  Copied unchanged files: {counts['copy-unchanged']}")
    print(f"  Output directory: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
