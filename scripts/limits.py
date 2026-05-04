#!/usr/bin/env python3

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import html
import json
import math
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


DEFAULT_DATACARD = REPO_ROOT / "datacards" / "datacard.txt"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "datacards" / "blind_limits"
DEFAULT_MASS = "125"
DEFAULT_QUANTILES = (0.16, 0.5, 0.84)
DEFAULT_POI_MIN = 0.0
DEFAULT_POI_MAX = 1_000_000.0
DEFAULT_POI_INITIAL = 1.0
QUICK_HYBRID_TOYS = 100
QUICK_CLS_ACC = 0.02
QUICK_R_REL_ACC = 0.10
QUICK_R_ABS_ACC = 10.0
WORKSPACE_OUTPUT_NAME = "multisignal_workspace.root"
LOCAL_WORKSPACE_NAME = "workspace.root"
LOCAL_DATACARD_NAME = "datacard.txt"


@dataclass(frozen=True)
class DatacardProcesses:
    signals: list[str]
    backgrounds: list[str]


@dataclass(frozen=True)
class PoiTarget:
    label: str
    poi: str
    processes: tuple[str, ...]


@dataclass(frozen=True)
class PoiScheme:
    name: str
    title: str
    description: str
    poi_maps: tuple[str, ...]
    targets: tuple[PoiTarget, ...]
    all_poi_names: tuple[str, ...] = ()

    @property
    def all_pois(self) -> tuple[str, ...]:
        if self.all_poi_names:
            return self.all_poi_names
        seen: set[str] = set()
        pois: list[str] = []
        for target in self.targets:
            if target.poi not in seen:
                pois.append(target.poi)
                seen.add(target.poi)
        return tuple(pois)


@dataclass(frozen=True)
class CommandJob:
    job_id: str
    kind: str
    method: str
    cwd: Path
    command: tuple[str, ...]
    output_patterns: tuple[str, ...]
    metadata: dict[str, Any] = field(default_factory=dict)
    inputs: tuple[dict[str, str], ...] = ()


def log(message: str) -> None:
    print(f"[limits] {message}")


def reset_directory(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def href_relative_to(base_dir: Path, target: Path) -> str:
    try:
        return str(target.resolve().relative_to(base_dir.resolve()))
    except ValueError:
        return str(target)


def ensure_file(path: Path, message: str) -> None:
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"{message}: {path}")


def format_number(value: float) -> str:
    if math.isfinite(value) and value.is_integer():
        return str(int(value))
    return f"{value:.12g}"


def parse_float_list(value: str) -> tuple[float, ...]:
    quantiles: list[float] = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        quantiles.append(float(item))
    if not quantiles:
        raise argparse.ArgumentTypeError("at least one quantile is required")
    return tuple(quantiles)


def quantile_key(value: float) -> str:
    return f"{value:.6g}"


def quantile_tag(value: float) -> str:
    return "q" + quantile_key(value).replace("-", "m").replace(".", "p")


def safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    value = value.strip("._")
    return value or "unnamed"


def is_int_token(value: str) -> bool:
    try:
        int(value)
    except ValueError:
        return False
    return True


def parse_datacard_processes(datacard_path: Path) -> DatacardProcesses:
    lines = datacard_path.read_text(encoding="ascii").splitlines()
    process_lines: list[tuple[int, list[str]]] = []
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        tokens = stripped.split()
        if tokens and tokens[0] == "process" and len(tokens) > 1:
            process_lines.append((line_number, tokens[1:]))

    for index in range(len(process_lines) - 1):
        _, names = process_lines[index]
        _, ids = process_lines[index + 1]
        if len(names) != len(ids):
            continue
        if all(is_int_token(item) for item in names):
            continue
        if not all(is_int_token(item) for item in ids):
            continue

        signals: list[str] = []
        backgrounds: list[str] = []
        for process_name, process_id_text in zip(names, ids):
            process_id = int(process_id_text)
            if process_id <= 0:
                signals.append(process_name)
            else:
                backgrounds.append(process_name)
        if not signals:
            raise RuntimeError(f"No signal processes found in {datacard_path}")
        return DatacardProcesses(signals=signals, backgrounds=backgrounds)

    raise RuntimeError(f"Could not identify process names and ids in {datacard_path}")


def process_state(process_name: str) -> str:
    parts = process_name.split("_", 1)
    if len(parts) != 2 or not parts[1]:
        raise ValueError(
            f"Cannot build grouped Upsilon POI for process {process_name!r}; expected names like H_1S or Z_1S."
        )
    return parts[1]


def process_boson(process_name: str) -> str:
    parts = process_name.split("_", 1)
    if len(parts) != 2 or not parts[0]:
        raise ValueError(
            f"Cannot identify boson for process {process_name!r}; expected names like H_1S or Z_1S."
        )
    return parts[0]


def process_specific_poi(process_name: str) -> str:
    return f"r_{process_name}"


def make_six_poi_scheme(
    signals: list[str],
    poi_initial: float,
    poi_min: float,
    poi_max: float,
) -> PoiScheme:
    maps: list[str] = []
    targets: list[PoiTarget] = []
    for process_name in signals:
        poi = process_specific_poi(process_name)
        maps.append(
            f"map=.*/{process_name}:{poi}[{format_number(poi_initial)},{format_number(poi_min)},{format_number(poi_max)}]"
        )
        targets.append(PoiTarget(label=process_name, poi=poi, processes=(process_name,)))
    return PoiScheme(
        name="six_poi",
        title="Six Independent Signal POIs",
        description=(
            "Each H/Z x Upsilon signal process has its own signal-strength POI. "
            "Limits are computed one POI at a time while the other signal POIs remain in the model."
        ),
        poi_maps=tuple(maps),
        targets=tuple(targets),
        all_poi_names=tuple(target.poi for target in targets),
    )


def make_boson_grouped_poi_scheme(
    signals: list[str],
    target_boson: str,
    poi_initial: float,
    poi_min: float,
    poi_max: float,
) -> PoiScheme:
    target_processes = [
        process_name for process_name in signals if process_boson(process_name) == target_boson
    ]
    profiled_processes = [
        process_name for process_name in signals if process_boson(process_name) != target_boson
    ]
    if not target_processes:
        raise RuntimeError(f"No {target_boson} signal processes found for grouped scheme")

    maps: list[str] = []
    all_pois: list[str] = []
    target_poi = f"r_{target_boson}_grouped"
    all_pois.append(target_poi)
    for index, process_name in enumerate(target_processes):
        if index == 0:
            maps.append(
                f"map=.*/{process_name}:{target_poi}[{format_number(poi_initial)},{format_number(poi_min)},{format_number(poi_max)}]"
            )
        else:
            maps.append(f"map=.*/{process_name}:{target_poi}")

    for process_name in profiled_processes:
        poi = process_specific_poi(process_name)
        all_pois.append(poi)
        maps.append(
            f"map=.*/{process_name}:{poi}[{format_number(poi_initial)},{format_number(poi_min)},{format_number(poi_max)}]"
        )

    profiled_bosons = sorted({process_boson(process_name) for process_name in profiled_processes})
    profiled_text = ", ".join(profiled_bosons) if profiled_bosons else "other"
    scheme_name = f"{target_boson.lower()}_grouped"

    return PoiScheme(
        name=scheme_name,
        title=f"{target_boson} Grouped Signal POI",
        description=(
            f"All {target_boson} signal processes share {target_poi}. "
            f"{profiled_text} signal strengths are mapped to individual profiled POIs."
        ),
        poi_maps=tuple(maps),
        targets=(
            PoiTarget(
                label=f"{target_boson}_grouped",
                poi=target_poi,
                processes=tuple(target_processes),
            ),
        ),
        all_poi_names=tuple(all_pois),
    )


def selected_schemes(args: argparse.Namespace, signals: list[str]) -> tuple[PoiScheme, ...]:
    schemes = {
        "six": make_six_poi_scheme(signals, args.poi_initial, args.poi_min, args.poi_max),
        "z_grouped": make_boson_grouped_poi_scheme(
            signals,
            "Z",
            args.poi_initial,
            args.poi_min,
            args.poi_max,
        ),
        "h_grouped": make_boson_grouped_poi_scheme(
            signals,
            "H",
            args.poi_initial,
            args.poi_min,
            args.poi_max,
        ),
    }
    if args.poi_scheme == "grouped":
        return (schemes["z_grouped"], schemes["h_grouped"])
    if args.poi_scheme == "both":
        return (schemes["six"], schemes["z_grouped"], schemes["h_grouped"])
    return (schemes[args.poi_scheme],)


def copy_file(src: Path, dst: Path) -> dict[str, str]:
    ensure_file(src, "Required input file does not exist")
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)
    return {
        "source": str(src.resolve()),
        "source_repo_relative": repo_relative(src),
        "staged": str(dst.resolve()),
        "staged_name": dst.name,
    }


def stage_datacard_inputs(datacard_path: Path, destination: Path) -> tuple[Path, list[dict[str, str]]]:
    destination.mkdir(parents=True, exist_ok=True)
    input_records: list[dict[str, str]] = []
    staged_shape_files: dict[str, Path] = {}
    staged_lines: list[str] = []

    for line in datacard_path.read_text(encoding="ascii").splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            tokens = stripped.split()
            if len(tokens) >= 5 and tokens[0] == "shapes":
                shape_token = tokens[3]
                if shape_token not in {"FAKE", "-"}:
                    shape_source = Path(shape_token)
                    if not shape_source.is_absolute():
                        shape_source = datacard_path.parent / shape_source
                    shape_source = shape_source.resolve()
                    staged_name = Path(shape_token).name
                    previous_source = staged_shape_files.get(staged_name)
                    if previous_source is not None and previous_source != shape_source:
                        raise RuntimeError(
                            f"Two datacard shape inputs share basename {staged_name!r}: {previous_source} and {shape_source}"
                        )
                    if previous_source is None:
                        staged_shape_files[staged_name] = shape_source
                        input_records.append(copy_file(shape_source, destination / staged_name))
                    tokens[3] = staged_name
                    line = " ".join(tokens)
        staged_lines.append(line)

    staged_datacard = destination / LOCAL_DATACARD_NAME
    staged_datacard.write_text("\n".join(staged_lines) + "\n", encoding="ascii")
    input_records.insert(
        0,
        {
            "source": str(datacard_path.resolve()),
            "source_repo_relative": repo_relative(datacard_path),
            "staged": str(staged_datacard.resolve()),
            "staged_name": staged_datacard.name,
        },
    )
    return staged_datacard, input_records


def set_parameters_argument(pois: tuple[str, ...], value: float = 0.0) -> str:
    return ",".join(f"{poi}={format_number(value)}" for poi in pois)


def parameter_ranges_argument(
    pois: tuple[str, ...],
    lower: float,
    upper: float,
) -> str:
    return ":".join(
        f"{poi}={format_number(lower)},{format_number(upper)}" for poi in pois
    )


def build_workspace_job(
    run_dir: Path,
    datacard_path: Path,
    mass: str,
    scheme: PoiScheme,
) -> CommandJob:
    cwd = run_dir / scheme.name / "workspace_build"
    reset_directory(cwd)
    _, inputs = stage_datacard_inputs(datacard_path, cwd)
    command: list[str] = [
        "text2workspace.py",
        LOCAL_DATACARD_NAME,
        "-m",
        mass,
        "-o",
        WORKSPACE_OUTPUT_NAME,
        "-P",
        "HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel",
    ]
    for poi_map in scheme.poi_maps:
        command.extend(["--PO", poi_map])
    return CommandJob(
        job_id=f"workspace_{scheme.name}",
        kind="workspace_build",
        method="text2workspace.py",
        cwd=cwd,
        command=tuple(command),
        output_patterns=(WORKSPACE_OUTPUT_NAME,),
        inputs=tuple(inputs),
        metadata={
            "scheme": scheme.name,
            "scheme_title": scheme.title,
            "scheme_description": scheme.description,
            "poi_maps": list(scheme.poi_maps),
            "pois": list(scheme.all_pois),
            "mass": mass,
        },
    )


def stage_common_combine_inputs(
    job_dir: Path,
    workspace_path: Path,
    datacard_path: Path,
) -> list[dict[str, str]]:
    records = [
        copy_file(workspace_path, job_dir / LOCAL_WORKSPACE_NAME),
        copy_file(datacard_path, job_dir / LOCAL_DATACARD_NAME),
    ]
    return records


def build_asymptotic_job(
    run_dir: Path,
    scheme: PoiScheme,
    target: PoiTarget,
    mass: str,
    workspace_path: Path,
    datacard_path: Path,
    poi_min: float,
    poi_max: float,
) -> CommandJob:
    cwd = run_dir / scheme.name / "combine" / "asymptotic" / safe_name(target.label)
    reset_directory(cwd)
    inputs = stage_common_combine_inputs(cwd, workspace_path, datacard_path)
    command = [
        "combine",
        "-M",
        "AsymptoticLimits",
        LOCAL_WORKSPACE_NAME,
        "--run",
        "blind",
        "--redefineSignalPOIs",
        target.poi,
        "--setParameters",
        set_parameters_argument(scheme.all_pois, 0.0),
        "--setParameterRanges",
        parameter_ranges_argument(scheme.all_pois, poi_min, poi_max),
        "-m",
        mass,
        "-n",
        f".asymptotic_blind.{safe_name(target.label)}",
    ]
    return CommandJob(
        job_id=f"asymptotic_{scheme.name}_{safe_name(target.label)}",
        kind="combine",
        method="AsymptoticLimits",
        cwd=cwd,
        command=tuple(command),
        output_patterns=("higgsCombine*.root",),
        inputs=tuple(inputs),
        metadata={
            "scheme": scheme.name,
            "target_label": target.label,
            "target_poi": target.poi,
            "target_processes": list(target.processes),
            "profiled_signal_pois": [poi for poi in scheme.all_pois if poi != target.poi],
            "blindness": [
                "AsymptoticLimits is run with --run blind.",
                "Combine therefore uses the pre-fit model state for the expected Asimov calculation instead of fitting observed data.",
                "Observed limits are not requested or extracted from this run.",
            ],
            "mass": mass,
        },
    )


def build_hybrid_job(
    run_dir: Path,
    scheme: PoiScheme,
    target: PoiTarget,
    quantile: float,
    mass: str,
    workspace_path: Path,
    datacard_path: Path,
    poi_min: float,
    poi_max: float,
    hybrid_toys: int | None,
    cls_acc: float | None,
    r_rel_acc: float | None,
    r_abs_acc: float | None,
    save_hybrid_result: bool,
) -> CommandJob:
    cwd = (
        run_dir
        / scheme.name
        / "combine"
        / "hybrid_lhc"
        / safe_name(target.label)
        / quantile_tag(quantile)
    )
    reset_directory(cwd)
    inputs = stage_common_combine_inputs(cwd, workspace_path, datacard_path)
    command = [
        "combine",
        "-M",
        "HybridNew",
        LOCAL_WORKSPACE_NAME,
        "--LHCmode",
        "LHC-limits",
        "--expectedFromGrid",
        quantile_key(quantile),
        "--redefineSignalPOIs",
        target.poi,
        "--setParameters",
        set_parameters_argument(scheme.all_pois, 0.0),
        "--setParameterRanges",
        parameter_ranges_argument(scheme.all_pois, poi_min, poi_max),
        "-m",
        mass,
        "-n",
        f".hybrid_blind.{safe_name(target.label)}.{quantile_tag(quantile)}",
    ]
    if save_hybrid_result:
        command.append("--saveHybridResult")
    if hybrid_toys is not None:
        command.extend(["-T", str(hybrid_toys)])
    if cls_acc is not None:
        command.extend(["--clsAcc", format_number(cls_acc)])
    if r_rel_acc is not None:
        command.extend(["--rRelAcc", format_number(r_rel_acc)])
    if r_abs_acc is not None:
        command.extend(["--rAbsAcc", format_number(r_abs_acc)])

    return CommandJob(
        job_id=f"hybrid_lhc_{scheme.name}_{safe_name(target.label)}_{quantile_tag(quantile)}",
        kind="combine",
        method="HybridNew",
        cwd=cwd,
        command=tuple(command),
        output_patterns=("higgsCombine*.root",),
        inputs=tuple(inputs),
        metadata={
            "scheme": scheme.name,
            "target_label": target.label,
            "target_poi": target.poi,
            "target_processes": list(target.processes),
            "profiled_signal_pois": [poi for poi in scheme.all_pois if poi != target.poi],
            "quantile": quantile,
            "hybrid_toys": hybrid_toys,
            "cls_acc": cls_acc,
            "r_rel_acc": r_rel_acc,
            "r_abs_acc": r_abs_acc,
            "blindness": [
                "HybridNew uses --LHCmode LHC-limits for LHC-style CLs limits.",
                "HybridNew is run directly with --expectedFromGrid for the requested quantile.",
                "No --dataset or --bypassFrequentistFit option is passed.",
            ],
            "mass": mass,
        },
    )


def format_duration(total_seconds: float) -> str:
    total_milliseconds = max(0, int(round(total_seconds * 1000.0)))
    total_seconds_int, milliseconds = divmod(total_milliseconds, 1000)
    total_minutes, seconds = divmod(total_seconds_int, 60)
    total_hours, minutes = divmod(total_minutes, 60)
    days, hours = divmod(total_hours, 24)
    return f"{days:02d}:{hours:02d}:{minutes:02d}:{seconds:02d}:{milliseconds:03d}"


def run_command_job(job: CommandJob) -> dict[str, Any]:
    job.cwd.mkdir(parents=True, exist_ok=True)
    started_at = dt.datetime.now(dt.timezone.utc)
    start_time = time.monotonic()
    stdout = ""
    stderr = ""
    returncode = 0

    try:
        completed = subprocess.run(
            list(job.command),
            cwd=job.cwd,
            capture_output=True,
            text=True,
        )
        stdout = completed.stdout or ""
        stderr = completed.stderr or ""
        returncode = int(completed.returncode)
    except FileNotFoundError as exc:
        returncode = 127
        stderr = str(exc)
    except Exception as exc:  # pragma: no cover - defensive runtime guard
        returncode = 1
        stderr = f"Unhandled exception while running command: {exc}"

    ended_at = dt.datetime.now(dt.timezone.utc)
    duration = time.monotonic() - start_time
    stdout_path = job.cwd / "stdout.txt"
    stderr_path = job.cwd / "stderr.txt"
    stdout_path.write_text(stdout, encoding="utf-8")
    stderr_path.write_text(stderr, encoding="utf-8")

    produced_files: list[str] = []
    for pattern in job.output_patterns:
        produced_files.extend(path.name for path in sorted(job.cwd.glob(pattern)))
    produced_files = sorted(set(produced_files))

    return {
        "schema_version": 1,
        "job_id": job.job_id,
        "kind": job.kind,
        "method": job.method,
        "cwd": str(job.cwd.resolve()),
        "cwd_repo_relative": repo_relative(job.cwd),
        "command": list(job.command),
        "command_string": shlex.join(job.command),
        "returncode": returncode,
        "status": "ok" if returncode == 0 else "failed",
        "started_at": started_at.isoformat(),
        "ended_at": ended_at.isoformat(),
        "duration_seconds": duration,
        "duration": format_duration(duration),
        "stdout": stdout,
        "stderr": stderr,
        "stdout_path": str(stdout_path.resolve()),
        "stderr_path": str(stderr_path.resolve()),
        "stdout_path_repo_relative": repo_relative(stdout_path),
        "stderr_path_repo_relative": repo_relative(stderr_path),
        "produced_root_files": produced_files,
        "inputs": list(job.inputs),
        "metadata": job.metadata,
    }


def job_manifest_label(job: CommandJob) -> str:
    details: list[str] = []
    metadata = job.metadata
    if metadata.get("scheme"):
        details.append(f"scheme={metadata['scheme']}")
    if metadata.get("target_poi"):
        details.append(f"poi={metadata['target_poi']}")
    if metadata.get("quantile") is not None:
        details.append(f"quantile={metadata['quantile']}")
    if metadata.get("target_processes"):
        details.append("processes=" + ",".join(metadata["target_processes"]))
    return "; ".join(details) if details else "preparation"


def print_job_manifest(wave_name: str, jobs: list[CommandJob], workers: int) -> None:
    if not jobs:
        log(f"No jobs planned for {wave_name}.")
        return
    max_workers = len(jobs) if workers <= 0 else min(workers, len(jobs))
    log(f"Planned jobs for {wave_name}: {len(jobs)} job(s), {max_workers} worker(s)")
    log("Each listed working directory has been cleared before input staging.")
    for index, job in enumerate(jobs, start=1):
        log(f"  {index}. {job.job_id} [{job.method}; {job_manifest_label(job)}]")
        log(f"     cwd: {repo_relative(job.cwd)}")
        log(f"     command: {shlex.join(job.command)}")


def run_jobs_parallel(
    jobs: list[CommandJob],
    workers: int,
    wave_name: str,
) -> list[dict[str, Any]]:
    print_job_manifest(wave_name, jobs, workers)
    if not jobs:
        return []
    max_workers = len(jobs) if workers <= 0 else min(workers, len(jobs))
    log(f"Running {len(jobs)} command(s) with {max_workers} worker(s)")
    results: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_job = {executor.submit(run_command_job, job): job for job in jobs}
        for future in concurrent.futures.as_completed(future_to_job):
            job = future_to_job[future]
            try:
                result = future.result()
            except Exception as exc:  # pragma: no cover - defensive runtime guard
                result = {
                    "schema_version": 1,
                    "job_id": job.job_id,
                    "kind": job.kind,
                    "method": job.method,
                    "cwd": str(job.cwd.resolve()),
                    "cwd_repo_relative": repo_relative(job.cwd),
                    "command": list(job.command),
                    "command_string": shlex.join(job.command),
                    "returncode": 1,
                    "status": "failed",
                    "duration_seconds": 0.0,
                    "duration": format_duration(0.0),
                    "stdout": "",
                    "stderr": f"Unhandled executor exception: {exc}",
                    "produced_root_files": [],
                    "inputs": list(job.inputs),
                    "metadata": job.metadata,
                }
            log(
                f"{result['status']}: {job.job_id} "
                f"in {result.get('duration', format_duration(float(result.get('duration_seconds', 0.0))))}"
            )
            results.append(result)
    return sorted(results, key=lambda item: str(item["job_id"]))


def leaf_value_to_json(value: float, type_name: str) -> int | float:
    integer_types = {
        "Char_t",
        "UChar_t",
        "Short_t",
        "UShort_t",
        "Int_t",
        "UInt_t",
        "Long_t",
        "ULong_t",
        "Long64_t",
        "ULong64_t",
        "Bool_t",
    }
    if type_name in integer_types:
        return int(round(value))
    return float(value)


def extract_root_file(root_path: Path) -> dict[str, Any]:
    try:
        from statistical_test_fit.root_runtime import configure_root

        ROOT = configure_root()
        root_file = ROOT.TFile.Open(str(root_path))
        if root_file is None or root_file.IsZombie():
            return {
                "file": str(root_path.resolve()),
                "file_name": root_path.name,
                "status": "failed",
                "message": "Could not open ROOT file.",
                "has_limit_tree": False,
                "has_toy_asimov": False,
                "limit_tree": None,
            }

        top_level_keys = [key.GetName() for key in list(root_file.GetListOfKeys())]
        toys_directory = root_file.Get("toys")
        toy_keys = (
            [key.GetName() for key in list(toys_directory.GetListOfKeys())]
            if toys_directory is not None and hasattr(toys_directory, "GetListOfKeys")
            else []
        )
        toy_asimov = root_file.Get("toys/toy_asimov")
        tree = root_file.Get("limit")
        if tree is None:
            root_file.Close()
            return {
                "file": str(root_path.resolve()),
                "file_name": root_path.name,
                "status": "ok",
                "top_level_keys": top_level_keys,
                "toy_keys": toy_keys,
                "has_limit_tree": False,
                "has_toy_asimov": toy_asimov is not None,
                "limit_tree": None,
            }

        branches = [branch.GetName() for branch in list(tree.GetListOfBranches())]
        branch_types: dict[str, str] = {}
        entries: list[dict[str, Any]] = []
        for entry_index in range(int(tree.GetEntries())):
            tree.GetEntry(entry_index)
            entry: dict[str, Any] = {"entry": entry_index}
            for branch_name in branches:
                leaf = tree.GetLeaf(branch_name)
                if leaf is None:
                    continue
                type_name = str(leaf.GetTypeName())
                branch_types[branch_name] = type_name
                leaf_len = max(1, int(leaf.GetLen()))
                values = [
                    leaf_value_to_json(float(leaf.GetValue(value_index)), type_name)
                    for value_index in range(leaf_len)
                ]
                entry[branch_name] = values[0] if leaf_len == 1 else values
            entries.append(entry)

        root_file.Close()
        return {
            "file": str(root_path.resolve()),
            "file_name": root_path.name,
            "status": "ok",
            "top_level_keys": top_level_keys,
            "toy_keys": toy_keys,
            "has_limit_tree": True,
            "has_toy_asimov": toy_asimov is not None,
            "limit_tree": {
                "branches": branches,
                "branch_types": branch_types,
                "entries": entries,
            },
        }
    except Exception as exc:
        return {
            "file": str(root_path.resolve()),
            "file_name": root_path.name,
            "status": "failed",
            "message": str(exc),
            "has_limit_tree": False,
            "has_toy_asimov": False,
            "limit_tree": None,
        }


def summarize_limits(root_outputs: list[dict[str, Any]]) -> dict[str, Any]:
    expected: dict[str, list[float]] = {}
    observed: list[float] = []
    generated: list[float] = []
    all_entries: list[dict[str, Any]] = []

    for root_output in root_outputs:
        limit_tree = root_output.get("limit_tree")
        if not limit_tree:
            continue
        for entry in limit_tree.get("entries", []):
            all_entries.append(entry)
            if "limit" not in entry or "quantileExpected" not in entry:
                continue
            limit_value = float(entry["limit"])
            quantile = float(entry["quantileExpected"])
            if math.isclose(quantile, -1.0, rel_tol=0.0, abs_tol=1e-5):
                observed.append(limit_value)
            elif math.isclose(quantile, -2.0, rel_tol=0.0, abs_tol=1e-5):
                generated.append(limit_value)
            elif quantile >= 0.0:
                expected.setdefault(quantile_key(quantile), []).append(limit_value)

    return {
        "expected": expected,
        "observed": observed,
        "generated": generated,
        "entries": all_entries,
    }


def attach_outputs(result: dict[str, Any]) -> dict[str, Any]:
    cwd = Path(result["cwd"])
    root_outputs = [
        extract_root_file(cwd / root_name) for root_name in result.get("produced_root_files", [])
    ]
    result["root_outputs"] = root_outputs
    result["limits"] = summarize_limits(root_outputs)
    result["json_path"] = str((cwd / "result.json").resolve())
    result["json_path_repo_relative"] = repo_relative(cwd / "result.json")
    result["summary_html_path"] = str((cwd / "summary.html").resolve())
    result["summary_html_path_repo_relative"] = repo_relative(cwd / "summary.html")
    (cwd / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (cwd / "summary.html").write_text(make_job_html_summary(result), encoding="utf-8")
    return result


def attach_and_write_outputs(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [attach_outputs(result) for result in results]


def html_escape(value: Any) -> str:
    return html.escape(str(value), quote=True)


def html_document(title: str, body: str) -> str:
    return f"""<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <title>{html_escape(title)}</title>
  <style>
    :root {{
      color-scheme: dark;
      --bg: #070b12;
      --panel: #101826;
      --panel-2: #0c1320;
      --text: #e6edf7;
      --muted: #9aa8bd;
      --line: #243247;
      --accent: #7dd3fc;
      --accent-2: #a78bfa;
      --ok: #34d399;
      --fail: #fb7185;
      --warn: #fbbf24;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: radial-gradient(circle at top left, rgba(125, 211, 252, 0.15), transparent 32rem),
        linear-gradient(180deg, #070b12 0%, #0a101c 100%);
      color: var(--text);
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      line-height: 1.5;
    }}
    main {{ max-width: 1200px; margin: 0 auto; padding: 2.2rem; }}
    h1 {{ margin: 0 0 0.5rem; font-size: clamp(1.9rem, 4vw, 3.1rem); letter-spacing: -0.04em; }}
    h2 {{ margin: 0 0 1rem; color: var(--accent); font-size: 1.05rem; text-transform: uppercase; letter-spacing: 0.12em; }}
    .subtitle {{ color: var(--muted); margin: 0 0 2rem; }}
    .grid {{ display: grid; gap: 1rem; grid-template-columns: repeat(auto-fit, minmax(18rem, 1fr)); }}
    .card {{
      background: linear-gradient(180deg, rgba(16, 24, 38, 0.96), rgba(12, 19, 32, 0.96));
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 1.2rem;
      box-shadow: 0 20px 60px rgba(0, 0, 0, 0.24);
      margin: 1rem 0;
      overflow: auto;
    }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.93rem; }}
    th, td {{ padding: 0.72rem 0.85rem; border-bottom: 1px solid var(--line); vertical-align: top; }}
    th {{ color: #dbeafe; background: rgba(125, 211, 252, 0.08); text-align: left; font-weight: 700; }}
    tr:hover td {{ background: rgba(125, 211, 252, 0.045); }}
    code, pre {{ font-family: "SFMono-Regular", Consolas, "Liberation Mono", monospace; }}
    pre {{
      white-space: pre-wrap;
      word-break: break-word;
      background: #050914;
      color: #dbeafe;
      border: 1px solid var(--line);
      border-radius: 14px;
      padding: 1rem;
      overflow: auto;
    }}
    a {{ color: var(--accent); text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .pill {{ display: inline-block; padding: 0.2rem 0.55rem; border-radius: 999px; font-size: 0.8rem; font-weight: 700; }}
    .ok {{ color: #052e1b; background: var(--ok); }}
    .failed {{ color: #3f0712; background: var(--fail); }}
    .muted {{ color: var(--muted); }}
    ul {{ margin: 0.4rem 0 0; padding-left: 1.2rem; }}
  </style>
</head>
<body>
  <main>
{body}
  </main>
</body>
</html>
"""


def html_table(headers: list[str], rows: list[list[Any]]) -> str:
    if not rows:
        rows = [["-" for _ in headers]]
    header_html = "".join(f"<th>{html_escape(header)}</th>" for header in headers)
    row_html = []
    for row in rows:
        row_html.append("<tr>" + "".join(f"<td>{html_escape(value)}</td>" for value in row) + "</tr>")
    return "<table><thead><tr>" + header_html + "</tr></thead><tbody>" + "".join(row_html) + "</tbody></table>"


def html_link(label: str, href: str) -> str:
    return f'<a href="{html_escape(href)}">{html_escape(label)}</a>'


def html_table_raw(headers: list[str], rows: list[list[Any]]) -> str:
    if not rows:
        rows = [["-" for _ in headers]]
    header_html = "".join(f"<th>{html_escape(header)}</th>" for header in headers)
    row_html = []
    for row in rows:
        row_html.append("<tr>" + "".join(f"<td>{value}</td>" for value in row) + "</tr>")
    return "<table><thead><tr>" + header_html + "</tr></thead><tbody>" + "".join(row_html) + "</tbody></table>"


def status_pill(status: Any) -> str:
    status_text = str(status)
    css = "ok" if status_text == "ok" else "failed" if status_text == "failed" else "muted"
    return f'<span class="pill {css}">{html_escape(status_text)}</span>'


def make_limit_entry_rows(root_outputs: list[dict[str, Any]]) -> list[list[Any]]:
    rows: list[list[Any]] = []
    for root_output in root_outputs:
        limit_tree = root_output.get("limit_tree")
        if not limit_tree:
            continue
        for entry in limit_tree.get("entries", []):
            rows.append(
                [
                    root_output.get("file_name", "-"),
                    entry.get("entry", "-"),
                    entry.get("quantileExpected", "-"),
                    entry.get("limit", "-"),
                    entry.get("limitErr", "-"),
                    entry.get("mh", "-"),
                    entry.get("iToy", "-"),
                    entry.get("iSeed", "-"),
                ]
            )
    return rows


def make_job_html_summary(result: dict[str, Any]) -> str:
    metadata = result.get("metadata", {})
    input_rows = [
        [item.get("staged_name", "-"), item.get("source_repo_relative", item.get("source", "-"))]
        for item in result.get("inputs", [])
    ]
    output_rows = [
        [
            output.get("file_name", "-"),
            "yes" if output.get("has_limit_tree") else "no",
            "yes" if output.get("has_toy_asimov") else "no",
            output.get("status", "-"),
        ]
        for output in result.get("root_outputs", [])
    ]
    blindness_lines = metadata.get("blindness", [])
    if blindness_lines:
        policy_html = "<ul>" + "".join(f"<li>{html_escape(line)}</li>" for line in blindness_lines) + "</ul>"
    else:
        policy_html = "<p class=\"muted\">Not a Combine limit call; this command prepares inputs.</p>"

    status_rows = [
        ["Kind", result.get("kind", "-")],
        ["Method", result.get("method", "-")],
        ["Status", result.get("status", "-")],
        ["Return code", result.get("returncode", "-")],
        ["Duration", result.get("duration", format_duration(float(result.get("duration_seconds", 0.0))))],
        ["Working directory", result.get("cwd_repo_relative", result.get("cwd", "-"))],
    ]
    logs_rows = [
        ["stdout", html_link("stdout.txt", "stdout.txt")],
        ["stderr", html_link("stderr.txt", "stderr.txt")],
        ["json", html_link("result.json", "result.json")],
    ]
    body = f"""
    <h1>{html_escape(result['job_id'])}</h1>
    <p class="subtitle">Per-job Combine execution summary</p>
    <section class="card"><h2>Status</h2>{html_table(['Field', 'Value'], status_rows)}</section>
    <section class="card"><h2>Command</h2><pre>{html_escape(result.get('command_string', ''))}</pre></section>
    <section class="card"><h2>Limit Policy</h2>{policy_html}</section>
    <section class="card"><h2>Inputs Staged In This Directory</h2>{html_table(['Local file', 'Source'], input_rows)}</section>
    <section class="card"><h2>ROOT Outputs</h2>{html_table(['File', 'Limit tree', 'toy_asimov', 'Status'], output_rows)}</section>
    <section class="card"><h2>Limit Tree Entries</h2>{html_table(['File', 'Entry', 'quantileExpected', 'limit', 'limitErr', 'mh', 'iToy', 'iSeed'], make_limit_entry_rows(result.get('root_outputs', [])))}</section>
    <section class="card"><h2>Logs</h2>{html_table_raw(['Stream', 'Path'], logs_rows)}</section>
    """
    return html_document(str(result["job_id"]), body)


def result_successful(result: dict[str, Any]) -> bool:
    return int(result.get("returncode", 1)) == 0


def result_output_path(result: dict[str, Any], expected_name: str | None = None) -> Path | None:
    cwd = Path(result["cwd"])
    root_files = [cwd / name for name in result.get("produced_root_files", [])]
    if expected_name is not None:
        for root_file in root_files:
            if root_file.name == expected_name:
                return root_file
    return root_files[0] if root_files else None


def short_result(result: dict[str, Any]) -> dict[str, Any]:
    metadata = result.get("metadata", {})
    return {
        "job_id": result.get("job_id"),
        "kind": result.get("kind"),
        "method": result.get("method"),
        "status": result.get("status"),
        "returncode": result.get("returncode"),
        "duration_seconds": result.get("duration_seconds"),
        "duration": result.get("duration"),
        "cwd": result.get("cwd_repo_relative", result.get("cwd")),
        "json": result.get("json_path_repo_relative", result.get("json_path")),
        "summary_html": result.get(
            "summary_html_path_repo_relative",
            result.get("summary_html_path"),
        ),
        "scheme": metadata.get("scheme"),
        "target_poi": metadata.get("target_poi"),
        "target_label": metadata.get("target_label"),
        "target_processes": metadata.get("target_processes"),
        "quantile": metadata.get("quantile"),
        "limits": result.get("limits"),
    }


def build_limit_summary(results: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for result in results:
        metadata = result.get("metadata", {})
        scheme = metadata.get("scheme")
        method = result.get("method")
        target_poi = metadata.get("target_poi")
        if not scheme or not target_poi or method not in {"AsymptoticLimits", "HybridNew"}:
            continue
        method_key = "hybrid_lhc" if method == "HybridNew" else "asymptotic"
        scheme_summary = summary.setdefault(scheme, {}).setdefault(method_key, {})
        poi_summary = scheme_summary.setdefault(
            target_poi,
            {
                "target_label": metadata.get("target_label"),
                "target_processes": metadata.get("target_processes", []),
                "jobs": [],
                "expected": {},
                "observed": [],
            },
        )
        poi_summary["jobs"].append(result.get("job_id"))
        limits = result.get("limits", {})
        for quantile, values in limits.get("expected", {}).items():
            poi_summary["expected"].setdefault(quantile, []).extend(values)
        poi_summary["observed"].extend(limits.get("observed", []))
    return summary


def write_run_summary(
    run_dir: Path,
    args: argparse.Namespace,
    processes: DatacardProcesses,
    schemes: tuple[PoiScheme, ...],
    results: list[dict[str, Any]],
) -> None:
    readme_html = run_dir / "README.html"
    payload = {
        "schema_version": 1,
        "run_dir": str(run_dir.resolve()),
        "run_dir_repo_relative": repo_relative(run_dir),
        "readme_html": str(readme_html.resolve()),
        "readme_html_repo_relative": repo_relative(readme_html),
        "datacard": str(Path(args.datacard).resolve()),
        "datacard_repo_relative": repo_relative(Path(args.datacard).resolve()),
        "mass": args.mass,
        "methods": args.methods,
        "quantiles": list(args.quantiles),
        "quick": args.quick,
        "poi_min": args.poi_min,
        "poi_max": args.poi_max,
        "hybrid_toys": args.hybrid_toys,
        "cls_acc": args.cls_acc,
        "r_rel_acc": args.r_rel_acc,
        "r_abs_acc": args.r_abs_acc,
        "signals": processes.signals,
        "backgrounds": processes.backgrounds,
        "schemes": [
            {
                "name": scheme.name,
                "title": scheme.title,
                "description": scheme.description,
                "pois": list(scheme.all_pois),
                "poi_maps": list(scheme.poi_maps),
                "targets": [
                    {
                        "label": target.label,
                        "poi": target.poi,
                        "processes": list(target.processes),
                    }
                    for target in scheme.targets
                ],
            }
            for scheme in schemes
        ],
        "jobs": [short_result(result) for result in results],
        "limits": build_limit_summary(results),
    }
    summary_json = run_dir / "blind_limits_summary.json"
    summary_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    job_rows = [
        [
            html_escape(result.get("job_id")),
            html_escape(result.get("method")),
            status_pill(result.get("status")),
            html_escape(result.get("metadata", {}).get("scheme", "-")),
            html_escape(result.get("metadata", {}).get("target_poi", "-")),
            html_escape(result.get("metadata", {}).get("quantile", "-")),
            html_escape(
                result.get(
                    "duration",
                    format_duration(float(result.get("duration_seconds", 0.0))),
                )
            ),
            html_link(
                "summary.html",
                href_relative_to(run_dir, Path(result.get("summary_html_path", ""))),
            )
            if result.get("summary_html_path")
            else "-",
        ]
        for result in results
    ]
    summary_rows = [
        ["Run directory", repo_relative(run_dir)],
        ["Datacard", repo_relative(Path(args.datacard).resolve())],
        ["Mass", args.mass],
        ["Methods", args.methods],
        ["Quick mode", "yes" if args.quick else "no"],
        ["POI scheme request", args.poi_scheme],
        ["POI range", f"{format_number(args.poi_min)} to {format_number(args.poi_max)}"],
        ["Quantiles", ", ".join(quantile_key(q) for q in args.quantiles)],
        ["Aggregate JSON", html_link("blind_limits_summary.json", "blind_limits_summary.json")],
    ]
    policy_html = """
      <ul>
        <li>AsymptoticLimits jobs use <code>--run blind</code>.</li>
        <li>HybridNew jobs use <code>--expectedFromGrid</code> for the requested expected quantiles.</li>
        <li>HybridNew jobs do not pass <code>--dataset</code> or <code>--bypassFrequentistFit</code>.</li>
      </ul>
    """
    body = f"""
    <h1>Blind Limit Run</h1>
    <p class="subtitle">Aggregate execution summary for Combine limit jobs.</p>
    <section class="card"><h2>Summary</h2>{html_table_raw(['Field', 'Value'], [[html_escape(k), v if str(v).startswith('<a ') else html_escape(v)] for k, v in summary_rows])}</section>
    <section class="card"><h2>Limit Policy</h2>{policy_html}</section>
    <section class="card"><h2>Jobs</h2>{html_table_raw(['Job', 'Method', 'Status', 'Scheme', 'Target POI', 'Quantile', 'Duration', 'Summary'], job_rows)}</section>
    """
    readme_html.write_text(html_document("Blind Limit Run", body), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run expected upper-limit workflows for the bundled simultaneous Combine card. "
            "Independent six-POI and H/Z grouped-POI schemes are supported."
        )
    )
    parser.add_argument(
        "--datacard",
        default=str(DEFAULT_DATACARD),
        help="Input text datacard produced by scripts/build_bundled_workspace.py.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Parent directory for blind-limit runs.",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help="Run subdirectory name. Defaults to a UTC timestamp.",
    )
    parser.add_argument("--mass", default=DEFAULT_MASS, help="Mass label passed to Combine.")
    parser.add_argument(
        "--poi-scheme",
        choices=("six", "z_grouped", "h_grouped", "grouped", "both"),
        default="both",
        help=(
            "POI scheme to run. `grouped` runs z_grouped and h_grouped; "
            "the default `both` runs six plus both grouped schemes."
        ),
    )
    parser.add_argument(
        "--methods",
        choices=("asymptotic", "hybrid", "both"),
        default="both",
        help="Limit method(s) to run.",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help=(
            "Run a fast non-default HybridNew configuration: "
            "--rRelAcc 0.10 --rAbsAcc 10 --clsAcc 0.02 -T 100."
        ),
    )
    parser.add_argument(
        "--quantiles",
        type=parse_float_list,
        default=DEFAULT_QUANTILES,
        help="Comma-separated HybridNew expected quantiles.",
    )
    parser.add_argument(
        "--poi-initial",
        type=float,
        default=DEFAULT_POI_INITIAL,
        help="Initial value used when defining POIs in text2workspace.py.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=0,
        help="Parallel subprocess workers per dependency wave. Use 0 to run all ready jobs at once.",
    )
    parser.add_argument(
        "--hybrid-toys",
        type=int,
        default=None,
        help="Optional -T value for HybridNew.",
    )
    parser.add_argument(
        "--cls-acc",
        type=float,
        default=None,
        help="Optional --clsAcc value for HybridNew.",
    )
    parser.add_argument(
        "--no-save-hybrid-result",
        action="store_true",
        help="Do not pass --saveHybridResult to HybridNew.",
    )
    return parser


def apply_quick_mode(args: argparse.Namespace) -> None:
    args.poi_min = DEFAULT_POI_MIN
    args.poi_max = DEFAULT_POI_MAX
    args.r_rel_acc = None
    args.r_abs_acc = None
    if not args.quick:
        return
    args.hybrid_toys = QUICK_HYBRID_TOYS
    args.cls_acc = QUICK_CLS_ACC
    args.r_rel_acc = QUICK_R_REL_ACC
    args.r_abs_acc = QUICK_R_ABS_ACC


def prepare_run_dir(output_dir: Path, run_name: str | None) -> Path:
    if run_name is None:
        run_name = dt.datetime.now(dt.timezone.utc).strftime("run_%Y%m%d_%H%M%S")
    run_dir = output_dir / run_name
    if run_dir.exists():
        log(f"Clearing run directory before starting: {repo_relative(run_dir)}")
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    apply_quick_mode(args)
    datacard_path = Path(args.datacard).resolve()
    args.datacard = str(datacard_path)
    ensure_file(datacard_path, "Input datacard not found")
    output_dir = Path(args.output_dir).resolve()
    run_dir = prepare_run_dir(output_dir, args.run_name)

    processes = parse_datacard_processes(datacard_path)
    schemes = selected_schemes(args, processes.signals)

    log(f"Run directory: {repo_relative(run_dir)}")
    log(f"Signals: {', '.join(processes.signals)}")
    log(f"Backgrounds: {', '.join(processes.backgrounds)}")
    log(f"POI schemes: {', '.join(scheme.name for scheme in schemes)}")
    if args.quick:
        log(
            "Quick mode enabled: HybridNew "
            "--rRelAcc 0.10 --rAbsAcc 10 --clsAcc 0.02 -T 100"
        )
    log("POI range: 0 to 1000000")

    all_results: list[dict[str, Any]] = []

    workspace_jobs = [
        build_workspace_job(run_dir, datacard_path, args.mass, scheme) for scheme in schemes
    ]
    workspace_results = attach_and_write_outputs(
        run_jobs_parallel(workspace_jobs, args.workers, "workspace builds")
    )
    all_results.extend(workspace_results)
    if not all(result_successful(result) for result in workspace_results):
        write_run_summary(run_dir, args, processes, schemes, all_results)
        return 1

    workspace_paths: dict[str, Path] = {}
    staged_datacards: dict[str, Path] = {}
    for result in workspace_results:
        scheme_name = result["metadata"]["scheme"]
        workspace_path = result_output_path(result, WORKSPACE_OUTPUT_NAME)
        if workspace_path is None:
            raise RuntimeError(f"No {WORKSPACE_OUTPUT_NAME} output found for {scheme_name}")
        workspace_paths[scheme_name] = workspace_path
        staged_datacards[scheme_name] = Path(result["cwd"]) / LOCAL_DATACARD_NAME

    run_hybrid = args.methods in {"hybrid", "both"}
    run_asymptotic = args.methods in {"asymptotic", "both"}

    limit_jobs: list[CommandJob] = []
    for scheme in schemes:
        if run_asymptotic:
            for target in scheme.targets:
                limit_jobs.append(
                    build_asymptotic_job(
                        run_dir,
                        scheme,
                        target,
                        args.mass,
                        workspace_paths[scheme.name],
                        staged_datacards[scheme.name],
                        args.poi_min,
                        args.poi_max,
                    )
                )
        if run_hybrid:
            for target in scheme.targets:
                for quantile in args.quantiles:
                    limit_jobs.append(
                        build_hybrid_job(
                            run_dir,
                            scheme,
                            target,
                            quantile,
                            args.mass,
                            workspace_paths[scheme.name],
                            staged_datacards[scheme.name],
                            args.poi_min,
                            args.poi_max,
                            args.hybrid_toys,
                            args.cls_acc,
                            args.r_rel_acc,
                            args.r_abs_acc,
                            not args.no_save_hybrid_result,
                        )
                    )

    limit_results = attach_and_write_outputs(
        run_jobs_parallel(limit_jobs, args.workers, "Combine limit jobs")
    )
    all_results.extend(limit_results)
    write_run_summary(run_dir, args, processes, schemes, all_results)

    failed = [result for result in all_results if not result_successful(result)]
    if failed:
        log(f"Completed with {len(failed)} failed command(s).")
        return 1

    log(f"Completed successfully. Summary: {repo_relative(run_dir / 'blind_limits_summary.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
