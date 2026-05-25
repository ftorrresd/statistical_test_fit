#!/usr/bin/env python3

from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import html
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.text import Text


REPO_ROOT = Path(__file__).resolve().parents[1]
CONSOLE = Console()
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from statistical_test_fit.display_names import (
    dataset_strategy_display,
    floating_pdf_display,
    infer_pdf_family_from_name,
    infer_pdf_order,
    normalize_pdf_family,
    pdf_state_display,
    poi_scheme_display,
    signal_target_display,
)


DEFAULT_DATACARD = REPO_ROOT / "datacards" / "datacard.txt"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "bias_study"
DEFAULT_MASS = "125"
DEFAULT_INJECTIONS = (0.0, 10.0, 100.0, 1000.0)
DEFAULT_TOYS = 1000
DEFAULT_DATASET_STRATEGY = "toys"
DATASET_STRATEGIES = ("toys", "asimov")
DEFAULT_PDF_TARGET_STRATEGY = "floating"
PDF_TARGET_STRATEGIES = ("fixed", "floating", "both")
DEFAULT_POI_INITIAL = 1.0
DEFAULT_POI_MIN = -1_000_000.0
DEFAULT_POI_MAX = 1_000_000.0
DEFAULT_PULL_RANGE = (-5.0, 5.0)
DEFAULT_PULL_BINS = 40
DEFAULT_MIN_FIT_ENTRIES = 3
DEFAULT_BIAS_PULL_THRESHOLD = 0.2
DEFAULT_PULL_WIDTH_THRESHOLD = 0.2
DEFAULT_SEED_BASE = 123450
DEFAULT_CMIN_DEFAULT_MINIMIZER_STRATEGY = 2
QUICK_CMIN_DEFAULT_MINIMIZER_STRATEGY = 0
SCATTER_PLOT_SUBPROCESS_TIMEOUT_SECONDS = 1800
HEATMAP_PLOT_SUBPROCESS_TIMEOUT_SECONDS = 1800
FIT_ALGO = "singles"
WORKSPACE_OUTPUT_NAME = "multisignal_workspace.root"
LOCAL_WORKSPACE_NAME = "workspace.root"
LOCAL_DATACARD_NAME = "datacard.txt"
LOCAL_TOYS_NAME = "toys.root"
NON_RESONANT_PROCESS_NAME = "non_resonant_bkg"
DEFAULT_PDF_INDEX_NAME = "pdfindex"


class StoreCminStrategy(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):  # type: ignore[override]
        setattr(namespace, self.dest, values)
        setattr(namespace, "cmin_strategy_explicit", True)


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
class TruthPdf:
    index: int
    name: str
    family: str
    label: str
    slug: str
    latex: str
    scan_order: int | None = None
    selection_role: str | None = None


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
    CONSOLE.print(Text("[bias_study]", style="dim cyan"), message)


def nproc() -> int:
    return max(1, os.cpu_count() or 1)


def resolve_worker_count(job_count: int, requested_workers: int) -> int:
    if job_count <= 0:
        return 0
    worker_cap = nproc()
    requested = worker_cap if requested_workers <= 0 else requested_workers
    return max(1, min(job_count, requested, worker_cap))


def worker_policy_text(requested_workers: int, job_count: int) -> str:
    worker_cap = nproc()
    effective = resolve_worker_count(job_count, requested_workers)
    if requested_workers <= 0:
        return f"requested=auto, nproc={worker_cap}, effective={effective}"
    if requested_workers > worker_cap:
        return f"requested={requested_workers}, capped_by_nproc={worker_cap}, effective={effective}"
    return f"requested={requested_workers}, nproc={worker_cap}, effective={effective}"


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
    values: list[float] = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        values.append(float(item))
    if not values:
        raise argparse.ArgumentTypeError("at least one value is required")
    return tuple(values)


def parse_scheme_injection_spec(value: str) -> tuple[str, tuple[float, ...]]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected SCHEME=R0,R1,...")
    scheme, injections = value.split("=", 1)
    scheme = scheme.strip()
    if not scheme:
        raise argparse.ArgumentTypeError("scheme name is required before '='")
    return scheme, parse_float_list(injections)


def parse_family_list(value: str) -> tuple[str, ...]:
    families: list[str] = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        family = normalize_pdf_family(item)
        families.append(family)
    if not families:
        raise argparse.ArgumentTypeError("at least one family is required")
    return tuple(families)


def safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    value = value.strip("._")
    return value or "unnamed"


def scheme_slug(name: str) -> str:
    return poi_scheme_display(name).slug


def scheme_label(name: str) -> str:
    return poi_scheme_display(name).text


def scheme_latex(name: str) -> str:
    return poi_scheme_display(name).latex


def target_slug(value: str) -> str:
    return signal_target_display(value).slug


def target_label(value: str) -> str:
    return signal_target_display(value).text


def target_latex(value: str) -> str:
    return signal_target_display(value).latex


def injection_output_slug(value: float) -> str:
    return "injected_r_" + format_number(value).replace("-", "minus_").replace(".", "p")


def injection_tag(value: float) -> str:
    return "r" + format_number(value).replace("-", "m").replace(".", "p")


def selected_dataset_strategies(args: argparse.Namespace) -> tuple[str, ...]:
    return (str(args.dataset_strategy),)


def normalize_pdf_target_strategy(value: Any) -> str:
    strategy = str(value or "floating").strip().lower()
    if strategy in {"free", "loose"}:
        return "floating"
    return strategy


def selected_pdf_target_strategies(args: argparse.Namespace) -> tuple[str, ...]:
    strategy = normalize_pdf_target_strategy(args.pdf_target_strategy)
    if strategy == "both":
        return ("floating", "fixed")
    return (strategy,)


def strategy_toys_argument(dataset_strategy: str, args: argparse.Namespace) -> int:
    return int(args.toys) if dataset_strategy == "toys" else -1


def strategy_dataset_count(dataset_strategy: str, args: argparse.Namespace) -> int:
    return int(args.toys) if dataset_strategy == "toys" else 1


def effective_poi_range(args: argparse.Namespace) -> tuple[float, float]:
    return float(args.poi_min), float(args.poi_max)


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


def make_three_poi_scheme(
    signals: list[str],
    target_boson: str,
    poi_initial: float,
    poi_min: float,
    poi_max: float,
) -> PoiScheme:
    maps: list[str] = []
    targets: list[PoiTarget] = []
    target_processes = [
        process_name for process_name in signals if process_boson(process_name) == target_boson
    ]
    if not target_processes:
        raise RuntimeError(f"No {target_boson} signal processes found for three-POI scheme")
    for process_name in signals:
        poi = process_specific_poi(process_name)
        maps.append(
            f"map=.*/{process_name}:{poi}[{format_number(poi_initial)},{format_number(poi_min)},{format_number(poi_max)}]"
        )
        if process_name in target_processes:
            targets.append(PoiTarget(label=process_name, poi=poi, processes=(process_name,)))
    scheme_name = f"three_poi_{target_boson.lower()}"
    return PoiScheme(
        name=scheme_name,
        title=f"Three Independent {target_boson} Signal POIs",
        description=(
            f"Each {target_boson} x Upsilon signal process has its own target signal-strength POI. "
            "The opposite-boson signal strengths remain in the model as profiled POIs."
        ),
        poi_maps=tuple(maps),
        targets=tuple(targets),
        all_poi_names=tuple(process_specific_poi(process_name) for process_name in signals),
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
        "three_poi_z": make_three_poi_scheme(
            signals,
            "Z",
            args.poi_initial,
            args.poi_min,
            args.poi_max,
        ),
        "three_poi_h": make_three_poi_scheme(
            signals,
            "H",
            args.poi_initial,
            args.poi_min,
            args.poi_max,
        ),
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
    if args.poi_scheme == "three_poi":
        return (schemes["three_poi_z"], schemes["three_poi_h"])
    if args.poi_scheme == "grouped":
        return (schemes["z_grouped"], schemes["h_grouped"])
    if args.poi_scheme == "both":
        return (schemes["three_poi_z"], schemes["three_poi_h"], schemes["z_grouped"], schemes["h_grouped"])
    return (schemes[args.poi_scheme],)


def injections_by_scheme(args: argparse.Namespace, schemes: tuple[PoiScheme, ...]) -> dict[str, tuple[float, ...]]:
    explicit = dict(getattr(args, "scheme_injections", []) or [])
    selected_names = {scheme.name for scheme in schemes}
    unknown = sorted(name for name in explicit if name not in selected_names)
    if unknown:
        raise ValueError(
            "--scheme-injections specified unselected scheme(s): " + ", ".join(unknown)
        )
    return {scheme.name: tuple(explicit.get(scheme.name, args.injections)) for scheme in schemes}


def format_injection_values(values: tuple[float, ...]) -> str:
    return ",".join(format_number(float(value)) for value in values)


def format_injections_by_scheme(mapping: dict[str, tuple[float, ...]]) -> str:
    return "; ".join(
        f"{scheme_label(scheme_name)}={format_injection_values(values)}"
        for scheme_name, values in mapping.items()
    )


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


def stage_common_combine_inputs(
    job_dir: Path,
    workspace_path: Path,
    datacard_path: Path,
) -> list[dict[str, str]]:
    return [
        copy_file(workspace_path, job_dir / LOCAL_WORKSPACE_NAME),
        copy_file(datacard_path, job_dir / LOCAL_DATACARD_NAME),
    ]


def set_parameters_argument(pois: tuple[str, ...], target_poi: str, value: float, mode: str) -> str:
    assignments: list[str] = []
    for poi in pois:
        assigned_value = value if (mode == "all-pois" or poi == target_poi) else 0.0
        assignments.append(f"{poi}={format_number(assigned_value)}")
    return ",".join(assignments)


def append_set_parameter(assignments: str, name: str, value: float | int) -> str:
    suffix = f"{name}={format_number(float(value))}"
    return f"{assignments},{suffix}" if assignments else suffix


def parameter_ranges_argument(pois: tuple[str, ...], lower: float, upper: float) -> str:
    return ":".join(f"{poi}={format_number(lower)},{format_number(upper)}" for poi in pois)


def build_workspace_job(run_dir: Path, datacard_path: Path, mass: str, scheme: PoiScheme) -> CommandJob:
    scheme_display = poi_scheme_display(scheme.name)
    cwd = run_dir / scheme_display.slug / "workspace_build"
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
        job_id=f"workspace_{scheme_display.slug}",
        kind="workspace_build",
        method="text2workspace.py",
        cwd=cwd,
        command=tuple(command),
        output_patterns=(WORKSPACE_OUTPUT_NAME,),
        inputs=tuple(inputs),
        metadata={
            "scheme": scheme.name,
            "scheme_slug": scheme_display.slug,
            "scheme_label": scheme_display.text,
            "scheme_latex": scheme_display.latex,
            "scheme_title": scheme.title,
            "scheme_description": scheme.description,
            "poi_maps": list(scheme.poi_maps),
            "pois": list(scheme.all_pois),
            "mass": mass,
        },
    )


def build_dataset_generation_job(
    run_dir: Path,
    scheme: PoiScheme,
    target: PoiTarget,
    truth_pdf: TruthPdf,
    injected_r: float,
    mass: str,
    workspace_path: Path,
    datacard_path: Path,
    args: argparse.Namespace,
    dataset_strategy: str,
    experiment_index: int,
) -> CommandJob:
    generation_dir_name = "toys" if dataset_strategy == "toys" else "asimov"
    scheme_display = poi_scheme_display(scheme.name)
    target_display = signal_target_display(target.label)
    dataset_display = dataset_strategy_display(dataset_strategy)
    truth_slug = truth_pdf.slug or safe_name(truth_pdf.label)
    cwd = (
        run_dir
        / scheme_display.slug
        / target_display.slug
        / truth_slug
        / injection_output_slug(injected_r)
        / generation_dir_name
    )
    reset_directory(cwd)
    inputs = stage_common_combine_inputs(cwd, workspace_path, datacard_path)
    toys_argument = strategy_toys_argument(dataset_strategy, args)
    poi_min, poi_max = effective_poi_range(args)
    set_parameters = set_parameters_argument(
        scheme.all_pois,
        target.poi,
        injected_r,
        args.injection_mode,
    )
    set_parameters = append_set_parameter(set_parameters, args.pdf_index_name, truth_pdf.index)
    freeze_parameters = [args.pdf_index_name]
    if args.injection_mode == "target-only":
        freeze_parameters.extend(poi for poi in scheme.all_pois if poi != target.poi)
    command = [
        "combine",
        "-M",
        "GenerateOnly",
        LOCAL_WORKSPACE_NAME,
        "-t",
        str(toys_argument),
        "--saveToys",
        "--toysFrequentist",
        "--bypassFrequentistFit",
        "--setParameters",
        set_parameters,
        "--freezeParameters",
        ",".join(freeze_parameters),
        "--setParameterRanges",
        parameter_ranges_argument(scheme.all_pois, poi_min, poi_max),
        "-m",
        mass,
        "-n",
        f".{dataset_display.slug}.{scheme_display.slug}.{target_display.slug}.{truth_slug}.{injection_output_slug(injected_r)}",
    ]
    seed = None
    if dataset_strategy == "toys":
        if args.random_seeds:
            seed = -1
            command.extend(["-s", str(seed)])
        elif args.seed_base is not None:
            seed = int(args.seed_base) + experiment_index
            command.extend(["-s", str(seed)])
    job_prefix = "toys" if dataset_strategy == "toys" else "asimov"
    return CommandJob(
        job_id=(
            f"{job_prefix}_{scheme_display.slug}_{target_display.slug}_"
            f"{truth_slug}_{injection_output_slug(injected_r)}"
        ),
        kind="toy_generation" if dataset_strategy == "toys" else "asimov_generation",
        method="GenerateOnly",
        cwd=cwd,
        command=tuple(command),
        output_patterns=("higgsCombine*.root",),
        inputs=tuple(inputs),
        metadata={
            "scheme": scheme.name,
            "scheme_slug": scheme_display.slug,
            "scheme_label": scheme_display.text,
            "scheme_latex": scheme_display.latex,
            "target_label": target.label,
            "target_slug": target_display.slug,
            "target_display_label": target_display.text,
            "target_latex": target_display.latex,
            "target_poi": target.poi,
            "target_processes": list(target.processes),
            "truth_pdf_index": truth_pdf.index,
            "truth_pdf": truth_pdf.name,
            "truth_pdf_label": truth_pdf.label,
            "truth_pdf_slug": truth_pdf.slug,
            "truth_pdf_latex": truth_pdf.latex,
            "truth_pdf_family": truth_pdf.family,
            "truth_pdf_scan_order": truth_pdf.scan_order,
            "truth_pdf_selection_role": truth_pdf.selection_role,
            "injected_r": injected_r,
            "injection_mode": args.injection_mode,
            "dataset_strategy": dataset_strategy,
            "frozen_parameters": freeze_parameters,
            "toys": toys_argument,
            "pseudo_datasets": strategy_dataset_count(dataset_strategy, args),
            "seed": seed,
            "requested_poi_min": args.poi_min,
            "requested_poi_max": args.poi_max,
            "poi_range_policy": "fixed requested range",
            "poi_min": poi_min,
            "poi_max": poi_max,
            "pdf_index_name": args.pdf_index_name,
            "mass": mass,
        },
    )


def build_fit_job(
    run_dir: Path,
    scheme: PoiScheme,
    target: PoiTarget,
    truth_pdf: TruthPdf,
    target_pdf: TruthPdf | None,
    injected_r: float,
    mass: str,
    workspace_path: Path,
    datacard_path: Path,
    toys_path: Path,
    args: argparse.Namespace,
    dataset_strategy: str,
) -> CommandJob:
    fit_dir_name = "fit" if dataset_strategy == "toys" else "asimov_fit"
    fit_pdf_mode = "floating" if target_pdf is None else "fixed"
    scheme_display = poi_scheme_display(scheme.name)
    target_display = signal_target_display(target.label)
    dataset_display = dataset_strategy_display(dataset_strategy)
    truth_slug = truth_pdf.slug or safe_name(truth_pdf.label)
    fit_pdf_display = floating_pdf_display() if target_pdf is None else pdf_state_display(
        index=target_pdf.index,
        family=target_pdf.family,
        order=target_pdf.scan_order,
        selection_role=target_pdf.selection_role,
        name=target_pdf.name,
    )
    fit_pdf_tag = fit_pdf_display.slug if target_pdf is None else f"fixed_{fit_pdf_display.slug}"
    cwd = (
        run_dir
        / scheme_display.slug
        / target_display.slug
        / truth_slug
        / injection_output_slug(injected_r)
        / fit_dir_name
        / fit_pdf_tag
    )
    reset_directory(cwd)
    inputs = stage_common_combine_inputs(cwd, workspace_path, datacard_path)
    inputs.append(copy_file(toys_path, cwd / LOCAL_TOYS_NAME))
    toys_argument = strategy_toys_argument(dataset_strategy, args)
    poi_min, poi_max = effective_poi_range(args)
    scanned_pois = (target.poi,)
    profiled_pois = tuple(poi for poi in scheme.all_pois if poi not in scanned_pois)
    set_parameters = set_parameters_argument(
        scheme.all_pois,
        target.poi,
        injected_r,
        args.injection_mode,
    )
    freeze_parameters: list[str] = []
    if target_pdf is not None:
        set_parameters = append_set_parameter(set_parameters, args.pdf_index_name, target_pdf.index)
        freeze_parameters.append(args.pdf_index_name)
    command = [
        "combine",
        "-M",
        "MultiDimFit",
        LOCAL_WORKSPACE_NAME,
        "--algo",
        FIT_ALGO,
        "-P",
        target.poi,
        "--floatOtherPOIs",
        "1",
        "--toysFile",
        LOCAL_TOYS_NAME,
        "-t",
        str(toys_argument),
        "--toysFrequentist",
        "--bypassFrequentistFit",
        "--redefineSignalPOIs",
        ",".join(scheme.all_pois),
        "--setParameters",
        set_parameters,
        "--setParameterRanges",
        parameter_ranges_argument(scheme.all_pois, poi_min, poi_max),
        "--cminDefaultMinimizerStrategy",
        str(args.cmin_default_minimizer_strategy),
        "--trackParameters",
        ",".join(scheme.all_pois),
        "--trackErrors",
        ",".join(scheme.all_pois),
        "--saveSpecifiedIndex",
        args.pdf_index_name,
        "-m",
        mass,
        "-n",
        f".{fit_dir_name}.{fit_pdf_tag}.{scheme_display.slug}.{target_display.slug}.{truth_slug}.{injection_output_slug(injected_r)}",
    ]
    if freeze_parameters:
        command.extend(["--freezeParameters", ",".join(freeze_parameters)])
    if args.profile_freeze_disassociated_params:
        command.extend(["--X-rtd", "MINIMIZER_freezeDisassociatedParams"])
    if args.robust_fit:
        command.extend(["--robustFit", "1"])

    job_prefix = "fit" if dataset_strategy == "toys" else "asimov_fit"
    return CommandJob(
        job_id=(
            f"{job_prefix}_{fit_pdf_tag}_{scheme_display.slug}_{target_display.slug}_"
            f"{truth_slug}_{injection_output_slug(injected_r)}"
        ),
        kind="toy_fit" if dataset_strategy == "toys" else "asimov_fit",
        method="MultiDimFit",
        cwd=cwd,
        command=tuple(command),
        output_patterns=("higgsCombine*.root",),
        inputs=tuple(inputs),
        metadata={
            "scheme": scheme.name,
            "scheme_slug": scheme_display.slug,
            "scheme_label": scheme_display.text,
            "scheme_latex": scheme_display.latex,
            "target_label": target.label,
            "target_slug": target_display.slug,
            "target_display_label": target_display.text,
            "target_latex": target_display.latex,
            "target_poi": target.poi,
            "all_pois": list(scheme.all_pois),
            "target_processes": list(target.processes),
            "truth_pdf_index": truth_pdf.index,
            "truth_pdf": truth_pdf.name,
            "truth_pdf_label": truth_pdf.label,
            "truth_pdf_slug": truth_pdf.slug,
            "truth_pdf_latex": truth_pdf.latex,
            "truth_pdf_family": truth_pdf.family,
            "truth_pdf_scan_order": truth_pdf.scan_order,
            "truth_pdf_selection_role": truth_pdf.selection_role,
            "fit_method": "MultiDimFit",
            "fit_algo": FIT_ALGO,
            "fit_pdf_mode": fit_pdf_mode,
            "fit_pdf_index": None if target_pdf is None else target_pdf.index,
            "fit_pdf": None if target_pdf is None else target_pdf.name,
            "fit_pdf_label": fit_pdf_display.text,
            "fit_pdf_slug": fit_pdf_display.slug,
            "fit_pdf_latex": fit_pdf_display.latex,
            "fit_pdf_family": None if target_pdf is None else target_pdf.family,
            "target_pdf_index": None if target_pdf is None else target_pdf.index,
            "target_pdf": None if target_pdf is None else target_pdf.name,
            "target_pdf_label": fit_pdf_display.text,
            "target_pdf_slug": fit_pdf_display.slug,
            "target_pdf_latex": fit_pdf_display.latex,
            "target_pdf_family": None if target_pdf is None else target_pdf.family,
            "free_floating_pois": list(scheme.all_pois),
            "scanned_pois": list(scanned_pois),
            "profiled_pois": list(profiled_pois),
            "poi_scan_mode": "target-only",
            "float_other_pois": True,
            "free_floating_pdf_indexes": [args.pdf_index_name] if target_pdf is None else [],
            "frozen_parameters": freeze_parameters,
            "injected_r": injected_r,
            "injection_mode": args.injection_mode,
            "dataset_strategy": dataset_strategy,
            "toys": toys_argument,
            "pseudo_datasets": strategy_dataset_count(dataset_strategy, args),
            "requested_poi_min": args.poi_min,
            "requested_poi_max": args.poi_max,
            "poi_range_policy": "fixed requested range",
            "poi_min": poi_min,
            "poi_max": poi_max,
            "quick": bool(args.quick),
            "cmin_default_minimizer_strategy": int(args.cmin_default_minimizer_strategy),
            "cmin_strategy_explicit": bool(args.cmin_strategy_explicit),
            "robust_fit": bool(args.robust_fit),
            "pdf_index_name": args.pdf_index_name,
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


def format_command_txt(command: list[str]) -> str:
    separator = " " + "\\" + "\n  "
    return separator.join(shlex.quote(part) for part in command)


def run_command_job(job: CommandJob) -> dict[str, Any]:
    job.cwd.mkdir(parents=True, exist_ok=True)
    command_string = shlex.join(job.command)
    command_path = job.cwd / "command.txt"
    command_path.write_text(
        format_command_txt(list(job.command)) + "\n",
        encoding="utf-8",
    )
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

    result = {
        "schema_version": 1,
        "job_id": job.job_id,
        "kind": job.kind,
        "method": job.method,
        "cwd": str(job.cwd.resolve()),
        "cwd_repo_relative": repo_relative(job.cwd),
        "command": list(job.command),
        "command_string": command_string,
        "command_path": str(command_path.resolve()),
        "command_path_repo_relative": repo_relative(command_path),
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
    write_job_summary(result)
    return result


def executor_exception_result(job: CommandJob, exc: Exception) -> dict[str, Any]:
    job.cwd.mkdir(parents=True, exist_ok=True)
    command_string = shlex.join(job.command)
    command_path = job.cwd / "command.txt"
    stdout_path = job.cwd / "stdout.txt"
    stderr_path = job.cwd / "stderr.txt"
    command_path.write_text(
        format_command_txt(list(job.command)) + "\n",
        encoding="utf-8",
    )
    stdout_path.write_text("", encoding="utf-8")
    stderr_text = f"Unhandled executor exception: {exc}"
    stderr_path.write_text(stderr_text, encoding="utf-8")
    now = dt.datetime.now(dt.timezone.utc).isoformat()
    result = {
        "schema_version": 1,
        "job_id": job.job_id,
        "kind": job.kind,
        "method": job.method,
        "cwd": str(job.cwd.resolve()),
        "cwd_repo_relative": repo_relative(job.cwd),
        "command": list(job.command),
        "command_string": command_string,
        "command_path": str(command_path.resolve()),
        "command_path_repo_relative": repo_relative(command_path),
        "returncode": 1,
        "status": "failed",
        "started_at": now,
        "ended_at": now,
        "duration_seconds": 0.0,
        "duration": format_duration(0.0),
        "stdout": "",
        "stderr": stderr_text,
        "stdout_path": str(stdout_path.resolve()),
        "stderr_path": str(stderr_path.resolve()),
        "stdout_path_repo_relative": repo_relative(stdout_path),
        "stderr_path_repo_relative": repo_relative(stderr_path),
        "produced_root_files": [],
        "inputs": list(job.inputs),
        "metadata": job.metadata,
    }
    write_job_summary(result)
    return result


def job_manifest_label(job: CommandJob) -> str:
    metadata = job.metadata
    details: list[str] = []
    for key in (
        "dataset_strategy",
        "scheme",
        "target_poi",
        "truth_pdf_label",
        "fit_pdf_mode",
        "fit_pdf_label",
        "injected_r",
    ):
        if metadata.get(key) is not None:
            details.append(f"{key}={metadata[key]}")
    return "; ".join(details) if details else "preparation"


def print_job_manifest(wave_name: str, jobs: list[CommandJob], workers: int) -> None:
    if not jobs:
        log(f"No jobs planned for {wave_name}.")
        return
    max_workers = resolve_worker_count(len(jobs), workers)
    log(f"Planned jobs for {wave_name}: {len(jobs)} job(s), {max_workers} worker(s)")
    log(f"Worker policy: {worker_policy_text(workers, len(jobs))}")
    log("Each listed working directory has been cleared before input staging.")
    for index, job in enumerate(jobs, start=1):
        log(f"  {index}. {job.job_id} [{job.method}; {job_manifest_label(job)}]")
        log(f"     cwd: {repo_relative(job.cwd)}")
        log(f"     command: {shlex.join(job.command)}")


def run_jobs_parallel(jobs: list[CommandJob], workers: int, wave_name: str) -> list[dict[str, Any]]:
    print_job_manifest(wave_name, jobs, workers)
    if not jobs:
        return []
    max_workers = resolve_worker_count(len(jobs), workers)
    log(f"Running {len(jobs)} command(s) with {max_workers} worker(s)")
    results: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_job = {executor.submit(run_command_job, job): job for job in jobs}
        pending_job_ids = {job.job_id for job in jobs}
        for future in concurrent.futures.as_completed(future_to_job):
            job = future_to_job[future]
            try:
                result = future.result()
            except Exception as exc:  # pragma: no cover - defensive runtime guard
                result = executor_exception_result(job, exc)
            results.append(result)
            pending_job_ids.discard(job.job_id)
            completed = len(results)
            remaining = len(jobs) - completed
            status_text = str(result["status"])
            status_style = "bold green" if status_text == "ok" else "bold red"
            CONSOLE.print(
                Text("[bias_study]", style="dim cyan"),
                Text(status_text, style=status_style),
                f": {job.job_id} "
                f"in {result.get('duration', format_duration(float(result.get('duration_seconds', 0.0))))}"
                f" [{completed}/{len(jobs)} done, {remaining} remaining]",
            )
            if remaining < 10:
                remaining_jobs = ", ".join(sorted(pending_job_ids)) or "none"
                log(f"remaining jobs ({remaining}): {remaining_jobs}")
    return sorted(results, key=lambda item: str(item["job_id"]))


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


def first_root_file(result: dict[str, Any], prefix: str) -> Path | None:
    cwd = Path(result["cwd"])
    for root_name in result.get("produced_root_files", []):
        if root_name.startswith(prefix):
            return cwd / root_name
    return None


def open_workspace_from_file(root_path: Path, preferred_names: tuple[str, ...] = ("w", "combined_workspace")) -> tuple[Any, Any, Any]:
    from statistical_test_fit.root_runtime import configure_root

    ROOT = configure_root()
    root_file = ROOT.TFile.Open(str(root_path))
    if root_file is None or root_file.IsZombie():
        raise RuntimeError(f"Could not open ROOT file {root_path}")

    for workspace_name in preferred_names:
        workspace = root_file.Get(workspace_name)
        if workspace is not None and workspace.InheritsFrom("RooWorkspace"):
            return ROOT, root_file, workspace

    for key in list(root_file.GetListOfKeys()):
        obj = root_file.Get(key.GetName())
        if obj is not None and obj.InheritsFrom("RooWorkspace"):
            return ROOT, root_file, obj

    root_file.Close()
    raise RuntimeError(f"No RooWorkspace found in {root_path}")


def infer_pdf_family(name: str) -> str:
    return infer_pdf_family_from_name(name)


def load_bundle_pdf_state_metadata(datacard_path: Path) -> dict[int, dict[str, Any]]:
    summary_path = datacard_path.parent / "bundle_summary.json"
    if not summary_path.exists():
        return {}
    try:
        payload = json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception as exc:
        log(f"Could not read bundle PDF metadata from {repo_relative(summary_path)}: {exc}")
        return {}
    states = payload.get("non_resonant_state_mappings", [])
    mapping: dict[int, dict[str, Any]] = {}
    for state in states:
        if not isinstance(state, dict) or state.get("index") is None:
            continue
        try:
            mapping[int(state["index"])] = dict(state)
        except (TypeError, ValueError):
            continue
    return mapping


def make_truth_pdf(index: int, pdf_name: str, state_metadata: dict[str, Any] | None = None) -> TruthPdf:
    state_metadata = state_metadata or {}
    family = normalize_pdf_family(
        state_metadata.get("pdf_family")
        or state_metadata.get("family")
        or infer_pdf_family_from_name(pdf_name, state_metadata.get("source_model_name"))
    )
    order = infer_pdf_order(
        state_metadata.get("scan_order"),
        pdf_name,
        state_metadata.get("source_model_name"),
        state_metadata.get("readable_label"),
    )
    selection_role = state_metadata.get("selection_role")
    display = pdf_state_display(
        index=index,
        family=family,
        order=order,
        selection_role=selection_role,
        name=pdf_name,
        source_model_name=state_metadata.get("source_model_name"),
        workspace_label=state_metadata.get("workspace_label"),
    )
    return TruthPdf(
        index=index,
        name=pdf_name,
        family=family,
        label=str(state_metadata.get("display_label") or display.text),
        slug=str(state_metadata.get("display_slug") or display.slug),
        latex=str(state_metadata.get("display_latex") or display.latex),
        scan_order=order,
        selection_role=None if selection_role is None else str(selection_role),
    )


def discover_truth_pdfs(
    workspace_path: Path,
    requested_families: tuple[str, ...],
    pdf_state_metadata: dict[int, dict[str, Any]] | None = None,
) -> tuple[list[TruthPdf], dict[str, Any]]:
    pdf_state_metadata = pdf_state_metadata or {}
    ROOT, root_file, workspace = open_workspace_from_file(workspace_path)
    try:
        multipdf_candidates = []
        for pdf in list(workspace.allPdfs()):
            if hasattr(pdf, "getNumPdfs") and hasattr(pdf, "getPdf"):
                try:
                    n_pdfs = int(pdf.getNumPdfs())
                except Exception:
                    continue
                if n_pdfs <= 0:
                    continue
                name = str(pdf.GetName())
                score = 0
                lower = name.lower()
                if NON_RESONANT_PROCESS_NAME in name:
                    score += 10
                if "non_resonant" in lower or "nonres" in lower:
                    score += 5
                multipdf_candidates.append((score, name, pdf, n_pdfs))
        if not multipdf_candidates:
            raise RuntimeError(
                f"Could not find a RooMultiPdf-like object in {workspace_path}. "
                "The bias study needs the non-resonant RooMultiPdf states."
            )
        multipdf_candidates.sort(key=lambda item: (-item[0], item[1]))
        score, multipdf_name, multipdf, n_pdfs = multipdf_candidates[0]
        available_pdfs: list[TruthPdf] = []
        truths: list[TruthPdf] = []
        for index in range(n_pdfs):
            pdf = multipdf.getPdf(index)
            pdf_name = str(pdf.GetName()) if pdf is not None else f"pdf_{index}"
            candidate = make_truth_pdf(index, pdf_name, pdf_state_metadata.get(index))
            available_pdfs.append(candidate)
            if requested_families == ("all",) or candidate.family in requested_families:
                truths.append(candidate)

        category_names = [cat.GetName() for cat in list(workspace.allCats())]
        metadata = {
            "workspace": workspace.GetName(),
            "multipdf": multipdf_name,
            "multipdf_score": score,
            "multipdf_states_total": n_pdfs,
            "categories": category_names,
            "available_target_pdfs": [target.__dict__ for target in available_pdfs],
            "selected_truth_pdfs": [truth.__dict__ for truth in truths],
        }
        if not truths:
            raise RuntimeError(
                f"No RooMultiPdf states matched requested families {requested_families}. "
                f"Available states are under {multipdf_name} with {n_pdfs} entries."
            )
        return truths, metadata
    finally:
        root_file.Close()


def leaf_value(tree: Any, name: str) -> float | None:
    leaf = tree.GetLeaf(name)
    try:
        if leaf is None or not leaf:
            return None
        return float(leaf.GetValue())
    except ReferenceError:
        return None


def extract_pull_values(
    fit_root_path: Path,
    poi: str,
    all_pois: list[str],
    injected_r: float,
    pdf_index_name: str,
    fit_algo: str,
    scanned_pois: list[str] | None = None,
) -> dict[str, Any]:
    from statistical_test_fit.root_runtime import configure_root

    ROOT = configure_root()
    root_file = ROOT.TFile.Open(str(fit_root_path))
    if root_file is None or root_file.IsZombie():
        return {
            "status": "failed",
            "message": f"Could not open {fit_root_path}",
            "pulls": [],
        }
    try:
        tree = root_file.Get("limit")
        if tree is None:
            return {
                "status": "failed",
                "message": "limit tree not found",
                "pulls": [],
            }
        branches = [branch.GetName() for branch in list(tree.GetListOfBranches())]
        required_branches = [poi, "quantileExpected"]
        missing = [branch for branch in required_branches if branch not in branches]
        if missing:
            return {
                "status": "failed",
                "message": "Missing required branches: " + ", ".join(missing),
                "branches": branches,
                "pulls": [],
            }

        all_pois = list(dict.fromkeys(str(poi_name) for poi_name in all_pois))
        if poi not in all_pois:
            all_pois = [poi]
        if scanned_pois is None:
            scanned_pois = list(all_pois)
        else:
            scanned_pois = list(dict.fromkeys(str(poi_name) for poi_name in scanned_pois))
        if poi not in scanned_pois:
            scanned_pois.insert(0, poi)
        tracked_pois = list(dict.fromkeys([*all_pois, *scanned_pois]))

        use_profile_endpoints = fit_algo != "none"
        row_entries: list[dict[str, Any]] = []
        for entry_index in range(int(tree.GetEntries())):
            tree.GetEntry(entry_index)
            row: dict[str, Any] = {
                "entry": entry_index,
                "quantileExpected": leaf_value(tree, "quantileExpected"),
                "iToy": leaf_value(tree, "iToy"),
                "iSeed": leaf_value(tree, "iSeed"),
                "deltaNLL": leaf_value(tree, "deltaNLL"),
                pdf_index_name: leaf_value(tree, pdf_index_name),
            }
            for poi_name in scanned_pois:
                row[poi_name] = leaf_value(tree, poi_name)
            for poi_name in tracked_pois:
                row[f"trackedParam_{poi_name}"] = leaf_value(tree, f"trackedParam_{poi_name}")
                row[f"trackedError_{poi_name}"] = leaf_value(tree, f"trackedError_{poi_name}")
            row_entries.append(row)

        grouped_rows: dict[tuple[int, int], list[dict[str, Any]]] = {}
        fallback_toy = 0
        for row in row_entries:
            i_toy_value = row.get("iToy")
            i_seed_value = row.get("iSeed")
            if i_toy_value is None:
                if row.get("quantileExpected") is not None and math.isclose(
                    float(row["quantileExpected"]), -1.0, rel_tol=0.0, abs_tol=1e-5
                ):
                    fallback_toy += 1
                i_toy = fallback_toy
            else:
                i_toy = int(round(float(i_toy_value)))
            i_seed = int(round(float(i_seed_value))) if i_seed_value is not None else 0
            grouped_rows.setdefault((i_seed, i_toy), []).append(row)

        pulls: list[float] = []
        entries: list[dict[str, Any]] = []
        skipped: list[dict[str, Any]] = []
        per_poi_pulls: dict[str, list[float]] = {poi_name: [] for poi_name in scanned_pois}
        per_poi_skipped: dict[str, int] = {poi_name: 0 for poi_name in scanned_pois}
        per_poi_sigma_sources: dict[str, dict[str, int]] = {poi_name: {} for poi_name in scanned_pois}
        for (i_seed, i_toy), rows in sorted(grouped_rows.items()):
            rows = sorted(rows, key=lambda item: int(item["entry"]))
            best_rows = [
                row
                for row in rows
                if row.get("quantileExpected") is not None
                and math.isclose(float(row["quantileExpected"]), -1.0, rel_tol=0.0, abs_tol=1e-5)
            ]
            if not best_rows:
                for poi_name in scanned_pois:
                    per_poi_skipped[poi_name] += 1
                skipped.append(
                    {
                        "iSeed": i_seed,
                        "iToy": i_toy,
                        "skip_reason": "missing_best_fit_entry",
                    }
                )
                continue
            best = best_rows[0]
            best_position = rows.index(best)
            for current_poi_index, current_poi in enumerate(scanned_pois):
                r_hat = best.get(current_poi)
                tracked_error = best.get(f"trackedError_{current_poi}")
                lower_row_index = best_position + 1 + 2 * current_poi_index
                upper_row_index = best_position + 2 + 2 * current_poi_index
                lower_value = None
                upper_value = None
                if use_profile_endpoints and upper_row_index < len(rows):
                    lower_value = rows[lower_row_index].get(current_poi)
                    upper_value = rows[upper_row_index].get(current_poi)
                record = {
                    "entry": best.get("entry"),
                    "iSeed": i_seed,
                    "iToy": i_toy,
                    "poi": current_poi,
                    "injected_r_reference": injected_r if current_poi == poi else 0.0,
                    "r_hat": r_hat,
                    "tracked_error": tracked_error,
                    "lower_endpoint": lower_value,
                    "upper_endpoint": upper_value,
                    "pdfindex_fit": best.get(pdf_index_name),
                    "deltaNLL": best.get("deltaNLL"),
                }
                if r_hat is None:
                    record["skip_reason"] = "missing_value"
                    per_poi_skipped[current_poi] += 1
                    if current_poi == poi:
                        skipped.append(record)
                    continue

                sigma_source = "tracked_error" if fit_algo == "none" else "profile_endpoints"
                sigma = None
                if use_profile_endpoints and lower_value is not None and upper_value is not None:
                    lo_err = abs(float(r_hat) - float(lower_value))
                    hi_err = abs(float(upper_value) - float(r_hat))
                    sigma = 0.5 * (lo_err + hi_err)
                    record["lo_err"] = lo_err
                    record["hi_err"] = hi_err
                if sigma is None or not math.isfinite(sigma) or sigma <= 0.0:
                    if tracked_error is not None:
                        sigma = abs(float(tracked_error))
                        sigma_source = "tracked_error"
                if sigma is None or not math.isfinite(sigma) or sigma <= 0.0:
                    record["skip_reason"] = "invalid_error"
                    record["sigma"] = sigma
                    per_poi_skipped[current_poi] += 1
                    if current_poi == poi:
                        skipped.append(record)
                    continue
                pull = (float(r_hat) - float(record["injected_r_reference"])) / sigma
                if not math.isfinite(pull):
                    record["skip_reason"] = "invalid_pull"
                    record["sigma"] = sigma
                    per_poi_skipped[current_poi] += 1
                    if current_poi == poi:
                        skipped.append(record)
                    continue
                record["sigma"] = sigma
                record["sigma_source"] = sigma_source
                record["pull"] = pull
                per_poi_pulls[current_poi].append(pull)
                per_poi_sigma_sources[current_poi][sigma_source] = (
                    per_poi_sigma_sources[current_poi].get(sigma_source, 0) + 1
                )
                if current_poi == poi:
                    pulls.append(pull)
                    entries.append(record)

        all_poi_pull_summaries = {
            poi_name: {
                "fit_algo": fit_algo,
                "reference_injected_r": injected_r if poi_name == poi else 0.0,
                "sample_mean": sample_mean(poi_pulls),
                "sample_sigma": sample_std(poi_pulls),
                "entries_used": len(poi_pulls),
                "entries_skipped": per_poi_skipped.get(poi_name, 0),
                "sigma_source": (
                    next(iter(per_poi_sigma_sources[poi_name]))
                    if len(per_poi_sigma_sources[poi_name]) == 1
                    else "mixed"
                    if per_poi_sigma_sources[poi_name]
                    else None
                ),
                "sigma_sources": per_poi_sigma_sources[poi_name],
            }
            for poi_name, poi_pulls in per_poi_pulls.items()
        }
        other_poi_pull_means = {
            poi_name: summary["sample_mean"]
            for poi_name, summary in all_poi_pull_summaries.items()
            if poi_name != poi
        }

        return {
            "status": "ok",
            "fit_root": str(fit_root_path.resolve()),
            "fit_root_repo_relative": repo_relative(fit_root_path),
            "tree": "limit",
            "branches": branches,
            "poi": poi,
            "all_pois": all_pois,
            "scanned_pois": scanned_pois,
            "fit_algo": fit_algo,
            "injected_r": injected_r,
            "pulls": pulls,
            "entries": entries,
            "skipped_entries": skipped,
            "all_poi_pull_summaries": all_poi_pull_summaries,
            "other_poi_pull_means": other_poi_pull_means,
            "entries_total": len(grouped_rows),
            "tree_entries_total": int(tree.GetEntries()),
            "entries_used": len(entries),
            "entries_skipped": len(skipped),
        }
    finally:
        root_file.Close()


def sample_mean(values: list[float]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def sample_std(values: list[float]) -> float | None:
    if len(values) < 2:
        return None
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
    return math.sqrt(max(0.0, variance))


def set_marker_alpha(graph: Any, color: int, alpha: float) -> None:
    if hasattr(graph, "SetMarkerColorAlpha"):
        graph.SetMarkerColorAlpha(color, alpha)
    else:
        graph.SetMarkerColor(color)


def make_pull_distribution_plot(
    pulls: list[float],
    output_dir: Path,
    title: str,
    pull_range: tuple[float, float],
    bins: int,
    min_fit_entries: int,
) -> dict[str, Any]:
    from statistical_test_fit.root_runtime import configure_root

    ROOT = configure_root()
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / "pull_distribution.png"
    pdf_path = output_dir / "pull_distribution.pdf"
    hist_name = "h_pull_" + safe_name(title)
    hist = ROOT.TH1F(hist_name, title, bins, pull_range[0], pull_range[1])
    hist.SetDirectory(0)
    hist.SetStats(0)
    ROOT.gStyle.SetOptStat(0)
    hist.SetLineColor(ROOT.kAzure + 1)
    hist.SetFillColorAlpha(ROOT.kAzure - 9, 0.45)
    hist.GetXaxis().SetTitle("pull")
    hist.GetYaxis().SetTitle("pseudo-datasets")
    for pull in pulls:
        hist.Fill(float(pull))

    gaussian_status = "skipped"
    gaussian_mean = sample_mean(pulls)
    gaussian_sigma = sample_std(pulls)
    gaussian_mean_error = None
    gaussian_sigma_error = None
    fit_status = None
    fit_curve = None
    if len(pulls) >= min_fit_entries and hist.GetEntries() > 0:
        initial_norm = max(1.0, float(hist.GetMaximum()))
        initial_mean = float(gaussian_mean if gaussian_mean is not None else 0.0)
        initial_sigma = abs(float(gaussian_sigma if gaussian_sigma is not None else 1.0))
        if initial_sigma <= 0.0 or not math.isfinite(initial_sigma):
            initial_sigma = 1.0
        fit_curve = ROOT.TF1(
            "f_pull_normal_" + safe_name(title),
            "gaus",
            pull_range[0],
            pull_range[1],
        )
        fit_curve.SetParameters(initial_norm, initial_mean, initial_sigma)
        fit_curve.SetParNames("normalization", "mean", "sigma")
        fit_result = hist.Fit(fit_curve, "QSN")
        fit_status = int(fit_result.Status()) if fit_result is not None else None
        if fit_curve:
            gaussian_mean = float(fit_curve.GetParameter(1))
            gaussian_sigma = abs(float(fit_curve.GetParameter(2)))
            gaussian_mean_error = float(fit_curve.GetParError(1))
            gaussian_sigma_error = float(fit_curve.GetParError(2))
        if fit_curve and (fit_status is None or fit_status == 0):
            gaussian_status = "ok"
        else:
            gaussian_status = "failed"

    canvas = ROOT.TCanvas("c_pull", "c_pull", 1100, 760)
    canvas.SetRightMargin(0.34)
    canvas.SetGridx(True)
    canvas.SetGridy(True)
    hist.Draw("HIST")
    if fit_curve:
        fit_curve.SetLineColor(ROOT.kRed + 1 if gaussian_status == "ok" else ROOT.kOrange + 7)
        fit_curve.SetLineStyle(1 if gaussian_status == "ok" else 2)
        fit_curve.SetLineWidth(3)
        fit_curve.Draw("SAME")

    legend = ROOT.TLegend(0.70, 0.74, 0.98, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.030)
    legend.AddEntry(hist, "Pull distribution", "f")
    if fit_curve:
        fit_label = (
            "Fitted normal (red)"
            if gaussian_status == "ok"
            else "Normal fit attempt (orange)"
        )
        legend.AddEntry(fit_curve, fit_label, "l")
    legend.Draw()

    stats = ROOT.TPaveText(0.70, 0.48, 0.98, 0.70, "NDC")
    stats.SetBorderSize(0)
    stats.SetFillColor(0)
    stats.SetFillStyle(0)
    stats.SetTextAlign(12)
    stats.SetTextSize(0.030)
    stats.AddText(f"N used = {len(pulls)}")
    stats.AddText(f"mean = {gaussian_mean:.4g}" if gaussian_mean is not None else "mean = n/a")
    stats.AddText(f"sigma = {gaussian_sigma:.4g}" if gaussian_sigma is not None else "sigma = n/a")
    stats.AddText(f"normal fit = {gaussian_status}")
    if gaussian_mean_error is not None:
        stats.AddText(f"mean err = {gaussian_mean_error:.3g}")
    if gaussian_sigma_error is not None:
        stats.AddText(f"sigma err = {gaussian_sigma_error:.3g}")
    stats.Draw()
    canvas.SaveAs(str(png_path))
    canvas.SaveAs(str(pdf_path))

    return {
        "plot_png": str(png_path.resolve()),
        "plot_png_repo_relative": repo_relative(png_path),
        "plot_pdf": str(pdf_path.resolve()),
        "plot_pdf_repo_relative": repo_relative(pdf_path),
        "hist_entries": int(hist.GetEntries()),
        "pull_range": list(pull_range),
        "pull_bins": bins,
        "sample_mean": sample_mean(pulls),
        "sample_sigma": sample_std(pulls),
        "gaussian_status": gaussian_status,
        "gaussian_fit_status": fit_status,
        "gaussian_mean": gaussian_mean,
        "gaussian_mean_error": gaussian_mean_error,
        "gaussian_sigma": gaussian_sigma,
        "gaussian_sigma_error": gaussian_sigma_error,
        "normal_fit_attempted": fit_curve is not None,
        "normal_fit_drawn": fit_curve is not None,
    }


def analyze_fit_result(result: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    metadata = result.get("metadata", {})
    dataset_strategy = str(metadata.get("dataset_strategy", "toys"))
    fit_root = first_root_file(result, "higgsCombine")
    if fit_root is None:
        analysis = {
            "status": "failed",
            "message": "No higgsCombine*.root file produced by MultiDimFit.",
            "pulls": [],
        }
    else:
        analysis = extract_pull_values(
            fit_root,
            str(metadata["target_poi"]),
            list(metadata.get("all_pois", [metadata["target_poi"]])),
            float(metadata["injected_r"]),
            str(metadata.get("pdf_index_name", args.pdf_index_name)),
            str(metadata.get("fit_algo", FIT_ALGO)),
            list(
                metadata.get("scanned_pois")
                or metadata.get("all_pois", [metadata["target_poi"]])
            ),
        )
    if analysis.get("status") == "ok":
        pulls = [float(value) for value in analysis.get("pulls", [])]
        entries = [entry for entry in analysis.get("entries", []) if isinstance(entry, dict)]
        r_hats = [
            float(value)
            for entry in entries
            for value in [finite_float_or_none(entry.get("r_hat"))]
            if value is not None
        ]
        fit_sigmas = [
            float(value)
            for entry in entries
            for value in [finite_float_or_none(entry.get("sigma"))]
            if value is not None
        ]
        r_hat_mean = sample_mean(r_hats)
        fit_sigma_mean = sample_mean(fit_sigmas)
        r_closure = r_hat_mean - float(metadata["injected_r"]) if r_hat_mean is not None else None
        analysis["r_hat"] = r_hat_mean
        analysis["r_hat_sample_sigma"] = sample_std(r_hats)
        analysis["fit_sigma"] = fit_sigma_mean
        analysis["fit_sigma_sample_sigma"] = sample_std(fit_sigmas)
        analysis["r_closure"] = r_closure
        if dataset_strategy == "asimov":
            closure_pull = sample_mean(pulls)
            analysis["bias_metric"] = "asimov_closure_pull"
            analysis["asimov_closure_pull"] = closure_pull
            analysis["asimov_r_hat"] = r_hat_mean
            analysis["asimov_fit_sigma"] = fit_sigma_mean
            analysis["asimov_r_closure"] = r_closure
            analysis["asimov_closure_status"] = "ok" if closure_pull is not None else "failed"
            analysis["bias_flag"] = (
                bool(abs(float(closure_pull)) >= args.bias_pull_threshold)
                if closure_pull is not None
                else None
            )
            analysis["pull_width_flag"] = None
        else:
            plot = make_pull_distribution_plot(
                pulls,
                Path(result["cwd"]),
                (
                    f"{metadata.get('scheme')} {metadata.get('target_poi')} "
                    f"truth={metadata.get('truth_pdf_label')} "
                    f"target={metadata.get('fit_pdf_label', 'floating')} "
                    f"r={metadata.get('injected_r')}"
                ),
                args.pull_range,
                args.pull_bins,
                args.min_fit_entries,
            )
            analysis["bias_metric"] = "toy_pull_distribution_mean"
            analysis["pull_distribution"] = plot
            mean = plot.get("gaussian_mean")
            sigma = plot.get("gaussian_sigma")
            if plot.get("gaussian_status") == "ok":
                analysis["bias_flag"] = bool(
                    mean is not None and abs(float(mean)) >= args.bias_pull_threshold
                )
                analysis["pull_width_flag"] = bool(
                    sigma is not None and abs(float(sigma) - 1.0) >= args.pull_width_threshold
                )
            else:
                analysis["bias_flag"] = None
                analysis["pull_width_flag"] = None
    else:
        analysis["bias_flag"] = None
        analysis["pull_width_flag"] = None
    result["bias_analysis"] = analysis
    write_job_summary(result)
    return result


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
      --ok: #34d399;
      --fail: #fb7185;
      --warn: #fbbf24;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: radial-gradient(circle at top left, rgba(125, 211, 252, 0.13), transparent 32rem),
        linear-gradient(180deg, #070b12 0%, #0a101c 100%);
      color: var(--text);
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      line-height: 1.5;
    }}
    main {{ width: 90vw; max-width: none; margin: 0 auto; padding: 2.2rem; }}
    h1 {{ margin: 0 0 0.5rem; font-size: clamp(1.9rem, 4vw, 3.1rem); letter-spacing: -0.04em; }}
    h2 {{ margin: 0 0 1rem; color: var(--accent); font-size: 1.05rem; text-transform: uppercase; letter-spacing: 0.12em; }}
    .subtitle {{ color: var(--muted); margin: 0 0 2rem; }}
    .card {{
      background: linear-gradient(180deg, rgba(16, 24, 38, 0.96), rgba(12, 19, 32, 0.96));
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 1.2rem;
      box-shadow: 0 20px 60px rgba(0, 0, 0, 0.24);
      margin: 1rem 0;
      overflow: auto;
    }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; }}
    th, td {{ padding: 0.68rem 0.78rem; border-bottom: 1px solid var(--line); vertical-align: top; }}
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
    .warn {{ color: #422006; background: var(--warn); }}
    .muted {{ color: var(--muted); }}
    .plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 1rem; }}
    .plot-tile {{ background: rgba(5, 9, 20, 0.45); border: 1px solid var(--line); border-radius: 14px; padding: 1rem; }}
    .plot-tile h3 {{ margin: 0 0 0.35rem; font-size: 1rem; color: #dbeafe; }}
    .plot-links {{ margin: 0.35rem 0 0.9rem; color: var(--muted); }}
    img {{ max-width: 100%; border: 1px solid var(--line); border-radius: 14px; background: #fff; }}
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


def html_table_raw(headers: list[str], rows: list[list[Any]]) -> str:
    if not rows:
        rows = [["-" for _ in headers]]
    header_html = "".join(f"<th>{html_escape(header)}</th>" for header in headers)
    row_html = []
    for row in rows:
        row_html.append("<tr>" + "".join(f"<td>{value}</td>" for value in row) + "</tr>")
    return "<table><thead><tr>" + header_html + "</tr></thead><tbody>" + "".join(row_html) + "</tbody></table>"


def html_link(label: str, href: str) -> str:
    return f'<a href="{html_escape(href)}">{html_escape(label)}</a>'


def status_pill(status: Any) -> str:
    status_text = str(status)
    css = "ok" if status_text == "ok" else "failed" if status_text == "failed" else "muted"
    return f'<span class="pill {css}">{html_escape(status_text)}</span>'


def flag_pill(flag: Any) -> str:
    if flag is True:
        return '<span class="pill warn">flagged</span>'
    if flag is False:
        return '<span class="pill ok">ok</span>'
    return '<span class="pill muted">n/a</span>'


def html_fixed_number_or_dash(value: Any, digits: int = 3) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "-"
    if not math.isfinite(numeric):
        return "-"
    return html_escape(f"{numeric:.{digits}f}")


def scatter_plot_tile_html(plot: dict[str, Any], run_dir: Path) -> str:
    label = html_escape(plot.get("label", "-"))
    status = str(plot.get("status", "-"))
    points_label = "cells" if plot.get("cells") is not None else "points"
    points = html_escape(plot.get("cells", plot.get("points", "-")))
    links = []
    if plot.get("plot_png"):
        links.append(html_link("png", href_relative_to(run_dir, Path(str(plot["plot_png"])))))
    if plot.get("plot_pdf"):
        links.append(html_link("pdf", href_relative_to(run_dir, Path(str(plot["plot_pdf"])))))
    if plot.get("experiment_map"):
        links.append(html_link("map", href_relative_to(run_dir, Path(str(plot["experiment_map"])))))
    if plot.get("values_json"):
        links.append(html_link("values", href_relative_to(run_dir, Path(str(plot["values_json"])))))
    link_html = " | ".join(links) if links else "-"
    message = html_escape(plot.get("message", ""))
    message_html = f'<p class="muted">{message}</p>' if message else ""
    image_html = ""
    if status == "ok" and plot.get("plot_png"):
        image_href = href_relative_to(run_dir, Path(str(plot["plot_png"])))
        image_html = f'<p><img src="{html_escape(image_href)}" alt="{label}"></p>'
    return f"""
      <div class="plot-tile">
        <h3>{label}</h3>
        <p class="muted">status: {html_escape(status)}; {html_escape(points_label)}: {points}</p>
        {message_html}
        <p class="plot-links">{link_html}</p>
        {image_html}
      </div>
    """


def scatter_plot_grid_section(title: str, plots: list[dict[str, Any]], run_dir: Path) -> str:
    if not plots:
        return ""
    tiles = "\n".join(scatter_plot_tile_html(plot, run_dir) for plot in plots)
    return f"""
    <section class="card"><h2>{html_escape(title)}</h2>
      <div class="plot-grid">{tiles}</div>
    </section>
    """


def make_job_html_summary(result: dict[str, Any]) -> str:
    metadata = result.get("metadata", {})
    input_rows = [
        [item.get("staged_name", "-"), item.get("source_repo_relative", item.get("source", "-"))]
        for item in result.get("inputs", [])
    ]
    output_rows = [[name] for name in result.get("produced_root_files", [])]
    status_rows = [
        ["Kind", result.get("kind", "-")],
        ["Method", result.get("method", "-")],
        ["Status", result.get("status", "-")],
        ["Return code", result.get("returncode", "-")],
        ["Duration", result.get("duration", "-")],
        ["Working directory", result.get("cwd_repo_relative", result.get("cwd", "-"))],
    ]
    metadata_rows = [[key, json.dumps(value) if isinstance(value, (list, dict)) else value] for key, value in sorted(metadata.items())]
    logs_rows = [
        ["command", html_link("command.txt", "command.txt")],
        ["stdout", html_link("stdout.txt", "stdout.txt")],
        ["stderr", html_link("stderr.txt", "stderr.txt")],
        ["json", html_link("summary.json", "summary.json")],
    ]
    analysis = result.get("bias_analysis") or {}
    analysis_section = ""
    if analysis:
        plot = (analysis.get("pull_distribution") or {}).get("plot_png_repo_relative")
        analysis_rows = [
            ["Analysis status", analysis.get("status", "-")],
            ["Bias metric", analysis.get("bias_metric", "-")],
            ["Entries used", analysis.get("entries_used", "-")],
            ["Entries skipped", analysis.get("entries_skipped", "-")],
            ["Fitted r", analysis.get("r_hat", "-")],
            ["Fitted r sigma", analysis.get("fit_sigma", "-")],
            ["r closure", analysis.get("r_closure", "-")],
            ["Asimov closure pull", analysis.get("asimov_closure_pull", "-")],
            ["Gaussian mean", (analysis.get("pull_distribution") or {}).get("gaussian_mean", "-")],
            ["Gaussian sigma", (analysis.get("pull_distribution") or {}).get("gaussian_sigma", "-")],
            ["Bias flag", "yes" if analysis.get("bias_flag") else "no" if analysis.get("bias_flag") is False else "n/a"],
            ["Pull-width flag", "yes" if analysis.get("pull_width_flag") else "no" if analysis.get("pull_width_flag") is False else "n/a"],
        ]
        image_html = f'<p><img src="{html_escape(Path(plot).name)}" alt="pull distribution"></p>' if plot else ""
        analysis_section = f"""
    <section class="card"><h2>Pull Analysis</h2>{html_table(['Field', 'Value'], analysis_rows)}{image_html}</section>
        """
    body = f"""
    <h1>{html_escape(result['job_id'])}</h1>
    <p class="subtitle">Per-job Combine execution summary</p>
    <section class="card"><h2>Status</h2>{html_table(['Field', 'Value'], status_rows)}</section>
    <section class="card"><h2>Command</h2><pre>{html_escape(result.get('command_string', ''))}</pre></section>
    <section class="card"><h2>Metadata</h2>{html_table(['Field', 'Value'], metadata_rows)}</section>
    {analysis_section}
    <section class="card"><h2>Inputs</h2>{html_table(['Local file', 'Source'], input_rows)}</section>
    <section class="card"><h2>ROOT Outputs</h2>{html_table(['File'], output_rows)}</section>
    <section class="card"><h2>Logs</h2>{html_table_raw(['Stream', 'Path'], logs_rows)}</section>
    """
    return html_document(str(result["job_id"]), body)


def write_job_summary(result: dict[str, Any]) -> None:
    cwd = Path(result["cwd"])
    result["summary_json_path"] = str((cwd / "summary.json").resolve())
    result["summary_json_path_repo_relative"] = repo_relative(cwd / "summary.json")
    result["summary_html_path"] = str((cwd / "summary.html").resolve())
    result["summary_html_path_repo_relative"] = repo_relative(cwd / "summary.html")
    (cwd / "summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (cwd / "summary.html").write_text(make_job_html_summary(result), encoding="utf-8")


def short_result(result: dict[str, Any]) -> dict[str, Any]:
    metadata = result.get("metadata", {})
    return {
        "job_id": result.get("job_id"),
        "kind": result.get("kind"),
        "method": result.get("method"),
        "status": result.get("status"),
        "returncode": result.get("returncode"),
        "duration": result.get("duration"),
        "cwd": result.get("cwd_repo_relative", result.get("cwd")),
        "command_file": result.get("command_path_repo_relative", result.get("command_path")),
        "summary_json": result.get("summary_json_path_repo_relative", result.get("summary_json_path")),
        "summary_html": result.get("summary_html_path_repo_relative", result.get("summary_html_path")),
        "scheme": metadata.get("scheme"),
        "scheme_slug": metadata.get("scheme_slug"),
        "scheme_label": metadata.get("scheme_label"),
        "scheme_latex": metadata.get("scheme_latex"),
        "target_poi": metadata.get("target_poi"),
        "target_label": metadata.get("target_label"),
        "target_slug": metadata.get("target_slug"),
        "target_display_label": metadata.get("target_display_label"),
        "target_latex": metadata.get("target_latex"),
        "scanned_pois": metadata.get("scanned_pois"),
        "profiled_pois": metadata.get("profiled_pois"),
        "poi_scan_mode": metadata.get("poi_scan_mode"),
        "float_other_pois": metadata.get("float_other_pois"),
        "truth_pdf_index": metadata.get("truth_pdf_index"),
        "truth_pdf_label": metadata.get("truth_pdf_label"),
        "truth_pdf_slug": metadata.get("truth_pdf_slug"),
        "truth_pdf_latex": metadata.get("truth_pdf_latex"),
        "truth_pdf_family": metadata.get("truth_pdf_family"),
        "truth_pdf_scan_order": metadata.get("truth_pdf_scan_order"),
        "truth_pdf_selection_role": metadata.get("truth_pdf_selection_role"),
        "fit_algo": metadata.get("fit_algo"),
        "fit_pdf_mode": metadata.get("fit_pdf_mode"),
        "fit_pdf_index": metadata.get("fit_pdf_index"),
        "fit_pdf_label": metadata.get("fit_pdf_label"),
        "fit_pdf_slug": metadata.get("fit_pdf_slug"),
        "fit_pdf_latex": metadata.get("fit_pdf_latex"),
        "fit_pdf_family": metadata.get("fit_pdf_family"),
        "target_pdf_index": metadata.get("target_pdf_index"),
        "target_pdf_label": metadata.get("target_pdf_label"),
        "target_pdf_slug": metadata.get("target_pdf_slug"),
        "target_pdf_latex": metadata.get("target_pdf_latex"),
        "target_pdf_family": metadata.get("target_pdf_family"),
        "injected_r": metadata.get("injected_r"),
        "dataset_strategy": metadata.get("dataset_strategy"),
        "pseudo_datasets": metadata.get("pseudo_datasets"),
        "requested_poi_min": metadata.get("requested_poi_min"),
        "requested_poi_max": metadata.get("requested_poi_max"),
        "poi_min": metadata.get("poi_min"),
        "poi_max": metadata.get("poi_max"),
        "quick": metadata.get("quick"),
        "cmin_default_minimizer_strategy": metadata.get("cmin_default_minimizer_strategy"),
        "cmin_strategy_explicit": metadata.get("cmin_strategy_explicit"),
        "robust_fit": metadata.get("robust_fit"),
        "bias_analysis": result.get("bias_analysis"),
    }


def experiment_summary_from_fit(result: dict[str, Any], index: int) -> dict[str, Any]:
    metadata = result.get("metadata", {})
    analysis = result.get("bias_analysis") or {}
    distribution = analysis.get("pull_distribution") or {}
    raw_scheme = str(metadata.get("scheme", ""))
    raw_target = str(metadata.get("target_label") or metadata.get("target_poi") or "")
    scheme_display = poi_scheme_display(raw_scheme)
    target_display = signal_target_display(raw_target)
    truth_display = pdf_state_display(
        index=metadata.get("truth_pdf_index"),
        family=metadata.get("truth_pdf_family"),
        order=metadata.get("truth_pdf_scan_order"),
        selection_role=metadata.get("truth_pdf_selection_role"),
        name=metadata.get("truth_pdf"),
    )
    return {
        "experiment_index": index,
        "job_id": result.get("job_id"),
        "status": result.get("status"),
        "scheme": metadata.get("scheme"),
        "scheme_slug": metadata.get("scheme_slug") or scheme_display.slug,
        "scheme_label": metadata.get("scheme_label") or scheme_display.text,
        "scheme_latex": metadata.get("scheme_latex") or scheme_display.latex,
        "target_label": metadata.get("target_label"),
        "target_slug": metadata.get("target_slug") or target_display.slug,
        "target_display_label": metadata.get("target_display_label") or target_display.text,
        "target_latex": metadata.get("target_latex") or target_display.latex,
        "target_poi": metadata.get("target_poi"),
        "scanned_pois": metadata.get("scanned_pois"),
        "profiled_pois": metadata.get("profiled_pois"),
        "poi_scan_mode": metadata.get("poi_scan_mode"),
        "float_other_pois": metadata.get("float_other_pois"),
        "truth_pdf_index": metadata.get("truth_pdf_index"),
        "truth_pdf": metadata.get("truth_pdf"),
        "truth_pdf_label": metadata.get("truth_pdf_label") or truth_display.text,
        "truth_pdf_slug": metadata.get("truth_pdf_slug") or truth_display.slug,
        "truth_pdf_latex": metadata.get("truth_pdf_latex") or truth_display.latex,
        "truth_pdf_family": metadata.get("truth_pdf_family"),
        "truth_pdf_scan_order": metadata.get("truth_pdf_scan_order"),
        "truth_pdf_selection_role": metadata.get("truth_pdf_selection_role"),
        "fit_algo": metadata.get("fit_algo"),
        "fit_pdf_mode": metadata.get("fit_pdf_mode"),
        "fit_pdf_index": metadata.get("fit_pdf_index"),
        "fit_pdf": metadata.get("fit_pdf"),
        "fit_pdf_label": metadata.get("fit_pdf_label"),
        "fit_pdf_slug": metadata.get("fit_pdf_slug"),
        "fit_pdf_latex": metadata.get("fit_pdf_latex"),
        "fit_pdf_family": metadata.get("fit_pdf_family"),
        "target_pdf_index": metadata.get("target_pdf_index"),
        "target_pdf": metadata.get("target_pdf"),
        "target_pdf_label": metadata.get("target_pdf_label"),
        "target_pdf_slug": metadata.get("target_pdf_slug"),
        "target_pdf_latex": metadata.get("target_pdf_latex"),
        "target_pdf_family": metadata.get("target_pdf_family"),
        "injected_r": metadata.get("injected_r"),
        "injection_mode": metadata.get("injection_mode"),
        "dataset_strategy": metadata.get("dataset_strategy", "toys"),
        "toys_requested": metadata.get("toys"),
        "pseudo_datasets_requested": metadata.get("pseudo_datasets"),
        "requested_poi_min": metadata.get("requested_poi_min"),
        "requested_poi_max": metadata.get("requested_poi_max"),
        "poi_min": metadata.get("poi_min"),
        "poi_max": metadata.get("poi_max"),
        "quick": metadata.get("quick"),
        "cmin_default_minimizer_strategy": metadata.get("cmin_default_minimizer_strategy"),
        "cmin_strategy_explicit": metadata.get("cmin_strategy_explicit"),
        "robust_fit": metadata.get("robust_fit"),
        "analysis_status": analysis.get("status"),
        "entries_total": analysis.get("entries_total"),
        "entries_used": analysis.get("entries_used"),
        "entries_skipped": analysis.get("entries_skipped"),
        "r_hat": analysis.get("r_hat"),
        "fit_sigma": analysis.get("fit_sigma"),
        "r_closure": analysis.get("r_closure"),
        "asimov_closure_pull": analysis.get("asimov_closure_pull"),
        "asimov_r_hat": analysis.get("asimov_r_hat"),
        "asimov_fit_sigma": analysis.get("asimov_fit_sigma"),
        "asimov_r_closure": analysis.get("asimov_r_closure"),
        "bias_metric": analysis.get("bias_metric"),
        "gaussian_mean": distribution.get("gaussian_mean"),
        "gaussian_mean_error": distribution.get("gaussian_mean_error"),
        "gaussian_sigma": distribution.get("gaussian_sigma"),
        "gaussian_sigma_error": distribution.get("gaussian_sigma_error"),
        "gaussian_status": distribution.get("gaussian_status"),
        "sample_mean": distribution.get("sample_mean"),
        "sample_sigma": distribution.get("sample_sigma"),
        "other_poi_pull_means": analysis.get("other_poi_pull_means", {}),
        "all_poi_pull_summaries": analysis.get("all_poi_pull_summaries", {}),
        "bias_flag": analysis.get("bias_flag"),
        "pull_width_flag": analysis.get("pull_width_flag"),
        "fit_summary_html": result.get("summary_html_path_repo_relative"),
        "pull_plot": distribution.get("plot_png_repo_relative"),
    }


def root_color_for_family(family: str) -> int:
    from statistical_test_fit.root_runtime import configure_root

    ROOT = configure_root()
    colors = {
        "johnson": ROOT.kAzure + 1,
        "bernstein": ROOT.kOrange + 7,
        "chebychev": ROOT.kGreen + 2,
        "other": ROOT.kMagenta + 1,
    }
    return int(colors.get(family, ROOT.kGray + 2))


def finite_float_or_none(value: Any) -> float | None:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def fit_target_pdf_label(experiment: dict[str, Any]) -> str:
    mode = normalize_pdf_target_strategy(experiment.get("fit_pdf_mode", "floating"))
    if mode == "fixed":
        label = experiment.get("target_pdf_label") or experiment.get("fit_pdf_label")
        index = experiment.get("target_pdf_index")
        if label:
            return str(label)
        if index is not None:
            return f"pdf{index}"
        return "fixed"
    return floating_pdf_display().text


def fit_target_pdf_plot_label(experiment: dict[str, Any]) -> str:
    mode = normalize_pdf_target_strategy(experiment.get("fit_pdf_mode", "floating"))
    if mode == "fixed":
        return str(
            experiment.get("target_pdf_latex")
            or experiment.get("fit_pdf_latex")
            or fit_target_pdf_label(experiment)
        )
    return floating_pdf_display().latex


def fit_target_pdf_sort_key(experiment: dict[str, Any]) -> tuple[int, int, str]:
    mode = normalize_pdf_target_strategy(experiment.get("fit_pdf_mode", "floating"))
    if mode == "fixed":
        index = experiment.get("target_pdf_index", experiment.get("fit_pdf_index"))
        try:
            pdf_index = int(index)
        except (TypeError, ValueError):
            pdf_index = 1_000_000
        return (0, pdf_index, fit_target_pdf_label(experiment))
    return (1, 1_000_000, "floating")


def scatter_point_values(experiment: dict[str, Any]) -> tuple[float, float, str] | None:
    if str(experiment.get("dataset_strategy", "toys")) == "asimov":
        closure_pull = finite_float_or_none(experiment.get("asimov_closure_pull"))
        if closure_pull is None:
            return None
        return closure_pull, 1.0, "asimov_closure"

    status = str(experiment.get("gaussian_status", ""))
    if status not in {"ok", "failed"}:
        return None
    mean = finite_float_or_none(experiment.get("gaussian_mean"))
    sigma = finite_float_or_none(experiment.get("gaussian_sigma"))
    source = "normal_fit_ok" if status == "ok" else "normal_fit_failed"
    if mean is None or sigma is None:
        return None
    return mean, abs(sigma), source


def experiment_full_label(experiment: dict[str, Any]) -> str:
    return (
        f"{experiment.get('target_display_label') or experiment.get('target_poi', '-')}\n"
        f"strategy = {experiment.get('dataset_strategy', 'toys')}\n"
        f"r = {format_number(float(experiment.get('injected_r', 0.0)))}\n"
        f"truth PDF = {experiment.get('truth_pdf_label', '?')}\n"
        f"target pdf = {fit_target_pdf_label(experiment)}"
    )


def make_scatter_plot(
    experiments: list[dict[str, Any]],
    run_dir: Path,
    bias_threshold: float,
    annotate_outliers: bool,
    plot_stem: str,
    label: str,
    split: dict[str, Any] | None = None,
    y_limits: tuple[float, float] | None = None,
) -> dict[str, Any]:
    plot_dir = run_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    png_path = plot_dir / f"{plot_stem}.png"
    pdf_path = plot_dir / f"{plot_stem}.pdf"
    map_path = plot_dir / f"{plot_stem}_experiments.json"
    split = dict(split or {})

    point_records: list[dict[str, Any]] = []
    map_records: list[dict[str, Any]] = []
    for fallback_index, experiment in enumerate(experiments, start=1):
        record = dict(experiment)
        record["experiment_index"] = int(record.get("experiment_index", fallback_index))
        values = scatter_point_values(experiment)
        if values is None:
            record["scatter_mean"] = None
            record["scatter_sigma"] = None
            record["scatter_source"] = None
            map_records.append(record)
            continue
        mean, sigma, source = values
        record["scatter_mean"] = mean
        record["scatter_sigma"] = sigma
        record["scatter_source"] = source
        map_records.append(record)
        point_records.append(record)
    map_path.write_text(json.dumps(map_records, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
    except Exception as exc:
        return {
            "status": "failed",
            "message": f"Could not import matplotlib for scatter plot: {exc}",
            "label": label,
            "split": split,
            "plot_stem": plot_stem,
            "points": len(point_records),
            "experiment_map": str(map_path.resolve()),
            "experiment_map_repo_relative": repo_relative(map_path),
        }

    n_points = len(point_records)
    if y_limits is None:
        y_values = [float(experiment["scatter_mean"]) for experiment in point_records]
        y_min = min(y_values + [-bias_threshold, 0.0])
        y_max = max(y_values + [bias_threshold, 0.0])
        margin = max(0.5, 0.15 * (y_max - y_min if y_max > y_min else 1.0))
        y_min -= margin
        y_max += margin
    else:
        y_min, y_max = y_limits

    pdf_index_colors = {
        0: "#1f77b4",
        1: "#ff7f0e",
        2: "#2ca02c",
        3: "#d62728",
        4: "#9467bd",
        5: "#8c564b",
        6: "#e377c2",
    }
    point_records.sort(
        key=lambda item: (
            str(item.get("scheme_slug") or item.get("scheme", "")),
            str(item.get("target_slug") or item.get("target_poi", "")),
            str(item.get("dataset_strategy", "toys")),
            float(item.get("injected_r", 0.0)),
            str(item.get("truth_pdf_slug") or item.get("truth_pdf_label", "")),
            fit_target_pdf_sort_key(item),
        )
    )
    for index, experiment in enumerate(point_records, start=1):
        experiment["scatter_x"] = index

    marker_cycle = ("o", "s", "^", "D", "P", "X", "v", "<", ">", "*", "h", "8")
    injection_values = sorted({float(experiment.get("injected_r", 0.0)) for experiment in point_records})
    injection_markers = {
        injected_r: marker_cycle[index % len(marker_cycle)]
        for index, injected_r in enumerate(injection_values)
    }

    side = min(max(15.0, 0.50 * max(1, n_points) + 8.0), 38.0)
    fig_width = side * 1.45
    fig, ax = plt.subplots(figsize=(fig_width, side), constrained_layout=False)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#f8fafc")
    ax.grid(True, axis="y", color="#cbd5e1", alpha=0.75, linewidth=0.8)
    ax.grid(True, axis="x", color="#e2e8f0", alpha=0.45, linewidth=0.5)

    if point_records:
        x_values = [int(experiment["scatter_x"]) for experiment in point_records]
        means = [float(experiment["scatter_mean"]) for experiment in point_records]
        sigmas = [abs(float(experiment["scatter_sigma"])) for experiment in point_records]
        colors = [
            pdf_index_colors.get(int(experiment.get("truth_pdf_index", -1)), "#7f7f7f")
            for experiment in point_records
        ]
        statuses = [str(experiment.get("scatter_source", "normal_fit_failed")) for experiment in point_records]
        injected_rs = [float(experiment.get("injected_r", 0.0)) for experiment in point_records]
        sizes = [max(180.0, min(3400.0, 220.0 + 420.0 * min(sigma, 7.0))) for sigma in sigmas]

        for status in ("normal_fit_ok", "normal_fit_failed", "asimov_closure"):
            for injected_r in injection_values:
                indices = [
                    index
                    for index, value in enumerate(statuses)
                    if value == status and injected_rs[index] == injected_r
                ]
                if not indices:
                    continue
                edge_color = "#0f172a"
                alpha = 0.55
                line_width = 0.85
                if status == "normal_fit_failed":
                    edge_color = "#dc2626"
                    alpha = 0.30
                    line_width = 3.0
                elif status == "asimov_closure":
                    edge_color = "#7c3aed"
                    alpha = 0.72
                    line_width = 2.4
                ax.scatter(
                    [x_values[index] for index in indices],
                    [means[index] for index in indices],
                    s=[sizes[index] for index in indices],
                    c=[colors[index] for index in indices],
                    marker=injection_markers[injected_r],
                    alpha=alpha,
                    edgecolors=edge_color,
                    linewidths=line_width,
                    zorder=3,
                )

        if annotate_outliers:
            for experiment in point_records:
                mean = float(experiment["scatter_mean"])
                if abs(mean) < bias_threshold:
                    continue
                ax.annotate(
                    experiment_full_label(experiment),
                    xy=(float(experiment["scatter_x"]), mean),
                    xytext=(14, 14 if mean >= 0.0 else -28),
                    textcoords="offset points",
                    fontsize=20,
                    color="#111827",
                    ha="left",
                    va="bottom" if mean >= 0.0 else "top",
                    arrowprops={"arrowstyle": "-", "color": "#64748b", "linewidth": 1.2},
                    bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "#cbd5e1", "alpha": 0.78},
                    zorder=5,
                )

    ax.axhline(0.0, color="#0f172a", linewidth=1.6, zorder=2)
    ax.axhline(bias_threshold, color="#dc2626", linestyle="--", linewidth=1.4, zorder=2)
    ax.axhline(-bias_threshold, color="#dc2626", linestyle="--", linewidth=1.4, zorder=2)
    ax.set_xlim(0.5, max(1, n_points) + 0.5)
    ax.set_yscale("symlog", linthresh=max(1e-3, bias_threshold), linscale=1.0)
    ax.set_ylim(y_min, y_max)
    ax.set_ylabel("Toy pull mean / Asimov closure pull", fontsize=34)
    ax.set_xlabel("", fontsize=34)
    ax.tick_params(axis="y", labelsize=28)
    ax.tick_params(axis="x", labelsize=24)

    max_x_ticks = 24
    if n_points <= max_x_ticks:
        x_ticks = list(range(1, n_points + 1))
    else:
        step = max(1, math.ceil((n_points - 1) / (max_x_ticks - 1)))
        x_ticks = list(range(1, n_points + 1, step))
        if x_ticks[-1] != n_points:
            x_ticks.append(n_points)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([str(index) for index in x_ticks], rotation=0, ha="center")

    group_start = 1
    groups: list[tuple[int, int, str]] = []
    while group_start <= n_points:
        current = point_records[group_start - 1]
        group_key = (current.get("scheme"), current.get("target_poi"))
        group_end = group_start
        while group_end <= n_points:
            candidate = point_records[group_end - 1]
            if (candidate.get("scheme"), candidate.get("target_poi")) != group_key:
                break
            group_end += 1
        end_index = group_end - 1
        group_label = str(current.get("target_latex") or current.get("target_display_label") or group_key[1])
        groups.append((group_start, end_index, group_label))
        group_start = group_end

    for group_index, (start, end, group_label) in enumerate(groups):
        if group_index % 2 == 0:
            ax.axvspan(start - 0.5, end + 0.5, color="#e2e8f0", alpha=0.32, zorder=0)
        if start > 1:
            ax.axvline(start - 0.5, color="#64748b", linewidth=1.0, alpha=0.65, zorder=1)
        midpoint = 0.5 * (start + end)
        ax.text(
            midpoint,
            1.015,
            group_label,
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=20,
            color="#334155",
            rotation=0,
            clip_on=False,
        )

    pdf_labels: dict[int, str] = {}
    for experiment in point_records:
        pdf_index = int(experiment.get("truth_pdf_index", -1))
        pdf_label = str(experiment.get("truth_pdf_latex") or experiment.get("truth_pdf_label", f"pdf{pdf_index}"))
        pdf_labels.setdefault(pdf_index, pdf_label)
    pdf_handles = [
        Patch(
            facecolor=pdf_index_colors.get(pdf_index, "#7f7f7f"),
            edgecolor="#0f172a",
            alpha=0.62,
            label=f"Truth PDF index {pdf_index}: {pdf_labels[pdf_index]}",
        )
        for pdf_index in sorted(pdf_labels)
    ]
    injection_handles = [
        Line2D(
            [0],
            [0],
            marker=injection_markers[injected_r],
            color="w",
            label=f"Injected r = {format_number(injected_r)}",
            markerfacecolor="#64748b",
            markeredgecolor="#0f172a",
            markersize=20,
        )
        for injected_r in injection_values
    ]
    status_handles = []
    if any(experiment.get("scatter_source") == "normal_fit_ok" for experiment in point_records):
        status_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label="Filled: normal fit converged; y = fitted mean, area = fitted sigma",
                markerfacecolor="#64748b",
                markeredgecolor="#0f172a",
                markersize=20,
            )
        )
    if any(experiment.get("scatter_source") == "normal_fit_failed" for experiment in point_records):
        status_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label="Red outline: normal fit failed; y/area from attempted fit",
                markerfacecolor="#64748b",
                markeredgecolor="#dc2626",
                markeredgewidth=3.0,
                markersize=20,
            )
        )
    if any(experiment.get("scatter_source") == "asimov_closure" for experiment in point_records):
        status_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label="Purple outline: Asimov closure pull; area is fixed",
                markerfacecolor="#64748b",
                markeredgecolor="#7c3aed",
                markeredgewidth=2.4,
                markersize=20,
            )
        )
    handles = pdf_handles + injection_handles + status_handles
    if handles:
        ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            frameon=False,
            borderaxespad=0.0,
            title="",
            title_fontsize=30,
            fontsize=20,
            labelspacing=0.9,
            handlelength=1.6,
        )

    if not point_records:
        ax.text(0.5, 0.5, "No valid pull summaries available", transform=ax.transAxes, ha="center", va="center", fontsize=28)

    bottom_margin = 0.20 if n_points <= 40 else 0.24
    fig.subplots_adjust(left=0.14, right=0.68, top=0.86, bottom=bottom_margin)
    fig.savefig(png_path, dpi=180)
    fig.savefig(pdf_path)
    plt.close(fig)

    map_path.write_text(json.dumps(map_records, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    annotation_note = (
        " Points above the bias threshold are annotated with the full name."
        if annotate_outliers
        else " Outlier annotations are disabled."
    )
    return {
        "status": "ok",
        "label": label,
        "split": split,
        "plot_stem": plot_stem,
        "points": n_points,
        "y_min": y_min,
        "y_max": y_max,
        "plot_png": str(png_path.resolve()),
        "plot_png_repo_relative": repo_relative(png_path),
        "plot_pdf": str(pdf_path.resolve()),
        "plot_pdf_repo_relative": repo_relative(pdf_path),
        "experiment_map": str(map_path.resolve()),
        "experiment_map_repo_relative": repo_relative(map_path),
        "annotate_outliers": annotate_outliers,
        "note": ("" if label == "Global" else f"Filtered to {label}. ") + "X axis shows index. Shaded bands and top labels group points by scheme and target POI. Point color shows truth PDF, marker shape shows injected signal strength, and point area shows fitted sigma for toy fits. Fixed and floating target PDF fits are separate points. Asimov closure points use fixed area." + annotation_note,
    }


def scatter_y_limits(
    experiments: list[dict[str, Any]],
    bias_threshold: float,
) -> tuple[float, float]:
    y_values: list[float] = []
    for experiment in experiments:
        values = scatter_point_values(experiment)
        if values is not None:
            y_values.append(float(values[0]))
    y_min = min(y_values + [-bias_threshold, 0.0])
    y_max = max(y_values + [bias_threshold, 0.0])
    margin = max(0.5, 0.15 * (y_max - y_min if y_max > y_min else 1.0))
    return y_min - margin, y_max + margin


def same_injected_r(experiment: dict[str, Any], injected_r: float) -> bool:
    value = finite_float_or_none(experiment.get("injected_r"))
    return value is not None and math.isclose(value, injected_r, rel_tol=0.0, abs_tol=1e-12)


def make_summary_scatter_plots(
    experiments: list[dict[str, Any]],
    run_dir: Path,
    bias_threshold: float,
    annotate_outliers: bool,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    y_limits = scatter_y_limits(experiments, bias_threshold)
    schemes = sorted(
        {str(experiment.get("scheme")) for experiment in experiments if experiment.get("scheme") is not None}
    )
    injected_rs = sorted(
        {
            value
            for experiment in experiments
            for value in [finite_float_or_none(experiment.get("injected_r"))]
            if value is not None
        }
    )
    total_plots = 1 + len(schemes) + len(injected_rs) + len(schemes) * len(injected_rs)
    completed_plots = 0

    log(
        f"Building {total_plots} summary scatter plot(s) from "
        f"{len(experiments)} fit summaries"
    )
    log(f"Summary scatter plot {completed_plots + 1}/{total_plots}: Global")
    global_plot = make_scatter_plot(
        experiments,
        run_dir,
        bias_threshold,
        annotate_outliers,
        "pull_mean_sigma_scatter",
        "Global",
        {},
        y_limits,
    )
    completed_plots += 1
    log(f"Finished summary scatter plot {completed_plots}/{total_plots}: Global")

    variations: dict[str, list[dict[str, Any]]] = {
        "by_poi_scheme": [],
        "by_injected_r": [],
        "by_poi_scheme_and_injected_r": [],
    }

    for scheme in schemes:
        display = poi_scheme_display(scheme)
        label = f"POI scheme {display.text}"
        log(f"Summary scatter plot {completed_plots + 1}/{total_plots}: {label}")
        variations["by_poi_scheme"].append(
            make_scatter_plot(
                [experiment for experiment in experiments if str(experiment.get("scheme")) == scheme],
                run_dir,
                bias_threshold,
                annotate_outliers,
                f"pull_mean_sigma_scatter__scheme_{display.slug}",
                label,
                {"scheme": scheme, "scheme_slug": display.slug, "scheme_label": display.text, "scheme_latex": display.latex},
                y_limits,
            )
        )
        completed_plots += 1
        log(f"Finished summary scatter plot {completed_plots}/{total_plots}: {label}")

    for injected_r in injected_rs:
        label = f"Injected r {format_number(injected_r)}"
        log(f"Summary scatter plot {completed_plots + 1}/{total_plots}: {label}")
        variations["by_injected_r"].append(
            make_scatter_plot(
                [experiment for experiment in experiments if same_injected_r(experiment, injected_r)],
                run_dir,
                bias_threshold,
                annotate_outliers,
                f"pull_mean_sigma_scatter__{injection_output_slug(injected_r)}",
                label,
                {"injected_r": injected_r},
                y_limits,
            )
        )
        completed_plots += 1
        log(f"Finished summary scatter plot {completed_plots}/{total_plots}: {label}")

    for scheme in schemes:
        display = poi_scheme_display(scheme)
        for injected_r in injected_rs:
            label = f"POI scheme {display.text}; injected r {format_number(injected_r)}"
            log(f"Summary scatter plot {completed_plots + 1}/{total_plots}: {label}")
            variations["by_poi_scheme_and_injected_r"].append(
                make_scatter_plot(
                    [
                        experiment
                        for experiment in experiments
                        if str(experiment.get("scheme")) == scheme
                        and same_injected_r(experiment, injected_r)
                    ],
                    run_dir,
                    bias_threshold,
                    annotate_outliers,
                    f"pull_mean_sigma_scatter__scheme_{display.slug}__{injection_output_slug(injected_r)}",
                    label,
                    {"scheme": scheme, "scheme_slug": display.slug, "scheme_label": display.text, "scheme_latex": display.latex, "injected_r": injected_r},
                    y_limits,
                )
            )
            completed_plots += 1
            log(f"Finished summary scatter plot {completed_plots}/{total_plots}: {label}")

    return global_plot, variations


def failed_summary_scatter_plot(
    run_dir: Path,
    message: str,
    returncode: int | None = None,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    plot_dir = run_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "status": "failed",
        "label": "Global",
        "split": {},
        "plot_stem": "pull_mean_sigma_scatter",
        "points": 0,
        "message": message,
        "stdout_path": str((plot_dir / "summary_scatter_stdout.txt").resolve()),
        "stdout_path_repo_relative": repo_relative(plot_dir / "summary_scatter_stdout.txt"),
        "stderr_path": str((plot_dir / "summary_scatter_stderr.txt").resolve()),
        "stderr_path_repo_relative": repo_relative(plot_dir / "summary_scatter_stderr.txt"),
    }
    if returncode is not None:
        payload["returncode"] = returncode
    return payload, {
        "by_poi_scheme": [],
        "by_injected_r": [],
        "by_poi_scheme_and_injected_r": [],
    }


def run_streamed_plot_subprocess(
    command: list[str],
    cwd: Path,
    stdout_path: Path,
    stderr_path: Path,
    timeout_seconds: int,
    label: str,
) -> dict[str, Any]:
    stdout_chunks: list[str] = []
    stderr_chunks: list[str] = []

    def read_stream(stream: Any, chunks: list[str], stream_label: str) -> None:
        if stream is None:
            return
        try:
            for line in stream:
                chunks.append(line)
                text = line.rstrip()
                if not text:
                    continue
                if stream_label == "stdout":
                    CONSOLE.print(Text(text))
                else:
                    CONSOLE.print(Text(f"[bias_study] {label} stderr:", style="yellow"), Text(text))
        finally:
            stream.close()

    log(f"{label}: starting plotting subprocess")
    process = subprocess.Popen(
        command,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )
    stdout_thread = threading.Thread(
        target=read_stream,
        args=(process.stdout, stdout_chunks, "stdout"),
        daemon=True,
    )
    stderr_thread = threading.Thread(
        target=read_stream,
        args=(process.stderr, stderr_chunks, "stderr"),
        daemon=True,
    )
    stdout_thread.start()
    stderr_thread.start()

    timed_out = False
    try:
        returncode = int(process.wait(timeout=timeout_seconds))
    except subprocess.TimeoutExpired:
        timed_out = True
        log(f"{label}: timed out after {timeout_seconds} seconds; terminating subprocess")
        process.kill()
        returncode = int(process.wait())

    stdout_thread.join(timeout=5.0)
    stderr_thread.join(timeout=5.0)
    if timed_out:
        stderr_chunks.append(f"Timed out after {timeout_seconds} seconds.\n")
    stdout_path.write_text("".join(stdout_chunks), encoding="utf-8")
    stderr_path.write_text("".join(stderr_chunks), encoding="utf-8")
    log(f"{label}: subprocess finished with return code {returncode}")
    return {
        "returncode": returncode,
        "timed_out": timed_out,
        "stdout": "".join(stdout_chunks),
        "stderr": "".join(stderr_chunks),
    }


def make_summary_scatter_plots_isolated(
    experiments: list[dict[str, Any]],
    run_dir: Path,
    bias_threshold: float,
    annotate_outliers: bool,
) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    plot_dir = run_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    input_path = plot_dir / "summary_scatter_input.json"
    output_path = plot_dir / "summary_scatter_output.json"
    stdout_path = plot_dir / "summary_scatter_stdout.txt"
    stderr_path = plot_dir / "summary_scatter_stderr.txt"
    input_payload = {
        "experiments": experiments,
        "run_dir": str(run_dir.resolve()),
        "bias_threshold": bias_threshold,
        "annotate_outliers": annotate_outliers,
    }
    input_path.write_text(json.dumps(input_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    code = r"""
import json
import sys
from pathlib import Path
import scripts.bias_study as bias_study

input_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
payload = json.loads(input_path.read_text(encoding='utf-8'))
scatter_plot, scatter_plot_variations = bias_study.make_summary_scatter_plots(
    payload['experiments'],
    Path(payload['run_dir']),
    float(payload['bias_threshold']),
    bool(payload['annotate_outliers']),
)
output_path.write_text(
    json.dumps(
        {
            'scatter_plot': scatter_plot,
            'scatter_plot_variations': scatter_plot_variations,
        },
        indent=2,
        sort_keys=True,
    ) + '\n',
    encoding='utf-8',
)
"""
    completed = run_streamed_plot_subprocess(
        [sys.executable, "-u", "-c", code, str(input_path), str(output_path)],
        REPO_ROOT,
        stdout_path,
        stderr_path,
        SCATTER_PLOT_SUBPROCESS_TIMEOUT_SECONDS,
        "summary scatter plots",
    )
    if completed["timed_out"]:
        return failed_summary_scatter_plot(
            run_dir,
            f"Summary scatter plotting timed out after {SCATTER_PLOT_SUBPROCESS_TIMEOUT_SECONDS} seconds.",
        )

    if int(completed["returncode"]) != 0:
        message = f"Summary scatter plotting subprocess failed with return code {completed['returncode']}."
        if int(completed["returncode"]) < 0:
            message += f" Signal {-int(completed['returncode'])} terminated the subprocess."
        return failed_summary_scatter_plot(run_dir, message, int(completed["returncode"]))
    if not output_path.exists():
        return failed_summary_scatter_plot(
            run_dir,
            "Summary scatter plotting subprocess finished without writing its output JSON.",
            int(completed["returncode"]),
        )

    try:
        output_payload = json.loads(output_path.read_text(encoding="utf-8"))
        return output_payload["scatter_plot"], output_payload["scatter_plot_variations"]
    except Exception as exc:
        return failed_summary_scatter_plot(
            run_dir,
            f"Could not read summary scatter plotting output JSON: {exc}",
            int(completed["returncode"]),
        )


def heatmap_metric_values(experiment: dict[str, Any]) -> tuple[float | None, float | None, str]:
    dataset_strategy = str(experiment.get("dataset_strategy", "toys"))
    if dataset_strategy == "asimov":
        metric = finite_float_or_none(experiment.get("asimov_closure_pull"))
        sigma = finite_float_or_none(experiment.get("asimov_fit_sigma"))
        if sigma is None:
            sigma = finite_float_or_none(experiment.get("fit_sigma"))
        return metric, sigma, "asimov_closure_pull"

    metric = finite_float_or_none(experiment.get("gaussian_mean"))
    sigma = finite_float_or_none(experiment.get("gaussian_sigma"))
    return metric, sigma, "pull_mean"


def experiment_fit_failed(experiment: dict[str, Any]) -> bool:
    return str(experiment.get("status", "ok")) != "ok" or str(experiment.get("analysis_status", "ok")) == "failed"


def truth_pdf_axis_key(experiment: dict[str, Any]) -> tuple[int, str]:
    try:
        index = int(experiment.get("truth_pdf_index"))
    except (TypeError, ValueError):
        index = 1_000_000
    label = str(experiment.get("truth_pdf_latex") or experiment.get("truth_pdf_label") or f"pdf{index}")
    return index, label


def target_pdf_axis_key(experiment: dict[str, Any]) -> tuple[int, int, str]:
    mode = normalize_pdf_target_strategy(experiment.get("fit_pdf_mode", "floating"))
    if mode == "fixed":
        index = experiment.get("target_pdf_index", experiment.get("fit_pdf_index"))
        try:
            pdf_index = int(index)
        except (TypeError, ValueError):
            pdf_index = 1_000_000
        return 0, pdf_index, fit_target_pdf_plot_label(experiment)
    return 1, 1_000_000, floating_pdf_display().latex


def summary_heatmap_row_key(experiment: dict[str, Any]) -> tuple[str, str, str, str, float, int, str]:
    truth_index, truth_label = truth_pdf_axis_key(experiment)
    injected_r = finite_float_or_none(experiment.get("injected_r"))
    return (
        str(experiment.get("dataset_strategy", "toys")),
        str(experiment.get("scheme", "-")),
        str(experiment.get("target_poi", "-")),
        str(experiment.get("target_latex") or experiment.get("target_display_label") or experiment.get("target_label") or experiment.get("target_poi", "-")),
        float(injected_r) if injected_r is not None else 0.0,
        truth_index,
        truth_label,
    )


def summary_heatmap_group_label(row_key: tuple[str, str, str, str, float, int, str]) -> str:
    dataset_strategy, scheme, _target_poi, target_label, injected_r, _truth_index, _truth_label = row_key
    return f"{dataset_strategy_display(dataset_strategy).text} | {scheme_latex(scheme)} | {target_label} | r={format_number(injected_r)}"


def make_pdf_target_summary_heatmap(
    experiments: list[dict[str, Any]],
    run_dir: Path,
    plt: Any,
    np: Any,
    TwoSlopeNorm: Any,
) -> dict[str, Any]:
    plot_dir = run_dir / "plots" / "pdf_target_heatmaps"
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_stem = "pdf_target_heatmap__all_tested_combinations"
    label = "All Tested Combinations"
    png_path = plot_dir / f"{plot_stem}.png"
    pdf_path = plot_dir / f"{plot_stem}.pdf"
    values_path = plot_dir / f"{plot_stem}_values.json"

    row_axis = sorted({summary_heatmap_row_key(experiment) for experiment in experiments})
    target_axis = sorted({target_pdf_axis_key(experiment) for experiment in experiments})
    if not row_axis or not target_axis:
        payload = {
            "status": "failed",
            "message": "No tested truth/target PDF combinations available.",
            "label": label,
            "split": {"summary": "all_tested_combinations"},
            "plot_stem": plot_stem,
            "points": 0,
        }
        values_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload

    row_positions = {axis_key: index for index, axis_key in enumerate(row_axis)}
    target_positions = {axis_key: index for index, axis_key in enumerate(target_axis)}
    matrix = np.full((len(row_axis), len(target_axis)), np.nan, dtype=float)
    sigma_matrix = np.full((len(row_axis), len(target_axis)), np.nan, dtype=float)
    failed_matrix = np.full((len(row_axis), len(target_axis)), False, dtype=bool)
    cell_records: list[dict[str, Any]] = []

    for experiment in sorted(
        experiments,
        key=lambda item: (
            summary_heatmap_row_key(item),
            target_pdf_axis_key(item),
            str(item.get("job_id", "")),
        ),
    ):
        row_key = summary_heatmap_row_key(experiment)
        target_key = target_pdf_axis_key(experiment)
        row = row_positions[row_key]
        column = target_positions[target_key]
        metric, sigma, metric_name = heatmap_metric_values(experiment)
        fit_failed = experiment_fit_failed(experiment)
        failed_matrix[row, column] = failed_matrix[row, column] or fit_failed
        if not fit_failed and metric is not None:
            matrix[row, column] = float(metric)
        if not fit_failed and sigma is not None:
            sigma_matrix[row, column] = abs(float(sigma))
        dataset_strategy, scheme, target_poi, target_label, injected_r, truth_index, truth_label = row_key
        cell_records.append(
            {
                "dataset_strategy": dataset_strategy,
                "scheme": scheme,
                "target_poi": target_poi,
                "target_label": target_label,
                "injected_r": injected_r,
                "truth_pdf_index": truth_index,
                "truth_pdf_label": truth_label,
                "target_pdf_index": None if target_key[0] == 1 else target_key[1],
                "target_pdf_label": target_key[2],
                "fit_pdf_mode": "floating" if target_key[0] == 1 else "fixed",
                "metric": metric,
                "metric_name": metric_name,
                "sigma": sigma,
                "fit_failed": fit_failed,
                "experiment_index": experiment.get("experiment_index"),
                "job_id": experiment.get("job_id"),
                "status": experiment.get("status"),
                "fit_summary_html": experiment.get("fit_summary_html"),
            }
        )

    color_scale_min = -1.0
    color_scale_max = 1.0
    from matplotlib.colors import ListedColormap

    plot_matrix = matrix.T
    plot_sigma_matrix = sigma_matrix.T
    plot_failed_matrix = failed_matrix.T
    masked_matrix = np.ma.masked_invalid(plot_matrix)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e5e7eb")
    norm = TwoSlopeNorm(vmin=color_scale_min, vcenter=0.0, vmax=color_scale_max)
    failed_color = "#000000"
    failed_cmap = ListedColormap([failed_color])
    failed_cmap.set_bad((1.0, 1.0, 1.0, 0.0))

    x_fontsize = 8 if len(row_axis) <= 40 else 6 if len(row_axis) <= 120 else 4
    y_fontsize = 10 if len(target_axis) <= 10 else 8
    fig_width = max(14.0, min(200.0, 0.52 * len(row_axis) + 8.0))
    fig_height = max(6.0, 0.80 * len(target_axis) + 4.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#f8fafc")
    image = ax.imshow(masked_matrix, cmap=cmap, norm=norm, aspect="auto", interpolation="none")
    if np.any(plot_failed_matrix):
        failed_overlay = np.ma.masked_where(
            ~plot_failed_matrix,
            np.ones(plot_failed_matrix.shape, dtype=float),
        )
        ax.imshow(failed_overlay, cmap=failed_cmap, aspect="auto", interpolation="none", zorder=3)

    column_labels: list[str] = []
    group_boundaries: list[int] = []
    groups: list[tuple[int, int, str]] = []
    previous_group = None
    group_start = 0
    for column_index, row_key in enumerate(row_axis):
        group_key = row_key[:5]
        group_label = summary_heatmap_group_label(row_key)
        truth_label = row_key[6]
        if previous_group is None or group_key != previous_group:
            if previous_group is not None:
                group_boundaries.append(column_index)
                groups.append((group_start, column_index - 1, summary_heatmap_group_label(row_axis[group_start])))
                group_start = column_index
            column_labels.append(f"{group_label} | {truth_label}")
            previous_group = group_key
        else:
            column_labels.append(truth_label)
    groups.append((group_start, len(row_axis) - 1, summary_heatmap_group_label(row_axis[group_start])))

    ax.set_xticks(range(len(row_axis)))
    ax.set_xticklabels(
        column_labels,
        rotation=90,
        ha="right",
        va="center",
        rotation_mode="anchor",
        fontsize=x_fontsize,
    )
    ax.set_yticks(range(len(target_axis)))
    ax.set_yticklabels([axis_key[2] for axis_key in target_axis], fontsize=y_fontsize)
    ax.set_xlabel("Grouped tested combination | truth PDF", fontsize=12)
    ax.set_ylabel("Target PDF", fontsize=12)
    ax.set_title(label, fontsize=13, pad=14)
    ax.set_xticks(np.arange(-0.5, len(row_axis), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(target_axis), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)
    for boundary in group_boundaries:
        ax.axvline(boundary - 0.5, color="#0f172a", linewidth=1.1, alpha=0.85)
    if len(groups) <= 40:
        for start, end, group_label in groups:
            midpoint = 0.5 * (start + end)
            ax.text(
                midpoint,
                1.01,
                group_label,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=max(5, x_fontsize),
                color="#334155",
                clip_on=False,
            )

    annotation_fontsize = 8 if len(row_axis) <= 80 else 6 if len(row_axis) <= 180 else 4
    annotate_cells = len(row_axis) * len(target_axis) <= 4500
    column_is_asimov = [row_key[0] == "asimov" for row_key in row_axis]
    has_asimov_columns = any(column_is_asimov)
    has_annotated_columns = any(not is_asimov for is_asimov in column_is_asimov)
    if annotate_cells:
        for row in range(len(target_axis)):
            for column in range(len(row_axis)):
                if column_is_asimov[column] or bool(plot_failed_matrix[row, column]):
                    continue
                metric_value = plot_matrix[row, column]
                sigma_value = plot_sigma_matrix[row, column]
                if np.isfinite(sigma_value):
                    text = f"{float(sigma_value):.2f}"
                else:
                    text = "-"
                text_color = "white" if np.isfinite(metric_value) and abs(float(metric_value)) > 0.55 * color_scale_max else "#111827"
                ax.text(column, row, text, ha="center", va="center", color=text_color, fontsize=annotation_fontsize, fontweight="bold")

    colorbar_label = "toy fitted pull mean / Asimov closure pull"
    colorbar = fig.colorbar(image, ax=ax, shrink=0.92, extend="both")
    colorbar.set_label(colorbar_label, fontsize=11)
    colorbar.set_ticks([color_scale_min, -0.5, 0.0, 0.5, color_scale_max])
    cell_text_note = None
    if annotate_cells and has_annotated_columns:
        cell_text_note = "Cell text: fitted sigma for toy columns only." if has_asimov_columns else "Cell text: fitted sigma."
    elif has_annotated_columns:
        cell_text_note = "Cell text omitted because the summary heatmap is dense; see values JSON for fitted sigma."
    if cell_text_note is not None and np.any(plot_failed_matrix):
        cell_text_note += f" Failed cells are {failed_color}."
    if cell_text_note is not None:
        fig.text(0.01, 0.01, cell_text_note, ha="left", va="bottom", fontsize=9, color="#475569")
    fig.savefig(png_path, dpi=180)
    fig.savefig(pdf_path)
    plt.close(fig)

    if not has_annotated_columns:
        annotation = "none"
    elif annotate_cells and has_asimov_columns:
        annotation = "fitted sigma for non-Asimov columns"
    elif annotate_cells:
        annotation = "fitted sigma"
    else:
        annotation = "see values_json"

    payload = {
        "status": "ok",
        "label": label,
        "split": {"summary": "all_tested_combinations"},
        "plot_stem": plot_stem,
        "points": len([record for record in cell_records if record.get("metric") is not None]),
        "cells": len(cell_records),
        "metric": colorbar_label,
        "annotation": annotation,
        "color_scale": {"min": color_scale_min, "max": color_scale_max},
        "failed_cell_color": failed_color,
        "plot_orientation": "target_pdf_rows_tested_combination_columns",
        "row_axis": [
            {
                "dataset_strategy": row_key[0],
                "scheme": row_key[1],
                "target_poi": row_key[2],
                "target_label": row_key[3],
                "injected_r": row_key[4],
                "truth_pdf_index": row_key[5],
                "truth_pdf_label": row_key[6],
            }
            for row_key in row_axis
        ],
        "target_axis": [
            {
                "mode": "floating" if mode_rank == 1 else "fixed",
                "index": None if mode_rank == 1 else pdf_index,
                "label": axis_label,
            }
            for mode_rank, pdf_index, axis_label in target_axis
        ],
        "values": cell_records,
        "plot_png": str(png_path.resolve()),
        "plot_png_repo_relative": repo_relative(png_path),
        "plot_pdf": str(pdf_path.resolve()),
        "plot_pdf_repo_relative": repo_relative(pdf_path),
        "values_json": str(values_path.resolve()),
        "values_json_repo_relative": repo_relative(values_path),
        "note": f"Single summary heatmap over all tested combinations. Columns are grouped by dataset strategy, POI scheme, target POI, injected r, and truth PDF; rows show target PDFs. Cell color uses a fixed -1 to 1 scale for toy fitted pull mean or Asimov closure pull. Failed cells are filled {failed_color} and are not part of the colorbar. Asimov columns have no cell text.",
    }
    values_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def make_pdf_target_heatmap(
    experiments: list[dict[str, Any]],
    run_dir: Path,
    plot_stem: str,
    label: str,
    split: dict[str, Any],
    plt: Any,
    np: Any,
    TwoSlopeNorm: Any,
) -> dict[str, Any]:
    plot_dir = run_dir / "plots" / "pdf_target_heatmaps"
    plot_dir.mkdir(parents=True, exist_ok=True)
    png_path = plot_dir / f"{plot_stem}.png"
    pdf_path = plot_dir / f"{plot_stem}.pdf"
    values_path = plot_dir / f"{plot_stem}_values.json"

    truth_axis = sorted({truth_pdf_axis_key(experiment) for experiment in experiments})
    target_axis = sorted({target_pdf_axis_key(experiment) for experiment in experiments})
    if not truth_axis or not target_axis:
        payload = {
            "status": "failed",
            "message": "No truth or target PDF axis entries available.",
            "label": label,
            "split": split,
            "plot_stem": plot_stem,
            "points": 0,
        }
        values_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload

    truth_positions = {axis_key: index for index, axis_key in enumerate(truth_axis)}
    target_positions = {axis_key: index for index, axis_key in enumerate(target_axis)}
    matrix = np.full((len(truth_axis), len(target_axis)), np.nan, dtype=float)
    sigma_matrix = np.full((len(truth_axis), len(target_axis)), np.nan, dtype=float)
    failed_matrix = np.full((len(truth_axis), len(target_axis)), False, dtype=bool)
    cell_records: list[dict[str, Any]] = []

    for experiment in sorted(
        experiments,
        key=lambda item: (
            truth_pdf_axis_key(item),
            target_pdf_axis_key(item),
            str(item.get("job_id", "")),
        ),
    ):
        truth_key = truth_pdf_axis_key(experiment)
        target_key = target_pdf_axis_key(experiment)
        row = truth_positions[truth_key]
        column = target_positions[target_key]
        metric, sigma, metric_name = heatmap_metric_values(experiment)
        fit_failed = experiment_fit_failed(experiment)
        failed_matrix[row, column] = failed_matrix[row, column] or fit_failed
        if not fit_failed and metric is not None:
            matrix[row, column] = float(metric)
        if not fit_failed and sigma is not None:
            sigma_matrix[row, column] = abs(float(sigma))
        cell_records.append(
            {
                "truth_pdf_index": truth_key[0],
                "truth_pdf_label": truth_key[1],
                "target_pdf_index": None if target_key[0] == 1 else target_key[1],
                "target_pdf_label": target_key[2],
                "fit_pdf_mode": "floating" if target_key[0] == 1 else "fixed",
                "metric": metric,
                "metric_name": metric_name,
                "sigma": sigma,
                "fit_failed": fit_failed,
                "experiment_index": experiment.get("experiment_index"),
                "job_id": experiment.get("job_id"),
                "status": experiment.get("status"),
                "fit_summary_html": experiment.get("fit_summary_html"),
            }
        )

    color_scale_min = -1.0
    color_scale_max = 1.0
    masked_matrix = np.ma.masked_invalid(matrix)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e5e7eb")
    norm = TwoSlopeNorm(vmin=color_scale_min, vcenter=0.0, vmax=color_scale_max)

    fig_width = max(8.5, 1.6 * len(target_axis) + 5.0)
    fig_height = max(5.5, 0.80 * len(truth_axis) + 3.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#f8fafc")
    image = ax.imshow(masked_matrix, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(range(len(target_axis)))
    ax.set_xticklabels([axis_key[2] for axis_key in target_axis], rotation=35, ha="right", fontsize=10)
    ax.set_yticks(range(len(truth_axis)))
    ax.set_yticklabels([axis_key[1] for axis_key in truth_axis], fontsize=10)
    ax.set_xlabel("Target PDF", fontsize=12)
    ax.set_ylabel("Truth PDF", fontsize=12)
    ax.set_title(label, fontsize=13, pad=14)
    ax.set_xticks(np.arange(-0.5, len(target_axis), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(truth_axis), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    dataset_strategy = str(split.get("dataset_strategy", "toys"))
    annotate_cells = dataset_strategy != "asimov"
    if annotate_cells:
        for row in range(len(truth_axis)):
            for column in range(len(target_axis)):
                metric_value = matrix[row, column]
                sigma_value = sigma_matrix[row, column]
                fit_failed = bool(failed_matrix[row, column])
                if fit_failed:
                    text = "FAILED"
                elif np.isfinite(sigma_value):
                    text = f"{float(sigma_value):.2f}"
                else:
                    text = "-"
                text_color = "#991b1b" if fit_failed else "white" if np.isfinite(metric_value) and abs(float(metric_value)) > 0.55 * color_scale_max else "#111827"
                ax.text(column, row, text, ha="center", va="center", color=text_color, fontsize=10, fontweight="bold")

    colorbar_label = "Asimov closure pull" if dataset_strategy == "asimov" else "fitted pull mean"
    colorbar = fig.colorbar(image, ax=ax, shrink=0.92, extend="both")
    colorbar.set_label(colorbar_label, fontsize=11)
    colorbar.set_ticks([color_scale_min, -0.5, 0.0, 0.5, color_scale_max])
    if annotate_cells:
        fig.text(0.01, 0.01, "Cell text: fitted sigma", ha="left", va="bottom", fontsize=9, color="#475569")
    fig.savefig(png_path, dpi=180)
    fig.savefig(pdf_path)
    plt.close(fig)

    note = (
        "Rows show truth PDFs, columns show target PDFs. Cell color is Asimov closure pull on a fixed -1 to 1 scale. Asimov heatmaps omit all cell text."
        if dataset_strategy == "asimov"
        else "Rows show truth PDFs, columns show target PDFs. Cell color is the fitted pull mean on a fixed -1 to 1 scale. Cell text is the fitted sigma rounded to two decimals; failed fits are marked FAILED."
    )

    payload = {
        "status": "ok",
        "label": label,
        "split": split,
        "plot_stem": plot_stem,
        "points": len([record for record in cell_records if record.get("metric") is not None]),
        "cells": len(cell_records),
        "metric": colorbar_label,
        "annotation": "none" if dataset_strategy == "asimov" else "fitted sigma",
        "color_scale": {"min": color_scale_min, "max": color_scale_max},
        "truth_axis": [{"index": index, "label": axis_label} for index, axis_label in truth_axis],
        "target_axis": [
            {
                "mode": "floating" if mode_rank == 1 else "fixed",
                "index": None if mode_rank == 1 else pdf_index,
                "label": axis_label,
            }
            for mode_rank, pdf_index, axis_label in target_axis
        ],
        "values": cell_records,
        "plot_png": str(png_path.resolve()),
        "plot_png_repo_relative": repo_relative(png_path),
        "plot_pdf": str(pdf_path.resolve()),
        "plot_pdf_repo_relative": repo_relative(pdf_path),
        "values_json": str(values_path.resolve()),
        "values_json_repo_relative": repo_relative(values_path),
        "note": note,
    }
    values_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def make_summary_heatmaps(
    experiments: list[dict[str, Any]],
    run_dir: Path,
) -> list[dict[str, Any]]:
    plot_dir = run_dir / "plots" / "pdf_target_heatmaps"
    plot_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.colors import TwoSlopeNorm
    except Exception as exc:
        return [
            {
                "status": "failed",
                "message": f"Could not import matplotlib/numpy for PDF target heatmaps: {exc}",
                "label": "PDF target heatmaps",
                "split": {},
                "plot_stem": "pdf_target_heatmaps",
                "points": 0,
            }
        ]

    heatmap_tasks: list[tuple[list[dict[str, Any]], str, str, dict[str, Any]]] = []
    dataset_strategies = sorted(
        {str(experiment.get("dataset_strategy", "toys")) for experiment in experiments}
    )
    schemes = sorted(
        {str(experiment.get("scheme")) for experiment in experiments if experiment.get("scheme") is not None}
    )
    for dataset_strategy in dataset_strategies:
        for scheme in schemes:
            target_pois = sorted(
                {
                    str(experiment.get("target_poi"))
                    for experiment in experiments
                    if str(experiment.get("dataset_strategy", "toys")) == dataset_strategy
                    and str(experiment.get("scheme")) == scheme
                    and experiment.get("target_poi") is not None
                }
            )
            injected_rs = sorted(
                {
                    value
                    for experiment in experiments
                    if str(experiment.get("dataset_strategy", "toys")) == dataset_strategy
                    and str(experiment.get("scheme")) == scheme
                    for value in [finite_float_or_none(experiment.get("injected_r"))]
                    if value is not None
                }
            )
            for target_poi in target_pois:
                target_label = next(
                    (
                        str(experiment.get("target_label"))
                        for experiment in experiments
                        if str(experiment.get("dataset_strategy", "toys")) == dataset_strategy
                        and str(experiment.get("scheme")) == scheme
                        and str(experiment.get("target_poi")) == target_poi
                        and experiment.get("target_label") is not None
                    ),
                    target_poi,
                )
                scheme_display = poi_scheme_display(scheme)
                target_display = signal_target_display(target_label)
                dataset_display = dataset_strategy_display(dataset_strategy)
                for injected_r in injected_rs:
                    selected = [
                        experiment
                        for experiment in experiments
                        if str(experiment.get("dataset_strategy", "toys")) == dataset_strategy
                        and str(experiment.get("scheme")) == scheme
                        and str(experiment.get("target_poi")) == target_poi
                        and same_injected_r(experiment, injected_r)
                    ]
                    if not selected:
                        continue
                    label = (
                        f"{dataset_display.text} | {scheme_display.text} | {target_display.text} | "
                        f"r={format_number(injected_r)}"
                    )
                    plot_stem = (
                        f"pdf_target_heatmap__{dataset_display.slug}"
                        f"__scheme_{scheme_display.slug}"
                        f"__target_{target_display.slug}"
                        f"__{injection_output_slug(injected_r)}"
                    )
                    heatmap_tasks.append(
                        (
                            selected,
                            plot_stem,
                            label,
                            {
                                "dataset_strategy": dataset_strategy,
                                "scheme": scheme,
                                "target_poi": target_poi,
                                "target_label": target_label,
                                "injected_r": injected_r,
                            },
                        )
                    )

    total_heatmaps = len(heatmap_tasks) + 1
    log(
        f"Building {total_heatmaps} PDF target heatmap(s) from "
        f"{len(experiments)} fit summaries"
    )
    heatmaps: list[dict[str, Any]] = []
    log(
        f"PDF target heatmap 1/{total_heatmaps}: All Tested Combinations "
        f"({len(experiments)} fit summaries)"
    )
    summary_heatmap = make_pdf_target_summary_heatmap(
        experiments,
        run_dir,
        plt,
        np,
        TwoSlopeNorm,
    )
    heatmaps.append(summary_heatmap)
    log(
        f"Finished PDF target heatmap 1/{total_heatmaps}: All Tested Combinations "
        f"status={summary_heatmap.get('status', '-')}"
    )
    for index, (selected, plot_stem, label, split) in enumerate(heatmap_tasks, start=2):
        log(
            f"PDF target heatmap {index}/{total_heatmaps}: {label} "
            f"({len(selected)} fit summaries)"
        )
        heatmap = make_pdf_target_heatmap(
            selected,
            run_dir,
            plot_stem,
            label,
            split,
            plt,
            np,
            TwoSlopeNorm,
        )
        heatmaps.append(heatmap)
        log(
            f"Finished PDF target heatmap {index}/{total_heatmaps}: {label} "
            f"status={heatmap.get('status', '-')}"
        )
    return heatmaps


def failed_summary_heatmaps(
    run_dir: Path,
    message: str,
    returncode: int | None = None,
) -> list[dict[str, Any]]:
    plot_dir = run_dir / "plots" / "pdf_target_heatmaps"
    plot_dir.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "status": "failed",
        "label": "PDF target heatmaps",
        "split": {},
        "plot_stem": "pdf_target_heatmaps",
        "points": 0,
        "message": message,
        "stdout_path": str((plot_dir / "summary_heatmaps_stdout.txt").resolve()),
        "stdout_path_repo_relative": repo_relative(plot_dir / "summary_heatmaps_stdout.txt"),
        "stderr_path": str((plot_dir / "summary_heatmaps_stderr.txt").resolve()),
        "stderr_path_repo_relative": repo_relative(plot_dir / "summary_heatmaps_stderr.txt"),
    }
    if returncode is not None:
        payload["returncode"] = returncode
    return [payload]


def make_summary_heatmaps_isolated(
    experiments: list[dict[str, Any]],
    run_dir: Path,
) -> list[dict[str, Any]]:
    plot_dir = run_dir / "plots" / "pdf_target_heatmaps"
    plot_dir.mkdir(parents=True, exist_ok=True)
    input_path = plot_dir / "summary_heatmaps_input.json"
    output_path = plot_dir / "summary_heatmaps_output.json"
    stdout_path = plot_dir / "summary_heatmaps_stdout.txt"
    stderr_path = plot_dir / "summary_heatmaps_stderr.txt"
    input_payload = {
        "experiments": experiments,
        "run_dir": str(run_dir.resolve()),
    }
    input_path.write_text(json.dumps(input_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    code = r"""
import json
import sys
from pathlib import Path
import scripts.bias_study as bias_study

input_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
payload = json.loads(input_path.read_text(encoding='utf-8'))
heatmaps = bias_study.make_summary_heatmaps(
    payload['experiments'],
    Path(payload['run_dir']),
)
output_path.write_text(
    json.dumps({'heatmaps': heatmaps}, indent=2, sort_keys=True) + '\n',
    encoding='utf-8',
)
"""
    completed = run_streamed_plot_subprocess(
        [sys.executable, "-u", "-c", code, str(input_path), str(output_path)],
        REPO_ROOT,
        stdout_path,
        stderr_path,
        HEATMAP_PLOT_SUBPROCESS_TIMEOUT_SECONDS,
        "PDF target heatmaps",
    )
    if completed["timed_out"]:
        return failed_summary_heatmaps(
            run_dir,
            f"Summary heatmap plotting timed out after {HEATMAP_PLOT_SUBPROCESS_TIMEOUT_SECONDS} seconds.",
        )

    if int(completed["returncode"]) != 0:
        message = f"Summary heatmap plotting subprocess failed with return code {completed['returncode']}."
        if int(completed["returncode"]) < 0:
            message += f" Signal {-int(completed['returncode'])} terminated the subprocess."
        return failed_summary_heatmaps(run_dir, message, int(completed["returncode"]))
    if not output_path.exists():
        return failed_summary_heatmaps(
            run_dir,
            "Summary heatmap plotting subprocess finished without writing its output JSON.",
            int(completed["returncode"]),
        )

    try:
        output_payload = json.loads(output_path.read_text(encoding="utf-8"))
        return list(output_payload.get("heatmaps", []))
    except Exception as exc:
        return failed_summary_heatmaps(
            run_dir,
            f"Could not read summary heatmap plotting output JSON: {exc}",
            int(completed["returncode"]),
        )


def resolve_plot_only_run_dir(args: argparse.Namespace, output_dir: Path) -> Path:
    plot_only_value = args.plot_only
    if plot_only_value:
        requested = Path(str(plot_only_value))
        if requested.is_absolute():
            return requested.resolve()
        candidates = [output_dir / requested, REPO_ROOT / requested]
        for candidate in candidates:
            if candidate.exists():
                return candidate.resolve()
        return candidates[0].resolve()

    if args.run_name:
        return (output_dir / str(args.run_name)).resolve()

    run_dirs = [path for path in output_dir.iterdir() if path.is_dir()]
    if not run_dirs:
        raise RuntimeError(f"No bias-study run directories found in {output_dir}")
    return max(run_dirs, key=lambda path: (path.stat().st_mtime, path.name)).resolve()


def load_existing_run_results(run_dir: Path) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for summary_path in sorted(run_dir.rglob("summary.json")):
        if summary_path.parent == run_dir:
            continue
        try:
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if isinstance(payload, dict) and payload.get("job_id") and payload.get("method"):
            results.append(payload)
    return results


def first_metadata_value(results: list[dict[str, Any]], key: str, default: Any = None) -> Any:
    for result in results:
        value = (result.get("metadata") or {}).get(key)
        if value is not None:
            return value
    return default


def ordered_unique(values: list[Any]) -> list[Any]:
    unique: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        if value in seen:
            continue
        unique.append(value)
        seen.add(value)
    return unique


def infer_plot_only_args(
    args: argparse.Namespace,
    run_dir: Path,
    results: list[dict[str, Any]],
    existing_summary: dict[str, Any] | None,
) -> None:
    existing_summary = existing_summary or {}
    fit_results = [result for result in results if result.get("method") == "MultiDimFit"]
    dataset_strategies = ordered_unique(
        [
            str((result.get("metadata") or {}).get("dataset_strategy", "toys"))
            for result in fit_results
        ]
    ) or list(existing_summary.get("dataset_strategies", [])) or [str(args.dataset_strategy)]
    args.dataset_strategies = tuple(dataset_strategies)
    args.dataset_strategy = dataset_strategies[0]

    injected_values = sorted(
        {
            float(value)
            for result in fit_results
            for value in [(result.get("metadata") or {}).get("injected_r")]
            if value is not None
        }
    )
    if injected_values:
        args.injections = tuple(injected_values)
    elif existing_summary.get("injections"):
        args.injections = tuple(float(value) for value in existing_summary["injections"])

    args.mass = str(first_metadata_value(results, "mass", existing_summary.get("mass", args.mass)))
    args.toys = int(first_metadata_value(fit_results, "toys", existing_summary.get("toys", args.toys)))
    args.injection_mode = str(
        first_metadata_value(fit_results, "injection_mode", existing_summary.get("injection_mode", args.injection_mode))
    )
    args.pdf_index_name = str(
        first_metadata_value(fit_results, "pdf_index_name", existing_summary.get("pdf_index_name", args.pdf_index_name))
    )
    fit_pdf_modes = ordered_unique(
        [
            normalize_pdf_target_strategy((result.get("metadata") or {}).get("fit_pdf_mode", "floating"))
            for result in fit_results
        ]
    )
    if set(fit_pdf_modes) >= {"fixed", "floating"}:
        args.pdf_target_strategy = "both"
    elif fit_pdf_modes:
        args.pdf_target_strategy = fit_pdf_modes[0]
    else:
        args.pdf_target_strategy = normalize_pdf_target_strategy(
            existing_summary.get("pdf_target_strategy", getattr(args, "pdf_target_strategy", DEFAULT_PDF_TARGET_STRATEGY))
        )
    args.pdf_target_strategies = selected_pdf_target_strategies(args)
    args.poi_min = float(
        first_metadata_value(fit_results, "requested_poi_min", existing_summary.get("requested_poi_min", args.poi_min))
    )
    args.poi_max = float(
        first_metadata_value(fit_results, "requested_poi_max", existing_summary.get("requested_poi_max", args.poi_max))
    )
    args.quick = bool(first_metadata_value(fit_results, "quick", existing_summary.get("quick", args.quick)))
    args.cmin_default_minimizer_strategy = int(
        first_metadata_value(
            fit_results,
            "cmin_default_minimizer_strategy",
            existing_summary.get("cmin_default_minimizer_strategy", args.cmin_default_minimizer_strategy),
        )
    )
    args.robust_fit = bool(
        first_metadata_value(fit_results, "robust_fit", existing_summary.get("robust_fit", args.robust_fit))
    )
    args.bias_pull_threshold = float(existing_summary.get("bias_pull_threshold", args.bias_pull_threshold))
    args.pull_width_threshold = float(existing_summary.get("pull_width_threshold", args.pull_width_threshold))
    args.annotate_scatter_outliers = bool(
        existing_summary.get("annotate_scatter_outliers", args.annotate_scatter_outliers)
    )

    datacard = existing_summary.get("datacard") or args.datacard
    args.datacard = str(Path(datacard).resolve())
    args.run_name = run_dir.name


def infer_processes_for_plot_only(
    args: argparse.Namespace,
    results: list[dict[str, Any]],
    existing_summary: dict[str, Any] | None,
) -> DatacardProcesses:
    existing_summary = existing_summary or {}
    datacard_path = Path(args.datacard)
    if datacard_path.exists():
        try:
            return parse_datacard_processes(datacard_path)
        except Exception:
            pass

    signals = list(existing_summary.get("signals", []))
    backgrounds = list(existing_summary.get("backgrounds", []))
    if not signals:
        for result in results:
            metadata = result.get("metadata") or {}
            for process_name in metadata.get("target_processes", []) or []:
                signals.append(str(process_name))
            for poi_name in metadata.get("all_pois", []) or []:
                poi_text = str(poi_name)
                if poi_text.startswith("r_") and not poi_text.endswith("_grouped"):
                    signals.append(poi_text[2:])
    return DatacardProcesses(signals=ordered_unique(signals), backgrounds=ordered_unique(backgrounds))


def infer_schemes_for_plot_only(
    results: list[dict[str, Any]],
    existing_summary: dict[str, Any] | None,
) -> tuple[PoiScheme, ...]:
    existing_summary = existing_summary or {}
    workspace_metadata_by_scheme = {
        str((result.get("metadata") or {}).get("scheme")): result.get("metadata") or {}
        for result in results
        if result.get("method") == "text2workspace.py" and (result.get("metadata") or {}).get("scheme")
    }
    existing_schemes = {
        str(scheme.get("name")): scheme
        for scheme in existing_summary.get("schemes", [])
        if isinstance(scheme, dict) and scheme.get("name")
    }
    scheme_names = ordered_unique(
        [
            str((result.get("metadata") or {}).get("scheme"))
            for result in results
            if (result.get("metadata") or {}).get("scheme")
        ]
    )
    schemes: list[PoiScheme] = []
    for scheme_name in scheme_names:
        workspace_metadata = workspace_metadata_by_scheme.get(scheme_name, {})
        existing_scheme = existing_schemes.get(scheme_name, {})
        all_pois = list(workspace_metadata.get("pois") or existing_scheme.get("pois") or [])
        poi_maps = tuple(workspace_metadata.get("poi_maps") or existing_scheme.get("poi_maps") or [])
        targets: list[PoiTarget] = []
        for result in results:
            metadata = result.get("metadata") or {}
            if metadata.get("scheme") != scheme_name or not metadata.get("target_poi"):
                continue
            target = PoiTarget(
                label=str(metadata.get("target_label", metadata["target_poi"])),
                poi=str(metadata["target_poi"]),
                processes=tuple(str(process) for process in metadata.get("target_processes", []) or []),
            )
            if target not in targets:
                targets.append(target)
            for poi_name in metadata.get("all_pois", []) or []:
                all_pois.append(str(poi_name))
        if not all_pois:
            all_pois = [target.poi for target in targets]
        schemes.append(
            PoiScheme(
                name=scheme_name,
                title=str(workspace_metadata.get("scheme_title", existing_scheme.get("title", scheme_name))),
                description=str(
                    workspace_metadata.get(
                        "scheme_description",
                        existing_scheme.get("description", "Reconstructed from existing fit summaries."),
                    )
                ),
                poi_maps=poi_maps,
                targets=tuple(targets),
                all_poi_names=tuple(ordered_unique(all_pois)),
            )
        )
    return tuple(schemes)


def infer_truth_pdf_metadata_for_plot_only(
    results: list[dict[str, Any]],
    existing_summary: dict[str, Any] | None,
) -> dict[str, Any]:
    existing_summary = existing_summary or {}
    if existing_summary.get("truth_pdf_metadata"):
        return dict(existing_summary["truth_pdf_metadata"])
    truths: dict[int, dict[str, Any]] = {}
    for result in results:
        metadata = result.get("metadata") or {}
        if metadata.get("truth_pdf_index") is None:
            continue
        index = int(metadata["truth_pdf_index"])
        display = pdf_state_display(
            index=index,
            family=metadata.get("truth_pdf_family"),
            order=metadata.get("truth_pdf_scan_order"),
            selection_role=metadata.get("truth_pdf_selection_role"),
            name=metadata.get("truth_pdf"),
        )
        truths[index] = {
            "index": index,
            "name": str(metadata.get("truth_pdf", "")),
            "family": str(metadata.get("truth_pdf_family", "")),
            "label": str(metadata.get("truth_pdf_label") or display.text),
            "slug": str(metadata.get("truth_pdf_slug") or display.slug),
            "latex": str(metadata.get("truth_pdf_latex") or display.latex),
            "scan_order": metadata.get("truth_pdf_scan_order"),
            "selection_role": metadata.get("truth_pdf_selection_role"),
        }
    return {"selected_truth_pdfs": [truths[index] for index in sorted(truths)]}


def infer_scheme_injection_map_for_plot_only(
    args: argparse.Namespace,
    schemes: tuple[PoiScheme, ...],
    fit_results: list[dict[str, Any]],
    existing_summary: dict[str, Any] | None,
) -> dict[str, tuple[float, ...]]:
    existing_summary = existing_summary or {}
    existing_map = existing_summary.get("injections_by_scheme")
    if not isinstance(existing_map, dict):
        existing_map = {}

    mapping: dict[str, tuple[float, ...]] = {}
    for scheme in schemes:
        injected_values = sorted(
            {
                float(value)
                for result in fit_results
                for metadata in [result.get("metadata") or {}]
                for value in [metadata.get("injected_r")]
                if metadata.get("scheme") == scheme.name and value is not None
            }
        )
        if not injected_values and scheme.name in existing_map:
            injected_values = [float(value) for value in existing_map[scheme.name]]
        mapping[scheme.name] = tuple(injected_values or args.injections)
    return mapping


def run_plot_only(args: argparse.Namespace, output_dir: Path) -> int:
    run_dir = resolve_plot_only_run_dir(args, output_dir)
    if not run_dir.exists() or not run_dir.is_dir():
        raise FileNotFoundError(f"Plot-only run directory does not exist: {run_dir}")
    top_summary_path = run_dir / "summary.json"
    existing_summary = None
    if top_summary_path.exists():
        existing_summary = json.loads(top_summary_path.read_text(encoding="utf-8"))
    results = load_existing_run_results(run_dir)
    fit_results = [result for result in results if result.get("method") == "MultiDimFit"]
    if not fit_results:
        raise RuntimeError(f"No existing MultiDimFit job summaries found under {run_dir}")

    infer_plot_only_args(args, run_dir, results, existing_summary)
    processes = infer_processes_for_plot_only(args, results, existing_summary)
    schemes = infer_schemes_for_plot_only(results, existing_summary)
    args.scheme_injection_map = infer_scheme_injection_map_for_plot_only(
        args,
        schemes,
        fit_results,
        existing_summary,
    )
    truth_pdf_metadata = infer_truth_pdf_metadata_for_plot_only(results, existing_summary)
    analyzed_fit_results = [analyze_fit_result(result, args) for result in fit_results]
    experiments = [
        experiment_summary_from_fit(result, index)
        for index, result in enumerate(analyzed_fit_results, start=1)
    ]
    scatter_plot, scatter_plot_variations = make_summary_scatter_plots_isolated(
        experiments,
        run_dir,
        args.bias_pull_threshold,
        args.annotate_scatter_outliers,
    )
    heatmap_plots = make_summary_heatmaps_isolated(experiments, run_dir)
    write_run_summary(
        run_dir,
        output_dir,
        args,
        processes,
        schemes,
        truth_pdf_metadata,
        results,
        scatter_plot,
        scatter_plot_variations,
        heatmap_plots,
    )
    log(f"Rebuilt plot-only summary: {repo_relative(run_dir / 'summary.json')}")
    log(f"Global summary: {repo_relative(output_dir / 'bias_study.html')}")
    return 0


def prepare_run_dir(output_dir: Path, run_name: str | None) -> Path:
    if run_name is None:
        run_name = dt.datetime.now(dt.timezone.utc).strftime("run_%Y%m%d_%H%M%S")
    run_dir = output_dir / run_name
    if run_dir.exists():
        log(f"Clearing run directory before starting: {repo_relative(run_dir)}")
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def write_run_summary(
    run_dir: Path,
    output_dir: Path,
    args: argparse.Namespace,
    processes: DatacardProcesses,
    schemes: tuple[PoiScheme, ...],
    truth_pdf_metadata: dict[str, Any],
    results: list[dict[str, Any]],
    scatter_plot: dict[str, Any] | None,
    scatter_plot_variations: dict[str, list[dict[str, Any]]] | None = None,
    heatmap_plots: list[dict[str, Any]] | None = None,
) -> None:
    fit_results = [result for result in results if result.get("method") == "MultiDimFit"]
    experiments = [experiment_summary_from_fit(result, index) for index, result in enumerate(fit_results, start=1)]
    bias_flags = [experiment for experiment in experiments if experiment.get("bias_flag") is True]
    width_flags = [experiment for experiment in experiments if experiment.get("pull_width_flag") is True]
    scheme_injection_map = getattr(args, "scheme_injection_map", None) or injections_by_scheme(args, schemes)
    payload = {
        "schema_version": 1,
        "run_dir": str(run_dir.resolve()),
        "run_dir_repo_relative": repo_relative(run_dir),
        "summary_json": str((run_dir / "summary.json").resolve()),
        "summary_html": str((run_dir / "summary.html").resolve()),
        "datacard": str(Path(args.datacard).resolve()),
        "datacard_repo_relative": repo_relative(Path(args.datacard).resolve()),
        "mass": args.mass,
        "dataset_strategy": args.dataset_strategy,
        "dataset_strategies": list(args.dataset_strategies),
        "toys": args.toys,
        "nproc": nproc(),
        "workers_requested": args.workers,
        "workers_effective_cap": nproc() if args.workers <= 0 else min(args.workers, nproc()),
        "worker_policy": "Each dependency wave uses min(ready_jobs, requested_workers_or_nproc, nproc).",
        "injections": list(args.injections),
        "default_injections": list(args.injections),
        "injections_by_scheme": {
            scheme_name: list(injected_values)
            for scheme_name, injected_values in scheme_injection_map.items()
        },
        "injection_mode": args.injection_mode,
        "fit_method": "MultiDimFit",
        "fit_algo": FIT_ALGO,
        "fit_poi_scan_mode": "target-only",
        "float_other_pois": True,
        "pdf_target_strategy": args.pdf_target_strategy,
        "pdf_target_strategies": list(args.pdf_target_strategies),
        "fit_pdf_mode": ",".join(args.pdf_target_strategies),
        "fit_pdf_index": None,
        "pdf_index_name": args.pdf_index_name,
        "poi_range_policy": "fixed requested range",
        "requested_poi_min": args.poi_min,
        "requested_poi_max": args.poi_max,
        "poi_min": args.poi_min,
        "poi_max": args.poi_max,
        "quick": bool(args.quick),
        "cmin_default_minimizer_strategy": int(args.cmin_default_minimizer_strategy),
        "cmin_strategy_explicit": bool(args.cmin_strategy_explicit),
        "robust_fit": bool(args.robust_fit),
        "bias_pull_threshold": args.bias_pull_threshold,
        "pull_width_threshold": args.pull_width_threshold,
        "annotate_scatter_outliers": args.annotate_scatter_outliers,
        "signals": processes.signals,
        "backgrounds": processes.backgrounds,
        "schemes": [
            {
                "name": scheme.name,
                "slug": scheme_slug(scheme.name),
                "label": scheme_label(scheme.name),
                "latex": scheme_latex(scheme.name),
                "title": scheme.title,
                "description": scheme.description,
                "pois": list(scheme.all_pois),
                "poi_maps": list(scheme.poi_maps),
                "targets": [
                    {
                        "label": target.label,
                        "slug": target_slug(target.label),
                        "display_label": target_label(target.label),
                        "latex": target_latex(target.label),
                        "poi": target.poi,
                        "processes": list(target.processes),
                    }
                    for target in scheme.targets
                ],
            }
            for scheme in schemes
        ],
        "truth_pdf_metadata": truth_pdf_metadata,
        "jobs": [short_result(result) for result in results],
        "experiments": experiments,
        "scatter_plot": scatter_plot,
        "scatter_plot_variations": scatter_plot_variations or {},
        "heatmap_plots": heatmap_plots or [],
        "bias_flags": bias_flags,
        "pull_width_flags": width_flags,
        "answer": {
            "question": "Did the choice of non-resonant background introduce a non-negligible bias?",
            "non_negligible_bias_found": bool(bias_flags),
            "criterion": f"toys: abs(fitted normal pull mean) >= {format_number(args.bias_pull_threshold)}; asimov: abs(closure pull) >= {format_number(args.bias_pull_threshold)}",
            "flagged_experiments": len(bias_flags),
            "pull_width_issue_found": bool(width_flags),
            "pull_width_criterion": f"abs(fitted normal pull sigma - 1) >= {format_number(args.pull_width_threshold)}",
        },
    }
    (run_dir / "summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    summary_rows = [
        ["Run directory", repo_relative(run_dir)],
        ["Datacard", repo_relative(Path(args.datacard).resolve())],
        ["Mass", args.mass],
        ["Dataset strategy", args.dataset_strategy],
        ["Worker request", "auto" if args.workers <= 0 else args.workers],
        ["System nproc", nproc()],
        ["Worker cap", nproc() if args.workers <= 0 else min(args.workers, nproc())],
        ["POI schemes", ", ".join(scheme_label(scheme.name) for scheme in schemes)],
        ["Default injections", format_injection_values(args.injections)],
        ["Injections by scheme", format_injections_by_scheme(scheme_injection_map)],
        ["Injection mode", args.injection_mode],
        ["Requested POI range", f"[{format_number(args.poi_min)}, {format_number(args.poi_max)}]"],
        ["Job POI range policy", "fixed requested range"],
        ["Fit method", f"MultiDimFit --algo {FIT_ALGO}"],
        ["Fit POI scan", "target POI only (-P target) with other signal POIs floating"],
        ["PDF target strategy", args.pdf_target_strategy],
        ["Fit PDF modes", ", ".join(args.pdf_target_strategies)],
        ["Quick mode", "enabled" if args.quick else "disabled"],
        ["Minimizer strategy", args.cmin_default_minimizer_strategy],
        ["Minimizer strategy explicit", "yes" if args.cmin_strategy_explicit else "no"],
        ["Robust fit", "enabled" if args.robust_fit else "disabled"],
        ["Bias threshold", args.bias_pull_threshold],
        ["Scatter outlier annotations", "enabled" if args.annotate_scatter_outliers else "disabled"],
        ["Aggregate JSON", html_link("summary.json", "summary.json")],
    ]
    if args.dataset_strategy == "toys":
        summary_rows.insert(4, ["Toys per experiment", args.toys])
        summary_rows.insert(-2, ["Pull-width threshold", args.pull_width_threshold])
    answer_rows = [
        ["Non-negligible bias found", "yes" if bias_flags else "no"],
        ["Flagged bias experiments", len(bias_flags)],
    ]
    if args.dataset_strategy == "toys":
        answer_rows.extend(
            [
                ["Pull-width issue found", "yes" if width_flags else "no"],
                ["Flagged pull-width experiments", len(width_flags)],
            ]
        )
    scatter_section = ""
    if scatter_plot and scatter_plot.get("status") == "ok":
        scatter_href = href_relative_to(run_dir, Path(scatter_plot["plot_png"]))
        scatter_section = f"""
    <section class="card"><h2>Global Pull Scatter</h2>
      <p class="muted">{html_escape(scatter_plot.get('note', ''))}</p>
      <p><img src="{html_escape(scatter_href)}" alt="pull mean and sigma scatter"></p>
    </section>
        """
    elif scatter_plot:
        scatter_section = f"""
    <section class="card"><h2>Global Pull Scatter</h2>
      <p class="muted">Plotting failed: {html_escape(scatter_plot.get('message', '-'))}</p>
    </section>
        """
    variation_plots = scatter_plot_variations or {}
    scatter_variation_sections = "".join(
        [
            scatter_plot_grid_section(
                "Pull Scatter By POI Scheme",
                variation_plots.get("by_poi_scheme", []),
                run_dir,
            ),
            scatter_plot_grid_section(
                "Pull Scatter By Injected R",
                variation_plots.get("by_injected_r", []),
                run_dir,
            ),
            scatter_plot_grid_section(
                "Pull Scatter By POI Scheme And Injected R",
                variation_plots.get("by_poi_scheme_and_injected_r", []),
                run_dir,
            ),
        ]
    )
    heatmap_section = scatter_plot_grid_section(
        "PDF Target Heatmaps",
        heatmap_plots or [],
        run_dir,
    )
    experiment_headers = [
        "Index",
        "Scheme",
        "Target POI",
        "Injected r",
        "Truth PDF",
        "Target PDF",
        "POI Range",
        "Entries Used",
    ]
    if args.dataset_strategy == "asimov":
        experiment_headers.extend(["Closure Pull", "Bias", "Fit Job"])
    else:
        experiment_headers.extend(["Pull Mean", "Pull Sigma", "Bias", "Width", "Pull Plot", "Fit Job"])

    experiment_rows = []
    for experiment in experiments:
        pull_plot = experiment.get("pull_plot")
        fit_html = experiment.get("fit_summary_html")
        poi_range = html_escape(
            f"[{format_number(float(experiment['poi_min']))}, {format_number(float(experiment['poi_max']))}]"
            if experiment.get("poi_min") is not None and experiment.get("poi_max") is not None
            else "-"
        )
        row = [
            html_escape(experiment.get("experiment_index", "-")),
            html_escape(experiment.get("scheme_label") or experiment.get("scheme", "-")),
            html_escape(experiment.get("target_display_label") or experiment.get("target_poi", "-")),
            html_escape(experiment.get("injected_r", "-")),
            html_escape(experiment.get("truth_pdf_label", "-")),
            html_escape(fit_target_pdf_label(experiment)),
            poi_range,
            html_escape(experiment.get("entries_used", "-")),
        ]
        if args.dataset_strategy == "asimov":
            row.extend(
                [
                    html_fixed_number_or_dash(experiment.get("asimov_closure_pull")),
                    flag_pill(experiment.get("bias_flag")),
                    html_link("fit", href_relative_to(run_dir, REPO_ROOT / fit_html)) if fit_html else "-",
                ]
            )
        else:
            row.extend(
                [
                    html_fixed_number_or_dash(experiment.get("gaussian_mean")),
                    html_fixed_number_or_dash(experiment.get("gaussian_sigma")),
                    flag_pill(experiment.get("bias_flag")),
                    flag_pill(experiment.get("pull_width_flag")),
                    html_link("pull", href_relative_to(run_dir, REPO_ROOT / pull_plot)) if pull_plot else "-",
                    html_link("fit", href_relative_to(run_dir, REPO_ROOT / fit_html)) if fit_html else "-",
                ]
            )
        experiment_rows.append(row)
    job_rows = [
        [
            html_escape(result.get("job_id", "-")),
            html_escape(result.get("method", "-")),
            status_pill(result.get("status", "-")),
            html_escape(result.get("metadata", {}).get("dataset_strategy", "-")),
            html_escape(result.get("metadata", {}).get("scheme_label", result.get("metadata", {}).get("scheme", "-"))),
            html_escape(result.get("metadata", {}).get("target_display_label", result.get("metadata", {}).get("target_poi", "-"))),
            html_escape(result.get("metadata", {}).get("truth_pdf_label", "-")),
            html_escape(result.get("metadata", {}).get("target_pdf_label", result.get("metadata", {}).get("fit_pdf_label", "-"))),
            html_escape(result.get("metadata", {}).get("injected_r", "-")),
            html_escape(result.get("duration", "-")),
            html_link("summary", href_relative_to(run_dir, Path(result.get("summary_html_path", ""))))
            if result.get("summary_html_path")
            else "-",
        ]
        for result in results
    ]
    policy_note = """
      <p>Pseudo-dataset generation freezes <code>pdfindex</code> to the requested truth PDF and uses <code>--bypassFrequentistFit</code>, so generated toys and Asimov datasets are based on pre-fit model values and do not use observed data to determine nuisance values. In target-only injection mode, non-tested POIs are set to zero and frozen during generation.</p>
      <p>The toy strategy uses <code>-t N</code>; the Asimov strategy uses <code>-t -1</code>. Fitting uses <code>MultiDimFit --algo singles -P &lt;target POI&gt; --floatOtherPOIs 1</code> with all scheme POIs in <code>--redefineSignalPOIs</code>. The <code>--pdf-target-strategy</code> option controls whether <code>pdfindex</code> floats/profiled, is fixed to each target PDF index, or both. Job-level POI ranges use the fixed requested range. Toy pull denominators use the target POI profile endpoint rows from <code>--algo singles</code>, with <code>trackedError_&lt;poi&gt;</code> as a fallback if endpoints are unusable. Asimov fits report r closure and the closure pull, and do not run a pull-width check. By default, fits use <code>--robustFit 1</code> and <code>--cminDefaultMinimizerStrategy 2</code>; <code>--quick</code> disables robust fitting and uses strategy 0 unless <code>--cmin-strategy</code> is explicitly set.</p>
    """
    body = f"""
    <h1>Bias Study Run</h1>
    <p class="subtitle">Pseudo-dataset generation and MultiDimFit pull study for non-resonant-background model bias.</p>
    <section class="card"><h2>Summary</h2>{html_table_raw(['Field', 'Value'], [[html_escape(k), v if str(v).startswith('<a ') else html_escape(v)] for k, v in summary_rows])}</section>
    <section class="card"><h2>Answer</h2>{html_table(['Field', 'Value'], answer_rows)}</section>
    <section class="card"><h2>Combine Prescription</h2>{policy_note}</section>
    {scatter_section}
    {scatter_variation_sections}
    {heatmap_section}
    <section class="card"><h2>Experiments</h2>{html_table_raw(experiment_headers, experiment_rows)}</section>
    <section class="card"><h2>Jobs</h2>{html_table_raw(['Job', 'Method', 'Status', 'Strategy', 'Scheme', 'Target POI', 'Truth PDF', 'Target PDF', 'Injected r', 'Duration', 'Summary'], job_rows)}</section>
    """
    (run_dir / "summary.html").write_text(html_document("Bias Study Run", body), encoding="utf-8")
    write_global_summary(output_dir)


def write_global_summary(output_dir: Path) -> None:
    summaries: list[dict[str, Any]] = []
    for summary_path in sorted(output_dir.glob("*/summary.json"), reverse=True):
        try:
            summaries.append(json.loads(summary_path.read_text(encoding="utf-8")))
        except Exception:
            continue
    rows = []
    for summary in summaries:
        run_dir = Path(summary.get("run_dir", ""))
        summary_injections = summary.get("injections_by_scheme")
        if isinstance(summary_injections, dict) and summary_injections:
            injections_text = "; ".join(
                f"{scheme_label(scheme)}={','.join(str(value) for value in values)}"
                for scheme, values in summary_injections.items()
            )
        else:
            injections_text = ", ".join(str(value) for value in summary.get("injections", []))
        pdf_target_text = summary.get("pdf_target_strategy") or summary.get("fit_pdf_mode", "-")
        rows.append(
            [
                html_escape(summary.get("run_dir_repo_relative", run_dir.name)),
                html_escape(", ".join(summary.get("dataset_strategies", [summary.get("dataset_strategy", "toys")]))),
                html_escape(summary.get("toys", "-")),
                html_escape(injections_text),
                html_escape(pdf_target_text),
                html_escape(summary.get("answer", {}).get("flagged_experiments", "-")),
                "yes" if summary.get("answer", {}).get("non_negligible_bias_found") else "no",
                html_link("summary.html", href_relative_to(output_dir, run_dir / "summary.html")),
                html_link("summary.json", href_relative_to(output_dir, run_dir / "summary.json")),
            ]
        )
    body = f"""
    <h1>Bias Study</h1>
    <p class="subtitle">Global index of bias-study runs.</p>
    <section class="card"><h2>Runs</h2>{html_table_raw(['Run', 'Strategy', 'Toys/Exp', 'Injections', 'PDF Target Strategy', 'Flagged Bias Experiments', 'Bias Found', 'HTML', 'JSON'], rows)}</section>
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "bias_study.html").write_text(html_document("Bias Study", body), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run a Combine pseudo-dataset bias study for non-resonant-background PDF choices "
            "using GenerateOnly and MultiDimFit."
        )
    )
    parser.set_defaults(cmin_strategy_explicit=False)
    parser.add_argument(
        "--datacard",
        default=str(DEFAULT_DATACARD),
        help="Input text datacard produced by scripts/build_bundled_workspace.py.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Parent directory for bias-study runs.",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help="Run subdirectory name. Defaults to a UTC timestamp.",
    )
    parser.add_argument(
        "--plot-only",
        nargs="?",
        const="",
        default=None,
        metavar="RUN_DIR",
        help=(
            "Rebuild only the top-level summary and scatter plots from existing per-job summaries. "
            "Optionally pass a run directory; without a value, --run-name is used, or the latest run under --output-dir."
        ),
    )
    parser.add_argument("--mass", default=DEFAULT_MASS, help="Mass label passed to Combine.")
    parser.add_argument(
        "--poi-scheme",
        choices=("three_poi_z", "three_poi_h", "three_poi", "z_grouped", "h_grouped", "grouped", "both"),
        default="both",
        help=(
            "POI scheme to run. `three_poi` runs three_poi_z and three_poi_h; "
            "`grouped` runs z_grouped and h_grouped; the default `both` runs all four schemes."
        ),
    )
    parser.add_argument(
        "--injections",
        type=parse_float_list,
        default=DEFAULT_INJECTIONS,
        help="Default comma-separated injected signal strengths for schemes without --scheme-injections.",
    )
    parser.add_argument(
        "--scheme-injections",
        action="append",
        type=parse_scheme_injection_spec,
        default=None,
        metavar="SCHEME=R0,R1,...",
        help=(
            "Override injected signal strengths for one selected scheme. "
            "May be repeated, e.g. --scheme-injections z_grouped=0,1,10 "
            "--scheme-injections h_grouped=0,10000,50000."
        ),
    )
    parser.add_argument(
        "--toys",
        type=int,
        default=DEFAULT_TOYS,
        help="Number of random toys per truth-PDF/injection/target experiment. Ignored for Asimov-only jobs.",
    )
    parser.add_argument(
        "--dataset-strategy",
        choices=DATASET_STRATEGIES,
        default=DEFAULT_DATASET_STRATEGY,
        help="Pseudo-dataset strategy: random toys or an Asimov dataset.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=0,
        help="Parallel subprocess workers per dependency wave. Use 0 for auto; the effective value is capped at nproc.",
    )
    parser.add_argument(
        "--poi-initial",
        type=float,
        default=DEFAULT_POI_INITIAL,
        help="Initial value used when defining POIs in text2workspace.py.",
    )
    parser.add_argument(
        "--poi-min",
        type=float,
        default=DEFAULT_POI_MIN,
        help="Lower bound used when defining POIs and running GenerateOnly/MultiDimFit jobs.",
    )
    parser.add_argument(
        "--poi-max",
        type=float,
        default=DEFAULT_POI_MAX,
        help="Upper bound used when defining POIs and running GenerateOnly/MultiDimFit jobs.",
    )
    parser.add_argument(
        "--injection-mode",
        choices=("target-only", "all-pois"),
        default="target-only",
        help="Inject the requested r into only the target POI, or into all POIs in the scheme.",
    )
    parser.add_argument(
        "--pdf-families",
        type=parse_family_list,
        default=("all",),
        help="Comma-separated truth PDF families to generate from, e.g. all,johnson,bernstein,chebychev,power_law,exponential.",
    )
    parser.add_argument(
        "--pdf-index-name",
        default=DEFAULT_PDF_INDEX_NAME,
        help="Name of the RooMultiPdf discrete category in the Combine workspace.",
    )
    parser.add_argument(
        "--pdf-target-strategy",
        choices=PDF_TARGET_STRATEGIES,
        default=DEFAULT_PDF_TARGET_STRATEGY,
        help=(
            "How to handle the target PDF during MultiDimFit: `floating` keeps pdfindex profiled, "
            "`fixed` runs one fit per fixed target pdfindex, and `both` runs both strategies. "
            "Default is `floating`."
        ),
    )
    parser.add_argument(
        "--profile-freeze-disassociated-params",
        action="store_true",
        help="Pass --X-rtd MINIMIZER_freezeDisassociatedParams to MultiDimFit jobs.",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help=(
            "Use faster, less robust MultiDimFit settings: disable --robustFit and use "
            f"--cminDefaultMinimizerStrategy {QUICK_CMIN_DEFAULT_MINIMIZER_STRATEGY} "
            "unless --cmin-strategy is explicitly set."
        ),
    )
    parser.add_argument(
        "--cmin-strategy",
        "--cmin-default-minimizer-strategy",
        dest="cmin_default_minimizer_strategy",
        type=int,
        choices=(0, 1, 2),
        default=DEFAULT_CMIN_DEFAULT_MINIMIZER_STRATEGY,
        action=StoreCminStrategy,
        help=(
            "Value passed as --cminDefaultMinimizerStrategy to MultiDimFit. "
            f"Default is {DEFAULT_CMIN_DEFAULT_MINIMIZER_STRATEGY}; --quick overrides this to "
            f"{QUICK_CMIN_DEFAULT_MINIMIZER_STRATEGY} unless this option is explicitly set."
        ),
    )
    parser.add_argument(
        "--robust-fit",
        action="store_true",
        default=True,
        help="Pass --robustFit 1 to MultiDimFit. Enabled by default; --quick disables it.",
    )
    parser.add_argument(
        "--seed-base",
        type=int,
        default=DEFAULT_SEED_BASE,
        help="Base seed for toy GenerateOnly jobs. Use --random-seeds to pass -s -1 instead of deterministic seeds. Ignored for Asimov jobs.",
    )
    parser.add_argument(
        "--random-seeds",
        action="store_true",
        help="Use random Combine seeds (-s -1) for toy generation instead of deterministic seeds. Ignored for Asimov jobs.",
    )
    parser.add_argument(
        "--pull-range",
        type=parse_float_list,
        default=DEFAULT_PULL_RANGE,
        help="Comma-separated pull histogram range, e.g. -5,5.",
    )
    parser.add_argument("--pull-bins", type=int, default=DEFAULT_PULL_BINS)
    parser.add_argument(
        "--min-fit-entries",
        type=int,
        default=DEFAULT_MIN_FIT_ENTRIES,
        help="Minimum usable toy entries required before fitting the pull histogram to a normal function.",
    )
    parser.add_argument(
        "--bias-pull-threshold",
        type=float,
        default=DEFAULT_BIAS_PULL_THRESHOLD,
        help="Flag toy experiments with abs(fitted normal pull mean), or Asimov experiments with abs(closure pull), at or above this value.",
    )
    parser.add_argument(
        "--pull-width-threshold",
        type=float,
        default=DEFAULT_PULL_WIDTH_THRESHOLD,
        help="Flag experiments with abs(fitted normal pull sigma - 1) at or above this value.",
    )
    parser.add_argument(
        "--annotate-scatter-outliers",
        action="store_true",
        help="Annotate points in pull_mean_sigma_scatter above the bias threshold. Disabled by default.",
    )
    parser.add_argument(
        "--skip-fits",
        action="store_true",
        help="Only build workspaces and generate pseudo-datasets; do not run MultiDimFit.",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    args.dataset_strategies = selected_dataset_strategies(args)
    args.pdf_target_strategy = normalize_pdf_target_strategy(args.pdf_target_strategy)
    args.pdf_target_strategies = selected_pdf_target_strategies(args)
    if args.quick:
        args.robust_fit = False
        if not args.cmin_strategy_explicit:
            args.cmin_default_minimizer_strategy = QUICK_CMIN_DEFAULT_MINIMIZER_STRATEGY
    if "toys" in args.dataset_strategies and args.toys <= 0:
        raise ValueError("--toys must be positive when --dataset-strategy includes toys")
    if args.workers < 0:
        raise ValueError("--workers must be non-negative")
    if args.cmin_default_minimizer_strategy < 0:
        raise ValueError("--cmin-strategy must be non-negative")
    if not args.injections:
        raise ValueError("--injections must contain at least one value")
    if any(not math.isfinite(float(value)) for value in args.injections):
        raise ValueError("--injections must contain only finite values")
    seen_scheme_injections: set[str] = set()
    args.scheme_injections = tuple(args.scheme_injections or ())
    for scheme_name, injected_values in args.scheme_injections:
        if scheme_name in seen_scheme_injections:
            raise ValueError(f"--scheme-injections was provided more than once for {scheme_name}")
        seen_scheme_injections.add(scheme_name)
        if not injected_values:
            raise ValueError(f"--scheme-injections for {scheme_name} must contain at least one value")
        if any(not math.isfinite(float(value)) for value in injected_values):
            raise ValueError(f"--scheme-injections for {scheme_name} must contain only finite values")
    if not math.isfinite(float(args.poi_min)) or not math.isfinite(float(args.poi_max)):
        raise ValueError("--poi-min and --poi-max must be finite")
    if args.poi_max <= args.poi_min:
        raise ValueError("--poi-max must be greater than --poi-min")
    if len(args.pull_range) != 2 or args.pull_range[1] <= args.pull_range[0]:
        raise ValueError("--pull-range must contain two increasing values")
    if args.pull_bins <= 0:
        raise ValueError("--pull-bins must be positive")
    if args.min_fit_entries <= 0:
        raise ValueError("--min-fit-entries must be positive")
    args.pull_range = (float(args.pull_range[0]), float(args.pull_range[1]))


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    args.datacard = str(Path(args.datacard).resolve())
    if args.plot_only is not None:
        return run_plot_only(args, output_dir)

    datacard_path = Path(args.datacard).resolve()
    ensure_file(datacard_path, "Input datacard not found")
    run_dir = prepare_run_dir(output_dir, args.run_name)
    processes = parse_datacard_processes(datacard_path)
    schemes = selected_schemes(args, processes.signals)
    scheme_injection_map = injections_by_scheme(args, schemes)
    args.scheme_injection_map = scheme_injection_map

    log(f"Run directory: {repo_relative(run_dir)}")
    log(f"Signals: {', '.join(processes.signals)}")
    log(f"Backgrounds: {', '.join(processes.backgrounds)}")
    log(f"POI schemes: {', '.join(scheme_label(scheme.name) for scheme in schemes)}")
    log(f"Default injections: {format_injection_values(args.injections)}")
    log(f"Injections by scheme: {format_injections_by_scheme(scheme_injection_map)}")
    log(f"Requested POI range: [{format_number(args.poi_min)}, {format_number(args.poi_max)}]")
    log("Job POI ranges use the fixed requested range.")
    log(f"Dataset strategies: {', '.join(args.dataset_strategies)}")
    log(f"PDF target strategies: {', '.join(args.pdf_target_strategies)}")
    if "toys" in args.dataset_strategies:
        log(f"Toys per toy experiment: {args.toys}")
    if "asimov" in args.dataset_strategies:
        log("Asimov experiments use Combine -t -1")
    log(f"Worker policy: {worker_policy_text(args.workers, max(1, len(schemes)))}")
    log(
        f"Fit method: MultiDimFit --algo {FIT_ALGO} -P target --floatOtherPOIs 1 "
        "with other POIs floating; pdfindex is either floating or fixed by --pdf-target-strategy"
    )
    log(f"Minimizer strategy: {args.cmin_default_minimizer_strategy}")
    if args.quick:
        if args.cmin_strategy_explicit:
            log(
                "Quick mode enabled: robust fitting disabled and explicit "
                f"minimizer strategy {args.cmin_default_minimizer_strategy} kept"
            )
        else:
            log("Quick mode enabled: robust fitting disabled and minimizer strategy set to 0")
    if args.robust_fit:
        log("Robust fit enabled: MultiDimFit commands include --robustFit 1")

    all_results: list[dict[str, Any]] = []
    workspace_jobs = [build_workspace_job(run_dir, datacard_path, args.mass, scheme) for scheme in schemes]
    workspace_results = run_jobs_parallel(workspace_jobs, args.workers, "workspace builds")
    all_results.extend(workspace_results)
    if not all(result_successful(result) for result in workspace_results):
        write_run_summary(run_dir, output_dir, args, processes, schemes, {}, all_results, None)
        return 1

    workspace_paths: dict[str, Path] = {}
    staged_datacards: dict[str, Path] = {}
    for result in workspace_results:
        scheme_name = str(result["metadata"]["scheme"])
        workspace_path = result_output_path(result, WORKSPACE_OUTPUT_NAME)
        if workspace_path is None:
            raise RuntimeError(f"No {WORKSPACE_OUTPUT_NAME} output found for {scheme_name}")
        workspace_paths[scheme_name] = workspace_path
        staged_datacards[scheme_name] = Path(result["cwd"]) / LOCAL_DATACARD_NAME

    first_scheme = schemes[0]
    pdf_state_metadata = load_bundle_pdf_state_metadata(datacard_path)
    truth_pdfs, truth_pdf_metadata = discover_truth_pdfs(
        workspace_paths[first_scheme.name],
        args.pdf_families,
        pdf_state_metadata,
    )
    target_pdfs = [
        make_truth_pdf(int(item["index"]), str(item["name"]), pdf_state_metadata.get(int(item["index"]), item))
        for item in truth_pdf_metadata.get("available_target_pdfs", [])
    ] or truth_pdfs
    log(
        "Truth PDFs: "
        + ", ".join(f"{truth.index}:{truth.label}({truth.name})" for truth in truth_pdfs)
    )
    log(
        "Target PDFs: "
        + ", ".join(f"{target.index}:{target.label}({target.name})" for target in target_pdfs)
    )

    dataset_jobs: list[CommandJob] = []
    experiment_index = 0
    for dataset_strategy in args.dataset_strategies:
        for scheme in schemes:
            for target in scheme.targets:
                for truth_pdf in truth_pdfs:
                    for injected_r in scheme_injection_map[scheme.name]:
                        experiment_index += 1
                        dataset_jobs.append(
                            build_dataset_generation_job(
                                run_dir,
                                scheme,
                                target,
                                truth_pdf,
                                injected_r,
                                args.mass,
                                workspace_paths[scheme.name],
                                staged_datacards[scheme.name],
                                args,
                                dataset_strategy,
                                experiment_index,
                            )
                        )

    dataset_results = run_jobs_parallel(dataset_jobs, args.workers, "pseudo-dataset generation jobs")
    all_results.extend(dataset_results)
    if not all(result_successful(result) for result in dataset_results):
        write_run_summary(run_dir, output_dir, args, processes, schemes, truth_pdf_metadata, all_results, None)
        return 1

    if args.skip_fits:
        write_run_summary(run_dir, output_dir, args, processes, schemes, truth_pdf_metadata, all_results, None)
        log(f"Completed pseudo-dataset generation only. Summary: {repo_relative(run_dir / 'summary.json')}")
        return 0

    dataset_result_by_key: dict[tuple[str, str, int, float, str], dict[str, Any]] = {}
    for result in dataset_results:
        metadata = result.get("metadata", {})
        key = (
            str(metadata["scheme"]),
            str(metadata["target_poi"]),
            int(metadata["truth_pdf_index"]),
            float(metadata["injected_r"]),
            str(metadata.get("dataset_strategy", "toys")),
        )
        dataset_result_by_key[key] = result

    fit_jobs: list[CommandJob] = []
    for dataset_strategy in args.dataset_strategies:
        for scheme in schemes:
            for target in scheme.targets:
                for truth_pdf in truth_pdfs:
                    for injected_r in scheme_injection_map[scheme.name]:
                        dataset_result = dataset_result_by_key[
                            (scheme.name, target.poi, truth_pdf.index, float(injected_r), dataset_strategy)
                        ]
                        toys_path = first_root_file(dataset_result, "higgsCombine")
                        if toys_path is None:
                            raise RuntimeError(f"No GenerateOnly ROOT file found for {dataset_result['job_id']}")
                        if "floating" in args.pdf_target_strategies:
                            fit_jobs.append(
                                build_fit_job(
                                    run_dir,
                                    scheme,
                                    target,
                                    truth_pdf,
                                    None,
                                    injected_r,
                                    args.mass,
                                    workspace_paths[scheme.name],
                                    staged_datacards[scheme.name],
                                    toys_path,
                                    args,
                                    dataset_strategy,
                                )
                            )
                        if "fixed" in args.pdf_target_strategies:
                            for target_pdf in target_pdfs:
                                fit_jobs.append(
                                    build_fit_job(
                                        run_dir,
                                        scheme,
                                        target,
                                        truth_pdf,
                                        target_pdf,
                                        injected_r,
                                        args.mass,
                                        workspace_paths[scheme.name],
                                        staged_datacards[scheme.name],
                                        toys_path,
                                        args,
                                        dataset_strategy,
                                    )
                                )

    fit_results = run_jobs_parallel(fit_jobs, args.workers, "MultiDimFit jobs")
    analyzed_fit_results = [analyze_fit_result(result, args) for result in fit_results]
    all_results.extend(analyzed_fit_results)
    experiments = [
        experiment_summary_from_fit(result, index)
        for index, result in enumerate(analyzed_fit_results, start=1)
    ]
    scatter_plot, scatter_plot_variations = make_summary_scatter_plots_isolated(
        experiments,
        run_dir,
        args.bias_pull_threshold,
        args.annotate_scatter_outliers,
    )
    heatmap_plots = make_summary_heatmaps_isolated(experiments, run_dir)
    write_run_summary(
        run_dir,
        output_dir,
        args,
        processes,
        schemes,
        truth_pdf_metadata,
        all_results,
        scatter_plot,
        scatter_plot_variations,
        heatmap_plots,
    )

    failed = [result for result in all_results if not result_successful(result)]
    if failed:
        log(f"Completed with {len(failed)} failed command(s).")
        return 1

    log(f"Completed successfully. Summary: {repo_relative(run_dir / 'summary.json')}")
    log(f"Global summary: {repo_relative(output_dir / 'bias_study.html')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
