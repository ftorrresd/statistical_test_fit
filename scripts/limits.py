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
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.text import Text


REPO_ROOT = Path(__file__).resolve().parents[1]
COMBINE_ROOT = REPO_ROOT.parent / "HiggsAnalysis" / "CombinedLimit"
PLOT_TEST_STAT_CLS_SCRIPT = COMBINE_ROOT / "test" / "plotTestStatCLs.py"
CONSOLE = Console()
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from statistical_test_fit.display_names import poi_scheme_display, signal_target_display


class StoreCminStrategy(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):  # type: ignore[override]
        setattr(namespace, self.dest, values)
        setattr(namespace, "cmin_strategy_explicit", True)


DEFAULT_DATACARD = REPO_ROOT / "datacards" / "datacard.txt"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "limits"
DEFAULT_MASS = "125"
DEFAULT_QUANTILES = (0.16, 0.5, 0.84)
DEFAULT_POI_MIN = 0.0
DEFAULT_POI_MAX = 1_000_000.0
DEFAULT_POI_INITIAL = 1.0
DEFAULT_CMIN_DEFAULT_MINIMIZER_STRATEGY = 2
HYBRID_RANGE_ASYMPTOTIC_SCALE = 10.0
HYBRID_RETRY_RANGE_SCALE = 10.0
DEFAULT_HYBRID_RANGE_MAX_RETRIES = 3
DEFAULT_HYBRID_GRID_TOYS = 2000
DEFAULT_HYBRID_GRID_CYCLES = 1
DEFAULT_HYBRID_GRID_UPPER_SCALE = 2.0
HYBRID_GRID_RMAX_HEADROOM_SCALE = 1.5
DEFAULT_HYBRID_GRID_COARSE_POINTS = 9
DEFAULT_HYBRID_GRID_SEED = 123456
HYBRID_GRID_RELATIVE_OFFSETS = (0.70, 0.80, 0.90, 0.95, 0.975, 1.0, 1.025, 1.05, 1.10, 1.20, 1.30)
RETRYABLE_HYBRID_WARNING_PATTERNS = (
    "[WARNING] Minimization finished with status",
    "covariance matrix forced positive definite",
)
QUICK_HYBRID_TOYS = 100
QUICK_CLS_ACC = 0.02
QUICK_R_REL_ACC = 0.10
QUICK_R_ABS_ACC = 10.0
WORKSPACE_OUTPUT_NAME = "multisignal_workspace.root"
LOCAL_WORKSPACE_NAME = "workspace.root"
LOCAL_DATACARD_NAME = "datacard.txt"
LOCAL_BLIND_ASIMOV_NAME = "blind_asimov.root"
LOCAL_HYBRID_GRID_NAME = "hybrid_grid.root"
BLIND_ASIMOV_DATASET = f"{LOCAL_BLIND_ASIMOV_NAME}:toys/toy_asimov"


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
    CONSOLE.print(Text("[limits]", style="dim"), message)


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


def point_tag(value: float) -> str:
    return safe_name("r_" + format_number(value).replace("-", "m").replace("+", "p"))


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


def common_minimizer_options(cmin_strategy: int) -> list[str]:
    return [
        "--cminDefaultMinimizerStrategy",
        str(cmin_strategy),
    ]


def use_common_minimizer_options(quick: bool, cmin_strategy_explicit: bool) -> bool:
    return (not quick) or cmin_strategy_explicit


def available_cpu_count() -> int:
    if hasattr(os, "sched_getaffinity"):
        try:
            return max(1, len(os.sched_getaffinity(0)))
        except OSError:
            pass
    return max(1, os.cpu_count() or 1)


def resolved_worker_count(requested_workers: int, job_count: int) -> int:
    if job_count <= 0:
        return 0
    if requested_workers <= 0:
        return min(job_count, available_cpu_count())
    return min(requested_workers, job_count)


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
    upper_by_poi: dict[str, float] | None = None,
) -> str:
    upper_by_poi = upper_by_poi or {}
    return ":".join(
        f"{poi}={format_number(lower)},{format_number(upper_by_poi.get(poi, upper))}"
        for poi in pois
    )


def hybrid_grid_rmax_option(grid_max: float) -> float:
    return HYBRID_GRID_RMAX_HEADROOM_SCALE * grid_max


def build_workspace_job(
    run_dir: Path,
    datacard_path: Path,
    mass: str,
    scheme: PoiScheme,
) -> CommandJob:
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


def stage_common_combine_inputs(
    job_dir: Path,
    workspace_path: Path,
    datacard_path: Path,
    blind_asimov_path: Path | None = None,
) -> list[dict[str, str]]:
    records = [
        copy_file(workspace_path, job_dir / LOCAL_WORKSPACE_NAME),
        copy_file(datacard_path, job_dir / LOCAL_DATACARD_NAME),
    ]
    if blind_asimov_path is not None:
        records.append(copy_file(blind_asimov_path, job_dir / LOCAL_BLIND_ASIMOV_NAME))
    return records


def build_blind_asimov_job(
    run_dir: Path,
    scheme: PoiScheme,
    mass: str,
    workspace_path: Path,
    datacard_path: Path,
    poi_min: float,
    poi_max: float,
    cmin_strategy: int,
) -> CommandJob:
    scheme_display = poi_scheme_display(scheme.name)
    cwd = run_dir / scheme_display.slug / "blind_asimov"
    reset_directory(cwd)
    inputs = stage_common_combine_inputs(cwd, workspace_path, datacard_path)
    command = [
        "combine",
        "-M",
        "GenerateOnly",
        LOCAL_WORKSPACE_NAME,
        "-t",
        "-1",
        "--expectSignal",
        "0",
        "--saveToys",
        "--toysFrequentist",
        "--bypassFrequentistFit",
        "--redefineSignalPOIs",
        ",".join(scheme.all_pois),
        "--setParameters",
        set_parameters_argument(scheme.all_pois, 0.0),
        "--setParameterRanges",
        parameter_ranges_argument(scheme.all_pois, poi_min, poi_max),
        "-m",
        mass,
        "-n",
        f".blind_asimov.{scheme_display.slug}",
    ]
    command.extend(common_minimizer_options(cmin_strategy))
    return CommandJob(
        job_id=f"blind_asimov_{scheme_display.slug}",
        kind="blind_asimov_generation",
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
            "pois": list(scheme.all_pois),
            "poi_min": poi_min,
            "poi_max": poi_max,
            "blind_asimov_dataset": "toys/toy_asimov",
            "cmin_default_minimizer_strategy": cmin_strategy,
            "cmin_minimizer_options_applied": True,
            "blindness": [
                "GenerateOnly uses -t -1 to create a pre-fit Asimov dataset.",
                "All signal POIs are set to zero before generation.",
                "--bypassFrequentistFit prevents fitting observed data before generation.",
            ],
            "mass": mass,
        },
    )


def build_asymptotic_job(
    run_dir: Path,
    scheme: PoiScheme,
    target: PoiTarget,
    mass: str,
    workspace_path: Path,
    datacard_path: Path,
    blind_asimov_path: Path,
    poi_min: float,
    poi_max: float,
    cmin_strategy: int,
) -> CommandJob:
    scheme_display = poi_scheme_display(scheme.name)
    target_display = signal_target_display(target.label)
    cwd = run_dir / scheme_display.slug / "combine" / "asymptotic" / target_display.slug
    reset_directory(cwd)
    inputs = stage_common_combine_inputs(cwd, workspace_path, datacard_path, blind_asimov_path)
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
        "--dataset",
        BLIND_ASIMOV_DATASET,
        "-m",
        mass,
        "-n",
        f".asymptotic_blind.{target_display.slug}",
    ]
    command.extend(common_minimizer_options(cmin_strategy))
    return CommandJob(
        job_id=f"asymptotic_{scheme_display.slug}_{target_display.slug}",
        kind="combine",
        method="AsymptoticLimits",
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
            "profiled_signal_pois": [poi for poi in scheme.all_pois if poi != target.poi],
            "poi_min": poi_min,
            "poi_max": poi_max,
            "range_source": "fixed",
            "dataset_override": BLIND_ASIMOV_DATASET,
            "cmin_default_minimizer_strategy": cmin_strategy,
            "cmin_minimizer_options_applied": True,
            "blindness": [
                "AsymptoticLimits is run with --run blind.",
                "Combine therefore uses the pre-fit model state for the expected Asimov calculation instead of fitting observed data.",
                f"The command also overrides data_obs with {BLIND_ASIMOV_DATASET} as a guardrail.",
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
    blind_asimov_path: Path,
    poi_min: float,
    poi_max: float,
    hybrid_toys: int | None,
    cls_acc: float | None,
    r_rel_acc: float | None,
    r_abs_acc: float | None,
    save_hybrid_result: bool,
    default_poi_max: float,
    *,
    cmin_strategy: int,
    cmin_strategy_explicit: bool = False,
    hint_method: str | None = None,
    range_source: str = "fixed",
    range_reference: dict[str, Any] | None = None,
    attempt: int = 0,
    retry_of: str | None = None,
    retry_reasons: tuple[str, ...] = (),
    quick: bool = False,
) -> CommandJob:
    scheme_display = poi_scheme_display(scheme.name)
    target_display = signal_target_display(target.label)
    base_cwd = (
        run_dir
        / scheme_display.slug
        / "combine"
        / "hybrid_lhc"
        / target_display.slug
        / quantile_tag(quantile)
    )
    cwd = base_cwd if attempt == 0 else base_cwd / f"retry_{attempt}"
    reset_directory(cwd)
    inputs = stage_common_combine_inputs(cwd, workspace_path, datacard_path, blind_asimov_path)
    retry_suffix = f".retry{attempt}" if attempt else ""
    command = [
        "combine",
        "-M",
        "HybridNew",
        LOCAL_WORKSPACE_NAME,
        "--LHCmode",
        "LHC-limits",
        "--expectedFromGrid",
        quantile_key(quantile),
        "--dataset",
        BLIND_ASIMOV_DATASET,
        "--bypassFrequentistFit",
        "--redefineSignalPOIs",
        target.poi,
        "--setParameters",
        set_parameters_argument(scheme.all_pois, 0.0),
        "--setParameterRanges",
        parameter_ranges_argument(
            scheme.all_pois,
            poi_min,
            default_poi_max,
            {target.poi: poi_max},
        ),
        "--rMin",
        format_number(poi_min),
        "--rMax",
        format_number(poi_max),
        "-m",
        mass,
        "-n",
        f".hybrid_blind.{target_display.slug}.{quantile_tag(quantile)}{retry_suffix}",
    ]
    if not quick:
        command.extend(common_minimizer_options(cmin_strategy))
    elif cmin_strategy_explicit:
        command.extend(common_minimizer_options(cmin_strategy))
    if hint_method is not None:
        command.extend(["-H", hint_method])
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
        job_id=f"hybrid_lhc_{scheme_display.slug}_{target_display.slug}_{quantile_tag(quantile)}"
        + (f"_retry{attempt}" if attempt else ""),
        kind="combine",
        method="HybridNew",
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
            "profiled_signal_pois": [poi for poi in scheme.all_pois if poi != target.poi],
            "quantile": quantile,
            "hybrid_mode": "adaptive",
            "poi_min": poi_min,
            "poi_max": poi_max,
            "default_poi_max": default_poi_max,
            "profiled_poi_max": default_poi_max,
            "range_source": range_source,
            "range_reference": range_reference,
            "attempt": attempt,
            "retry_of": retry_of,
            "retry_reasons": list(retry_reasons),
            "hybrid_toys": hybrid_toys,
            "cls_acc": cls_acc,
            "r_rel_acc": r_rel_acc,
            "r_abs_acc": r_abs_acc,
            "hint_method": hint_method,
            "dataset_override": BLIND_ASIMOV_DATASET,
            "bypass_frequentist_fit": True,
            "r_min_option": poi_min,
            "r_max_option": poi_max,
            "cmin_default_minimizer_strategy": cmin_strategy
            if use_common_minimizer_options(quick, cmin_strategy_explicit)
            else None,
            "cmin_strategy_explicit": cmin_strategy_explicit,
            "cmin_minimizer_options_applied": use_common_minimizer_options(
                quick, cmin_strategy_explicit
            ),
            "blindness": [
                "HybridNew uses --LHCmode LHC-limits for LHC-style CLs limits.",
                "HybridNew is run directly with --expectedFromGrid for the requested quantile.",
                f"data_obs is overridden with the blind pre-fit Asimov dataset {BLIND_ASIMOV_DATASET}.",
                "--bypassFrequentistFit prevents fitting observed data before toy generation.",
            ],
            "mass": mass,
        },
    )


def input_record(path: Path) -> dict[str, str]:
    return {
        "source": str(path.resolve()),
        "source_repo_relative": repo_relative(path),
        "staged": str(path.resolve()),
        "staged_name": path.name,
    }


def build_hybrid_grid_point_job(
    run_dir: Path,
    scheme: PoiScheme,
    target: PoiTarget,
    point_value: float,
    point_index: int,
    cycle_index: int,
    seed: int,
    mass: str,
    workspace_path: Path,
    datacard_path: Path,
    blind_asimov_path: Path,
    poi_min: float,
    poi_max: float,
    default_poi_max: float,
    hybrid_toys: int,
    *,
    cmin_strategy: int,
    cmin_strategy_explicit: bool = False,
    quick: bool = False,
    range_reference: dict[str, Any] | None = None,
) -> CommandJob:
    scheme_display = poi_scheme_display(scheme.name)
    target_display = signal_target_display(target.label)
    cwd = (
        run_dir
        / scheme_display.slug
        / "combine"
        / "hybrid_grid"
        / target_display.slug
        / "points"
        / f"{point_index:04d}_{point_tag(point_value)}"
        / f"cycle_{cycle_index:03d}"
    )
    reset_directory(cwd)
    inputs = stage_common_combine_inputs(cwd, workspace_path, datacard_path, blind_asimov_path)
    r_max_option = hybrid_grid_rmax_option(poi_max)
    command = [
        "combine",
        "-M",
        "HybridNew",
        LOCAL_WORKSPACE_NAME,
        "--LHCmode",
        "LHC-limits",
        "--singlePoint",
        f"{target.poi}={format_number(point_value)}",
        "--saveHybridResult",
        "--clsAcc",
        "0",
        "-T",
        str(hybrid_toys),
        "-s",
        str(seed),
        "--dataset",
        BLIND_ASIMOV_DATASET,
        "--bypassFrequentistFit",
        "--redefineSignalPOIs",
        target.poi,
        "--setParameters",
        set_parameters_argument(scheme.all_pois, 0.0),
        "--setParameterRanges",
        parameter_ranges_argument(
            scheme.all_pois,
            poi_min,
            default_poi_max,
            {target.poi: poi_max},
        ),
        "--rMin",
        format_number(poi_min),
        "--rMax",
        format_number(r_max_option),
        "-m",
        mass,
        "-n",
        f".hybrid_grid_blind.{target_display.slug}.{point_tag(point_value)}.cycle{cycle_index:03d}",
    ]
    if not quick:
        command.extend(common_minimizer_options(cmin_strategy))
    elif cmin_strategy_explicit:
        command.extend(common_minimizer_options(cmin_strategy))

    return CommandJob(
        job_id=(
            f"hybrid_grid_point_{scheme_display.slug}_{target_display.slug}_"
            f"{point_index:04d}_cycle{cycle_index:03d}"
        ),
        kind="combine_grid_point",
        method="HybridNewGridPoint",
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
            "profiled_signal_pois": [poi for poi in scheme.all_pois if poi != target.poi],
            "hybrid_mode": "grid",
            "grid_point": point_value,
            "grid_point_index": point_index,
            "grid_cycle": cycle_index,
            "grid_seed": seed,
            "poi_min": poi_min,
            "poi_max": poi_max,
            "default_poi_max": default_poi_max,
            "profiled_poi_max": default_poi_max,
            "range_source": "asymptotic_grid",
            "range_reference": range_reference,
            "hybrid_toys": hybrid_toys,
            "grid_cls_acc": 0.0,
            "dataset_override": BLIND_ASIMOV_DATASET,
            "bypass_frequentist_fit": True,
            "r_min_option": poi_min,
            "r_max_option": r_max_option,
            "r_max_headroom_scale": HYBRID_GRID_RMAX_HEADROOM_SCALE,
            "cmin_default_minimizer_strategy": cmin_strategy
            if use_common_minimizer_options(quick, cmin_strategy_explicit)
            else None,
            "cmin_strategy_explicit": cmin_strategy_explicit,
            "cmin_minimizer_options_applied": use_common_minimizer_options(
                quick, cmin_strategy_explicit
            ),
            "blindness": [
                "HybridNew grid points use --singlePoint and --saveHybridResult.",
                "Grid generation sets --clsAcc 0 so the requested -T toys define the point precision.",
                f"data_obs is overridden with the blind pre-fit Asimov dataset {BLIND_ASIMOV_DATASET}.",
                "--bypassFrequentistFit prevents fitting observed data before toy generation.",
                "--rMax is set to 1.5 times the asymptotic-derived grid maximum so large POI points are away from the upper boundary.",
            ],
            "mass": mass,
        },
    )


def build_hybrid_grid_merge_job(
    run_dir: Path,
    scheme: PoiScheme,
    target: PoiTarget,
    root_paths: list[Path],
    range_reference: dict[str, Any],
) -> CommandJob:
    scheme_display = poi_scheme_display(scheme.name)
    target_display = signal_target_display(target.label)
    cwd = run_dir / scheme_display.slug / "combine" / "hybrid_grid" / target_display.slug / "merged"
    reset_directory(cwd)
    command = ["hadd", "-f", LOCAL_HYBRID_GRID_NAME]
    command.extend(str(path.resolve()) for path in root_paths)
    return CommandJob(
        job_id=f"hybrid_grid_merge_{scheme_display.slug}_{target_display.slug}",
        kind="hybrid_grid_merge",
        method="hadd",
        cwd=cwd,
        command=tuple(command),
        output_patterns=(LOCAL_HYBRID_GRID_NAME,),
        inputs=tuple(input_record(path) for path in root_paths),
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
            "hybrid_mode": "grid",
            "grid_input_files": len(root_paths),
            "range_reference": range_reference,
        },
    )


def build_hybrid_grid_distribution_plot_job(
    run_dir: Path,
    scheme: PoiScheme,
    target: PoiTarget,
    mass: str,
    grid_path: Path,
    range_reference: dict[str, Any],
) -> CommandJob:
    ensure_file(PLOT_TEST_STAT_CLS_SCRIPT, "Required Combine diagnostic plotting script does not exist")
    scheme_display = poi_scheme_display(scheme.name)
    target_display = signal_target_display(target.label)
    cwd = run_dir / scheme_display.slug / "combine" / "hybrid_grid" / target_display.slug / "diagnostics"
    reset_directory(cwd)
    inputs = [copy_file(grid_path, cwd / LOCAL_HYBRID_GRID_NAME)]
    output_name = f"test_stat_distributions_{target_display.slug}.root"
    command = [
        "python3",
        str(PLOT_TEST_STAT_CLS_SCRIPT),
        "--input",
        LOCAL_HYBRID_GRID_NAME,
        "--poi",
        target.poi,
        "--val",
        "all",
        "--mass",
        mass,
        "--output",
        output_name,
        "--save-as-pdf",
    ]
    return CommandJob(
        job_id=f"hybrid_grid_diagnostics_{scheme_display.slug}_{target_display.slug}",
        kind="hybrid_grid_diagnostics",
        method="plotTestStatCLs.py",
        cwd=cwd,
        command=tuple(command),
        output_patterns=(),
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
            "hybrid_mode": "grid",
            "grid_file": str(grid_path.resolve()),
            "grid_file_repo_relative": repo_relative(grid_path),
            "diagnostic_output_patterns": [output_name, f"{output_name}_*.pdf"],
            "range_reference": range_reference,
            "mass": mass,
            "blindness": [
                "plotTestStatCLs.py reads the merged HybridNew grid and plots test-statistic distributions at each grid point.",
            ],
        },
    )


def build_hybrid_grid_read_job(
    run_dir: Path,
    scheme: PoiScheme,
    target: PoiTarget,
    quantile: float,
    mass: str,
    workspace_path: Path,
    datacard_path: Path,
    blind_asimov_path: Path,
    grid_path: Path,
    poi_min: float,
    poi_max: float,
    default_poi_max: float,
    hybrid_toys: int,
    grid_cycles: int,
    grid_points: tuple[float, ...],
    *,
    cmin_strategy: int,
    cmin_strategy_explicit: bool = False,
    quick: bool = False,
    range_reference: dict[str, Any] | None = None,
) -> CommandJob:
    scheme_display = poi_scheme_display(scheme.name)
    target_display = signal_target_display(target.label)
    cwd = (
        run_dir
        / scheme_display.slug
        / "combine"
        / "hybrid_grid"
        / target_display.slug
        / "read"
        / quantile_tag(quantile)
    )
    reset_directory(cwd)
    inputs = stage_common_combine_inputs(cwd, workspace_path, datacard_path, blind_asimov_path)
    inputs.append(copy_file(grid_path, cwd / LOCAL_HYBRID_GRID_NAME))
    limit_scan_plot = f"limit_scan_{target_display.slug}_{quantile_tag(quantile)}.png"
    r_max_option = hybrid_grid_rmax_option(poi_max)
    command = [
        "combine",
        "-M",
        "HybridNew",
        LOCAL_WORKSPACE_NAME,
        "--LHCmode",
        "LHC-limits",
        "--readHybridResults",
        "--grid",
        LOCAL_HYBRID_GRID_NAME,
        "--expectedFromGrid",
        quantile_key(quantile),
        "--dataset",
        BLIND_ASIMOV_DATASET,
        "--bypassFrequentistFit",
        "--redefineSignalPOIs",
        target.poi,
        "--setParameters",
        set_parameters_argument(scheme.all_pois, 0.0),
        "--setParameterRanges",
        parameter_ranges_argument(
            scheme.all_pois,
            poi_min,
            default_poi_max,
            {target.poi: poi_max},
        ),
        "--rMin",
        format_number(poi_min),
        "--rMax",
        format_number(r_max_option),
        "--plot",
        limit_scan_plot,
        "-m",
        mass,
        "-n",
        f".hybrid_grid_blind.{target_display.slug}.{quantile_tag(quantile)}",
    ]
    if not quick:
        command.extend(common_minimizer_options(cmin_strategy))
    elif cmin_strategy_explicit:
        command.extend(common_minimizer_options(cmin_strategy))

    return CommandJob(
        job_id=f"hybrid_grid_lhc_{scheme_display.slug}_{target_display.slug}_{quantile_tag(quantile)}",
        kind="combine",
        method="HybridNew",
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
            "profiled_signal_pois": [poi for poi in scheme.all_pois if poi != target.poi],
            "quantile": quantile,
            "hybrid_mode": "grid",
            "poi_min": poi_min,
            "poi_max": poi_max,
            "default_poi_max": default_poi_max,
            "profiled_poi_max": default_poi_max,
            "range_source": "asymptotic_grid",
            "range_reference": range_reference,
            "hybrid_toys": hybrid_toys,
            "grid_cycles": grid_cycles,
            "grid_point_count": len(grid_points),
            "grid_points": [format_number(value) for value in grid_points],
            "grid_file": str(grid_path.resolve()),
            "grid_file_repo_relative": repo_relative(grid_path),
            "diagnostic_output_patterns": [limit_scan_plot],
            "limit_scan_plot": limit_scan_plot,
            "grid_cls_acc": 0.0,
            "dataset_override": BLIND_ASIMOV_DATASET,
            "bypass_frequentist_fit": True,
            "r_min_option": poi_min,
            "r_max_option": r_max_option,
            "r_max_headroom_scale": HYBRID_GRID_RMAX_HEADROOM_SCALE,
            "cmin_default_minimizer_strategy": cmin_strategy
            if use_common_minimizer_options(quick, cmin_strategy_explicit)
            else None,
            "cmin_strategy_explicit": cmin_strategy_explicit,
            "cmin_minimizer_options_applied": use_common_minimizer_options(
                quick, cmin_strategy_explicit
            ),
            "blindness": [
                "HybridNew uses --readHybridResults with a merged --singlePoint toy grid.",
                "The expected quantile is extracted with --expectedFromGrid from the precomputed grid.",
                f"data_obs is overridden with the blind pre-fit Asimov dataset {BLIND_ASIMOV_DATASET}.",
                "--bypassFrequentistFit prevents fitting observed data before any required test-statistic update.",
                "--rMax is set to 1.5 times the asymptotic-derived grid maximum so large POI limits are away from the upper boundary.",
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


def format_command_txt(
    command: list[str],
    extra_shell_parts: list[str] | None = None,
) -> str:
    parts = [shlex.quote(part) for part in command]
    if extra_shell_parts:
        parts.extend(extra_shell_parts)
    separator = " " + "\\" + "\n  "
    return separator.join(parts)


def run_command_job(job: CommandJob) -> dict[str, Any]:
    job.cwd.mkdir(parents=True, exist_ok=True)
    command_string = shlex.join(job.command)
    stdout_path = job.cwd / "stdout.txt"
    stderr_path = job.cwd / "stderr.txt"
    shell_command = (
        f"{command_string} > {shlex.quote(stdout_path.name)} "
        f"2> {shlex.quote(stderr_path.name)}"
    )
    command_path = job.cwd / "command.txt"
    command_path.write_text(
        format_command_txt(
            list(job.command),
            [
                f"> {shlex.quote(stdout_path.name)}",
                f"2> {shlex.quote(stderr_path.name)}",
            ],
        )
        + "\n",
        encoding="utf-8",
    )
    started_at = dt.datetime.now(dt.timezone.utc)
    start_time = time.monotonic()
    stdout = ""
    stderr = ""
    returncode = 0

    try:
        process = subprocess.Popen(shell_command, cwd=job.cwd, shell=True)
        returncode = int(process.wait())
    except FileNotFoundError as exc:
        returncode = 127
        stderr = str(exc)
        stderr_path.write_text(stderr, encoding="utf-8")
    except Exception as exc:  # pragma: no cover - defensive runtime guard
        returncode = 1
        stderr = f"Unhandled exception while running command: {exc}"
        stderr_path.write_text(stderr, encoding="utf-8")

    ended_at = dt.datetime.now(dt.timezone.utc)
    duration = time.monotonic() - start_time
    if stdout_path.exists():
        stdout = stdout_path.read_text(encoding="utf-8")
    if stderr_path.exists():
        stderr = stderr_path.read_text(encoding="utf-8")

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
        "command_string": command_string,
        "shell_command": shell_command,
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


def job_manifest_label(job: CommandJob) -> str:
    details: list[str] = []
    metadata = job.metadata
    if metadata.get("scheme"):
        details.append(f"scheme={metadata.get('scheme_label', metadata['scheme'])}")
    if metadata.get("target_poi"):
        details.append(f"target={metadata.get('target_display_label', metadata['target_poi'])}")
    if metadata.get("quantile") is not None:
        details.append(f"quantile={metadata['quantile']}")
    if metadata.get("grid_point") is not None:
        details.append(f"grid_point={format_number(float(metadata['grid_point']))}")
    if metadata.get("grid_cycle") is not None:
        details.append(f"cycle={metadata['grid_cycle']}")
    if metadata.get("poi_max") is not None:
        details.append(f"target_range=0,{format_number(float(metadata['poi_max']))}")
    if (
        metadata.get("profiled_poi_max") is not None
        and metadata.get("profiled_poi_max") != metadata.get("poi_max")
    ):
        details.append(
            f"profiled_range=0,{format_number(float(metadata['profiled_poi_max']))}"
        )
    if metadata.get("target_processes"):
        details.append("processes=" + ",".join(metadata["target_processes"]))
    return "; ".join(details) if details else "preparation"


def print_job_manifest(wave_name: str, jobs: list[CommandJob], workers: int) -> None:
    if not jobs:
        log(f"No jobs planned for {wave_name}.")
        return
    max_workers = resolved_worker_count(workers, len(jobs))
    log(f"Planned jobs for {wave_name}: {len(jobs)} job(s), {max_workers} worker(s)")
    log("Each listed working directory has been cleared before input staging.")
    for index, job in enumerate(jobs, start=1):
        log(f"  {index}. {job.job_id} [{job.method}; {job_manifest_label(job)}]")
        log(f"     cwd: {repo_relative(job.cwd)}")
        log(f"     command: {shlex.join(job.command)}")


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
    return {
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


def run_jobs_parallel(
    jobs: list[CommandJob],
    workers: int,
    wave_name: str,
) -> list[dict[str, Any]]:
    print_job_manifest(wave_name, jobs, workers)
    if not jobs:
        return []
    max_workers = resolved_worker_count(workers, len(jobs))
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
            result = attach_outputs(result)
            results.append(result)
            pending_job_ids.discard(job.job_id)
            completed = len(results)
            remaining = len(jobs) - completed
            log(
                f"{result['status']}: {job.job_id} "
                f"in {result.get('duration', format_duration(float(result.get('duration_seconds', 0.0))))}"
                f" [{completed}/{len(jobs)} done, {remaining} remaining]"
            )
            if result.get("status") == "ok" and remaining < 10:
                remaining_jobs = ", ".join(sorted(pending_job_ids)) or "none"
                log(f"remaining jobs ({remaining}): {remaining_jobs}")
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


def collect_diagnostic_outputs(cwd: Path, metadata: dict[str, Any]) -> list[dict[str, str]]:
    outputs: list[dict[str, str]] = []
    seen: set[Path] = set()
    patterns = metadata.get("diagnostic_output_patterns", [])
    if not isinstance(patterns, list):
        return outputs
    for pattern in patterns:
        if not isinstance(pattern, str) or not pattern:
            continue
        for path in sorted(cwd.glob(pattern)):
            if not path.is_file() or path in seen:
                continue
            seen.add(path)
            outputs.append(
                {
                    "file": str(path.resolve()),
                    "file_name": path.name,
                    "file_repo_relative": repo_relative(path),
                }
            )
    return outputs


def attach_outputs(result: dict[str, Any]) -> dict[str, Any]:
    cwd = Path(result["cwd"])
    root_outputs = [
        extract_root_file(cwd / root_name) for root_name in result.get("produced_root_files", [])
    ]
    result["root_outputs"] = root_outputs
    result["limits"] = summarize_limits(root_outputs)
    result["diagnostic_outputs"] = collect_diagnostic_outputs(cwd, result.get("metadata", {}))
    result["json_path"] = str((cwd / "result.json").resolve())
    result["json_path_repo_relative"] = repo_relative(cwd / "result.json")
    result["summary_html_path"] = str((cwd / "summary.html").resolve())
    result["summary_html_path_repo_relative"] = repo_relative(cwd / "summary.html")
    write_attached_result_outputs(result)
    return result


def write_attached_result_outputs(result: dict[str, Any]) -> None:
    cwd = Path(result["cwd"])
    (cwd / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (cwd / "summary.html").write_text(make_job_html_summary(result), encoding="utf-8")


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
    main {{ width: 100%; max-width: none; margin: 0; padding: 2.2rem; }}
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


def diagnostic_links_html(result: dict[str, Any], base_dir: Path) -> str:
    outputs = result.get("diagnostic_outputs", [])
    if not outputs:
        return "-"
    links = []
    for output in outputs:
        file_path = output.get("file")
        file_name = output.get("file_name", "diagnostic")
        if file_path:
            links.append(html_link(file_name, href_relative_to(base_dir, Path(file_path))))
        else:
            links.append(html_escape(file_name))
    return "<br>".join(links)


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
    diagnostic_rows = [
        [html_link(output.get("file_name", "-"), output.get("file_name", "-"))]
        for output in result.get("diagnostic_outputs", [])
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
    if result.get("method") in {"HybridNew", "HybridNewGridPoint"}:
        profiled_poi_max = metadata.get("profiled_poi_max", metadata.get("poi_max", DEFAULT_POI_MAX))
        status_rows.extend(
            [
                ["Hybrid mode", metadata.get("hybrid_mode") or "adaptive"],
                ["Retry", retry_summary_text(result)],
                ["Target POI range", f"0 to {format_number(float(metadata.get('poi_max', DEFAULT_POI_MAX)))}"],
                ["Profiled POI range", f"0 to {format_number(float(profiled_poi_max))}"],
                ["rMax option", metadata.get("r_max_option", "-")],
                ["Hint method", metadata.get("hint_method") or "-"],
            ]
        )
        if metadata.get("grid_point") is not None:
            status_rows.extend(
                [
                    ["Grid point", format_number(float(metadata["grid_point"]))],
                    ["Grid cycle", metadata.get("grid_cycle", "-")],
                    ["Grid seed", metadata.get("grid_seed", "-")],
                    ["Grid toys", metadata.get("hybrid_toys", "-")],
                ]
            )
        if metadata.get("grid_file") is not None:
            status_rows.append(["Grid file", metadata.get("grid_file_repo_relative", metadata.get("grid_file", "-"))])
    logs_rows = [
        ["command", html_link("command.txt", "command.txt")],
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
    <section class="card"><h2>Diagnostic Outputs</h2>{html_table_raw(['File'], diagnostic_rows)}</section>
    <section class="card"><h2>Limit Tree Entries</h2>{html_table(['File', 'Entry', 'quantileExpected', 'limit', 'limitErr', 'mh', 'iToy', 'iSeed'], make_limit_entry_rows(result.get('root_outputs', [])))}</section>
    <section class="card"><h2>Logs</h2>{html_table_raw(['Stream', 'Path'], logs_rows)}</section>
    """
    return html_document(str(result["job_id"]), body)


def result_successful(result: dict[str, Any]) -> bool:
    return int(result.get("returncode", 1)) == 0


def result_superseded(result: dict[str, Any]) -> bool:
    return bool(result.get("metadata", {}).get("superseded_by"))


def asymptotic_results_by_target(
    results: list[dict[str, Any]],
) -> dict[tuple[str, str], dict[str, Any]]:
    lookup: dict[tuple[str, str], dict[str, Any]] = {}
    for result in results:
        if result.get("method") != "AsymptoticLimits":
            continue
        metadata = result.get("metadata", {})
        scheme = metadata.get("scheme")
        target_poi = metadata.get("target_poi")
        if scheme and target_poi:
            lookup[(str(scheme), str(target_poi))] = result
    return lookup


def first_positive_expected_limit(
    result: dict[str, Any],
    quantile: float,
) -> float | None:
    values = result.get("limits", {}).get("expected", {}).get(quantile_key(quantile), [])
    if isinstance(values, (int, float)):
        values = [values]
    if not isinstance(values, list):
        return None
    for value in values:
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(numeric) and numeric > 0.0:
            return numeric
    return None


def hybrid_range_from_asymptotic_result(
    asymptotic_lookup: dict[tuple[str, str], dict[str, Any]],
    scheme: PoiScheme,
    target: PoiTarget,
    quantile: float,
) -> tuple[float, dict[str, Any]]:
    result = asymptotic_lookup.get((scheme.name, target.poi))
    if result is None:
        raise RuntimeError(
            f"No asymptotic result found for {scheme.name}/{target.poi}; "
            "cannot derive HybridNew r range."
        )

    asymptotic_limit = first_positive_expected_limit(result, quantile)
    if asymptotic_limit is None:
        raise RuntimeError(
            f"Asymptotic result {result.get('job_id')} does not contain a positive "
            f"expected limit for quantile {quantile_key(quantile)}; "
            "cannot derive HybridNew r range."
        )

    poi_max = min(DEFAULT_POI_MAX, HYBRID_RANGE_ASYMPTOTIC_SCALE * asymptotic_limit)
    if poi_max <= DEFAULT_POI_MIN:
        raise RuntimeError(
            f"Derived invalid HybridNew range for {scheme.name}/{target.poi}/"
            f"{quantile_key(quantile)}: 0,{format_number(poi_max)}"
        )

    return poi_max, {
        "strategy": "target POI max = min(DEFAULT_POI_MAX, 10 * asymptotic_expected_limit)",
        "applies_to": "redefined target POI only",
        "asymptotic_job": result.get("job_id"),
        "asymptotic_quantile": quantile_key(quantile),
        "asymptotic_limit": asymptotic_limit,
        "scale": HYBRID_RANGE_ASYMPTOTIC_SCALE,
        "cap": DEFAULT_POI_MAX,
    }


def hybrid_grid_toys(args: argparse.Namespace) -> int:
    if args.hybrid_grid_toys is not None:
        return int(args.hybrid_grid_toys)
    if args.hybrid_toys is not None:
        return int(args.hybrid_toys)
    return DEFAULT_HYBRID_GRID_TOYS


def uses_hybrid_grid(args: argparse.Namespace) -> bool:
    return args.methods == "HybridGrid"


def hybrid_mode_label(args: argparse.Namespace) -> str:
    return "grid" if uses_hybrid_grid(args) else "adaptive"


def unique_sorted_points(values: list[float], lower: float, upper: float) -> tuple[float, ...]:
    by_name: dict[str, float] = {}
    for value in values:
        if not math.isfinite(value):
            continue
        clipped = min(upper, max(lower, value))
        by_name[format_number(clipped)] = clipped
    return tuple(sorted(by_name.values()))


def hybrid_grid_from_asymptotic_result(
    asymptotic_lookup: dict[tuple[str, str], dict[str, Any]],
    scheme: PoiScheme,
    target: PoiTarget,
    quantiles: tuple[float, ...],
    poi_min: float,
    poi_cap: float,
    upper_scale: float,
    coarse_points: int,
) -> tuple[tuple[float, ...], float, dict[str, Any]]:
    result = asymptotic_lookup.get((scheme.name, target.poi))
    if result is None:
        raise RuntimeError(
            f"No asymptotic result found for {scheme.name}/{target.poi}; "
            "cannot derive HybridNew grid."
        )

    asymptotic_limits: dict[str, float] = {}
    for quantile in quantiles:
        asymptotic_limit = first_positive_expected_limit(result, quantile)
        if asymptotic_limit is None:
            raise RuntimeError(
                f"Asymptotic result {result.get('job_id')} does not contain a positive "
                f"expected limit for quantile {quantile_key(quantile)}; "
                "cannot derive HybridNew grid."
            )
        asymptotic_limits[quantile_key(quantile)] = asymptotic_limit

    largest_limit = max(asymptotic_limits.values())
    poi_max = min(poi_cap, upper_scale * largest_limit)
    if poi_max <= poi_min:
        raise RuntimeError(
            f"Derived invalid HybridNew grid range for {scheme.name}/{target.poi}: "
            f"{format_number(poi_min)},{format_number(poi_max)}"
        )

    grid_values: list[float] = [poi_min, poi_max]
    for index in range(coarse_points):
        fraction = index / max(1, coarse_points - 1)
        grid_values.append(poi_min + fraction * (poi_max - poi_min))
    for asymptotic_limit in asymptotic_limits.values():
        for relative_offset in HYBRID_GRID_RELATIVE_OFFSETS:
            grid_values.append(relative_offset * asymptotic_limit)

    grid_points = unique_sorted_points(grid_values, poi_min, poi_max)
    if len(grid_points) < 2:
        raise RuntimeError(
            f"Derived too few HybridNew grid points for {scheme.name}/{target.poi}: "
            f"{', '.join(format_number(point) for point in grid_points)}"
        )

    return grid_points, poi_max, {
        "strategy": "grid max = min(DEFAULT_POI_MAX, upper_scale * max(requested asymptotic expected limits)); grid clustered around each asymptotic expected limit",
        "applies_to": "redefined target POI only",
        "asymptotic_job": result.get("job_id"),
        "asymptotic_limits": asymptotic_limits,
        "largest_asymptotic_limit": largest_limit,
        "upper_scale": upper_scale,
        "coarse_points": coarse_points,
        "relative_offsets": list(HYBRID_GRID_RELATIVE_OFFSETS),
        "cap": poi_cap,
        "grid_min": poi_min,
        "grid_max": poi_max,
        "grid_points": [format_number(point) for point in grid_points],
    }


def hybrid_retry_reasons(result: dict[str, Any]) -> list[str]:
    if result.get("method") != "HybridNew":
        return []
    log_text = "\n".join(
        str(result.get(stream, "")) for stream in ("stderr", "stdout")
    )
    return [
        pattern
        for pattern in RETRYABLE_HYBRID_WARNING_PATTERNS
        if pattern in log_text
    ]


def update_attached_result_metadata(result: dict[str, Any], updates: dict[str, Any]) -> None:
    result.setdefault("metadata", {}).update(updates)
    if result.get("json_path") and result.get("summary_html_path"):
        write_attached_result_outputs(result)


def retry_alert_from_result(
    result: dict[str, Any],
    reasons: list[str],
    message: str,
) -> dict[str, Any]:
    metadata = result.get("metadata", {})
    return {
        "job_id": result.get("job_id"),
        "status": result.get("status"),
        "returncode": result.get("returncode"),
        "scheme": metadata.get("scheme"),
        "target_poi": metadata.get("target_poi"),
        "target_label": metadata.get("target_label"),
        "quantile": metadata.get("quantile"),
        "attempt": metadata.get("attempt", 0),
        "poi_max": metadata.get("poi_max"),
        "reasons": reasons,
        "message": message,
        "cwd": result.get("cwd_repo_relative", result.get("cwd")),
    }


def log_large_retry_alert(alert: dict[str, Any]) -> None:
    CONSOLE.rule(style="bright_red")
    CONSOLE.print("[bold bright_red]HYBRIDNEW RETRY WARNING EXHAUSTED[/bold bright_red]")
    CONSOLE.print(alert.get("message", "HybridNew retry warning remains."), style="red")
    CONSOLE.print(
        "job={job} scheme={scheme} poi={poi} quantile={quantile} attempt={attempt} range=0,{poi_max}".format(
            job=alert.get("job_id"),
            scheme=alert.get("scheme"),
            poi=alert.get("target_poi"),
            quantile=alert.get("quantile"),
            attempt=alert.get("attempt"),
            poi_max=format_number(float(alert.get("poi_max", DEFAULT_POI_MAX))),
        ),
        style="red",
    )
    CONSOLE.print("Return code 0 jobs remain successful.", style="dim red")
    CONSOLE.rule(style="bright_red")


def hybrid_retry_alerts(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    alerts: list[dict[str, Any]] = []
    for result in results:
        metadata = result.get("metadata", {})
        if metadata.get("retry_exhausted"):
            alerts.append(
                retry_alert_from_result(
                    result,
                    list(metadata.get("retry_reasons", [])),
                    str(metadata.get("retry_message", "HybridNew warning persisted after retries.")),
                )
            )
    return alerts


def retry_summary_text(result: dict[str, Any]) -> str:
    metadata = result.get("metadata", {})
    parts: list[str] = []
    if result.get("method") == "HybridNew":
        parts.append(f"attempt {metadata.get('attempt', 0)}")
    if metadata.get("retry_of"):
        parts.append(f"retry of {metadata['retry_of']}")
    if metadata.get("superseded_by"):
        parts.append(f"superseded by {metadata['superseded_by']}")
    if metadata.get("retry_exhausted"):
        parts.append("WARNING exhausted")
    return "; ".join(parts) if parts else "-"


def make_hybrid_retry_job(
    result: dict[str, Any],
    args: argparse.Namespace,
    run_dir: Path,
    scheme_by_name: dict[str, PoiScheme],
    target_by_key: dict[tuple[str, str], PoiTarget],
    workspace_paths: dict[str, Path],
    staged_datacards: dict[str, Path],
    blind_asimov_paths: dict[str, Path],
) -> CommandJob | None:
    reasons = hybrid_retry_reasons(result)
    if not reasons:
        return None

    metadata = result.get("metadata", {})
    attempt = int(metadata.get("attempt", 0))
    previous_poi_max = float(metadata.get("poi_max", DEFAULT_POI_MAX))

    if attempt >= args.hybrid_range_max_retries:
        message = (
            f"HybridNew warning persisted after {attempt} retry attempt(s); "
            "no further retries will be scheduled."
        )
        update_attached_result_metadata(
            result,
            {
                "retry_exhausted": True,
                "retry_status": "warning_after_max_retries",
                "retry_reasons": reasons,
                "retry_message": message,
            },
        )
        alert = retry_alert_from_result(result, reasons, message)
        log_large_retry_alert(alert)
        return None

    next_poi_max = HYBRID_RETRY_RANGE_SCALE * previous_poi_max
    if next_poi_max <= previous_poi_max:
        message = (
            "HybridNew warning persisted, but the computed retry r range did not "
            "increase; no further retries will be scheduled."
        )
        update_attached_result_metadata(
            result,
            {
                "retry_exhausted": True,
                "retry_status": "retry_range_not_increased",
                "retry_reasons": reasons,
                "retry_message": message,
            },
        )
        alert = retry_alert_from_result(result, reasons, message)
        log_large_retry_alert(alert)
        return None

    scheme_name = str(metadata["scheme"])
    target_poi = str(metadata["target_poi"])
    quantile = float(metadata["quantile"])
    scheme = scheme_by_name[scheme_name]
    target = target_by_key[(scheme_name, target_poi)]
    minimizer_options_applied = use_common_minimizer_options(
        args.quick,
        args.cmin_strategy_explicit,
    )
    range_reference = dict(metadata.get("range_reference") or {})
    range_reference.setdefault("retry_history", [])
    range_reference["retry_history"].append(
        {
            "source_job": result.get("job_id"),
            "attempt": attempt,
            "previous_poi_max": previous_poi_max,
            "next_poi_max": next_poi_max,
            "scale": HYBRID_RETRY_RANGE_SCALE,
            "cmin_default_minimizer_strategy": args.cmin_strategy
            if minimizer_options_applied
            else None,
            "cmin_strategy_explicit": args.cmin_strategy_explicit,
            "cmin_minimizer_options_applied": minimizer_options_applied,
            "reasons": reasons,
        }
    )

    retry_job = build_hybrid_job(
        run_dir,
        scheme,
        target,
        quantile,
        args.mass,
        workspace_paths[scheme_name],
        staged_datacards[scheme_name],
        blind_asimov_paths[scheme_name],
        args.poi_min,
        next_poi_max,
        args.hybrid_toys,
        args.cls_acc,
        args.r_rel_acc,
        args.r_abs_acc,
        not args.no_save_hybrid_result,
        args.poi_max,
        cmin_strategy=args.cmin_strategy,
        cmin_strategy_explicit=args.cmin_strategy_explicit,
        quick=args.quick,
        range_source=str(metadata.get("range_source", "asymptotic")),
        range_reference=range_reference,
        attempt=attempt + 1,
        retry_of=str(result.get("job_id")),
        retry_reasons=tuple(reasons),
    )
    update_attached_result_metadata(
        result,
        {
            "superseded_by": retry_job.job_id,
            "retry_status": "retry_scheduled",
            "retry_reasons": reasons,
            "retry_next_poi_max": next_poi_max,
            "retry_range_scale": HYBRID_RETRY_RANGE_SCALE,
            "retry_cmin_default_minimizer_strategy": args.cmin_strategy
            if minimizer_options_applied
            else None,
            "retry_cmin_minimizer_options_applied": minimizer_options_applied,
        },
    )
    cmin_log = (
        f"--cminDefaultMinimizerStrategy {args.cmin_strategy}"
        if minimizer_options_applied
        else "no --cminDefaultMinimizerStrategy in quick mode"
    )
    log(
        "Submitting immediate HybridNew retry after retryable warning: "
        f"{result.get('job_id')} -> {retry_job.job_id}; "
        f"range 0,{format_number(previous_poi_max)} -> 0,{format_number(next_poi_max)}; "
        f"{cmin_log}"
    )
    log(f"     cwd: {repo_relative(retry_job.cwd)}")
    log(f"     command: {shlex.join(retry_job.command)}")
    return retry_job


def run_hybrid_jobs_with_adaptive_retries(
    initial_jobs: list[CommandJob],
    args: argparse.Namespace,
    run_dir: Path,
    schemes: tuple[PoiScheme, ...],
    workspace_paths: dict[str, Path],
    staged_datacards: dict[str, Path],
    blind_asimov_paths: dict[str, Path],
) -> list[dict[str, Any]]:
    scheme_by_name = {scheme.name: scheme for scheme in schemes}
    target_by_key = {
        (scheme.name, target.poi): target
        for scheme in schemes
        for target in scheme.targets
    }
    all_results: list[dict[str, Any]] = []

    print_job_manifest("HybridNew limit jobs", initial_jobs, args.workers)
    if not initial_jobs:
        return []
    max_workers = resolved_worker_count(args.workers, len(initial_jobs))
    log(f"Running {len(initial_jobs)} initial HybridNew command(s) with {max_workers} worker(s)")
    log("Adaptive HybridNew retries are submitted as soon as each job finishes.")
    submitted_count = 0
    completed_count = 0
    future_to_job: dict[concurrent.futures.Future[dict[str, Any]], CommandJob] = {}
    pending_job_ids: set[str] = set()

    def submit_job(
        executor: concurrent.futures.ThreadPoolExecutor,
        job: CommandJob,
    ) -> None:
        nonlocal submitted_count
        future_to_job[executor.submit(run_command_job, job)] = job
        pending_job_ids.add(job.job_id)
        submitted_count += 1

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        for job in initial_jobs:
            submit_job(executor, job)

        while future_to_job:
            done_futures, _ = concurrent.futures.wait(
                set(future_to_job),
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            for future in done_futures:
                job = future_to_job.pop(future)
                try:
                    result = future.result()
                except Exception as exc:  # pragma: no cover - defensive runtime guard
                    result = executor_exception_result(job, exc)
                result = attach_outputs(result)
                all_results.append(result)
                pending_job_ids.discard(job.job_id)
                completed_count += 1
                log(
                    f"{result['status']}: {job.job_id} "
                    f"in {result.get('duration', format_duration(float(result.get('duration_seconds', 0.0))))}"
                    f" [{completed_count}/{submitted_count} done, {len(future_to_job)} active]"
                )

                retry_job = make_hybrid_retry_job(
                    result,
                    args,
                    run_dir,
                    scheme_by_name,
                    target_by_key,
                    workspace_paths,
                    staged_datacards,
                    blind_asimov_paths,
                )
                if retry_job is not None:
                    submit_job(executor, retry_job)

                if len(future_to_job) < 10:
                    active_jobs = ", ".join(sorted(pending_job_ids)) or "none"
                    log(f"active jobs ({len(future_to_job)}): {active_jobs}")

    return all_results


def run_hybrid_grid_workflow(
    args: argparse.Namespace,
    run_dir: Path,
    schemes: tuple[PoiScheme, ...],
    asymptotic_results: list[dict[str, Any]],
    workspace_paths: dict[str, Path],
    staged_datacards: dict[str, Path],
    blind_asimov_paths: dict[str, Path],
) -> list[dict[str, Any]]:
    asymptotic_lookup = asymptotic_results_by_target(asymptotic_results)
    grid_toys = hybrid_grid_toys(args)
    planned: dict[tuple[str, str], dict[str, Any]] = {}
    point_jobs: list[CommandJob] = []
    next_seed = int(args.hybrid_grid_seed)

    for scheme in schemes:
        for target in scheme.targets:
            grid_points, poi_max, range_reference = hybrid_grid_from_asymptotic_result(
                asymptotic_lookup,
                scheme,
                target,
                args.quantiles,
                args.poi_min,
                args.poi_max,
                args.hybrid_grid_upper_scale,
                args.hybrid_grid_coarse_points,
            )
            planned[(scheme.name, target.poi)] = {
                "scheme": scheme,
                "target": target,
                "grid_points": grid_points,
                "poi_max": poi_max,
                "range_reference": range_reference,
            }
            log(
                "HybridNew grid from asymptotic: "
                f"{scheme.name}/{target.poi} has {len(grid_points)} point(s), "
                f"r range {format_number(args.poi_min)} to {format_number(poi_max)}, "
                f"-T {grid_toys}, cycles {args.hybrid_grid_cycles}, "
                f"profiled POIs keep 0 to {format_number(args.poi_max)}"
            )
            for point_index, point_value in enumerate(grid_points):
                for cycle_index in range(args.hybrid_grid_cycles):
                    point_jobs.append(
                        build_hybrid_grid_point_job(
                            run_dir,
                            scheme,
                            target,
                            point_value,
                            point_index,
                            cycle_index,
                            next_seed,
                            args.mass,
                            workspace_paths[scheme.name],
                            staged_datacards[scheme.name],
                            blind_asimov_paths[scheme.name],
                            args.poi_min,
                            poi_max,
                            args.poi_max,
                            grid_toys,
                            cmin_strategy=args.cmin_strategy,
                            cmin_strategy_explicit=args.cmin_strategy_explicit,
                            quick=args.quick,
                            range_reference=range_reference,
                        )
                    )
                    next_seed += 1

    results: list[dict[str, Any]] = []
    point_results = run_jobs_parallel(point_jobs, args.workers, "HybridNew grid point jobs")
    results.extend(point_results)
    if not all(result_successful(result) for result in point_results):
        return results

    root_paths_by_target: dict[tuple[str, str], list[Path]] = {key: [] for key in planned}
    for result in point_results:
        metadata = result.get("metadata", {})
        key = (str(metadata.get("scheme")), str(metadata.get("target_poi")))
        root_path = result_output_path(result)
        if root_path is None:
            raise RuntimeError(f"No HybridNew grid point ROOT output found for {result.get('job_id')}")
        root_paths_by_target[key].append(root_path)

    merge_jobs = [
        build_hybrid_grid_merge_job(
            run_dir,
            plan["scheme"],
            plan["target"],
            root_paths_by_target[key],
            plan["range_reference"],
        )
        for key, plan in planned.items()
    ]
    merge_results = run_jobs_parallel(merge_jobs, args.workers, "HybridNew grid merge jobs")
    results.extend(merge_results)
    if not all(result_successful(result) for result in merge_results):
        return results

    grid_paths: dict[tuple[str, str], Path] = {}
    for result in merge_results:
        metadata = result.get("metadata", {})
        key = (str(metadata.get("scheme")), str(metadata.get("target_poi")))
        grid_path = result_output_path(result, LOCAL_HYBRID_GRID_NAME)
        if grid_path is None:
            raise RuntimeError(f"No merged HybridNew grid output found for {result.get('job_id')}")
        grid_paths[key] = grid_path

    diagnostic_jobs = [
        build_hybrid_grid_distribution_plot_job(
            run_dir,
            plan["scheme"],
            plan["target"],
            args.mass,
            grid_paths[key],
            plan["range_reference"],
        )
        for key, plan in planned.items()
    ]
    diagnostic_results = run_jobs_parallel(
        diagnostic_jobs,
        args.workers,
        "HybridNew grid diagnostic plot jobs",
    )
    results.extend(diagnostic_results)

    read_jobs: list[CommandJob] = []
    for key, plan in planned.items():
        scheme = plan["scheme"]
        target = plan["target"]
        for quantile in args.quantiles:
            read_jobs.append(
                build_hybrid_grid_read_job(
                    run_dir,
                    scheme,
                    target,
                    quantile,
                    args.mass,
                    workspace_paths[scheme.name],
                    staged_datacards[scheme.name],
                    blind_asimov_paths[scheme.name],
                    grid_paths[key],
                    args.poi_min,
                    plan["poi_max"],
                    args.poi_max,
                    grid_toys,
                    args.hybrid_grid_cycles,
                    plan["grid_points"],
                    cmin_strategy=args.cmin_strategy,
                    cmin_strategy_explicit=args.cmin_strategy_explicit,
                    quick=args.quick,
                    range_reference=plan["range_reference"],
                )
            )

    read_results = run_jobs_parallel(read_jobs, args.workers, "HybridNew grid readback jobs")
    results.extend(read_results)
    return results


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
        "command_file": result.get(
            "command_path_repo_relative",
            result.get("command_path"),
        ),
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
        "poi_min": metadata.get("poi_min"),
        "poi_max": metadata.get("poi_max"),
        "default_poi_max": metadata.get("default_poi_max"),
        "profiled_poi_max": metadata.get("profiled_poi_max"),
        "hint_method": metadata.get("hint_method"),
        "hybrid_mode": metadata.get("hybrid_mode"),
        "grid_point": metadata.get("grid_point"),
        "grid_point_index": metadata.get("grid_point_index"),
        "grid_cycle": metadata.get("grid_cycle"),
        "grid_seed": metadata.get("grid_seed"),
        "grid_point_count": metadata.get("grid_point_count"),
        "grid_cycles": metadata.get("grid_cycles"),
        "grid_file": metadata.get("grid_file_repo_relative", metadata.get("grid_file")),
        "grid_cls_acc": metadata.get("grid_cls_acc"),
        "r_max_option": metadata.get("r_max_option"),
        "dataset_override": metadata.get("dataset_override"),
        "bypass_frequentist_fit": metadata.get("bypass_frequentist_fit"),
        "range_source": metadata.get("range_source"),
        "range_reference": metadata.get("range_reference"),
        "attempt": metadata.get("attempt"),
        "retry_of": metadata.get("retry_of"),
        "superseded_by": metadata.get("superseded_by"),
        "retry_status": metadata.get("retry_status"),
        "retry_reasons": metadata.get("retry_reasons"),
        "retry_exhausted": metadata.get("retry_exhausted"),
        "retry_message": metadata.get("retry_message"),
        "diagnostic_outputs": result.get("diagnostic_outputs", []),
        "limits": result.get("limits"),
    }


def build_limit_summary(results: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for result in results:
        if result_superseded(result):
            continue
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
                "target_slug": metadata.get("target_slug"),
                "target_display_label": metadata.get("target_display_label"),
                "target_latex": metadata.get("target_latex"),
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
    retry_alerts = hybrid_retry_alerts(results)
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
        "hybrid_mode": hybrid_mode_label(args),
        "hybrid_range_from_asymptotic": args.hybrid_range_from_asymptotic,
        "hybrid_range_asymptotic_scale": HYBRID_RANGE_ASYMPTOTIC_SCALE,
        "hybrid_retry_range_scale": HYBRID_RETRY_RANGE_SCALE,
        "hybrid_grid_toys": hybrid_grid_toys(args),
        "hybrid_grid_cycles": args.hybrid_grid_cycles,
        "hybrid_grid_upper_scale": args.hybrid_grid_upper_scale,
        "hybrid_grid_rmax_headroom_scale": HYBRID_GRID_RMAX_HEADROOM_SCALE,
        "hybrid_grid_coarse_points": args.hybrid_grid_coarse_points,
        "hybrid_grid_seed": args.hybrid_grid_seed,
        "hybrid_grid_relative_offsets": list(HYBRID_GRID_RELATIVE_OFFSETS),
        "fixed_range_hybrid_hint_method": "AsymptoticLimits",
        "blind_asimov_dataset": BLIND_ASIMOV_DATASET,
        "cmin_default_minimizer_strategy": args.cmin_strategy,
        "cmin_strategy_explicit": args.cmin_strategy_explicit,
        "hybrid_range_max_retries": args.hybrid_range_max_retries,
        "hybrid_retry_alerts": retry_alerts,
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
            html_escape(result.get("metadata", {}).get("scheme_label", result.get("metadata", {}).get("scheme", "-"))),
            html_escape(result.get("metadata", {}).get("target_display_label", result.get("metadata", {}).get("target_poi", "-"))),
            html_escape(result.get("metadata", {}).get("quantile", "-")),
            html_escape(
                f"0 to {format_number(float(result.get('metadata', {}).get('poi_max')))}"
                if result.get("metadata", {}).get("poi_max") is not None
                else "-"
            ),
            html_escape(
                f"0 to {format_number(float(result.get('metadata', {}).get('profiled_poi_max')))}"
                if result.get("metadata", {}).get("profiled_poi_max") is not None
                else "-"
            ),
            html_escape(retry_summary_text(result)),
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
            diagnostic_links_html(result, run_dir),
        ]
        for result in results
    ]
    summary_rows = [
        ["Run directory", repo_relative(run_dir)],
        ["Datacard", repo_relative(Path(args.datacard).resolve())],
        ["Mass", args.mass],
        ["Methods", args.methods],
        ["Quick mode", "yes" if args.quick else "no"],
        ["Hybrid mode", hybrid_mode_label(args)],
        [
            "Hybrid range strategy",
            "asymptotic-derived target POI grid"
            if uses_hybrid_grid(args)
            else "asymptotic-derived target POI per quantile"
            if args.hybrid_range_from_asymptotic
            else "fixed default range",
        ],
        ["Hybrid grid toys per cycle", hybrid_grid_toys(args) if uses_hybrid_grid(args) else "-"],
        ["Hybrid grid cycles", args.hybrid_grid_cycles if uses_hybrid_grid(args) else "-"],
        ["Hybrid grid upper scale", args.hybrid_grid_upper_scale if uses_hybrid_grid(args) else "-"],
        ["Hybrid grid rMax headroom", HYBRID_GRID_RMAX_HEADROOM_SCALE if uses_hybrid_grid(args) else "-"],
        ["Hybrid range max retries", args.hybrid_range_max_retries],
        ["Hybrid retry range scale", HYBRID_RETRY_RANGE_SCALE],
        ["Blind Asimov dataset", BLIND_ASIMOV_DATASET],
        ["Minimizer strategy", args.cmin_strategy],
        ["Minimizer strategy explicit", "yes" if args.cmin_strategy_explicit else "no"],
        ["POI scheme request", args.poi_scheme],
        ["Selected POI schemes", ", ".join(scheme_label(scheme.name) for scheme in schemes)],
        ["Default POI range", f"{format_number(args.poi_min)} to {format_number(args.poi_max)}"],
        ["Quantiles", ", ".join(quantile_key(q) for q in args.quantiles)],
        ["Aggregate JSON", html_link("blind_limits_summary.json", "blind_limits_summary.json")],
    ]
    policy_html = f"""
      <ul>
        <li>AsymptoticLimits jobs use <code>--run blind</code>.</li>
        <li>A pre-fit background-only Asimov dataset is generated per POI scheme with <code>GenerateOnly -t -1 --bypassFrequentistFit</code>.</li>
        <li>All limit jobs override <code>data_obs</code> with <code>{BLIND_ASIMOV_DATASET}</code>.</li>
        <li>HybridNew jobs use <code>--expectedFromGrid</code> for the requested expected quantiles.</li>
        <li>HybridNew jobs pass <code>--bypassFrequentistFit</code>.</li>
        <li>Fixed-range HybridNew jobs pass <code>-H AsymptoticLimits</code>; asymptotic-derived-range HybridNew jobs do not.</li>
        <li>Non-quick Combine limit jobs pass <code>--cminDefaultMinimizerStrategy {args.cmin_strategy}</code>.</li>
        <li>Quick-mode HybridNew jobs pass <code>--cminDefaultMinimizerStrategy</code> only when <code>--cmin-strategy</code> is explicitly set.</li>
        <li>When <code>--hybrid-range-from-asymptotic</code> is used, only the <code>--redefineSignalPOIs</code> target gets <code>0,min(1000000,10*r_asymp)</code>; profiled POIs keep the default range.</li>
        <li>When <code>--methods HybridGrid</code> is used, HybridNew point jobs run <code>--singlePoint</code> with <code>--clsAcc 0</code>, merge with <code>hadd</code>, and read back expected limits with <code>--readHybridResults --grid</code>.</li>
        <li>HybridNew grid jobs set <code>--rMax</code> to <code>1.5*grid_rmax</code> so large POI limits are not clipped by Combine's default maximum or placed exactly at the upper boundary.</li>
        <li>Retryable HybridNew minimization warnings scale the r range by 10 and immediately re-enter the job queue when that job finishes, until the retry cap is reached.</li>
      </ul>
    """
    retry_alert_section = ""
    if retry_alerts:
        alert_rows = [
            [
                alert.get("job_id", "-"),
                alert.get("scheme", "-"),
                alert.get("target_poi", "-"),
                alert.get("quantile", "-"),
                alert.get("attempt", "-"),
                f"0 to {format_number(float(alert.get('poi_max', DEFAULT_POI_MAX)))}",
                "; ".join(str(reason) for reason in alert.get("reasons", [])),
                alert.get("message", "-"),
            ]
            for alert in retry_alerts
        ]
        retry_alert_section = f"""
    <section class="card"><h2>Large HybridNew Retry Error Notification</h2>
      <p><strong>WARNING:</strong> one or more HybridNew jobs still emitted retryable minimization warnings after the retry policy stopped. Return code 0 jobs are reported but remain successful.</p>
      {html_table(['Job', 'Scheme', 'Target POI', 'Quantile', 'Attempt', 'POI Range', 'Reasons', 'Message'], alert_rows)}
    </section>
        """
    body = f"""
    <h1>Blind Limit Run</h1>
    <p class="subtitle">Aggregate execution summary for Combine limit jobs.</p>
    <section class="card"><h2>Summary</h2>{html_table_raw(['Field', 'Value'], [[html_escape(k), v if str(v).startswith('<a ') else html_escape(v)] for k, v in summary_rows])}</section>
    <section class="card"><h2>Limit Policy</h2>{policy_html}</section>
    {retry_alert_section}
    <section class="card"><h2>Jobs</h2>{html_table_raw(['Job', 'Method', 'Status', 'Scheme', 'Target POI', 'Quantile', 'Target POI Range', 'Profiled POI Range', 'Retry', 'Duration', 'Summary', 'Diagnostics'], job_rows)}</section>
    """
    readme_html.write_text(html_document("Blind Limit Run", body), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run expected upper-limit workflows for the bundled simultaneous Combine card. "
            "Independent three-POI H/Z and grouped H/Z schemes are supported."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        help="Parent directory for limit runs.",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help="Run subdirectory name. Defaults to a UTC timestamp.",
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
        "--methods",
        choices=("asymptotic", "hybrid", "HybridGrid", "both"),
        default="both",
        help=(
            "Limit method(s) to run. `HybridGrid` runs asymptotic limits first, "
            "then parallel HybridNew --singlePoint toy grids, then readback with --readHybridResults."
        ),
    )
    parser.add_argument(
        "--hybrid-range-from-asymptotic",
        action="store_true",
        help=(
            "Run asymptotic limits before HybridNew and set each HybridNew r range to "
            "0,min(1000000,10*r_asymp) for the matching target and quantile. "
            "Default behavior keeps the fixed 0,1000000 range."
        ),
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
        "--cmin-strategy",
        type=int,
        choices=(0, 1, 2),
        default=DEFAULT_CMIN_DEFAULT_MINIMIZER_STRATEGY,
        action=StoreCminStrategy,
        help=(
            "Value passed to Combine as --cminDefaultMinimizerStrategy. "
            "Quick-mode HybridNew jobs keep their current behavior unless this option is explicitly set."
        ),
    )
    parser.add_argument(
        "--hybrid-grid-toys",
        type=int,
        default=None,
        help=(
            "Toys per HybridNew grid point cycle. Defaults to --hybrid-toys when set "
            f"(including --quick), otherwise {DEFAULT_HYBRID_GRID_TOYS}."
        ),
    )
    parser.add_argument(
        "--hybrid-grid-cycles",
        type=int,
        default=DEFAULT_HYBRID_GRID_CYCLES,
        help="Independent seed cycles per HybridNew grid point; outputs are merged before readback.",
    )
    parser.add_argument(
        "--hybrid-grid-upper-scale",
        type=float,
        default=DEFAULT_HYBRID_GRID_UPPER_SCALE,
        help="Set grid rMax to this factor times the largest requested asymptotic expected limit, capped by the default POI maximum.",
    )
    parser.add_argument(
        "--hybrid-grid-coarse-points",
        type=int,
        default=DEFAULT_HYBRID_GRID_COARSE_POINTS,
        help="Number of coarse linear grid points between rMin and the asymptotic-derived grid rMax.",
    )
    parser.add_argument(
        "--hybrid-grid-seed",
        type=int,
        default=DEFAULT_HYBRID_GRID_SEED,
        help="Base random seed for deterministic HybridNew grid point cycles.",
    )
    parser.add_argument(
        "--hybrid-range-max-retries",
        type=int,
        default=DEFAULT_HYBRID_RANGE_MAX_RETRIES,
        help=(
            "Maximum number of 10x-r-range HybridNew retries after retryable "
            "minimization warnings in --hybrid-range-from-asymptotic mode. "
            "Non-quick limit jobs pass the configured --cmin-strategy."
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
        help="Parallel subprocess workers per dependency wave. Use 0 to cap at the available CPU count; adaptive HybridNew retries reuse this pool immediately.",
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
    if args.hybrid_range_max_retries < 0:
        raise ValueError("--hybrid-range-max-retries must be non-negative")
    apply_quick_mode(args)
    if uses_hybrid_grid(args) and args.no_save_hybrid_result:
        raise ValueError("--no-save-hybrid-result is incompatible with --methods HybridGrid")
    if args.hybrid_grid_cycles <= 0:
        raise ValueError("--hybrid-grid-cycles must be positive")
    if args.hybrid_grid_upper_scale <= 0.0:
        raise ValueError("--hybrid-grid-upper-scale must be positive")
    if args.hybrid_grid_coarse_points < 2:
        raise ValueError("--hybrid-grid-coarse-points must be at least 2")
    if hybrid_grid_toys(args) <= 0:
        raise ValueError("Hybrid grid toys must be positive")
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
    log(f"POI schemes: {', '.join(scheme_label(scheme.name) for scheme in schemes)}")
    if args.quick:
        if uses_hybrid_grid(args):
            log("Quick mode enabled: HybridNew grid point cycles use -T 100; grid points keep --clsAcc 0")
        else:
            log(
                "Quick mode enabled: HybridNew "
                "--rRelAcc 0.10 --rAbsAcc 10 --clsAcc 0.02 -T 100"
            )
        if args.cmin_strategy_explicit:
            log(
                "Quick mode explicit --cmin-strategy enabled: HybridNew jobs pass "
                f"--cminDefaultMinimizerStrategy {args.cmin_strategy}"
            )
    log("Default POI range: 0 to 1000000")
    if args.hybrid_range_from_asymptotic:
        log("Optional HybridNew range mode enabled: target POI ranges come from asymptotic limits")
    if uses_hybrid_grid(args):
        log("HybridGrid method enabled: asymptotic limits will seed parallel --singlePoint grids")

    all_results: list[dict[str, Any]] = []

    workspace_jobs = [
        build_workspace_job(run_dir, datacard_path, args.mass, scheme) for scheme in schemes
    ]
    workspace_results = run_jobs_parallel(workspace_jobs, args.workers, "workspace builds")
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

    blind_asimov_jobs = [
        build_blind_asimov_job(
            run_dir,
            scheme,
            args.mass,
            workspace_paths[scheme.name],
            staged_datacards[scheme.name],
            args.poi_min,
            args.poi_max,
            args.cmin_strategy,
        )
        for scheme in schemes
    ]
    blind_asimov_results = run_jobs_parallel(
        blind_asimov_jobs,
        args.workers,
        "blind Asimov generation jobs",
    )
    all_results.extend(blind_asimov_results)
    if not all(result_successful(result) for result in blind_asimov_results):
        write_run_summary(run_dir, args, processes, schemes, all_results)
        return 1

    blind_asimov_paths: dict[str, Path] = {}
    for result in blind_asimov_results:
        scheme_name = result["metadata"]["scheme"]
        asimov_path = result_output_path(result)
        if asimov_path is None:
            raise RuntimeError(f"No blind Asimov output found for {scheme_name}")
        blind_asimov_paths[scheme_name] = asimov_path

    run_hybrid_grid = uses_hybrid_grid(args)
    run_hybrid = args.methods in {"hybrid", "both"} or run_hybrid_grid
    run_requested_asymptotic = args.methods in {"asymptotic", "both"}
    use_asymptotic_hybrid_ranges = (args.hybrid_range_from_asymptotic or run_hybrid_grid) and run_hybrid
    run_asymptotic = run_requested_asymptotic or use_asymptotic_hybrid_ranges
    if args.hybrid_range_from_asymptotic and not run_hybrid:
        log("Hybrid range from asymptotic requested, but HybridNew is not selected; using normal asymptotic-only flow.")

    asymptotic_jobs: list[CommandJob] = []
    if run_asymptotic:
        for scheme in schemes:
            for target in scheme.targets:
                asymptotic_jobs.append(
                    build_asymptotic_job(
                        run_dir,
                        scheme,
                        target,
                        args.mass,
                        workspace_paths[scheme.name],
                        staged_datacards[scheme.name],
                        blind_asimov_paths[scheme.name],
                        args.poi_min,
                        args.poi_max,
                        args.cmin_strategy,
                    )
                )

    if use_asymptotic_hybrid_ranges:
        asymptotic_results = run_jobs_parallel(
            asymptotic_jobs,
            args.workers,
            "Asymptotic limit jobs",
        )
        all_results.extend(asymptotic_results)
        if not all(result_successful(result) for result in asymptotic_results):
            write_run_summary(run_dir, args, processes, schemes, all_results)
            return 1

        if run_hybrid_grid:
            try:
                hybrid_results = run_hybrid_grid_workflow(
                    args,
                    run_dir,
                    schemes,
                    asymptotic_results,
                    workspace_paths,
                    staged_datacards,
                    blind_asimov_paths,
                )
            except RuntimeError as exc:
                log(str(exc))
                write_run_summary(run_dir, args, processes, schemes, all_results)
                return 1
            all_results.extend(hybrid_results)
        else:
            asymptotic_lookup = asymptotic_results_by_target(asymptotic_results)
            hybrid_jobs: list[CommandJob] = []
            try:
                for scheme in schemes:
                    for target in scheme.targets:
                        for quantile in args.quantiles:
                            poi_max, range_reference = hybrid_range_from_asymptotic_result(
                                asymptotic_lookup,
                                scheme,
                                target,
                                quantile,
                            )
                            log(
                                "HybridNew target range from asymptotic: "
                                f"{scheme.name}/{target.poi}/{quantile_key(quantile)} "
                                f"0 to {format_number(poi_max)} "
                                f"from r_asymp={format_number(float(range_reference['asymptotic_limit']))}; "
                                f"profiled POIs keep 0 to {format_number(args.poi_max)}"
                            )
                            hybrid_jobs.append(
                                build_hybrid_job(
                                    run_dir,
                                    scheme,
                                    target,
                                    quantile,
                                    args.mass,
                                    workspace_paths[scheme.name],
                                    staged_datacards[scheme.name],
                                    blind_asimov_paths[scheme.name],
                                    args.poi_min,
                                    poi_max,
                                    args.hybrid_toys,
                                    args.cls_acc,
                                    args.r_rel_acc,
                                    args.r_abs_acc,
                                    not args.no_save_hybrid_result,
                                    args.poi_max,
                                    cmin_strategy=args.cmin_strategy,
                                    cmin_strategy_explicit=args.cmin_strategy_explicit,
                                    quick=args.quick,
                                    range_source="asymptotic",
                                    range_reference=range_reference,
                                )
                            )
            except RuntimeError as exc:
                log(str(exc))
                write_run_summary(run_dir, args, processes, schemes, all_results)
                return 1
            hybrid_results = run_hybrid_jobs_with_adaptive_retries(
                hybrid_jobs,
                args,
                run_dir,
                schemes,
                workspace_paths,
                staged_datacards,
                blind_asimov_paths,
            )
            all_results.extend(hybrid_results)
    else:
        limit_jobs: list[CommandJob] = list(asymptotic_jobs if run_requested_asymptotic else [])
        if run_hybrid:
            for scheme in schemes:
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
                                blind_asimov_paths[scheme.name],
                                args.poi_min,
                                args.poi_max,
                                args.hybrid_toys,
                                args.cls_acc,
                                args.r_rel_acc,
                                args.r_abs_acc,
                                not args.no_save_hybrid_result,
                                args.poi_max,
                                cmin_strategy=args.cmin_strategy,
                                cmin_strategy_explicit=args.cmin_strategy_explicit,
                                quick=args.quick,
                                hint_method="AsymptoticLimits",
                            )
                        )

        limit_results = run_jobs_parallel(limit_jobs, args.workers, "Combine limit jobs")
        all_results.extend(limit_results)
    write_run_summary(run_dir, args, processes, schemes, all_results)

    failed = [
        result
        for result in all_results
        if not result_successful(result) and not result_superseded(result)
    ]
    if failed:
        log(f"Completed with {len(failed)} failed command(s).")
        return 1

    retry_alert_count = len(hybrid_retry_alerts(all_results))
    if retry_alert_count:
        log(
            f"Completed with {retry_alert_count} HybridNew retry warning alert(s); "
            "returning 0 because final command return codes were successful."
        )

    log(f"Completed successfully. Summary: {repo_relative(run_dir / 'blind_limits_summary.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
