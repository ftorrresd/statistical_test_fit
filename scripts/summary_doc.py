#!/usr/bin/env python3

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IMAGE = "docker://ghcr.io/xu-cheng/texlive-full:latest"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "summary-docs"
DEFAULT_OUTPUT_NAME = "summary_doc"

PLOT_EXTENSIONS = {".pdf"}
TEXT_ARTIFACT_NAMES: set[str] = set()
PULL_DISTRIBUTION_NAME = "pull_distribution.pdf"
LIMIT_METHOD_LABELS = {
    "asymptotic": "asymptotic CLs",
    "hybrid_lhc": "HybridNew CLs",
    "hybrid_grid": "HybridNew grid CLs",
}
LIMIT_SCHEME_LABELS = {
    "h_grouped": "grouped H",
    "z_grouped": "grouped Z",
    "three_poi_h": "individual H",
    "three_poi_z": "individual Z",
}

@dataclass(frozen=True)
class Artifact:
    source: Path
    kind: str
    caption: str
    reason: str


@dataclass(frozen=True)
class RenderedArtifact:
    source: Path
    kind: str
    caption: str
    reason: str
    asset_relative: str


@dataclass(frozen=True)
class Section:
    key: str
    title: str
    source_dir: Path
    main_artifacts: list[Artifact]
    appendix_artifacts: list[Artifact]
    grid_groups: list["GridGroup"]
    notes: list[str]
    warnings: list[str]
    description_tex: str = ""


@dataclass(frozen=True)
class RenderedSection:
    key: str
    title: str
    source_dir: Path
    main_artifacts: list[RenderedArtifact]
    appendix_artifacts: list[RenderedArtifact]
    grid_groups: list["RenderedGridGroup"]
    notes: list[str]
    warnings: list[str]
    description_tex: str = ""


@dataclass(frozen=True)
class GridItem:
    artifact: Artifact
    label: str
    sort_key: tuple[Any, ...]


@dataclass(frozen=True)
class GridGroup:
    title: str
    caption: str
    sort_key: tuple[Any, ...]
    items: list[GridItem]
    placement: str = "appendix"
    columns: int = 3
    max_items_per_figure: int = 12
    image_height: str = r"0.145\textheight"
    figure_spec: str = "p"


@dataclass(frozen=True)
class RenderedGridItem:
    artifact: RenderedArtifact
    label: str
    sort_key: tuple[Any, ...]


@dataclass(frozen=True)
class RenderedGridGroup:
    title: str
    caption: str
    sort_key: tuple[Any, ...]
    items: list[RenderedGridItem]
    placement: str
    columns: int
    max_items_per_figure: int
    image_height: str
    figure_spec: str


def log(message: str) -> None:
    print(f"[summary_doc] {message}")


def repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def is_relative_to(path: Path, base: Path) -> bool:
    try:
        path.resolve().relative_to(base.resolve())
        return True
    except ValueError:
        return False


def resolve_input_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path.resolve()


def latex_escape_text(value: Any) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in text)


def latex_escape_caption(value: Any) -> str:
    text = str(value)
    patterns = [
        (r"m\s*_?\s*\{\s*#mu\s*#mu\s*#gamma\s*\}", lambda m: r"\(m_{\mu\mu\gamma}\)"),
        (r"m\s*_?\s*\{\s*#mu\s*#mu\s*\}", lambda m: r"\(m_{\mu\mu}\)"),
        (r"H\s*/\s*Z\s*->\s*Upsilon\s*\(\s*nS\s*\)\s*gamma", lambda m: r"$H/Z \to \Upsilon(nS)\gamma$"),
        (
            r"\b(H|Z)\s*->\s*Upsilon\s*\(\s*(\d+S|nS)\s*\)\s*gamma\b",
            lambda m: rf"${m.group(1)} \to \Upsilon({m.group(2)})\gamma$",
        ),
    ]
    placeholders: dict[str, str] = {}
    for index, (pattern, repl_fn) in enumerate(patterns):
        placeholder = f"@@ROOTMATH{index}@@"
        placeholders[placeholder] = None

        def _repl(m: re.Match[str], fn=repl_fn, token: str = placeholder) -> str:
            placeholders[token] = fn(m)
            return token

        text = re.sub(pattern, _repl, text)

    escaped = latex_escape_text(text)
    for placeholder, replacement in placeholders.items():
        if replacement is not None:
            escaped = escaped.replace(placeholder, replacement)
    return escaped


def latex_path(path: str) -> str:
    return rf"\detokenize{{{path}}}"


def slugify(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("._") or "artifact"


def title_from_stem(stem: str) -> str:
    text = stem.replace("__", " / ").replace("_", " ").replace("-", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text[:1].upper() + text[1:] if text else "Artifact"


def format_number(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not number == number:
        return "nan"
    if number.is_integer():
        return str(int(number))
    return f"{number:.4g}"


def human_join(items: Iterable[str]) -> str:
    values = [item for item in items if item]
    if not values:
        return ""
    if len(values) == 1:
        return values[0]
    if len(values) == 2:
        return f"{values[0]} and {values[1]}"
    return ", ".join(values[:-1]) + f", and {values[-1]}"


def mass_projection_from_name(name: str) -> str:
    if "m_mumugamma" in name or "_boson_" in name:
        return "boson-mass"
    if "m_mumu" in name or "_upsilon_" in name:
        return "dimuon-mass"
    return "mass"


def signal_caption_from_path(path: Path) -> str:
    projection = mass_projection_from_name(path.name)
    process = "H" if "HToUps" in path.name or "ggH" in path.name else "Z"
    state_match = re.search(r"Ups(?:ilon)?([123])S", path.name)
    state = f"Upsilon({state_match.group(1)}S)" if state_match else "Upsilon(nS)"
    return (
        f"Fitted signal-shape {projection} projection for {process} -> {state} gamma; "
        "the simulated Run 2 signal sample is compared with the parametric PDF used in the workspace"
    )


def resonant_caption_from_path(path: Path) -> str:
    name = path.name
    if name.startswith("normalization_extrapolation"):
        return "Control-region normalization extrapolation used to scale resonant background yields into the signal region"
    cr_match = re.search(r"MC_data_(CR[1-4])_fit", name)
    if cr_match:
        return (
            f"Resonant-background control-region fit in {cr_match.group(1)}; "
            "the data/MC comparison constrains the normalization transfer used for the peaking components"
        )
    if "HiggsDalitz" in name:
        projection = mass_projection_from_name(name)
        return f"Higgs Dalitz resonant-background {projection} projection used for the peaking H-background model"
    if "ZG" in name:
        return "Z gamma resonant-background boson-mass fit used for the peaking Z-background model"
    if "_Z_MC_" in name:
        projection = mass_projection_from_name(name)
        return f"Z resonant-background {projection} projection used for the peaking Z-background model"
    return title_from_stem(path.stem)


def non_resonant_projection_label(projection_key: str) -> str:
    if projection_key == "mumugamma":
        return "boson-mass"
    if projection_key == "mumu":
        return "dimuon-mass"
    return projection_key


def generic_caption_from_path(path: Path) -> str:
    name = path.name
    if name.startswith("pull_mean_sigma_scatter"):
        return "Bias-study pull mean and width summary; each point summarizes one valid toy pull fit"
    if name == PULL_DISTRIBUTION_NAME:
        return "Toy pull distribution for one injected-signal and truth-PDF configuration"
    if "branching_fraction_limits" in name:
        return (
            "H/Z -> Upsilon(nS) gamma branching-fraction upper-limit tables; "
            "columns report raw signal-strength limits and theory-scaled branching-fraction limits"
        )
    return title_from_stem(path.stem)


def read_json(path: Path, warnings: list[str]) -> dict[str, Any] | None:
    if not path.exists() or not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        warnings.append(f"Could not read JSON {repo_relative(path)}: {exc}")
        return None


def limit_table_metadata_candidates(path: Path) -> list[Path]:
    name = path.name
    candidates: list[Path] = []
    if name.endswith(".table.tex"):
        candidates.append(path.with_name(name.removesuffix(".table.tex") + ".json"))
    if name.endswith(".table.pdf"):
        candidates.append(path.with_name(name.removesuffix(".table.pdf") + ".json"))
    if path.suffix.lower() == ".pdf":
        candidates.append(path.with_suffix(".json"))
    if path.suffix.lower() == ".tex":
        candidates.append(path.with_suffix(".json"))
    return list(dict.fromkeys(candidates))


def limit_method_label(method: Any) -> str:
    method_key = str(method)
    return LIMIT_METHOD_LABELS.get(method_key, title_from_stem(method_key))


def limit_scheme_label(scheme: Any) -> str:
    scheme_key = str(scheme)
    return LIMIT_SCHEME_LABELS.get(scheme_key, title_from_stem(scheme_key))


def ordered_limit_values(values: list[Any], known_order: dict[str, str]) -> list[Any]:
    order = {key: index for index, key in enumerate(known_order)}
    return sorted(values, key=lambda value: (order.get(str(value), 999), str(value)))


def limit_table_has_observed_values(metadata: dict[str, Any]) -> bool:
    tables = metadata.get("tables", [])
    if not isinstance(tables, list):
        return False
    for table in tables:
        if not isinstance(table, dict):
            continue
        rows = table.get("rows", [])
        if not isinstance(rows, list):
            continue
        for row in rows:
            if not isinstance(row, dict):
                continue
            if (
                row.get("observed_strength_limit") is not None
                or row.get("observed_bf_limit") is not None
            ):
                return True
    return False


def limit_caption_from_path(path: Path, warnings: list[str]) -> str:
    metadata = None
    for metadata_path in limit_table_metadata_candidates(path):
        metadata = read_json(metadata_path, warnings)
        if metadata is not None:
            break

    if not metadata:
        return generic_caption_from_path(path)

    methods = metadata.get("methods", [])
    if not isinstance(methods, list):
        methods = []
    schemes = metadata.get("schemes", [])
    if not isinstance(schemes, list):
        schemes = []

    method_text = human_join(
        limit_method_label(method)
        for method in ordered_limit_values(methods, LIMIT_METHOD_LABELS)
    )
    scheme_text = human_join(
        limit_scheme_label(scheme)
        for scheme in ordered_limit_values(schemes, LIMIT_SCHEME_LABELS)
    )
    observed_text = "observed and expected" if limit_table_has_observed_values(metadata) else "blind expected"

    details = [
        f"{observed_text.capitalize()} 95% CL upper limits for H/Z -> Upsilon(nS) gamma",
    ]
    if method_text:
        details.append(f"computed with {method_text}")
    if scheme_text:
        details.append(f"for {scheme_text} POI schemes")
    caption = "; ".join(details)
    return (
        f"{caption}. Tables list raw Combine signal-strength limits r and the corresponding "
        "branching-fraction limits after scaling by the theory branching fractions; expected entries quote "
        "the median with the +/-1 sigma band"
    )


def is_latex_fragment(path: Path) -> bool:
    try:
        sample = path.read_text(encoding="utf-8", errors="replace")[:4096]
    except OSError:
        return False
    return "\\documentclass" not in sample and "\\begin{document}" not in sample


def classify_artifact(path: Path) -> str | None:
    suffix = path.suffix.lower()
    if suffix in PLOT_EXTENSIONS:
        return "plot"
    if suffix == ".tex" and is_latex_fragment(path):
        return "table"
    if path.name in TEXT_ARTIFACT_NAMES:
        return "text"
    return None


def artifact_from_path(
    path: Path,
    *,
    caption: str | None = None,
    reason: str = "matched artifact",
) -> Artifact | None:
    if not path.exists() or not path.is_file():
        return None
    kind = classify_artifact(path)
    if kind is None:
        return None
    return Artifact(
        source=path.resolve(),
        kind=kind,
        caption=caption or generic_caption_from_path(path),
        reason=reason,
    )


def artifact_group_key(artifact: Artifact) -> str:
    if artifact.kind == "plot":
        return str(artifact.source.with_suffix("").resolve())
    return str(artifact.source.resolve())


def plot_extension_priority(path: Path) -> int:
    return {".pdf": 0}.get(path.suffix.lower(), 99)


def dedupe_artifacts(artifacts: Iterable[Artifact]) -> list[Artifact]:
    selected: dict[str, Artifact] = {}
    for artifact in artifacts:
        key = artifact_group_key(artifact)
        existing = selected.get(key)
        if existing is None:
            selected[key] = artifact
            continue
        if artifact.kind == "plot" and plot_extension_priority(artifact.source) < plot_extension_priority(existing.source):
            selected[key] = artifact
    return sorted(selected.values(), key=lambda item: (item.kind, str(item.source)))


def collect_supported_artifacts(root: Path, output_dir: Path) -> list[Artifact]:
    """Deprecated — appendix is no longer auto-populated from directory scans.
    Kept for reference; no longer called in the document build pipeline."""
    del output_dir


def resolve_payload_path(value: Any, section_dir: Path) -> Path | None:
    if not isinstance(value, str) or not value.strip():
        return None
    path = Path(value).expanduser()
    candidates = []
    if path.is_absolute():
        candidates.append(path)
    else:
        candidates.append(section_dir / path)
        candidates.append(REPO_ROOT / path)
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return candidates[-1].resolve() if candidates else None


def sorted_globs(root: Path, patterns: Iterable[str]) -> list[Path]:
    paths: list[Path] = []
    if not root.exists() or not root.is_dir():
        return paths
    for pattern in patterns:
        paths.extend(sorted(root.glob(pattern)))
    return paths


def select_signal_artifacts(root: Path, warnings: list[str]) -> list[Artifact]:
    del warnings
    patterns = [
        "signal_fit_boson_*.pdf",
        "signal_fit_upsilon_*.pdf",
    ]
    return dedupe_artifacts(
        artifact
        for path in sorted_globs(root, patterns)
        for artifact in [
            artifact_from_path(path, caption=signal_caption_from_path(path), reason="signal projection plot")
        ]
        if artifact is not None
    )


def signal_pair_process(path: Path) -> tuple[str, int, int]:
    name = path.name
    if "HToUps" in name or "ggH" in name:
        process = "H"
        boson_order = 1
    else:
        process = "Z"
        boson_order = 0
    state_match = re.search(r"Ups(?:ilon)?([123])S", name)
    state_order = int(state_match.group(1)) if state_match else 4
    return process, boson_order, state_order


def signal_pair_stem(path: Path) -> str:
    stem = path.stem
    for marker in ("_boson_", "_upsilon_"):
        if marker in stem:
            return stem.replace(marker, "_PAIR_")
    return stem


def build_signal_grid_groups(root: Path, warnings: list[str]) -> list[GridGroup]:
    del warnings
    boson_paths = sorted(filter(bool, [
        artifact_from_path(p, caption=signal_caption_from_path(p), reason="signal grid plot")
        for p in sorted_globs(root, ["signal_fit_boson_*.pdf"])
    ]), key=lambda a: signal_pair_process(a.source))
    upsilon_paths = sorted(filter(bool, [
        artifact_from_path(p, caption=signal_caption_from_path(p), reason="signal grid plot")
        for p in sorted_globs(root, ["signal_fit_upsilon_*.pdf"])
    ]), key=lambda a: signal_pair_process(a.source))

    upsilon_by_stem: dict[str, Artifact] = {
        signal_pair_stem(a.source): a for a in upsilon_paths
    }

    items: list[GridItem] = []
    index = 0
    for boson in boson_paths:
        stem = signal_pair_stem(boson.source)
        upsilon = upsilon_by_stem.get(stem)
        if upsilon is None:
            continue
        process, _boson_order, state_order = signal_pair_process(boson.source)
        state_label = f"Upsilon({state_order}S)" if state_order != 4 else "Upsilon(nS)"
        pair_label = f"{process} -> {state_label} gamma"
        items.append(GridItem(artifact=boson, label=f"{pair_label} boson mass", sort_key=(index,)))
        index += 1
        items.append(GridItem(artifact=upsilon, label=f"{pair_label} dimuon mass", sort_key=(index,)))
        index += 1
    if not items:
        return []
    return [
        GridGroup(
            title="Signal Model Projections",
            caption="Fitted signal-shape boson-mass and dimuon-mass projections for all H/Z -> Upsilon(nS) gamma processes",
            sort_key=(0,),
            items=items,
            placement="main",
            columns=2,
            max_items_per_figure=2,
            image_height=r"0.38\textheight",
            figure_spec="H",
        )
    ]


def select_resonant_artifacts(root: Path, warnings: list[str]) -> list[Artifact]:
    del warnings
    patterns = [
        "resonant_background_fit_HiggsDalitz_MC_m_mumugamma.*",
        "resonant_background_fit_HiggsDalitz_MC_m_mumu.*",
        "resonant_background_fit_Z_MC_m_mumugamma.*",
        "resonant_background_fit_Z_MC_m_mumu.*",
        "resonant_background_fit_ZG_MC_m_mumugamma.*",
        "normalization_extrapolation.*",
    ]
    return dedupe_artifacts(
        artifact
        for path in sorted_globs(root, patterns)
        for artifact in [
            artifact_from_path(path, caption=resonant_caption_from_path(path), reason="resonant-background summary plot")
        ]
        if artifact is not None
    )


def control_region_sort_key(path: Path) -> tuple[Any, ...]:
    match = re.search(r"CR([1-4])", path.name)
    return (int(match.group(1)) if match else 99, path.name)


def build_resonant_grid_groups(root: Path, warnings: list[str]) -> list[GridGroup]:
    del warnings
    cr_items: list[GridItem] = []

    for path in sorted_globs(root, ["resonant_background_MC_data_CR*_fit_m_mumugamma.pdf"]):
        artifact = artifact_from_path(
            path,
            caption=resonant_caption_from_path(path),
            reason="resonant-background control-region grid plot",
        )
        if artifact is None:
            continue
        match = re.search(r"CR([1-4])", path.name)
        label = f"CR{match.group(1)}" if match else title_from_stem(path.stem)
        cr_items.append(GridItem(artifact=artifact, label=label, sort_key=control_region_sort_key(path)))

    z_items: list[GridItem] = []
    z_index = 0
    for tag, pair_label in [("Z", "Z Resonant Background")]:
        boson_artifact = None
        dimuon_artifact = None
        for candidate in sorted_globs(root, [f"resonant_background_fit_{tag}_MC_m_mumugamma.pdf"]):
            boson_artifact = artifact_from_path(
                candidate,
                caption=resonant_caption_from_path(candidate),
                reason="resonant-background projection grid plot",
            )
            if boson_artifact is not None:
                break
        for candidate in sorted_globs(root, [f"resonant_background_fit_{tag}_MC_m_mumu.pdf"]):
            dimuon_artifact = artifact_from_path(
                candidate,
                caption=resonant_caption_from_path(candidate),
                reason="resonant-background projection grid plot",
            )
            if dimuon_artifact is not None:
                break
        if boson_artifact is None or dimuon_artifact is None:
            continue
        z_items.append(GridItem(artifact=boson_artifact, label=f"{pair_label} boson mass", sort_key=(z_index,)))
        z_index += 1
        z_items.append(GridItem(artifact=dimuon_artifact, label=f"{pair_label} dimuon mass", sort_key=(z_index,)))
        z_index += 1

    hd_items: list[GridItem] = []
    hd_index = 0
    for tag, pair_label in [("HiggsDalitz", "Higgs Dalitz Resonant Background")]:
        boson_artifact = None
        dimuon_artifact = None
        for candidate in sorted_globs(root, [f"resonant_background_fit_{tag}_MC_m_mumugamma.pdf"]):
            boson_artifact = artifact_from_path(
                candidate,
                caption=resonant_caption_from_path(candidate),
                reason="resonant-background projection grid plot",
            )
            if boson_artifact is not None:
                break
        for candidate in sorted_globs(root, [f"resonant_background_fit_{tag}_MC_m_mumu.pdf"]):
            dimuon_artifact = artifact_from_path(
                candidate,
                caption=resonant_caption_from_path(candidate),
                reason="resonant-background projection grid plot",
            )
            if dimuon_artifact is not None:
                break
        if boson_artifact is None or dimuon_artifact is None:
            continue
        hd_items.append(GridItem(artifact=boson_artifact, label=f"{pair_label} boson mass", sort_key=(hd_index,)))
        hd_index += 1
        hd_items.append(GridItem(artifact=dimuon_artifact, label=f"{pair_label} dimuon mass", sort_key=(hd_index,)))
        hd_index += 1

    groups: list[GridGroup] = []
    if cr_items:
        groups.append(
            GridGroup(
                title="Resonant Background Control Regions",
                caption="Control-region fits CR1-CR4 for the resonant background normalization transfer",
                sort_key=(0,),
                items=cr_items,
                placement="main",
                columns=2,
                max_items_per_figure=4,
                image_height=r"0.30\textheight",
                figure_spec="H",
            )
        )
    if z_items:
        groups.append(
            GridGroup(
                title="Z Resonant Background",
                caption="Z resonant-background boson-mass and dimuon-mass projection plots used for the peaking Z-background model",
                sort_key=(1,),
                items=z_items,
                placement="main",
                columns=2,
                max_items_per_figure=2,
                image_height=r"0.38\textheight",
                figure_spec="H",
            )
        )
    if hd_items:
        groups.append(
            GridGroup(
                title="Higgs Dalitz Resonant Background",
                caption="Higgs Dalitz resonant-background boson-mass and dimuon-mass projection plots used for the peaking H-background model",
                sort_key=(1,),
                items=hd_items,
                placement="main",
                columns=2,
                max_items_per_figure=2,
                image_height=r"0.38\textheight",
                figure_spec="H",
            )
        )
    return groups


def select_non_resonant_artifacts(root: Path, warnings: list[str]) -> list[Artifact]:
    summary = read_json(root / "non_resonant_fit_summary.json", warnings)
    artifacts: list[Artifact] = []
    if summary:
        families = summary.get("families", {})
        if isinstance(families, dict):
            for family_name, family_payload in sorted(families.items()):
                if not isinstance(family_payload, dict):
                    continue
                family_label = str(family_payload.get("family_label") or family_name)
                plots = family_payload.get("plots", {})
                if not isinstance(plots, dict):
                    continue
                for projection_key in ("mumugamma", "mumu"):
                    path = resolve_payload_path(plots.get(projection_key), root)
                    if path is None:
                        continue
                    artifact = artifact_from_path(
                        path,
                        caption=(
                            f"{family_label} non-resonant sideband {non_resonant_projection_label(projection_key)} "
                            "projection used to validate the candidate 2D background shape"
                        ),
                        reason="non-resonant family summary plot",
                    )
                    if artifact is not None:
                        artifacts.append(artifact)
    if artifacts:
        dimuon_fit = artifact_from_path(
            root / "dimuon_mass_fit_Run2.pdf",
            caption=(
                "Run 2 dimuon-mass fit that defines the Upsilon component of the non-resonant "
                "background model before the 2D sideband fit"
            ),
            reason="non-resonant dimuon fit summary plot",
        )
        if dimuon_fit is not None:
            artifacts.append(dimuon_fit)
        return dedupe_artifacts(artifacts)

    patterns = ["mumugamma_*.pdf", "mumu_*.pdf", "dimuon_mass_fit_Run2.pdf"]
    return dedupe_artifacts(
        artifact
        for path in sorted_globs(root, patterns)
        for artifact in [artifact_from_path(path, reason="non-resonant projection plot")]
        if artifact is not None
    )


def non_resonant_family_label(family_name: str) -> str:
    mapping = {
        "Bernstein": "Bernstein",
        "Chebychev": "Chebyshev",
        "Johnson": "Johnson",
    }
    return mapping.get(family_name, family_name)


def build_non_resonant_grid_groups(root: Path, warnings: list[str]) -> list[GridGroup]:
    del warnings
    items: list[GridItem] = []
    index = 0
    for family in ("Bernstein", "Chebychev", "Johnson"):
        dimuon = artifact_from_path(
            root / f"mumu_bkg_only_{family}.pdf",
            caption=f"{family} non-resonant sideband dimuon-mass projection",
            reason="non-resonant family grid plot",
        )
        boson = artifact_from_path(
            root / f"mumugamma_bkg_only_{family}.pdf",
            caption=f"{family} non-resonant sideband boson-mass projection",
            reason="non-resonant family grid plot",
        )
        if dimuon is None or boson is None:
            continue
        family_display = non_resonant_family_label(family)
        items.append(GridItem(artifact=dimuon, label=f"{family_display} dimuon mass", sort_key=(index,)))
        index += 1
        items.append(GridItem(artifact=boson, label=f"{family_display} boson mass", sort_key=(index,)))
        index += 1
    if not items:
        return []
    return [
        GridGroup(
            title="Non-Resonant Background Projections",
            caption="Non-resonant 2D sideband projection plots used to validate each candidate background model family",
            sort_key=(0,),
            items=items,
            placement="main",
            columns=2,
            max_items_per_figure=2,
            image_height=r"0.38\textheight",
            figure_spec="H",
        )
    ]


def bundled_projection_caption(key: str, plot_info: dict[str, Any] | None = None) -> str:
    if key == "mumu":
        return (
            "Projection on m_{#mu#mu} (linear scale); data points from the blinded Z and H "
            "boson-mass regions are excluded"
        )
    if key == "mumugamma_log":
        return "Projection on m_{#mu#mu#gamma} (log scale)"
    if plot_info is not None:
        return str(plot_info.get("caption") or title_from_stem(key))
    return title_from_stem(key)


def select_statistical_model_artifacts(root: Path, warnings: list[str]) -> list[Artifact]:
    artifacts: list[Artifact] = []
    selected_projection_keys = ("mumu", "mumugamma_log")
    summary = read_json(root / "bundle_summary.json", warnings)
    if summary:
        projection_plots = summary.get("projection_plots", {}).get("plots", {})
        if isinstance(projection_plots, dict):
            for key in selected_projection_keys:
                plot_info = projection_plots.get(key)
                if not isinstance(plot_info, dict):
                    continue
                path = None
                for path_key in ("pdf_path", "path"):
                    path = resolve_payload_path(plot_info.get(path_key), root)
                    if path is not None and path.exists():
                        break
                if path is None:
                    continue
                artifact = artifact_from_path(
                    path,
                    caption=bundled_projection_caption(key, plot_info),
                    reason="bundled model projection plot",
                )
                if artifact is not None:
                    artifacts.append(artifact)

    if not artifacts:
        fallback_paths = [
            ("mumu", root / "bundle_model_projection_mumu.pdf"),
            ("mumugamma_log", root / "bundle_model_projection_mumugamma_log.pdf"),
        ]
        for key, path in fallback_paths:
            artifact = artifact_from_path(
                path,
                caption=bundled_projection_caption(key),
                reason="bundled model projection plot",
            )
            if artifact is not None:
                artifacts.append(artifact)

    for path in sorted(root.glob("*.table.tex")) if root.exists() else []:
        artifact = artifact_from_path(path, reason="statistical-model table fragment")
        if artifact is not None:
            artifacts.append(artifact)

    return dedupe_artifacts(artifacts)


def build_statistical_model_grid_groups(root: Path, warnings: list[str]) -> list[GridGroup]:
    summary = read_json(root / "bundle_summary.json", warnings)
    items: list[GridItem] = []
    pairs = [
        ("mumu", "bundle_model_projection_mumu.pdf", "dimuon mass (linear)",
         bundled_projection_caption("mumu")),
        ("mumugamma_log", "bundle_model_projection_mumugamma_log.pdf", "boson mass (log)",
         bundled_projection_caption("mumugamma_log")),
    ]
    for key, fallback_name, label, caption in pairs:
        path = None
        if summary:
            plot_info = summary.get("projection_plots", {}).get("plots", {}).get(key)
            if isinstance(plot_info, dict):
                for path_key in ("pdf_path", "path"):
                    path = resolve_payload_path(plot_info.get(path_key), root)
                    if path is not None and path.exists():
                        break
        if path is None:
            path = root / fallback_name
        artifact = artifact_from_path(
            path, caption=caption,
            reason="bundled model projection grid plot",
        )
        if artifact is not None:
            items.append(GridItem(artifact=artifact, label=label, sort_key=(len(items),)))
    if not items:
        return []
    return [
        GridGroup(
            title="Bundled Workspace Projections",
            caption="Pre-fit and post-fit projections of the bundled Combine workspace on the dimuon and boson masses",
            sort_key=(0,),
            items=items,
            placement="main",
            columns=2,
            max_items_per_figure=2,
            image_height=r"0.40\textheight",
            figure_spec="H",
        )
    ]


def select_limit_artifacts(root: Path, warnings: list[str]) -> list[Artifact]:
    artifacts: list[Artifact] = []
    table_dir = root / "tables"
    table_patterns = ["*.table.tex", "*.pdf"]
    for path in sorted_globs(table_dir, table_patterns):
        if path.suffix.lower() == ".tex" or ".table" in path.name:
            artifact = artifact_from_path(
                path,
                caption=limit_caption_from_path(path, warnings),
                reason="limit table artifact",
            )
            if artifact is not None:
                artifacts.append(artifact)

    if artifacts:
        return dedupe_artifacts(artifacts)

    summary_candidates: list[Path] = []
    if root.exists():
        summary_candidates.extend(sorted(root.glob("*summary*.pdf")))
    return dedupe_artifacts(
        artifact
        for path in summary_candidates
        for artifact in [
            artifact_from_path(
                path,
                caption=limit_caption_from_path(path, warnings),
                reason="limit summary plot",
            )
        ]
        if artifact is not None
    )


def best_plot_from_payload(payload: dict[str, Any], root: Path) -> Path | None:
    for key in (
        "plot_pdf",
        "plot_pdf_repo_relative",
        "pdf_path",
        "pdf_path_repo_relative",
        "path",
    ):
        path = resolve_payload_path(payload.get(key), root)
        if path is not None and path.exists() and path.suffix.lower() == ".pdf":
            return path
    return None


def has_valid_pull_points(payload: dict[str, Any]) -> bool:
    if payload.get("status") != "ok":
        return False
    try:
        return int(payload.get("points", 0)) > 0
    except (TypeError, ValueError):
        return False


def iter_bias_scatter_payloads(summary: dict[str, Any]) -> Iterable[dict[str, Any]]:
    scatter_plot = summary.get("scatter_plot")
    if isinstance(scatter_plot, dict):
        yield scatter_plot
    variations = summary.get("scatter_plot_variations", {})
    if not isinstance(variations, dict):
        return
    for plots in variations.values():
        if not isinstance(plots, list):
            continue
        for plot in plots:
            if isinstance(plot, dict):
                yield plot


def valid_bias_scatter_paths(root: Path, warnings: list[str]) -> set[Path]:
    summary = read_json(root / "summary.json", warnings)
    if not summary:
        return set()
    paths: set[Path] = set()
    for plot in iter_bias_scatter_payloads(summary):
        if not has_valid_pull_points(plot):
            continue
        path = best_plot_from_payload(plot, root)
        if path is not None:
            paths.add(path.resolve())
    return paths


def is_bias_scatter_artifact(artifact: Artifact) -> bool:
    return artifact.kind == "plot" and artifact.source.name.startswith("pull_mean_sigma_scatter")


def pull_distribution_pdf_from_experiment(experiment: dict[str, Any], root: Path) -> Path | None:
    path = resolve_payload_path(experiment.get("pull_plot"), root)
    candidates: list[Path] = []
    if path is not None:
        candidates.append(path.with_suffix(".pdf") if path.suffix.lower() != ".pdf" else path)
    fit_summary = experiment.get("fit_summary_html")
    fit_summary_path = resolve_payload_path(fit_summary, root)
    if fit_summary_path is not None:
        candidates.append(fit_summary_path.parent / PULL_DISTRIBUTION_NAME)
    for candidate in candidates:
        if candidate.exists() and candidate.is_file() and candidate.suffix.lower() == ".pdf":
            return candidate.resolve()
    return None


def is_valid_pull_experiment(experiment: dict[str, Any]) -> bool:
    if experiment.get("dataset_strategy") != "toys":
        return False
    if experiment.get("status") != "ok":
        return False
    try:
        return int(experiment.get("entries_used", 0)) > 0
    except (TypeError, ValueError):
        return False


def poi_group_metadata(experiment: dict[str, Any]) -> tuple[str, str, str, tuple[Any, ...]]:
    target_poi = str(experiment.get("target_poi") or experiment.get("target_label") or "")
    target_label = str(experiment.get("target_label") or "")
    scheme = str(experiment.get("scheme") or "")
    label_text = " ".join(
        str(experiment.get(key) or "")
        for key in ("target_display_label", "scheme_label", "target_slug")
    )
    haystack = " ".join([target_poi, target_label, scheme, label_text])
    boson = "Z" if re.search(r"\bZ\b|r_Z|Z_", haystack) else "H"

    grouped = "grouped" in scheme or "combined" in label_text.lower() or target_poi.endswith("grouped")
    if grouped:
        state = "combined"
        mode = "grouped"
    else:
        mode = "separate"
        state_match = re.search(r"([123])S", haystack, flags=re.IGNORECASE)
        state = f"{state_match.group(1)}S" if state_match else "state"

    boson_order = {"Z": 0, "H": 1}.get(boson, 99)
    mode_order = {"separate": 0, "grouped": 1}.get(mode, 99)
    state_order = {"1S": 0, "2S": 1, "3S": 2, "combined": 3}.get(state, 99)
    return boson, mode, state, (boson_order, mode_order, state_order)


def poi_group_title(boson: str, mode: str, state: str) -> str:
    if mode == "grouped":
        return f"{boson} grouped Upsilon(nS) POI"
    return f"{boson} separate Upsilon({state}) POI"


def pull_grid_item_label(experiment: dict[str, Any]) -> str:
    truth_pdf_index = experiment.get("truth_pdf_index")
    truth_pdf = (
        f"pdf {truth_pdf_index}"
        if truth_pdf_index is not None
        else str(experiment.get("truth_pdf_slug") or "truth PDF")
    )
    fit_pdf = str(experiment.get("fit_pdf_mode") or experiment.get("target_pdf_slug") or "fit")
    parts = [
        f"r={format_number(experiment.get('injected_r', '-'))}",
        truth_pdf,
    ]
    if fit_pdf:
        parts.append(fit_pdf)
    mean = experiment.get("gaussian_mean")
    sigma = experiment.get("gaussian_sigma")
    if experiment.get("gaussian_status") == "ok" and mean is not None and sigma is not None:
        parts.append(f"mu={format_number(mean)}, sig={format_number(sigma)}")
    elif experiment.get("gaussian_status"):
        parts.append(f"pull summary: {experiment['gaussian_status']}")
    return "; ".join(parts)


def pull_grid_sort_key(experiment: dict[str, Any]) -> tuple[Any, ...]:
    return (
        float(experiment.get("injected_r", 0.0) or 0.0),
        int(experiment.get("truth_pdf_index", 999) if experiment.get("truth_pdf_index") is not None else 999),
        str(experiment.get("fit_pdf_mode") or experiment.get("target_pdf_slug") or ""),
        int(experiment.get("experiment_index", 0) or 0),
    )


def build_bias_pull_grid_groups(root: Path, warnings: list[str]) -> list[GridGroup]:
    summary = read_json(root / "summary.json", warnings)
    if not summary:
        return []
    strategies = [str(item) for item in (summary.get("dataset_strategies") or [summary.get("dataset_strategy")])]
    if "toys" not in strategies:
        return []
    groups: dict[tuple[str, str, str], dict[str, Any]] = {}
    experiments = summary.get("experiments", [])
    if not isinstance(experiments, list):
        return []
    for experiment in experiments:
        if not isinstance(experiment, dict) or not is_valid_pull_experiment(experiment):
            continue
        path = pull_distribution_pdf_from_experiment(experiment, root)
        if path is None:
            continue
        boson, mode, state, sort_key = poi_group_metadata(experiment)
        group_key = (boson, mode, state)
        group_payload = groups.setdefault(
            group_key,
            {
                "title": poi_group_title(boson, mode, state),
                "sort_key": sort_key,
                "items": [],
            },
        )
        artifact = artifact_from_path(
            path,
            caption="Toy pull distribution",
            reason="bias-study pull distribution",
        )
        if artifact is None:
            continue
        group_payload["items"].append(
            GridItem(
                artifact=artifact,
                label=pull_grid_item_label(experiment),
                sort_key=pull_grid_sort_key(experiment),
            )
        )
    grid_groups: list[GridGroup] = []
    for payload in groups.values():
        items = sorted(payload["items"], key=lambda item: item.sort_key)
        if not items:
            continue
        title = str(payload["title"])
        grid_groups.append(
            GridGroup(
                title=title,
                caption=(
                    f"Toy pull distributions for {title}. Each panel shows the fitted pull distribution for one "
                    "truth-PDF, injected-signal, and fit-PDF configuration."
                ),
                sort_key=tuple(payload["sort_key"]),
                items=items,
            )
        )
    return sorted(grid_groups, key=lambda group: group.sort_key)


def is_bias_heatmap_path(path: Path) -> bool:
    return "pdf_target_heatmap" in path.name or "heatmap" in path.name.lower()


def is_all_tested_heatmap(plot: dict[str, Any], path: Path | None = None) -> bool:
    fields = [
        str(plot.get("label") or ""),
        str(plot.get("plot_stem") or ""),
    ]
    if path is not None:
        fields.append(path.stem)
    haystack = " ".join(fields).lower().replace("_", " ")
    return "all tested combinations" in haystack


def build_bias_heatmap_grid_groups(root: Path, warnings: list[str]) -> list[GridGroup]:
    summary = read_json(root / "summary.json", warnings)
    if not summary:
        return []
    heatmap_plots = summary.get("heatmap_plots", [])
    if not isinstance(heatmap_plots, list):
        return []

    items: list[GridItem] = []
    for index, plot in enumerate(heatmap_plots):
        if not isinstance(plot, dict) or plot.get("status", "ok") != "ok":
            continue
        path = best_plot_from_payload(plot, root)
        if path is None or is_all_tested_heatmap(plot, path):
            continue
        artifact = artifact_from_path(
            path,
            caption=str(plot.get("label") or title_from_stem(path.stem)),
            reason="bias-study heatmap grid plot",
        )
        if artifact is None:
            continue
        items.append(
            GridItem(
                artifact=artifact,
                label=str(plot.get("label") or title_from_stem(path.stem)),
                sort_key=(index,),
            )
        )
    if not items:
        return []
    return [
        GridGroup(
            title="Bias Study Heatmaps",
            caption=(
                "PDF-target heatmaps for the bias study, grouped four per page. Each panel summarizes "
                "the fitted pull mean for one POI target and injected signal strength; the global aggregate "
                "heatmap is skipped"
            ),
            sort_key=(-1,),
            items=items,
            placement="appendix",
            columns=2,
            max_items_per_figure=4,
            image_height=r"0.31\textheight",
            figure_spec="p",
        )
    ]


def build_bias_scatter_grid_groups(root: Path, warnings: list[str]) -> list[GridGroup]:
    paths = sorted(valid_bias_scatter_paths(root, warnings))
    if not paths:
        return []
    global_item = None
    appendix_items: list[GridItem] = []
    for path in paths:
        artifact = artifact_from_path(
            path,
            caption=generic_caption_from_path(path),
            reason="bias-study scatter grid plot",
        )
        if artifact is None:
            continue
        label = title_from_stem(path.stem)
        item = GridItem(artifact=artifact, label=label, sort_key=(len(appendix_items),))
        if path.stem == "pull_mean_sigma_scatter":
            global_item = item
        else:
            appendix_items.append(item)
    groups: list[GridGroup] = []
    if global_item is not None:
        groups.append(
            GridGroup(
                title="Bias Study Pull Scatter",
                caption="Global pull mean and width summary across all valid toy pull summaries",
                sort_key=(0,),
                items=[global_item],
                placement="main",
                columns=1,
                max_items_per_figure=1,
                image_height=r"0.78\textheight",
                figure_spec="H",
            )
        )
    if appendix_items:
        groups.append(
            GridGroup(
                title="Bias Study Pull Scatter Plots",
                caption="Pull mean and width scatter plots from the bias study, two per page",
                sort_key=(1,),
                items=appendix_items,
                placement="appendix",
                columns=1,
                max_items_per_figure=2,
                image_height=r"0.40\textheight",
                figure_spec="p",
            )
        )
    return groups


def build_bias_grid_groups(root: Path, warnings: list[str]) -> list[GridGroup]:
    return [
        *build_bias_scatter_grid_groups(root, warnings),
        *build_bias_heatmap_grid_groups(root, warnings),
        *build_bias_pull_grid_groups(root, warnings),
    ]


def filter_bias_appendix_artifacts(
    root: Path,
    artifacts: list[Artifact],
    grid_groups: list[GridGroup],
    warnings: list[str],
) -> list[Artifact]:
    valid_scatter_paths = valid_bias_scatter_paths(root, warnings)
    grid_artifact_keys = {
        artifact_group_key(item.artifact)
        for group in grid_groups
        for item in group.items
    }
    filtered: list[Artifact] = []
    for artifact in artifacts:
        if artifact.source.name == PULL_DISTRIBUTION_NAME:
            continue
        if is_bias_heatmap_path(artifact.source):
            continue
        if is_bias_scatter_artifact(artifact):
            continue
        if artifact_group_key(artifact) in grid_artifact_keys:
            continue
        filtered.append(artifact)
    return filtered


def filter_statistical_model_appendix_artifacts(artifacts: list[Artifact]) -> list[Artifact]:
    return [
        artifact
        for artifact in artifacts
        if not (
            artifact.kind == "plot"
            and artifact.source.name.startswith("bundle_model_projection_")
        )
    ]


def select_bias_artifacts(root: Path, warnings: list[str]) -> list[Artifact]:
    artifacts: list[Artifact] = []
    summary = read_json(root / "summary.json", warnings)
    if summary:
        scatter_plot = summary.get("scatter_plot")
        if isinstance(scatter_plot, dict) and has_valid_pull_points(scatter_plot):
            path = best_plot_from_payload(scatter_plot, root)
            if path is not None:
                artifact = artifact_from_path(
                    path,
                    caption=(
                        "Global pull mean and width summary across valid toy pull summaries"
                    ),
                    reason="bias-study global scatter plot",
                )
                if artifact is not None:
                    artifacts.append(artifact)

    if artifacts:
        return dedupe_artifacts(artifacts)

    fallback_patterns = [
        "plots/pull_mean_sigma_scatter.pdf",
    ]
    return dedupe_artifacts(
        artifact
        for path in sorted_globs(root, fallback_patterns)
        for artifact in [artifact_from_path(path, reason="bias-study summary plot")]
        if artifact is not None
    )


def make_section(
    *,
    key: str,
    title: str,
    source_dir: Path,
    output_dir: Path,
    selector: Any,
    grid_builder: Any | None = None,
    description_tex: str = "",
) -> Section:
    warnings: list[str] = []
    notes: list[str] = []
    if not source_dir.exists() or not source_dir.is_dir():
        warnings.append(f"Input directory does not exist: {repo_relative(source_dir)}")
        return Section(key, title, source_dir, [], [], [], notes, warnings)

    main_artifacts = dedupe_artifacts(selector(source_dir, warnings))
    grid_groups = grid_builder(source_dir, warnings) if grid_builder is not None else []
    main_keys = {artifact_group_key(artifact) for artifact in main_artifacts}
    grid_keys = {
        artifact_group_key(item.artifact)
        for group in grid_groups
        for item in group.items
    }
    main_artifacts = [
        artifact
        for artifact in main_artifacts
        if artifact_group_key(artifact) not in grid_keys
    ]
    main_grid_groups = [group for group in grid_groups if group.placement == "main"]
    appendix_grid_groups = [group for group in grid_groups if group.placement != "main"]
    appendix_artifacts: list[Artifact] = []
    if not main_artifacts and not main_grid_groups:
        notes.append("No summary-level supported plots or table fragments were found.")
    if not appendix_artifacts and not appendix_grid_groups:
        notes.append("No additional supported appendix plots or table fragments were found.")
    return Section(key, title, source_dir, main_artifacts, appendix_artifacts, grid_groups, notes, warnings, description_tex)


def latest_limits_run(warnings: list[str]) -> Path:
    summaries = sorted(
        (REPO_ROOT / "limits").glob("*/blind_limits_summary.json"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if summaries:
        return summaries[0].parent.resolve()
    warnings.append("No default limits run found under limits/*/blind_limits_summary.json.")
    return (REPO_ROOT / "limits" / "missing_limits_run").resolve()


def latest_bias_run(strategy: str, warnings: list[str]) -> Path:
    candidates: list[Path] = []
    for summary_path in sorted((REPO_ROOT / "bias_study").glob("*/summary.json")):
        try:
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        strategies = payload.get("dataset_strategies") or [payload.get("dataset_strategy")]
        if strategy in [str(item) for item in strategies if item is not None]:
            candidates.append(summary_path)
    candidates.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    if candidates:
        return candidates[0].parent.resolve()
    warnings.append(f"No default {strategy} bias-study run found under bias_study/*/summary.json.")
    return (REPO_ROOT / "bias_study" / f"missing_{strategy}_run").resolve()


MAX_ASSET_NAME_LENGTH = 120


def asset_name(source: Path, source_root: Path, used_names: set[str]) -> str:
    try:
        rel = source.resolve().relative_to(source_root.resolve())
        base = "__".join(rel.parts)
    except ValueError:
        base = source.name
    suffix = "".join(source.suffixes) if source.suffixes else source.suffix
    suffix_lower = suffix.lower()
    suffix_len = len(suffix_lower)
    max_stem_len = MAX_ASSET_NAME_LENGTH - suffix_len
    if suffix and base.endswith(suffix):
        raw_stem = base[: -len(suffix)]
    else:
        raw_stem = source.stem
    stem = slugify(raw_stem)
    if len(stem) > max_stem_len:
        import hashlib
        fingerprint = hashlib.shake_128(raw_stem.encode(), usedforsecurity=False).hexdigest(4)
        keep = max(8, max_stem_len - len(fingerprint) - 1)
        stem = f"{stem[:keep]}_{fingerprint}"
    candidate = f"{stem}{suffix_lower}"
    if candidate not in used_names:
        used_names.add(candidate)
        return candidate
    counter = 2
    while True:
        numbered = f"{slugify(raw_stem)}_{counter}{suffix_lower}"
        if len(numbered) > MAX_ASSET_NAME_LENGTH:
            numbered = f"{stem}_{counter}{suffix_lower}"
        if numbered not in used_names:
            used_names.add(numbered)
            return numbered
        counter += 1


MAX_FILES_PER_DIR = 150


class AssetDir:
    def __init__(self, output_dir: Path):
        self._output_dir = output_dir
        self._counter = 0
        self._dir_index = 1
        self._current_dir: Path | None = None

    def next_path(self, name: str) -> Path:
        if self._current_dir is None or self._counter >= MAX_FILES_PER_DIR:
            self._dir_index += 1 if self._current_dir is not None else 0
            self._current_dir = self._output_dir / f"assets_summary_docs_{self._dir_index}"
            self._current_dir.mkdir(parents=True, exist_ok=True)
            self._counter = 0
        self._counter += 1
        return self._current_dir / name


def materialize_section(section: Section, output_dir: Path, asset_dir: AssetDir) -> RenderedSection:
    used_names: set[str] = set()

    def materialize(artifact: Artifact) -> RenderedArtifact:
        name = asset_name(artifact.source, section.source_dir, used_names)
        target = asset_dir.next_path(name)
        shutil.copy2(artifact.source, target)
        return RenderedArtifact(
            source=artifact.source,
            kind=artifact.kind,
            caption=artifact.caption,
            reason=artifact.reason,
            asset_relative=str(target.relative_to(output_dir)),
        )

    def materialize_grid_group(group: GridGroup) -> RenderedGridGroup:
        return RenderedGridGroup(
            title=group.title,
            caption=group.caption,
            sort_key=group.sort_key,
            items=[
                RenderedGridItem(
                    artifact=materialize(item.artifact),
                    label=item.label,
                    sort_key=item.sort_key,
                )
                for item in group.items
            ],
            placement=group.placement,
            columns=group.columns,
            max_items_per_figure=group.max_items_per_figure,
            image_height=group.image_height,
            figure_spec=group.figure_spec,
        )

    return RenderedSection(
        key=section.key,
        title=section.title,
        source_dir=section.source_dir,
        main_artifacts=[materialize(artifact) for artifact in section.main_artifacts],
        appendix_artifacts=[materialize(artifact) for artifact in section.appendix_artifacts],
        grid_groups=[materialize_grid_group(group) for group in section.grid_groups],
        notes=section.notes,
        warnings=section.warnings,
        description_tex=section.description_tex,
    )


def sentence_text(value: str) -> str:
    value = value.strip()
    if not value or value.endswith(('.', '!', '?')):
        return value
    return f"{value}."


def table_heading_from_caption(caption: str) -> str:
    heading = re.split(r"[.;:]", caption, maxsplit=1)[0].strip()
    return heading or "Table"


def render_artifact(artifact: RenderedArtifact) -> str:
    caption = latex_escape_caption(sentence_text(artifact.caption))
    if artifact.kind == "plot":
        return rf"""
\begin{{figure}}[H]
\centering
\includegraphics[width=0.96\textwidth,height=0.78\textheight,keepaspectratio]{{{latex_path(artifact.asset_relative)}}}
\caption{{{caption}}}
\end{{figure}}
\FloatBarrier
"""
    if artifact.kind == "table":
        heading = latex_escape_caption(table_heading_from_caption(artifact.caption))
        description = latex_escape_caption(sentence_text(artifact.caption))
        return rf"""
\subsection*{{{heading}}}
\noindent\textit{{{description}}}
\par\medskip
\input{{{latex_path(artifact.asset_relative)}}}
\FloatBarrier
"""
    if artifact.kind == "text":
        return rf"""
\subsection*{{{caption}}}
\VerbatimInput[fontsize=\scriptsize,breaklines=true,breakanywhere=true]{{{latex_path(artifact.asset_relative)}}}
\FloatBarrier
"""
    raise ValueError(f"Unknown artifact kind: {artifact.kind}")


def chunks(items: list[RenderedGridItem], size: int) -> Iterable[list[RenderedGridItem]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def render_grid_group(group: RenderedGridGroup) -> str:
    columns = max(int(group.columns), 1)
    max_items = max(int(group.max_items_per_figure), 1)
    cell_width = max(0.1, min(0.96 / columns - 0.01, 0.96))
    body = [rf"\subsection{{{latex_escape_caption(group.title)}}}"]
    for index, group_items in enumerate(chunks(sorted(group.items, key=lambda item: item.sort_key), max_items), start=1):
        if len(group_items) == 2 and columns == 2:
            grid_tex = (
                rf"\includegraphics[width=0.45\textwidth]{{{latex_path(group_items[0].artifact.asset_relative)}}}"
                + rf"\hspace*{{1.cm}}"
                + rf"\includegraphics[width=0.45\textwidth]{{{latex_path(group_items[1].artifact.asset_relative)}}}"
            )
        else:
            row_tex = []
            for row_items in chunks(group_items, columns):
                row_cells = []
                for item in row_items:
                    row_cells.append(
                        rf"""\begin{{subfigure}}[t]{{{cell_width:.2f}\textwidth}}
\centering
\includegraphics[width=\linewidth,height={group.image_height},keepaspectratio]{{{latex_path(item.artifact.asset_relative)}}}
\caption{{{latex_escape_caption(item.label)}}}
\end{{subfigure}}"""
                    )
                row_tex.append("%\n\\hfill\n".join(row_cells))
            grid_tex = "\n\\par\\medskip\n".join(row_tex)
        caption_text = group.caption if index == 1 else f"{group.caption} (continued {index})"
        body.append(
            rf"""
\begin{{figure}}[{group.figure_spec}]
\begin{{center}}
{grid_tex}
\end{{center}}
\caption{{{latex_escape_caption(sentence_text(caption_text))}}}
\end{{figure}}
\FloatBarrier
"""
        )
    return "\n".join(body)


def render_notes(section: RenderedSection) -> str:
    rows = []
    for note in section.notes:
        rows.append(rf"\item {latex_escape_text(note)}")
    for warning in section.warnings:
        rows.append(rf"\item Warning: {latex_escape_text(warning)}")
    return "\n".join(rows)


def render_section(section: RenderedSection) -> str:
    body = [rf"\section{{{latex_escape_text(section.title)}}}"]
    if section.description_tex:
        body.append(section.description_tex)
    main_grid_groups = [group for group in section.grid_groups if group.placement == "main"]
    notes = render_notes(section)
    if notes:
        body.append(r"\begin{itemize}")
        body.append(notes)
        body.append(r"\end{itemize}")
    if section.main_artifacts:
        body.extend(render_artifact(artifact) for artifact in section.main_artifacts)
    if main_grid_groups:
        body.extend(render_grid_group(group) for group in main_grid_groups)
    if not section.main_artifacts and not main_grid_groups:
        body.append(r"\noindent Placeholder for this section's narrative and summary artifacts.")
    return "\n".join(body)


def render_appendix_section(section: RenderedSection) -> str:
    body = [rf"\section{{{latex_escape_text(section.title)}}}"]
    body.append(r"\noindent Additional plots and table fragments for this analysis section.")
    appendix_grid_groups = [group for group in section.grid_groups if group.placement != "main"]
    if appendix_grid_groups:
        body.append(r"\subsection*{Additional Plot Grids}")
        body.extend(render_grid_group(group) for group in appendix_grid_groups)
    if section.appendix_artifacts:
        body.extend(render_artifact(artifact) for artifact in section.appendix_artifacts)
    elif not appendix_grid_groups:
        body.append(r"\par\noindent No additional supported plots or table fragments were found for this section.")
    return "\n".join(body)


def build_latex_document(sections: list[RenderedSection], generated_at: str) -> str:
    section_tex = "\n\n".join(render_section(section) for section in sections)
    appendix_tex = "\n\n\\clearpage\n\n".join(render_appendix_section(section) for section in sections)
    return rf"""\documentclass[11pt]{{article}}
\usepackage[margin=0.75in]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{float}}
\usepackage{{placeins}}
\usepackage{{booktabs}}
\usepackage{{longtable}}
\usepackage{{array}}
\usepackage{{caption}}
\usepackage{{subcaption}}
\captionsetup[subfigure]{{font=scriptsize}}
\usepackage{{fvextra}}
\usepackage[hidelinks]{{hyperref}}

\title{{Statistical Test Fit Summary}}
\author{{}}
\date{{Generated {latex_escape_text(generated_at)}}}

\begin{{document}}
\maketitle
\tableofcontents
\clearpage

{section_tex}

\clearpage
\appendix
\section*{{Appendix}}
\addcontentsline{{toc}}{{section}}{{Appendix}}
The appendix contains supported plots and table fragments found in each input directory that were not selected for the main summary sections.

{appendix_tex}

\end{{document}}
"""


def write_metadata(
    output_dir: Path,
    output_name: str,
    sections: list[RenderedSection],
    generated_at: str,
    document_path: Path,
    global_warnings: list[str],
) -> Path:
    metadata_path = output_dir / f"{output_name}.json"
    payload = {
        "schema_version": 1,
        "generated_at": generated_at,
        "document": str(document_path.resolve()),
        "document_repo_relative": repo_relative(document_path),
        "sections": [
            {
                "key": section.key,
                "title": section.title,
                "source_dir": str(section.source_dir),
                "source_dir_repo_relative": repo_relative(section.source_dir),
                "notes": section.notes,
                "warnings": section.warnings,
                "main_artifacts": [artifact_metadata(artifact) for artifact in section.main_artifacts],
                "appendix_artifacts": [artifact_metadata(artifact) for artifact in section.appendix_artifacts],
                "grid_groups": [grid_group_metadata(group) for group in section.grid_groups],
            }
            for section in sections
        ],
        "warnings": global_warnings,
    }
    metadata_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return metadata_path


def artifact_metadata(artifact: RenderedArtifact) -> dict[str, Any]:
    return {
        "source": str(artifact.source),
        "source_repo_relative": repo_relative(artifact.source),
        "kind": artifact.kind,
        "caption": artifact.caption,
        "reason": artifact.reason,
        "asset_relative": artifact.asset_relative,
    }


def grid_group_metadata(group: RenderedGridGroup) -> dict[str, Any]:
    return {
        "title": group.title,
        "caption": group.caption,
        "placement": group.placement,
        "columns": group.columns,
        "max_items_per_figure": group.max_items_per_figure,
        "image_height": group.image_height,
        "figure_spec": group.figure_spec,
        "items": [
            {
                "label": item.label,
                "artifact": artifact_metadata(item.artifact),
            }
            for item in group.items
        ],
    }


def compile_with_apptainer(
    document_path: Path,
    image: str,
    apptainer_bin: str,
) -> Path:
    output_dir = document_path.parent.resolve()
    command = [
        apptainer_bin,
        "exec",
        "--cleanenv",
        "--bind",
        f"{output_dir}:/work",
        "--pwd",
        "/work",
        image,
        "latexmk",
        "-pdf",
        "-interaction=nonstopmode",
        "-halt-on-error",
        document_path.name,
    ]
    stdout_path = output_dir / "latex_compile_stdout.txt"
    stderr_path = output_dir / "latex_compile_stderr.txt"
    try:
        completed = subprocess.run(command, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Could not execute {apptainer_bin!r}. Install Apptainer or pass --skip-compile-pdf."
        ) from exc
    stdout_path.write_text(completed.stdout or "", encoding="utf-8")
    stderr_path.write_text(completed.stderr or "", encoding="utf-8")
    compile_json = output_dir / "latex_compile.json"
    compile_json.write_text(
        json.dumps(
            {
                "command": command,
                "command_string": shlex.join(command),
                "returncode": completed.returncode,
                "stdout": str(stdout_path.resolve()),
                "stderr": str(stderr_path.resolve()),
                "image": image,
                "engine": apptainer_bin,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"LaTeX compilation failed with return code {completed.returncode}. See {stdout_path} and {stderr_path}."
        )
    pdf_path = output_dir / f"{document_path.stem}.pdf"
    if not pdf_path.exists():
        raise RuntimeError(f"LaTeX command succeeded but PDF was not found: {pdf_path}")
    return pdf_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a LaTeX/PDF summary document from produced analysis plots and tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--signal-modeling",
        default="plots/signal_fit",
        help="Directory with signal-modeling plots.",
    )
    parser.add_argument(
        "--resonant-background-modeling",
        default="plots/resonant_background",
        help="Directory with resonant-background plots.",
    )
    parser.add_argument(
        "--non-resonant-background-model",
        default="plots/fit_2d_data",
        help="Directory with non-resonant-background plots and summary JSON.",
    )
    parser.add_argument(
        "--statistical-model",
        default="datacards",
        help="Directory with bundled workspace summary and model plots.",
    )
    parser.add_argument(
        "--limits",
        default=None,
        help="Limits run directory. Defaults to the newest limits/* run with blind_limits_summary.json.",
    )
    parser.add_argument(
        "--bias-study-toys",
        default=None,
        help="Bias-study run directory for the toys dataset strategy. Defaults to the newest matching bias_study/* run.",
    )
    parser.add_argument(
        "--output-dir",
        default="summary-docs",
        help="Output directory for the LaTeX document, staged assets/inputs, metadata, and optional PDF.",
    )
    parser.add_argument(
        "--output-name",
        default=DEFAULT_OUTPUT_NAME,
        help="Base filename for .tex, .json, and optional .pdf outputs.",
    )
    parser.add_argument(
        "--skip-compile-pdf",
        action="store_true",
        help="Write LaTeX and metadata only; do not compile the PDF.",
    )
    parser.add_argument(
        "--latex-image",
        default=DEFAULT_IMAGE,
        help="Docker image URI passed to Apptainer.",
    )
    parser.add_argument(
        "--apptainer-bin",
        default="apptainer",
        help="Apptainer executable name or path. Pass 'singularity' if that is the available compatible runtime.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    output_dir = resolve_input_path(args.output_dir)
    if not args.skip_compile_pdf and not args.latex_image.startswith("docker://"):
        raise ValueError("--latex-image must be a docker:// URI when compiling through Apptainer")

    global_warnings: list[str] = []
    limits_dir = resolve_input_path(args.limits) if args.limits else latest_limits_run(global_warnings)
    bias_toys_dir = (
        resolve_input_path(args.bias_study_toys)
        if args.bias_study_toys
        else latest_bias_run("toys", global_warnings)
    )

    sections = [
        make_section(
            key="signal-modeling",
            title="Signal Modeling",
            source_dir=resolve_input_path(args.signal_modeling),
            output_dir=output_dir,
            selector=select_signal_artifacts,
            grid_builder=build_signal_grid_groups,
            description_tex=(
                r"The signal processes $H\to\Upsilon(nS)\gamma$ and $Z\to\Upsilon(nS)\gamma$ are "
                r"modeled using simulated Monte Carlo samples. For each signal process and Upsilon state "
                r"($1S$, $2S$, $3S$), parametric shapes are extracted from fits to the dimuon invariant "
                r"mass $m_{\mu\mu}$ and the boson invariant mass $m_{\mu\mu\gamma}$ distributions. "
                r"The $m_{\mu\mu}$ distribution is parameterized with a Double Crystal Ball (DCB) "
                r"function to capture the detector resolution and the radiative tail. The "
                r"$m_{\mu\mu\gamma}$ distribution is modeled with a Crystal Ball convolved with a "
                r"Gaussian resolution function. \\[4pt]"
                r"\noindent Signal samples are generated at next-to-leading order (NLO) in QCD: "
                r"$Z\to\Upsilon\gamma$ events are simulated with \textsc{MadGraph5\_aMC@NLO}, and "
                r"$gg\to H\to\Upsilon\gamma$ events with \textsc{Powheg} + \textsc{Pythia8}, "
                r"using the NNPDF\,3.1 parton distribution functions and the CUETP8M1 tune for the "
                r"underlying event. All signal shape parameters are extracted from the simulated "
                r"samples and fixed in the final simultaneous fit to collision data."
            ),
        ),
        make_section(
            key="resonant-background-modeling",
            title="Resonant Background Modeling",
            source_dir=resolve_input_path(args.resonant_background_modeling),
            output_dir=output_dir,
            selector=select_resonant_artifacts,
            grid_builder=build_resonant_grid_groups,
            description_tex=(
                r"Resonant backgrounds arise from standard-model processes that produce a genuine "
                r"$\Upsilon(nS)$ meson in association with a photon, but whose production mechanism "
                r"differs from the Higgs or Z boson decays under study. The dominant resonant "
                r"backgrounds are $Z\gamma$ production where the Z boson decays to a muon pair "
                r"($Z\to\mu\mu\gamma$, with a small tail from $\tau\tau\gamma$), and Higgs boson "
                r"Dalitz decays $H\to\Upsilon\gamma$ via an intermediate virtual photon. "
                r"\\[4pt]"
                r"\noindent The resonant $Z$ background shape is modeled using a parametric fit to "
                r"simulated $Z\gamma$ events in the $m_{\mu\mu\gamma}$ and $m_{\mu\mu}$ distributions, "
                r"stored as a RooKeysPdf template in the workspace. The Higgs Dalitz background is "
                r"similarly modeled from simulated $H\to\Upsilon\gamma$ samples and stored as a "
                r"RooKeysPdf template. The $Z$ boson plus photon ($Z\gamma$) contribution provides "
                r"the largest resonant background, while the Dalitz contribution is subdominant but "
                r"non-negligible in the H boson mass window. "
                r"\\[4pt]"
                r"\noindent Since both resonant backgrounds peak in the same boson-mass region "
                r"as the signal, their normalizations are constrained using four dedicated control "
                r"regions (CR1--CR4) enriched in $Z\to\mu\mu$ events. These control regions are "
                r"defined in the $m_{\mu\mu}$--$m_{\mu\mu\gamma}$ plane around the Z boson peak, "
                r"selecting dimuon candidates compatible with the Z resonance but with varying "
                r"photon-quality or kinematic requirements, providing orthogonal samples that map "
                r"the resonant background composition across the signal-selection phase space. "
                r"A simultaneous fit to the observed data and the simulated $Z\gamma$ template "
                r"in CR1--CR4 determines the normalization transfer factors that extrapolate "
                r"the resonant background yields from the control regions into the signal region. "
                r"The resonant background shape parameters are fixed from simulation and "
                r"only the overall normalizations are allowed to float, subject to the "
                r"control-region constraints, in the final simultaneous fit to data."
            ),
        ),
        make_section(
            key="non-resonant-background-model",
            title="Non-Resonant Background Model",
            source_dir=resolve_input_path(args.non_resonant_background_model),
            output_dir=output_dir,
            selector=select_non_resonant_artifacts,
            grid_builder=build_non_resonant_grid_groups,
            description_tex=(
                r"The non-resonant background comprises all standard-model processes that do not "
                r"produce a genuine $\Upsilon(nS)$ meson, primarily continuum quarkonium production, "
                r"Drell--Yan dimuon events, and multijet QCD processes where jets are misidentified "
                r"as photons. This background is the dominant component in the analysis and is "
                r"modeled directly from collision data using a two-dimensional (2D) parametric "
                r"function in the $m_{\mu\mu}$--$m_{\mu\mu\gamma}$ plane. "
                r"\\[4pt]"
                r"\noindent The modeling strategy proceeds in two stages. First, a dedicated fit to "
                r"the dimuon invariant mass distribution in a sideband region away from the H and Z "
                r"boson mass windows determines the $\Upsilon(nS)$ resonance component of the "
                r"background: three Crystal Ball functions (one per $\Upsilon(1S,2S,3S)$ state) on "
                r"top of a smooth continuum. The fitted $\Upsilon(nS)$ yields and resolution "
                r"parameters are then propagated into the 2D background model. "
                r"\\[4pt]"
                r"\noindent In the second stage, the 2D sideband distribution is described by a "
                r"factorized parametric function"
                r"\\[4pt]"
                r"\begin{equation}"
                r"B(m_{\mu\mu},\,m_{\mu\mu\gamma}) = P(m_{\mu\mu},\,m_{\mu\mu\gamma}) "
                r"+ \sum_{n=1}^{3} \mathrm{CB}_n(m_{\mu\mu})\,G_n(m_{\mu\mu\gamma})"
                r"\end{equation}"
                r"\\[4pt]"
                r"\noindent where $P$ is a smooth 2D polynomial describing the continuum and $\mathrm{CB}_n$ "
                r"($G_n$) is a Crystal Ball (Gaussian) function for the $n$-th Upsilon state in the "
                r"dimuon (boson) mass projection. Three families of 2D polynomial functions are "
                r"tested as candidates for $P$: the Johnson $S_U$ distribution, Bernstein "
                r"polynomials, and Chebyshev polynomials of the first kind. "
                r"\\[4pt]"
                r"\noindent Within each polynomial family, the optimal order is determined by "
                r"fitting candidate models of increasing polynomial degree to the 2D sideband "
                r"data. Each candidate fit must pass a set of quality checks before being "
                r"considered further: the fit must converge (status $=0$), the covariance matrix "
                r"must be accurately computed (covQual $\geq 3$), the estimated distance to the "
                r"minimum must satisfy EDM $\leq 10^{-3}$, and the negative log-likelihood "
                r"must be finite. "
                r"\\[4pt]"
                r"\noindent A two-step procedure selects the optimal order. First, the "
                r"goodness-of-fit of each candidate is assessed via a binned 2D chi-squared "
                r"test on the $m_{\mu\mu\gamma}$ projection over the sideband regions, yielding "
                r"a $\chi^2$ and its associated $p$-value. Starting from the lowest tested "
                r"order, the first model whose $\chi^2$ $p$-value exceeds 0.05 is designated "
                r"as the starting point. Second, a likelihood ratio test (LRT) is performed "
                r"between consecutive surviving models. For models $M_d$ and $M_{d+1}$ of "
                r"successive orders $d$ and $d+1$, the test statistic"
                r"\\[4pt]"
                r"\begin{equation}"
                r"q = 2\left[\mathrm{NLL}(M_d) - \mathrm{NLL}(M_{d+1})\right]"
                r"\end{equation}"
                r"\\[4pt]"
                r"\noindent is distributed as a $\chi^2$ with "
                r"$\Delta k = k_{d+1} - k_d$ degrees of freedom under the null hypothesis "
                r"that $M_d$ is sufficient. If the LRT $p$-value exceeds 0.05, the simpler "
                r"model $M_d$ is selected as the winner; otherwise the procedure advances to "
                r"the next pair, and the most complex model becomes the winner if all "
                r"successive tests reject the simpler alternative. The Johnson $S_U$, "
                r"Bernstein, and Chebyshev families are independently optimized following "
                r"this procedure, with a relaxed fallback that selects the best available "
                r"model should no candidate pass the chi-squared threshold."
                r"\\[4pt]"
                r"\noindent The complete set of non-resonant background PDFs carried into the final "
                r"simultaneous fit comprises seven configurations spanning three families: the "
                r"Johnson $S_U$ (nominal model); the Bernstein polynomial of order~7 (selected "
                r"model) flanked by orders~6 and~8 (below and above the selected order); and the "
                r"Chebyshev polynomial of order~5 (selected model) flanked by orders~4 and~6 "
                r"(below and above). The Johnson $S_U$ serves as the nominal background "
                r"parameterization, while the Bernstein and Chebyshev families provide systematic "
                r"envelope variations that probe the dependence of the result on the choice of "
                r"background functional form. The non-resonant background normalization is left "
                r"floating in the final simultaneous fit to data."
            ),
        ),
        make_section(
            key="statistical-model",
            title="Statistical Model",
            source_dir=resolve_input_path(args.statistical_model),
            output_dir=output_dir,
            selector=select_statistical_model_artifacts,
            grid_builder=build_statistical_model_grid_groups,
        ),
        make_section(
            key="limits",
            title="Limits",
            source_dir=limits_dir,
            output_dir=output_dir,
            selector=select_limit_artifacts,
        ),
        make_section(
            key="bias-study-toys",
            title="Bias Study - Toys",
            source_dir=bias_toys_dir,
            output_dir=output_dir,
            selector=select_bias_artifacts,
            grid_builder=build_bias_grid_groups,
        ),
    ]

    output_dir.mkdir(parents=True, exist_ok=True)
    for stale in sorted(output_dir.glob("assets_summary_docs_*")):
        if stale.is_dir():
            shutil.rmtree(stale)
    asset_dir = AssetDir(output_dir)
    rendered_sections = [materialize_section(section, output_dir, asset_dir) for section in sections]

    generated_at = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    document_path = output_dir / f"{args.output_name}.tex"
    document_path.write_text(build_latex_document(rendered_sections, generated_at), encoding="utf-8")
    metadata_path = write_metadata(
        output_dir,
        args.output_name,
        rendered_sections,
        generated_at,
        document_path,
        global_warnings,
    )

    log(f"Wrote LaTeX: {repo_relative(document_path)}")
    log(f"Wrote metadata JSON: {repo_relative(metadata_path)}")
    for warning in global_warnings:
        log(f"warning: {warning}")
    for section in rendered_sections:
        log(
            f"{section.title}: {len(section.main_artifacts)} main artifact(s), "
            f"{len(section.appendix_artifacts)} appendix artifact(s)"
        )
        for warning in section.warnings:
            log(f"warning: {section.title}: {warning}")

    if not args.skip_compile_pdf:
        pdf_path = compile_with_apptainer(
            document_path,
            image=args.latex_image,
            apptainer_bin=args.apptainer_bin,
        )
        log(f"Wrote PDF: {repo_relative(pdf_path)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
