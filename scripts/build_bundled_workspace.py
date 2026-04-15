#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from statistical_test_fit.root_runtime import configure_root

ROOT = configure_root()

from ROOT import (  # type: ignore
    RooArgList,  # type: ignore
    RooArgSet,  # type: ignore
    RooCategory,  # type: ignore
    RooDataSet,  # type: ignore
    RooFit,  # type: ignore
    RooMultiPdf,  # type: ignore
    RooRealVar,  # type: ignore
    RooWorkspace,  # type: ignore
    TFile,  # type: ignore
)

from statistical_test_fit.mass_ranges import (
    BOSON_MASS_LOWER,
    BOSON_MASS_UPPER,
    LEFT_SIDEBAND_LOWER,
    LEFT_SIDEBAND_UPPER,
    MIDDLE_SIDEBAND_LOWER,
    MIDDLE_SIDEBAND_UPPER,
    RIGHT_SIDEBAND_LOWER,
    RIGHT_SIDEBAND_UPPER,
    UPSILON_MASS_LOWER,
    UPSILON_MASS_UPPER,
)
from statistical_test_fit.signal_modeling import SIGNAL_SAMPLES
from statistical_test_fit.ws_helper import freeze_pdf_params


DEFAULT_OUTPUT_DIR = REPO_ROOT / "datacards"
WORKSPACE_NAME = "combined_workspace"
CENTRAL_WORKSPACE_FILENAME = "workspace.root"
NON_RESONANT_WORKSPACE_PATH = REPO_ROOT / "non_resonant_background_workspace.root"
NON_RESONANT_WORKSPACE_NAME = "non_resonant_background_ws"
NON_RESONANT_SUMMARY_PATH = (
    REPO_ROOT / "plots" / "fit_2d_data" / "non_resonant_fit_summary.json"
)
RESONANT_SUMMARY_PATH = (
    REPO_ROOT / "plots" / "resonant_background" / "resonant_background_summary.json"
)
LEGACY_SYSTEMATICS_SOURCE = (
    REPO_ROOT
    / "from_mauricio"
    / "HZUpsilonPhotonRun2Statistics"
    / "HZUpsilonPhotonRun2Statistics"
    / "combine_helpers.py"
)

SIGNAL_WORKSPACE_NAME = "ws"
SIGNAL_SOURCE_PDF_NAME = "signal_model"
SIGNAL_SOURCE_DATA_NAME = "signal_data"

RESONANT_CONFIG = {
    "H": {
        "root_path": REPO_ROOT / "resonant_background_fit_HiggsDalitz.root",
        "workspace_name": "resonant_background_Higgs_ws",
        "source_pdf": "resonant_background_model_Higgs",
        "source_data": "resonant_background_data",
    },
    "Z": {
        "root_path": REPO_ROOT / "resonant_background_fit_ZGamma.root",
        "workspace_name": "resonant_background_Z_ws",
        "source_pdf": "resonant_background_model_Z",
        "source_data": "resonant_background_data",
    },
}

COPIED_SYSTEMATICS: dict[str, list[tuple[str, str, str, str, str]]] = {
    "H": [
        ("lumi", "lnN", "1.025", "-", "1.025"),
        ("HZ_xs_sc", "lnN", "0.933/1.046", "-", "0.933/1.046"),
        ("HZ_xs_pdf", "lnN", "1.032", "-", "1.032"),
        ("br_peak", "lnN", "-", "-", "1.06"),
        ("pu_r", "lnN", "1.006", "-", "1.009"),
        ("trg", "lnN", "1.055", "-", "1.061"),
        ("muon_id", "lnN", "1.0435", "-", "1.0435"),
        ("ph_id", "lnN", "1.012", "-", "1.012"),
        ("ele_veto", "lnN", "1.0104", "-", "1.0104"),
    ],
    "Z": [
        ("lumi", "lnN", "1.025", "-", "1.025"),
        ("HZ_xs_sc", "lnN", "1.033", "-", "1.05"),
        ("HZ_xs_pdf", "lnN", "1.0173", "-", "1.05"),
        ("br_peak", "lnN", "-", "-", "1.0"),
        ("pu_r", "lnN", "1.0065", "-", "1.0062"),
        ("trg", "lnN", "1.045", "-", "1.047"),
        ("muon_id", "lnN", "1.048", "-", "1.045"),
        ("ph_id", "lnN", "1.011", "-", "1.011"),
        ("ele_veto", "lnN", "1.0102", "-", "1.0102"),
    ],
}


@dataclass(frozen=True)
class ChannelSpec:
    name: str
    process: str
    state: str
    signal_workspace_path: Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Bundle persisted H/Z -> Upsilon(nS)+Photon workspaces into a single "
            "Combine-ready RooWorkspace and write one datacard per analysis channel."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory where datacards, copied workspace ROOT files, and README are written.",
    )
    parser.add_argument(
        "--workspace-name",
        default=WORKSPACE_NAME,
        help="Name of the bundled RooWorkspace stored in the ROOT file.",
    )
    parser.add_argument(
        "--workspace-file-name",
        default=CENTRAL_WORKSPACE_FILENAME,
        help="Filename of the bundled ROOT workspace file.",
    )
    parser.add_argument(
        "--strict-mode",
        action="store_true",
        help=(
            "Require strict non-resonant family selection. Relaxed mode is the default."
        ),
    )
    parser.add_argument(
        "--skip-validation",
        action="store_true",
        help="Skip text2workspace.py and combine validation commands.",
    )
    parser.add_argument(
        "--signal-mass-label",
        default="125",
        help="Mass label passed to text2workspace.py and combine during validation.",
    )
    return parser


def log(message: str) -> None:
    print(f"[bundler] {message}")


def log_kv(key: str, value: Any) -> None:
    print(f"[bundler]   - {key}: {value}")


def relative_to_repo(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def markdown_table(headers: list[str], rows: Iterable[Iterable[Any]]) -> str:
    def stringify(value: Any) -> str:
        return str(value).replace("|", r"\|")

    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(stringify(value) for value in row) + " |")
    return "\n".join(lines)


def ensure_file(path: Path, hint: str | None = None) -> None:
    if path.exists():
        return
    message = f"Required input not found: {relative_to_repo(path)}"
    if hint:
        message += f". {hint}"
    raise FileNotFoundError(message)


def boson_obs_name(channel: str) -> str:
    return f"{channel}_boson_mass"


def upsilon_obs_name(channel: str) -> str:
    return f"{channel}_upsilon_mass"


def make_channel_specs() -> list[ChannelSpec]:
    signal_paths: dict[tuple[str, str], Path] = {}
    for sample in SIGNAL_SAMPLES:
        signal_paths[(sample.process, sample.state)] = (
            REPO_ROOT / f"signal_workspace_{sample.inner_file_name}.root"
        )

    channels: list[ChannelSpec] = []
    for process in ("H", "Z"):
        for state in ("1S", "2S", "3S"):
            channel_name = (
                f"HToUpsilon{state}Photon"
                if process == "H"
                else f"ZToUpsilon{state}Photon"
            )
            channels.append(
                ChannelSpec(
                    name=channel_name,
                    process=process,
                    state=state,
                    signal_workspace_path=signal_paths[(process, state)],
                )
            )
    return channels


def load_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        return json.load(handle)


def open_workspace(path: Path, workspace_name: str) -> tuple[Any, Any]:
    root_file = TFile.Open(str(path))
    if root_file is None or root_file.IsZombie():
        raise RuntimeError(f"Could not open ROOT file {relative_to_repo(path)}")
    workspace = root_file.Get(workspace_name)
    if workspace is None:
        raise RuntimeError(
            f"Could not find workspace {workspace_name!r} in {relative_to_repo(path)}"
        )
    return root_file, workspace


def ensure_channel_observables(workspace: RooWorkspace, channel: str) -> None:
    boson_name = boson_obs_name(channel)
    upsilon_name = upsilon_obs_name(channel)
    if not workspace.var(boson_name):
        workspace.factory(f"{boson_name}[{BOSON_MASS_LOWER},{BOSON_MASS_UPPER}]")
    if not workspace.var(upsilon_name):
        workspace.factory(f"{upsilon_name}[{UPSILON_MASS_LOWER},{UPSILON_MASS_UPPER}]")

    boson_mass = workspace.var(boson_name)
    upsilon_mass = workspace.var(upsilon_name)
    if not boson_mass or not upsilon_mass:
        raise RuntimeError(f"Failed to create observables for {channel}")

    boson_mass.SetTitle("m_{#mu#mu#gamma}")
    boson_mass.setUnit("GeV")
    boson_mass.setBins(60)
    boson_mass.setRange("LEFT", LEFT_SIDEBAND_LOWER, LEFT_SIDEBAND_UPPER)
    boson_mass.setRange("MIDDLE", MIDDLE_SIDEBAND_LOWER, MIDDLE_SIDEBAND_UPPER)
    boson_mass.setRange("RIGHT", RIGHT_SIDEBAND_LOWER, RIGHT_SIDEBAND_UPPER)
    boson_mass.setRange("FULL", BOSON_MASS_LOWER, BOSON_MASS_UPPER)

    upsilon_mass.SetTitle("m_{#mu#mu}")
    upsilon_mass.setUnit("GeV")
    upsilon_mass.setBins(60)


def pdf_names(workspace: RooWorkspace) -> set[str]:
    return {pdf.GetName() for pdf in list(workspace.allPdfs())}


def import_pdf_from_workspace(
    target_workspace: RooWorkspace,
    source_workspace,
    source_pdf_name: str,
    target_pdf_name: str,
    channel: str,
    role_tag: str,
) -> str:
    source_pdf = source_workspace.pdf(source_pdf_name)
    if source_pdf is None:
        raise RuntimeError(
            f"Could not find PDF {source_pdf_name!r} in source workspace {source_workspace.GetName()}"
        )

    before_names = pdf_names(target_workspace)
    clone = source_pdf.cloneTree()
    clone.SetName(target_pdf_name)
    getattr(target_workspace, "import")(
        clone,
        ROOT.RooFit.RenameAllNodes(role_tag),
        ROOT.RooFit.RenameAllVariablesExcept(role_tag, "boson_mass,upsilon_mass"),
        ROOT.RooFit.RenameVariable("boson_mass", boson_obs_name(channel)),
        ROOT.RooFit.RenameVariable("upsilon_mass", upsilon_obs_name(channel)),
    )
    after_names = pdf_names(target_workspace)

    imported_names = sorted(
        name
        for name in (after_names - before_names)
        if name.startswith(target_pdf_name)
    )
    if not imported_names:
        if not target_workspace.pdf(target_pdf_name):
            raise RuntimeError(
                f"Unable to identify imported top-level PDF for {target_pdf_name!r}"
            )
        return target_pdf_name

    imported_name = imported_names[0]
    if imported_name != target_pdf_name:
        imported_pdf = target_workspace.pdf(imported_name)
        if not imported_pdf:
            raise RuntimeError(
                f"Imported PDF {imported_name!r} was not found for final rename"
            )
        imported_pdf.SetName(target_pdf_name)
    return target_pdf_name


def import_dataset_from_workspace(
    target_workspace: RooWorkspace,
    source_workspace,
    source_dataset_name: str,
    target_dataset_name: str,
    channel: str,
    dataset_role: str,
) -> str:
    source_data = source_workspace.data(source_dataset_name)
    if source_data is None:
        raise RuntimeError(
            f"Could not find dataset {source_dataset_name!r} in source workspace {source_workspace.GetName()}"
        )

    clone = source_data.Clone(f"{target_dataset_name}_clone")
    import_options = [
        RooFit.Rename(target_dataset_name),
        RooFit.RenameVariable("boson_mass", boson_obs_name(channel)),
        RooFit.RenameVariable("upsilon_mass", upsilon_obs_name(channel)),
    ]
    if source_data.get().find("weight"):
        import_options.append(
            RooFit.RenameVariable("weight", f"{channel}_{dataset_role}_weight")
        )
    getattr(target_workspace, "import")(clone, *import_options)
    return target_dataset_name


def signal_norm_from_workspace(signal_workspace_path: Path) -> float:
    root_file, workspace = open_workspace(signal_workspace_path, SIGNAL_WORKSPACE_NAME)
    data = workspace.data(SIGNAL_SOURCE_DATA_NAME)
    if data is None:
        raise RuntimeError(
            f"Could not find dataset {SIGNAL_SOURCE_DATA_NAME!r} in {relative_to_repo(signal_workspace_path)}"
        )
    value = float(data.sumEntries())
    root_file.Close()
    return value


def build_nonres_candidate_selection(
    nonres_summary: dict[str, Any],
    relaxed_mode: bool,
) -> dict[str, Any]:
    mode_key = "relaxed" if relaxed_mode else "strict"
    families = nonres_summary.get("families", {})
    if "johnson" not in families:
        raise RuntimeError("Non-resonant summary is missing the Johnson family.")

    selected_candidates: list[dict[str, Any]] = []
    seen_model_names: set[str] = set()
    selection_rows: list[list[Any]] = []

    johnson_candidate = dict(families["johnson"]["candidates"][0])
    johnson_candidate["selection_role"] = "nominal"
    johnson_candidate["workspace_state_index"] = len(selected_candidates)
    johnson_candidate["workspace_state_label"] = (
        f"_pdf{johnson_candidate['workspace_state_index']}"
    )
    johnson_candidate["readable_state_label"] = "johnson_nominal"
    selected_candidates.append(johnson_candidate)
    seen_model_names.add(johnson_candidate["model_name"])
    selection_rows.append(
        [
            "johnson",
            mode_key,
            "nominal",
            johnson_candidate["model_name"],
            johnson_candidate["scan_order"],
            f"{johnson_candidate['initial_norm_estimate']:.6f}",
            johnson_candidate["workspace_state_label"],
            johnson_candidate["readable_state_label"],
        ]
    )

    for family_name in ("bernstein", "chebychev"):
        family_info = families.get(family_name)
        if family_info is None:
            raise RuntimeError(
                f"Non-resonant summary is missing the {family_name!r} family."
            )
        selection = family_info.get("selections", {}).get(mode_key)
        if selection is None:
            raise RuntimeError(
                f"Non-resonant summary does not contain {mode_key!r} selection for {family_name}."
            )
        if selection.get("status") != "ok":
            raise RuntimeError(
                f"{mode_key.title()} selection failed for {family_name}: {selection.get('message', 'unknown error')}"
            )

        winner_index = int(selection["winner_index"])
        candidates = family_info["candidates"]
        for candidate_index, selection_role in (
            (winner_index - 1, "below"),
            (winner_index, "winner"),
            (winner_index + 1, "above"),
        ):
            if candidate_index < 0 or candidate_index >= len(candidates):
                continue
            candidate = dict(candidates[candidate_index])
            candidate["selection_role"] = selection_role
            model_name = candidate["model_name"]
            if model_name in seen_model_names:
                continue
            candidate["workspace_state_index"] = len(selected_candidates)
            candidate["workspace_state_label"] = (
                f"_pdf{candidate['workspace_state_index']}"
            )
            candidate["readable_state_label"] = (
                f"{family_name}_{selection_role}_order{candidate['scan_order']}"
            )
            seen_model_names.add(model_name)
            selected_candidates.append(candidate)
            selection_rows.append(
                [
                    family_name,
                    mode_key,
                    selection_role,
                    model_name,
                    candidate["scan_order"],
                    f"{candidate['initial_norm_estimate']:.6f}",
                    candidate["workspace_state_label"],
                    candidate["readable_state_label"],
                ]
            )

    return {
        "mode_key": mode_key,
        "selected_candidates": selected_candidates,
        "selection_rows": selection_rows,
        "johnson_initial_norm": float(johnson_candidate["initial_norm_estimate"]),
        "sideband_entries": int(nonres_summary["entries"]["sidebands"]),
        "full_entries": int(nonres_summary["entries"]["full"]),
    }


def nonres_alias_name(channel: str, candidate: dict[str, Any]) -> str:
    family = candidate["pdf_family"]
    if family == "johnson":
        return f"{channel}_non_resonant_background_johnson_nominal_pdf"
    order = candidate["scan_order"]
    return f"{channel}_non_resonant_background_{family}_order{order}_pdf"


def make_signal_dataset_name(channel: str) -> str:
    return f"{channel}_signal_data"


def make_resonant_dataset_name(channel: str) -> str:
    return f"{channel}_resonant_background_data"


def make_nonres_sideband_dataset_name(channel: str) -> str:
    return f"{channel}_non_resonant_sidebands"


def import_signal_artifacts(
    target_workspace: RooWorkspace,
    channel: ChannelSpec,
) -> dict[str, Any]:
    log(f"Importing signal artifacts for {channel.name}")
    root_file, source_workspace = open_workspace(
        channel.signal_workspace_path,
        SIGNAL_WORKSPACE_NAME,
    )
    ensure_channel_observables(target_workspace, channel.name)
    import_pdf_from_workspace(
        target_workspace,
        source_workspace,
        SIGNAL_SOURCE_PDF_NAME,
        f"{channel.name}_signal_pdf",
        channel.name,
        f"__{channel.name}__signal",
    )
    import_dataset_from_workspace(
        target_workspace,
        source_workspace,
        SIGNAL_SOURCE_DATA_NAME,
        make_signal_dataset_name(channel.name),
        channel.name,
        "signal",
    )

    signal_norm = RooRealVar(
        f"{channel.name}_signal_pdf_norm",
        f"{channel.name}_signal_pdf_norm",
        signal_norm_from_workspace(channel.signal_workspace_path),
    )
    signal_norm.setConstant(True)
    getattr(target_workspace, "import")(signal_norm)
    root_file.Close()

    return {
        "source_workspace": relative_to_repo(channel.signal_workspace_path),
        "pdf": f"{channel.name}_signal_pdf",
        "dataset": make_signal_dataset_name(channel.name),
        "norm": float(signal_norm.getVal()),
    }


def import_resonant_artifacts(
    target_workspace: RooWorkspace,
    channel: ChannelSpec,
    resonant_summary: dict[str, Any],
) -> dict[str, Any]:
    log(f"Importing resonant-background artifacts for {channel.name}")
    source_info = RESONANT_CONFIG[channel.process]
    root_file, source_workspace = open_workspace(
        source_info["root_path"],
        source_info["workspace_name"],
    )
    ensure_channel_observables(target_workspace, channel.name)
    import_pdf_from_workspace(
        target_workspace,
        source_workspace,
        source_info["source_pdf"],
        f"{channel.name}_resonant_background_pdf",
        channel.name,
        f"__{channel.name}__resonant",
    )
    import_dataset_from_workspace(
        target_workspace,
        source_workspace,
        source_info["source_data"],
        make_resonant_dataset_name(channel.name),
        channel.name,
        "resonant",
    )

    initial_norm = float(resonant_summary["process_initial_norms"][channel.process])
    resonant_norm = RooRealVar(
        f"{channel.name}_resonant_background_pdf_norm",
        f"{channel.name}_resonant_background_pdf_norm",
        initial_norm,
        0.0,
        max(initial_norm * 10.0, 10.0),
    )
    resonant_norm.setConstant(False)
    getattr(target_workspace, "import")(resonant_norm)
    root_file.Close()

    return {
        "source_workspace": relative_to_repo(source_info["root_path"]),
        "pdf": f"{channel.name}_resonant_background_pdf",
        "dataset": make_resonant_dataset_name(channel.name),
        "norm": float(resonant_norm.getVal()),
    }


def import_non_resonant_artifacts(
    target_workspace: RooWorkspace,
    channel: ChannelSpec,
    nonres_workspace,
    nonres_selection: dict[str, Any],
) -> dict[str, Any]:
    log(f"Importing non-resonant artifacts for {channel.name}")
    import_dataset_from_workspace(
        target_workspace,
        nonres_workspace,
        "data_obs",
        f"{channel.name}_data_obs",
        channel.name,
        "data_obs",
    )
    import_dataset_from_workspace(
        target_workspace,
        nonres_workspace,
        "data_sidebands",
        make_nonres_sideband_dataset_name(channel.name),
        channel.name,
        "non_resonant_sidebands",
    )
    ensure_channel_observables(target_workspace, channel.name)

    category = RooCategory(f"{channel.name}_pdfindex", f"{channel.name}_pdfindex")
    getattr(target_workspace, "import")(category)

    multipdf_inputs = RooArgList()
    candidate_aliases: list[str] = []
    state_mappings: list[dict[str, Any]] = []
    for candidate in nonres_selection["selected_candidates"]:
        alias_name = nonres_alias_name(channel.name, candidate)
        import_pdf_from_workspace(
            target_workspace,
            nonres_workspace,
            candidate["model_name"],
            alias_name,
            channel.name,
            f"__{channel.name}__{candidate['model_name']}",
        )
        freeze_pdf_params(
            target_workspace.pdf(alias_name),
            RooArgSet(
                target_workspace.var(upsilon_obs_name(channel.name)),
                target_workspace.var(boson_obs_name(channel.name)),
            ),
        )
        multipdf_inputs.add(target_workspace.pdf(alias_name))
        candidate_aliases.append(alias_name)
        state_mappings.append(
            {
                "index": int(candidate["workspace_state_index"]),
                "workspace_label": candidate["workspace_state_label"],
                "readable_label": candidate["readable_state_label"],
                "pdf_name": alias_name,
                "source_model_name": candidate["model_name"],
                "selection_role": candidate["selection_role"],
                "scan_order": candidate["scan_order"],
                "pdf_family": candidate["pdf_family"],
            }
        )

    if not candidate_aliases:
        raise RuntimeError(
            f"No non-resonant candidates were selected for {channel.name}"
        )
    multipdf = RooMultiPdf(
        f"{channel.name}_non_resonant_background_pdf",
        f"{channel.name}_non_resonant_background_pdf",
        target_workspace.cat(f"{channel.name}_pdfindex"),
        multipdf_inputs,
    )
    getattr(target_workspace, "import")(multipdf)
    target_workspace.cat(f"{channel.name}_pdfindex").setIndex(0)

    nonres_norm = RooRealVar(
        f"{channel.name}_non_resonant_background_pdf_norm",
        f"{channel.name}_non_resonant_background_pdf_norm",
        float(nonres_selection["johnson_initial_norm"]),
        0.0,
        max(float(nonres_selection["full_entries"]) * 10.0, 10.0),
    )
    nonres_norm.setConstant(False)
    getattr(target_workspace, "import")(nonres_norm)

    return {
        "source_workspace": relative_to_repo(NON_RESONANT_WORKSPACE_PATH),
        "pdf": f"{channel.name}_non_resonant_background_pdf",
        "dataset": f"{channel.name}_data_obs",
        "sideband_dataset": make_nonres_sideband_dataset_name(channel.name),
        "norm": float(nonres_norm.getVal()),
        "mode": nonres_selection["mode_key"],
        "candidates": candidate_aliases,
        "state_mappings": state_mappings,
    }


def write_datacard(channel_dir: Path, channel: ChannelSpec) -> Path:
    systematics = COPIED_SYSTEMATICS[channel.process]
    lines = [
        f"# Datacard for {channel.name}",
        f"# Process: {channel.process}, Upsilon state: {channel.state}",
        "imax 1 number of channels",
        "jmax 2 number of backgrounds",
        "kmax * number of nuisance parameters",
        "------------",
        f"shapes signal {channel.name} {CENTRAL_WORKSPACE_FILENAME} {WORKSPACE_NAME}:{channel.name}_signal_pdf",
        f"shapes non_resonant_bkg {channel.name} {CENTRAL_WORKSPACE_FILENAME} {WORKSPACE_NAME}:{channel.name}_non_resonant_background_pdf",
        f"shapes resonant_bkg {channel.name} {CENTRAL_WORKSPACE_FILENAME} {WORKSPACE_NAME}:{channel.name}_resonant_background_pdf",
        f"shapes data_obs {channel.name} {CENTRAL_WORKSPACE_FILENAME} {WORKSPACE_NAME}:{channel.name}_data_obs",
        "------------",
        f"bin {channel.name}",
        "observation -1",
        "------------",
        f"bin {channel.name} {channel.name} {channel.name}",
        "process signal non_resonant_bkg resonant_bkg",
        "process 0 1 2",
        "rate 1 1 1",
        "------------",
    ]
    for (
        nuisance_name,
        nuisance_pdf,
        signal_value,
        nonres_value,
        resonant_value,
    ) in systematics:
        lines.append(
            f"{nuisance_name} {nuisance_pdf} {signal_value} {nonres_value} {resonant_value}"
        )
    lines.extend(
        [
            "------------",
            f"{channel.name}_pdfindex discrete",
            "",
        ]
    )
    datacard_path = channel_dir / f"{channel.name}.txt"
    datacard_path.write_text("\n".join(lines), encoding="ascii")
    return datacard_path


def validate_channel(
    channel_dir: Path, channel_name: str, mass_label: str
) -> dict[str, Any]:
    datacard_path = channel_dir / f"{channel_name}.txt"
    text2workspace_cmd = [
        "text2workspace.py",
        datacard_path.name,
        "-m",
        mass_label,
        "-o",
        "validation_workspace.root",
    ]
    text2workspace_res = subprocess.run(
        text2workspace_cmd,
        cwd=channel_dir,
        capture_output=True,
        text=True,
    )

    combine_status = "blocked"
    combine_cmd: list[str] | None = None
    combine_output = "combine step skipped because text2workspace.py failed."
    combine_returncode = None
    if text2workspace_res.returncode == 0:
        combine_cmd = [
            "combine",
            "-M",
            "AsymptoticLimits",
            "validation_workspace.root",
            "--run",
            "blind",
            "-m",
            mass_label,
            "-n",
            f".{channel_name}_validate",
        ]
        combine_res = subprocess.run(
            combine_cmd,
            cwd=channel_dir,
            capture_output=True,
            text=True,
        )
        combine_status = "ok" if combine_res.returncode == 0 else "failed"
        combine_output = (combine_res.stdout or "") + (combine_res.stderr or "")
        combine_returncode = combine_res.returncode

    log_lines = [
        f"# Validation log for {channel_name}",
        "",
        "## text2workspace.py",
        f"Command: {' '.join(text2workspace_cmd)}",
        f"Return code: {text2workspace_res.returncode}",
        "",
        "```text",
        (text2workspace_res.stdout or "") + (text2workspace_res.stderr or ""),
        "```",
        "",
        "## combine -M AsymptoticLimits",
        f"Command: {' '.join(combine_cmd) if combine_cmd else '(skipped)'}",
        f"Status: {combine_status}",
        "",
        "```text",
        combine_output,
        "```",
        "",
    ]
    log_path = channel_dir / "validation.log"
    log_path.write_text("\n".join(log_lines), encoding="ascii")
    return {
        "text2workspace_status": "ok"
        if text2workspace_res.returncode == 0
        else "failed",
        "combine_status": combine_status,
        "text2workspace_returncode": text2workspace_res.returncode,
        "combine_returncode": combine_returncode,
        "log": relative_to_repo(log_path),
    }


def collect_workspace_objects(workspace: RooWorkspace) -> dict[str, list[str]]:
    return {
        "datasets": sorted(obj.GetName() for obj in list(workspace.allData())),
        "pdfs": sorted(obj.GetName() for obj in list(workspace.allPdfs())),
        "functions": sorted(obj.GetName() for obj in list(workspace.allFunctions())),
        "variables": sorted(obj.GetName() for obj in list(workspace.allVars())),
        "categories": sorted(obj.GetName() for obj in list(workspace.allCats())),
    }


def collect_parameter_rows(workspace: RooWorkspace) -> list[list[Any]]:
    rows: list[list[Any]] = []
    for variable in sorted(list(workspace.allVars()), key=lambda item: item.GetName()):
        name = variable.GetName()
        if name.endswith("_boson_mass") or name.endswith("_upsilon_mass"):
            status = "observable"
        else:
            status = "fixed" if variable.isConstant() else "floating"
        rows.append(
            [
                name,
                type(variable).__name__,
                status,
                f"{variable.getVal():.8g}",
                f"{variable.getMin():.8g}",
                f"{variable.getMax():.8g}",
            ]
        )
    for category in sorted(list(workspace.allCats()), key=lambda item: item.GetName()):
        rows.append(
            [
                category.GetName(),
                type(category).__name__,
                "fixed" if category.isConstant() else "floating",
                category.getLabel(),
                "-",
                "-",
            ]
        )
    return rows


def collect_channel_parameter_rows(
    parameter_rows: list[list[Any]],
    channels: list[ChannelSpec],
) -> dict[str, list[list[Any]]]:
    return {
        channel.name: [row for row in parameter_rows if channel.name in str(row[0])]
        for channel in channels
    }


def write_summary_readme(
    output_dir: Path,
    workspace_name: str,
    workspace_file_name: str,
    workspace: RooWorkspace,
    channels: list[ChannelSpec],
    channel_reports: dict[str, dict[str, Any]],
    nonres_selection: dict[str, Any],
    resonant_summary: dict[str, Any],
    validation_results: dict[str, dict[str, Any]],
    produced_files: list[Path],
    placeholders: list[str],
    missing_inputs: list[str],
) -> None:
    object_lists = collect_workspace_objects(workspace)
    parameter_rows = collect_parameter_rows(workspace)
    channel_parameter_rows = collect_channel_parameter_rows(parameter_rows, channels)

    channel_rows = []
    state_mapping_rows = []
    for channel in channels:
        report = channel_reports[channel.name]
        state_summary = ", ".join(
            f"{state['workspace_label']}={state['readable_label']}"
            for state in report["non_resonant"]["state_mappings"]
        )
        channel_rows.append(
            [
                channel.name,
                channel.process,
                channel.state,
                f"{report['signal']['norm']:.8g}",
                f"{report['resonant']['norm']:.8g}",
                f"{report['non_resonant']['norm']:.8g}",
                report["non_resonant"]["mode"],
                state_summary,
                ", ".join(
                    Path(name).stem for name in report["non_resonant"]["candidates"]
                ),
            ]
        )
        for state in report["non_resonant"]["state_mappings"]:
            state_mapping_rows.append(
                [
                    channel.name,
                    state["index"],
                    state["workspace_label"],
                    state["readable_label"],
                    state["pdf_name"],
                ]
            )

    validation_rows = []
    for channel in channels:
        if channel.name in validation_results:
            result = validation_results[channel.name]
            validation_rows.append(
                [
                    channel.name,
                    result["text2workspace_status"],
                    result["combine_status"],
                    result["log"],
                ]
            )
        else:
            validation_rows.append([channel.name, "skipped", "skipped", "-"])

    systematics_rows = []
    for process in ("H", "Z"):
        for (
            nuisance_name,
            nuisance_pdf,
            signal_value,
            nonres_value,
            resonant_value,
        ) in COPIED_SYSTEMATICS[process]:
            systematics_rows.append(
                [
                    process,
                    nuisance_name,
                    nuisance_pdf,
                    signal_value,
                    nonres_value,
                    resonant_value,
                    relative_to_repo(LEGACY_SYSTEMATICS_SOURCE),
                ]
            )

    lines = [
        "# Bundled Datacards",
        "",
        "## Summary",
        f"- Bundled ROOT file: `{workspace_file_name}`",
        f"- RooWorkspace name: `{workspace_name}`",
        "- Observable policy: full-range in every channel",
        f"- `upsilon_mass`: `{UPSILON_MASS_LOWER}` to `{UPSILON_MASS_UPPER}` GeV",
        f"- `boson_mass`: `{BOSON_MASS_LOWER}` to `{BOSON_MASS_UPPER}` GeV",
        f"- Non-resonant selection mode used by the bundler: `{nonres_selection['mode_key']}`",
        "- `pseudodata.py` is intentionally excluded from production datacard inputs and treated as a toy-study workflow only.",
        "",
        "## Naming Convention Chosen",
        "- Channels use the required names: `HToUpsilon1SPhoton`, `HToUpsilon2SPhoton`, `HToUpsilon3SPhoton`, `ZToUpsilon1SPhoton`, `ZToUpsilon2SPhoton`, `ZToUpsilon3SPhoton`.",
        "- Public Combine-facing PDFs are named `CHANNEL_signal_pdf`, `CHANNEL_resonant_background_pdf`, and `CHANNEL_non_resonant_background_pdf`.",
        "- Public normalizations are named `CHANNEL_signal_pdf_norm`, `CHANNEL_resonant_background_pdf_norm`, and `CHANNEL_non_resonant_background_pdf_norm`.",
        "- Public observed datasets are named `CHANNEL_data_obs`.",
        "- Public discrete categories are named `CHANNEL_pdfindex`.",
        "- `RooMultiPdf` auto-generates category state labels as `_pdfN`; the readable state names are recorded below.",
        "",
        "## Collision-Avoidance Strategy",
        "- The bundler imports cloned source PDFs and datasets from persisted workspaces, not refitted or reconstructed models.",
        "- Imported dependency graphs are renamed deterministically with channel/role tags to avoid collisions in one shared workspace.",
        "- If a required source object or unique target name is missing, the bundler fails loudly.",
        "",
        "## Workspace Layout",
        "- One shared workspace contains all six channels.",
        "- Each channel has its own observable pair, observed dataset, signal dataset, resonant dataset, signal PDF, resonant PDF, non-resonant candidate PDFs, and final non-resonant `RooMultiPdf`.",
        "- The final non-resonant `RooMultiPdf` contains Johnson nominal plus the selected winner and immediate neighbor orders for the Bernstein and Chebychev families.",
        "",
        "## Channel Inputs",
        markdown_table(
            [
                "Channel",
                "Process",
                "State",
                "Signal norm",
                "Resonant init",
                "Non-res init",
                "Non-res mode",
                "RooMultiPdf states",
                "Non-res candidates",
            ],
            channel_rows,
        ),
        "",
        "## Non-Resonant Candidate Selection",
        markdown_table(
            [
                "Family",
                "Mode",
                "Role",
                "Model",
                "Order",
                "Initial norm estimate",
                "Workspace state",
                "Readable state",
            ],
            nonres_selection["selection_rows"],
        ),
        "",
        "## RooMultiPdf State Mapping",
        markdown_table(
            ["Channel", "Index", "Workspace label", "Readable label", "Imported PDF"],
            state_mapping_rows,
        ),
        "",
        "## Imported Objects",
        "### Datasets",
        "```text",
        *object_lists["datasets"],
        "```",
        "",
        "### PDFs",
        "```text",
        *object_lists["pdfs"],
        "```",
        "",
        "### Functions",
        "```text",
        *object_lists["functions"],
        "```",
        "",
        "### Variables",
        "```text",
        *object_lists["variables"],
        "```",
        "",
        "### Categories",
        "```text",
        *object_lists["categories"],
        "```",
        "",
        "## Parameters",
        "",
        "## Copied Systematics And Sources",
        markdown_table(
            [
                "Process",
                "Nuisance",
                "Type",
                "Signal",
                "Non-resonant",
                "Resonant",
                "Source",
            ],
            systematics_rows,
        ),
        "",
        "## What Was Copied From Legacy Inputs",
        f"- Systematic names and values were copied from `{relative_to_repo(LEGACY_SYSTEMATICS_SOURCE)}`.",
        "- Signal normalizations come from the weighted `signal_data.sumEntries()` values stored in the signal workspaces produced by `signal.py`.",
        f"- Resonant initial normalization seeds come from the produced summary `{relative_to_repo(RESONANT_SUMMARY_PATH)}`: `H` uses the stored Higgs resonant dataset sum of weights and `Z` uses the stored control-region extrapolation `y0`.",
        f"- Non-resonant candidate selection and initial normalization estimates come from `{relative_to_repo(NON_RESONANT_SUMMARY_PATH)}` and the persisted non-resonant workspace `{relative_to_repo(NON_RESONANT_WORKSPACE_PATH)}`.",
        "",
        "## Fixed Vs Floating Identification",
        "- `RooRealVar.isConstant()` is used for all real-valued parameters stored in the bundled workspace.",
        "- `RooCategory.isConstant()` is used for discrete parameters such as `CHANNEL_pdfindex`.",
        "- Observables are reported as `observable` rather than `floating`.",
        "",
        "## Produced Files",
        "```text",
        *sorted(relative_to_repo(path) for path in produced_files),
        "```",
        "",
        "## Detected Placeholders",
    ]

    parameter_headers = ["Name", "Type", "Status", "Value/Label", "Min", "Max"]
    parameter_insert_index = lines.index("## Copied Systematics And Sources")
    parameter_sections: list[str] = []
    for channel in channels:
        parameter_sections.extend(
            [
                f"### {channel.name}",
                markdown_table(
                    parameter_headers,
                    channel_parameter_rows[channel.name],
                ),
                "",
            ]
        )
    lines[parameter_insert_index:parameter_insert_index] = parameter_sections

    if placeholders:
        lines.extend([f"- {item}" for item in placeholders])
    else:
        lines.append("- None.")

    lines.extend(["", "## Missing Inputs"])
    if missing_inputs:
        lines.extend([f"- {item}" for item in missing_inputs])
    else:
        lines.append("- None.")

    lines.extend(
        [
            "",
            "## Validation",
            markdown_table(
                ["Channel", "text2workspace.py", "combine -M AsymptoticLimits", "Log"],
                validation_rows,
            ),
            "",
            "## What Remains Unresolved",
            "- The bundler intentionally excludes `pseudodata.py` from production datacard inputs.",
            "- Systematic correlations still follow the legacy Mauricio helper conventions because this repo does not yet expose an independent native nuisance table.",
        ]
    )

    if validation_results and any(
        result["text2workspace_status"] != "ok"
        for result in validation_results.values()
    ):
        lines.append(
            "- Validation failures, if any, are channel-local and captured in the `validation.log` files listed above."
        )

    lines.append("")
    (output_dir / "README.md").write_text("\n".join(lines), encoding="ascii")


def verify_required_outputs(channels: list[ChannelSpec]) -> None:
    log("Checking persisted upstream outputs")
    for channel in channels:
        ensure_file(
            channel.signal_workspace_path,
            "Run `python3 scripts/signal.py` first.",
        )
        log_kv("signal workspace", relative_to_repo(channel.signal_workspace_path))
    ensure_file(
        RESONANT_CONFIG["H"]["root_path"],
        "Run `python3 scripts/resonant_background.py` first.",
    )
    ensure_file(
        RESONANT_CONFIG["Z"]["root_path"],
        "Run `python3 scripts/resonant_background.py` first.",
    )
    ensure_file(
        RESONANT_SUMMARY_PATH,
        "Run `python3 scripts/resonant_background.py` with the updated code so the summary is written.",
    )
    ensure_file(
        NON_RESONANT_WORKSPACE_PATH,
        "Run `python3 scripts/non_resonant_background.py` with the updated code so the non-resonant workspace is written.",
    )
    ensure_file(
        NON_RESONANT_SUMMARY_PATH,
        "Run `python3 scripts/non_resonant_background.py` with the updated code so the non-resonant summary is written.",
    )
    ensure_file(LEGACY_SYSTEMATICS_SOURCE)
    log_kv(
        "resonant workspace (H)", relative_to_repo(RESONANT_CONFIG["H"]["root_path"])
    )
    log_kv(
        "resonant workspace (Z)", relative_to_repo(RESONANT_CONFIG["Z"]["root_path"])
    )
    log_kv("resonant summary", relative_to_repo(RESONANT_SUMMARY_PATH))
    log_kv("non-resonant workspace", relative_to_repo(NON_RESONANT_WORKSPACE_PATH))
    log_kv("non-resonant summary", relative_to_repo(NON_RESONANT_SUMMARY_PATH))


def validate_summary_shapes(
    nonres_summary: dict[str, Any],
    resonant_summary: dict[str, Any],
) -> None:
    for family in ("johnson", "bernstein", "chebychev"):
        family_info = nonres_summary.get("families", {}).get(family)
        if family_info is None:
            raise RuntimeError(
                f"Non-resonant summary is missing family {family!r}. Re-run `python3 non_resonant_background.py`."
            )
        if "selections" not in family_info:
            raise RuntimeError(
                "Non-resonant summary does not contain strict/relaxed selections. "
                "Re-run `python3 scripts/non_resonant_background.py` with the updated code."
            )
        if family == "johnson" and not family_info.get("candidates"):
            raise RuntimeError(
                "Johnson family has no saved candidates in the non-resonant summary."
            )
        for candidate in family_info.get("candidates", []):
            if "initial_norm_estimate" not in candidate:
                raise RuntimeError(
                    "Non-resonant summary is missing candidate normalization estimates. "
                    "Re-run `python3 scripts/non_resonant_background.py` with the updated code."
                )

    if "process_initial_norms" not in resonant_summary:
        raise RuntimeError(
            "Resonant summary is missing `process_initial_norms`. "
            "Re-run `python3 scripts/resonant_background.py` with the updated code."
        )


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    channels = make_channel_specs()
    produced_files: list[Path] = []
    placeholders: list[str] = []
    missing_inputs: list[str] = []

    log("Starting bundled workspace build")
    log_kv("output directory", relative_to_repo(output_dir))
    log_kv("workspace name", args.workspace_name)
    log_kv("workspace file name", args.workspace_file_name)
    log_kv("selection mode", "strict" if args.strict_mode else "relaxed")
    log_kv("validation", "disabled" if args.skip_validation else "enabled")

    verify_required_outputs(channels)
    nonres_summary = load_json(NON_RESONANT_SUMMARY_PATH)
    resonant_summary = load_json(RESONANT_SUMMARY_PATH)
    validate_summary_shapes(nonres_summary, resonant_summary)
    nonres_selection = build_nonres_candidate_selection(
        nonres_summary, relaxed_mode=not args.strict_mode
    )

    log("Selected non-resonant candidates for the final RooMultiPdf")
    for row in nonres_selection["selection_rows"]:
        log_kv(
            f"{row[0]} {row[2]}",
            f"model={row[3]}, order={row[4]}, init_norm={row[5]}",
        )

    bundled_workspace = RooWorkspace(args.workspace_name)
    nonres_root_file, nonres_workspace = open_workspace(
        NON_RESONANT_WORKSPACE_PATH,
        NON_RESONANT_WORKSPACE_NAME,
    )

    channel_reports: dict[str, dict[str, Any]] = {}
    for channel in channels:
        log(f"Building channel {channel.name}")
        ensure_channel_observables(bundled_workspace, channel.name)
        nonres_report = import_non_resonant_artifacts(
            bundled_workspace,
            channel,
            nonres_workspace,
            nonres_selection,
        )
        signal_report = import_signal_artifacts(bundled_workspace, channel)
        resonant_report = import_resonant_artifacts(
            bundled_workspace,
            channel,
            resonant_summary,
        )
        channel_reports[channel.name] = {
            "signal": signal_report,
            "resonant": resonant_report,
            "non_resonant": nonres_report,
        }
        log_kv("signal pdf", signal_report["pdf"])
        log_kv("resonant pdf", resonant_report["pdf"])
        log_kv("non-resonant pdf", nonres_report["pdf"])
        log_kv("data_obs", nonres_report["dataset"])

    nonres_root_file.Close()

    central_workspace_path = output_dir / args.workspace_file_name
    bundled_workspace.writeToFile(str(central_workspace_path))
    produced_files.append(central_workspace_path)
    log(f"Wrote bundled workspace: {relative_to_repo(central_workspace_path)}")

    for channel in channels:
        channel_dir = output_dir / channel.name
        channel_dir.mkdir(parents=True, exist_ok=True)
        copied_workspace_path = channel_dir / "workspace.root"
        shutil.copy2(central_workspace_path, copied_workspace_path)
        datacard_path = write_datacard(channel_dir, channel)
        produced_files.extend([copied_workspace_path, datacard_path])
        log(f"Wrote channel outputs for {channel.name}")
        log_kv("datacard", relative_to_repo(datacard_path))
        log_kv("workspace copy", relative_to_repo(copied_workspace_path))

    validation_results: dict[str, dict[str, Any]] = {}
    if not args.skip_validation:
        log("Running validation commands")
        for channel in channels:
            channel_dir = output_dir / channel.name
            log(f"Validating {channel.name}")
            validation_results[channel.name] = validate_channel(
                channel_dir,
                channel.name,
                args.signal_mass_label,
            )
            produced_files.append(channel_dir / "validation.log")
            log_kv(
                "text2workspace.py",
                validation_results[channel.name]["text2workspace_status"],
            )
            log_kv(
                "combine -M AsymptoticLimits",
                validation_results[channel.name]["combine_status"],
            )

    readme_path = output_dir / "README.md"
    produced_files.append(readme_path)
    write_summary_readme(
        output_dir=output_dir,
        workspace_name=args.workspace_name,
        workspace_file_name=args.workspace_file_name,
        workspace=bundled_workspace,
        channels=channels,
        channel_reports=channel_reports,
        nonres_selection=nonres_selection,
        resonant_summary=resonant_summary,
        validation_results=validation_results,
        produced_files=produced_files,
        placeholders=placeholders,
        missing_inputs=missing_inputs,
    )
    log(f"Wrote summary report: {relative_to_repo(readme_path)}")

    log("Produced artifacts")
    for path in sorted(produced_files):
        log_kv("file", relative_to_repo(path))

    if validation_results:
        log("Validation summary")
        for channel_name, result in validation_results.items():
            log_kv(
                channel_name,
                f"text2workspace={result['text2workspace_status']}, combine={result['combine_status']}",
            )

    log("Bundler setup complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
