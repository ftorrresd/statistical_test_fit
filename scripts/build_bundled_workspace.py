#!/usr/bin/env python3

from __future__ import annotations

import argparse
import html
import json
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
DEFAULT_WORKSPACE_NAME = "combined_workspace"
DEFAULT_WORKSPACE_FILENAME = "workspace.root"
DEFAULT_DATACARD_FILENAME = "datacard.txt"
ANALYSIS_BIN_NAME = "cat1"
NON_RESONANT_PROCESS_NAME = "non_resonant_bkg"
RESONANT_H_PROCESS_NAME = "resonant_H_bkg"
RESONANT_Z_PROCESS_NAME = "resonant_Z_bkg"
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
        "public_process": RESONANT_H_PROCESS_NAME,
    },
    "Z": {
        "root_path": REPO_ROOT / "resonant_background_fit_ZGamma.root",
        "workspace_name": "resonant_background_Z_ws",
        "source_pdf": "resonant_background_model_Z",
        "public_process": RESONANT_Z_PROCESS_NAME,
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
class SignalProcessSpec:
    process: str
    state: str
    signal_workspace_path: Path

    @property
    def combine_name(self) -> str:
        return f"{self.process}_{self.state}"

    @property
    def channel_label(self) -> str:
        return (
            f"HToUpsilon{self.state}Photon"
            if self.process == "H"
            else f"ZToUpsilon{self.state}Photon"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Bundle the persisted H/Z -> Upsilon(nS)+Photon workspaces into one "
            "Combine-ready RooWorkspace and write a single simultaneous parametric datacard."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory where the bundled workspace, datacard, and bundle_summary.html are written.",
    )
    parser.add_argument(
        "--workspace-name",
        default=DEFAULT_WORKSPACE_NAME,
        help="Name of the bundled RooWorkspace stored in the ROOT file.",
    )
    parser.add_argument(
        "--workspace-file-name",
        default=DEFAULT_WORKSPACE_FILENAME,
        help="Filename of the bundled ROOT workspace file.",
    )
    parser.add_argument(
        "--datacard-file-name",
        default=DEFAULT_DATACARD_FILENAME,
        help="Filename of the simultaneous datacard written next to the workspace.",
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


def html_escape(value: Any) -> str:
    escaped = html.escape(str(value), quote=True)
    return escaped.encode("ascii", "xmlcharrefreplace").decode("ascii")


def html_table(headers: list[str], rows: Iterable[Iterable[Any]]) -> str:
    body_rows = [list(row) for row in rows]
    header_html = "".join(f"<th>{html_escape(header)}</th>" for header in headers)
    if body_rows:
        body_html = "\n".join(
            "<tr>"
            + "".join(f"<td>{html_escape(value)}</td>" for value in row)
            + "</tr>"
            for row in body_rows
        )
    else:
        body_html = f'<tr><td class="empty" colspan="{len(headers)}">No rows</td></tr>'
    return f"""<div class="table-wrap">
<table>
<thead><tr>{header_html}</tr></thead>
<tbody>
{body_html}
</tbody>
</table>
</div>"""


def html_list(items: Iterable[Any]) -> str:
    return "<ul>" + "".join(f"<li>{html_escape(item)}</li>" for item in items) + "</ul>"


def html_code_block(lines: Iterable[Any] | str) -> str:
    if isinstance(lines, str):
        text = lines
    else:
        text = "\n".join(str(line) for line in lines)
    return f"<pre><code>{html_escape(text)}</code></pre>"


def html_metric_grid(items: Iterable[tuple[str, Any]]) -> str:
    cards = "\n".join(
        f"""<div class="metric-card">
<div class="metric-label">{html_escape(label)}</div>
<div class="metric-value">{html_escape(value)}</div>
</div>"""
        for label, value in items
    )
    return f'<div class="metric-grid">\n{cards}\n</div>'


def html_section(title: str, body: str, intro: str | None = None) -> str:
    intro_html = f'<p class="section-intro">{html_escape(intro)}</p>' if intro else ""
    return f"""<section class="section-card">
<div class="section-heading">
<h2>{html_escape(title)}</h2>
{intro_html}
</div>
{body}
</section>"""


def html_page(title: str, subtitle: str, body: str) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html_escape(title)}</title>
<style>
:root {{
  color-scheme: light;
  --bg: #f5f7fb;
  --bg-accent: #e7efff;
  --card: #ffffff;
  --text: #172033;
  --muted: #667085;
  --line: #dde5f2;
  --accent: #3b82f6;
  --shadow: 0 20px 55px rgba(15, 23, 42, 0.11);
  --mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
  --sans: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  background:
    radial-gradient(circle at top left, var(--bg-accent), transparent 34rem),
    linear-gradient(180deg, #ffffff 0%, var(--bg) 28rem);
  color: var(--text);
  font-family: var(--sans);
}}
.page {{ width: min(1480px, calc(100% - 48px)); margin: 0 auto; padding: 36px 0 64px; }}
.hero {{
  border-radius: 28px;
  padding: 34px;
  color: #ffffff;
  background:
    linear-gradient(135deg, rgba(29, 78, 216, 0.96), rgba(14, 116, 144, 0.9)),
    radial-gradient(circle at top right, rgba(255, 255, 255, 0.32), transparent 26rem);
  box-shadow: var(--shadow);
}}
.eyebrow {{ margin: 0 0 10px; font-size: 0.77rem; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; opacity: 0.82; }}
h1 {{ margin: 0; font-size: clamp(2rem, 4vw, 3.8rem); line-height: 1.02; letter-spacing: -0.045em; }}
.subtitle {{ max-width: 860px; margin: 16px 0 0; color: rgba(255, 255, 255, 0.86); font-size: 1.04rem; line-height: 1.65; }}
.metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 14px; margin-top: 18px; }}
.metric-grid + ul {{ margin-top: 18px; }}
.metric-card {{ border: 1px solid var(--line); border-radius: 18px; background: #ffffff; padding: 16px 18px; box-shadow: 0 10px 24px rgba(15, 23, 42, 0.06); }}
.metric-label {{ color: var(--muted); font-size: 0.74rem; font-weight: 800; letter-spacing: 0.09em; text-transform: uppercase; }}
.metric-value {{ margin-top: 8px; color: var(--text); font-size: 1.08rem; font-weight: 750; overflow-wrap: anywhere; }}
.section-card {{ margin-top: 22px; padding: 24px; background: rgba(255, 255, 255, 0.92); border: 1px solid var(--line); border-radius: 24px; box-shadow: 0 12px 32px rgba(15, 23, 42, 0.07); }}
.section-heading {{ display: flex; align-items: baseline; justify-content: space-between; gap: 18px; flex-wrap: wrap; margin-bottom: 16px; }}
h2 {{ margin: 0; font-size: 1.22rem; letter-spacing: -0.02em; }}
h3 {{ margin: 22px 0 10px; font-size: 1rem; color: #344054; }}
.section-intro {{ margin: 0; color: var(--muted); line-height: 1.6; max-width: 820px; }}
ul {{ margin: 0; padding-left: 1.2rem; color: #344054; line-height: 1.7; }}
.table-wrap {{ width: 100%; overflow: auto; border: 1px solid var(--line); border-radius: 16px; background: var(--card); }}
table {{ width: 100%; min-width: 760px; border-collapse: separate; border-spacing: 0; font-size: 0.9rem; }}
th {{ position: sticky; top: 0; z-index: 1; background: #f8fbff; color: #475467; font-size: 0.72rem; font-weight: 800; letter-spacing: 0.08em; text-align: left; text-transform: uppercase; }}
th, td {{ padding: 11px 13px; border-bottom: 1px solid var(--line); vertical-align: top; }}
td {{ color: #263247; overflow-wrap: anywhere; }}
tbody tr:hover td {{ background: #f8fbff; }}
tbody tr:last-child td {{ border-bottom: 0; }}
.empty {{ color: var(--muted); text-align: center; }}
code {{ font-family: var(--mono); font-size: 0.88em; }}
pre {{ margin: 0; max-height: 540px; overflow: auto; border-radius: 16px; border: 1px solid #172033; background: #0f172a; color: #dbeafe; padding: 18px; line-height: 1.55; box-shadow: inset 0 1px 0 rgba(255,255,255,0.06); }}
.split-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 18px; }}
.subcard {{ border: 1px solid var(--line); border-radius: 18px; background: #ffffff; padding: 18px; }}
.footer {{ color: var(--muted); font-size: 0.84rem; margin-top: 28px; text-align: center; }}
@media (max-width: 720px) {{
  .page {{ width: min(100% - 28px, 1480px); padding-top: 18px; }}
  .hero, .section-card {{ border-radius: 20px; padding: 20px; }}
  table {{ min-width: 680px; }}
}}
</style>
</head>
<body>
<main class="page">
<header class="hero">
<p class="eyebrow">Combine workspace bundle</p>
<h1>{html_escape(title)}</h1>
<p class="subtitle">{html_escape(subtitle)}</p>
</header>
{body}
<p class="footer">Generated by scripts/build_bundled_workspace.py</p>
</main>
</body>
</html>
"""


def ensure_file(path: Path, hint: str | None = None) -> None:
    if path.exists():
        return
    message = f"Required input not found: {relative_to_repo(path)}"
    if hint:
        message += f". {hint}"
    raise FileNotFoundError(message)


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


def ensure_shared_observables(workspace: RooWorkspace) -> None:
    if not workspace.var("boson_mass"):
        workspace.factory(f"boson_mass[{BOSON_MASS_LOWER},{BOSON_MASS_UPPER}]")
    if not workspace.var("upsilon_mass"):
        workspace.factory(f"upsilon_mass[{UPSILON_MASS_LOWER},{UPSILON_MASS_UPPER}]")

    boson_mass = workspace.var("boson_mass")
    upsilon_mass = workspace.var("upsilon_mass")
    if not boson_mass or not upsilon_mass:
        raise RuntimeError("Failed to create shared observables")

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


def import_pdf_with_shared_observables(
    target_workspace: RooWorkspace,
    source_workspace,
    source_pdf_name: str,
    target_pdf_name: str,
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
        ROOT.RooFit.RecycleConflictNodes(),
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


def import_dataset_with_shared_observables(
    target_workspace: RooWorkspace,
    source_workspace,
    source_dataset_name: str,
    target_dataset_name: str,
) -> str:
    source_data = source_workspace.data(source_dataset_name)
    if source_data is None:
        raise RuntimeError(
            f"Could not find dataset {source_dataset_name!r} in source workspace {source_workspace.GetName()}"
        )

    clone = source_data.Clone(f"{target_dataset_name}_clone")
    import_options = [RooFit.Rename(target_dataset_name)]
    if source_data.get().find("weight"):
        import_options.append(
            RooFit.RenameVariable("weight", f"{target_dataset_name}_weight")
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


def make_signal_specs() -> list[SignalProcessSpec]:
    return [
        SignalProcessSpec(
            process=sample.process,
            state=sample.state,
            signal_workspace_path=REPO_ROOT
            / f"signal_workspace_{sample.inner_file_name}.root",
        )
        for sample in SIGNAL_SAMPLES
    ]


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


def nonres_candidate_alias_name(candidate: dict[str, Any]) -> str:
    family = candidate["pdf_family"]
    if family == "johnson":
        return f"{NON_RESONANT_PROCESS_NAME}_johnson_nominal_pdf"
    order = candidate["scan_order"]
    return f"{NON_RESONANT_PROCESS_NAME}_{family}_order{order}_pdf"


def import_signal_artifacts(
    target_workspace: RooWorkspace,
    signal_spec: SignalProcessSpec,
) -> dict[str, Any]:
    log(f"Importing signal PDF for {signal_spec.combine_name}")
    root_file, source_workspace = open_workspace(
        signal_spec.signal_workspace_path,
        SIGNAL_WORKSPACE_NAME,
    )
    ensure_shared_observables(target_workspace)
    import_pdf_with_shared_observables(
        target_workspace,
        source_workspace,
        SIGNAL_SOURCE_PDF_NAME,
        signal_spec.combine_name,
        f"__{signal_spec.combine_name}__signal",
    )

    signal_norm = RooRealVar(
        f"{signal_spec.combine_name}_norm",
        f"{signal_spec.combine_name}_norm",
        signal_norm_from_workspace(signal_spec.signal_workspace_path),
    )
    signal_norm.setConstant(True)
    getattr(target_workspace, "import")(signal_norm)
    root_file.Close()

    return {
        "process": signal_spec.combine_name,
        "channel_label": signal_spec.channel_label,
        "source_workspace": relative_to_repo(signal_spec.signal_workspace_path),
        "pdf": signal_spec.combine_name,
        "norm_name": signal_norm.GetName(),
        "norm": float(signal_norm.getVal()),
        "norm_status": "fixed",
    }


def import_resonant_artifacts(
    target_workspace: RooWorkspace,
    process: str,
    resonant_summary: dict[str, Any],
) -> dict[str, Any]:
    log(f"Importing resonant background PDF for {process}")
    source_info = RESONANT_CONFIG[process]
    root_file, source_workspace = open_workspace(
        source_info["root_path"],
        source_info["workspace_name"],
    )
    ensure_shared_observables(target_workspace)
    public_process = str(source_info["public_process"])
    import_pdf_with_shared_observables(
        target_workspace,
        source_workspace,
        str(source_info["source_pdf"]),
        public_process,
        f"__{public_process}__resonant",
    )

    initial_norm = float(resonant_summary["process_initial_norms"][process])
    resonant_norm = RooRealVar(
        f"{public_process}_norm",
        f"{public_process}_norm",
        initial_norm,
        0.0,
        max(initial_norm * 10.0, 10.0),
    )
    resonant_norm.setConstant(False)
    getattr(target_workspace, "import")(resonant_norm)
    root_file.Close()

    return {
        "process": public_process,
        "source_workspace": relative_to_repo(source_info["root_path"]),
        "pdf": public_process,
        "norm_name": resonant_norm.GetName(),
        "norm": float(resonant_norm.getVal()),
        "norm_status": "floating",
    }


def import_non_resonant_artifacts(
    target_workspace: RooWorkspace,
    nonres_workspace,
    nonres_selection: dict[str, Any],
) -> dict[str, Any]:
    log("Importing observed data and non-resonant background PDFs")
    ensure_shared_observables(target_workspace)
    import_dataset_with_shared_observables(
        target_workspace,
        nonres_workspace,
        "data_obs",
        "data_obs",
    )

    category = RooCategory("pdfindex", "pdfindex")
    getattr(target_workspace, "import")(category)

    multipdf_inputs = RooArgList()
    candidate_aliases: list[str] = []
    state_mappings: list[dict[str, Any]] = []
    observables = RooArgSet(
        target_workspace.var("upsilon_mass"),
        target_workspace.var("boson_mass"),
    )
    for candidate in nonres_selection["selected_candidates"]:
        alias_name = nonres_candidate_alias_name(candidate)
        import_pdf_with_shared_observables(
            target_workspace,
            nonres_workspace,
            candidate["model_name"],
            alias_name,
            f"__{alias_name}",
        )
        freeze_pdf_params(target_workspace.pdf(alias_name), observables)
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
        raise RuntimeError("No non-resonant candidates were selected")

    multipdf = RooMultiPdf(
        NON_RESONANT_PROCESS_NAME,
        NON_RESONANT_PROCESS_NAME,
        target_workspace.cat("pdfindex"),
        multipdf_inputs,
    )
    getattr(target_workspace, "import")(multipdf)
    target_workspace.cat("pdfindex").setIndex(0)

    nonres_norm = RooRealVar(
        f"{NON_RESONANT_PROCESS_NAME}_norm",
        f"{NON_RESONANT_PROCESS_NAME}_norm",
        float(nonres_selection["johnson_initial_norm"]),
        0.0,
        max(float(nonres_selection["full_entries"]) * 10.0, 10.0),
    )
    nonres_norm.setConstant(False)
    getattr(target_workspace, "import")(nonres_norm)

    return {
        "process": NON_RESONANT_PROCESS_NAME,
        "source_workspace": relative_to_repo(NON_RESONANT_WORKSPACE_PATH),
        "dataset": "data_obs",
        "pdf": NON_RESONANT_PROCESS_NAME,
        "norm_name": nonres_norm.GetName(),
        "norm": float(nonres_norm.getVal()),
        "norm_status": "floating",
        "mode": nonres_selection["mode_key"],
        "candidates": candidate_aliases,
        "state_mappings": state_mappings,
    }


def combined_process_names(signal_specs: list[SignalProcessSpec]) -> list[str]:
    return [spec.combine_name for spec in signal_specs] + [
        NON_RESONANT_PROCESS_NAME,
        RESONANT_H_PROCESS_NAME,
        RESONANT_Z_PROCESS_NAME,
    ]


def combined_process_ids(signal_specs: list[SignalProcessSpec]) -> list[int]:
    signal_ids = [0] + [-index for index in range(1, len(signal_specs))]
    return signal_ids + [1, 2, 3]


def build_systematics_rows(
    signal_specs: list[SignalProcessSpec],
) -> list[tuple[str, str, list[str]]]:
    h_signals = [spec.combine_name for spec in signal_specs if spec.process == "H"]
    z_signals = [spec.combine_name for spec in signal_specs if spec.process == "Z"]
    process_order = combined_process_names(signal_specs)
    h_systematics = {
        name: (pdf_type, signal_value, nonres_value, resonant_value)
        for name, pdf_type, signal_value, nonres_value, resonant_value in COPIED_SYSTEMATICS[
            "H"
        ]
    }
    z_systematics = {
        name: (pdf_type, signal_value, nonres_value, resonant_value)
        for name, pdf_type, signal_value, nonres_value, resonant_value in COPIED_SYSTEMATICS[
            "Z"
        ]
    }

    rows: list[tuple[str, str, list[str]]] = []
    for (
        nuisance_name,
        h_pdf_type,
        h_signal_value,
        _,
        h_resonant_value,
    ) in COPIED_SYSTEMATICS["H"]:
        if nuisance_name not in z_systematics:
            raise RuntimeError(f"Missing Z systematic row for {nuisance_name!r}")
        z_pdf_type, z_signal_value, _, z_resonant_value = z_systematics[nuisance_name]
        if z_pdf_type != h_pdf_type:
            raise RuntimeError(
                f"Systematic {nuisance_name!r} has inconsistent types between H and Z"
            )

        if (h_signal_value, h_resonant_value) == (z_signal_value, z_resonant_value):
            values = []
            for process_name in process_order:
                if process_name in h_signals or process_name in z_signals:
                    values.append(h_signal_value)
                elif process_name == NON_RESONANT_PROCESS_NAME:
                    values.append("-")
                elif process_name == RESONANT_H_PROCESS_NAME:
                    values.append(h_resonant_value)
                elif process_name == RESONANT_Z_PROCESS_NAME:
                    values.append(z_resonant_value)
                else:
                    raise RuntimeError(f"Unsupported process name {process_name!r}")
            rows.append((nuisance_name, h_pdf_type, values))
            continue

        h_values = []
        z_values = []
        for process_name in process_order:
            if process_name in h_signals:
                h_values.append(h_signal_value)
                z_values.append("-")
            elif process_name in z_signals:
                h_values.append("-")
                z_values.append(z_signal_value)
            elif process_name == NON_RESONANT_PROCESS_NAME:
                h_values.append("-")
                z_values.append("-")
            elif process_name == RESONANT_H_PROCESS_NAME:
                h_values.append(h_resonant_value)
                z_values.append("-")
            elif process_name == RESONANT_Z_PROCESS_NAME:
                h_values.append("-")
                z_values.append(z_resonant_value)
            else:
                raise RuntimeError(f"Unsupported process name {process_name!r}")

        rows.append((f"{nuisance_name}_H", h_pdf_type, h_values))
        rows.append((f"{nuisance_name}_Z", z_pdf_type, z_values))

    return rows


def write_datacard(
    output_dir: Path,
    datacard_file_name: str,
    workspace_file_name: str,
    workspace_name: str,
    signal_specs: list[SignalProcessSpec],
) -> tuple[Path, list[str]]:
    process_names = combined_process_names(signal_specs)
    process_ids = combined_process_ids(signal_specs)
    systematics_rows = build_systematics_rows(signal_specs)

    lines = [
        "# Single simultaneous datacard for H/Z -> Upsilon(nS) + photon",
        f"# Bin: {ANALYSIS_BIN_NAME}",
        "imax 1 number of channels",
        f"jmax {len(process_names) - 1} number of processes minus 1",
        "kmax * number of nuisance parameters",
        "------------",
        f"shapes data_obs {ANALYSIS_BIN_NAME} {workspace_file_name} {workspace_name}:data_obs",
    ]
    for signal_spec in signal_specs:
        lines.append(
            f"shapes {signal_spec.combine_name} {ANALYSIS_BIN_NAME} {workspace_file_name} {workspace_name}:{signal_spec.combine_name}"
        )
    lines.extend(
        [
            f"shapes {NON_RESONANT_PROCESS_NAME} {ANALYSIS_BIN_NAME} {workspace_file_name} {workspace_name}:{NON_RESONANT_PROCESS_NAME}",
            f"shapes {RESONANT_H_PROCESS_NAME} {ANALYSIS_BIN_NAME} {workspace_file_name} {workspace_name}:{RESONANT_H_PROCESS_NAME}",
            f"shapes {RESONANT_Z_PROCESS_NAME} {ANALYSIS_BIN_NAME} {workspace_file_name} {workspace_name}:{RESONANT_Z_PROCESS_NAME}",
            "------------",
            f"bin {ANALYSIS_BIN_NAME}",
            "observation -1",
            "------------",
            "bin " + " ".join([ANALYSIS_BIN_NAME] * len(process_names)),
            "process " + " ".join(process_names),
            "process " + " ".join(str(process_id) for process_id in process_ids),
            "rate " + " ".join(["1"] * len(process_names)),
            "------------",
        ]
    )
    for nuisance_name, nuisance_pdf, values in systematics_rows:
        lines.append(f"{nuisance_name} {nuisance_pdf} " + " ".join(values))
    lines.extend(["------------", "pdfindex discrete", ""])

    datacard_path = output_dir / datacard_file_name
    datacard_path.write_text("\n".join(lines), encoding="ascii")
    return datacard_path, lines


def validate_datacard(
    output_dir: Path,
    datacard_file_name: str,
    mass_label: str,
) -> dict[str, Any]:
    datacard_path = output_dir / datacard_file_name
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
        cwd=output_dir,
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
            ".validate",
        ]
        combine_res = subprocess.run(
            combine_cmd,
            cwd=output_dir,
            capture_output=True,
            text=True,
        )
        combine_status = "ok" if combine_res.returncode == 0 else "failed"
        combine_output = (combine_res.stdout or "") + (combine_res.stderr or "")
        combine_returncode = combine_res.returncode

    log_lines = [
        "# Validation log for single simultaneous datacard",
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
    log_path = output_dir / "validation.log"
    log_path.write_text("\n".join(log_lines), encoding="ascii")

    produced_files = [log_path]
    validation_workspace_path = output_dir / "validation_workspace.root"
    if validation_workspace_path.exists():
        produced_files.append(validation_workspace_path)
    combine_root_path = (
        output_dir / f"higgsCombine.validate.AsymptoticLimits.mH{mass_label}.root"
    )
    if combine_root_path.exists():
        produced_files.append(combine_root_path)

    return {
        "text2workspace_status": "ok"
        if text2workspace_res.returncode == 0
        else "failed",
        "combine_status": combine_status,
        "text2workspace_returncode": text2workspace_res.returncode,
        "combine_returncode": combine_returncode,
        "log": relative_to_repo(log_path),
        "produced_files": produced_files,
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
        if name in {"boson_mass", "upsilon_mass"}:
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


def write_summary_report(
    output_dir: Path,
    workspace_name: str,
    workspace_file_name: str,
    datacard_file_name: str,
    workspace: RooWorkspace,
    signal_specs: list[SignalProcessSpec],
    signal_reports: list[dict[str, Any]],
    resonant_reports: list[dict[str, Any]],
    nonres_report: dict[str, Any],
    nonres_selection: dict[str, Any],
    datacard_lines: list[str],
    validation_result: dict[str, Any] | None,
    produced_files: list[Path],
) -> None:
    object_lists = collect_workspace_objects(workspace)
    parameter_rows = collect_parameter_rows(workspace)
    systematics_rows = build_systematics_rows(signal_specs)

    signal_rows = [
        [
            report["process"],
            report["channel_label"],
            report["pdf"],
            report["norm_name"],
            f"{report['norm']:.8g}",
            report["norm_status"],
            report["source_workspace"],
        ]
        for report in signal_reports
    ]
    background_rows = [
        [
            nonres_report["process"],
            nonres_report["pdf"],
            nonres_report["norm_name"],
            f"{nonres_report['norm']:.8g}",
            nonres_report["norm_status"],
            nonres_report["source_workspace"],
        ]
    ] + [
        [
            report["process"],
            report["pdf"],
            report["norm_name"],
            f"{report['norm']:.8g}",
            report["norm_status"],
            report["source_workspace"],
        ]
        for report in resonant_reports
    ]
    state_mapping_rows = [
        [
            state["index"],
            state["workspace_label"],
            state["readable_label"],
            state["pdf_name"],
            state["source_model_name"],
        ]
        for state in nonres_report["state_mappings"]
    ]
    systematics_table_rows = [
        [name, pdf_type, *values] for name, pdf_type, values in systematics_rows
    ]
    validation_rows = (
        [
            [
                validation_result["text2workspace_status"],
                validation_result["combine_status"],
                validation_result["log"],
            ]
        ]
        if validation_result
        else [["skipped", "skipped", "-"]]
    )

    summary_metrics = [
        ("Bundled ROOT", workspace_file_name),
        ("Workspace", workspace_name),
        ("Datacard", datacard_file_name),
        ("Analysis bin", ANALYSIS_BIN_NAME),
        ("Signal PDFs", len(signal_reports)),
        ("Background PDFs", len(background_rows)),
        ("Parameters", len(parameter_rows)),
        ("Non-resonant mode", nonres_selection["mode_key"]),
    ]
    summary_notes = [
        "One shared observed dataset is used for the full simultaneous model.",
        "The six public signal processes are H_1S, H_2S, H_3S, Z_1S, Z_2S, and Z_3S.",
        "Background processes are non_resonant_bkg, resonant_H_bkg, and resonant_Z_bkg.",
        "Signal yields use fixed workspace-side *_norm objects; background yields use floating *_norm objects.",
    ]
    observable_rows = [
        ["upsilon_mass", UPSILON_MASS_LOWER, UPSILON_MASS_UPPER, "GeV", "-"],
        ["boson_mass", BOSON_MASS_LOWER, BOSON_MASS_UPPER, "GeV", "LEFT, MIDDLE, RIGHT, FULL"],
    ]
    imported_objects_html = "<div class=\"split-grid\">" + "\n".join(
        f"""<div class="subcard">
<h3>{html_escape(title)}</h3>
{html_code_block(object_lists[key])}
</div>"""
        for title, key in (
            ("Datasets", "datasets"),
            ("PDFs", "pdfs"),
            ("Functions", "functions"),
            ("Variables", "variables"),
            ("Categories", "categories"),
        )
    ) + "\n</div>"

    sections = [
        html_section(
            "Summary",
            html_metric_grid(summary_metrics) + html_list(summary_notes),
            "Single-card Combine package produced from the persisted upstream fits.",
        ),
        html_section(
            "Observable Layout",
            html_table(["Observable", "Lower", "Upper", "Unit", "Named ranges"], observable_rows),
        ),
        html_section(
            "Signal Processes",
            html_table(
                [
                    "Process",
                    "Legacy channel label",
                    "PDF",
                    "Norm",
                    "Norm value",
                    "Status",
                    "Source workspace",
                ],
                signal_rows,
            ),
        ),
        html_section(
            "Background Processes",
            html_table(
                ["Process", "PDF", "Norm", "Norm value", "Status", "Source workspace"],
                background_rows,
            ),
        ),
        html_section(
            "Non-Resonant Candidate Selection",
            html_table(
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
        ),
        html_section(
            "RooMultiPdf State Mapping",
            html_table(
                [
                    "Index",
                    "Workspace label",
                    "Readable label",
                    "Imported PDF",
                    "Source model",
                ],
                state_mapping_rows,
            ),
        ),
        html_section(
            "Systematics Written To The Single Card",
            html_table(
                [
                    "Nuisance",
                    "Type",
                    "H_1S",
                    "H_2S",
                    "H_3S",
                    "Z_1S",
                    "Z_2S",
                    "Z_3S",
                    NON_RESONANT_PROCESS_NAME,
                    RESONANT_H_PROCESS_NAME,
                    RESONANT_Z_PROCESS_NAME,
                ],
                systematics_table_rows,
            ),
        ),
        html_section("Datacard Contents", html_code_block(datacard_lines)),
        html_section("Imported Objects", imported_objects_html),
        html_section(
            "Parameters",
            html_table(
                ["Name", "Type", "Status", "Value/Label", "Min", "Max"],
                parameter_rows,
            ),
        ),
        html_section(
            "Source Inputs",
            html_list(
                [
                    "Signal workspaces come from the outputs of python3 scripts/signal.py.",
                    f"Resonant background workspaces come from python3 scripts/resonant_background.py and summary {relative_to_repo(RESONANT_SUMMARY_PATH)}.",
                    f"The shared observed dataset and non-resonant candidate PDFs come from {relative_to_repo(NON_RESONANT_WORKSPACE_PATH)} and summary {relative_to_repo(NON_RESONANT_SUMMARY_PATH)}.",
                    f"Legacy lnN values were copied from {relative_to_repo(LEGACY_SYSTEMATICS_SOURCE)} and expanded onto the single-card process layout.",
                ]
            ),
        ),
        html_section(
            "Produced Files",
            html_code_block(sorted(relative_to_repo(path) for path in produced_files)),
        ),
        html_section(
            "Notes",
            html_list(
                [
                    "Validation runs text2workspace.py and then a blind combine -M AsymptoticLimits smoke test by default.",
                    "Use --skip-validation to write only the workspace and parametric datacard.",
                    "The six public signal process names are ready to be mapped to six independent POIs in a later Combine step.",
                ]
            ),
        ),
        html_section(
            "Validation",
            html_table(["text2workspace.py", "combine -M AsymptoticLimits", "Log"], validation_rows),
        ),
    ]
    report_html = html_page(
        "Bundled Single Datacard",
        "A clean overview of the simultaneous H/Z to Upsilon(nS)+photon Combine workspace, datacard, imported objects, nuisance layout, and validation status.",
        "\n".join(sections),
    )
    (output_dir / "bundle_summary.html").write_text(report_html, encoding="ascii")


def verify_required_outputs(signal_specs: list[SignalProcessSpec]) -> None:
    log("Checking persisted upstream outputs")
    for signal_spec in signal_specs:
        ensure_file(
            signal_spec.signal_workspace_path,
            "Run `python3 scripts/signal.py` first.",
        )
        log_kv("signal workspace", relative_to_repo(signal_spec.signal_workspace_path))
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
                f"Non-resonant summary is missing family {family!r}. Re-run `python3 scripts/non_resonant_background.py`."
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
    signal_specs = make_signal_specs()
    produced_files: list[Path] = []

    log("Starting single-card workspace build")
    log_kv("output directory", relative_to_repo(output_dir))
    log_kv("workspace name", args.workspace_name)
    log_kv("workspace file name", args.workspace_file_name)
    log_kv("datacard file name", args.datacard_file_name)
    log_kv("selection mode", "strict" if args.strict_mode else "relaxed")
    log_kv("validation", "disabled" if args.skip_validation else "enabled")

    verify_required_outputs(signal_specs)
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
    ensure_shared_observables(bundled_workspace)

    nonres_root_file, nonres_workspace = open_workspace(
        NON_RESONANT_WORKSPACE_PATH,
        NON_RESONANT_WORKSPACE_NAME,
    )
    nonres_report = import_non_resonant_artifacts(
        bundled_workspace,
        nonres_workspace,
        nonres_selection,
    )
    nonres_root_file.Close()
    log_kv("data_obs", nonres_report["dataset"])
    log_kv("non-resonant pdf", nonres_report["pdf"])

    signal_reports: list[dict[str, Any]] = []
    for signal_spec in signal_specs:
        report = import_signal_artifacts(bundled_workspace, signal_spec)
        signal_reports.append(report)
        log_kv("signal pdf", report["pdf"])
        log_kv("signal norm", f"{report['norm_name']}={report['norm']:.8g}")

    resonant_reports: list[dict[str, Any]] = []
    for process in ("H", "Z"):
        report = import_resonant_artifacts(
            bundled_workspace,
            process,
            resonant_summary,
        )
        resonant_reports.append(report)
        log_kv("resonant pdf", report["pdf"])
        log_kv("resonant norm", f"{report['norm_name']}={report['norm']:.8g}")

    workspace_path = output_dir / args.workspace_file_name
    bundled_workspace.writeToFile(str(workspace_path))
    produced_files.append(workspace_path)
    log(f"Wrote bundled workspace: {relative_to_repo(workspace_path)}")

    datacard_path, datacard_lines = write_datacard(
        output_dir=output_dir,
        datacard_file_name=args.datacard_file_name,
        workspace_file_name=args.workspace_file_name,
        workspace_name=args.workspace_name,
        signal_specs=signal_specs,
    )
    produced_files.append(datacard_path)
    log(f"Wrote simultaneous datacard: {relative_to_repo(datacard_path)}")

    validation_result: dict[str, Any] | None = None
    if not args.skip_validation:
        log("Running validation commands")
        validation_result = validate_datacard(
            output_dir=output_dir,
            datacard_file_name=args.datacard_file_name,
            mass_label=args.signal_mass_label,
        )
        produced_files.extend(validation_result["produced_files"])
        log_kv("text2workspace.py", validation_result["text2workspace_status"])
        log_kv("combine -M AsymptoticLimits", validation_result["combine_status"])

    report_path = output_dir / "bundle_summary.html"
    produced_files.append(report_path)
    write_summary_report(
        output_dir=output_dir,
        workspace_name=args.workspace_name,
        workspace_file_name=args.workspace_file_name,
        datacard_file_name=args.datacard_file_name,
        workspace=bundled_workspace,
        signal_specs=signal_specs,
        signal_reports=signal_reports,
        resonant_reports=resonant_reports,
        nonres_report=nonres_report,
        nonres_selection=nonres_selection,
        datacard_lines=datacard_lines,
        validation_result=validation_result,
        produced_files=produced_files,
    )
    log(f"Wrote summary report: {relative_to_repo(report_path)}")

    log("Produced artifacts")
    for path in sorted(produced_files):
        log_kv("file", relative_to_repo(path))

    log("Bundler setup complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
