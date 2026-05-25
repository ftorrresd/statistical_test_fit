#!/usr/bin/env python3

from __future__ import annotations

import argparse
import html
import json
import math
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
    RooAddPdf,  # type: ignore
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
from statistical_test_fit.display_names import pdf_state_display
from statistical_test_fit.fastplot import fastplot
from statistical_test_fit.signal_modeling import SIGNAL_SAMPLES
from statistical_test_fit.ws_helper import freeze_pdf_params


DEFAULT_OUTPUT_DIR = REPO_ROOT / "datacards"
DEFAULT_WORKSPACE_NAME = "combined_workspace"
DEFAULT_WORKSPACE_FILENAME = "workspace.root"
DEFAULT_DATACARD_FILENAME = "datacard.txt"
DEFAULT_SUMMARY_JSON_FILENAME = "bundle_summary.json"
MUMUGAMMA_PROJECTION_NBINS = 72
MUMU_PROJECTION_NBINS = 80
ANALYSIS_BIN_NAME = "cat1"
NON_RESONANT_PROCESS_NAME = "non_resonant_bkg"
NON_RESONANT_PROJECTION_PDF_NAME = f"{NON_RESONANT_PROCESS_NAME}_johnson_nominal_pdf"
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
YIELD_SYSTEMATICS_PATH = REPO_ROOT / "inputs" / "yields_nevents.json"
MASS_SYSTEMATICS_PATH = REPO_ROOT / "inputs" / "mass_systematics_summary.json"
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

HARDCODED_SYSTEMATICS: dict[str, list[tuple[str, str, str, str, str]]] = {
    "H": [
        ("lumi", "lnN", "1.025", "-", "1.025"),
        ("HZ_xs_sc", "lnN", "0.933/1.046", "-", "0.933/1.046"),
    ],
    "Z": [
        ("lumi", "lnN", "1.025", "-", "1.025"),
        ("HZ_xs_sc", "lnN", "1.033", "-", "1.05"),
    ],
}

YIELD_VARIATION_TO_NUISANCE: tuple[tuple[str, str], ...] = (
    ("pileup", "pu_r"),
    ("trigger_sf", "trg"),
    ("muon_id", "muon_id"),
    ("muon_iso", "muon_iso"),
    ("photon_id", "ph_id"),
    ("photon_electron_veto", "ele_veto"),
    ("pdf_alpha_s_weight", "pdf_alpha_s_weight"),
    ("l1_prefiring", "l1_prefiring"),
)

SIGNAL_MASS_SYSTEMATICS: tuple[tuple[str, str], ...] = (
    ("muon_cor", "CMS_sig_mmg_muon_cor"),
    ("photon_E_scale", "CMS_sig_mmg_photon_E_scale"),
    ("photon_E_smearing", "CMS_sig_mmg_photon_E_smearing"),
)


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
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        help="Skip the text2workspace.py validation command.",
    )
    parser.add_argument(
        "--signal-mass-label",
        default="125",
        help="Mass label passed to text2workspace.py during validation.",
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


def html_plot_card(title: str, png_file_name: str, pdf_file_name: str, caption: str) -> str:
    safe_png_file_name = html_escape(png_file_name)
    safe_pdf_file_name = html_escape(pdf_file_name)
    return f"""<div class="subcard">
<h3>{html_escape(title)}</h3>
<a href="{safe_pdf_file_name}"><img class="plot-object" src="{safe_png_file_name}" alt="{html_escape(title)}"></a>
<p><a class="plot-link" href="{safe_pdf_file_name}">Open PDF</a></p>
<p class="section-intro">{html_escape(caption)}</p>
</div>"""


def html_page(title: str, subtitle: str, body: str) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html_escape(title)}</title>
<style>
:root {{
  color-scheme: dark;
  --bg: #08111f;
  --bg-accent: #12305c;
  --card: #101b2d;
  --card-strong: #14233a;
  --text: #e6edf7;
  --muted: #99a8bc;
  --line: #263750;
  --accent: #7dd3fc;
  --accent-strong: #38bdf8;
  --shadow: 0 24px 70px rgba(0, 0, 0, 0.42);
  --mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
  --sans: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  background:
    radial-gradient(circle at top left, rgba(56, 189, 248, 0.20), transparent 34rem),
    radial-gradient(circle at top right, rgba(99, 102, 241, 0.16), transparent 30rem),
    linear-gradient(180deg, #0b1324 0%, var(--bg) 32rem);
  color: var(--text);
  font-family: var(--sans);
}}
.page {{ width: min(1480px, calc(100% - 48px)); margin: 0 auto; padding: 36px 0 64px; }}
.hero {{
  border-radius: 28px;
  padding: 34px;
  color: var(--text);
  background:
    linear-gradient(135deg, rgba(14, 30, 55, 0.96), rgba(12, 74, 110, 0.82)),
    radial-gradient(circle at top right, rgba(125, 211, 252, 0.28), transparent 26rem);
  border: 1px solid rgba(125, 211, 252, 0.22);
  box-shadow: var(--shadow);
}}
.eyebrow {{ margin: 0 0 10px; font-size: 0.77rem; font-weight: 800; letter-spacing: 0.14em; text-transform: uppercase; opacity: 0.82; }}
h1 {{ margin: 0; font-size: clamp(2rem, 4vw, 3.8rem); line-height: 1.02; letter-spacing: -0.045em; }}
.subtitle {{ max-width: 860px; margin: 16px 0 0; color: #c8d6e8; font-size: 1.04rem; line-height: 1.65; }}
.metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 14px; margin-top: 18px; }}
.metric-grid + ul {{ margin-top: 18px; }}
.metric-card {{ border: 1px solid var(--line); border-radius: 18px; background: linear-gradient(180deg, var(--card-strong), var(--card)); padding: 16px 18px; box-shadow: 0 14px 35px rgba(0, 0, 0, 0.22); }}
.metric-label {{ color: var(--muted); font-size: 0.74rem; font-weight: 800; letter-spacing: 0.09em; text-transform: uppercase; }}
.metric-value {{ margin-top: 8px; color: var(--text); font-size: 1.08rem; font-weight: 750; overflow-wrap: anywhere; }}
.section-card {{ margin-top: 22px; padding: 24px; background: rgba(16, 27, 45, 0.88); border: 1px solid var(--line); border-radius: 24px; box-shadow: 0 18px 46px rgba(0, 0, 0, 0.28); backdrop-filter: blur(10px); }}
.section-heading {{ display: flex; align-items: baseline; justify-content: space-between; gap: 18px; flex-wrap: wrap; margin-bottom: 16px; }}
h2 {{ margin: 0; font-size: 1.22rem; letter-spacing: -0.02em; }}
h3 {{ margin: 22px 0 10px; font-size: 1rem; color: #c8d6e8; }}
.section-intro {{ margin: 0; color: var(--muted); line-height: 1.6; max-width: 820px; }}
ul {{ margin: 0; padding-left: 1.2rem; color: #c8d6e8; line-height: 1.7; }}
.table-wrap {{ width: 100%; overflow: auto; border: 1px solid var(--line); border-radius: 16px; background: var(--card); }}
table {{ width: 100%; min-width: 760px; border-collapse: separate; border-spacing: 0; font-size: 0.9rem; }}
th {{ position: sticky; top: 0; z-index: 1; background: #152238; color: #b6c5d9; font-size: 0.72rem; font-weight: 800; letter-spacing: 0.08em; text-align: left; text-transform: uppercase; }}
th, td {{ padding: 11px 13px; border-bottom: 1px solid var(--line); vertical-align: top; }}
td {{ color: #d8e2f0; overflow-wrap: anywhere; }}
tbody tr:hover td {{ background: rgba(125, 211, 252, 0.08); }}
tbody tr:last-child td {{ border-bottom: 0; }}
.empty {{ color: var(--muted); text-align: center; }}
code {{ font-family: var(--mono); font-size: 0.88em; }}
pre {{ margin: 0; max-height: 540px; overflow: auto; border-radius: 16px; border: 1px solid #314563; background: #050b15; color: #dbeafe; padding: 18px; line-height: 1.55; box-shadow: inset 0 1px 0 rgba(255,255,255,0.06); }}
.split-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 18px; }}
.subcard {{ border: 1px solid var(--line); border-radius: 18px; background: var(--card); padding: 18px; }}
.plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 18px; }}
.plot-object {{ display: block; width: 100%; height: auto; border: 1px solid var(--line); border-radius: 14px; background: #050b15; }}
.plot-link {{ color: var(--accent); font-weight: 700; }}
.footer {{ color: var(--muted); font-size: 0.84rem; margin-top: 28px; text-align: center; }}
@media (max-width: 720px) {{
  .page {{ width: min(100% - 28px, 1480px); padding-top: 18px; }}
  .hero, .section-card {{ border-radius: 20px; padding: 20px; }}
  .plot-grid {{ grid-template-columns: 1fr; }}
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


def clear_default_datacards_dir(output_dir: Path) -> bool:
    if output_dir != DEFAULT_OUTPUT_DIR.resolve():
        return False
    if DEFAULT_OUTPUT_DIR.exists() and DEFAULT_OUTPUT_DIR.is_symlink():
        raise RuntimeError(
            f"Refusing to clear symlinked output directory: {relative_to_repo(DEFAULT_OUTPUT_DIR)}"
        )
    DEFAULT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for child in DEFAULT_OUTPUT_DIR.iterdir():
        if child.is_symlink() or not child.is_dir():
            child.unlink()
        else:
            shutil.rmtree(child)
    return True


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


def signal_mass_state_key(signal_spec: SignalProcessSpec) -> str:
    return signal_spec.state.lower()


def signal_mass_shape_nuisance_name(
    source_name: str,
    signal_spec: SignalProcessSpec,
) -> str:
    for configured_source, nuisance_name in SIGNAL_MASS_SYSTEMATICS:
        if source_name == configured_source:
            return f"{nuisance_name}_{signal_spec.process}"
    raise RuntimeError(f"Unsupported signal mass systematic source {source_name!r}")


def signal_mass_shape_nuisance_rows(
    signal_specs: list[SignalProcessSpec],
) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for signal_spec in signal_specs:
        for source_name, _ in SIGNAL_MASS_SYSTEMATICS:
            nuisance_name = signal_mass_shape_nuisance_name(source_name, signal_spec)
            key = (signal_spec.process, source_name, nuisance_name)
            if key in seen:
                continue
            seen.add(key)
            rows.append((signal_spec.process, source_name, nuisance_name))
    return rows


def validate_signal_mass_systematics(
    mass_summary: dict[str, Any],
    signal_specs: list[SignalProcessSpec],
) -> None:
    for signal_spec in signal_specs:
        process_info = mass_summary.get(signal_spec.process)
        if not isinstance(process_info, dict):
            raise RuntimeError(
                f"Mass systematic summary is missing process {signal_spec.process!r}"
            )
        state_key = signal_mass_state_key(signal_spec)
        state_info = process_info.get(state_key)
        if not isinstance(state_info, dict):
            raise RuntimeError(
                f"Mass systematic summary is missing state {signal_spec.process}/{state_key}"
            )
        for source_name, _ in SIGNAL_MASS_SYSTEMATICS:
            source_info = state_info.get(source_name)
            if not isinstance(source_info, dict):
                raise RuntimeError(
                    f"Mass systematic summary is missing source {signal_spec.process}/{state_key}/{source_name}"
                )
            for key in ("mean_diff_percent", "sigma_diff_percent"):
                value = float(source_info.get(key, float("nan")))
                if not math.isfinite(value) or value < 0.0:
                    raise RuntimeError(
                        f"Mass systematic {signal_spec.process}/{state_key}/{source_name}/{key} is invalid: {value!r}"
                    )


def ensure_signal_mass_nuisance_vars(
    workspace: RooWorkspace,
    signal_spec: SignalProcessSpec,
) -> None:
    for source_name, _ in SIGNAL_MASS_SYSTEMATICS:
        nuisance_name = signal_mass_shape_nuisance_name(source_name, signal_spec)
        nuisance = workspace.var(nuisance_name)
        if not nuisance:
            workspace.factory(f"{nuisance_name}[0,-5,5]")
            nuisance = workspace.var(nuisance_name)
        if not nuisance:
            raise RuntimeError(f"Failed to create signal mass nuisance {nuisance_name}")
        nuisance.setConstant(False)


def make_relative_signal_mass_terms(
    mass_summary: dict[str, Any],
    signal_spec: SignalProcessSpec,
    field_name: str,
) -> list[tuple[str, str, float]]:
    state_info = mass_summary[signal_spec.process][signal_mass_state_key(signal_spec)]
    terms = []
    for source_name, _ in SIGNAL_MASS_SYSTEMATICS:
        nuisance_name = signal_mass_shape_nuisance_name(source_name, signal_spec)
        percent = float(state_info[source_name][field_name])
        terms.append((source_name, nuisance_name, percent / 100.0))
    return terms


def formula_expression_for_relative_shifts(terms: list[tuple[str, str, float]]) -> str:
    expression_terms = []
    for index, (_, _, coefficient) in enumerate(terms, start=1):
        sign = "+" if coefficient >= 0.0 else "-"
        expression_terms.append(f" {sign} {abs(coefficient):.12g}*@{index}")
    return "@0*(1" + "".join(expression_terms) + ")"


def create_signal_mass_syst_function(
    workspace: RooWorkspace,
    function_name: str,
    nominal_var_name: str,
    terms: list[tuple[str, str, float]],
) -> str:
    nominal_var = workspace.var(nominal_var_name)
    if nominal_var is None:
        raise RuntimeError(f"Could not find signal mass parameter {nominal_var_name!r}")

    args = [nominal_var_name] + [nuisance_name for _, nuisance_name, _ in terms]
    expression = formula_expression_for_relative_shifts(terms)
    workspace.factory(f"expr::{function_name}(\"{expression}\", {', '.join(args)})")
    if workspace.function(function_name) is None:
        raise RuntimeError(f"Failed to create signal mass function {function_name}")
    return function_name


def apply_signal_mass_shape_systematics(
    workspace: RooWorkspace,
    signal_spec: SignalProcessSpec,
    raw_pdf_name: str,
    mass_summary: dict[str, Any],
) -> dict[str, Any]:
    ensure_signal_mass_nuisance_vars(workspace, signal_spec)
    process_name = signal_spec.combine_name
    role_suffix = f"___{process_name}__signal"
    mean_var_name = f"mean_boson{role_suffix}"
    sigma_var_name = f"sigma_boson{role_suffix}"
    mean_function_name = f"mean_boson_syst{role_suffix}"
    sigma_function_name = f"sigma_boson_syst{role_suffix}"

    mean_terms = [
        term
        for term in make_relative_signal_mass_terms(
            mass_summary,
            signal_spec,
            "mean_diff_percent",
        )
        if term[0] != "photon_E_smearing"
    ]
    sigma_terms = make_relative_signal_mass_terms(
        mass_summary,
        signal_spec,
        "sigma_diff_percent",
    )
    create_signal_mass_syst_function(
        workspace,
        mean_function_name,
        mean_var_name,
        mean_terms,
    )
    create_signal_mass_syst_function(
        workspace,
        sigma_function_name,
        sigma_var_name,
        sigma_terms,
    )

    edited_pdf = workspace.factory(
        f"EDIT::{process_name}("
        f"{raw_pdf_name},"
        f"{mean_var_name}={mean_function_name},"
        f"{sigma_var_name}={sigma_function_name})"
    )
    if edited_pdf is None or workspace.pdf(process_name) is None:
        raise RuntimeError(f"Failed to create signal PDF {process_name} with mass systematics")

    return {
        "source": relative_to_repo(MASS_SYSTEMATICS_PATH),
        "raw_pdf": raw_pdf_name,
        "mean_nominal_var": mean_var_name,
        "sigma_nominal_var": sigma_var_name,
        "mean_function": mean_function_name,
        "sigma_function": sigma_function_name,
        "mean_terms": [
            {
                "source": source_name,
                "nuisance": nuisance_name,
                "relative_uncertainty": coefficient,
            }
            for source_name, nuisance_name, coefficient in mean_terms
        ],
        "sigma_terms": [
            {
                "source": source_name,
                "nuisance": nuisance_name,
                "relative_uncertainty": coefficient,
            }
            for source_name, nuisance_name, coefficient in sigma_terms
        ],
    }


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
    display = pdf_state_display(
        index=johnson_candidate["workspace_state_index"],
        family=johnson_candidate.get("pdf_family", "johnson"),
        order=johnson_candidate.get("scan_order"),
        selection_role=johnson_candidate.get("selection_role"),
        name=johnson_candidate.get("model_name"),
        workspace_label=johnson_candidate.get("workspace_state_label"),
    )
    johnson_candidate.update(display.metadata())
    johnson_candidate["readable_state_label"] = display.text
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
            display = pdf_state_display(
                index=candidate["workspace_state_index"],
                family=candidate.get("pdf_family", family_name),
                order=candidate.get("scan_order"),
                selection_role=selection_role,
                name=model_name,
                workspace_label=candidate.get("workspace_state_label"),
            )
            candidate.update(display.metadata())
            candidate["readable_state_label"] = display.text
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
    mass_summary: dict[str, Any],
) -> dict[str, Any]:
    log(f"Importing signal PDF for {signal_spec.combine_name}")
    root_file, source_workspace = open_workspace(
        signal_spec.signal_workspace_path,
        SIGNAL_WORKSPACE_NAME,
    )
    ensure_shared_observables(target_workspace)
    raw_pdf_name = f"{signal_spec.combine_name}_nominal_pdf"
    import_pdf_with_shared_observables(
        target_workspace,
        source_workspace,
        SIGNAL_SOURCE_PDF_NAME,
        raw_pdf_name,
        f"__{signal_spec.combine_name}__signal",
    )
    mass_systematics = apply_signal_mass_shape_systematics(
        target_workspace,
        signal_spec,
        raw_pdf_name,
        mass_summary,
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
        "raw_pdf": raw_pdf_name,
        "mass_systematics": mass_systematics,
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
        display = pdf_state_display(
            index=candidate.get("workspace_state_index"),
            family=candidate.get("pdf_family"),
            order=candidate.get("scan_order"),
            selection_role=candidate.get("selection_role"),
            name=alias_name,
            source_model_name=candidate.get("model_name"),
            workspace_label=candidate.get("workspace_state_label"),
        )
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
                "readable_label": candidate.get("display_label", display.text),
                "pdf_name": alias_name,
                "source_model_name": candidate["model_name"],
                "selection_role": candidate["selection_role"],
                "scan_order": candidate["scan_order"],
                "pdf_family": candidate["pdf_family"],
                **display.metadata(),
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


def projection_pdf_name(process_name: str) -> str:
    if process_name == NON_RESONANT_PROCESS_NAME:
        return NON_RESONANT_PROJECTION_PDF_NAME
    return process_name


def edge_aligned_bin_edges(
    boundaries: list[float],
    nbins: int,
) -> list[float]:
    if len(boundaries) < 2:
        raise ValueError("At least two binning boundaries are required")
    if nbins < len(boundaries) - 1:
        raise ValueError("Number of bins must be at least the number of intervals")
    for low, high in zip(boundaries[:-1], boundaries[1:]):
        if low >= high:
            raise ValueError("Binning boundaries must be strictly increasing")

    widths = [high - low for low, high in zip(boundaries[:-1], boundaries[1:])]
    total_width = sum(widths)
    ideal_bins = [nbins * width / total_width for width in widths]
    bins_per_interval = [max(1, int(value)) for value in ideal_bins]

    while sum(bins_per_interval) > nbins:
        candidates = [idx for idx, bins in enumerate(bins_per_interval) if bins > 1]
        idx = min(candidates, key=lambda item: ideal_bins[item] - bins_per_interval[item])
        bins_per_interval[idx] -= 1

    deficit = nbins - sum(bins_per_interval)
    remainder_order = sorted(
        range(len(bins_per_interval)),
        key=lambda idx: ideal_bins[idx] - int(ideal_bins[idx]),
        reverse=True,
    )
    for idx in remainder_order[:deficit]:
        bins_per_interval[idx] += 1

    edges = [boundaries[0]]
    for low, high, interval_bins in zip(
        boundaries[:-1], boundaries[1:], bins_per_interval
    ):
        step = (high - low) / float(interval_bins)
        for bin_idx in range(1, interval_bins + 1):
            edges.append(high if bin_idx == interval_bins else low + step * bin_idx)

    return edges


def projection_process_sideband_norm(
    workspace: RooWorkspace,
    process_name: str,
) -> float:
    pdf_name = projection_pdf_name(process_name)
    pdf = workspace.pdf(pdf_name)
    norm_name = f"{process_name}_norm"
    norm = workspace.var(norm_name)
    boson_mass = workspace.var("boson_mass")
    if pdf is None or norm is None or boson_mass is None:
        raise RuntimeError(
            f"Cannot compute sideband norm for {process_name!r}: missing PDF, norm, or boson_mass"
        )

    norm_value = float(norm.getVal())
    obs_set = RooArgSet(boson_mass)

    sideband_int = pdf.createIntegral(obs_set, RooFit.Range("LEFT,MIDDLE,RIGHT")).getVal()
    full_int = pdf.createIntegral(obs_set, RooFit.Range("FULL")).getVal()
    if full_int <= 0.0:
        raise RuntimeError(f"Full-range integral is non-positive for {pdf_name!r}")
    return norm_value * float(sideband_int) / float(full_int)


def build_bundled_projection_model(
    workspace: RooWorkspace,
    signal_specs: list[SignalProcessSpec],
    projection_range: str | None = None,
) -> tuple[Any, float, list[dict[str, Any]]]:
    pdfs = RooArgList()
    norms = RooArgList()
    local_norms = []
    norm_rows: list[dict[str, Any]] = []
    total_norm = 0.0

    for process_name in combined_process_names(signal_specs):
        pdf_name = projection_pdf_name(process_name)
        pdf = workspace.pdf(pdf_name)
        if pdf is None:
            raise RuntimeError(f"Cannot build projection model: missing PDF {pdf_name!r}")
        norm_name = f"{process_name}_norm"
        norm = workspace.var(norm_name)
        if norm is None:
            raise RuntimeError(f"Cannot build projection model: missing norm {norm_name!r}")
        workspace_norm_value = float(norm.getVal())

        if projection_range is None:
            plot_norm_value = workspace_norm_value
            plot_norm = norm
        else:
            plot_norm_value = projection_process_sideband_norm(workspace, process_name)
            plot_norm = RooRealVar(
                f"projection_{process_name}_norm",
                f"projection_{process_name}_norm",
                plot_norm_value,
            )
            plot_norm.setConstant(True)
            local_norms.append(plot_norm)

        pdfs.add(pdf)
        norms.add(plot_norm)
        total_norm += plot_norm_value
        norm_rows.append(
            {
                "process": process_name,
                "pdf": pdf.GetName(),
                "norm": norm_name,
                "norm_value": workspace_norm_value,
                "projection_norm_value": plot_norm_value,
                "norm_status": "fixed" if norm.isConstant() else "floating",
            }
        )

    model = RooAddPdf(
        "bundled_projection_model",
        "Bundled projection model",
        pdfs,
        norms,
    )
    model._keepalive = {"pdfs": pdfs, "norms": norms, "local_norms": local_norms}
    return model, total_norm, norm_rows


def bundled_projection_components(
    signal_specs: list[SignalProcessSpec],
    total_norm: float,
) -> list[dict[str, Any]]:
    h_signals = [spec.combine_name for spec in signal_specs if spec.process == "H"]
    z_signals = [spec.combine_name for spec in signal_specs if spec.process == "Z"]
    components: list[dict[str, Any]] = [
        {
            "selector": f"{NON_RESONANT_PROJECTION_PDF_NAME},{RESONANT_H_PROCESS_NAME},{RESONANT_Z_PROCESS_NAME}",
            "name": "stack_nonres_resH_resZ",
            "label": "Z resonant #rightarrow #mu#mu#gamma",
            "fill_color": ROOT.kAzure + 5,
            "line_color": ROOT.kAzure + 2,
            "fill_style": 1001,
            "fill_line_width": 0,
            "legend_option": "F",
            "line": False,
        },
        {
            "selector": f"{NON_RESONANT_PROJECTION_PDF_NAME},{RESONANT_H_PROCESS_NAME}",
            "name": "stack_nonres_resH",
            "label": "H resonant #rightarrow #mu#mu#gamma",
            "fill_color": ROOT.kGreen - 2,
            "line_color": ROOT.kGreen + 2,
            "fill_style": 1001,
            "fill_line_width": 0,
            "legend_option": "F",
            "line": False,
        },
        {
            "selector": NON_RESONANT_PROJECTION_PDF_NAME,
            "name": "stack_nonres",
            "label": "Non-resonant",
            "fill_color": ROOT.kOrange - 2,
            "line_color": ROOT.kOrange + 7,
            "fill_style": 1001,
            "fill_line_width": 0,
            "legend_option": "F",
            "line": False,
        },
        {
            "selector": f"{NON_RESONANT_PROJECTION_PDF_NAME},{RESONANT_H_PROCESS_NAME},{RESONANT_Z_PROCESS_NAME}",
            "name": "line_total_background",
            "label": "Total Bkg Model",
            "line_color": ROOT.kBlack,
            "line_style": 1,
            "line_width": 3,
            "fill": False,
            "line": True,
        },
    ]

    if h_signals:
        components.append(
            {
                "selector": ",".join(h_signals),
                "name": "line_H_signal",
                "label": "H #rightarrow #Upsilon + #gamma (#times 10^{5})",
                "line_color": ROOT.kRed + 1,
                "line_style": 2,
                "line_width": 3,
                "fill": False,
                "line": True,
                "normalization": total_norm * 100000.0,
                "draw_after_model": True,
            }
        )
    if z_signals:
        components.append(
            {
                "selector": ",".join(z_signals),
                "name": "line_Z_signal",
                "label": "Z #rightarrow #Upsilon + #gamma (#times 60)",
                "line_color": ROOT.kBlue + 1,
                "line_style": 2,
                "line_width": 3,
                "fill": False,
                "line": True,
                "normalization": total_norm * 60.0,
                "draw_after_model": True,
            }
        )

    return components


def write_bundled_projection_plots(
    workspace: RooWorkspace,
    output_dir: Path,
    signal_specs: list[SignalProcessSpec],
) -> tuple[list[Path], dict[str, Any]]:
    data = workspace.data("data_obs")
    if data is None:
        raise RuntimeError("Cannot draw bundled projections: missing data_obs")
    upsilon_mass = workspace.var("upsilon_mass")
    boson_mass = workspace.var("boson_mass")
    if upsilon_mass is None or boson_mass is None:
        raise RuntimeError("Cannot draw bundled projections: missing shared observables")

    full_model, full_total_norm, full_norm_rows = build_bundled_projection_model(
        workspace,
        signal_specs,
    )
    mumugamma_bin_edges = edge_aligned_bin_edges(
        [
            BOSON_MASS_LOWER,
            LEFT_SIDEBAND_UPPER,
            MIDDLE_SIDEBAND_LOWER,
            MIDDLE_SIDEBAND_UPPER,
            RIGHT_SIDEBAND_LOWER,
            BOSON_MASS_UPPER,
        ],
        MUMUGAMMA_PROJECTION_NBINS,
    )
    mumu_bin_edges = edge_aligned_bin_edges(
        [UPSILON_MASS_LOWER, UPSILON_MASS_UPPER],
        MUMU_PROJECTION_NBINS,
    )
    blind_regions = [
        (LEFT_SIDEBAND_UPPER, MIDDLE_SIDEBAND_LOWER),
        (MIDDLE_SIDEBAND_UPPER, RIGHT_SIDEBAND_LOWER),
    ]
    data_blind_cut = (
        f"((boson_mass<{LEFT_SIDEBAND_UPPER}) || "
        f"(boson_mass>{RIGHT_SIDEBAND_LOWER}) || "
        f"((boson_mass>{MIDDLE_SIDEBAND_LOWER}) && "
        f"(boson_mass<{MIDDLE_SIDEBAND_UPPER})))"
    )
    blinded_data = data.reduce(RooFit.Cut(data_blind_cut))
    if blinded_data is None:
        raise RuntimeError("Could not build blinded data projection dataset")
    data_norm = float(blinded_data.sumEntries())
    plot_specs = [
        {
            "key": "mumugamma",
            "model": full_model,
            "total_norm": full_total_norm,
            "norm_rows": full_norm_rows,
            "components": bundled_projection_components(signal_specs, full_total_norm),
            "projection_range": None,
            "observable": boson_mass,
            "pdf_file_name": "bundle_model_projection_mumugamma.pdf",
            "png_file_name": "bundle_model_projection_mumugamma.png",
            "bin_edges": mumugamma_bin_edges,
            "legend": [0.65, 0.48, 0.92, 0.92],
            "blind_regions": blind_regions,
            "log_y": False,
            "y_headroom_factor": 1.2,
            "caption": "Projection on m_{#mu#mu#gamma}",
        },
        {
            "key": "mumugamma_log",
            "model": full_model,
            "total_norm": full_total_norm,
            "norm_rows": full_norm_rows,
            "components": bundled_projection_components(signal_specs, full_total_norm),
            "projection_range": None,
            "observable": boson_mass,
            "pdf_file_name": "bundle_model_projection_mumugamma_log.pdf",
            "png_file_name": "bundle_model_projection_mumugamma_log.png",
            "bin_edges": mumugamma_bin_edges,
            "legend": [0.65, 0.48, 0.92, 0.92],
            "blind_regions": blind_regions,
            "log_y": True,
            "y_headroom_factor": 1.0,
            "caption": "Projection on m_{#mu#mu#gamma} (log scale)",
        },
        {
            "key": "mumu",
            "model": full_model,
            "total_norm": data_norm,
            "norm_rows": full_norm_rows,
            "components": bundled_projection_components(signal_specs, data_norm),
            "projection_range": None,
            "observable": upsilon_mass,
            "pdf_file_name": "bundle_model_projection_mumu.pdf",
            "png_file_name": "bundle_model_projection_mumu.png",
            "bin_edges": mumu_bin_edges,
            "legend": [0.65, 0.48, 0.92, 0.92],
            "blind_regions": None,
            "log_y": False,
            "y_headroom_factor": 1.5,
            "caption": "Projection on m_{#mu#mu}; data points from the blinded Z and H boson-mass regions are excluded.",
        },
        {
            "key": "mumu_log",
            "model": full_model,
            "total_norm": data_norm,
            "norm_rows": full_norm_rows,
            "components": bundled_projection_components(signal_specs, data_norm),
            "projection_range": None,
            "observable": upsilon_mass,
            "pdf_file_name": "bundle_model_projection_mumu_log.pdf",
            "png_file_name": "bundle_model_projection_mumu_log.png",
            "bin_edges": mumu_bin_edges,
            "legend": [0.65, 0.48, 0.92, 0.92],
            "blind_regions": None,
            "log_y": True,
            "y_headroom_factor": 1.6,
            "caption": "Projection on m_{#mu#mu}; log scale.",
        },
    ]

    produced_paths: list[Path] = []
    plot_report: dict[str, Any] = {
        "model": full_model.GetName(),
        "normalization_source": "workspace *_norm factors for m_{#mu#mu#gamma}; data.sumEntries() for m_{#mu#mu}",
        "total_norm": full_total_norm,
        "data_norm": data_norm,
        "norms": full_norm_rows,
        "plots": {},
    }
    for spec in plot_specs:
        pdf_path = output_dir / str(spec["pdf_file_name"])
        png_path = output_dir / str(spec["png_file_name"])
        fastplot(
            spec["model"],
            blinded_data,
            spec["observable"],
            str(pdf_path),
            extra_filenames=[str(png_path)],
            components=spec["components"],
            nbins=60,
            plot_bin_edges=spec["bin_edges"],
            legend=spec["legend"],
            legend_columns=1,
            legend_text_size=0.022,
            normalization=spec["total_norm"],
            projection_range=spec["projection_range"],
            blind_regions=spec["blind_regions"],
            show_model_curve=False,
            y_headroom_factor=spec.get("y_headroom_factor", 1.5),
            log_y=spec["log_y"],
            show_chi2_ndf=False,
            residual_y_range=(-1.0, 1.0),
        )
        produced_paths.extend([pdf_path, png_path])
        plot_report["plots"][spec["key"]] = {
            "path": relative_to_repo(pdf_path),
            "pdf_path": relative_to_repo(pdf_path),
            "png_path": relative_to_repo(png_path),
            "pdf_file_name": pdf_path.name,
            "png_file_name": png_path.name,
            "observable": spec["observable"].GetName(),
            "caption": spec["caption"],
            "blind_regions": spec["blind_regions"],
            "data_blind_cut": data_blind_cut,
            "bin_edges": spec["bin_edges"],
            "model_norm": spec["total_norm"],
            "data_sum_entries": float(blinded_data.sumEntries()),
            "projection_range": spec["projection_range"],
        }

    return produced_paths, plot_report


def build_hardcoded_systematics_rows(
    signal_specs: list[SignalProcessSpec],
) -> list[tuple[str, str, list[str]]]:
    h_signals = [spec.combine_name for spec in signal_specs if spec.process == "H"]
    z_signals = [spec.combine_name for spec in signal_specs if spec.process == "Z"]
    process_order = combined_process_names(signal_specs)
    h_systematics = {
        name: (pdf_type, signal_value, nonres_value, resonant_value)
        for name, pdf_type, signal_value, nonres_value, resonant_value in HARDCODED_SYSTEMATICS[
            "H"
        ]
    }
    z_systematics = {
        name: (pdf_type, signal_value, nonres_value, resonant_value)
        for name, pdf_type, signal_value, nonres_value, resonant_value in HARDCODED_SYSTEMATICS[
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
    ) in HARDCODED_SYSTEMATICS["H"]:
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


def find_yield_sample(
    yield_summary: dict[str, Any],
    process_name: str,
    prefixes: tuple[str, ...],
) -> str:
    matches = [
        sample_name
        for sample_name in yield_summary
        if any(sample_name.startswith(prefix) for prefix in prefixes)
    ]
    if len(matches) != 1:
        prefix_text = ", ".join(prefixes)
        raise RuntimeError(
            f"Expected exactly one yield entry for {process_name} matching {prefix_text}; found {len(matches)}."
        )
    return matches[0]


def build_yield_process_sample_map(
    yield_summary: dict[str, Any],
    signal_specs: list[SignalProcessSpec],
) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for signal_spec in signal_specs:
        if signal_spec.process == "H":
            prefixes = (f"selected_ggH_HToUps{signal_spec.state}G",)
        elif signal_spec.process == "Z":
            prefixes = (f"selected_ZToUpsilon{signal_spec.state}Gamma",)
        else:
            raise RuntimeError(f"Unsupported signal process {signal_spec.process!r}")
        mapping[signal_spec.combine_name] = find_yield_sample(
            yield_summary,
            signal_spec.combine_name,
            prefixes,
        )

    mapping[RESONANT_H_PROCESS_NAME] = find_yield_sample(
        yield_summary,
        RESONANT_H_PROCESS_NAME,
        ("selected_GluGluHToMuMuG",),
    )
    mapping[RESONANT_Z_PROCESS_NAME] = find_yield_sample(
        yield_summary,
        RESONANT_Z_PROCESS_NAME,
        ("selected_ZGTo2MuG",),
    )
    return mapping


def format_lnn_ratio(value: float) -> str:
    if not math.isfinite(value) or value <= 0.0:
        raise RuntimeError(f"Invalid lnN ratio from yield variations: {value!r}")
    return f"{value:.6g}"


def yield_variation_lnn(
    yield_summary: dict[str, Any],
    sample_name: str,
    variation_name: str,
) -> str:
    sample_info = yield_summary.get(sample_name)
    if not isinstance(sample_info, dict):
        raise RuntimeError(f"Yield entry {sample_name!r} is missing or malformed")

    nominal = float(sample_info.get("n_wgt", 0.0))
    if not math.isfinite(nominal) or nominal <= 0.0:
        raise RuntimeError(f"Yield entry {sample_name!r} has invalid nominal n_wgt={nominal!r}")

    variations = sample_info.get("variations")
    if not isinstance(variations, dict) or variation_name not in variations:
        raise RuntimeError(f"Yield entry {sample_name!r} is missing variation {variation_name!r}")
    variation = variations[variation_name]
    if not isinstance(variation, dict) or "up" not in variation or "down" not in variation:
        raise RuntimeError(
            f"Yield variation {variation_name!r} for {sample_name!r} must contain up/down values"
        )

    down_ratio = float(variation["down"]) / nominal
    up_ratio = float(variation["up"]) / nominal
    return f"{format_lnn_ratio(down_ratio)}/{format_lnn_ratio(up_ratio)}"


def build_yield_systematics_rows(
    signal_specs: list[SignalProcessSpec],
    yield_summary: dict[str, Any],
) -> tuple[list[tuple[str, str, list[str]]], dict[str, Any]]:
    process_order = combined_process_names(signal_specs)
    process_sample_map = build_yield_process_sample_map(yield_summary, signal_specs)
    rows: list[tuple[str, str, list[str]]] = []
    row_reports: list[dict[str, Any]] = []

    for variation_name, nuisance_name in YIELD_VARIATION_TO_NUISANCE:
        values: list[str] = []
        for process_name in process_order:
            if process_name == NON_RESONANT_PROCESS_NAME:
                values.append("-")
                continue
            sample_name = process_sample_map[process_name]
            values.append(yield_variation_lnn(yield_summary, sample_name, variation_name))

        rows.append((nuisance_name, "lnN", values))
        row_reports.append(
            {
                "nuisance": nuisance_name,
                "variation": variation_name,
                "values": dict(zip(process_order, values)),
            }
        )

    return rows, {
        "source": relative_to_repo(YIELD_SYSTEMATICS_PATH),
        "process_samples": process_sample_map,
        "variation_mappings": [
            {"variation": variation_name, "nuisance": nuisance_name}
            for variation_name, nuisance_name in YIELD_VARIATION_TO_NUISANCE
        ],
        "rows": row_reports,
    }


def build_systematics_rows(
    signal_specs: list[SignalProcessSpec],
    yield_summary: dict[str, Any],
) -> tuple[list[tuple[str, str, list[str]]], dict[str, Any]]:
    hardcoded_rows = build_hardcoded_systematics_rows(signal_specs)
    yield_rows, yield_report = build_yield_systematics_rows(signal_specs, yield_summary)
    mass_param_rows = [
        (nuisance_name, "param", ["0", "1"])
        for _, _, nuisance_name in signal_mass_shape_nuisance_rows(signal_specs)
    ]
    return hardcoded_rows + yield_rows + mass_param_rows, yield_report


def format_aligned_datacard_rows(rows: Iterable[Iterable[Any]]) -> list[str]:
    normalized_rows = [[str(cell) for cell in row] for row in rows]
    if not normalized_rows:
        return []

    column_count = max(len(row) for row in normalized_rows)
    widths = [0] * column_count
    for row in normalized_rows:
        for index in range(column_count):
            cell = row[index] if index < len(row) else ""
            widths[index] = max(widths[index], len(cell))

    formatted_rows = []
    for row in normalized_rows:
        padded_row = row + [""] * (column_count - len(row))
        formatted_rows.append(
            "  ".join(
                cell.ljust(widths[index])
                for index, cell in enumerate(padded_row)
            ).rstrip()
        )
    return formatted_rows


def write_datacard(
    output_dir: Path,
    datacard_file_name: str,
    workspace_file_name: str,
    workspace_name: str,
    signal_specs: list[SignalProcessSpec],
    systematics_rows: list[tuple[str, str, list[str]]],
) -> tuple[Path, list[str]]:
    process_names = combined_process_names(signal_specs)
    process_ids = combined_process_ids(signal_specs)

    header_rows = [
        ["imax", "1", "number of channels"],
        ["jmax", str(len(process_names) - 1), "number of processes minus 1"],
        ["kmax", "*", "number of nuisance parameters"],
    ]
    shape_rows = [
        [
            "shapes",
            "data_obs",
            ANALYSIS_BIN_NAME,
            workspace_file_name,
            f"{workspace_name}:data_obs",
        ]
    ]
    for signal_spec in signal_specs:
        shape_rows.append(
            [
                "shapes",
                signal_spec.combine_name,
                ANALYSIS_BIN_NAME,
                workspace_file_name,
                f"{workspace_name}:{signal_spec.combine_name}",
            ]
        )
    shape_rows.extend(
        [
            [
                "shapes",
                NON_RESONANT_PROCESS_NAME,
                ANALYSIS_BIN_NAME,
                workspace_file_name,
                f"{workspace_name}:{NON_RESONANT_PROCESS_NAME}",
            ],
            [
                "shapes",
                RESONANT_H_PROCESS_NAME,
                ANALYSIS_BIN_NAME,
                workspace_file_name,
                f"{workspace_name}:{RESONANT_H_PROCESS_NAME}",
            ],
            [
                "shapes",
                RESONANT_Z_PROCESS_NAME,
                ANALYSIS_BIN_NAME,
                workspace_file_name,
                f"{workspace_name}:{RESONANT_Z_PROCESS_NAME}",
            ],
        ]
    )
    observation_rows = [
        ["bin", ANALYSIS_BIN_NAME],
        ["observation", "-1"],
    ]
    model_rows = [
        ["bin", "", *([ANALYSIS_BIN_NAME] * len(process_names))],
        ["process", "", *process_names],
        ["process", "", *(str(process_id) for process_id in process_ids)],
        ["rate", "", *("1" for _ in process_names)],
        *(
            [nuisance_name, nuisance_pdf, *values]
            for nuisance_name, nuisance_pdf, values in systematics_rows
        ),
        ["pdfindex", "discrete"],
    ]
    formatted_model_rows = format_aligned_datacard_rows(model_rows)
    rate_row_count = 4

    lines = [
        "# Single simultaneous datacard for H/Z -> Upsilon(nS) + photon",
        f"# Bin: {ANALYSIS_BIN_NAME}",
        *format_aligned_datacard_rows(header_rows),
        "------------",
        *format_aligned_datacard_rows(shape_rows),
        "------------",
        *format_aligned_datacard_rows(observation_rows),
        "------------",
        *formatted_model_rows[:rate_row_count],
        "------------",
        *formatted_model_rows[rate_row_count:-1],
        "------------",
        formatted_model_rows[-1],
        "",
    ]

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

    log_lines = [
        "# text2workspace.py validation log for single simultaneous datacard",
        "",
        "## text2workspace.py",
        f"Command: {' '.join(text2workspace_cmd)}",
        f"Return code: {text2workspace_res.returncode}",
        "",
        "```text",
        (text2workspace_res.stdout or "") + (text2workspace_res.stderr or ""),
        "```",
        "",
    ]
    log_path = output_dir / "validation.log"
    log_path.write_text("\n".join(log_lines), encoding="ascii")

    produced_files = [log_path]
    validation_workspace_path = output_dir / "validation_workspace.root"
    if validation_workspace_path.exists():
        produced_files.append(validation_workspace_path)
    return {
        "text2workspace_status": "ok"
        if text2workspace_res.returncode == 0
        else "failed",
        "text2workspace_returncode": text2workspace_res.returncode,
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


def read_datacard_lines(datacard_path: Path) -> list[str]:
    ensure_file(datacard_path)
    return datacard_path.read_text(encoding="ascii").splitlines()


def is_integer_token(value: str) -> bool:
    try:
        int(value)
    except ValueError:
        return False
    return True


def parse_shape_reference(reference: str) -> tuple[str | None, str]:
    if ":" not in reference:
        return None, reference
    workspace_name, object_name = reference.split(":", 1)
    return workspace_name, object_name


def parse_datacard(datacard_lines: list[str]) -> dict[str, Any]:
    shapes: dict[str, dict[str, Any]] = {}
    process_rows: list[list[str]] = []
    process_bins: list[str] = []
    observation_bin = ""
    observation: list[str] = []
    rates: list[str] = []
    nuisance_rows: list[dict[str, Any]] = []
    discrete_rows: list[dict[str, Any]] = []

    for line_number, line in enumerate(datacard_lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or set(stripped) == {"-"}:
            continue
        parts = stripped.split()
        if not parts:
            continue

        keyword = parts[0]
        if keyword == "shapes" and len(parts) >= 5:
            workspace_ref, object_name = parse_shape_reference(parts[4])
            shapes[parts[1]] = {
                "process": parts[1],
                "bin": parts[2],
                "file": parts[3],
                "workspace": workspace_ref,
                "object": object_name,
                "reference": parts[4],
                "line": line_number,
            }
            continue

        if keyword == "bin":
            if len(parts) == 2 and not observation:
                observation_bin = parts[1]
            elif len(parts) > 2:
                process_bins = parts[1:]
            continue

        if keyword == "observation":
            observation = parts[1:]
            continue

        if keyword == "process":
            process_rows.append(parts[1:])
            continue

        if keyword == "rate":
            rates = parts[1:]
            continue

        if len(parts) >= 2 and parts[1] in {"lnN", "param", "discrete"}:
            row = {
                "name": parts[0],
                "type": parts[1],
                "values": parts[2:],
                "line": line_number,
                "raw": stripped,
            }
            if parts[1] == "discrete":
                discrete_rows.append(row)
            else:
                nuisance_rows.append(row)

    process_names: list[str] = []
    process_ids: list[str] = []
    for row in process_rows:
        if row and all(is_integer_token(value) for value in row):
            process_ids = row
        elif row:
            process_names = row

    return {
        "shapes": shapes,
        "observation_bin": observation_bin,
        "observation": observation,
        "process_bins": process_bins,
        "process_names": process_names,
        "process_ids": process_ids,
        "rates": rates,
        "nuisance_rows": nuisance_rows,
        "lnn_rows": [row for row in nuisance_rows if row["type"] == "lnN"],
        "param_rows": [row for row in nuisance_rows if row["type"] == "param"],
        "discrete_rows": discrete_rows,
    }


def workspace_has_shape_object(
    workspace: RooWorkspace,
    process_name: str,
    object_name: str,
) -> bool:
    if process_name == "data_obs":
        return bool(workspace.data(object_name))
    return bool(workspace.pdf(object_name))


def validation_row(check: str, ok: bool, detail: str) -> list[str]:
    return [check, "ok" if ok else "failed", detail]


def validate_persisted_bundle(
    workspace: RooWorkspace,
    workspace_path: Path,
    datacard_path: Path,
    workspace_name: str,
    datacard_info: dict[str, Any],
) -> tuple[list[list[str]], str]:
    rows: list[list[str]] = []
    shapes = datacard_info["shapes"]
    process_names = datacard_info["process_names"]

    rows.append(
        validation_row(
            "workspace file",
            workspace_path.exists(),
            relative_to_repo(workspace_path),
        )
    )
    rows.append(
        validation_row(
            "datacard file",
            datacard_path.exists(),
            relative_to_repo(datacard_path),
        )
    )
    rows.append(
        validation_row(
            "workspace name",
            workspace.GetName() == workspace_name,
            workspace.GetName(),
        )
    )

    expected_shape_processes = ["data_obs", *process_names]
    for process_name in expected_shape_processes:
        shape = shapes.get(process_name)
        if shape is None:
            rows.append(validation_row(f"shape {process_name}", False, "missing from datacard"))
            continue
        file_ok = shape["file"] == workspace_path.name
        workspace_ok = shape["workspace"] == workspace_name
        object_ok = workspace_has_shape_object(workspace, process_name, shape["object"])
        rows.append(
            validation_row(
                f"shape {process_name}",
                file_ok and workspace_ok and object_ok,
                f"{shape['file']} {shape['reference']}",
            )
        )

    for process_name in process_names:
        norm_name = f"{process_name}_norm"
        norm = workspace.var(norm_name)
        rows.append(
            validation_row(
                f"norm {process_name}",
                bool(norm),
                norm_name if norm else "missing",
            )
        )

    for row in datacard_info["param_rows"]:
        param = workspace.var(row["name"])
        rows.append(
            validation_row(
                f"param {row['name']}",
                bool(param),
                "workspace variable found" if param else "missing from workspace",
            )
        )

    for row in datacard_info["discrete_rows"]:
        category = workspace.cat(row["name"])
        rows.append(
            validation_row(
                f"discrete {row['name']}",
                bool(category),
                "workspace category found" if category else "missing from workspace",
            )
        )

    row_lengths_ok = True
    for row_name, values in (
        ("process bins", datacard_info["process_bins"]),
        ("process ids", datacard_info["process_ids"]),
        ("rates", datacard_info["rates"]),
    ):
        ok = len(values) == len(process_names)
        row_lengths_ok = row_lengths_ok and ok
        rows.append(
            validation_row(
                row_name,
                ok,
                f"{len(values)} values for {len(process_names)} processes",
            )
        )

    status = "ok" if all(row[1] == "ok" for row in rows) and row_lengths_ok else "failed"
    return rows, status


def write_summary_report(
    output_dir: Path,
    workspace_name: str,
    workspace_path: Path,
    datacard_path: Path,
    signal_specs: list[SignalProcessSpec],
    signal_reports: list[dict[str, Any]],
    resonant_reports: list[dict[str, Any]],
    nonres_report: dict[str, Any],
    nonres_selection: dict[str, Any],
    yield_systematics_report: dict[str, Any],
    validation_result: dict[str, Any] | None,
    projection_plot_report: dict[str, Any],
    produced_files: list[Path],
) -> None:
    workspace_file_name = workspace_path.name
    datacard_file_name = datacard_path.name
    datacard_lines = read_datacard_lines(datacard_path)
    datacard_info = parse_datacard(datacard_lines)
    root_file, persisted_workspace = open_workspace(workspace_path, workspace_name)
    object_lists = collect_workspace_objects(persisted_workspace)
    parameter_rows = collect_parameter_rows(persisted_workspace)
    bundle_validation_rows, bundle_validation_status = validate_persisted_bundle(
        persisted_workspace,
        workspace_path,
        datacard_path,
        workspace_name,
        datacard_info,
    )

    signal_report_by_process = {report["process"]: report for report in signal_reports}
    background_report_by_process = {
        nonres_report["process"]: nonres_report,
        **{report["process"]: report for report in resonant_reports},
    }

    def norm_row_details(process_name: str) -> tuple[str, str, str]:
        norm_name = f"{process_name}_norm"
        norm = persisted_workspace.var(norm_name)
        if not norm:
            return norm_name, "missing", "missing"
        return (
            norm_name,
            f"{norm.getVal():.8g}",
            "fixed" if norm.isConstant() else "floating",
        )

    signal_rows = [
        [
            process_name,
            signal_report_by_process[process_name]["channel_label"],
            datacard_info["shapes"].get(process_name, {}).get("object", "missing"),
            *norm_row_details(process_name),
            signal_report_by_process[process_name]["source_workspace"],
        ]
        for process_name in datacard_info["process_names"]
        if process_name in signal_report_by_process
    ]
    signal_mass_systematics_rows = []
    for report in signal_reports:
        mass_report = report["mass_systematics"]
        mean_terms = {item["source"]: item for item in mass_report["mean_terms"]}
        sigma_terms = {item["source"]: item for item in mass_report["sigma_terms"]}
        for source_name, _ in SIGNAL_MASS_SYSTEMATICS:
            mean_term = mean_terms.get(source_name)
            sigma_term = sigma_terms.get(source_name)
            if mean_term is None and sigma_term is None:
                continue
            nuisance_name = (mean_term or sigma_term)["nuisance"]
            signal_mass_systematics_rows.append(
                [
                    report["process"],
                    source_name,
                    nuisance_name,
                    "-"
                    if mean_term is None
                    else f"{100.0 * mean_term['relative_uncertainty']:.6g}",
                    "-"
                    if sigma_term is None
                    else f"{100.0 * sigma_term['relative_uncertainty']:.6g}",
                    mass_report["mean_function"],
                    mass_report["sigma_function"],
                ]
            )
    background_rows = [
        [
            process_name,
            datacard_info["shapes"].get(process_name, {}).get("object", "missing"),
            *norm_row_details(process_name),
            background_report_by_process.get(process_name, {}).get("source_workspace", "-"),
        ]
        for process_name in datacard_info["process_names"]
        if process_name not in signal_report_by_process
    ]
    root_file.Close()

    state_mapping_rows = [
        [
            state["index"],
            state["workspace_label"],
            state.get("display_label", state["readable_label"]),
            state.get("display_latex", "-"),
            state["pdf_name"],
            state["source_model_name"],
        ]
        for state in nonres_report["state_mappings"]
    ]
    systematics_table_rows = [
        [row["name"], row["type"], *row["values"]]
        for row in datacard_info["lnn_rows"]
    ]
    param_table_rows = [
        [row["name"], row["type"], *row["values"]]
        for row in datacard_info["param_rows"]
    ]
    discrete_table_rows = [
        [row["name"], row["type"], *row["values"]]
        for row in datacard_info["discrete_rows"]
    ]
    process_order = datacard_info["process_names"]
    yield_process_rows = [
        [process_name, yield_systematics_report["process_samples"][process_name]]
        for process_name in process_order
        if process_name in yield_systematics_report["process_samples"]
    ]
    yield_mapping_rows = [
        [item["variation"], item["nuisance"]]
        for item in yield_systematics_report["variation_mappings"]
    ]
    validation_rows = (
        [
            [
                validation_result["text2workspace_status"],
                validation_result["log"],
            ]
        ]
        if validation_result
        else [["skipped", "-"]]
    )
    projection_norm_rows = [
        [
            row["process"],
            row["pdf"],
            row["norm"],
            f"{row['norm_value']:.8g}",
            row["norm_status"],
        ]
        for row in projection_plot_report.get("norms", [])
    ]
    projection_plot_cards = []
    for key, title in (
        ("mumugamma", "m_{#mu#mu#gamma} Projection"),
        ("mumu", "m_{#mu#mu} Projection"),
    ):
        plot_info = projection_plot_report.get("plots", {}).get(key)
        if not plot_info:
            continue
        projection_plot_cards.append(
            html_plot_card(
                title,
                plot_info["png_file_name"],
                plot_info["pdf_file_name"],
                plot_info["caption"],
            )
        )
    projection_plots_html = (
        html_metric_grid(
            [
                ("Projection model", projection_plot_report.get("model", "-")),
                ("Total *_norm", f"{projection_plot_report.get('total_norm', 0.0):.8g}"),
                ("Plot normalizations", projection_plot_report.get("normalization_source", "-")),
            ]
        )
        + '<div class="plot-grid">'
        + "\n".join(projection_plot_cards)
        + "\n</div>"
        + "<h3>Norm factors used</h3>"
        + html_table(["Process", "PDF", "Norm", "Norm value", "Status"], projection_norm_rows)
    )

    summary_json_payload = {
        "schema_version": 1,
        "summary_source": "persisted_workspace_and_datacard",
        "workspace_name": workspace_name,
        "workspace_file_name": workspace_file_name,
        "datacard_file_name": datacard_file_name,
        "analysis_bin": ANALYSIS_BIN_NAME,
        "non_resonant_mode": nonres_selection["mode_key"],
        "datacard": datacard_info,
        "workspace_objects": object_lists,
        "bundle_validation": {
            "status": bundle_validation_status,
            "checks": [
                {"check": row[0], "status": row[1], "detail": row[2]}
                for row in bundle_validation_rows
            ],
        },
        "signals": signal_reports,
        "backgrounds": [nonres_report, *resonant_reports],
        "non_resonant_selection": {
            key: value
            for key, value in nonres_selection.items()
            if key != "selected_candidates"
        },
        "non_resonant_state_mappings": nonres_report["state_mappings"],
        "yield_systematics": yield_systematics_report,
        "signal_mass_systematics": {
            "source": relative_to_repo(MASS_SYSTEMATICS_PATH),
            "nuisance_parameters": [
                {
                    "process": process_name,
                    "source": source_name,
                    "nuisance": nuisance_name,
                }
                for process_name, source_name, nuisance_name in signal_mass_shape_nuisance_rows(
                    signal_specs
                )
            ],
            "signals": {
                report["process"]: report["mass_systematics"]
                for report in signal_reports
            },
        },
        "projection_plots": projection_plot_report,
        "validation": (
            {
                **validation_result,
                "produced_files": [
                    relative_to_repo(path) for path in validation_result.get("produced_files", [])
                ],
            }
            if validation_result
            else None
        ),
        "produced_files": sorted(relative_to_repo(path) for path in produced_files),
    }
    (output_dir / DEFAULT_SUMMARY_JSON_FILENAME).write_text(
        json.dumps(summary_json_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    summary_metrics = [
        ("Bundled ROOT", workspace_file_name),
        ("Workspace", workspace_name),
        ("Datacard", datacard_file_name),
        ("Analysis bin", ANALYSIS_BIN_NAME),
        ("Signal PDFs", len(signal_rows)),
        ("Background PDFs", len(background_rows)),
        ("Parameters", len(parameter_rows)),
        ("Bundle validation", bundle_validation_status),
        ("Non-resonant mode", nonres_selection["mode_key"]),
    ]
    summary_notes = [
        "This report reopens the written ROOT workspace and rereads the written datacard before summarizing them.",
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
            "Single-card Combine package inspected from persisted artifacts after generation.",
        ),
        html_section(
            "Bundle Validation",
            html_table(["Check", "Status", "Detail"], bundle_validation_rows),
            "These checks compare datacard references with objects found in the written RooWorkspace.",
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
                    "PDF from datacard",
                    "Norm",
                    "Norm value",
                    "Status",
                    "Source workspace",
                ],
                signal_rows,
            ),
        ),
        html_section(
            "Signal Mass Shape Systematics",
            html_table(
                [
                    "Process",
                    "Source",
                    "Nuisance",
                    "Mean rel. unc. (%)",
                    "Sigma rel. unc. (%)",
                    "Mean function",
                    "Sigma function",
                ],
                signal_mass_systematics_rows,
            ),
            "Gaussian-constrained param nuisances modify the imported signal boson DCB mean and sigma through RooFormulaVar replacements.",
        ),
        html_section(
            "Background Processes",
            html_table(
                ["Process", "PDF from datacard", "Norm", "Norm value", "Status", "Source workspace"],
                background_rows,
            ),
        ),
        html_section(
            "Model Projections",
            projection_plots_html,
            "Fastplot projections of a temporary RooAddPdf assembled from the bundled workspace PDFs and their *_norm factors.",
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
                    "LaTeX label",
                    "Imported PDF",
                    "Source model",
                ],
                state_mapping_rows,
            ),
        ),
        html_section(
            "lnN Nuisances From Datacard",
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
        html_section(
            "Param Nuisances From Datacard",
            html_table(["Nuisance", "Type", "Mean", "Sigma"], param_table_rows),
        ),
        html_section(
            "Discrete Nuisances From Datacard",
            html_table(["Nuisance", "Type"], discrete_table_rows),
        ),
        html_section(
            "Yield Systematics Input",
            html_metric_grid(
                [
                    ("Source", yield_systematics_report["source"]),
                    ("Derived nuisances", len(yield_systematics_report["variation_mappings"])),
                ]
            )
            + "<h3>Process mapping</h3>"
            + html_table(["Process", "Yield sample"], yield_process_rows)
            + "<h3>Variation mapping</h3>"
            + html_table(["Yield variation", "Datacard nuisance"], yield_mapping_rows),
            "Asymmetric lnN values are formatted as (variation_down / nominal)/(variation_up / nominal).",
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
                    f"Yield-based asymmetric lnN values come from {relative_to_repo(YIELD_SYSTEMATICS_PATH)}; lumi and HZ_xs_sc remain embedded in this script.",
                    f"Signal mean/sigma shape params come from {relative_to_repo(MASS_SYSTEMATICS_PATH)} and are constrained by datacard param rows.",
                    f"Original legacy lnN reference path: {relative_to_repo(LEGACY_SYSTEMATICS_SOURCE)}.",
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
                    "Validation runs text2workspace.py by default.",
                    "Use --skip-validation to write only the workspace and parametric datacard.",
                    "The six public signal process names are ready to be mapped to six independent POIs in a later Combine step.",
                ]
            ),
        ),
        html_section(
            "Validation",
            html_table(["text2workspace.py", "Log"], validation_rows),
        ),
    ]
    report_html = html_page(
        "Bundled Single Datacard",
        "A persisted-artifact validation report for the simultaneous H/Z to Upsilon(nS)+photon Combine workspace, datacard, imported objects, nuisance layout, and validation status.",
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
    ensure_file(
        YIELD_SYSTEMATICS_PATH,
        "Provide selected-event nominal/variation yields for datacard lnN nuisances.",
    )
    ensure_file(
        MASS_SYSTEMATICS_PATH,
        "Provide signal mass mean/sigma percent variations for datacard param nuisances.",
    )
    log_kv(
        "resonant workspace (H)", relative_to_repo(RESONANT_CONFIG["H"]["root_path"])
    )
    log_kv(
        "resonant workspace (Z)", relative_to_repo(RESONANT_CONFIG["Z"]["root_path"])
    )
    log_kv("resonant summary", relative_to_repo(RESONANT_SUMMARY_PATH))
    log_kv("non-resonant workspace", relative_to_repo(NON_RESONANT_WORKSPACE_PATH))
    log_kv("non-resonant summary", relative_to_repo(NON_RESONANT_SUMMARY_PATH))
    log_kv("yield systematics", relative_to_repo(YIELD_SYSTEMATICS_PATH))
    log_kv("signal mass systematics", relative_to_repo(MASS_SYSTEMATICS_PATH))


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
    signal_specs = make_signal_specs()
    produced_files: list[Path] = []

    log("Starting single-card workspace build")
    if clear_default_datacards_dir(output_dir):
        log(f"Cleared generated output directory: {relative_to_repo(DEFAULT_OUTPUT_DIR)}")
    output_dir.mkdir(parents=True, exist_ok=True)
    log_kv("output directory", relative_to_repo(output_dir))
    log_kv("workspace name", args.workspace_name)
    log_kv("workspace file name", args.workspace_file_name)
    log_kv("datacard file name", args.datacard_file_name)
    log_kv("selection mode", "strict" if args.strict_mode else "relaxed")
    log_kv("validation", "disabled" if args.skip_validation else "enabled")

    verify_required_outputs(signal_specs)
    nonres_summary = load_json(NON_RESONANT_SUMMARY_PATH)
    resonant_summary = load_json(RESONANT_SUMMARY_PATH)
    yield_summary = load_json(YIELD_SYSTEMATICS_PATH)
    mass_summary = load_json(MASS_SYSTEMATICS_PATH)
    validate_summary_shapes(nonres_summary, resonant_summary)
    validate_signal_mass_systematics(mass_summary, signal_specs)
    nonres_selection = build_nonres_candidate_selection(
        nonres_summary, relaxed_mode=not args.strict_mode
    )
    systematics_rows, yield_systematics_report = build_systematics_rows(
        signal_specs,
        yield_summary,
    )

    log("Selected non-resonant candidates for the final RooMultiPdf")
    for row in nonres_selection["selection_rows"]:
        log_kv(
            f"{row[0]} {row[2]}",
            f"model={row[3]}, order={row[4]}, init_norm={row[5]}",
        )
    log("Built datacard lnN nuisance rows")
    for variation_mapping in yield_systematics_report["variation_mappings"]:
        log_kv(
            variation_mapping["nuisance"],
            f"from {variation_mapping['variation']} in {relative_to_repo(YIELD_SYSTEMATICS_PATH)}",
        )
    log("Configured signal mass shape nuisance parameters")
    for process_name, source_name, nuisance_name in signal_mass_shape_nuisance_rows(
        signal_specs
    ):
        log_kv(
            nuisance_name,
            f"{process_name} from {source_name} in {relative_to_repo(MASS_SYSTEMATICS_PATH)}",
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
        report = import_signal_artifacts(bundled_workspace, signal_spec, mass_summary)
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
        systematics_rows=systematics_rows,
    )
    produced_files.append(datacard_path)
    log(f"Wrote simultaneous datacard: {relative_to_repo(datacard_path)}")

    log("Drawing bundled model projections")
    projection_plot_paths, projection_plot_report = write_bundled_projection_plots(
        bundled_workspace,
        output_dir,
        signal_specs,
    )
    produced_files.extend(projection_plot_paths)
    for plot_path in projection_plot_paths:
        log_kv("projection plot", relative_to_repo(plot_path))

    validation_result: dict[str, Any] | None = None
    if not args.skip_validation:
        log("Running text2workspace.py validation")
        validation_result = validate_datacard(
            output_dir=output_dir,
            datacard_file_name=args.datacard_file_name,
            mass_label=args.signal_mass_label,
        )
        produced_files.extend(validation_result["produced_files"])
        log_kv("text2workspace.py", validation_result["text2workspace_status"])

    summary_json_path = output_dir / DEFAULT_SUMMARY_JSON_FILENAME
    produced_files.append(summary_json_path)
    report_path = output_dir / "bundle_summary.html"
    produced_files.append(report_path)
    write_summary_report(
        output_dir=output_dir,
        workspace_name=args.workspace_name,
        workspace_path=workspace_path,
        datacard_path=datacard_path,
        signal_specs=signal_specs,
        signal_reports=signal_reports,
        resonant_reports=resonant_reports,
        nonres_report=nonres_report,
        nonres_selection=nonres_selection,
        yield_systematics_report=yield_systematics_report,
        validation_result=validation_result,
        projection_plot_report=projection_plot_report,
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
