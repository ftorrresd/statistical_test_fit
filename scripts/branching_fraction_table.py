#!/usr/bin/env python3

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


DEFAULT_LIMITS_DIR = REPO_ROOT / "datacards" / "blind_limits"
DEFAULT_IMAGE = "docker://ghcr.io/xu-cheng/texlive-full:latest"
DEFAULT_OUTPUT_NAME = "branching_fraction_limits"
DEFAULT_SIG_FIGS = 3
QUANTILE_LABELS = (
    ("0.16", r"Expected $-1\sigma$"),
    ("0.5", r"Expected median"),
    ("0.84", r"Expected $+1\sigma$"),
)
PROCESS_ORDER = ("H_1S", "H_2S", "H_3S", "Z_1S", "Z_2S", "Z_3S")
SCHEME_PROCESS_ORDER = {
    "six_poi": PROCESS_ORDER,
    "h_grouped": tuple(process for process in PROCESS_ORDER if process.startswith("H_")),
    "z_grouped": tuple(process for process in PROCESS_ORDER if process.startswith("Z_")),
}
THEORY_BRANCHING_FRACTIONS = {
    "H_1S": 5.22e-9,
    "H_2S": 1.42e-9,
    "H_3S": 0.91e-9,
    "Z_1S": 4.8e-8,
    "Z_2S": 2.44e-8,
    "Z_3S": 1.88e-8,
}


@dataclass(frozen=True)
class TableRow:
    process: str
    poi: str
    theory_bf: float
    observed_bf_limit: float | None
    expected_bf_limits: dict[str, float | None]


METHOD_LABELS = {
    "asymptotic": "Asymptotic Limits",
    "hybrid_lhc": "HybridNew Limits",
}
DEFAULT_METHODS = ("asymptotic", "hybrid_lhc")


@dataclass(frozen=True)
class TableSection:
    scheme: str
    method: str
    title: str
    caption: str
    label: str
    rows: list[TableRow]


SECTION_CONFIG = {
    "h_grouped": {
        "title": "Grouped H signal-strength limit",
        "caption": "Grouped H signal-strength limits translated to branching fractions.",
        "label_base": "tab:bf_limits_h_grouped",
    },
    "z_grouped": {
        "title": "Grouped Z signal-strength limit",
        "caption": "Grouped Z signal-strength limits translated to branching fractions.",
        "label_base": "tab:bf_limits_z_grouped",
    },
    "six_poi": {
        "title": "Individual signal-strength limits",
        "caption": "Individual signal-strength limits translated to branching fractions.",
        "label_base": "tab:bf_limits_six_poi",
    },
}
DEFAULT_TABLE_SCHEMES = ("h_grouped", "z_grouped", "six_poi")


def log(message: str) -> None:
    print(f"[branching_fraction_table] {message}")


def repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def ensure_file(path: Path, message: str) -> None:
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"{message}: {path}")


def latest_summary_path() -> Path:
    summaries = sorted(
        DEFAULT_LIMITS_DIR.glob("*/blind_limits_summary.json"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if not summaries:
        raise FileNotFoundError(
            f"No blind limit summary found under {DEFAULT_LIMITS_DIR}. Run scripts/limits.py first or pass --summary."
        )
    return summaries[0]


def normalize_quantile_key(value: str | float) -> str:
    return f"{float(value):.6g}"


def first_value(values: Any, warnings: list[str], context: str) -> float | None:
    if values is None:
        return None
    if isinstance(values, (int, float)):
        return float(values)
    if not isinstance(values, list) or not values:
        return None
    if len(values) > 1:
        warnings.append(f"{context}: found {len(values)} values; using the first one.")
    return float(values[0])


def process_latex(process: str) -> str:
    boson, state = process.split("_", 1)
    state_label = state.replace("S", r"S")
    return rf"$\mathcal{{B}}({boson}\to\Upsilon({state_label})\gamma)$"


def latex_escape_text(value: str) -> str:
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
    return "".join(replacements.get(char, char) for char in value)


def format_latex_sci(value: float | None, sig_figs: int) -> str:
    if value is None:
        return r"--"
    if value == 0.0:
        return r"$0$"
    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10**exponent)
    decimals = max(sig_figs - 1, 0)
    mantissa_text = f"{mantissa:.{decimals}f}".rstrip("0").rstrip(".")
    if mantissa_text in {"10", "10."}:
        mantissa_text = "1"
        exponent += 1
    return rf"${mantissa_text}\times 10^{{{exponent}}}$"


def limit_summary_for_method(
    summary: dict[str, Any],
    scheme: str,
    method: str,
) -> dict[str, Any]:
    limits = summary.get("limits", {})
    if scheme not in limits:
        available = ", ".join(sorted(limits)) or "none"
        raise KeyError(f"Scheme {scheme!r} not found in summary. Available schemes: {available}")
    scheme_limits = limits[scheme]
    if method not in scheme_limits:
        available = ", ".join(sorted(scheme_limits)) or "none"
        raise KeyError(
            f"Method {method!r} not found for scheme {scheme!r}. Available methods: {available}"
        )
    return scheme_limits[method]


def build_process_to_poi(method_limits: dict[str, Any]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for poi, payload in method_limits.items():
        target_processes = payload.get("target_processes") or []
        target_label = payload.get("target_label")
        if not target_processes and isinstance(target_label, str):
            target_processes = [target_label]
        for process in target_processes:
            if process in THEORY_BRANCHING_FRACTIONS:
                mapping[process] = poi
    return mapping


def empty_rows_for_scheme(scheme: str) -> list[TableRow]:
    return [
        TableRow(
            process=process,
            poi="-",
            theory_bf=THEORY_BRANCHING_FRACTIONS[process],
            observed_bf_limit=None,
            expected_bf_limits={key: None for key, _ in QUANTILE_LABELS},
        )
        for process in SCHEME_PROCESS_ORDER.get(scheme, PROCESS_ORDER)
    ]


def build_rows(
    summary: dict[str, Any],
    scheme: str,
    method: str,
    warnings: list[str],
) -> list[TableRow]:
    process_order = SCHEME_PROCESS_ORDER.get(scheme, PROCESS_ORDER)
    try:
        method_limits = limit_summary_for_method(summary, scheme, method)
    except KeyError as exc:
        warnings.append(str(exc))
        return empty_rows_for_scheme(scheme)
    process_to_poi = build_process_to_poi(method_limits)
    rows: list[TableRow] = []

    for process in process_order:
        if process not in process_to_poi:
            warnings.append(f"No limit entry found for {process} in scheme {scheme}, method {method}.")
            rows.append(
                TableRow(
                    process=process,
                    poi="-",
                    theory_bf=THEORY_BRANCHING_FRACTIONS[process],
                    observed_bf_limit=None,
                    expected_bf_limits={key: None for key, _ in QUANTILE_LABELS},
                )
            )
            continue

        poi = process_to_poi[process]
        payload = method_limits[poi]
        theory = THEORY_BRANCHING_FRACTIONS[process]
        observed_strength = first_value(
            payload.get("observed"),
            warnings,
            f"{scheme}/{method}/{poi}/observed",
        )
        expected_payload = payload.get("expected", {})
        expected_bf_limits: dict[str, float | None] = {}
        for quantile_key, _ in QUANTILE_LABELS:
            value = first_value(
                expected_payload.get(quantile_key),
                warnings,
                f"{scheme}/{method}/{poi}/expected/{quantile_key}",
            )
            expected_bf_limits[quantile_key] = None if value is None else value * theory

        rows.append(
            TableRow(
                process=process,
                poi=poi,
                theory_bf=theory,
                observed_bf_limit=None
                if observed_strength is None
                else observed_strength * theory,
                expected_bf_limits=expected_bf_limits,
            )
        )
    return rows


def format_scaled(value: float, sig_figs: int) -> str:
    decimals = max(sig_figs - 1, 0)
    return f"{value:.{decimals}f}".rstrip("0").rstrip(".")


def format_latex_sci_with_delta(
    median: float | None,
    low: float | None,
    high: float | None,
    sig_figs: int,
) -> str:
    if median is None:
        return r"--"
    if low is None or high is None:
        return format_latex_sci(median, sig_figs)
    if median == 0.0:
        plus = max(high - median, 0.0)
        minus = max(median - low, 0.0)
        return rf"$0^{{+{format_scaled(plus, sig_figs)}}}_{{-{format_scaled(minus, sig_figs)}}}$"
    exponent = int(math.floor(math.log10(abs(median))))
    scale = 10**exponent
    mantissa = median / scale
    plus = max(high - median, 0.0) / scale
    minus = max(median - low, 0.0) / scale
    return (
        rf"${format_scaled(mantissa, sig_figs)}"
        rf"^{{+{format_scaled(plus, sig_figs)}}}_{{-{format_scaled(minus, sig_figs)}}}"
        rf"\times 10^{{{exponent}}}$"
    )


def format_expected_limit(row: TableRow, sig_figs: int) -> str:
    return format_latex_sci_with_delta(
        row.expected_bf_limits.get("0.5"),
        row.expected_bf_limits.get("0.16"),
        row.expected_bf_limits.get("0.84"),
        sig_figs,
    )


def make_tabular(rows: list[TableRow], sig_figs: int) -> str:
    header_cells = [
        r"Decay mode",
        r"POI",
        r"Theory",
        r"Observed",
        r"Expected",
    ]
    lines = [
        r"\begin{tabular}{llccc}",
        r"\toprule",
        " & ".join(header_cells) + r" \\",
        r"\midrule",
    ]
    for row in rows:
        cells = [
            process_latex(row.process),
            rf"\texttt{{{latex_escape_text(row.poi)}}}" if row.poi != "-" else r"--",
            format_latex_sci(row.theory_bf, sig_figs),
            format_latex_sci(row.observed_bf_limit, sig_figs),
            format_expected_limit(row, sig_figs),
        ]
        lines.append(" & ".join(cells) + r" \\")
    lines.extend([r"\bottomrule", r"\end{tabular}", ""])
    return "\n".join(lines)


def make_table_environment(section: TableSection, sig_figs: int) -> str:
    tabular = make_tabular(section.rows, sig_figs)
    return rf"""\begin{{table}}[htbp]
\centering
\caption{{{latex_escape_text(section.caption)}}}
\label{{{section.label}}}
\resizebox{{\textwidth}}{{!}}{{%
{tabular}%
}}
\end{{table}}
"""


def make_standalone_document(
    sections: list[TableSection],
    source_summary: Path,
    sig_figs: int,
) -> str:
    methods_seen: list[str] = []
    for section in sections:
        if section.method not in methods_seen:
            methods_seen.append(section.method)

    grouped: list[tuple[str, list[TableSection]]] = []
    seen: set[str] = set()
    for section in sections:
        if section.method not in seen:
            seen.add(section.method)
            grouped.append((section.method, [section]))
        else:
            for method_name, group in grouped:
                if method_name == section.method:
                    group.append(section)
                    break

    body_parts: list[str] = []
    for method_name, group in grouped:
        if len(methods_seen) > 1:
            body_parts.append(rf"\section*{{{METHOD_LABELS[method_name]}}}")
        body_parts.append(
            "\n".join(make_table_environment(section, sig_figs) for section in group)
        )

    tables_tex = "\n".join(body_parts)

    notes = [
        r"Limits are quoted on branching fractions and are obtained by multiplying the Combine signal-strength limit by the theory branching fraction in the Theory column.",
        r"Expected limits are shown as the median value with the one-standard-deviation band written as superscript/subscript deltas.",
        rf"Source summary: \texttt{{{latex_escape_text(repo_relative(source_summary))}}}.",
    ]
    schemes_seen = sorted({section.scheme for section in sections})
    if any(scheme in {"h_grouped", "z_grouped"} for scheme in schemes_seen):
        notes.append(
            r"For a grouped boson POI scheme, one common signal-strength POI scales all signal processes for the selected boson while the other boson signal strengths are profiled individually."
        )

    method_labels = ", ".join(METHOD_LABELS[m] for m in methods_seen)
    notes_text = "\n".join(rf"\item {note}" for note in notes)
    return rf"""\documentclass[11pt]{{article}}
\usepackage[margin=0.7in]{{geometry}}
\usepackage{{booktabs}}
\usepackage{{caption}}
\usepackage{{graphicx}}
\usepackage{{amsmath}}
\usepackage{{array}}
\begin{{document}}
{tables_tex}

\noindent\textbf{{Configuration:}} schemes \texttt{{{latex_escape_text(', '.join(schemes_seen))}}}, methods: {method_labels}.
\begin{{itemize}}
{notes_text}
\end{{itemize}}
\end{{document}}
"""


def write_outputs(
    sections: list[TableSection],
    summary_path: Path,
    output_dir: Path,
    output_name: str,
    sig_figs: int,
    warnings: list[str],
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    fragment_path = output_dir / f"{output_name}.table.tex"
    document_path = output_dir / f"{output_name}.tex"
    metadata_path = output_dir / f"{output_name}.json"

    all_tables_tex = "\n".join(
        make_table_environment(section, sig_figs) for section in sections
    )
    fragment_path.write_text(all_tables_tex, encoding="ascii")
    document_path.write_text(
        make_standalone_document(
            sections=sections,
            source_summary=summary_path,
            sig_figs=sig_figs,
        ),
        encoding="ascii",
    )

    methods_seen = sorted({section.method for section in sections})
    metadata = {
        "schema_version": 1,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_summary": str(summary_path.resolve()),
        "source_summary_repo_relative": repo_relative(summary_path),
        "methods": methods_seen,
        "schemes": sorted({section.scheme for section in sections}),
        "theory_branching_fractions": THEORY_BRANCHING_FRACTIONS,
        "tables": [
            {
                "method": section.method,
                "scheme": section.scheme,
                "title": section.title,
                "caption": section.caption,
                "label": section.label,
                "rows": [
                    {
                        "process": row.process,
                        "poi": row.poi,
                        "theory_bf": row.theory_bf,
                        "observed_bf_limit": row.observed_bf_limit,
                        "expected_bf_limits": row.expected_bf_limits,
                    }
                    for row in section.rows
                ],
            }
            for section in sections
        ],
        "warnings": warnings,
        "outputs": {
            "fragment_tex": str(fragment_path.resolve()),
            "standalone_tex": str(document_path.resolve()),
        },
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "fragment": fragment_path,
        "document": document_path,
        "metadata": metadata_path,
    }


def compile_with_singularity(
    document_path: Path,
    image: str,
    singularity_bin: str,
) -> Path:
    output_dir = document_path.parent.resolve()
    command = [
        singularity_bin,
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
    completed = subprocess.run(command, capture_output=True, text=True)
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
                "engine": singularity_bin,
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
        description=(
            "Create a LaTeX table of theory branching fractions and observed/expected branching-fraction limits "
            "from scripts/limits.py JSON output."
        )
    )
    parser.add_argument(
        "--summary",
        default=None,
        help="Path to blind_limits_summary.json. Defaults to the newest summary under datacards/blind_limits.",
    )
    parser.add_argument(
        "--scheme",
        choices=("all", "six_poi", "z_grouped", "h_grouped"),
        default="all",
        help="POI scheme to tabulate. Default writes grouped H, grouped Z, and six-POI tables.",
    )
    parser.add_argument(
        "--method",
        choices=("hybrid_lhc", "asymptotic", "both"),
        default="both",
        help="Limit method to tabulate. 'both' writes asymptotic and HybridNew tables in one document.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory. Defaults to a tables/ directory next to the summary JSON.",
    )
    parser.add_argument(
        "--output-name",
        default=DEFAULT_OUTPUT_NAME,
        help="Base filename for .tex, .table.tex, .json, and optional .pdf outputs.",
    )
    parser.add_argument(
        "--sig-figs",
        type=int,
        default=DEFAULT_SIG_FIGS,
        help="Significant figures used in scientific-notation table cells.",
    )
    parser.add_argument(
        "--compile-pdf",
        action="store_true",
        help="Compile the standalone .tex to PDF using Singularity with a Docker LaTeX image.",
    )
    parser.add_argument(
        "--latex-image",
        default=DEFAULT_IMAGE,
        help="Docker image URI passed to Singularity, e.g. docker://...",
    )
    parser.add_argument(
        "--singularity-bin",
        default="singularity",
        help="Singularity executable name or path.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    summary_path = Path(args.summary).resolve() if args.summary else latest_summary_path().resolve()
    ensure_file(summary_path, "Blind-limit summary JSON not found")
    output_dir = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else summary_path.parent / "tables"
    )
    if args.sig_figs < 1:
        raise ValueError("--sig-figs must be at least 1")
    if args.compile_pdf and not args.latex_image.startswith("docker://"):
        raise ValueError("--latex-image must be a docker:// URI when compiling through Singularity")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    warnings: list[str] = []
    schemes = list(DEFAULT_TABLE_SCHEMES if args.scheme == "all" else (args.scheme,))
    methods = list(DEFAULT_METHODS if args.method == "both" else (args.method,))
    sections: list[TableSection] = []
    for method in methods:
        for scheme in schemes:
            config = SECTION_CONFIG[scheme]
            sections.append(
                TableSection(
                    scheme=scheme,
                    method=method,
                    title=config["title"],
                    caption=config["caption"],
                    label=f"{config['label_base']}_{method}",
                    rows=build_rows(summary, scheme, method, warnings),
                )
            )
    paths = write_outputs(
        sections=sections,
        summary_path=summary_path,
        output_dir=output_dir,
        output_name=args.output_name,
        sig_figs=args.sig_figs,
        warnings=warnings,
    )

    log(f"Wrote table fragment: {repo_relative(paths['fragment'])}")
    log(f"Wrote standalone LaTeX: {repo_relative(paths['document'])}")
    log(f"Wrote metadata JSON: {repo_relative(paths['metadata'])}")
    for warning in warnings:
        log(f"warning: {warning}")

    if args.compile_pdf:
        pdf_path = compile_with_singularity(
            paths["document"],
            image=args.latex_image,
            singularity_bin=args.singularity_bin,
        )
        log(f"Wrote PDF: {repo_relative(pdf_path)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
