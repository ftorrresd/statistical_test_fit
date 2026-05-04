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
            f"No blind limit summary found under {DEFAULT_LIMITS_DIR}. Run scripts/blind_limits.py first or pass --summary."
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


def build_rows(
    summary: dict[str, Any],
    scheme: str,
    method: str,
    warnings: list[str],
) -> list[TableRow]:
    method_limits = limit_summary_for_method(summary, scheme, method)
    process_to_poi = build_process_to_poi(method_limits)
    rows: list[TableRow] = []

    for process in PROCESS_ORDER:
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


def make_tabular(rows: list[TableRow], sig_figs: int) -> str:
    header_cells = [
        r"Decay mode",
        r"POI",
        r"Theory",
        r"Observed",
        *[label for _, label in QUANTILE_LABELS],
    ]
    lines = [
        r"\begin{tabular}{llccccc}",
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
        ]
        cells.extend(
            format_latex_sci(row.expected_bf_limits.get(quantile_key), sig_figs)
            for quantile_key, _ in QUANTILE_LABELS
        )
        lines.append(" & ".join(cells) + r" \\")
    lines.extend([r"\bottomrule", r"\end{tabular}", ""])
    return "\n".join(lines)


def make_standalone_document(
    tabular: str,
    caption: str,
    label: str,
    source_summary: Path,
    scheme: str,
    method: str,
    grouped_note: bool,
) -> str:
    notes = [
        r"Limits are quoted on branching fractions and are obtained by multiplying the Combine signal-strength limit by the theory branching fraction in the Theory column.",
        rf"Source summary: \texttt{{{latex_escape_text(repo_relative(source_summary))}}}.",
    ]
    if grouped_note:
        notes.append(
            r"For the grouped Upsilon POI scheme, one common signal-strength POI scales the H and Z processes with the same Upsilon state; the H and Z branching-fraction limits in each state are therefore common-strength translations, not independent one-process limits."
        )

    notes_text = "\n".join(rf"\item {note}" for note in notes)
    return rf"""\documentclass[11pt]{{article}}
\usepackage[margin=0.7in]{{geometry}}
\usepackage{{booktabs}}
\usepackage{{caption}}
\usepackage{{graphicx}}
\usepackage{{amsmath}}
\usepackage{{array}}
\begin{{document}}
\begin{{table}}[htbp]
\centering
\caption{{{caption}}}
\label{{{label}}}
\resizebox{{\textwidth}}{{!}}{{%
{tabular}%
}}
\end{{table}}

\noindent\textbf{{Configuration:}} scheme \texttt{{{latex_escape_text(scheme)}}}, method \texttt{{{latex_escape_text(method)}}}.
\begin{{itemize}}
{notes_text}
\end{{itemize}}
\end{{document}}
"""


def write_outputs(
    rows: list[TableRow],
    summary_path: Path,
    output_dir: Path,
    output_name: str,
    caption: str,
    label: str,
    scheme: str,
    method: str,
    sig_figs: int,
    warnings: list[str],
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    tabular = make_tabular(rows, sig_figs)
    fragment_path = output_dir / f"{output_name}.table.tex"
    document_path = output_dir / f"{output_name}.tex"
    metadata_path = output_dir / f"{output_name}.json"
    fragment_path.write_text(tabular, encoding="ascii")
    document_path.write_text(
        make_standalone_document(
            tabular=tabular,
            caption=caption,
            label=label,
            source_summary=summary_path,
            scheme=scheme,
            method=method,
            grouped_note=scheme == "grouped_upsilon_poi",
        ),
        encoding="ascii",
    )

    metadata = {
        "schema_version": 1,
        "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_summary": str(summary_path.resolve()),
        "source_summary_repo_relative": repo_relative(summary_path),
        "scheme": scheme,
        "method": method,
        "theory_branching_fractions": THEORY_BRANCHING_FRACTIONS,
        "rows": [
            {
                "process": row.process,
                "poi": row.poi,
                "theory_bf": row.theory_bf,
                "observed_bf_limit": row.observed_bf_limit,
                "expected_bf_limits": row.expected_bf_limits,
            }
            for row in rows
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
            "from scripts/blind_limits.py JSON output."
        )
    )
    parser.add_argument(
        "--summary",
        default=None,
        help="Path to blind_limits_summary.json. Defaults to the newest summary under datacards/blind_limits.",
    )
    parser.add_argument(
        "--scheme",
        choices=("six_poi", "grouped_upsilon_poi"),
        default="six_poi",
        help="POI scheme to tabulate.",
    )
    parser.add_argument(
        "--method",
        choices=("hybrid_lhc", "asymptotic"),
        default="hybrid_lhc",
        help="Limit method to tabulate.",
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
        "--caption",
        default="Theory branching fractions and expected upper limits on branching fractions.",
        help="LaTeX table caption.",
    )
    parser.add_argument(
        "--label",
        default="tab:branching_fraction_limits",
        help="LaTeX label for the standalone table.",
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
    rows = build_rows(summary, args.scheme, args.method, warnings)
    paths = write_outputs(
        rows=rows,
        summary_path=summary_path,
        output_dir=output_dir,
        output_name=args.output_name,
        caption=args.caption,
        label=args.label,
        scheme=args.scheme,
        method=args.method,
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
