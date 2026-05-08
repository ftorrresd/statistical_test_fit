from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class DisplayName:
    slug: str
    text: str
    latex: str

    def metadata(self, prefix: str = "display") -> dict[str, str]:
        return {
            f"{prefix}_slug": self.slug,
            f"{prefix}_label": self.text,
            f"{prefix}_latex": self.latex,
        }


def safe_slug(value: Any) -> str:
    text = str(value or "").strip().lower()
    text = text.replace("->", " to ")
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = text.strip("_")
    return text or "unnamed"


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


def normalize_pdf_family(value: Any) -> str:
    family = safe_slug(value)
    aliases = {
        "cheb": "chebychev",
        "chebyshev": "chebychev",
        "cheby": "chebychev",
        "johnson_su": "johnson",
        "power": "power_law",
        "powerlaw": "power_law",
        "power_law_pdf": "power_law",
        "exp": "exponential",
        "expo": "exponential",
    }
    return aliases.get(family, family)


def infer_pdf_family_from_name(*names: Any) -> str:
    lower = " ".join(str(name or "").lower() for name in names)
    if "johnson" in lower:
        return "johnson"
    if "bernstein" in lower:
        return "bernstein"
    if "chebychev" in lower or "chebyshev" in lower or "cheb" in lower:
        return "chebychev"
    if "power_law" in lower or "powerlaw" in lower or "power-law" in lower:
        return "power_law"
    if "exponential" in lower or "rooexppoly" in lower:
        return "exponential"
    return "other"


def _int_or_none(value: Any) -> int | None:
    if value is None or value == "":
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not numeric.is_integer():
        return None
    return int(numeric)


def infer_pdf_order(*values: Any) -> int | None:
    for value in values:
        direct = _int_or_none(value)
        if direct is not None:
            return direct

    text = " ".join(str(value or "") for value in values).lower()
    patterns = (
        r"(?:order|degree|deg)[_-]?(\d+)",
        r"(?:bernstein|chebychev|chebyshev|power_law|powerlaw|exponential)(?:_2d)?[_-](\d+)(?:[_\W]|$)",
        r"background_(?:bernstein|chebychev|chebyshev|power_law|powerlaw|exponential)_2d_(\d+)",
    )
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return int(match.group(1))
    return None


def _role_display(selection_role: Any) -> tuple[str, str] | None:
    role = safe_slug(selection_role)
    if role in {"", "unnamed", "nominal"}:
        return None
    return {
        "winner": ("selected_model", "selected model"),
        "selected": ("selected_model", "selected model"),
        "below": ("below_selected_order", "below selected order"),
        "above": ("above_selected_order", "above selected order"),
    }.get(role, (role, str(selection_role).replace("_", " ")))


def _family_base_display(family: str, order: int | None) -> tuple[str, str, str]:
    if family == "johnson":
        return "johnson_su_nominal", "Johnson SU nominal", "Johnson SU nominal"

    family_templates = {
        "bernstein": ("bernstein_polynomial", "Bernstein polynomial"),
        "chebychev": ("chebyshev_polynomial", "Chebyshev polynomial"),
        "power_law": ("power_law_model", "Power-law model"),
        "exponential": ("exponential_model", "Exponential model"),
    }
    slug_base, label_base = family_templates.get(
        family,
        (safe_slug(family), str(family).replace("_", " ").title()),
    )
    if order is None:
        return slug_base, label_base, latex_escape_text(label_base)
    return (
        f"{slug_base}_order{order}",
        f"{label_base}, order {order}",
        f"{latex_escape_text(label_base)}, order ${order}$",
    )


def pdf_state_display(
    *,
    index: Any = None,
    family: Any = None,
    order: Any = None,
    selection_role: Any = None,
    name: Any = None,
    source_model_name: Any = None,
    workspace_label: Any = None,
) -> DisplayName:
    family_name = normalize_pdf_family(family) if family not in (None, "") else ""
    if not family_name:
        family_name = infer_pdf_family_from_name(name, source_model_name, workspace_label)
    order_value = infer_pdf_order(order, name, source_model_name, workspace_label)
    index_value = _int_or_none(index)

    slug, text, latex = _family_base_display(family_name, order_value)
    role = _role_display(selection_role)
    qualifiers: list[str] = []
    latex_qualifiers: list[str] = []
    if role is not None:
        role_slug, role_text = role
        slug = f"{slug}_{role_slug}"
        qualifiers.append(role_text)
        latex_qualifiers.append(latex_escape_text(role_text))
    if index_value is not None:
        slug = f"pdfindex{index_value}_{slug}"
        qualifiers.append(f"pdfindex {index_value}")
        latex_qualifiers.append(rf"$\mathrm{{pdfindex}}={index_value}$")

    if qualifiers:
        text = f"{text} (" + ", ".join(qualifiers) + ")"
        latex = f"{latex} (" + ", ".join(latex_qualifiers) + ")"
    return DisplayName(slug=safe_slug(slug), text=text, latex=latex)


def floating_pdf_display() -> DisplayName:
    return DisplayName(
        slug="floating_pdfindex",
        text="floating pdfindex",
        latex=r"floating $\mathrm{pdfindex}$",
    )


def poi_scheme_display(name: Any) -> DisplayName:
    raw = str(name or "")
    normalized = safe_slug(raw)
    aliases = {
        "grouped_z": "z_grouped",
        "grouped_h": "h_grouped",
    }
    normalized = aliases.get(normalized, normalized)
    mapping = {
        "three_poi_z": DisplayName(
            "z_separate_upsilon_state_pois",
            "Z: separate Upsilon(1S,2S,3S) POIs",
            r"$Z$: separate $\Upsilon(1S,2S,3S)$ POIs",
        ),
        "three_poi_h": DisplayName(
            "h_separate_upsilon_state_pois",
            "H: separate Upsilon(1S,2S,3S) POIs",
            r"$H$: separate $\Upsilon(1S,2S,3S)$ POIs",
        ),
        "z_grouped": DisplayName(
            "z_combined_upsilon_state_poi",
            "Z: combined Upsilon(nS) POI",
            r"$Z$: combined $\Upsilon(nS)$ POI",
        ),
        "h_grouped": DisplayName(
            "h_combined_upsilon_state_poi",
            "H: combined Upsilon(nS) POI",
            r"$H$: combined $\Upsilon(nS)$ POI",
        ),
    }
    return mapping.get(
        normalized,
        DisplayName(safe_slug(raw), raw or "unnamed scheme", latex_escape_text(raw or "unnamed scheme")),
    )


def signal_target_display(value: Any) -> DisplayName:
    raw = str(value or "").strip()
    token = raw[2:] if raw.startswith("r_") else raw
    parts = token.split("_", 1)
    if len(parts) == 2 and parts[0] in {"H", "Z"}:
        boson, state = parts
        boson_slug = boson.lower()
        if state.lower() == "grouped":
            return DisplayName(
                f"{boson_slug}_combined_upsilon",
                f"{boson} -> Upsilon(nS) gamma combined",
                rf"${boson} \to \Upsilon(nS)\gamma$ combined",
            )
        state_latex = state
        return DisplayName(
            f"{boson_slug}_upsilon_{state.lower()}",
            f"{boson} -> Upsilon({state}) gamma",
            rf"${boson} \to \Upsilon({state_latex})\gamma$",
        )
    return DisplayName(safe_slug(raw), raw or "unnamed target", latex_escape_text(raw or "unnamed target"))


def dataset_strategy_display(value: Any) -> DisplayName:
    strategy = safe_slug(value)
    if strategy == "toys":
        return DisplayName("toy_datasets", "toy datasets", "toy datasets")
    if strategy == "asimov":
        return DisplayName("asimov_dataset", "Asimov dataset", "Asimov dataset")
    return DisplayName(strategy, str(value), latex_escape_text(value))
