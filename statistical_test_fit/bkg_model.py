import math
from pprint import pprint
from dataclasses import dataclass
from typing import Any, Optional

from ROOT import (  # type: ignore
    RooArgList,  # type: ignore
    RooAddPdf,  # type: ignore
    RooArgSet,  # type: ignore
    RooBernstein,  # type: ignore
    RooChebychev,  # type: ignore
    RooJohnson,  # type: ignore
    kTRUE,  # type: ignore
    RooGenericPdf,  # type: ignore
    RooProdPdf,  # type: ignore
    RooRealVar,  # type: ignore
    TMath,  # type: ignore
)

try:
    from ROOT import RooExpPoly  # type: ignore
except ImportError:
    RooExpPoly = None

from .bkg_pdf_families import BkgPdfFamily
from .chi2_test import ChiSquareResult
from .dimuon_non_correlated import build_upsilon_model
from .ws_helper import set_pdf_parameters, get_pdf_parameters
from .resonant_bkg_modeling import build_resonant_background_modeling_Z


@dataclass
class BkgModel:
    pdf_family: BkgPdfFamily
    model: Any
    scan_order: Optional[int] = None
    n_params: Optional[int] = None
    n_float_params: Optional[int] = None
    chi_square_res: Optional[ChiSquareResult] = None
    NLL: Optional[float] = None
    fit_res: Optional[Any] = None
    fit_status: Optional[int] = None
    fit_cov_qual: Optional[int] = None
    fit_edm: Optional[float] = None
    fit_ok: Optional[bool] = None
    fit_quality_reason: Optional[str] = None

    def __post_init__(self):
        self.n_params = self.model.getParameters(RooArgSet()).getSize()

    def is_complete(self):
        return (
            (self.chi_square_res is not None)
            and (self.NLL is not None)
            and (self.n_params is not None)
            and (self.n_float_params is not None)
            and (self.fit_ok is not None)
            and (
                (self.fit_res is not None)
                or ((self.fit_status is not None) and (self.fit_cov_qual is not None))
            )
        )

    def __hash__(self) -> int:
        return hash(
            (
                self.pdf_family,
                self.scan_order,
                self.chi_square_res,
                self.n_params,
                self.n_float_params,
                self.NLL,
            )
        )

    def __str__(self) -> str:
        self.model.Print()
        string_repr = []
        string_repr.append(f"Polinomial family: {self.pdf_family}")
        if self.scan_order is not None:
            string_repr.append(f"Scanned order: {self.scan_order}")
        string_repr.append(f"N params: {self.n_params}")
        if self.n_float_params is not None:
            string_repr.append(f"N floated params: {self.n_float_params}")
        if self.fit_ok is not None:
            fit_quality = "OK" if self.fit_ok else "FAILED"
            string_repr.append(f"Fit quality: {fit_quality}")
        if self.fit_quality_reason is not None:
            string_repr.append(f"Fit note: {self.fit_quality_reason}")
        if self.fit_res is not None:
            string_repr.append(
                f"Status: {self.fit_res.status()}  (0=OK)   CovQual: {self.fit_res.covQual()}  (3=good)"
            )
            pars = self.fit_res.floatParsFinal()
            for i in range(pars.getSize()):
                p = pars[i]
                string_repr.append(
                    f" {p.GetName():>12s} = {p.getVal():9.4f} ± {p.getError():7.4f}"
                )

            string_repr.append(f"Chi2 test results: {self.chi_square_res}")
            string_repr.append(f"NLL: {self.NLL}")
        elif self.fit_status is not None and self.fit_cov_qual is not None:
            string_repr.append(
                f"Status: {self.fit_status}  (0=OK)   CovQual: {self.fit_cov_qual}  (3=good)"
            )
            if self.fit_edm is not None:
                string_repr.append(f" EDM: {self.fit_edm:.3e}")
            if self.chi_square_res is not None:
                string_repr.append(f"Chi2 test results: {self.chi_square_res}")
            if self.NLL is not None:
                string_repr.append(f"NLL: {self.NLL}")

        return "\n".join(string_repr)


def _print_big_warning(title: str, details: Optional[list[str]] = None) -> None:
    banner = "!" * 90
    print()
    print(banner)
    print(f"!!! BIG WARNING: {title}")
    if details is not None:
        for detail in details:
            print(f"!!! {detail}")
    print(banner)
    print()


def _assess_fit_result(
    fit_res,
    require_covqual: int = 3,
    max_edm: float = 1e-3,
) -> tuple[bool, str]:
    reasons = []

    status = fit_res.status()
    if status != 0:
        reasons.append(f"status={status} (expected 0)")

    cov_qual = fit_res.covQual()
    if cov_qual < require_covqual:
        reasons.append(f"covQual={cov_qual} (expected >= {require_covqual})")

    if hasattr(fit_res, "edm"):
        edm = fit_res.edm()
        if not math.isfinite(edm):
            reasons.append("edm is not finite")
        elif edm > max_edm:
            reasons.append(f"edm={edm:.3e} (expected <= {max_edm:.3e})")

    if hasattr(fit_res, "minNll") and not math.isfinite(fit_res.minNll()):
        reasons.append("minNll is not finite")

    if reasons:
        return False, "; ".join(reasons)

    return True, f"status=0, covQual>={require_covqual}, edm<={max_edm:.3e}"


def update_fit_quality(test_bkg_model: BkgModel) -> None:
    if test_bkg_model.fit_res is None:
        raise RuntimeError(
            f"Cannot assess fit quality for {test_bkg_model.model.GetName()} without a fit result."
        )

    test_bkg_model.n_float_params = test_bkg_model.fit_res.floatParsFinal().getSize()
    test_bkg_model.fit_status = test_bkg_model.fit_res.status()
    test_bkg_model.fit_cov_qual = test_bkg_model.fit_res.covQual()
    if hasattr(test_bkg_model.fit_res, "edm"):
        test_bkg_model.fit_edm = test_bkg_model.fit_res.edm()

    fit_ok, fit_quality_reason = _assess_fit_result(test_bkg_model.fit_res)
    if test_bkg_model.NLL is None or not math.isfinite(test_bkg_model.NLL):
        fit_ok = False
        if test_bkg_model.NLL is None:
            fit_quality_reason = f"{fit_quality_reason}; NLL was not stored"
        else:
            fit_quality_reason = (
                f"{fit_quality_reason}; NLL={test_bkg_model.NLL} is not finite"
            )

    test_bkg_model.fit_ok = fit_ok
    test_bkg_model.fit_quality_reason = fit_quality_reason

    if not test_bkg_model.fit_ok:
        details = [
            f"Family: {test_bkg_model.pdf_family}",
            f"Model: {test_bkg_model.model.GetName()}",
            f"Status: {test_bkg_model.fit_res.status()}",
            f"CovQual: {test_bkg_model.fit_res.covQual()}",
            f"N floated params: {test_bkg_model.n_float_params}",
        ]
        if hasattr(test_bkg_model.fit_res, "edm"):
            details.append(f"EDM: {test_bkg_model.fit_res.edm():.3e}")
        if test_bkg_model.NLL is not None:
            details.append(f"NLL: {test_bkg_model.NLL:.6f}")
        details.append(f"Reason: {test_bkg_model.fit_quality_reason}")
        _print_big_warning(
            f"Fit failed quality checks for {test_bkg_model.model.GetName()}",
            details,
        )


def _format_candidate_summary(idx: int, test_bkg_model: BkgModel) -> str:
    summary = [f"[{idx}] {test_bkg_model.model.GetName()}"]

    if test_bkg_model.n_params is not None:
        summary.append(f"n_params={test_bkg_model.n_params}")
    if test_bkg_model.scan_order is not None:
        summary.append(f"order={test_bkg_model.scan_order}")
    if test_bkg_model.n_float_params is not None:
        summary.append(f"n_float={test_bkg_model.n_float_params}")
    if test_bkg_model.fit_ok is True:
        summary.append("fit=OK")
    elif test_bkg_model.fit_ok is False:
        summary.append(f"fit=FAILED ({test_bkg_model.fit_quality_reason})")

    if test_bkg_model.chi_square_res is not None:
        summary.append(f"chi2_p={test_bkg_model.chi_square_res.pvalue:.4f}")
    if test_bkg_model.NLL is not None:
        summary.append(f"NLL={test_bkg_model.NLL:.6f}")

    return " | ".join(summary)


def _raise_family_strict_selection_error(
    test_bkg_models: list[BkgModel],
    message: str,
) -> None:
    family = "unknown"
    if len(test_bkg_models) > 0:
        family = str(test_bkg_models[0].pdf_family)

    details = [
        f"Family-strict model selection aborted for {family}.",
        message,
    ]
    details.extend(
        _format_candidate_summary(idx, test_bkg_model)
        for idx, test_bkg_model in enumerate(test_bkg_models)
    )
    _print_big_warning(f"Family-strict selection failed for {family}", details)
    raise RuntimeError(message)


def _select_relaxed_winner(
    test_bkg_models: list[BkgModel],
    failure_message: str,
) -> tuple[int, int]:
    family = str(test_bkg_models[0].pdf_family)
    valid_models = [
        (idx, test_bkg_model)
        for idx, test_bkg_model in enumerate(test_bkg_models)
        if test_bkg_model.fit_ok
    ]

    details = [
        f"Relaxed selection enabled for {family}.",
        failure_message,
    ]

    if len(valid_models) == 0:
        _print_big_warning(
            f"No fit-quality-passing candidates available for {family}",
            details
            + [
                _format_candidate_summary(idx, test_bkg_model)
                for idx, test_bkg_model in enumerate(test_bkg_models)
            ],
        )
        fallback_idx = len(test_bkg_models) - 1
        return fallback_idx, fallback_idx

    start_idx = None
    winner_idx = None
    best_pvalue = -1.0
    for valid_idx, (_, test_bkg_model) in enumerate(valid_models):
        assert test_bkg_model.chi_square_res is not None
        pvalue = test_bkg_model.chi_square_res.pvalue
        if pvalue > 0.05:
            start_idx = valid_idx
            winner_idx = valid_idx
            break
        if pvalue > best_pvalue:
            best_pvalue = pvalue
            start_idx = valid_idx
            winner_idx = valid_idx

    assert start_idx is not None and winner_idx is not None

    details.append(
        "Falling back to the fit-quality-passing candidate with the best available chi-square p-value."
    )
    details.extend(
        _format_candidate_summary(idx, test_bkg_model)
        for idx, test_bkg_model in enumerate(test_bkg_models)
    )
    _print_big_warning(f"Relaxed model selection fallback for {family}", details)

    return valid_models[start_idx][0], valid_models[winner_idx][0]


def _compute_lrt_pvalue(
    simple_model: BkgModel,
    complex_model: BkgModel,
    tol: float = 1e-6,
) -> tuple[float, int, float]:
    if simple_model.NLL is None or complex_model.NLL is None:
        raise ValueError("Both models must have a finite stored NLL before the LRT.")

    if simple_model.n_float_params is None or complex_model.n_float_params is None:
        raise ValueError("Both models must store the number of floated parameters.")

    delta_n_float_params = complex_model.n_float_params - simple_model.n_float_params
    if delta_n_float_params <= 0:
        raise ValueError(
            "The complex model must have more floated parameters than the simple one. "
            f"Received df={delta_n_float_params}."
        )

    q = 2.0 * (simple_model.NLL - complex_model.NLL)
    if not math.isfinite(q):
        raise ValueError(f"2*DeltaNLL is not finite: {q}")
    if q < -tol:
        raise ValueError(
            f"Encountered a negative 2*DeltaNLL beyond tolerance. Received {q:.6f}."
        )

    q = max(q, 0.0)
    pval = TMath.Prob(q, delta_n_float_params)
    return q, delta_n_float_params, pval


def compute_lrt_summary(
    simple_model: BkgModel,
    complex_model: BkgModel,
    tol: float = 1e-6,
) -> dict[str, Any]:
    q, delta_n_float_params, pval = _compute_lrt_pvalue(
        simple_model,
        complex_model,
        tol=tol,
    )
    assert simple_model.NLL is not None and complex_model.NLL is not None
    return {
        "delta_nll": simple_model.NLL - complex_model.NLL,
        "two_delta_nll": q,
        "delta_n_float_params": delta_n_float_params,
        "pvalue": pval,
    }


def compute_winner_and_start_indexes(
    test_bkg_models: list[BkgModel],
    strict_mode: bool = True,
) -> tuple[int, int]:
    if len(test_bkg_models) == 0:
        raise RuntimeError("Cannot select a winner from an empty family.")

    valid_models = []
    for idx, test_bkg_model in enumerate(test_bkg_models):
        if not test_bkg_model.is_complete():
            message = f"Model {idx} is missing fit diagnostics required by selection."
            if strict_mode:
                _raise_family_strict_selection_error(test_bkg_models, message)
            return _select_relaxed_winner(test_bkg_models, message)

        if test_bkg_model.fit_ok:
            valid_models.append((idx, test_bkg_model))

    if len(valid_models) == 0:
        message = "No candidate in this family passed the fit-quality requirements."
        if strict_mode:
            _raise_family_strict_selection_error(test_bkg_models, message)
        return _select_relaxed_winner(test_bkg_models, message)

    start_idx = None
    winner_idx = None
    for valid_idx, (_, test_bkg_model) in enumerate(valid_models):
        assert test_bkg_model.chi_square_res is not None
        if test_bkg_model.chi_square_res.pvalue > 0.05:
            start_idx = valid_idx
            winner_idx = valid_idx
            break

    if start_idx is None or winner_idx is None:
        message = "No candidate in this family passed the chi-square p-value threshold after fit-quality filtering."
        if strict_mode:
            _raise_family_strict_selection_error(test_bkg_models, message)
        return _select_relaxed_winner(test_bkg_models, message)

    for valid_idx in range(start_idx + 1, len(valid_models)):
        _, previous_model = valid_models[valid_idx - 1]
        _, current_model = valid_models[valid_idx]

        try:
            _, _, pval = _compute_lrt_pvalue(previous_model, current_model)
        except ValueError as exc:
            message = (
                "Invalid likelihood-ratio comparison between "
                f"{previous_model.model.GetName()} and {current_model.model.GetName()}: {exc}"
            )
            if strict_mode:
                _raise_family_strict_selection_error(test_bkg_models, message)
            return _select_relaxed_winner(test_bkg_models, message)

        if pval > 0.05:
            winner_idx = valid_idx - 1
            break

        winner_idx = valid_idx

    return valid_models[start_idx][0], valid_models[winner_idx][0]


def _get_power_law_formula(exponents):
    # build formula: @0 is x, then a0,b0,a1,b1,...
    # term i uses indices (1+2*i) for a_i and (2+2*i) for b_i
    formula = ""
    if len(exponents) == 1:
        formula = "pow(x, -a0)"
    elif len(exponents) == 2:
        formula = "pow(x, -a0) * pow(1 - x, -a1)"
    elif len(exponents) == 3:
        formula = "pow(x, -(a0 + a2*log(x))) * pow(1 - x, -a1)"
    elif len(exponents) == 4:
        formula = "pow(x, -(a0 + a2*log(x) + a3*log(x))) * pow(1 - x, -a1)"
    elif len(exponents) == 5:
        formula = (
            "pow(x, -(a0 + a2*log(x) + a3*log(x)+a4*pow(log(x),2))) * pow(1 - x, -a1)"
        )
    elif len(exponents) == 6:
        formula = "pow(x, -(a0 + a2*log(x) + a3*log(x)+a3*pow(log(x),2)+a4*pow(log(x),3) )) * pow(1 - x, -a1)"
    elif len(exponents) == 7:
        formula = "pow(x, -(a0 + a2*log(x) + a3*log(x)+a3*pow(log(x),2)+a4*pow(log(x),3)+a5*pow(log(x),4) )) * pow(1 - x, -a1)"
    elif len(exponents) == 8:
        formula = "pow(x, -(a0 + a2*log(x) + a3*log(x)+a3*pow(log(x),2)+a4*pow(log(x),3)+a5*pow(log(x),4)+a6*pow(log(x),5) )) * pow(1 - x, -a1)"
    elif len(exponents) == 9:
        formula = "pow(x, -(a0 + a2*log(x) + a3*log(x)+a3*pow(log(x),2)+a4*pow(log(x),3)+a5*pow(log(x),4)+a6*pow(log(x),5)+a7*pow(log(x),6)+a8*pow(log(x),7) )) * pow(1 - x, -a1)"
    else:
        raise ValueError("Number of exponents > 8")

    formula = formula.replace("x", "@0")
    for i in range(len(exponents)):
        formula = formula.replace(f"a{i}", f"@{i + 1}")

    return formula


def build_background_cheb(x, cheb_coeffs, name="background_chebychev"):
    coeff_vars = [
        RooRealVar(f"c{i}", f"Chebychev c{i}", float(v), -1.0, 1.0)
        for i, v in enumerate(cheb_coeffs, start=1)
    ]

    coeff_list = RooArgList()
    for v in coeff_vars:
        coeff_list.add(v)

    bkg = RooChebychev(name, "Chebychev background", x, coeff_list)
    bkg._keepalive = {"coeff_vars": coeff_vars, "coeff_list": coeff_list}

    return bkg


def build_background_bernstein(
    x, bern_coeffs, name="background_bernstein", force_positive=True
):
    """
    Build a RooBernstein PDF:
      - len(bern_coeffs) = n+1 gives a degree-n Bernstein polynomial.
      - By default, coefficients are constrained to be >= 0 for positivity.
    """
    lo = 0.0 if force_positive else -10.0
    hi = 10.0

    coeff_vars = [
        RooRealVar(f"b{i}", f"Bernstein b{i}", float(v), lo, hi)
        for i, v in enumerate(bern_coeffs)
    ]

    coeff_list = RooArgList()
    for v in coeff_vars:
        coeff_list.add(v)

    bkg = RooBernstein(name, "Bernstein background", x, coeff_list)
    bkg._keepalive = {"coeff_vars": coeff_vars, "coeff_list": coeff_list}
    return bkg


def build_background_power_law(
    x,
    exponents,
    name="background_power_law",
    title="Power-sum background",
    exponents_bounds=(-10.0, +10.0),
):
    if not (len(exponents) >= 1):
        raise ValueError("At least one exponent should be provided")

    # parameters
    a_vars = []
    for i, v in enumerate(exponents):
        a_vars.append(RooRealVar(f"a{i}", f"a{i}", float(v), *exponents_bounds))

    formula = _get_power_law_formula(exponents)

    args = RooArgList(x)
    for i in range(len(exponents)):
        args.add(a_vars[i])

    pdf = RooGenericPdf(name, title, formula, args)
    pdf._keepalive = {"a_vars": a_vars, "args": args}

    return pdf


def build_background_exponential(
    x,
    coeffs,
    name="background_exponential",
    title="Exponential background",
    coeffs_bounds=(-1000000000000.0, 0.0),
):
    if RooExpPoly is None:
        raise ImportError("RooExpPoly is not available in this ROOT build.")

    if not (len(coeffs) >= 1):
        raise ValueError("At least one coeffs should be provided")

    # parameters
    a_vars = []
    for i, v in enumerate(coeffs):
        a_vars.append(RooRealVar(f"a{i}", f"a{i}", float(v), *coeffs_bounds))

    args = RooArgList(x)
    for i in range(len(coeffs)):
        args.add(a_vars[i])

    pdf = RooExpPoly(name, title, x, args)
    pdf._keepalive = {"a_vars": a_vars, "args": args}

    return pdf


def _build_x_background_component(
    x,
    name: str,
    coeff_vars_x,
    coeff_list_x,
    upsilon_params=None,
):
    bkg_x_lin = RooChebychev(
        f"{name}_chebychev_x", "Chebychev background 2D - x", x, coeff_list_x
    )
    x_keepalive = {
        "coeff_vars_x": coeff_vars_x,
        "coeff_list_x": coeff_list_x,
        "bkg_x_lin": bkg_x_lin,
    }

    if upsilon_params is None:
        return bkg_x_lin, x_keepalive

    upsilon_frac = RooRealVar(f"{name}_upsilon_frac", f"{name}_upsilon_frac", 0.2, 0, 1)
    upsilon_model, _, _, _ = build_upsilon_model(x, sufix=name)

    set_pdf_parameters(
        upsilon_model,
        upsilon_params,
        x,
        make_constant=True,
        sufix=name,
    )

    bkg_x = RooAddPdf(
        f"{name}_x",
        "Chebychev + upsilon model background 2D - x",
        RooArgList(upsilon_model, bkg_x_lin),
        RooArgList(upsilon_frac),
        kTRUE,
    )
    x_keepalive.update(
        {
            "bkg_x": bkg_x,
            "upsilon_frac": upsilon_frac,
            "upsilon_model": upsilon_model,
        }
    )
    return bkg_x, x_keepalive


def build_background_cheb_2d(
    x,
    y,
    cheb_coeffs,
    upsilon_params=None,
    name="background_chebychev_2d",
    x_coeffs=None,
):
    name = f"{name}_{len(cheb_coeffs)}"
    if x_coeffs is None:
        x_coeffs = [0.0]
    coeff_vars_x = [
        RooRealVar(
            f"c_{name}_x{i}",
            f"Chebychev c_x{i}",
            float(v),
            -1.0,
            1.0,
        )
        for i, v in enumerate(x_coeffs, start=1)
    ]
    coeff_list_x = RooArgList()
    for v in coeff_vars_x:
        coeff_list_x.add(v)

    bkg_x, x_keepalive = _build_x_background_component(
        x,
        name,
        coeff_vars_x,
        coeff_list_x,
        upsilon_params,
    )

    coeff_vars_y = [
        RooRealVar(
            f"c_{name}_{i}",
            f"Chebychev c{i}",
            float(v),
            -1.0,
            1.0,
        )
        for i, v in enumerate(cheb_coeffs, start=1)
    ]
    coeff_list_y = RooArgList()
    for v in coeff_vars_y:
        coeff_list_y.add(v)

    bkg_y = RooChebychev(
        f"{name}_y",
        "Chebychev background 2D - y",
        y,
        coeff_list_y,
    )

    bkg = RooProdPdf(
        f"{name}",
        f"{name} background 2D",
        RooArgList(bkg_x, bkg_y),
    )

    bkg._keepalive = {
        **x_keepalive,
        "coeff_vars_y": coeff_vars_y,
        "coeff_list_y": coeff_list_y,
        "bkg_y": bkg_y,
    }

    return bkg


def build_background_bernstein_2d(
    x,
    y,
    bern_coeffs,
    upsilon_params=None,
    name="background_bernstein_2d",
    force_positive=True,
    x_coeffs=None,
):
    name = f"{name}_{len(bern_coeffs)}"
    lo = 0.0 if force_positive else -10.0
    hi = 10.0
    if x_coeffs is None:
        x_coeffs = [0.0]

    coeff_vars_x = [
        RooRealVar(
            f"c_{name}_x{i}",
            f"Chebychev c_x{i}",
            float(v),
            -1.0,
            1.0,
        )
        for i, v in enumerate(x_coeffs, start=1)
    ]
    coeff_list_x = RooArgList()
    for v in coeff_vars_x:
        coeff_list_x.add(v)

    bkg_x, x_keepalive = _build_x_background_component(
        x,
        name,
        coeff_vars_x,
        coeff_list_x,
        upsilon_params,
    )

    coeff_vars_y = [
        RooRealVar(
            f"b_{name}_{i}",
            f"Bernstein b{i}",
            float(v),
            lo,
            hi,
        )
        for i, v in enumerate(bern_coeffs)
    ]

    coeff_list_y = RooArgList()
    for v in coeff_vars_y:
        coeff_list_y.add(v)

    bkg_y = RooBernstein(
        f"{name}_y",
        "Bernstein background 2D - y",
        y,
        coeff_list_y,
    )

    bkg = RooProdPdf(
        f"{name}",
        f"{name} background 2D",
        RooArgList(bkg_x, bkg_y),
    )

    bkg._keepalive = {
        **x_keepalive,
        "coeff_vars_y": coeff_vars_y,
        "coeff_list_y": coeff_list_y,
        "bkg_y": bkg_y,
    }

    return bkg


def build_background_johnson_2d(
    x,
    y,
    dummy_coeffs,
    upsilon_params=None,
    name="background_johnson_2d",
    x_coeffs=None,
):
    name = name
    if x_coeffs is None:
        x_coeffs = [0.0]

    coeff_vars_x = [
        RooRealVar(
            f"c_{name}_x{i}",
            f"Chebychev c_x{i}",
            float(v),
            -1.0,
            1.0,
        )
        for i, v in enumerate(x_coeffs, start=1)
    ]
    coeff_list_x = RooArgList()
    for v in coeff_vars_x:
        coeff_list_x.add(v)

    bkg_x, x_keepalive = _build_x_background_component(
        x,
        name,
        coeff_vars_x,
        coeff_list_x,
        upsilon_params,
    )

    # Create parameters for the Johnson PDF
    mu = RooRealVar(f"mean_{name}", f"mean_{name}", 9.96, 0.01, 100)
    lam = RooRealVar(f"lambda_{name}", f"lambda_{name}", 24.3, 0.01, 100)  # > 0
    gamma = RooRealVar(f"gamma_{name}", f"gamma_{name}", -2.16, -10.0, 10.0)
    delta = RooRealVar(f"delta_{name}", f"delta_{name}", 1.42, 0.01, 20.0)  # > 0

    bkg_y = RooJohnson(
        f"{name}_y",
        "Johnson background 2D - y",
        y,
        mu,
        lam,
        gamma,
        delta,
    )

    bkg = RooProdPdf(
        f"{name}",
        f"{name} background 2D",
        RooArgList(bkg_x, bkg_y),
    )

    bkg._keepalive = {
        **x_keepalive,
        "mu": mu,
        "lam": lam,
        "gamma": gamma,
        "delta": delta,
        "bkg_y": bkg_y,
    }

    return bkg


def build_background_power_law_2d(
    x,
    y,
    exponents,
    upsilon_params=None,
    name="background_power_law_2d",
    exponents_bounds=(-10.0, +10.0),
    x_coeffs=None,
):
    name = f"{name}_{len(exponents)}"
    if not (len(exponents) >= 1):
        raise ValueError("At least one exponent should be provided")
    if x_coeffs is None:
        x_coeffs = [0.0]

    coeff_vars_x = [
        RooRealVar(
            f"c_{name}_x{i}",
            f"Chebychev c_x{i}",
            float(v),
            -1.0,
            1.0,
        )
        for i, v in enumerate(x_coeffs, start=1)
    ]
    coeff_list_x = RooArgList()
    for v in coeff_vars_x:
        coeff_list_x.add(v)

    bkg_x, x_keepalive = _build_x_background_component(
        x,
        name,
        coeff_vars_x,
        coeff_list_x,
        upsilon_params,
    )

    formula = _get_power_law_formula(exponents).replace("x", "boson_mass")
    print(formula)

    # parameters
    a_vars_y = []
    for i, v in enumerate(exponents):
        a_vars_y.append(
            RooRealVar(
                f"a_{name}_{i}",
                f"Power law a{i}",
                float(v),
                *exponents_bounds,
            )
        )

    args_y = RooArgList(y)
    for i in range(len(exponents)):
        args_y.add(a_vars_y[i])

    bkg_y = RooGenericPdf(
        f"{name}_y",
        # title,
        "Power-law background 2D - y",
        formula,
        args_y,
    )

    bkg = RooProdPdf(
        f"{name}",
        f"{name} background 2D",
        RooArgList(bkg_x, bkg_y),
    )

    bkg._keepalive = {
        **x_keepalive,
        "a_vars_y": a_vars_y,
        "args_y": args_y,
        "bkg_y": bkg_y,
    }
    return bkg


def build_background_models(
    pdf_family: BkgPdfFamily,
    x,
    n_coeffs: int = 6,
    min_n_coeffs: int = 1,
    initial_coeff=None,
):
    test_bkg_pdfs = {}
    test_bkg_pdfs[pdf_family] = []

    model_builder = None

    if pdf_family == BkgPdfFamily.CHEBYCHEV:
        model_builder = build_background_cheb

    if pdf_family == BkgPdfFamily.BERNSTEIN:
        model_builder = build_background_bernstein

    if pdf_family == BkgPdfFamily.POWER_LAW:
        model_builder = build_background_power_law

    if pdf_family == BkgPdfFamily.EXPONENTIAL:
        model_builder = build_background_exponential

    if model_builder == None:
        raise ValueError("Invalid PDF Family")

    if initial_coeff == None:
        initial_coeff = 0.0

    for i in range(min_n_coeffs, n_coeffs):
        test_bkg_pdfs[pdf_family].append(
            BkgModel(
                model=model_builder(
                    x,
                    [initial_coeff] * i,
                ),
                pdf_family=pdf_family,
                scan_order=i,
            )
        )

    return test_bkg_pdfs


def build_background_models_2d(
    pdf_family: BkgPdfFamily,
    x,
    y,
    upsilon_params=None,
    n_coeffs: int = 8,
    min_n_coeffs: int = 1,
    initial_coeff=None,
    x_initial_coeff: float = 0.0,
):
    test_bkg_pdfs = {}
    test_bkg_pdfs[pdf_family] = []

    model_builder = None

    if pdf_family == BkgPdfFamily.CHEBYCHEV:
        model_builder = build_background_cheb_2d

    if pdf_family == BkgPdfFamily.BERNSTEIN:
        model_builder = build_background_bernstein_2d

    if pdf_family == BkgPdfFamily.POWER_LAW:
        model_builder = build_background_power_law_2d

    if pdf_family == BkgPdfFamily.JOHNSON:
        model_builder = build_background_johnson_2d
    #
    # if pdf_family == BkgPdfFamily.EXPONENTIAL:
    #     model_builder = build_background_exponential_2d

    if model_builder == None:
        raise ValueError("Invalid PDF Family")

    if initial_coeff == None:
        initial_coeff = 0.0

    if pdf_family != BkgPdfFamily.JOHNSON:
        for i in range(min_n_coeffs, n_coeffs):
            test_bkg_pdfs[pdf_family].append(
                BkgModel(
                    model=model_builder(
                        x,
                        y,
                        [initial_coeff] * i,
                        upsilon_params,
                        x_coeffs=[x_initial_coeff],
                    ),
                    pdf_family=pdf_family,
                    scan_order=i,
                )
            )
    else:
        test_bkg_pdfs[pdf_family].append(
            BkgModel(
                model=model_builder(
                    x,
                    y,
                    [0.0],
                    upsilon_params,
                    x_coeffs=[x_initial_coeff],
                ),
                pdf_family=BkgPdfFamily.JOHNSON,
            )
        )

    return test_bkg_pdfs


def build_background_model_2d_candidate(
    pdf_family: BkgPdfFamily,
    x,
    y,
    *,
    scan_order: Optional[int],
    upsilon_params=None,
    initial_coeff: float = 0.0,
    x_initial_coeff: float = 0.0,
):
    if pdf_family == BkgPdfFamily.JOHNSON:
        return build_background_johnson_2d(
            x,
            y,
            [0.0],
            upsilon_params,
            x_coeffs=[x_initial_coeff],
        )

    if scan_order is None:
        raise ValueError(f"scan_order is required for {pdf_family}")

    if pdf_family == BkgPdfFamily.CHEBYCHEV:
        return build_background_cheb_2d(
            x,
            y,
            [initial_coeff] * scan_order,
            upsilon_params,
            x_coeffs=[x_initial_coeff],
        )

    if pdf_family == BkgPdfFamily.BERNSTEIN:
        return build_background_bernstein_2d(
            x,
            y,
            [initial_coeff] * scan_order,
            upsilon_params,
            x_coeffs=[x_initial_coeff],
        )

    if pdf_family == BkgPdfFamily.POWER_LAW:
        return build_background_power_law_2d(
            x,
            y,
            [initial_coeff] * scan_order,
            upsilon_params,
            x_coeffs=[x_initial_coeff],
        )

    raise ValueError(f"Invalid PDF Family: {pdf_family}")
