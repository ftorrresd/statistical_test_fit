from dataclasses import dataclass
import tempfile
import os
from typing import Any, Optional

from ROOT import (  # type: ignore
    RooArgSet,  # type: ignore
    RooChebychev,  # type: ignore
    RooRealVar,  # type: ignore
    RooArgList,  # type: ignore
    RooProdPdf,  # type: ignore
    RooExpPoly,  # type: ignore
    RooGenericPdf,  # type: ignore
    RooBernstein,  # type: ignore
    TMath,  # type: ignore
)

from .bkg_pdf_families import BkgPdfFamily
from .chi2_test import ChiSquareResult


@dataclass
class BkgModel:
    pdf_family: BkgPdfFamily
    model: Any
    n_params: Optional[int] = None
    chi_square_res: Optional[ChiSquareResult] = None
    NLL: Optional[float] = None
    fit_res: Optional[Any] = None

    def __post_init__(self):
        self.n_params = self.model.getParameters(RooArgSet()).getSize()

    def is_complete(self):
        return (
            (self.chi_square_res is not None)
            and (self.NLL is not None)
            and (self.n_params is not None)
            and (self.fit_res is not None)
        )

    def __hash__(self) -> int:
        return hash((self.pdf_family, self.chi_square_res, self.n_params, self.NLL))

    def __str__(self) -> str:
        self.model.Print()
        string_repr = []
        string_repr.append(f"Polinomial family: {self.pdf_family}")
        string_repr.append(f"N params: {self.n_params}")
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

        return "\n".join(string_repr)


def compute_winner_and_start_indexes(
    test_bkg_models: list[BkgModel],
) -> tuple[Optional[int], Optional[int]]:
    start = None
    winner = None
    for idx, test_bkg_model in enumerate(test_bkg_models):
        assert test_bkg_model.is_complete()
        if start is None:
            assert test_bkg_model.chi_square_res is not None
            if test_bkg_model.chi_square_res.pvalue > 0.05:
                start = idx
                winner = idx

        if start is not None:
            if idx > 0 and idx > start:
                previous_model = test_bkg_models[idx - 1]

                assert previous_model.NLL is not None and test_bkg_model.NLL is not None
                delta_NLL = abs(test_bkg_model.NLL - previous_model.NLL)

                assert (
                    previous_model.n_params is not None
                    and test_bkg_model.n_params is not None
                )
                delta_n_params = test_bkg_model.n_params - previous_model.n_params

                pval = TMath.Prob(delta_NLL, delta_n_params)
                # print(f"DEBUG: {delta_NLL, delta_n_params, pval}")
                if pval > 0.05:
                    winner = idx - 1
                    break

    if start is not None and winner is not None:
        return start, winner

    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("WARNING: Could not find best polynomial order!")
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!")

    return (
        len(test_bkg_models) - 1,
        len(test_bkg_models) - 1,
    )


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


def build_background_cheb_2d(x, y, cheb_coeffs, name="background_chebychev_2d"):
    coeff_vars_x = [
        RooRealVar(f"c_x{i}", f"Chebychev c_x{i}", float(v), -1.0, 1.0)
        for i, v in enumerate([0.3], start=1)
    ]
    coeff_list_x = RooArgList()
    for v in coeff_vars_x:
        coeff_list_x.add(v)
    bkg_x = RooChebychev(f"{name}_x", "Chebychev background 2D - x", x, coeff_list_x)

    coeff_vars_y = [
        RooRealVar(f"c{i}", f"Chebychev c{i}", float(v), -1.0, 1.0)
        for i, v in enumerate(cheb_coeffs, start=1)
    ]
    coeff_list_y = RooArgList()
    for v in coeff_vars_y:
        coeff_list_y.add(v)
    bkg_y = RooChebychev(f"{name}_y", "Chebychev background 2D - y", y, coeff_list_y)

    bkg = RooProdPdf("model", "Chebychev background 2D", RooArgList(bkg_x, bkg_y))

    bkg._keepalive = {
        "coeff_vars_x": coeff_vars_x,
        "coeff_vars_y": coeff_vars_y,
        "coeff_list_x": coeff_list_x,
        "coeff_list_y": coeff_list_y,
        "bkg_x": bkg_x,
        "bkg_y": bkg_y,
    }

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
            )
        )

    return test_bkg_pdfs


def build_background_models_2d(
    pdf_family: BkgPdfFamily,
    x,
    y,
    n_coeffs: int = 6,
    min_n_coeffs: int = 1,
    initial_coeff=None,
):
    test_bkg_pdfs = {}
    test_bkg_pdfs[pdf_family] = []

    model_builder = None

    if pdf_family == BkgPdfFamily.CHEBYCHEV:
        model_builder = build_background_cheb_2d

    # if pdf_family == BkgPdfFamily.BERNSTEIN:
    #     model_builder = build_background_bernstein
    #
    # if pdf_family == BkgPdfFamily.POWER_LAW:
    #     model_builder = build_background_power_law
    #
    # if pdf_family == BkgPdfFamily.EXPONENTIAL:
    #     model_builder = build_background_exponential

    if model_builder == None:
        raise ValueError("Invalid PDF Family")

    if initial_coeff == None:
        initial_coeff = 0.0

    for i in range(min_n_coeffs, n_coeffs):
        test_bkg_pdfs[pdf_family].append(
            BkgModel(
                model=model_builder(
                    x,
                    y,
                    [initial_coeff] * i,
                ),
                pdf_family=pdf_family,
            )
        )

    return test_bkg_pdfs
