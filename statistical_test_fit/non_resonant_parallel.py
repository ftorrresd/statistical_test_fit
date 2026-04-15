from dataclasses import dataclass
from typing import Dict, List, Optional

from .bkg_pdf_families import BkgPdfFamily
from .parallel_utils import ParallelJob


@dataclass(frozen=True)
class NonResonantCandidateSpec:
    pdf_family: BkgPdfFamily
    scan_order: Optional[int]
    upsilon_params: Dict[str, float]
    nbins: int
    input_file: str
    upsilon_mass_lower: float
    upsilon_mass_upper: float
    boson_mass_lower: float
    boson_mass_upper: float
    left_lower: float
    left_upper: float
    middle_lower: float
    middle_upper: float
    right_lower: float
    right_upper: float
    initial_coeff: float = 0.0
    x_initial_coeff: float = 0.0

    @property
    def key(self) -> str:
        order = "single" if self.scan_order is None else str(self.scan_order)
        return f"{self.pdf_family.value}_{order}"

    @property
    def label(self) -> str:
        if self.scan_order is None:
            return str(self.pdf_family)
        return f"{self.pdf_family} order {self.scan_order}"


@dataclass(frozen=True)
class NonResonantCandidateResult:
    spec: NonResonantCandidateSpec
    model_name: str
    params_dict: Dict[str, float]
    chi_square_result: Dict[str, float]
    NLL: float
    n_params: int
    n_float_params: int
    fit_status: int
    fit_cov_qual: int
    fit_edm: Optional[float]
    fit_ok: bool
    fit_quality_reason: str


def build_non_resonant_candidate_specs(
    *,
    upsilon_params: Dict[str, float],
    nbins: int,
    input_file: str,
    upsilon_mass_lower: float,
    upsilon_mass_upper: float,
    boson_mass_lower: float,
    boson_mass_upper: float,
    left_lower: float,
    left_upper: float,
    middle_lower: float,
    middle_upper: float,
    right_lower: float,
    right_upper: float,
) -> List[NonResonantCandidateSpec]:
    specs = [
        NonResonantCandidateSpec(
            pdf_family=BkgPdfFamily.JOHNSON,
            scan_order=None,
            upsilon_params=upsilon_params,
            nbins=nbins,
            input_file=input_file,
            upsilon_mass_lower=upsilon_mass_lower,
            upsilon_mass_upper=upsilon_mass_upper,
            boson_mass_lower=boson_mass_lower,
            boson_mass_upper=boson_mass_upper,
            left_lower=left_lower,
            left_upper=left_upper,
            middle_lower=middle_lower,
            middle_upper=middle_upper,
            right_lower=right_lower,
            right_upper=right_upper,
        )
    ]

    for family in (BkgPdfFamily.CHEBYCHEV, BkgPdfFamily.BERNSTEIN):
        for order in range(3, 10):
            specs.append(
                NonResonantCandidateSpec(
                    pdf_family=family,
                    scan_order=order,
                    upsilon_params=upsilon_params,
                    nbins=nbins,
                    input_file=input_file,
                    upsilon_mass_lower=upsilon_mass_lower,
                    upsilon_mass_upper=upsilon_mass_upper,
                    boson_mass_lower=boson_mass_lower,
                    boson_mass_upper=boson_mass_upper,
                    left_lower=left_lower,
                    left_upper=left_upper,
                    middle_lower=middle_lower,
                    middle_upper=middle_upper,
                    right_lower=right_lower,
                    right_upper=right_upper,
                )
            )

    return specs


def build_non_resonant_candidate_jobs(
    specs: List[NonResonantCandidateSpec],
) -> List[ParallelJob]:
    return [ParallelJob(key=spec.key, label=spec.label, payload=spec) for spec in specs]


def fit_non_resonant_candidate_job(
    spec: NonResonantCandidateSpec,
) -> NonResonantCandidateResult:
    from .root_runtime import configure_root

    ROOT = configure_root()

    from ROOT import (  # type: ignore
        RooArgList,  # type: ignore
        RooArgSet,  # type: ignore
        RooDataSet,  # type: ignore
        RooEffProd,  # type: ignore
        RooFit,  # type: ignore
        RooFormulaVar,  # type: ignore
        RooRealVar,  # type: ignore
        TFile,  # type: ignore
    )

    from .bkg_model import (
        BkgModel,
        build_background_model_2d_candidate,
        update_fit_quality,
    )
    from .chi2_test import ChiSquareResult
    from .ws_helper import get_pdf_parameters

    upsilon_mass = RooRealVar(
        "upsilon_mass",
        "upsilon_mass",
        spec.upsilon_mass_lower,
        spec.upsilon_mass_upper,
    )
    boson_mass = RooRealVar(
        "boson_mass",
        "boson_mass",
        spec.boson_mass_lower,
        spec.boson_mass_upper,
    )
    boson_mass.setRange("LEFT", spec.left_lower, spec.left_upper)
    boson_mass.setRange("MIDDLE", spec.middle_lower, spec.middle_upper)
    boson_mass.setRange("RIGHT", spec.right_lower, spec.right_upper)
    boson_mass.setRange("FULL", spec.left_lower, spec.right_upper)

    root_file = TFile.Open(spec.input_file)
    if root_file is None or root_file.IsZombie():
        raise RuntimeError(f"Could not open input file {spec.input_file}")

    data_full = RooDataSet(
        "data_obs",
        "data_obs",
        RooArgSet(upsilon_mass, boson_mass),
        RooFit.Import(root_file.Events),
        RooFit.Cut(
            f"(upsilon_mass < {spec.upsilon_mass_upper} && upsilon_mass >= {spec.upsilon_mass_lower} && boson_mass>={spec.left_lower} && boson_mass<{spec.right_upper})"
        ),
    )
    data_sb = data_full.reduce(
        RooFit.Cut(
            f"((boson_mass<{spec.left_upper}) || (boson_mass>{spec.right_lower}) || ((boson_mass>{spec.middle_lower}) && (boson_mass<{spec.middle_upper})))"
        )
    )

    model = build_background_model_2d_candidate(
        spec.pdf_family,
        upsilon_mass,
        boson_mass,
        scan_order=spec.scan_order,
        upsilon_params=spec.upsilon_params,
        initial_coeff=spec.initial_coeff,
        x_initial_coeff=spec.x_initial_coeff,
    )
    test_bkg_pdf = BkgModel(
        model=model,
        pdf_family=spec.pdf_family,
        scan_order=spec.scan_order,
    )

    lu = RooRealVar("lu", "left_upper", spec.left_upper)
    rl = RooRealVar("rl", "right_lower", spec.right_lower)
    ml = RooRealVar("ml", "middle_lower", spec.middle_lower)
    mu = RooRealVar("mu", "middle_upper", spec.middle_upper)
    for var in (lu, rl, ml, mu):
        var.setConstant(True)

    acc = RooFormulaVar(
        "acc",
        "(boson_mass < lu) || (boson_mass > rl) || ((boson_mass > ml) && (boson_mass < mu))",
        RooArgList(boson_mass, lu, rl, ml, mu),
    )
    eff_pdf = RooEffProd("eff_pdf", "pdf * acceptance", test_bkg_pdf.model, acc)

    test_bkg_pdf.fit_res = eff_pdf.fitTo(
        data_sb,
        RooFit.Save(True),
        RooFit.PrintLevel(-1),
        RooFit.Verbose(False),
        RooFit.Minimizer("Minuit2", "Migrad"),
    )
    test_bkg_pdf.n_float_params = test_bkg_pdf.fit_res.floatParsFinal().getSize()

    control_dir = (
        f"plots/fit_2d_data/control/{str(test_bkg_pdf.pdf_family).replace(' ', '_')}"
    )
    import os

    os.makedirs(control_dir, exist_ok=True)

    chi_square_res = ChiSquareResult.compute_chi_square_2d(
        test_bkg_pdf.model,
        data_sb,
        upsilon_mass,
        boson_mass,
        outprefix=spec.key,
        pdf_family=spec.pdf_family,
        nbins=spec.nbins,
        is_data=True,
        nfloatpars=test_bkg_pdf.n_float_params,
    )
    test_bkg_pdf.chi_square_res = chi_square_res
    test_bkg_pdf.NLL = eff_pdf.createNLL(data_sb).getVal()
    update_fit_quality(test_bkg_pdf)

    params_dict = get_pdf_parameters(
        test_bkg_pdf.model,
        RooArgSet(upsilon_mass, boson_mass),
    )

    root_file.Close()

    return NonResonantCandidateResult(
        spec=spec,
        model_name=test_bkg_pdf.model.GetName(),
        params_dict=params_dict,
        chi_square_result={
            "chi2": chi_square_res.chi2,
            "pvalue": chi_square_res.pvalue,
            "ndf": chi_square_res.ndf,
            "chi2_ndf": chi_square_res.chi2_ndf,
        },
        NLL=test_bkg_pdf.NLL,
        n_params=test_bkg_pdf.n_params if test_bkg_pdf.n_params is not None else 0,
        n_float_params=test_bkg_pdf.n_float_params,
        fit_status=test_bkg_pdf.fit_status
        if test_bkg_pdf.fit_status is not None
        else -1,
        fit_cov_qual=test_bkg_pdf.fit_cov_qual
        if test_bkg_pdf.fit_cov_qual is not None
        else -1,
        fit_edm=test_bkg_pdf.fit_edm,
        fit_ok=bool(test_bkg_pdf.fit_ok),
        fit_quality_reason=test_bkg_pdf.fit_quality_reason or "",
    )


def reconstruct_non_resonant_candidate(
    result: NonResonantCandidateResult,
    upsilon_mass,
    boson_mass,
):
    from ROOT import RooArgSet  # type: ignore

    from .bkg_model import BkgModel, build_background_model_2d_candidate
    from .chi2_test import ChiSquareResult
    from .ws_helper import set_pdf_parameters

    model = build_background_model_2d_candidate(
        result.spec.pdf_family,
        upsilon_mass,
        boson_mass,
        scan_order=result.spec.scan_order,
        upsilon_params=result.spec.upsilon_params,
        initial_coeff=result.spec.initial_coeff,
        x_initial_coeff=result.spec.x_initial_coeff,
    )
    set_pdf_parameters(
        model,
        result.params_dict,
        RooArgSet(upsilon_mass, boson_mass),
        make_constant=True,
    )

    reconstructed = BkgModel(
        model=model,
        pdf_family=result.spec.pdf_family,
        scan_order=result.spec.scan_order,
    )
    reconstructed.n_float_params = result.n_float_params
    reconstructed.chi_square_res = ChiSquareResult(**result.chi_square_result)
    reconstructed.NLL = result.NLL
    reconstructed.fit_status = result.fit_status
    reconstructed.fit_cov_qual = result.fit_cov_qual
    reconstructed.fit_edm = result.fit_edm
    reconstructed.fit_ok = result.fit_ok
    reconstructed.fit_quality_reason = result.fit_quality_reason
    return reconstructed
