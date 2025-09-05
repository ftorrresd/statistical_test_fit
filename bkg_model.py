from dataclasses import dataclass
from typing import Any, Optional

from ROOT import (
    RooArgSet,  # type: ignore
    TMath,  # type: ignore
)

from bkg_pdf_families import BkgPdfFamily
from chi2_test import ChiSquareResult


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
        string_repr = []
        string_repr.append("=======")
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
        string_repr.append("=======")

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
