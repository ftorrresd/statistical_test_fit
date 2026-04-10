from pprint import pprint
import os
import time
from functools import wraps
from argparse import Namespace

from ROOT import (  # type: ignore
    RooArgList,  # type: ignore
    RooArgSet,  # type: ignore
    RooCategory,  # type: ignore
    RooDataSet,  # type: ignore
    RooEffProd,  # type: ignore
    RooFit,  # type: ignore
    RooFormulaVar,  # type: ignore
    RooMultiPdf,  # type: ignore
    RooRealVar,  # type: ignore
    RooWorkspace,  # type: ignore
    TFile,  # type: ignore
)

from .bkg_model import (
    BkgModel,
    build_background_cheb_2d,
    build_background_models_2d,
    compute_winner_and_start_indexes,
    update_fit_quality,
)
from .bkg_pdf_families import BkgPdfFamily
from .chi2_test import ChiSquareResult
from .dimuon_non_correlated import dimuon_non_correlated
from .mass_ranges import (
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
from .make_plots import DataType, ProjDim, make_plots_2d


def execution_time(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        formatted_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
        print(f"Execution time of {func.__name__}: {formatted_time}")
        return result

    return wrapper


@execution_time
def run_fit_2d_data(args: Namespace):
    w = RooWorkspace("ws")

    # Limits
    upsilon_mass_lower = UPSILON_MASS_LOWER
    upsilon_mass_upper = UPSILON_MASS_UPPER

    left_lower = LEFT_SIDEBAND_LOWER
    left_upper = LEFT_SIDEBAND_UPPER

    middle_lower = MIDDLE_SIDEBAND_LOWER
    middle_upper = MIDDLE_SIDEBAND_UPPER

    right_lower = RIGHT_SIDEBAND_LOWER
    right_upper = RIGHT_SIDEBAND_UPPER

    outprefix = "bkg_only"

    LOAD_FROM_CACHE = args.use_cache
    if not LOAD_FROM_CACHE:
        os.system("rm -rf *.json")

    # build upsilon models
    upsilon_params = dimuon_non_correlated(
        upsilon_mass_lower, upsilon_mass_upper, load_from_cache=LOAD_FROM_CACHE
    )
    print("\n\nUpsilon Parameters:")
    pprint(upsilon_params)
    print("\n\n")

    # Observable
    upsilon_mass = RooRealVar(
        "upsilon_mass", "upsilon_mass", upsilon_mass_lower, upsilon_mass_upper
    )
    upsilon_mass.SetTitle("m_{#mu#mu}")  # LaTeX-style title
    upsilon_mass.setUnit("GeV")  # physical unit

    boson_mass = RooRealVar(
        "boson_mass", "boson_mass", BOSON_MASS_LOWER, BOSON_MASS_UPPER
    )
    boson_mass.SetTitle("m_{#mu#mu#gamma}")  # LaTeX-style title
    boson_mass.setUnit("GeV")  # physical unit

    f = TFile.Open("inputs/mass_Run2.root")
    data_full = RooDataSet(
        "data_obs",
        "data_obs",
        RooArgSet(
            upsilon_mass,
            boson_mass,
        ),
        RooFit.Import(f.Events),
        RooFit.Cut(
            f"(upsilon_mass < {upsilon_mass_upper} && upsilon_mass >= {upsilon_mass_lower} && boson_mass>={left_lower} && boson_mass<{right_upper})"
        ),
    )
    getattr(w, "import")(data_full)

    # Named ranges for sidebands
    boson_mass.setRange("LEFT", left_lower, left_upper)
    boson_mass.setRange("MIDDLE", middle_lower, middle_upper)
    boson_mass.setRange("RIGHT", right_lower, right_upper)
    boson_mass.setRange("FULL", left_lower, right_upper)

    # Sideband-only dataset (useful for plotting)
    cut_expr = f"((boson_mass<{left_upper}) || (boson_mass>{right_lower}) || ((boson_mass>{middle_lower}) && (boson_mass<{middle_upper})))"
    # data_sb = data_full.reduce(RooFit.Cut(cut_expr))
    data_sb = data_full.reduce(RooFit.Cut(cut_expr))
    print(f"Sideband entries: {data_sb.numEntries()} (out of {data_full.numEntries()})")

    test_bkg_pdfs: dict[BkgPdfFamily, list[BkgModel]] = {}
    for family in BkgPdfFamily:
        test_bkg_pdfs[family] = []

    # build the many bkg models
    test_bkg_pdfs |= build_background_models_2d(
        BkgPdfFamily.JOHNSON,
        upsilon_mass,
        boson_mass,
        # min_n_coeffs=1,
        # n_coeffs=3,
        upsilon_params=upsilon_params,
    )
    test_bkg_pdfs |= build_background_models_2d(
        BkgPdfFamily.CHEBYCHEV,
        upsilon_mass,
        boson_mass,
        min_n_coeffs=3,
        n_coeffs=10,
        upsilon_params=upsilon_params,
    )
    test_bkg_pdfs |= build_background_models_2d(
        BkgPdfFamily.BERNSTEIN,
        upsilon_mass,
        boson_mass,
        min_n_coeffs=3,
        n_coeffs=10,
        upsilon_params=upsilon_params,
    )
    # test_bkg_pdfs |= build_background_models_2d(
    #     BkgPdfFamily.POWER_LAW,
    #     upsilon_mass,
    #     boson_mass,
    #     min_n_coeffs=1,
    #     n_coeffs=2,
    #     upsilon_params=upsilon_params,
    #     initial_coeff=4.0,
    # )
    # test_bkg_pdfs |= build_background_models(
    #     BkgPdfFamily.EXPONENTIAL, x, initial_coeff=-150_000
    # )

    winners = {}
    for family in BkgPdfFamily:
        if len(test_bkg_pdfs[family]) > 0:
            print()
            print("##############################################")
            print("##############################################")
            print(f"################## {str(family).upper()} ########################")
            print("##############################################")
            print("##############################################")
            print()
            for test_bkg_pdf in test_bkg_pdfs[family]:
                # sideband acceptance (depends only on boson_mass)
                lu = RooRealVar("lu", "left_upper", left_upper)
                rl = RooRealVar("rl", "right_lower", right_lower)
                ml = RooRealVar("ml", "middle_lower", middle_lower)
                mu = RooRealVar("mu", "middle_upper", middle_upper)
                for v in (lu, rl, ml, mu):
                    v.setConstant(True)

                # acceptance: 1 in sidebands, 0 in the blinded window
                acc = RooFormulaVar(
                    "acc",
                    "(boson_mass < lu) || (boson_mass > rl) || "
                    "((boson_mass > ml) && (boson_mass < mu))",
                    RooArgList(boson_mass, lu, rl, ml, mu),
                )

                # Apply acceptance to your 2D pdf (shape depends on upsilon_mass, boson_mass)
                # eff_pdf(upsilon_mass, boson_mass) ∝ pdf(upsilon_mass, boson_mass) * acc(boson_mass)
                eff_pdf = RooEffProd(
                    "eff_pdf", "pdf * acceptance", test_bkg_pdf.model, acc
                )

                print(f"Will fit PDF: {test_bkg_pdf.model.GetName()}")
                test_bkg_pdf.fit_res = eff_pdf.fitTo(
                    data_sb,
                    RooFit.Save(True),
                    RooFit.PrintLevel(-1),
                    RooFit.Verbose(False),
                    RooFit.Minimizer("Minuit2", "Migrad"),
                )
                test_bkg_pdf.fit_res.Print("v")
                print("... done fitting.")
                test_bkg_pdf.n_float_params = (
                    test_bkg_pdf.fit_res.floatParsFinal().getSize()
                )

                os.system(
                    f"mkdir -p plots/fit_2d_data/control/{str(test_bkg_pdf.pdf_family).replace(' ', '_')}"
                )

                test_bkg_pdf.chi_square_res = ChiSquareResult.compute_chi_square_2d(
                    test_bkg_pdf.model,
                    data_sb,
                    upsilon_mass,
                    boson_mass,
                    outprefix=f"test_bkg_pdf_{test_bkg_pdf.n_params}",
                    pdf_family=family,
                    nbins=args.nbins,
                    is_data=True,
                    nfloatpars=test_bkg_pdf.n_float_params,
                )

                # test_bkg_pdf.NLL = test_bkg_pdf.model.createNLL(data_sb).getVal()
                test_bkg_pdf.NLL = eff_pdf.createNLL(data_sb).getVal()
                update_fit_quality(test_bkg_pdf)

            print("\n\n=== Test Background-only fit (sidebands) ===")
            for i, test_bkg_pdf in enumerate(test_bkg_pdfs[family]):
                assert test_bkg_pdf.is_complete()
                if i != 0:
                    print("")

                print(test_bkg_pdf)

            # compute winner function
            start, winner = compute_winner_and_start_indexes(
                test_bkg_pdfs[family],
                strict_mode=args.strict_mode,
            )
            winners[family] = winner

            # Plots
            plot_file_name = make_plots_2d(
                ProjDim.Y,
                upsilon_mass,
                boson_mass,
                # data_full,
                None,
                data_sb,
                None,
                test_bkg_pdfs[family],
                family,
                left_lower,
                left_upper,
                middle_lower,
                middle_upper,
                right_lower,
                right_upper,
                f"{outprefix}_{family}",
                nbins=args.nbins,
                data_type=DataType.REAL,
                start=start,
                winner=winner,
                # components=[("Z Resonant BKG", test_bkg_pdfs[family][0].model)],
            )
            print(f"\nPlot saved to: {plot_file_name}")

            plot_file_name = make_plots_2d(
                ProjDim.X,
                upsilon_mass,
                boson_mass,
                # data_full,
                None,
                data_sb,
                None,
                test_bkg_pdfs[family],
                family,
                left_lower,
                left_upper,
                middle_lower,
                middle_upper,
                right_lower,
                right_upper,
                f"{outprefix}_{family}",
                nbins=args.nbins,
                data_type=DataType.REAL,
                start=start,
                winner=winner,
            )
            print(f"\nPlot saved to: {plot_file_name}")

            # Optional: estimate expected background count in the blinded window (from fitted pdf)
            # Compute integral of background over blind window, normalized over full x.
            boson_mass.setRange("BLIND_LEFT", left_upper, middle_lower)
            boson_mass.setRange("BLIND_RIGHT", middle_upper, right_lower)

    pdf_index = RooCategory("pdfIndex", "pdfIndex")
    getattr(w, "import")(pdf_index)
    pdf_list = RooArgList()
    for family in BkgPdfFamily:
        if len(test_bkg_pdfs[family]) > 0:
            w.cat("pdfIndex").defineType(family.value)
            print(f"--> Adding winner for {family}")
            getattr(w, "import")(test_bkg_pdfs[family][winners[family]].model)
            # pdf_list.add(w.pdf(family.value))
            pdf_list.add(w.pdf(test_bkg_pdfs[family][winners[family]].model.GetName()))

    multi_pdf = RooMultiPdf(
        "multiPdf",
        "multiPdf",
        w.cat("pdfIndex"),
        pdf_list,
    )
    getattr(w, "import")(multi_pdf)

    w.Print("v")
    print("pdfIndex: ", w.cat("pdfIndex").states())
