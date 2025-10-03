import os
import sys
from argparse import Namespace

from ROOT import (  # type: ignore
    RooAddPdf,  # type: ignore
    RooArgList,  # type: ignore
    RooArgSet,  # type: ignore
    RooCategory,  # type: ignore
    RooDataSet,  # type: ignore
    RooEffProd,  # type: ignore
    RooFit,  # type: ignore
    RooFormulaVar,  # type: ignore
    RooGaussian,  # type: ignore
    RooMultiPdf,  # type: ignore
    RooProdPdf,  # type: ignore
    RooRealVar,  # type: ignore
    RooWorkspace,  # type: ignore
    TFile,  # type: ignore
)

from .bkg_model import (
    BkgModel,
    build_background_cheb_2d,
    build_background_models_2d,
    compute_winner_and_start_indexes,
)
from .bkg_pdf_families import BkgPdfFamily
from .chi2_test import ChiSquareResult
from .dimuon_non_correlated import dimuon_non_correlated
from .make_plots import DataType, ProjDim, make_plots_2d
from .normalization_fit import fit_and_plot
from .resonant_bkg_modeling import (
    ControlRegion,
    get_normalization_from_CR,
    resonant_background_modeling_Higgs,
    resonant_background_modeling_Z,
)


def build_signal(x, y, mean_x, sigma_x, mean_y, sigma_y, name: str):
    mu_x = RooRealVar(
        "mu_x",
        "signal mean - x",
        mean_x,
        x.getMin(),
        x.getMax(),
    )
    sig_x = RooRealVar(
        "sigma_x", "signal sigma - x", sigma_x, sigma_x * (0.7), sigma_x * (1.3)
    )
    gaus_x = RooGaussian(f"{name}_x", "Gaussian signal - x", x, mu_x, sig_x)

    mu_y = RooRealVar(
        "mu_y",
        "signal mean - x",
        mean_y,
        y.getMin(),
        y.getMax(),
    )
    sig_y = RooRealVar(
        "sigma_y", "signal sigma - x", sigma_y, sigma_y * (0.7), sigma_y * (1.3)
    )
    gaus_y = RooGaussian(f"{name}_y", "Gaussian signal - y", y, mu_y, sig_y)

    sig = RooProdPdf(f"{name}_model", "Gaussian signal 2D", RooArgList(gaus_x, gaus_y))

    sig._keepalive = {
        "gaus_x": gaus_x,
        "gaus_y": gaus_y,
        "mu_x": mu_x,
        "sig_x": sig_x,
        "mu_y": mu_y,
        "sig_y": sig_y,
    }

    return sig


def build_gen_model_unextended(
    z_model,
    higgs_model,
    bkg: BkgModel,
    f_Zsig: float,
    f_Hsig: float,
):
    # Mixture model for GENERATION ONLY (fraction, not extended)
    fs_Z = RooRealVar("fzsig_gen", "Z signal fraction (gen)", float(f_Zsig), 0.0, 1.0)
    fs_H = RooRealVar("fhsig_gen", "H signal fraction (gen)", float(f_Hsig), 0.0, 1.0)
    model = RooAddPdf(
        "model_gen",
        "signal+background (gen)",
        RooArgList(z_model, higgs_model, bkg),
        RooArgList(f_Zsig, f_Hsig),
    )
    model._keepalive = {
        "z_model": z_model,
        "higgs_model": higgs_model,
        "bkg": bkg,
        "fs_Z": fs_Z,
        "fs_H": fs_H,
    }
    return model


def run_fit_2d_data(args: Namespace):
    w = RooWorkspace("ws")

    # Limits
    upsilon_mass_lower = 8.0
    upsilon_mass_upper = 12.0

    left_lower = 70.0
    left_upper = 80.0

    middle_lower = 100.0
    middle_upper = 115.0

    right_lower = 135.0
    right_upper = 200.0

    z_sigfrac = 0.05
    h_sigfrac = 0.05 / 2.0

    outprefix = "bkg_only"

    # # build upsilon models
    # upsilon_model = dimuon_non_correlated(upsilon_mass_lower, upsilon_mass_upper)
    #
    # # build higss resonant bkg
    # higgs_resonant_bkg_ws = resonant_background_modeling_Higgs()
    # higgs_resonant_bkg_ws.pdf("resonant_background_model").Print("v")
    #
    # Z_resonant_bkg_ws, Z_resonant_bkg_parameters = resonant_background_modeling_Z()
    # Z_resonant_bkg_ws.pdf("resonant_background_model").Print("v")
    #
    # normalizations_from_CR = []
    # normalizations_from_CR.append(
    #     get_normalization_from_CR(
    #         Z_resonant_bkg_parameters,
    #         ControlRegion.CR1,
    #     )
    # )
    # normalizations_from_CR.append(
    #     get_normalization_from_CR(
    #         Z_resonant_bkg_parameters,
    #         ControlRegion.CR2,
    #     )
    # )
    # normalizations_from_CR.append(
    #     get_normalization_from_CR(
    #         Z_resonant_bkg_parameters,
    #         ControlRegion.CR3,
    #     )
    # )
    # normalizations_from_CR.append(
    #     get_normalization_from_CR(
    #         Z_resonant_bkg_parameters,
    #         ControlRegion.CR4,
    #     )
    # )
    #
    # normalization_extrapolation = fit_and_plot(
    #     [c.value.midpoint for c in ControlRegion],
    #     [r["normalization"] for r in normalizations_from_CR],
    #     [r["normalization_unc"] for r in normalizations_from_CR],
    #     x0=10.0,
    #     output_pdf="plots/fit_2d_data/normalization_extrapolation.pdf",
    #     x_lines=[
    #         ControlRegion.CR1.value.upper,
    #         ControlRegion.CR2.value.lower,
    #         ControlRegion.CR2.value.upper,
    #         ControlRegion.CR3.value.upper,
    #     ],
    #     point_labels=["CR1", "CR2", "CR3", "CR4"],
    # )
    # print(f"Extrapolated normalization: {normalization_extrapolation}")
    #
    # Observable
    upsilon_mass = RooRealVar(
        "upsilon_mass", "upsilon_mass", upsilon_mass_lower, upsilon_mass_upper
    )
    upsilon_mass.SetTitle("m_{#mu#mu}")  # LaTeX-style title
    upsilon_mass.setUnit("GeV")  # physical unit

    boson_mass = RooRealVar("boson_mass", "boson_mass", left_lower, right_upper)
    boson_mass.SetTitle("m_{#mu#mu#gamma}")  # LaTeX-style title
    boson_mass.setUnit("GeV")  # physical unit

    observables = RooArgSet(upsilon_mass, boson_mass)

    f = TFile.Open("inputs/mass_Run2.root")
    data_full = RooDataSet(
        "data_obs",
        "data_obs",
        RooArgSet(boson_mass, upsilon_mass),
        RooFit.Import(f.Events),
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
        BkgPdfFamily.CHEBYCHEV,
        upsilon_mass,
        boson_mass,
        # min_n_coeffs=1,
        # n_coeffs=2,
    )
    test_bkg_pdfs |= build_background_models_2d(
        BkgPdfFamily.BERNSTEIN,
        upsilon_mass,
        boson_mass,
        # min_n_coeffs=1,
        # n_coeffs=2,
    )
    test_bkg_pdfs |= build_background_models_2d(
        BkgPdfFamily.POWER_LAW,
        upsilon_mass,
        boson_mass,
        # min_n_coeffs=1,
        # n_coeffs=2,
    )
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

                test_bkg_pdf.fit_res = eff_pdf.fitTo(
                    data_sb,
                    RooFit.Save(True),
                    RooFit.PrintLevel(-1),
                    RooFit.Verbose(False),
                    RooFit.Minimizer("Minuit2", "Migrad"),
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
                )

                # test_bkg_pdf.NLL = test_bkg_pdf.model.createNLL(data_sb).getVal()
                test_bkg_pdf.NLL = eff_pdf.createNLL(data_sb).getVal()

            print("\n\n=== Test Background-only fit (sidebands) ===")
            for i, test_bkg_pdf in enumerate(test_bkg_pdfs[family]):
                assert test_bkg_pdf.is_complete()
                if i != 0:
                    print("")

                print(test_bkg_pdf)

            # compute winner function
            start, winner = compute_winner_and_start_indexes(test_bkg_pdfs[family])
            winners[family] = winner

            # Plots
            plot_file_name = make_plots_2d(
                ProjDim.Y,
                upsilon_mass,
                boson_mass,
                data_full,
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

            plot_file_name = make_plots_2d(
                ProjDim.X,
                upsilon_mass,
                boson_mass,
                data_full,
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
            getattr(w, "import")(test_bkg_pdfs[family][winners[family]].model)
            pdf_list.add(w.pdf(family.value))

    multi_pdf = RooMultiPdf(
        "multiPdf",
        "multiPdf",
        w.cat("pdfIndex"),
        pdf_list,
    )
    getattr(w, "import")(multi_pdf)

    w.Print("v")
    print("pdfIndex: ", w.cat("pdfIndex").states())
