import sys
from argparse import Namespace

from ROOT import (  # type: ignore
    RooAddPdf,  # type: ignore
    RooArgList,  # type: ignore
    RooArgSet,  # type: ignore
    RooFit,  # type: ignore
    RooGaussian,  # type: ignore
    RooProdPdf,  # type: ignore
    RooRandom,  # type: ignore
    RooRealVar,  # type: ignore
)

from .bkg_model import (
    BkgModel,
    build_background_cheb_2d,
    build_background_models_2d,
    compute_winner_and_start_indexes,
)
from .bkg_pdf_families import BkgPdfFamily
from .chi2_test import ChiSquareResult
from .make_plots import ProjDim, make_plots_2d


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


def run_fit_2d(args: Namespace):
    # Limits
    m_mumu_lower = 8.0
    m_mumu_upper = 12.0

    left_lower = 70.0
    left_upper = 80.0

    middle_lower = 100.0
    middle_upper = 115.0

    right_lower = 135.0
    right_upper = 200.0

    z_sigfrac = 0.05
    h_sigfrac = 0.05 / 2.0

    outprefix = "bkg_only"

    # Observable
    m_mumu = RooRealVar("m_mumu", "m_mumu", m_mumu_lower, m_mumu_upper)
    m_mumugamma = RooRealVar("m_mumugamma", "m_mumugamma", left_lower, right_upper)

    observables = RooArgSet(m_mumu, m_mumugamma)

    # Background pdf
    cheb_coeffs = [float(v.strip()) for v in args.cheb.split(",") if v.strip() != ""]
    if not cheb_coeffs:
        print("Provide at least one Chebychev coefficient via --cheb", file=sys.stderr)
        sys.exit(1)
    bkg_pdf = BkgModel(
        model=build_background_cheb_2d(m_mumu, m_mumugamma, cheb_coeffs),
        pdf_family=BkgPdfFamily.CHEBYCHEV,
    )

    # Generation model (S+B, unextended)
    z_model = build_signal(m_mumu, m_mumugamma, 10.0, 0.2, 91.0, 2.0, "z_signal")
    higgs_model = build_signal(m_mumu, m_mumugamma, 10.0, 0.02, 125.0, 0.3, "h_signal")
    gen_model = build_gen_model_unextended(
        z_model,
        higgs_model,
        bkg_pdf.model,
        z_sigfrac,
        h_sigfrac,
    )

    # Random seed
    if args.seed:
        RooRandom.randomGenerator().SetSeed(args.seed)
    else:
        print("INFO: No seed was provided.")

    # Generate exactly N events (UNEXTENDED)
    data_full = gen_model.generate(observables, int(args.events))
    print(
        f"Generated dataset with {data_full.numEntries()} events (true z_fsig={z_sigfrac:.3f} and true h_fsig={h_sigfrac:.3f})"
    )

    # Named ranges for sidebands
    m_mumugamma.setRange("left", left_lower, left_upper)
    m_mumugamma.setRange("middle", middle_lower, middle_upper)
    m_mumugamma.setRange("right", right_lower, right_upper)
    m_mumu.setRange("left", m_mumu_lower, m_mumu_upper)
    m_mumu.setRange("middle", m_mumu_lower, m_mumu_upper)
    m_mumu.setRange("right", m_mumu_lower, m_mumu_upper)

    m_mumugamma.setRange("full", left_lower, right_upper)
    m_mumu.setRange("full", m_mumu_lower, m_mumu_upper)

    # Sideband-only dataset (useful for plotting)
    cut_expr = f"((m_mumugamma<{left_upper}) || (m_mumugamma>{right_lower}) || ((m_mumugamma>{middle_lower}) && (m_mumugamma<{middle_upper})))"
    data_sb = data_full.reduce(RooFit.Cut(cut_expr))
    print(f"Sideband entries: {data_sb.numEntries()} (out of {data_full.numEntries()})")

    # Print a short summary
    print("\n\n=== Background-only fit (sidebands) ===")
    print(bkg_pdf)

    test_bkg_pdfs: dict[BkgPdfFamily, list[BkgModel]] = {}
    for family in BkgPdfFamily:
        test_bkg_pdfs[family] = []

    # build the many bkg models
    test_bkg_pdfs |= build_background_models_2d(
        BkgPdfFamily.CHEBYCHEV,
        m_mumu,
        m_mumugamma,
        # min_n_coeffs=1,
        # n_coeffs=2,
    )
    # test_bkg_pdfs |= build_background_models(BkgPdfFamily.BERNSTEIN, x)
    # test_bkg_pdfs |= build_background_models(BkgPdfFamily.POWER_LAW, x)
    # test_bkg_pdfs|=build_background_models(BkgPdfFamily.EXPONENTIAL, x, initial_coeff=-150_000)

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
                test_bkg_pdf.fit_res = test_bkg_pdf.model.fitTo(
                    data_sb,
                    RooFit.Range("middle"),
                    # data_full,
                    # RooFit.Range(
                    #     # "left,middle,right"
                    #     "middle"
                    # ),  # ...but restrict the likelihood to sidebands
                    RooFit.Save(True),
                    RooFit.PrintLevel(-1),
                    RooFit.Verbose(False),
                    RooFit.Minimizer("Minuit2", "Migrad"),
                )

                test_bkg_pdf.chi_square_res = ChiSquareResult.compute_chi_square_2d(
                    test_bkg_pdf.model,
                    data_sb,
                    m_mumu,
                    m_mumugamma,
                    outprefix=f"test_bkg_pdf_{test_bkg_pdf.n_params}",
                    pdf_family=family,
                    nbins=args.nbins,
                )

                test_bkg_pdf.NLL = test_bkg_pdf.model.createNLL(
                    data_sb,
                ).getVal()

            print("\n\n=== Test Background-only fit (sidebands) ===")
            for i, test_bkg_pdf in enumerate(test_bkg_pdfs[family]):
                assert test_bkg_pdf.is_complete()
                if i != 0:
                    print("")

                print(test_bkg_pdf)

            # compute winner function
            start, winner = compute_winner_and_start_indexes(test_bkg_pdfs[family])

            # Plots
            plot_file_name = make_plots_2d(
                ProjDim.Y,
                m_mumu,
                m_mumugamma,
                data_full,
                data_sb,
                bkg_pdf,
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
                start=start,
                winner=winner,
            )
            print(f"\nPlot saved to: {plot_file_name}")

            plot_file_name = make_plots_2d(
                ProjDim.X,
                m_mumu,
                m_mumugamma,
                data_full,
                data_sb,
                bkg_pdf,
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
                start=start,
                winner=winner,
            )
            print(f"\nPlot saved to: {plot_file_name}")

            # Optional: estimate expected background count in the blinded window (from fitted pdf)
            # Compute integral of background over blind window, normalized over full x.
            m_mumugamma.setRange("blind_left", left_upper, middle_lower)
            m_mumugamma.setRange("blind_right", middle_upper, right_lower)

            # For an expected count, scale by the number of observed sideband events with
            # the pdf’s fraction in blind vs sideband ranges.
            i_blind = bkg_pdf.model.createIntegral(
                RooArgSet(m_mumugamma), RooFit.Range("blind_left,blind_right")
            ).getVal()

            i_full = bkg_pdf.model.createIntegral(
                RooArgSet(m_mumugamma),
                RooFit.Range("left,blind_left,middle,blind_right,right"),
            ).getVal()

            frac_blind = i_blind / i_full if i_full > 0 else 0.0
            n_est = frac_blind * data_full.numEntries()
            print(
                f"\nBackground estimate in blind window (from fitted pdf): "
                f"N_bkg(blind) ≈ {n_est:.1f} (fraction={frac_blind:.4f})"
            )
