import sys
from argparse import Namespace

from ROOT import (  # type: ignore
    RooArgSet,  # type: ignore
    RooFit,  # type: ignore
    RooRandom,  # type: ignore
    RooRealVar,  # type: ignore
    RooGaussian,  # type: ignore
    RooAddPdf,  # type: ignore
    RooArgList,  # type: ignore
)

from .bkg_model import (
    BkgModel,
    compute_winner_and_start_indexes,
    build_background_cheb,
    build_background_models,
)
from .bkg_pdf_families import BkgPdfFamily
from .chi2_test import ChiSquareResult
from .make_plots import make_plots_2d


def build_signal(x, mean, sigma, name: str):
    mu = RooRealVar(
        "mu",
        "signal mean",
        mean,
        x.getMin() + 0.05 * (x.getMax() - x.getMin()),
        x.getMax() - 0.05 * (x.getMax() - x.getMin()),
    )
    sig = RooRealVar(
        "sigma",
        "signal sigma",
        sigma,
        0.02 * (x.getMax() - x.getMin()),
        0.5 * (x.getMax() - x.getMin()),
    )
    gaus = RooGaussian(name, "Gaussian signal", x, mu, sig)
    gaus._keepalive = {"mu": mu, "sig": sig}
    return gaus


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
    x = RooRealVar("x", "Observable x", left_lower, right_upper)

    # Background pdf
    cheb_coeffs = [float(v.strip()) for v in args.cheb.split(",") if v.strip() != ""]
    if not cheb_coeffs:
        print("Provide at least one Chebychev coefficient via --cheb", file=sys.stderr)
        sys.exit(1)
    bkg_pdf = BkgModel(
        model=build_background_cheb(x, cheb_coeffs),
        pdf_family=BkgPdfFamily.CHEBYCHEV,
    )

    # Generation model (S+B, unextended)
    z_model = build_signal(x, 91.0, 2.0, "z_signal")
    higgs_model = build_signal(x, 125.0, 0.3, "h_signal")
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
    data_full = gen_model.generate(RooArgSet(x), int(args.events))
    print(
        f"Generated dataset with {data_full.numEntries()} events (true z_fsig={z_sigfrac:.3f} and true h_fsig={h_sigfrac:.3f})"
    )

    # Named ranges for sidebands
    x.setRange("left", left_lower, left_upper)
    x.setRange("middle", middle_lower, middle_upper)
    x.setRange("right", right_lower, right_upper)

    # Sideband-only dataset (useful for plotting)
    cut_expr = f"((x<{left_upper}) || (x>{right_lower}) || ((x>{middle_lower}) && (x<{middle_upper})))"
    data_sb = data_full.reduce(RooFit.Cut(cut_expr))
    print(f"Sideband entries: {data_sb.numEntries()} (out of {data_full.numEntries()})")

    # Print a short summary
    print("\n\n=== Background-only fit (sidebands) ===")
    print(bkg_pdf)

    test_bkg_pdfs: dict[BkgPdfFamily, list[BkgModel]] = {}
    for family in BkgPdfFamily:
        test_bkg_pdfs[family] = []

    # build the many bkg models
    test_bkg_pdfs |= build_background_models(BkgPdfFamily.CHEBYCHEV, x)
    test_bkg_pdfs |= build_background_models(BkgPdfFamily.BERNSTEIN, x)
    test_bkg_pdfs |= build_background_models(BkgPdfFamily.POWER_LAW, x)
    # test_bkg_pdfs|=build_background_models(BkgPdfFamily.EXPONENTIAL, x, initial_coeff=-150_000)

    for family in BkgPdfFamily:
        print()
        print("##############################################")
        print("##############################################")
        print(f"################## {str(family).upper()} ########################")
        print("##############################################")
        print("##############################################")
        print()
        for test_bkg_pdf in test_bkg_pdfs[family]:
            # if test_bkg_pdf.pdf_family == BkgPdfFamily.EXPONENTIAL:
            #     test_bkg_pdf.model.Print("all")

            test_bkg_pdf.fit_res = test_bkg_pdf.model.fitTo(
                data_full,  # we can pass the full dataset...
                RooFit.Range(
                    "left,middle,right"
                ),  # ...but restrict the likelihood to sidebands
                RooFit.Save(True),
                RooFit.PrintLevel(-1),
                RooFit.Verbose(False),
                RooFit.Minimizer("Minuit2", "Migrad"),
            )

            test_bkg_pdf.chi_square_res = ChiSquareResult.compute_chi_square_2d(
                test_bkg_pdf.model,
                data_full,
                x,
                outprefix=f"test_bkg_pdf_{test_bkg_pdf.n_params}",
                pdf_family=family,
                nbins=args.nbins,
            )

            test_bkg_pdf.NLL = test_bkg_pdf.model.createNLL(
                data_full,
                RooFit.Range("left,middle,right"),
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
            x,
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
        x.setRange("blind_left", left_upper, middle_lower)
        x.setRange("blind_right", middle_upper, right_lower)

        # For an expected count, scale by the number of observed sideband events with
        # the pdf’s fraction in blind vs sideband ranges.
        i_blind = bkg_pdf.model.createIntegral(
            RooArgSet(x), RooFit.Range("blind_left,blind_right")
        ).getVal()

        i_full = bkg_pdf.model.createIntegral(
            RooArgSet(x), RooFit.Range("left,blind_left,middle,blind_right,right")
        ).getVal()

        frac_blind = i_blind / i_full if i_full > 0 else 0.0
        n_est = frac_blind * data_full.numEntries()
        print(
            f"\nBackground estimate in blind window (from fitted pdf): "
            f"N_bkg(blind) ≈ {n_est:.1f} (fraction={frac_blind:.4f})"
        )
