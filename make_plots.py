from itertools import cycle

from ROOT import (
    RooArgSet,  # type: ignore
    RooFit,  # type: ignore
    TBox,  # type: ignore
    TCanvas,  # type: ignore
    TLegend,  # type: ignore
    gPad,  # type: ignore
    kDashed,  # type: ignore
    kGray,  # type: ignore
)

from bkg_model import BkgModel, BkgPdfFamily

COLORS = cycle([1, 2, 3, 4, 5, 6, 7, 8, 9])


def shade_blind_region(bl_lo, bl_hi, ymin, ymax, color=kGray, alpha=0.35):
    """Return a TBox that shades the blinded region."""
    box = TBox(bl_lo, ymin, bl_hi, ymax)
    box.SetFillColorAlpha(color, alpha)
    box.SetLineColor(color)
    return box


def make_plots(
    x,
    data_full,
    data_sb,
    bkg_pdf: BkgModel,
    test_bkg_pdfs: list[BkgModel],
    bkg_pdf_family: BkgPdfFamily,
    left_lower,
    left_upper,
    middle_lower,
    middle_upper,
    right_lower,
    right_upper,
    outprefix,
    nbins,
    start: int,
    winner: int,
) -> str:
    # Named ranges for sidebands
    x.setRange("left", left_lower, left_upper)
    x.setRange("middle", middle_lower, middle_upper)
    x.setRange("right", right_lower, right_upper)

    # Frame
    frame = x.frame(
        RooFit.Bins(nbins),
        RooFit.Title(""),
    )

    # Plot only sideband data points
    data_full.plotOn(frame, RooFit.Name("data_sb"))
    # data_sb.plotOn(frame, RooFit.Name("data_sb"))

    # Draw fitted background on sidebands only (solid), normalized to sidebands
    bkg_pdf.model.plotOn(
        frame,
        RooFit.Range("left,middle,right"),
        RooFit.NormRange("left,middle,right"),
        RooFit.Name("bkg_sidebands"),
    )

    for test_bkg_pdf in test_bkg_pdfs:
        c = next(COLORS)
        test_bkg_pdf.model.plotOn(
            frame,
            RooFit.Range("left,middle,right"),
            RooFit.NormRange("left,middle,right"),
            # RooFit.Name("bkg_sidebands"),
            RooFit.LineColor(c),
            RooFit.LineStyle(kDashed),
        )

        test_bkg_pdf.model.plotOn(
            frame,
            RooFit.NormRange("left,middle,right"),
            RooFit.Name(
                f"test_bkg_{test_bkg_pdf.model.getParameters(RooArgSet()).getSize()}_params"
            ),
            RooFit.LineColor(c),
            RooFit.LineStyle(kDashed),
        )

    frame.SetMaximum(1.5 * frame.GetMaximum())  # leave  headroom

    # Pulls (computed vs the full-range curve). For clarity we compute pulls only where points exist (sidebands).
    pull_hist = frame.pullHist("data_sb", "bkg_sidebands")
    pull_frame = x.frame(RooFit.Title("Pulls (sidebands only)"))
    pull_frame.addPlotable(pull_hist, "P")

    # Canvas with shaded blinded region
    can = TCanvas("c", "c", 820, 920)
    can.Divide(1, 2)

    # Top pad: fit
    can.cd(1)
    gPad.SetPad(0, 0.30, 1, 1)
    frame.Draw()

    # Shade blinded region
    ymax = frame.GetMaximum()
    box_left = shade_blind_region(left_upper, middle_lower, 0, 1.02 * ymax)
    box_left.Draw("same")
    box_right = shade_blind_region(middle_upper, right_lower, 0, 1.02 * ymax)
    box_right.Draw("same")

    # Add a legend
    leg = TLegend(0.4, 0.5, 0.88, 0.88)
    # leg.SetTextAlign(13)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(frame.findObject("data_sb"), "Data", "lep")
    leg.AddEntry(
        frame.findObject("bkg_sidebands"),
        f"True BKG ({bkg_pdf.pdf_family} - {bkg_pdf.n_params} params)",
        "l",
    )

    for idx, test_bkg_pdf in enumerate(test_bkg_pdfs):
        assert test_bkg_pdf.chi_square_res is not None

        start_winner_str = ""
        if idx == start and idx == winner:
            start_winner_str = " (Start - Winner)"
        elif idx == start:
            start_winner_str = " (Start)"
        elif idx == winner:
            start_winner_str = " (Winner)"

        chi2_pval = test_bkg_pdf.chi_square_res.pvalue
        leg.AddEntry(
            frame.findObject(
                f"test_bkg_{test_bkg_pdf.model.getParameters(RooArgSet()).getSize()}_params"
            ),
            f"{bkg_pdf_family}: {test_bkg_pdf.model.getParameters(RooArgSet()).getSize()} params - #chi^{{2}} p-value: {chi2_pval:.3f}{start_winner_str}",
            "l",
        )

    leg.Draw()

    # Plot only sideband data points - again
    data_full.plotOn(frame, RooFit.Name("data_sb"))
    # data_sb.plotOn(frame, RooFit.Name("data_sb"))

    # Bottom pad: pulls
    can.cd(2)
    gPad.SetPad(0, 0.00, 1, 0.30)
    gPad.SetGridy(True)
    pull_frame.GetYaxis().SetNdivisions(505)
    pull_frame.GetYaxis().SetRangeUser(-5, 5)
    pull_frame.Draw()

    plot_file_name = f"plots/{outprefix}.pdf"
    can.SaveAs(plot_file_name)

    return plot_file_name
