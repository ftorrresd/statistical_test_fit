from array import array
from enum import Enum, auto
from itertools import cycle
import os
from typing import Optional, List, Any

from ROOT import (
    RooArgSet,  # type: ignore
    RooBinning,  # type: ignore
    RooFit,  # type: ignore
    TBox,  # type: ignore
    TCanvas,  # type: ignore
    TLegend,  # type: ignore
    gPad,  # type: ignore
    gStyle,  # type: ignore
    kDashed,  # type: ignore
    kGray,  # type: ignore
    kRed,  # type: ignore
    kSolid,  # type: ignore
)
from typing_extensions import ReadOnly

from .bkg_model import BkgModel, BkgPdfFamily

try:
    import sys as _sys

    _cmsstyle_src = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "cmsstyle",
        "src",
    )
    _cmsstyle_src = os.path.abspath(_cmsstyle_src)
    if _cmsstyle_src not in _sys.path:
        _sys.path.insert(0, _cmsstyle_src)
    import cmsstyle as CMS  # noqa: E402

    _has_cmsstyle = True
except ImportError:
    _has_cmsstyle = False


_cmsstyle_labels_done = False


def _init_cmsstyle_labels():
    global _cmsstyle_labels_done
    if not _has_cmsstyle or _cmsstyle_labels_done:
        return
    _cmsstyle_labels_done = True
    CMS.setCMSStyle()
    CMS.SetExtraText("Private work (CMS data)")
    CMS.SetEnergy(13)
    CMS.SetLumi(138)


def _cms_style_setup():
    _init_cmsstyle_labels()


def _cms_add_labels(can):
    if not _has_cmsstyle:
        return
    CMS.UpdatePad(can)
    CMS.CMS_lumi(can, iPosX=0, scaleLumi=0.8)


def _curve_name_for_model(test_bkg_pdf: BkgModel, idx: int) -> str:
    return f"test_bkg_{test_bkg_pdf.pdf_family.value}_{idx}"


def _format_model_label(test_bkg_pdf: BkgModel) -> str:
    label = str(test_bkg_pdf.pdf_family)
    if test_bkg_pdf.scan_order is not None:
        label += f" order {test_bkg_pdf.scan_order}"
    if test_bkg_pdf.n_float_params is not None:
        label += f" ({test_bkg_pdf.n_float_params} floated)"
    return label


def _plot_order(n_candidates: int, winner: Optional[int]) -> list[int]:
    order = [idx for idx in range(n_candidates) if idx != winner]
    if winner is not None and 0 <= winner < n_candidates:
        order.append(winner)
    return order


def get_ROOT_colors(data_type):
    if data_type == DataType.PSEUDO:
        return cycle([1, 2, 3, 5, 6, 7, 8, 9, 4])
    if data_type == DataType.REAL:
        return cycle([4, 1, 2, 3, 6, 7, 8, 9, 5, 4])


def shade_blind_region(bl_lo, bl_hi, ymin, ymax, color=kGray, alpha=0.35):
    """Return a TBox that shades the blinded region."""
    box = TBox(bl_lo, ymin, bl_hi, ymax)
    box.SetFillColorAlpha(color, alpha)
    box.SetLineColor(color)
    return box


def _build_edge_aligned_binning(nbins: int, boundaries: list[float], name: str):
    """Build non-uniform plot binning with exact user-supplied boundaries."""
    if nbins < len(boundaries) - 1:
        return None

    widths = [high - low for low, high in zip(boundaries[:-1], boundaries[1:])]
    total_width = sum(widths)
    ideal = [nbins * width / total_width for width in widths]
    bins_per_interval = [max(1, int(value)) for value in ideal]

    while sum(bins_per_interval) > nbins:
        candidates = [idx for idx, bins in enumerate(bins_per_interval) if bins > 1]
        idx = min(candidates, key=lambda i: ideal[i] - bins_per_interval[i])
        bins_per_interval[idx] -= 1

    deficit = nbins - sum(bins_per_interval)
    remainder_order = sorted(
        range(len(bins_per_interval)),
        key=lambda idx: ideal[idx] - int(ideal[idx]),
        reverse=True,
    )
    for idx in remainder_order[:deficit]:
        bins_per_interval[idx] += 1

    edges = [boundaries[0]]
    for low, high, interval_bins in zip(
        boundaries[:-1], boundaries[1:], bins_per_interval
    ):
        step = (high - low) / interval_bins
        for bin_idx in range(1, interval_bins + 1):
            edge = high if bin_idx == interval_bins else low + step * bin_idx
            edges.append(edge)

    return RooBinning(len(edges) - 1, array("d", edges), name)


def make_plots_1d(
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
    start: Optional[int],
    winner: Optional[int],
) -> str:
    _init_cmsstyle_labels()
    # Frame
    frame = x.frame(
        RooFit.Bins(nbins),
        RooFit.Title(""),
    )

    # Plot only sideband data points
    # data_full.plotOn(frame, RooFit.Name("data_sb"))
    data_sb.plotOn(frame, RooFit.Name("data_sb"))

    # Draw fitted background on sidebands only (solid), normalized to sidebands
    bkg_pdf.model.plotOn(
        frame,
        # RooFit.Range("left,middle,right"),
        RooFit.NormRange("left,middle,right"),
        RooFit.Name("bkg_sidebands"),
    )

    ROOT_COLORS = get_ROOT_colors(DataType.PSEUDO)
    for idx, test_bkg_pdf in enumerate(test_bkg_pdfs):
        c = next(ROOT_COLORS)
        test_bkg_pdf.model.plotOn(
            frame,
            RooFit.NormRange("left,middle,right"),
            RooFit.Name(_curve_name_for_model(test_bkg_pdf, idx)),
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
    can.SetLeftMargin(0.18)
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
    leg = TLegend(0.26, 0.58, 0.54, 0.90)
    # leg.SetTextAlign(13)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.022)
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
            frame.findObject(_curve_name_for_model(test_bkg_pdf, idx)),
            f"{_format_model_label(test_bkg_pdf)} - #chi^{{2}} p-value: {chi2_pval:.3f}{start_winner_str}",
            "l",
        )

    leg.Draw()

    # Bottom pad: pulls
    can.cd(2)
    gPad.SetPad(0, 0.00, 1, 0.30)
    gPad.SetGridy(True)
    pull_frame.GetYaxis().SetNdivisions(505)
    pull_frame.GetYaxis().SetRangeUser(-5, 5)
    pull_frame.Draw()

    plot_file_name = f"plots/fit_1d/{outprefix}.pdf"

    if _has_cmsstyle:
        can.cd(1)
        _cms_add_labels(can)
    can.SaveAs(plot_file_name)

    return plot_file_name


class ProjDim(Enum):
    X = "x"
    Y = "y"


class DataType(Enum):
    REAL = auto()
    PSEUDO = auto()


def make_plots_2d(
    proj_dim: ProjDim,
    x,
    y,
    data_full,
    data_sb,
    bkg_pdf: Optional[BkgModel],
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
    data_type: DataType,
    start: Optional[int],
    winner: Optional[int],
    components: Optional[List[Any]] = None,
) -> str:
    _init_cmsstyle_labels()
    data_postfix = ""
    if data_type == DataType.REAL:
        data_postfix = "_data"

    # Frame
    plot_binning = None
    if proj_dim == ProjDim.Y:
        if data_type == DataType.REAL:
            plot_binning = _build_edge_aligned_binning(
                nbins,
                [
                    left_lower,
                    left_upper,
                    middle_lower,
                    middle_upper,
                    right_lower,
                    right_upper,
                ],
                "mumugamma_plot_binning",
            )
        frame = y.frame(
            RooFit.Bins(nbins),
            RooFit.Title(""),
        )
    else:
        frame = x.frame(
            RooFit.Bins(nbins),
            RooFit.Title(""),
        )

    frame.SetTitle("")

    # Plot only sideband data points
    data_plot_args = [RooFit.Name("data_sb")]
    if plot_binning is not None:
        data_plot_args.append(RooFit.Binning(plot_binning))
    data_sb.plotOn(frame, *data_plot_args)
    if data_full is not None and proj_dim == ProjDim.Y:
        full_data_plot_args = [
            RooFit.Name("data_full"),
            RooFit.MarkerStyle(20),
            RooFit.MarkerColor(kRed),
        ]
        if plot_binning is not None:
            full_data_plot_args.append(RooFit.Binning(plot_binning))
        data_full.plotOn(frame, *full_data_plot_args)

    # Draw fitted background on sidebands only (solid), normalized to sidebands
    if bkg_pdf is not None:
        if proj_dim == ProjDim.Y:
            bkg_pdf.model.plotOn(
                frame,
                # RooFit.Range("LEFT,MIDDLE,RIGHT"),
                RooFit.NormRange("LEFT,MIDDLE,RIGHT"),
                RooFit.Name("bkg_sidebands"),
            )
        else:
            bkg_pdf.model.plotOn(
                frame,
                RooFit.Name("bkg_sidebands"),
            )

    ROOT_COLORS = get_ROOT_colors(data_type)
    color_by_idx = {idx: next(ROOT_COLORS) for idx in range(len(test_bkg_pdfs))}
    for idx in _plot_order(len(test_bkg_pdfs), winner):
        test_bkg_pdf = test_bkg_pdfs[idx]
        line_style = kSolid if idx == winner else kDashed
        if proj_dim == ProjDim.Y:
            test_bkg_pdf.model.plotOn(
                frame,
                RooFit.NormRange("LEFT,MIDDLE,RIGHT"),
                RooFit.Name(_curve_name_for_model(test_bkg_pdf, idx)),
                RooFit.LineColor(color_by_idx[idx]),
                RooFit.LineStyle(line_style),
            )
        else:
            test_bkg_pdf.model.plotOn(
                frame,
                RooFit.Name(_curve_name_for_model(test_bkg_pdf, idx)),
                RooFit.LineColor(color_by_idx[idx]),
                RooFit.LineStyle(line_style),
            )

    leg = TLegend(0.26, 0.58, 0.54, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.022)
    if components is not None:
        for n, c in components:
            test_bkg_pdfs[0].model.plotOn(
                frame,
                RooFit.Components(c.GetName()),
                RooFit.LineColor(kRed - 7),
                RooFit.FillColor(kRed - 7),
                RooFit.FillStyle(1001),
                RooFit.Name(c.GetName()),
                RooFit.DrawOption("F"),
            )

            leg.AddEntry(frame.findObject(c.GetName()), c.getTitle().Data())

    # if proj_dim == ProjDim.X:
    frame.SetMaximum(2.0 * frame.GetMaximum())  # leave  headroom

    # # Pulls (computed vs the full-range curve). For clarity we compute pulls only where points exist (sidebands).
    # pull_hist = frame.pullHist("data_sb", "bkg_sidebands")
    # pull_frame = x.frame(RooFit.Title("Pulls (sidebands only)"))
    # pull_frame.addPlotable(pull_hist, "P")

    # Canvas with shaded blinded region
    can = TCanvas("c", "c", 820, 920)
    can.SetLeftMargin(0.18)
    # can.Divide(1, 2)
    # Top pad: fit
    # can.cd(1)
    # gPad.SetPad(0, 0.30, 1, 1)
    frame.Draw()

    # Shade blinded region
    ymax = frame.GetMaximum()
    box_left = shade_blind_region(left_upper, middle_lower, 0, 1.02 * ymax)
    box_left.Draw("same")
    box_right = shade_blind_region(middle_upper, right_lower, 0, 1.02 * ymax)
    box_right.Draw("same")

    # Add a legend
    # leg.SetTextAlign(13)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    if data_type == DataType.REAL:
        leg.AddEntry(frame.findObject("data_sb"), "Data", "lep")
    if data_type == DataType.PSEUDO:
        leg.AddEntry(frame.findObject("data_sb"), "Sideband pseudodata", "lep")
        if data_full is not None and proj_dim == ProjDim.Y:
            leg.AddEntry(
                frame.findObject("data_full"),
                "Full pseudodata (incl. blinded region)",
                "lep",
            )
    if bkg_pdf is not None:
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
            frame.findObject(_curve_name_for_model(test_bkg_pdf, idx)),
            f"{_format_model_label(test_bkg_pdf)} - #chi^{{2}} p-value: {chi2_pval:.3f}{start_winner_str}",
            "l",
        )

    leg.Draw()

    # Plot only sideband data points - again
    data_sb.plotOn(frame, *data_plot_args)

    # # Bottom pad: pulls
    # can.cd(2)
    # gPad.SetPad(0, 0.00, 1, 0.30)
    # gPad.SetGridy(True)
    # pull_frame.GetYaxis().SetNdivisions(505)
    # pull_frame.GetYaxis().SetRangeUser(-5, 5)
    # pull_frame.Draw()

    if proj_dim == ProjDim.X:
        plot_file_name = f"plots/fit_2d{data_postfix}/mumu_{outprefix}.pdf"
    else:
        plot_file_name = f"plots/fit_2d{data_postfix}/mumugamma_{outprefix}.pdf"
    if _has_cmsstyle:
        _cms_add_labels(can)
    can.SaveAs(plot_file_name)

    return plot_file_name
