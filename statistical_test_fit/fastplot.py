# REF: https://www.desy.de/~swehle/pyroofit/_modules/pyroofit/plotting.html#fast_plot

from __future__ import print_function

import math
from enum import Enum

import ROOT

DEFAULT_PALETTE = [
    1,
    ROOT.kRed - 7,
    ROOT.kAzure + 5,
    ROOT.kGreen - 2,
    ROOT.kMagenta + 1,
    ROOT.kYellow,
]
DEFAULT_STYLES = [0, 1001, 3004, 3005, 3009, 3006]


class ResidualMode(Enum):
    PULLS = "pulls"
    DIFF = "diff"
    RELATIVE_DIFF = "relative_diff"


def set_root_style(font_scale=1.0, label_scale=1.0):
    """Setting a general style that one can look at plots without getting eye-cancer.

    Args:
        font_scale (float): Scale of the fonts
        label_scale (float): Scale of the labels

    Todo:
        * Absolute font size

    """
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetLabelSize(0.04 * label_scale, "xy")
    ROOT.gStyle.SetLabelOffset(0.006, "y")
    ROOT.gStyle.SetTitleSize(0.06 * font_scale, "xy")
    ROOT.gStyle.SetTitleOffset(0.9, "x")
    ROOT.gStyle.SetTitleOffset(1.15, "y")
    ROOT.gStyle.SetNdivisions(505, "x")

    ROOT.gStyle.SetPadLeftMargin(0.14)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadTopMargin(0.05)

    ROOT.gStyle.SetFillColor(0)
    ROOT.gStyle.SetMarkerSize(0.8)
    ROOT.gStyle.SetLineColor(ROOT.kBlack)
    ROOT.gStyle.SetLineWidth(1)

    ROOT.gStyle.SetLegendBorderSize(0)


def round_to_1(x):
    from math import floor, log10

    return round(x, -int(floor(log10(abs(x)))))


def fastplot(
    model,
    data,
    observable,
    filename,
    components=None,
    nbins=None,
    extra_info=None,
    size=1280,
    average=True,
    pi_label=False,
    font_scale=1.0,
    label_scale=1.0,
    legend=False,
    extra_text=None,
    round_bins=5,
    tick_len=30,
    color_cycle=DEFAULT_PALETTE,
    fill_cycle=DEFAULT_STYLES,
    lw=2,
    line_shade=0,
    data_range=None,
    plot_range=None,
    y_min=0,
    y_max=-999,
    is_data=True,
    model_legend_name="Fit",
    show_residuals=True,
    residual_mode=ResidualMode.RELATIVE_DIFF,
):
    """Generic plot function

    Args:
        model (RooAbsPDF):
            Fit model to be drawn
        data (RooDataSet):
            Dataset to be plotted
        observable (RooAbsVar):
            Observable to be drawn
        filename (str):
            Name of the output file. Suffix determines file type
        components (list of tuples):
            Normalisation and ROOT.RooAbsPDF to be drawn searately
        nbins (int):
            Number of bins
        extra_info (list or TPaveText):
        lw (int):
            Width of the line of the total fit model
        size (int):
            Plot size in pixels
        average (bool):
            Average bin content for calculating the pull distribution, if false, take central value
        pi_label (bool):
            Calculate the bin count in radians
        font_scale (float):
            Set relative font scale
        label_scale (float):
            Set relative lable scale
        color_cycle (list of ROOT.TColor):
            Overwrite for the default color cycle
        fill_cycle (list of ROOT.TAttrFill):
            Overwrite for the default fill cycle
        line_shade (int):
            Integer to add to the color cycle for the fill color
        legend (list):
            Vector with four coordinates for the TLegend position
        extra_text (list of ROOT.TPaveText or ROOT.TPaveText):
            Extra text to be drawn on the plot
        round_bins (int) :
            magic to for automatically choosing the bin numbers
        tick_len (int) :
            Sets the length of the bins, EQUALLY (yes root. this is possible.), choose between 0-100
        data_range (str) :
            Range to plot
        show_residuals (bool):
            Draw a residual panel below the main plot
        residual_mode (ResidualMode):
            Residual strategy for the lower panel

    Todo:
        * Change or remove extra_info
    """
    set_root_style(font_scale, label_scale)
    if nbins is None:
        nbins = 60

    if not isinstance(residual_mode, ResidualMode):
        raise TypeError("residual_mode must be a ResidualMode enum value")

    # Use a slightly finer display binning for plotting only.
    plot_nbins = max(nbins + 10, int(round(1.5 * nbins)))

    x_min = float(plot_range[0]) if plot_range is not None else observable.getMin()
    x_max = float(plot_range[1]) if plot_range is not None else observable.getMax()
    bin_half_width = 0.5 * (x_max - x_min) / float(plot_nbins)

    def _build_frame(title="Fit Result"):
        frame_args = [ROOT.RooFit.Title(title), ROOT.RooFit.Bins(plot_nbins)]
        if plot_range is not None:
            frame_args.append(
                ROOT.RooFit.Range(float(plot_range[0]), float(plot_range[1]))
            )
        return observable.frame(*frame_args)

    def _plot_data(target_frame, name):
        data.plotOn(
            target_frame,
            ROOT.RooFit.Name(name),
            ROOT.RooFit.Binning(plot_nbins),
            ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),
        )

    def _plot_model(target_frame, name):
        if data_range is not None:
            model.plotOn(
                target_frame,
                ROOT.RooFit.Name(name),
                ROOT.RooFit.LineColor(color_cycle[0]),
                ROOT.RooFit.Range("full"),
                ROOT.RooFit.NormRange(data_range),
                ROOT.RooFit.ProjectionRange("full"),
            )
        else:
            model.plotOn(
                target_frame,
                ROOT.RooFit.Name(name),
                ROOT.RooFit.LineColor(color_cycle[0]),
            )

    def _max_abs_graph_y(graph):
        max_abs = 0.0
        y_values = graph.GetY()
        for i in range(graph.GetN()):
            y = float(y_values[i])
            if not math.isfinite(y):
                continue
            max_abs = max(max_abs, abs(y))
        return max_abs

    def _build_residual_hist(residual_source):
        if residual_mode == ResidualMode.PULLS:
            return residual_source.pullHist(
                "DataResidual",
                "ModelResidual",
                average,
            )

        residual_hist = residual_source.residHist(
            "DataResidual",
            "ModelResidual",
            False,
            average,
        )
        if residual_mode == ResidualMode.DIFF:
            return residual_hist

        model_curve = residual_source.getCurve("ModelResidual")
        if model_curve is None:
            raise RuntimeError("Could not find model curve for residual panel")

        curve_scale = _max_abs_graph_y(model_curve)
        min_fit_abs = max(1e-12, 1e-9 * curve_scale)
        for i in range(residual_hist.GetN() - 1, -1, -1):
            x = float(residual_hist.GetX()[i])
            residual_y = float(residual_hist.GetY()[i])
            fit_y = float(model_curve.Eval(x))
            if (
                not math.isfinite(residual_y)
                or not math.isfinite(fit_y)
                or abs(fit_y) <= min_fit_abs
            ):
                residual_hist.RemovePoint(i)
                continue

            residual_hist.SetPoint(i, x, residual_y / fit_y)
            residual_hist.SetPointEYlow(
                i,
                residual_hist.GetErrorYlow(i) / abs(fit_y),
            )
            residual_hist.SetPointEYhigh(
                i,
                residual_hist.GetErrorYhigh(i) / abs(fit_y),
            )

        return residual_hist

    def _residual_axis_title():
        if residual_mode == ResidualMode.PULLS:
            return "(Data - Fit) / #sigma"
        if residual_mode == ResidualMode.DIFF:
            return "Data - Fit"
        return "(Data - Fit) / Fit"

    def _default_residual_extent():
        if residual_mode == ResidualMode.PULLS:
            return 3.5
        return 1.0

    # nbins = get_optimal_bin_size(data.numEntries(), round_bins) if nbins is None else nbins
    # if isinstance(data, ROOT.RooDataHist):
    #     nbins = observable.getBins()
    frame = _build_frame()

    if isinstance(legend, list):
        assert len(legend) == 4, "Please provide four coordinates for the legend"
        leg = ROOT.TLegend(*legend)
    else:
        leg = ROOT.TLegend(0.7, 0.78, 0.93, 0.92)

    _plot_data(frame, "Data")

    data_legend_name = "Data"
    if not is_data:
        data_legend_name = "MC"

    leg.AddEntry(frame.findObject("Data"), data_legend_name, "LEP")

    _plot_model(frame, "Model")

    leg.AddEntry(frame.findObject("Model"), model_legend_name, "L")

    if components is not None:
        n_col = 1
        for c, component_label in components:
            if data_range != None:  # here it should be !=
                model.plotOn(
                    frame,
                    ROOT.RooFit.Components(c.GetName()),
                    ROOT.RooFit.LineColor(color_cycle[n_col] + line_shade),
                    # ROOT.RooFit.Normalization(ni, 2),
                    ROOT.RooFit.FillColor(color_cycle[n_col]),
                    ROOT.RooFit.FillStyle(fill_cycle[n_col]),
                    ROOT.RooFit.Name(c.GetName()),
                    ROOT.RooFit.DrawOption("F"),
                    ROOT.RooFit.NormRange(data_range),
                    #                     ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange(data_range), ROOT.RooFit.ProjectionRange("full")
                )
            else:
                model.plotOn(
                    frame,
                    ROOT.RooFit.Components(c.GetName()),
                    ROOT.RooFit.LineColor(color_cycle[n_col] + line_shade),
                    # ROOT.RooFit.Normalization(ni, 2),
                    ROOT.RooFit.FillColor(color_cycle[n_col]),
                    ROOT.RooFit.FillStyle(fill_cycle[n_col]),
                    ROOT.RooFit.Name(c.GetName()),
                    ROOT.RooFit.DrawOption("F"),
                )

            if not isinstance(component_label, str):
                component_label = c.getTitle().Data()

            leg.AddEntry(frame.findObject(c.GetName()), component_label)

            if data_range != None:  # here it should be !=
                model.plotOn(
                    frame,
                    ROOT.RooFit.Components(c.GetName()),
                    ROOT.RooFit.LineColor(color_cycle[n_col] + line_shade),
                    #  ROOT.RooFit.Normalization(ni, 2),
                    ROOT.RooFit.FillColor(color_cycle[n_col]),
                    ROOT.RooFit.LineWidth(lw),
                    ROOT.RooFit.NormRange(data_range),
                    #                       ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange(data_range), ROOT.RooFit.ProjectionRange("full"),
                    # ROOT.RooFit.DrawOption("F"),
                )  # ROOT.RooFit.DrawOption("F")) #4050
            else:
                model.plotOn(
                    frame,
                    ROOT.RooFit.Components(c.GetName()),
                    ROOT.RooFit.LineColor(color_cycle[n_col] + line_shade),
                    #  ROOT.RooFit.Normalization(ni, 2),
                    ROOT.RooFit.FillColor(color_cycle[n_col]),
                    ROOT.RooFit.LineWidth(lw),
                    # ROOT.RooFit.DrawOption("BC"),
                )  # ROOT.RooFit.DrawOption("F")) #4050
            n_col += 1

    _plot_model(frame, "Model")

    _plot_data(frame, "Data")

    # Print num of entries
    print("sumEntries for {0}: {1}".format(filename, data.sumEntries()))

    residual_frame = None
    zero_line = None
    if show_residuals:
        residual_source = _build_frame()
        _plot_data(residual_source, "DataResidual")
        _plot_model(residual_source, "ModelResidual")

        residual_hist = _build_residual_hist(residual_source)
        for i in range(residual_hist.GetN()):
            residual_hist.SetPointEXlow(i, bin_half_width)
            residual_hist.SetPointEXhigh(i, bin_half_width)

        residual_frame = _build_frame("")
        residual_frame.addPlotable(residual_hist, "PE1")
        residual_frame.SetTitle("")
        residual_frame.SetYTitle(_residual_axis_title())
        residual_frame.GetYaxis().CenterTitle()
        residual_frame.GetYaxis().SetNdivisions(505)

        max_abs_residual = _max_abs_graph_y(residual_hist)
        residual_extent = (
            max(_default_residual_extent(), 1.15 * max_abs_residual)
            if max_abs_residual > 0
            else _default_residual_extent()
        )
        residual_frame.SetMinimum(-residual_extent)
        residual_frame.SetMaximum(residual_extent)

        residual_frame.GetXaxis().SetTitleSize(0.12 * font_scale)
        residual_frame.GetXaxis().SetLabelSize(0.10 * label_scale)
        residual_frame.GetXaxis().SetTitleOffset(1.0)
        residual_frame.GetYaxis().SetTitleSize(0.10 * font_scale)
        residual_frame.GetYaxis().SetLabelSize(0.09 * label_scale)
        residual_frame.GetYaxis().SetTitleOffset(0.6)

        zero_line = ROOT.TLine(x_min, 0.0, x_max, 0.0)
        zero_line.SetLineColor(ROOT.kGray + 2)
        zero_line.SetLineStyle(2)
        zero_line.SetLineWidth(1)

    # Create Canvas
    canvas = ROOT.TCanvas("plot", "plot", size, size)
    if show_residuals:
        top_pad = ROOT.TPad("top_pad", "top_pad", 0.0, 0.28, 1.0, 1.0)
        bottom_pad = ROOT.TPad("bottom_pad", "bottom_pad", 0.0, 0.0, 1.0, 0.28)

        top_pad.SetBottomMargin(0.03)
        top_pad.SetLeftMargin(0.14)
        top_pad.SetRightMargin(0.05)
        top_pad.SetTopMargin(0.05)
        top_pad.SetTicks(1, 1)

        bottom_pad.SetTopMargin(0.03)
        bottom_pad.SetBottomMargin(0.35)
        bottom_pad.SetLeftMargin(0.14)
        bottom_pad.SetRightMargin(0.05)
        bottom_pad.SetTicks(1, 1)
        bottom_pad.SetGridy(1)

        top_pad.Draw()
        bottom_pad.Draw()
        main_pad = top_pad
    else:
        canvas.SetRightMargin(0.05)
        canvas.SetTicks(1, 1)
        main_pad = canvas

    # Pi label because of...
    if pi_label:
        pifactor = 1 if observable.getMax() > 1.9 else 2
        ylabel = "Events / ( %.2f #pi rad )" % (1.0 / float(pifactor * plot_nbins))
        frame.SetYTitle(ylabel)
    else:
        obs_range = round_to_1(x_max - x_min)  # stupid overflow artefacts
        # obs_range = (observable.getMax() - observable.getMin())  # stupid overflow artefacts
        div = round(plot_nbins / obs_range)
        # div = (float(nbins)/obs_range)
        # print(div,obs_range,numbins)
        unit = observable.getUnit()
        if unit is not None or unit != "":
            ylabel = "Events / ( %s / %d )" % (observable.getUnit(), div)
            frame.SetYTitle(ylabel)

    if show_residuals:
        frame.GetXaxis().SetLabelSize(0)
        frame.GetXaxis().SetTitleSize(0)
        frame.GetXaxis().SetTickLength(0)

    # Draw All The Stuff
    main_pad.cd()
    if y_max != -999:
        yaxis = frame.GetYaxis()
        yaxis.SetRangeUser(y_min, y_max)
    else:
        frame.SetMinimum(y_min)
        headroom_factor = 1.35 if legend is not False else 1.20
        frame.SetMaximum(headroom_factor * frame.GetMaximum())
    frame.Draw()
    if legend is not False:
        leg.Draw("same")

    if show_residuals:
        bottom_pad.cd()
        residual_frame.Draw()
        zero_line.Draw("same")

    if extra_text is not None:
        main_pad.cd()
        if isinstance(extra_text, ROOT.TPaveText):
            extra_info.Draw("Same")
        if isinstance(extra_text, list):
            for txt in extra_text:
                assert isinstance(txt, ROOT.TPaveText), (
                    "Please provide extra_txt with a list or ROOT.TPaveText"
                )
                txt.Draw("Same")

    if extra_info is not None:
        main_pad.cd()
        if isinstance(extra_info, ROOT.TPaveText):
            extra_info.Draw("Same")
        else:
            assert isinstance(extra_info, list), (
                "Please provide extra_info with a list or ROOT.TPaveText"
            )
            box = ROOT.TPaveText(0.2, 0.75, 0.4, 0.9, "NDC")
            box.SetFillColor(10)
            box.SetBorderSize(0)
            box.SetTextAlign(12)
            box.SetTextSize(0.04)
            box.SetFillStyle(1001)
            box.SetFillColor(10)
            for info in extra_info:
                try:
                    if not isinstance(info, list):
                        if isinstance(info, ROOT.TPaveText):
                            info.Draw("same")
                        else:
                            info = [info]
                    if len(info) == 1:
                        box.AddText(info[0])
                    elif len(info) == 3:
                        box.AddText(info[0] + " = %.2f #pm %.2f" % (info[1], info[2]))
                    else:
                        print("Could not add to legend ", info)
                except IndexError:
                    print("Something went wrong in plotting")

            box.Draw("same")

    canvas.SaveAs(filename)
