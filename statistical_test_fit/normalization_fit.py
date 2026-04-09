from array import array
from collections import namedtuple
from typing import List, Optional

import numpy as np
from numpy.typing import ArrayLike
from ROOT import (  # type: ignore
    TCanvas,  # type: ignore
    TGraph,  # type: ignore
    TGraphErrors,  # type: ignore
    TLegend,  # type: ignore
    TLatex,  # type: ignore
    TLine,  # type: ignore
    kBlue,  # type: ignore
    kGreen,  # type: ignore
    kRed,  # type: ignore
)

# Define result type
FitResult = namedtuple("FitResult", ["coeffs", "cov", "chi2", "dof", "y0", "sy0"])


def _to_double_array(values) -> array:
    return array("d", [float(v) for v in values])


def fit_and_plot(
    x: ArrayLike,
    y: ArrayLike,
    sigma_y: ArrayLike,
    x0: float,
    output_pdf: str = "fit_result.pdf",
    x_lines: Optional[List[float]] = None,
    point_labels: Optional[List[str]] = None,
) -> FitResult:
    """
    Perform weighted quadratic fit, predict at x0 with uncertainty,
    plot everything, save to PDF, and close all figures.
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sy = np.asarray(sigma_y, dtype=float)

    if np.any(sy <= 0):
        raise ValueError("All sigma_y must be > 0.")

    X = np.column_stack([x**2, x, np.ones_like(x)])
    w = 1.0 / sy**2
    WX = X * w[:, None]
    XT_W_X = X.T @ WX
    XT_W_y = X.T @ (w * y)

    coeffs = np.linalg.solve(XT_W_X, XT_W_y)

    residuals = y - X @ coeffs
    chi2 = float(np.sum((residuals / sy) ** 2))
    dof = max(len(x) - 3, 1)
    cov = np.linalg.inv(XT_W_X) * (chi2 / dof)

    v = np.array([x0**2, x0, 1.0])
    y0 = float(v @ coeffs)
    var_y0 = float(v @ cov @ v)
    sy0 = float(np.sqrt(var_y0)) if var_y0 > 0 else 0.0

    xx = np.linspace(min(x) - 0.5, max(x) + 1.0, 300)
    yy = coeffs[0] * xx**2 + coeffs[1] * xx + coeffs[2]

    canvas = TCanvas("normalization_fit", "normalization_fit", 900, 700)
    canvas.SetLeftMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.06)

    graph = TGraphErrors(
        len(x),
        _to_double_array(x),
        _to_double_array(y),
        _to_double_array(np.zeros_like(x)),
        _to_double_array(sy),
    )
    graph.SetTitle("")
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.0)
    graph.SetLineWidth(2)
    graph.GetXaxis().SetTitle("CR Midpoint")
    graph.GetYaxis().SetTitle("Normalization")

    y_values = np.concatenate([y - sy, y + sy, yy, np.array([y0 - sy0, y0 + sy0])])
    y_min = float(np.min(y_values))
    y_max = float(np.max(y_values))
    y_span = y_max - y_min
    if y_span <= 0:
        y_span = max(abs(y_max), 1.0)
    graph.SetMinimum(y_min - 0.15 * y_span)
    graph.SetMaximum(y_max + 0.20 * y_span)
    graph.Draw("AP")

    fit_graph = TGraph(len(xx), _to_double_array(xx), _to_double_array(yy))
    fit_graph.SetLineColor(kRed)
    fit_graph.SetLineWidth(2)
    fit_graph.Draw("L SAME")

    pred_graph = TGraphErrors(
        1,
        _to_double_array([x0]),
        _to_double_array([y0]),
        _to_double_array([0.0]),
        _to_double_array([sy0]),
    )
    pred_graph.SetMarkerStyle(21)
    pred_graph.SetMarkerSize(1.1)
    pred_graph.SetMarkerColor(kGreen + 2)
    pred_graph.SetLineColor(kGreen + 2)
    pred_graph.SetLineWidth(2)
    pred_graph.Draw("P SAME")

    line_bottom = graph.GetMinimum()
    line_top = graph.GetMaximum()
    lines = []
    if x_lines is not None:
        for xv in x_lines:
            line = TLine(float(xv), line_bottom, float(xv), line_top)
            line.SetLineColor(kBlue)
            line.SetLineStyle(2)
            line.Draw("SAME")
            lines.append(line)

    labels = []
    if point_labels is not None:
        if len(point_labels) != len(x):
            raise ValueError("Length of point_labels must match length of x.")
        label_offset = 0.05 * (line_top - line_bottom)
        for xi, yi, lbl in zip(x, y, point_labels):
            text = TLatex(float(xi), float(yi + label_offset), lbl)
            text.SetTextAlign(22)
            text.SetTextSize(0.03)
            text.Draw()
            labels.append(text)

    legend = TLegend(0.52, 0.70, 0.90, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(graph, "CR Normalization", "lep")
    legend.AddEntry(fit_graph, "Quadratic fit", "l")
    legend.AddEntry(
        pred_graph,
        f"Extrapolation for SR (x={x0:.1f}): {y0:.2f} +/- {sy0:.2f}",
        "lep",
    )
    legend.Draw()

    canvas.SetGrid(1, 1)
    canvas.SaveAs(output_pdf)
    canvas.Close()

    return FitResult(coeffs, cov, chi2, dof, y0, sy0)


if __name__ == "__main__":
    x = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
    y = np.array([1.1, 1.5, 2.2, 3.2, 4.1])
    sy = np.array([0.10, 0.12, 0.10, 0.15, 0.20])

    result = fit_and_plot(x, y, sy, x0=2.5, output_pdf="quadratic_fit.pdf")

    print("Coefficients [a, b, c]:", result.coeffs)
    print("Covariance matrix:\n", result.cov)
    print(f"chi2/dof = {result.chi2:.2f}/{result.dof} = {result.chi2 / result.dof:.3f}")
    print(f"Prediction at x0=2.5: {result.y0:.4f} ± {result.sy0:.4f}")
