from collections import namedtuple
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import ArrayLike

# Define result type
FitResult = namedtuple("FitResult", ["coeffs", "cov", "chi2", "dof", "y0", "sy0"])


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

    Parameters
    ----------
    x, y, sigma_y : ArrayLike
        Data points and uncertainties.
    x0 : float
        Extrapolation point.
    output_pdf : str, default="fit_result.pdf"
        Path of the PDF file to save the plot.

    Returns
    -------
    FitResult : namedtuple
        Fields:
        - coeffs : np.ndarray, shape (3,), fitted coefficients [a, b, c].
        - cov    : np.ndarray, shape (3,3), covariance matrix of coefficients.
        - chi2   : float, chi-square of the fit.
        - dof    : int, degrees of freedom.
        - y0     : float, predicted value at x0.
        - sy0    : float, uncertainty on prediction at x0.
    """

    # --- Prepare arrays ---
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sy = np.asarray(sigma_y, dtype=float)

    if np.any(sy <= 0):
        raise ValueError("All sigma_y must be > 0.")

    # --- Weighted fit ---
    X = np.column_stack([x**2, x, np.ones_like(x)])
    w = 1.0 / sy**2
    WX = X * w[:, None]
    XT_W_X = X.T @ WX
    XT_W_y = X.T @ (w * y)

    coeffs = np.linalg.solve(XT_W_X, XT_W_y)

    r = y - X @ coeffs
    chi2 = float(np.sum((r / sy) ** 2))
    dof = max(len(x) - 3, 1)
    cov = np.linalg.inv(XT_W_X) * (chi2 / dof)

    # --- Prediction at x0 ---
    v = np.array([x0**2, x0, 1.0])
    y0 = float(v @ coeffs)
    var_y0 = float(v @ cov @ v)
    sy0 = float(np.sqrt(var_y0)) if var_y0 > 0 else 0.0

    # --- Plot ---
    fig, ax = plt.subplots()
    ax.errorbar(
        x, y, yerr=sy, fmt="o", markersize=3, capsize=3, label="CR Normalization"
    )

    xx = np.linspace(min(x) - 0.5, max(x) + 1.0, 300)
    yy = coeffs[0] * xx**2 + coeffs[1] * xx + coeffs[2]
    ax.plot(xx, yy, "r-", label="Quadratic fit")

    ax.errorbar(
        [x0],
        [y0],
        yerr=[sy0],
        fmt="o",
        color="green",
        capsize=3,
        markersize=3,
        label=f"Extrapolation for SR (x={x0}): {y0:.2f} +/- {sy0:.2f}",
    )
    # Draw vertical lines if provided
    if x_lines is not None:
        for xv in x_lines:
            ax.axvline(x=xv, color="blue", linestyle="--", alpha=0.5)

    # Annotate data points if labels are provided
    if point_labels is not None:
        if len(point_labels) != len(x):
            raise ValueError("Length of point_labels must match length of x.")
        for xi, yi, lbl in zip(x, y, point_labels):
            ax.text(
                xi,
                yi + 0.07 * (max(y) - min(y)),  # small offset above point
                lbl,
                ha="center",
                va="bottom",
                fontsize=9,
                rotation=0,
            )

    ax.set_xlabel("CR Midpoint")
    ax.set_ylabel("Normalization")
    ax.legend()
    ax.grid(True, ls="--", alpha=0.5)

    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close("all")

    return FitResult(coeffs, cov, chi2, dof, y0, sy0)


# --- Example usage ---
if __name__ == "__main__":
    x = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
    y = np.array([1.1, 1.5, 2.2, 3.2, 4.1])
    sy = np.array([0.10, 0.12, 0.10, 0.15, 0.20])

    result = fit_and_plot(x, y, sy, x0=2.5, output_pdf="quadratic_fit.pdf")

    print("Coefficients [a, b, c]:", result.coeffs)
    print("Covariance matrix:\n", result.cov)
    print(f"chi2/dof = {result.chi2:.2f}/{result.dof} = {result.chi2 / result.dof:.3f}")
    print(f"Prediction at x0=2.5: {result.y0:.4f} ± {result.sy0:.4f}")
