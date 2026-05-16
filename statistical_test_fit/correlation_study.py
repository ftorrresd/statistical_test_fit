from __future__ import annotations

from argparse import Namespace
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path

from .mass_ranges import (
    UPSILON_MASS_LOWER,
    UPSILON_MASS_UPPER,
    get_signal_boson_plot_range,
)
from .parallel_utils import ParallelJob, run_parallel_jobs
from .signal_modeling import SIGNAL_SAMPLES, SignalSample


PLOT_DIR = "plots/correlation_study"
RESULTS_JSON = f"{PLOT_DIR}/correlation_study.json"
SUMMARY_PLOT = f"{PLOT_DIR}/correlation_study_summary.pdf"
SIGNAL_PROCESSES = ("H", "Z")
DEFAULT_M_MUMU_LOWER = 9.0
DEFAULT_M_MUMU_UPPER = 10.6


@dataclass(frozen=True)
class MumugammaWindow:
    index: int
    low: float
    high: float
    include_high: bool

    @property
    def label(self) -> str:
        return f"window_{self.index:02d}_{self.low:.3f}_{self.high:.3f}"

    @property
    def cut(self) -> str:
        upper_operator = "<=" if self.include_high else "<"
        return (
            f"upsilon_mass >= {self.low:.12g} && "
            f"upsilon_mass {upper_operator} {self.high:.12g}"
        )


@dataclass(frozen=True)
class MumugammaWindowFitJob:
    process: str
    samples: tuple[SignalSample, ...]
    window: MumugammaWindow
    nbins: int
    plot_dir: str

    @property
    def key(self) -> str:
        return f"{self.process}_{self.window.label}"

    @property
    def label(self) -> str:
        return (
            f"{self.process} m_mumu [{self.window.low:.3f}, "
            f"{self.window.high:.3f}{']' if self.window.include_high else ')'}"
        )


def _build_windows(n_windows: int, lower: float, upper: float) -> list[MumugammaWindow]:
    if n_windows <= 0:
        raise ValueError("Number of windows must be positive")
    if lower >= upper:
        raise ValueError("m_mumu lower bound must be smaller than upper bound")

    width = (upper - lower) / n_windows
    return [
        MumugammaWindow(
            index=idx,
            low=lower + idx * width,
            high=lower + (idx + 1) * width,
            include_high=idx == n_windows - 1,
        )
        for idx in range(n_windows)
    ]


def _samples_by_process(processes: tuple[str, ...]) -> dict[str, tuple[SignalSample, ...]]:
    requested = set(processes)
    unknown = requested - set(SIGNAL_PROCESSES)
    if unknown:
        raise ValueError(f"Unsupported signal process(es): {', '.join(sorted(unknown))}")

    samples_by_process: dict[str, tuple[SignalSample, ...]] = {}
    for process in processes:
        samples_by_process[process] = tuple(
            sample for sample in SIGNAL_SAMPLES if sample.process == process
        )
        if len(samples_by_process[process]) == 0:
            raise ValueError(f"No signal samples configured for process {process}")

    return samples_by_process


def _build_mumugamma_model(w, process: str) -> None:
    if process == "Z":
        w.factory(
            "RooCBShape::signal_model_boson_cb("
            "boson_mass[70,120], "
            "mean_boson[91.1876, 70, 120], "
            "sigma_boson[2, 0.6, 2], "
            "alpha_boson[3, 0, 3],"
            "n_boson[0.5, 0.1, 3]"
            ")"
        )
        w.factory(
            "Gaussian::signal_model_boson_gauss("
            "boson_mass,"
            "mean_boson,"
            "sigma_boson_gauss[2, 0.6, 2]"
            ")"
        )
    elif process == "H":
        w.factory(
            "RooCBShape::signal_model_boson_cb("
            "boson_mass[100,150], "
            "mean_boson[125, 100, 150], "
            "sigma_boson[1.5, 1, 5], "
            "alpha_boson[1.5, 1, 3],"
            "n_boson[5, 4, 6]"
            ")"
        )
        w.factory(
            "Gaussian::signal_model_boson_gauss("
            "boson_mass,"
            "mean_boson,"
            "sigma_boson_gauss[0.5, 0, 3]"
            ")"
        )
    else:
        raise ValueError(f"Unsupported process {process!r}")

    w.factory(
        "RSUM::signal_model_mumugamma("
        "signal_cb_frac[0.5,0,1]*signal_model_boson_cb,"
        "signal_model_boson_gauss)"
    )
    w.factory(f"upsilon_mass[{UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}]")
    w.factory("weight[-100,100]")


def _safe_float(value) -> float | None:
    value = float(value)
    if not math.isfinite(value):
        return None
    return value


def _parameter_dict(fit_result) -> dict[str, dict[str, float | None]]:
    params = {}
    arg_list = fit_result.floatParsFinal()
    for idx in range(arg_list.getSize()):
        param = arg_list[idx]
        params[param.GetName()] = {
            "value": _safe_float(param.getVal()),
            "error": _safe_float(param.getError()),
        }
    return params


def _compute_chi2_summary(
    model,
    data,
    observable,
    nbins: int,
    nfloatpars: int,
    plot_range: tuple[float, float],
):
    from ROOT import RooAbsData, RooFit  # type: ignore

    plot_nbins = max(nbins + 10, int(round(1.5 * nbins)))
    frame = observable.frame(
        RooFit.Bins(plot_nbins),
        RooFit.Range(float(plot_range[0]), float(plot_range[1])),
    )
    data.plotOn(
        frame,
        RooFit.Name("Data"),
        RooFit.Binning(plot_nbins),
        RooFit.DataError(RooAbsData.SumW2),
    )
    model.plotOn(frame, RooFit.Name("Model"))

    data_hist = frame.findObject("Data")
    npoints = int(data_hist.GetN()) if data_hist is not None else 0
    if npoints <= 0:
        return {
            "chi2": None,
            "ndf": 0,
            "chi2_ndf": None,
            "npoints": npoints,
            "nbins": plot_nbins,
        }

    effective_nfloatpars = min(max(0, int(nfloatpars)), max(0, npoints - 1))
    ndf = max(0, npoints - effective_nfloatpars)
    if ndf <= 0:
        return {
            "chi2": None,
            "ndf": ndf,
            "chi2_ndf": None,
            "npoints": npoints,
            "nbins": plot_nbins,
        }

    chi2_ndf = frame.chiSquare("Model", "Data", effective_nfloatpars)
    chi2_ndf = _safe_float(chi2_ndf)
    chi2 = None if chi2_ndf is None else chi2_ndf * ndf
    return {
        "chi2": _safe_float(chi2) if chi2 is not None else None,
        "ndf": ndf,
        "chi2_ndf": chi2_ndf,
        "npoints": npoints,
        "nbins": plot_nbins,
    }


def _load_window_dataset(job: MumugammaWindowFitJob, w):
    from ROOT import RooArgSet, RooDataSet, RooFit, TFile  # type: ignore

    observables = RooArgSet(w.var("boson_mass"), w.var("upsilon_mass"), w.var("weight"))
    data = RooDataSet(
        "signal_data",
        "signal_data",
        observables,
        RooFit.WeightVar(w.var("weight")),
    )

    for sample in job.samples:
        root_file = TFile.Open(sample.input_file)
        if root_file is None or root_file.IsZombie():
            raise RuntimeError(f"Could not open signal input {sample.input_file}")

        events = root_file.Get("Events")
        if events is None:
            root_file.Close()
            raise RuntimeError(f"Could not find 'Events' tree in {sample.input_file}")

        sample_data = RooDataSet(
            f"signal_data_{sample.state}",
            f"signal_data_{sample.state}",
            observables,
            RooFit.Import(events),
            RooFit.Cut(job.window.cut),
            RooFit.WeightVar(w.var("weight")),
        )
        data.append(sample_data)
        root_file.Close()

    return data


def _empty_result(job: MumugammaWindowFitJob, data) -> dict:
    return {
        "process": job.process,
        "window_index": job.window.index,
        "m_mumu_low": job.window.low,
        "m_mumu_high": job.window.high,
        "m_mumu_high_inclusive": job.window.include_high,
        "m_mumu_cut": job.window.cut,
        "n_unweighted_events": int(data.numEntries()),
        "sum_weights": _safe_float(data.sumEntries()),
        "fit_status": None,
        "cov_qual": None,
        "edm": None,
        "min_nll": None,
        "chi2": None,
        "ndf": None,
        "chi2_ndf": None,
        "chi2_npoints": None,
        "chi2_nbins": None,
        "parameters": {},
        "plot": None,
        "m_mumugamma_plot_range": list(get_signal_boson_plot_range(job.process)),
        "samples": [sample.input_file for sample in job.samples],
        "status": "skipped_empty",
    }


def _fit_mumugamma_window(job: MumugammaWindowFitJob) -> dict:
    from .root_runtime import configure_root

    configure_root()

    from ROOT import RooFit, RooWorkspace  # type: ignore

    from .fastplot import fastplot

    print("\n\n##############################################")
    print(f"######## SIGNAL {job.process} {job.window.label} ########")
    print("##############################################")

    w = RooWorkspace("ws")
    _build_mumugamma_model(w, job.process)
    data = _load_window_dataset(job, w)
    getattr(w, "import")(data)

    if data.numEntries() == 0:
        return _empty_result(job, data)

    model = w.pdf("signal_model_mumugamma")
    fit_result = model.fitTo(data, RooFit.Save(), RooFit.SumW2Error(True))
    nfloatpars = int(fit_result.floatParsFinal().getSize())

    boson_mass = w.var("boson_mass")
    boson_mass.SetTitle(r"m_{#mu#mu#gamma}")
    boson_mass.setUnit(r"GeV")
    w.var("upsilon_mass").SetTitle(r"m_{#mu#mu}")
    w.var("upsilon_mass").setUnit(r"GeV")
    w.pdf("signal_model_boson_gauss").SetTitle(r"Gaussian Component")
    w.pdf("signal_model_boson_cb").SetTitle(r"CB Component")

    plot_range = get_signal_boson_plot_range(job.process)
    chi2_summary = _compute_chi2_summary(
        model,
        data,
        boson_mass,
        job.nbins,
        nfloatpars,
        plot_range,
    )

    plot_file = f"{job.plot_dir}/signal_mumugamma_fit_{job.key}.pdf"
    fastplot(
        model,
        data,
        boson_mass,
        plot_file,
        components=[
            (w.pdf("signal_model_boson_cb"), "CB Component"),
            (w.pdf("signal_model_boson_gauss"), "Gaussian Component"),
        ],
        nbins=job.nbins,
        legend=[0.2, 0.6, 0.5, 0.92]
        if job.process == "H"
        else [0.6, 0.6, 0.93, 0.92],
        is_data=False,
        plot_range=plot_range,
        residual_y_range=(-2.0, 2.0),
        chi2_nfloatpars=nfloatpars,
    )

    print("\n\n--> Fit parameters")
    fit_result.Print("v")
    print("data.numEntries(): ", data.numEntries())
    print("data.sumEntries(): ", data.sumEntries())

    return {
        "process": job.process,
        "window_index": job.window.index,
        "m_mumu_low": job.window.low,
        "m_mumu_high": job.window.high,
        "m_mumu_high_inclusive": job.window.include_high,
        "m_mumu_cut": job.window.cut,
        "n_unweighted_events": int(data.numEntries()),
        "sum_weights": _safe_float(data.sumEntries()),
        "fit_status": int(fit_result.status()),
        "cov_qual": int(fit_result.covQual()),
        "edm": _safe_float(fit_result.edm()),
        "min_nll": _safe_float(fit_result.minNll()),
        "chi2": chi2_summary["chi2"],
        "ndf": chi2_summary["ndf"],
        "chi2_ndf": chi2_summary["chi2_ndf"],
        "chi2_npoints": chi2_summary["npoints"],
        "chi2_nbins": chi2_summary["nbins"],
        "parameters": _parameter_dict(fit_result),
        "plot": plot_file,
        "m_mumugamma_plot_range": list(plot_range),
        "samples": [sample.input_file for sample in job.samples],
        "status": "fit",
    }


def _write_results_json(output_json: str, payload: dict) -> None:
    output_path = Path(output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as output_file:
        json.dump(payload, output_file, indent=2, sort_keys=True)
        output_file.write("\n")


def _build_m_mumu_histograms(
    processes: tuple[str, ...],
    samples_by_process: dict[str, tuple[SignalSample, ...]],
    lower: float,
    upper: float,
    hist_bins: int,
) -> dict[str, dict]:
    from ROOT import TFile  # type: ignore

    if hist_bins <= 0:
        raise ValueError("m_mumu histogram bin count must be positive")

    width = (upper - lower) / hist_bins
    bin_edges = [lower + idx * width for idx in range(hist_bins + 1)]
    histograms: dict[str, dict] = {}
    for process in processes:
        counts = [0.0 for _ in range(hist_bins)]
        for sample in samples_by_process[process]:
            root_file = TFile.Open(sample.input_file)
            if root_file is None or root_file.IsZombie():
                raise RuntimeError(f"Could not open signal input {sample.input_file}")

            events = root_file.Get("Events")
            if events is None:
                root_file.Close()
                raise RuntimeError(f"Could not find 'Events' tree in {sample.input_file}")

            for event in events:
                upsilon_mass = float(event.upsilon_mass)
                if upsilon_mass < lower or upsilon_mass > upper:
                    continue
                bin_index = int((upsilon_mass - lower) / width)
                if bin_index == hist_bins and upsilon_mass == upper:
                    bin_index = hist_bins - 1
                if 0 <= bin_index < hist_bins:
                    counts[bin_index] += float(event.weight)

            root_file.Close()

        histograms[process] = {
            "bin_edges": bin_edges,
            "counts": counts,
            "hist_bins": hist_bins,
            "source": "weighted_signal_events",
        }

    return histograms


def _summary_value(result: dict, parameter_name: str, field: str) -> float | None:
    parameter = result.get("parameters", {}).get(parameter_name, {})
    value = parameter.get(field)
    if value is None:
        return None
    return _safe_float(value)


def _normalize_excluded_windows(excluded_windows) -> tuple[int, ...]:
    if excluded_windows is None:
        return ()
    normalized = tuple(sorted(set(int(window) for window in excluded_windows)))
    if any(window < 0 for window in normalized):
        raise ValueError("Excluded summary window indices must be non-negative")
    return normalized


def _summary_rows(
    payload: dict,
    processes: tuple[str, ...] | None,
    excluded_windows: tuple[int, ...] | None = None,
) -> dict[str, list[dict]]:
    requested = set(processes) if processes is not None else None
    excluded_window_set = set(_normalize_excluded_windows(excluded_windows))
    rows_by_process: dict[str, list[dict]] = {}
    for result in payload.get("results", []):
        process = result.get("process")
        if process is None or (requested is not None and process not in requested):
            continue
        if result.get("status") != "fit":
            continue
        window_index = int(result.get("window_index", -1))
        if window_index in excluded_window_set:
            continue

        mean = _summary_value(result, "mean_boson", "value")
        sigma = _summary_value(result, "sigma_boson", "value")
        if mean is None or sigma is None:
            continue

        low = _safe_float(result["m_mumu_low"])
        high = _safe_float(result["m_mumu_high"])
        rows_by_process.setdefault(process, []).append(
            {
                "window_index": window_index,
                "low": low,
                "high": high,
                "midpoint": 0.5 * (low + high),
                "n_unweighted_events": int(
                    result.get("n_unweighted_events", result.get("n_events", 0))
                ),
                "mean": mean,
                "mean_error": _summary_value(result, "mean_boson", "error"),
                "sigma": sigma,
                "sigma_error": _summary_value(result, "sigma_boson", "error"),
            }
        )

    return {
        process: sorted(rows, key=lambda row: row["midpoint"])
        for process, rows in rows_by_process.items()
        if len(rows) > 0
    }


def _finite_values(values) -> list[float]:
    return [float(value) for value in values if value is not None and math.isfinite(value)]


def _apply_headroom(axis, values, errors=None, *, zero_floor=False) -> None:
    finite_values = _finite_values(values)
    if not finite_values:
        return

    if errors is None:
        lows = finite_values
        highs = finite_values
    else:
        finite_errors = [0.0 if error is None else float(error) for error in errors]
        lows = [value - error for value, error in zip(finite_values, finite_errors)]
        highs = [value + error for value, error in zip(finite_values, finite_errors)]

    low = min(lows)
    high = max(highs)
    if zero_floor:
        low = 0.0
    span = high - low
    if span <= 0.0:
        span = max(abs(high), 1.0) * 0.1

    bottom = 0.0 if zero_floor else low - 0.08 * span
    axis.set_ylim(bottom, high + 0.5 * span)


def _window_edges(rows: list[dict]) -> list[float]:
    if not rows:
        return []
    edges = [rows[0]["low"]]
    edges.extend(row["high"] for row in rows)
    return edges


def _histogram_line(payload: dict, process: str, rows: list[dict]) -> tuple[list[float], list[float], str]:
    histograms = payload.get("m_mumu_histograms", {})
    histogram = histograms.get(process)
    if histogram is not None:
        bin_edges = histogram.get("bin_edges", [])
        counts = histogram.get("counts", [])
        centers = [0.5 * (bin_edges[idx] + bin_edges[idx + 1]) for idx in range(len(counts))]
        source = histogram.get("source")
        label = (
            r"weighted $m_{\mu\mu}$ histogram"
            if source == "weighted_signal_events"
            else r"$m_{\mu\mu}$ histogram"
        )
        return centers, counts, label

    return (
        [row["midpoint"] for row in rows],
        [row["n_unweighted_events"] for row in rows],
        r"$m_{\mu\mu}$ window counts",
    )


def make_correlation_summary_plot(
    input_json: str = RESULTS_JSON,
    output_plot: str = SUMMARY_PLOT,
    processes: tuple[str, ...] | None = None,
    excluded_windows: tuple[int, ...] | None = None,
) -> str:
    input_path = Path(input_json)
    if not input_path.exists():
        raise FileNotFoundError(f"Could not find correlation-study JSON: {input_json}")

    with input_path.open("r", encoding="utf-8") as input_file:
        payload = json.load(input_file)

    excluded_windows = _normalize_excluded_windows(excluded_windows)
    rows_by_process = _summary_rows(payload, processes, excluded_windows)
    if not rows_by_process:
        raise RuntimeError(f"No fitted correlation-study rows found in {input_json}")
    all_rows_by_process = _summary_rows(payload, processes)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selected_processes = [
        process
        for process in (processes or tuple(payload.get("config", {}).get("processes", ())))
        if process in rows_by_process
    ]
    if not selected_processes:
        selected_processes = sorted(rows_by_process)

    fig, axes = plt.subplots(
        len(selected_processes),
        1,
        figsize=(10.5, 4.2 * len(selected_processes)),
        squeeze=False,
    )
    colors = {
        "histogram": "#1f77b4",
        "mean": "#d62728",
        "sigma": "#2ca02c",
    }

    for axis_row, process in zip(axes, selected_processes):
        ax_hist = axis_row[0]
        ax_mean = ax_hist.twinx()
        ax_sigma = ax_hist.twinx()
        ax_sigma.spines["right"].set_position(("axes", 1.14))
        ax_sigma.spines["right"].set_visible(True)

        rows = rows_by_process[process]
        x_values = [row["midpoint"] for row in rows]
        means = [row["mean"] for row in rows]
        mean_errors = [row["mean_error"] or 0.0 for row in rows]
        sigmas = [row["sigma"] for row in rows]
        sigma_errors = [row["sigma_error"] or 0.0 for row in rows]

        edges = _window_edges(all_rows_by_process.get(process, rows))
        for idx, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
            if idx % 2 == 1:
                ax_hist.axvspan(
                    low,
                    high,
                    color="0.5",
                    alpha=0.10,
                    linewidth=0,
                    zorder=0,
                )
            ax_hist.axvline(
                low,
                color="0.55",
                linestyle=":",
                linewidth=0.6,
                alpha=0.25,
                zorder=1,
            )
        if edges:
            ax_hist.axvline(
                edges[-1],
                color="0.55",
                linestyle=":",
                linewidth=0.6,
                alpha=0.25,
                zorder=1,
            )

        hist_x, hist_counts, hist_label = _histogram_line(payload, process, rows)
        hist_handle = ax_hist.plot(
            hist_x,
            hist_counts,
            color=colors["histogram"],
            drawstyle="steps-mid",
            alpha=0.50,
            linewidth=1.7,
            label=hist_label,
            zorder=2,
        )[0]

        mean_handle = ax_mean.errorbar(
            x_values,
            means,
            yerr=mean_errors,
            color=colors["mean"],
            marker="s",
            capsize=3,
            capthick=1,
            label="CB mean",
        ).lines[0]
        sigma_handle = ax_sigma.errorbar(
            x_values,
            sigmas,
            yerr=sigma_errors,
            color=colors["sigma"],
            marker="^",
            capsize=3,
            capthick=1,
            label="CB sigma",
        ).lines[0]

        ax_hist.set_title(f"{process} signal m_mumugamma fit vs m_mumu window")
        ax_hist.set_xlabel(r"$m_{\mu\mu}$ window midpoint [GeV]")
        ax_hist.set_ylabel(r"weighted $m_{\mu\mu}$ histogram", color=colors["histogram"])
        ax_mean.set_ylabel("CB mean [GeV]", color=colors["mean"])
        ax_sigma.set_ylabel("CB sigma [GeV]", color=colors["sigma"])
        ax_hist.tick_params(axis="y", labelcolor=colors["histogram"])
        ax_mean.tick_params(axis="y", labelcolor=colors["mean"])
        ax_sigma.tick_params(axis="y", labelcolor=colors["sigma"])
        ax_hist.grid(True, axis="both", alpha=0.25)
        _apply_headroom(ax_hist, hist_counts, zero_floor=True)
        _apply_headroom(ax_mean, means, mean_errors)
        _apply_headroom(ax_sigma, sigmas, sigma_errors)
        ax_hist.legend(
            [hist_handle, mean_handle, sigma_handle],
            [hist_label, "CB mean", "CB sigma"],
            loc="upper left",
            ncol=3,
            framealpha=0.92,
        )

    fig.subplots_adjust(right=0.78, hspace=0.35)
    output_path = Path(output_plot)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)
    print(f"Wrote correlation summary plot to {output_plot}")
    return str(output_path)


def run_correlation_study(args: Namespace):
    plot_dir = getattr(args, "plot_dir", PLOT_DIR)
    output_json = getattr(args, "output_json", RESULTS_JSON)
    summary_plot = getattr(args, "summary_plot", SUMMARY_PLOT)
    nbins = getattr(args, "nbins", 60)
    workers = getattr(args, "workers", None)
    n_windows = getattr(args, "windows", 10)
    lower = getattr(args, "m_mumu_min", DEFAULT_M_MUMU_LOWER)
    upper = getattr(args, "m_mumu_max", DEFAULT_M_MUMU_UPPER)
    processes = tuple(getattr(args, "processes", SIGNAL_PROCESSES))
    excluded_summary_windows = _normalize_excluded_windows(
        getattr(args, "exclude_summary_windows", ())
    )
    hist_bins = getattr(args, "m_mumu_hist_bins", None)
    if hist_bins is None:
        hist_bins = max(50, 10 * n_windows)

    os.makedirs(plot_dir, exist_ok=True)
    windows = _build_windows(n_windows, lower, upper)
    samples_by_process = _samples_by_process(processes)

    jobs = []
    for process in processes:
        for window in windows:
            payload = MumugammaWindowFitJob(
                process=process,
                samples=samples_by_process[process],
                window=window,
                nbins=nbins,
                plot_dir=plot_dir,
            )
            jobs.append(
                ParallelJob(
                    key=payload.key,
                    label=payload.label,
                    payload=payload,
                )
            )

    fit_results_by_key = run_parallel_jobs(
        "Signal m_mumugamma window fits",
        jobs,
        _fit_mumugamma_window,
        workers=workers,
    )
    fit_results = [fit_results_by_key[job.key] for job in jobs]

    payload = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "config": {
            "nbins": nbins,
            "windows": n_windows,
            "m_mumu_min": lower,
            "m_mumu_max": upper,
            "processes": list(processes),
            "plot_dir": plot_dir,
            "summary_plot": summary_plot,
            "m_mumu_hist_bins": hist_bins,
            "exclude_summary_windows": list(excluded_summary_windows),
        },
        "samples": {
            process: [sample.input_file for sample in samples]
            for process, samples in samples_by_process.items()
        },
        "m_mumu_histograms": _build_m_mumu_histograms(
            processes,
            samples_by_process,
            lower,
            upper,
            hist_bins,
        ),
        "results": fit_results,
    }
    _write_results_json(output_json, payload)
    print(f"Wrote aggregated results to {output_json}")
    make_correlation_summary_plot(
        output_json,
        summary_plot,
        processes,
        excluded_summary_windows,
    )
    return payload
