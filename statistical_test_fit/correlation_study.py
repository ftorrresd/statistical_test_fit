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
SIGNAL_PROCESSES = ("H", "Z")


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
        "n_events": int(data.numEntries()),
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
        "workspace": None,
        "m_mumugamma_plot_range": list(get_signal_boson_plot_range(job.process)),
        "samples": [sample.input_file for sample in job.samples],
        "status": "skipped_empty",
    }


def _fit_mumugamma_window(job: MumugammaWindowFitJob) -> dict:
    from .root_runtime import configure_root

    configure_root()

    from ROOT import RooFit, RooWorkspace  # type: ignore

    from .fastplot import fastplot

    from .ws_helper import set_constant

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

    w = set_constant(w)
    workspace_file = f"{job.plot_dir}/signal_mumugamma_workspace_{job.key}.root"
    w.writeToFile(workspace_file)

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
        "n_events": int(data.numEntries()),
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
        "workspace": workspace_file,
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


def run_correlation_study(args: Namespace):
    plot_dir = getattr(args, "plot_dir", PLOT_DIR)
    output_json = getattr(args, "output_json", RESULTS_JSON)
    nbins = getattr(args, "nbins", 60)
    workers = getattr(args, "workers", None)
    n_windows = getattr(args, "windows", 10)
    lower = getattr(args, "m_mumu_min", UPSILON_MASS_LOWER)
    upper = getattr(args, "m_mumu_max", UPSILON_MASS_UPPER)
    processes = tuple(getattr(args, "processes", SIGNAL_PROCESSES))

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
        },
        "samples": {
            process: [sample.input_file for sample in samples]
            for process, samples in samples_by_process.items()
        },
        "results": fit_results,
    }
    _write_results_json(output_json, payload)
    print(f"Wrote aggregated results to {output_json}")
    return payload
