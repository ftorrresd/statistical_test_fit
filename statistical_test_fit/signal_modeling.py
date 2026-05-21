from argparse import Namespace
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import math
import multiprocessing as mp
import os

from .mass_ranges import (
    UPSILON_MASS_LOWER,
    UPSILON_MASS_SEEDS,
    UPSILON_MASS_UPPER,
    get_signal_boson_plot_range,
    get_signal_upsilon_plot_range,
)
from .parallel_utils import ParallelJob, resolve_worker_count, run_parallel_jobs


PLOT_DIR = "plots/signal_fit"


@dataclass(frozen=True)
class SignalSample:
    process: str
    state: str
    input_file: str

    @property
    def inner_file_name(self) -> str:
        return os.path.splitext(os.path.basename(self.input_file))[0].removeprefix("mass_")


@dataclass(frozen=True)
class SignalFitJob:
    sample: SignalSample
    nbins: int
    plot_dir: str
    high_fit_effort: bool = False
    high_fit_starts: int = 24
    start_workers: int = 1

    @property
    def key(self) -> str:
        return self.sample.inner_file_name

    @property
    def label(self) -> str:
        return f"{self.sample.process} {self.sample.state}"


@dataclass(frozen=True)
class SignalStartFitJob:
    sample: SignalSample
    start_index: int
    start_values: dict[str, float]


SIGNAL_SAMPLES = [
    SignalSample(
        "H",
        "1S",
        "inputs/selected_ggH_HToUps1SG_M125_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Run2.root",
    ),
    SignalSample(
        "H",
        "2S",
        "inputs/selected_ggH_HToUps2SG_M125_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Run2.root",
    ),
    SignalSample(
        "H",
        "3S",
        "inputs/selected_ggH_HToUps3SG_M125_NNPDF31_TuneCP5_13TeV-powheg-pythia8_Run2.root",
    ),
    SignalSample(
        "Z",
        "1S",
        "inputs/selected_ZToUpsilon1SGamma_TuneCP5_13TeV-amcatnloFXFX-pythia8_Run2.root",
    ),
    SignalSample(
        "Z",
        "2S",
        "inputs/selected_ZToUpsilon2SGamma_TuneCP5_13TeV-amcatnloFXFX-pythia8_Run2.root",
    ),
    SignalSample(
        "Z",
        "3S",
        "inputs/selected_ZToUpsilon3SGamma_TuneCP5_13TeV-amcatnloFXFX-pythia8_Run2.root",
    ),
]


HIGH_EFFORT_STARTS = (
    {},
    {
        "sigma_boson": 0.7,
        "sigma_boson_gauss": 0.4,
        "sigma_upsilon": 0.05,
        "sigma_upsilon_gauss": 0.04,
        "signal_dcb_frac": 0.9,
        "signal_upsilon_dcb_frac": 0.9,
    },
    {
        "sigma_boson": 1.2,
        "sigma_boson_gauss": 0.8,
        "sigma_upsilon": 0.08,
        "sigma_upsilon_gauss": 0.08,
        "signal_dcb_frac": 0.8,
        "signal_upsilon_dcb_frac": 0.8,
    },
    {
        "sigma_boson": 2.5,
        "sigma_boson_gauss": 1.5,
        "sigma_upsilon": 0.2,
        "sigma_upsilon_gauss": 0.15,
        "signal_dcb_frac": 0.7,
        "signal_upsilon_dcb_frac": 0.7,
    },
    {
        "sigma_boson": 4.0,
        "sigma_boson_gauss": 2.5,
        "sigma_upsilon": 0.4,
        "sigma_upsilon_gauss": 0.3,
        "signal_dcb_frac": 0.5,
        "signal_upsilon_dcb_frac": 0.5,
    },
    {
        "alpha1_boson": 0.4,
        "n1_boson": 0.6,
        "alpha2_boson": 0.4,
        "n2_boson": 0.6,
        "alpha1_upsilon": 0.4,
        "n1_upsilon": 0.6,
        "alpha2_upsilon": 0.4,
        "n2_upsilon": 0.6,
        "signal_dcb_frac": 0.9,
        "signal_upsilon_dcb_frac": 0.9,
    },
    {
        "alpha1_boson": 1.0,
        "n1_boson": 1.0,
        "alpha2_boson": 1.0,
        "n2_boson": 1.0,
        "alpha1_upsilon": 1.0,
        "n1_upsilon": 1.0,
        "alpha2_upsilon": 1.0,
        "n2_upsilon": 1.0,
        "signal_dcb_frac": 0.8,
        "signal_upsilon_dcb_frac": 0.8,
    },
    {
        "alpha1_boson": 5.0,
        "n1_boson": 20.0,
        "alpha2_boson": 5.0,
        "n2_boson": 20.0,
        "alpha1_upsilon": 5.0,
        "n1_upsilon": 20.0,
        "alpha2_upsilon": 5.0,
        "n2_upsilon": 20.0,
        "signal_dcb_frac": 0.5,
        "signal_upsilon_dcb_frac": 0.5,
    },
    {
        "alpha1_boson": 0.3,
        "n1_boson": 0.6,
        "alpha2_boson": 3.0,
        "n2_boson": 10.0,
        "alpha1_upsilon": 0.3,
        "n1_upsilon": 0.6,
        "alpha2_upsilon": 3.0,
        "n2_upsilon": 10.0,
    },
    {
        "alpha1_boson": 3.0,
        "n1_boson": 10.0,
        "alpha2_boson": 0.3,
        "n2_boson": 0.6,
        "alpha1_upsilon": 3.0,
        "n1_upsilon": 10.0,
        "alpha2_upsilon": 0.3,
        "n2_upsilon": 0.6,
    },
    {
        "signal_dcb_frac": 0.5,
        "signal_upsilon_dcb_frac": 0.5,
        "sigma_boson_gauss": 1.5,
        "sigma_upsilon_gauss": 0.2,
    },
    {
        "signal_dcb_frac": 0.95,
        "signal_upsilon_dcb_frac": 0.95,
        "sigma_boson": 3.0,
        "sigma_upsilon": 0.2,
    },
    {
        "sigma_boson": 0.3,
        "sigma_boson_gauss": 0.2,
        "sigma_upsilon": 0.03,
        "sigma_upsilon_gauss": 0.03,
        "signal_dcb_frac": 0.7,
        "signal_upsilon_dcb_frac": 0.7,
    },
    {
        "sigma_boson": 0.5,
        "sigma_boson_gauss": 1.0,
        "sigma_upsilon": 0.1,
        "sigma_upsilon_gauss": 0.05,
        "signal_dcb_frac": 0.6,
        "signal_upsilon_dcb_frac": 0.9,
    },
    {
        "sigma_boson": 1.0,
        "sigma_boson_gauss": 0.3,
        "sigma_upsilon": 0.15,
        "sigma_upsilon_gauss": 0.25,
        "signal_dcb_frac": 0.9,
        "signal_upsilon_dcb_frac": 0.6,
    },
    {
        "sigma_boson": 1.8,
        "sigma_boson_gauss": 3.0,
        "sigma_upsilon": 0.06,
        "sigma_upsilon_gauss": 0.4,
        "signal_dcb_frac": 0.75,
        "signal_upsilon_dcb_frac": 0.75,
    },
    {
        "sigma_boson": 3.5,
        "sigma_boson_gauss": 0.7,
        "sigma_upsilon": 0.3,
        "sigma_upsilon_gauss": 0.1,
        "signal_dcb_frac": 0.55,
        "signal_upsilon_dcb_frac": 0.85,
    },
    {
        "sigma_boson": 5.0,
        "sigma_boson_gauss": 4.5,
        "sigma_upsilon": 0.6,
        "sigma_upsilon_gauss": 0.5,
        "signal_dcb_frac": 0.85,
        "signal_upsilon_dcb_frac": 0.55,
    },
    {
        "alpha1_boson": 0.7,
        "n1_boson": 2.0,
        "alpha2_boson": 4.0,
        "n2_boson": 15.0,
        "signal_dcb_frac": 0.65,
    },
    {
        "alpha1_boson": 4.0,
        "n1_boson": 15.0,
        "alpha2_boson": 0.7,
        "n2_boson": 2.0,
        "signal_dcb_frac": 0.85,
    },
    {
        "alpha1_upsilon": 0.7,
        "n1_upsilon": 2.0,
        "alpha2_upsilon": 4.0,
        "n2_upsilon": 15.0,
        "signal_upsilon_dcb_frac": 0.65,
    },
    {
        "alpha1_upsilon": 4.0,
        "n1_upsilon": 15.0,
        "alpha2_upsilon": 0.7,
        "n2_upsilon": 2.0,
        "signal_upsilon_dcb_frac": 0.85,
    },
    {
        "sigma_boson": 0.8,
        "alpha1_boson": 2.5,
        "alpha2_boson": 2.5,
        "sigma_upsilon": 0.08,
        "alpha1_upsilon": 2.5,
        "alpha2_upsilon": 2.5,
    },
    {
        "sigma_boson": 2.0,
        "alpha1_boson": 0.5,
        "alpha2_boson": 5.0,
        "sigma_upsilon": 0.25,
        "alpha1_upsilon": 5.0,
        "alpha2_upsilon": 0.5,
    },
)


Z_FIT_INITIALS = {
    "1S": {
        "mean_boson": 91.1251,
        "sigma_boson": 1.28845,
        "alpha1_boson": 1.13806,
        "n1_boson": 9.12668,
        "alpha2_boson": 2.30193,
        "n2_boson": 4.28881,
        "sigma_boson_gauss": 0.627102,
        "signal_dcb_frac": 0.8392,
        "mean_upsilon": 9.45798,
        "sigma_upsilon": 0.0898103,
        "alpha1_upsilon": 2.00974,
        "n1_upsilon": 1.10676,
        "alpha2_upsilon": 3.17943,
        "n2_upsilon": 1.16456,
        "sigma_upsilon_gauss": 0.206566,
        "signal_upsilon_dcb_frac": 0.768843,
    },
    "2S": {
        "mean_boson": 91.1665,
        "sigma_boson": 1.06662,
        "alpha1_boson": 1.01452,
        "n1_boson": 8.88621,
        "alpha2_boson": 3.11102,
        "n2_boson": 1.15361,
        "sigma_boson_gauss": 2.23579,
        "signal_dcb_frac": 0.852226,
        "mean_upsilon": 10.0179,
        "sigma_upsilon": 0.09,
        "alpha1_upsilon": 2.37887,
        "n1_upsilon": 0.545057,
        "alpha2_upsilon": 2.86702,
        "n2_upsilon": 1.59121,
        "sigma_upsilon_gauss": 0.186,
        "signal_upsilon_dcb_frac": 0.65,
    },
    "3S": {
        "mean_boson": 91.1776,
        "sigma_boson": 1.02036,
        "alpha1_boson": 1.05157,
        "n1_boson": 5.25265,
        "alpha2_boson": 1.99113,
        "n2_boson": 4.30759,
        "sigma_boson_gauss": 1.43633,
        "signal_dcb_frac": 0.885992,
        "mean_upsilon": 10.3499,
        "sigma_upsilon": 0.0889776,
        "alpha1_upsilon": 2.31292,
        "n1_upsilon": 0.73828,
        "alpha2_upsilon": 2.97709,
        "n2_upsilon": 0.973556,
        "sigma_upsilon_gauss": 0.192096,
        "signal_upsilon_dcb_frac": 0.603921,
    },
}


H_FIT_INITIALS = {
    "1S": {
        "mean_boson": 125.02214561,
        "sigma_boson": 1.19655568795,
        "alpha1_boson": 1.52713636172,
        "n1_boson": 2.10600153652,
        "alpha2_boson": 1.53158750437,
        "n2_boson": 4.96311685292,
        "sigma_boson_gauss": 1.96765449134,
        "signal_dcb_frac": 0.66975578095,
        "mean_upsilon": 9.45981940878,
        "sigma_upsilon": 0.0790037498453,
        "alpha1_upsilon": 1.02754631305,
        "n1_upsilon": 4.19311827063,
        "alpha2_upsilon": 0.992182341674,
        "n2_upsilon": 4.96706598492,
        "sigma_upsilon_gauss": 0.0832827965007,
        "signal_upsilon_dcb_frac": 0.78370990067,
    },
    "2S": {
        "mean_boson": 124.995851389,
        "sigma_boson": 1.18103644201,
        "alpha1_boson": 1.51820514651,
        "n1_boson": 2.10864239688,
        "alpha2_boson": 1.52045540453,
        "n2_boson": 4.9651067451,
        "sigma_boson_gauss": 1.84923754962,
        "signal_dcb_frac": 0.696926079,
        "mean_upsilon": 10.0197510483,
        "sigma_upsilon": 0.0750684153274,
        "alpha1_upsilon": 1.14468956752,
        "n1_upsilon": 1.63792767564,
        "alpha2_upsilon": 0.974685799505,
        "n2_upsilon": 4.86009580752,
        "sigma_upsilon_gauss": 0.0904366087853,
        "signal_upsilon_dcb_frac": 0.762266051197,
    },
    "3S": {
        "mean_boson": 124.983611204,
        "sigma_boson": 1.17973712738,
        "alpha1_boson": 1.52421099342,
        "n1_boson": 2.11289370425,
        "alpha2_boson": 1.52109683895,
        "n2_boson": 4.95952000293,
        "sigma_boson_gauss": 1.90656350359,
        "signal_dcb_frac": 0.691539431083,
        "mean_upsilon": 10.3500961844,
        "sigma_upsilon": 0.076391166293,
        "alpha1_upsilon": 1.0993598888,
        "n1_upsilon": 2.22799732628,
        "alpha2_upsilon": 0.983087470533,
        "n2_upsilon": 4.8815918708,
        "sigma_upsilon_gauss": 0.0875083601584,
        "signal_upsilon_dcb_frac": 0.769236387607,
    },
}


def _build_signal_model(w, sample: SignalSample) -> None:
    boson_dcb_fraction = "signal_dcb_frac[0.9,0.5,1]"
    upsilon_dcb_fraction = "signal_upsilon_dcb_frac[0.9,0.5,1]"

    if sample.process == "Z":
        z_initials = Z_FIT_INITIALS[sample.state]
        boson_dcb_fraction = (
            f"signal_dcb_frac[{z_initials['signal_dcb_frac']},0.5,1]"
        )
        upsilon_dcb_fraction = (
            f"signal_upsilon_dcb_frac[{z_initials['signal_upsilon_dcb_frac']},0.6,1]"
        )
        w.factory(
            "RooDoubleCB::signal_model_boson_dcb("
            "boson_mass[70,120], "
            f"mean_boson[{z_initials['mean_boson']}, 70, 120], "
            f"sigma_boson[{z_initials['sigma_boson']}, 1E-3, 6], "
            f"alpha1_boson[{z_initials['alpha1_boson']}, 0.1, 20],"
            f"n1_boson[{z_initials['n1_boson']}, 0.5, 30],"
            f"alpha2_boson[{z_initials['alpha2_boson']}, 0.1, 20],"
            f"n2_boson[{z_initials['n2_boson']}, 0.5, 30]"
            ")"
        )
        w.factory(
            "Gaussian::signal_model_boson_gauss("
            "boson_mass,"
            "mean_boson,"
            f"sigma_boson_gauss[{z_initials['sigma_boson_gauss']}, 1E-3, 6]"
            ")"
        )
        w.factory(
            f"RooDoubleCB::signal_model_upsilon_dcb("
            f"upsilon_mass[{UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}],"
            f"mean_upsilon[{z_initials['mean_upsilon']}, {UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}],"
            f"sigma_upsilon[{z_initials['sigma_upsilon']}, 1E-3, 6],"
            f"alpha1_upsilon[{z_initials['alpha1_upsilon']}, 0.1, 20],"
            f"n1_upsilon[{z_initials['n1_upsilon']}, 0.5, 30],"
            f"alpha2_upsilon[{z_initials['alpha2_upsilon']}, 0.1, 20],"
            f"n2_upsilon[{z_initials['n2_upsilon']}, 0.5, 30]"
            ")"
        )
        w.factory(
            "Gaussian::signal_model_upsilon_gauss("
            "upsilon_mass,"
            "mean_upsilon,"
            f"sigma_upsilon_gauss[{z_initials['sigma_upsilon_gauss']}, 0.15, 6]"
            ")"
        )
    elif sample.process == "H":
        h_initials = H_FIT_INITIALS[sample.state]
        boson_dcb_fraction = (
            f"signal_dcb_frac[{h_initials['signal_dcb_frac']},0.5,1]"
        )
        upsilon_dcb_fraction = (
            f"signal_upsilon_dcb_frac[{h_initials['signal_upsilon_dcb_frac']},0.6,1]"
        )
        w.factory(
            "RooDoubleCB::signal_model_boson_dcb("
            "boson_mass[100,150], "
            f"mean_boson[{h_initials['mean_boson']}, 100, 150], "
            f"sigma_boson[{h_initials['sigma_boson']}, 1E-3, 6], "
            f"alpha1_boson[{h_initials['alpha1_boson']}, 0.1, 20],"
            f"n1_boson[{h_initials['n1_boson']}, 0.5, 30],"
            f"alpha2_boson[{h_initials['alpha2_boson']}, 0.1, 20],"
            f"n2_boson[{h_initials['n2_boson']}, 0.5, 30]"
            ")"
        )
        w.factory(
            "Gaussian::signal_model_boson_gauss("
            "boson_mass,"
            "mean_boson,"
            f"sigma_boson_gauss[{h_initials['sigma_boson_gauss']}, 1E-3, 6]"
            ")"
        )
        w.factory(
            f"RooDoubleCB::signal_model_upsilon_dcb("
            f"upsilon_mass[{UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}],"
            f"mean_upsilon[{h_initials['mean_upsilon']}, {UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}],"
            f"sigma_upsilon[{h_initials['sigma_upsilon']}, 1E-3, 6],"
            f"alpha1_upsilon[{h_initials['alpha1_upsilon']}, 0.1, 20],"
            f"n1_upsilon[{h_initials['n1_upsilon']}, 0.5, 30],"
            f"alpha2_upsilon[{h_initials['alpha2_upsilon']}, 0.1, 20],"
            f"n2_upsilon[{h_initials['n2_upsilon']}, 0.5, 30]"
            ")"
        )
        w.factory(
            "Gaussian::signal_model_upsilon_gauss("
            "upsilon_mass,"
            "mean_upsilon,"
            f"sigma_upsilon_gauss[{h_initials['sigma_upsilon_gauss']}, 1E-3, 6]"
            ")"
        )
    else:
        raise ValueError(f"Unsupported process {sample.process!r}")

    w.factory(
        "RSUM::signal_model_boson("
        f"{boson_dcb_fraction}*signal_model_boson_dcb,"
        "signal_model_boson_gauss)"
    )
    w.factory(
        "RSUM::signal_model_upsilon("
        f"{upsilon_dcb_fraction}*signal_model_upsilon_dcb,"
        "signal_model_upsilon_gauss)"
    )
    w.factory("PROD::signal_model(signal_model_boson,signal_model_upsilon)")
    w.factory("weight[-100,100]")


def _configure_high_fit_effort(ROOT) -> None:
    ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2", "Migrad")
    ROOT.Math.MinimizerOptions.SetDefaultStrategy(2)
    ROOT.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(100000)
    ROOT.Math.MinimizerOptions.SetDefaultMaxIterations(100000)
    ROOT.Math.MinimizerOptions.SetDefaultTolerance(1e-4)


def _signal_parameter_names(w, data) -> list[str]:
    return [p.GetName() for p in w.pdf("signal_model").getParameters(data.get())]


def _capture_parameter_state(w, parameter_names: list[str]) -> dict[str, tuple[float, float]]:
    state = {}
    for name in parameter_names:
        var = w.var(name)
        if var is not None:
            state[name] = (float(var.getVal()), float(var.getError()))
    return state


def _restore_parameter_state(
    w,
    state: dict[str, tuple[float, float]],
    *,
    make_float: bool = False,
) -> None:
    for name, (value, error) in state.items():
        var = w.var(name)
        if var is None:
            continue
        if make_float:
            var.setConstant(False)
        var.setVal(value)
        var.setError(error)


def _apply_start_values(w, start_values: dict[str, float]) -> None:
    for name, value in start_values.items():
        var = w.var(name)
        if var is None:
            continue
        lower = float(var.getMin())
        upper = float(var.getMax())
        var.setVal(min(max(float(value), lower), upper))


def _fit_with_default_effort(w, data, RooFit):
    return w.pdf("signal_model").fitTo(data, RooFit.Save())


def _run_high_effort_start(payload: SignalStartFitJob) -> dict:
    from .root_runtime import configure_root

    configure_root()

    import ROOT  # type: ignore

    from ROOT import (  # type: ignore
        RooArgSet,  # type: ignore
        RooDataSet,  # type: ignore
        RooFit,  # type: ignore
        RooWorkspace,  # type: ignore
        TFile,  # type: ignore
    )

    _configure_high_fit_effort(ROOT)

    w = RooWorkspace("ws")
    _build_signal_model(w, payload.sample)

    root_file = TFile.Open(payload.sample.input_file)
    if root_file is None or root_file.IsZombie():
        raise RuntimeError(f"Could not open signal input {payload.sample.input_file}")

    try:
        events = root_file.Get("Events")
        if events is None:
            raise RuntimeError(
                f"Could not find 'Events' tree in {payload.sample.input_file}"
            )

        data = RooDataSet(
            "signal_data",
            "signal_data",
            RooArgSet(w.var("boson_mass"), w.var("upsilon_mass"), w.var("weight")),
            RooFit.Import(events),
            RooFit.WeightVar(w.var("weight")),
        )
        getattr(w, "import")(data)

        parameter_names = _signal_parameter_names(w, data)
        _apply_start_values(w, payload.start_values)
        fit_result = w.pdf("signal_model").fitTo(
            data,
            RooFit.Save(True),
            RooFit.SumW2Error(True),
            RooFit.Strategy(2),
            RooFit.Minimizer("Minuit2", "Migrad"),
            RooFit.Hesse(True),
            RooFit.PrintLevel(-1),
            RooFit.Verbose(False),
        )
        if fit_result is None:
            return {
                "start_index": payload.start_index,
                "min_nll": math.inf,
                "status": -1,
                "cov_qual": -1,
                "state": None,
                "message": "fit returned None",
            }

        return {
            "start_index": payload.start_index,
            "min_nll": float(fit_result.minNll()),
            "status": int(fit_result.status()),
            "cov_qual": int(fit_result.covQual()),
            "state": _capture_parameter_state(w, parameter_names),
            "message": "ok",
        }
    finally:
        root_file.Close()


def _run_high_effort_starts(
    sample: SignalSample,
    start_workers: int,
    high_fit_starts: int,
) -> list[dict]:
    if high_fit_starts < 1:
        raise ValueError("--high-fit-starts must be at least 1")
    if high_fit_starts > len(HIGH_EFFORT_STARTS):
        raise ValueError(
            f"--high-fit-starts={high_fit_starts} exceeds the "
            f"{len(HIGH_EFFORT_STARTS)} configured deterministic starts"
        )
    start_jobs = [
        SignalStartFitJob(
            sample=sample,
            start_index=index,
            start_values=dict(start_values),
        )
        for index, start_values in enumerate(
            HIGH_EFFORT_STARTS[:high_fit_starts],
            start=1,
        )
    ]
    max_workers = max(1, min(int(start_workers), len(start_jobs)))
    print(
        f"--> High fit effort: trying {len(start_jobs)} start points "
        f"with {max_workers} start worker(s)"
    )

    if max_workers == 1:
        results = []
        for start_job in start_jobs:
            try:
                results.append(_run_high_effort_start(start_job))
            except Exception as exc:
                results.append(
                    {
                        "start_index": start_job.start_index,
                        "min_nll": math.inf,
                        "status": -1,
                        "cov_qual": -1,
                        "state": None,
                        "message": repr(exc),
                    }
                )
        return results

    results = []
    context = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=max_workers, mp_context=context) as executor:
        future_to_start = {
            executor.submit(_run_high_effort_start, start_job): start_job
            for start_job in start_jobs
        }
        for future in as_completed(future_to_start):
            start_job = future_to_start[future]
            try:
                results.append(future.result())
            except Exception as exc:
                results.append(
                    {
                        "start_index": start_job.start_index,
                        "min_nll": math.inf,
                        "status": -1,
                        "cov_qual": -1,
                        "state": None,
                        "message": repr(exc),
                    }
                )

    return sorted(results, key=lambda result: int(result["start_index"]))


def _final_high_effort_fit(w, data, RooFit):
    return w.pdf("signal_model").fitTo(
        data,
        RooFit.Save(True),
        RooFit.SumW2Error(True),
        RooFit.Strategy(2),
        RooFit.Minimizer("Minuit2", "Migrad"),
        RooFit.Hesse(True),
        RooFit.PrintLevel(-1),
        RooFit.Verbose(False),
    )


def _fit_with_high_effort(w, data, ROOT, RooFit, job: SignalFitJob):
    _configure_high_fit_effort(ROOT)

    parameter_names = _signal_parameter_names(w, data)
    nominal_state = _capture_parameter_state(w, parameter_names)
    best_result = None
    best_state = None
    best_nll = math.inf

    for result in _run_high_effort_starts(
        job.sample,
        job.start_workers,
        job.high_fit_starts,
    ):
        start_index = int(result["start_index"])
        if result["state"] is None:
            print(f"    start {start_index}: failed ({result['message']})")
            continue

        min_nll = float(result["min_nll"])
        status = int(result["status"])
        cov_qual = int(result["cov_qual"])
        print(
            f"    start {start_index}: minNll={min_nll:.6g}, "
            f"status={status}, covQual={cov_qual}"
        )
        if math.isfinite(min_nll) and min_nll < best_nll:
            best_nll = min_nll
            best_result = result
            best_state = result["state"]

    if best_state is None or best_result is None:
        print("--> High fit effort did not find a valid multistart result; keeping nominal fit")
        _restore_parameter_state(w, nominal_state, make_float=True)
        return _final_high_effort_fit(w, data, RooFit)

    _restore_parameter_state(w, best_state, make_float=True)
    print(f"--> High fit effort selected minNll={best_nll:.6g}")
    fit_result = _final_high_effort_fit(w, data, RooFit)
    if fit_result is None:
        _restore_parameter_state(w, best_state, make_float=True)
    return fit_result


def _fit_signal_sample(job: SignalFitJob) -> dict[str, str]:
    from .root_runtime import configure_root

    configure_root()

    import ROOT  # type: ignore

    from ROOT import (  # type: ignore
        RooArgSet,  # type: ignore
        RooDataSet,  # type: ignore
        RooFit,  # type: ignore
        RooWorkspace,  # type: ignore
        TFile,  # type: ignore
    )

    from .fastplot import fastplot
    from .ws_helper import set_constant

    sample = job.sample

    print("\n\n##############################################")
    print(
        f"################## SIGNAL {sample.process} {sample.state} ##################"
    )
    print("##############################################")

    w = RooWorkspace("ws")
    _build_signal_model(w, sample)

    root_file = TFile.Open(sample.input_file)
    if root_file is None or root_file.IsZombie():
        raise RuntimeError(f"Could not open signal input {sample.input_file}")

    events = root_file.Get("Events")
    if events is None:
        raise RuntimeError(f"Could not find 'Events' tree in {sample.input_file}")

    data = RooDataSet(
        "signal_data",
        "signal_data",
        RooArgSet(w.var("boson_mass"), w.var("upsilon_mass"), w.var("weight")),
        RooFit.Import(events),
        RooFit.WeightVar(w.var("weight")),
    )
    getattr(w, "import")(data)

    if job.high_fit_effort:
        fit_result = _fit_with_high_effort(w, data, ROOT, RooFit, job)
    else:
        fit_result = _fit_with_default_effort(w, data, RooFit)

    w.var("boson_mass").SetTitle(r"m_{#mu#mu#gamma}")
    w.var("boson_mass").setUnit(r"GeV")
    w.var("upsilon_mass").SetTitle(r"m_{#mu#mu}")
    w.var("upsilon_mass").setUnit(r"GeV")
    w.pdf("signal_model_boson_dcb").SetTitle(r"Double-sided CB Component")
    w.pdf("signal_model_boson_gauss").SetTitle(r"Gaussian Component")
    w.pdf("signal_model_upsilon_dcb").SetTitle(r"Double-sided CB Component")
    w.pdf("signal_model_upsilon_gauss").SetTitle(r"Gaussian Component")

    boson_plot = f"{job.plot_dir}/signal_fit_boson_{sample.inner_file_name}.pdf"
    fastplot(
        w.pdf("signal_model"),
        w.data("signal_data"),
        w.var("boson_mass"),
        boson_plot,
        components=[
            (w.pdf("signal_model_boson_dcb"), "Double-sided CB Component"),
            (w.pdf("signal_model_boson_gauss"), "Gaussian Component"),
        ],
        nbins=job.nbins,
        legend=[0.2, 0.6, 0.5, 0.92]
        if sample.process == "H"
        else [0.6, 0.6, 0.93, 0.92],
        is_data=False,
        model_legend_name="Double-sided CB + Gaussian",
        plot_range=get_signal_boson_plot_range(sample.process),
        residual_y_range=(-2.0, 2.0),
    )

    upsilon_plot = f"{job.plot_dir}/signal_fit_upsilon_{sample.inner_file_name}.pdf"
    fastplot(
        w.pdf("signal_model"),
        w.data("signal_data"),
        w.var("upsilon_mass"),
        upsilon_plot,
        components=[
            (w.pdf("signal_model_upsilon_dcb"), "Double-sided CB Component"),
            (w.pdf("signal_model_upsilon_gauss"), "Gaussian Component"),
        ],
        nbins=job.nbins,
        legend=[0.65, 0.7, 0.9, 0.92],
        is_data=False,
        model_legend_name="Double-sided CB + Gaussian",
        plot_range=get_signal_upsilon_plot_range(sample.state),
        residual_y_range=(-2.0, 2.0),
    )

    w = set_constant(w)

    print("\n\n--> Fit parameters")
    if fit_result is not None:
        fit_result.Print("v")
    else:
        print("Fit result unavailable; workspace parameters reflect the selected start state.")
    print("data.sumEntries(): ", data.sumEntries())
    print("model Integral: ", w.pdf("signal_model").createIntegral(data.get()).getVal())

    workspace_file = f"signal_workspace_{sample.inner_file_name}.root"
    w.writeToFile(workspace_file)
    root_file.Close()

    return {
        "workspace": workspace_file,
        "boson_plot": boson_plot,
        "upsilon_plot": upsilon_plot,
    }


def run_signal_modeling(args: Namespace):
    os.makedirs(PLOT_DIR, exist_ok=True)
    nbins = getattr(args, "nbins", 60)
    workers = getattr(args, "workers", None)
    high_fit_effort = bool(getattr(args, "high_fit_effort", False))
    high_fit_starts = int(getattr(args, "high_fit_starts", 24))
    selected_process = getattr(args, "process", "all")
    if selected_process == "all":
        selected_samples = SIGNAL_SAMPLES
    else:
        selected_samples = [
            sample for sample in SIGNAL_SAMPLES if sample.process == selected_process
        ]

    if not selected_samples:
        raise ValueError(f"No signal samples selected for process {selected_process!r}")

    outer_workers = resolve_worker_count(workers, len(selected_samples))
    start_workers = 1
    if high_fit_effort:
        if high_fit_starts < 1:
            raise ValueError("--high-fit-starts must be at least 1")
        if high_fit_starts > len(HIGH_EFFORT_STARTS):
            raise ValueError(
                f"--high-fit-starts={high_fit_starts} exceeds the "
                f"{len(HIGH_EFFORT_STARTS)} configured deterministic starts"
            )
        cpu_count = os.cpu_count() or 1
        start_workers = max(
            1,
            min(high_fit_starts, math.ceil(cpu_count / max(1, outer_workers))),
        )
        print(
            "High fit effort enabled: "
            f"{high_fit_starts} start point(s), "
            f"up to {outer_workers} sample worker(s), "
            f"{start_workers} start worker(s) per sample."
        )
    jobs = [
        ParallelJob(
            key=sample.inner_file_name,
            label=f"{sample.process} {sample.state}",
            payload=SignalFitJob(
                sample=sample,
                nbins=nbins,
                plot_dir=PLOT_DIR,
                high_fit_effort=high_fit_effort,
                high_fit_starts=high_fit_starts,
                start_workers=start_workers,
            ),
        )
        for sample in selected_samples
    ]
    results = run_parallel_jobs(
        "Signal fits",
        jobs,
        _fit_signal_sample,
        workers=workers,
    )
    return [results[sample.inner_file_name] for sample in selected_samples]
