from argparse import Namespace
from dataclasses import dataclass
import os

from ROOT import (  # type: ignore
    RooArgSet,  # type: ignore
    RooDataSet,  # type: ignore
    RooFit,  # type: ignore
    RooWorkspace,  # type: ignore
    TFile,  # type: ignore
)

from .fastplot import fastplot
from .mass_ranges import (
    BOSON_MASS_LOWER,
    BOSON_MASS_UPPER,
    UPSILON_MASS_LOWER,
    UPSILON_MASS_SEEDS,
    UPSILON_MASS_UPPER,
    get_signal_boson_plot_range,
    get_signal_upsilon_plot_range,
)
from .ws_helper import set_constant


PLOT_DIR = "plots/signal_fit"


@dataclass(frozen=True)
class SignalSample:
    process: str
    state: str
    input_file: str

    @property
    def inner_file_name(self) -> str:
        return self.input_file.replace("inputs/mass_", "").replace(".root", "")


SIGNAL_SAMPLES = [
    SignalSample("H", "1S", "inputs/mass_H_HToUps1SG_Run2.root"),
    SignalSample("H", "2S", "inputs/mass_H_HToUps2SG_Run2.root"),
    SignalSample("H", "3S", "inputs/mass_H_HToUps3SG_Run2.root"),
    SignalSample("Z", "1S", "inputs/mass_Z_ZToUpsilon1SGamma_Run2.root"),
    SignalSample("Z", "2S", "inputs/mass_Z_ZToUpsilon2SGamma_Run2.root"),
    SignalSample("Z", "3S", "inputs/mass_Z_ZToUpsilon3SGamma_Run2.root"),
]


def _build_signal_model(w: RooWorkspace, sample: SignalSample) -> None:
    if sample.process == "Z":
        w.factory(
            "RooCBShape::signal_model_boson_cb("
            f"boson_mass[70,120], "
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
        w.factory(
            f"RooDoubleCB::signal_model_upsilon("
            f"upsilon_mass[{UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}],"
            f"mean_upsilon[{UPSILON_MASS_SEEDS[sample.state]}, {UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}],"
            "sigma_upsilon[0.1, 0., 1],"
            "alpha1_upsilon[3, 0, 10],"
            "n1_upsilon[1, 0, 30],"
            "alpha2_upsilon[3, 0, 10],"
            "n2_upsilon[1, 0, 30]"
            ")"
        )
    elif sample.process == "H":
        w.factory(
            "RooCBShape::signal_model_boson_cb("
            f"boson_mass[100,150], "
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
        w.factory(
            f"RooDoubleCB::signal_model_upsilon("
            f"upsilon_mass[{UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}],"
            f"mean_upsilon[{UPSILON_MASS_SEEDS[sample.state]}, {UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}],"
            "sigma_upsilon[0.1, 0., 0.3],"
            "alpha1_upsilon[1, 0, 1.5],"
            "n1_upsilon[5, 3.5, 6],"
            "alpha2_upsilon[1, 0, 10],"
            "n2_upsilon[5, 3.5, 6]"
            ")"
        )
    else:
        raise ValueError(f"Unsupported process {sample.process!r}")

    w.factory(
        "RSUM::signal_model_boson("
        "signal_cb_fracs[0.5,0,1]*signal_model_boson_cb,"
        "signal_model_boson_gauss)"
    )
    w.factory("PROD::signal_model(signal_model_boson,signal_model_upsilon)")
    w.factory("weight[-100,100]")


def _fit_signal_sample(sample: SignalSample, nbins: int) -> dict[str, str]:
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

    fit_result = w.pdf("signal_model").fitTo(data, RooFit.Save())

    w.var("boson_mass").SetTitle(r"m_{#mu#mu#gamma}")
    w.var("boson_mass").setUnit(r"GeV")
    w.var("upsilon_mass").SetTitle(r"m_{#mu#mu}")
    w.var("upsilon_mass").setUnit(r"GeV")
    w.pdf("signal_model_boson_gauss").SetTitle(r"Gaussian Component")
    w.pdf("signal_model_boson_cb").SetTitle(r"CB Component")

    boson_plot = f"{PLOT_DIR}/signal_fit_boson_{sample.inner_file_name}.pdf"
    fastplot(
        w.pdf("signal_model"),
        w.data("signal_data"),
        w.var("boson_mass"),
        boson_plot,
        components=[
            (w.pdf("signal_model_boson_cb"), "CB Component"),
            (w.pdf("signal_model_boson_gauss"), "Gaussian Component"),
        ],
        nbins=nbins,
        legend=[0.2, 0.6, 0.5, 0.92]
        if sample.process == "H"
        else [0.6, 0.6, 0.93, 0.92],
        is_data=False,
        plot_range=get_signal_boson_plot_range(sample.process),
    )

    upsilon_plot = f"{PLOT_DIR}/signal_fit_upsilon_{sample.inner_file_name}.pdf"
    fastplot(
        w.pdf("signal_model"),
        w.data("signal_data"),
        w.var("upsilon_mass"),
        upsilon_plot,
        nbins=nbins,
        legend=[0.65, 0.7, 0.9, 0.92],
        is_data=False,
        plot_range=get_signal_upsilon_plot_range(sample.state),
    )

    w = set_constant(w)

    print("\n\n--> Fit parameters")
    fit_result.Print("v")
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
    return [_fit_signal_sample(sample, nbins=nbins) for sample in SIGNAL_SAMPLES]
