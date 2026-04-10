from argparse import Namespace
from dataclasses import dataclass, field
import json
import os
from enum import Enum
from pprint import pprint

from ROOT import (  # type: ignore
    RooArgSet,  # type: ignore
    RooArgList,  # type: ignore
    RooDataSet,  # type: ignore
    RooDoubleCB,  # type: ignore
    RooBernstein,  # type: ignore
    RooProdPdf,  # type: ignore
    RooFit,  # type: ignore
    RooWorkspace,  # type: ignore
    TFile,  # type: ignore
    TMath,  # type: ignore
)

from .fastplot import fastplot
from .mass_ranges import (
    BOSON_MASS_LOWER,
    BOSON_MASS_UPPER,
    RESONANT_CR_UPSILON_MASS_LOWER,
    RESONANT_CR_UPSILON_MASS_UPPER,
    UPSILON_MASS_LOWER,
    UPSILON_MASS_UPPER,
    get_signal_boson_plot_range,
)
from .normalization_fit import fit_and_plot
from .ws_helper import *


PLOT_DIR = "plots/resonant_background"


def build_resonant_background_Higgs_ws(plot_dir=PLOT_DIR, nbins=80):
    """
    Resonant Background Modeling
    """

    input_file = "inputs/mass_H_GluGluHToMuMuG_M125_MLL-0To60_Dalitz_012j_Run2.root"

    w = RooWorkspace("resonant_background_Higgs_ws")

    ### boson model

    w.factory(
        "RooDoubleCB::resonant_background_model_Higgs_boson("
        f"boson_mass[{BOSON_MASS_LOWER}, {BOSON_MASS_UPPER}], "
        "mean_boson[125, 57, 200], "
        "sigma_boson[2, 0.5, 4],"
        "alpha1_upsilon[3, 0, 10],"
        "n1_upsilon[0.5, 0.1, 57],"
        "alpha2_upsilon[3, 0, 10],"
        "n2_upsilon[0.5, 0.1, 57]"
        ")"
    )

    # upsilon model
    w.factory(
        f"RooBernstein::resonant_background_model_Higgs_upsilon(upsilon_mass[{UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}],{{1,p1[5, 0, 10]}})"
    )

    # 2D model
    w.factory(
        "PROD::resonant_background_model_Higgs(resonant_background_model_Higgs_boson,resonant_background_model_Higgs_upsilon)"
    )

    w.factory("weight[-100,100]")

    # load data
    f = TFile.Open(input_file)
    data = RooDataSet(
        "resonant_background_data",
        "resonant_background_data",
        RooArgSet(w.var("boson_mass"), w.var("upsilon_mass"), w.var("weight")),
        RooFit.Import(f.Events),
        RooFit.WeightVar(w.var("weight")),
    )
    # getattr(w, "import")(data)

    # fit to data
    fit_result = w.pdf("resonant_background_model_Higgs").fitTo(
        data, RooFit.Save(), RooFit.SumW2Error(True)
    )

    # plot data the and the pdf
    print("\n\n--> Saving plot ")
    w.var("boson_mass").SetTitle(r"m_{#mu#mu#gamma}")
    w.var("boson_mass").setUnit(r"GeV")
    w.var("upsilon_mass").SetTitle(r"m_{#mu#mu}")
    w.var("upsilon_mass").setUnit(r"GeV")
    # w.pdf("resonant_background_model_boson_gauss").SetTitle(r"Gaussian Component")
    # w.pdf("resonant_background_model_boson_cb").SetTitle(r"CB Component")
    # w.pdf("signal_model_upsilon").SetTitle(r"CB Component")

    # plot boson mass fit
    fastplot(
        w.pdf("resonant_background_model_Higgs"),
        data,
        w.var("boson_mass"),
        f"{plot_dir}/resonant_background_fit_HiggsDalitz_MC_m_mumugamma.pdf",
        nbins=nbins,
        legend=[0.6, 0.6, 0.9, 0.92],  # legend=[0.2, 0.6, 0.5, 0.92], #type: ignore
        is_data=False,
        plot_range=get_signal_boson_plot_range("H"),
    )

    # plot upsilon mass fit
    fastplot(
        w.pdf("resonant_background_model_Higgs"),
        data,
        w.var("upsilon_mass"),
        f"{plot_dir}/resonant_background_fit_HiggsDalitz_MC_m_mumu.pdf",
        # components=[
        #             (upsilon_1S, 10),
        #             (upsilon_2S, 10),
        #             (upsilon_3S, 10),
        #             (background, 20)
        #             ],
        nbins=nbins,
        legend=[0.2, 0.7, 0.53, 0.9],  # legend=[0.6, 0.2, 0.93, 0.4], #type: ignore
        is_data=False,
    )

    # setting all var as constants
    w = set_constant(w)

    print("\n\n--> Fit parameters ")
    fit_result.Print("v")

    print("data.sumEntries(): ", data.sumEntries())

    w.SaveAs("resonant_background_fit_HiggsDalitz.root")


def build_resonant_background_Z_ws(plot_dir=PLOT_DIR, nbins=80):
    """
    Resonant Background Modeling
    """

    input_file = "inputs/mass_Z_ZGTo2MuG_MMuMu-2To15_Run2.root"

    w = RooWorkspace("resonant_background_Z_ws")

    w.factory("weight[-100,100]")
    w.factory(f"upsilon_mass[{UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}]")
    w.factory(f"boson_mass[{BOSON_MASS_LOWER}, {BOSON_MASS_UPPER}]")

    resonant_background_model_Z = build_resonant_background_modeling_Z(
        w.var("boson_mass"),
    )
    getattr(w, "import")(resonant_background_model_Z)

    # upsilon model
    w.factory(
        f"RooBernstein::resonant_background_model_Z_upsilon(upsilon_mass[{UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}],{{1,p1[5, 0, 10]}})"
    )

    # 2D model
    w.factory(
        "PROD::resonant_background_model_Z(resonant_background_model_ZG_boson,resonant_background_model_Z_upsilon)"
    )

    # load data
    f = TFile.Open(input_file)
    data = RooDataSet(
        "resonant_background_data",
        "resonant_background_data",
        RooArgSet(w.var("boson_mass"), w.var("upsilon_mass"), w.var("weight")),
        RooFit.Import(f.Events),
        RooFit.WeightVar(w.var("weight")),
    )
    # getattr(w, "import")(data)

    # fit to data
    fit_result = w.pdf("resonant_background_model_Z").fitTo(
        data, RooFit.Save(), RooFit.SumW2Error(True)
    )

    # plot data the and the pdf
    print("\n\n--> Saving plot ")
    w.var("boson_mass").SetTitle(r"m_{#mu#mu#gamma}")
    w.var("boson_mass").setUnit(r"GeV")
    w.var("upsilon_mass").SetTitle(r"m_{#mu#mu}")
    w.var("upsilon_mass").setUnit(r"GeV")
    # w.pdf("resonant_background_model_boson_gauss").SetTitle(r"Gaussian Component")
    # w.pdf("resonant_background_model_boson_cb").SetTitle(r"CB Component")
    # w.pdf("signal_model_upsilon").SetTitle(r"CB Component")

    # plot boson mass fit
    fastplot(
        w.pdf("resonant_background_model_Z"),
        data,
        w.var("boson_mass"),
        f"{plot_dir}/resonant_background_fit_Z_MC_m_mumugamma.pdf",
        nbins=nbins,
        legend=[0.6, 0.6, 0.9, 0.92],  # legend=[0.2, 0.6, 0.5, 0.92], #type: ignore
        is_data=False,
        plot_range=get_signal_boson_plot_range("Z"),
    )

    # plot upsilon mass fit
    fastplot(
        w.pdf("resonant_background_model_Z"),
        data,
        w.var("upsilon_mass"),
        f"{plot_dir}/resonant_background_fit_Z_MC_m_mumu.pdf",
        # components=[
        #             (upsilon_1S, 10),
        #             (upsilon_2S, 10),
        #             (upsilon_3S, 10),
        #             (background, 20)
        #             ],
        nbins=nbins,
        legend=[0.2, 0.7, 0.53, 0.9],  # legend=[0.6, 0.2, 0.93, 0.4], #type: ignore
        is_data=False,
    )

    # setting all var as constants
    w = set_constant(w)

    print("\n\n--> Fit parameters ")
    fit_result.Print("v")

    print("data.sumEntries(): ", data.sumEntries())

    w.SaveAs("resonant_background_fit_ZGamma.root")


def build_resonant_background_modeling_Z(boson_mass, sufix=None):
    if sufix is not None:
        sufix = f"_{sufix}"
    else:
        sufix = ""

    # === Boson (Double Crystal Ball) ===
    mean_boson = RooRealVar(
        f"resonant_background_model_ZG_mean_boson{sufix}",
        "mean{sufix}",
        91.1876,
        57.0,
        157.0,
        "GeV",
    )
    sigma_boson = RooRealVar(
        f"resonant_background_model_ZG_sigma_boson{sufix}",
        "sigma{sufix}",
        2.0,
        0.5,
        4.0,
        "GeV",
    )

    alpha1_upsilon = RooRealVar(
        f"resonant_background_model_ZG_alpha1_boson{sufix}",
        "alpha1{sufix}",
        3.0,
        0.0,
        10.0,
    )
    n1_upsilon = RooRealVar(
        f"resonant_background_model_ZG_n1_boson{sufix}", "n1{sufix}", 0.5, 0.1, 57.0
    )
    alpha2_upsilon = RooRealVar(
        f"resonant_background_model_ZG_alpha2_boson{sufix}",
        "alpha2{sufix}",
        3.0,
        0.0,
        10.0,
    )
    n2_upsilon = RooRealVar(
        f"resonant_background_model_ZG_n2_boson{sufix}", "n2{sufix}", 0.5, 0.1, 57.0
    )

    resonant_background_model_ZG_boson = RooDoubleCB(
        f"resonant_background_model_ZG_boson{sufix}",
        "Boson DoubleCB {sufix}",
        boson_mass,
        mean_boson,
        sigma_boson,
        alpha1_upsilon,
        n1_upsilon,
        alpha2_upsilon,
        n2_upsilon,
    )

    resonant_background_model_ZG_boson._keepalive = [
        mean_boson,
        sigma_boson,
        alpha1_upsilon,
        n1_upsilon,
        alpha2_upsilon,
        n2_upsilon,
    ]

    return resonant_background_model_ZG_boson


def resonant_background_modeling_Z(load_from_cache=False, plot_dir=PLOT_DIR, nbins=80):
    """
    Resonant Background Modeling
    """
    cached_resonant_background_model_Z_params = None

    if load_from_cache:
        with open("resonant_background_model_Z_params.json", "r") as f:
            cached_resonant_background_model_Z_params = json.load(f)

    input_file = "inputs/mass_Z_ZGTo2MuG_MMuMu-2To15_Run2.root"

    w = RooWorkspace("resonant_background_ws")

    w.factory("weight[-100,100]")
    w.factory(f"upsilon_mass[{UPSILON_MASS_LOWER}, {UPSILON_MASS_UPPER}]")
    w.factory(f"boson_mass[{BOSON_MASS_LOWER}, {BOSON_MASS_UPPER}]")

    # load data
    f = TFile.Open(input_file)
    data = RooDataSet(
        "resonant_background_data",
        "resonant_background_data",
        RooArgSet(w.var("boson_mass"), w.var("weight")),
        RooFit.Import(f.Events),
        RooFit.WeightVar(w.var("weight")),
    )
    getattr(w, "import")(data)

    resonant_background_model_Z = build_resonant_background_modeling_Z(
        w.var("boson_mass"),
    )
    getattr(w, "import")(resonant_background_model_Z)

    if cached_resonant_background_model_Z_params is not None:
        set_pdf_parameters(
            w.pdf("resonant_background_model_ZG_boson"),
            cached_resonant_background_model_Z_params,
            RooArgSet(w.var("boson_mass")),
            make_constant=True,
        )
    else:
        # fit to data
        fit_result = w.pdf("resonant_background_model_ZG_boson").fitTo(
            data, RooFit.Save(), RooFit.SumW2Error(True)
        )
        fit_result.Print("v")

    # plot data the and the pdf
    print("\n\n--> Saving plot ")
    w.var("boson_mass").SetTitle(r"m_{#mu#mu#gamma}")
    w.var("boson_mass").setUnit(r"GeV")
    w.var("upsilon_mass").SetTitle(r"m_{#mu#mu}")
    w.var("upsilon_mass").setUnit(r"GeV")

    # plot boson mass fit
    fastplot(
        w.pdf("resonant_background_model_ZG_boson"),
        w.data("resonant_background_data"),
        w.var("boson_mass"),
        f"{plot_dir}/resonant_background_fit_ZG_MC_m_mumugamma.pdf",
        # components=[
        #             (w.pdf("resonant_background_model_boson_cb"), 10),
        #             (w.pdf("resonant_background_model_boson_gauss"), 10),
        #             ],
        nbins=nbins,
        legend=[0.6, 0.6, 0.9, 0.92],  # legend=[0.6, 0.6, 0.93, 0.92], #type:ignore
        is_data=False,
        plot_range=get_signal_boson_plot_range("Z"),
    )

    print("data.sumEntries(): ", data.sumEntries())

    if cached_resonant_background_model_Z_params is not None:
        print("\n\n--> Loaded resonant Z parameters from cache")
        pprint(cached_resonant_background_model_Z_params)
        return cached_resonant_background_model_Z_params

    resonant_background_model_Z_params = get_pdf_parameters(
        w.pdf("resonant_background_model_ZG_boson"),
        RooArgSet(
            w.var("upsilon_mass"),
            w.var("boson_mass"),
        ),
    )
    print("\n\n--> Fit parameters ")
    pprint(resonant_background_model_Z_params)

    with open("resonant_background_model_Z_params.json", "w") as f:
        json.dump(resonant_background_model_Z_params, f, indent=4)

    return resonant_background_model_Z_params


@dataclass
class ControlRegionProps:
    name: str
    lower: float
    upper: float
    width: float = field(init=False)
    midpoint: float = field(init=False)

    def __post_init__(self):
        if self.upper <= self.lower:
            raise ValueError(
                "Invalid Control Region: Upper should be greater than lower"
            )

        self.width = self.upper - self.lower
        self.midpoint = (self.upper + self.lower) / 2.0


class ControlRegion(Enum):
    CR1 = ControlRegionProps(name="CR1", lower=4.0, upper=8.0)
    CR2 = ControlRegionProps(name="CR2", lower=12.0, upper=16.0)
    CR3 = ControlRegionProps(name="CR3", lower=16.0, upper=20.0)
    CR4 = ControlRegionProps(name="CR4", lower=20.0, upper=24.0)


def get_normalization_from_CR(
    boson_parameters,
    control_region: ControlRegion,
    load_from_cache=False,
    plot_dir=PLOT_DIR,
    nbins=40,
):
    """Resonant Background Modeling"""
    cached_norm_para = None
    if load_from_cache:
        with open(f"NormParams_{control_region.value.name}.json", "r") as f:
            cached_norm_para = json.load(f)

    input_file = "inputs/selected_Run2_resonant_background_modeling_MC_data_.root"

    w = RooWorkspace("resonant_background_ws")

    # load upsilon and boson fit parameters from json
    print("###########------------>>>>>>>>>>>  boson_parmeters:  ", boson_parameters)
    ### resonant model
    w.factory(
        "RooDoubleCB::resonant_background_model_res("
        f"boson_mass[{BOSON_MASS_LOWER}, {BOSON_MASS_UPPER}], "
        "" + str(boson_parameters["resonant_background_model_ZG_mean_boson"]) + ", "
        "" + str(boson_parameters["resonant_background_model_ZG_sigma_boson"]) + ", "
        "" + str(boson_parameters["resonant_background_model_ZG_alpha1_boson"]) + ", "
        "" + str(boson_parameters["resonant_background_model_ZG_n1_boson"]) + ", "
        "" + str(boson_parameters["resonant_background_model_ZG_alpha2_boson"]) + ", "
        "" + str(boson_parameters["resonant_background_model_ZG_n2_boson"]) + ""
        ")"
    )

    # Create parameters for the Johnson PDF
    meanJ = w.factory("mean[9.96, 0.01, 100]")
    sigmaJ = w.factory("lambda[24.3, 0.01, 100]")
    gammaJ = w.factory("gamma[-2.16, -10.0, 10.0]")
    deltaJ = w.factory("delta[1.42, 0.01, 20.0]")
    w.factory(
        "RooJohnson::resonant_background_model_non_res(boson_mass, mean, lambda, gamma, delta)"
    )

    # 1D model
    w.factory(
        "SUM::resonant_background_model(res_frac[0.5, 0., 1.]*resonant_background_model_res,resonant_background_model_non_res)"
    )

    w.factory("weight[-100,100]")

    w.factory(
        f"upsilon_mass[{RESONANT_CR_UPSILON_MASS_LOWER},{RESONANT_CR_UPSILON_MASS_UPPER}]"
    )

    # load data
    f = TFile.Open(input_file)
    data_ = RooDataSet(
        "resonant_background_data",
        "resonant_background_data",
        RooArgSet(w.var("boson_mass"), w.var("upsilon_mass"), w.var("weight")),
        RooFit.Import(f.Events),
        RooFit.WeightVar(w.var("weight")),
    )

    data = data_.reduce(
        f"(upsilon_mass < {control_region.value.upper} && upsilon_mass >= {control_region.value.lower} && boson_mass>{BOSON_MASS_LOWER} && boson_mass<{BOSON_MASS_UPPER})"
    )
    # data.Print()
    getattr(w, "import")(data)

    # fit to data
    _ = w.pdf("resonant_background_model").fitTo(data, RooFit.Save())

    # Get normalization
    resonant_background_model_integral = (
        w.pdf("resonant_background_model").createIntegral(data.get()).getVal()
    )
    print("#data: ", data.numEntries())
    print("#model: ", resonant_background_model_integral)
    print("res_frac: ", w.var("res_frac").getVal())
    normalization = w.var("res_frac").getVal() * (
        data.numEntries() / control_region.value.width
    )
    print("normalization: ", normalization)
    NormPara = {}
    NormPara["normalization"] = normalization
    NormPara["normalization_unc"] = (normalization) * TMath.Sqrt(
        (w.var("res_frac").getError() / w.var("res_frac").getVal()) ** 2
        + (TMath.Sqrt(data.numEntries()) / data.numEntries()) ** 2
    )

    print(f"-- > Normalization: {NormPara}")

    # plot data the and the pdf
    print("\n\n--> Saving plot for ")
    w.var("boson_mass").SetTitle(r"m_{#mu#mu#gamma}")
    w.var("boson_mass").setUnit(r"GeV")

    # plot boson mass fit
    fastplot(
        w.pdf("resonant_background_model"),
        w.data("resonant_background_data"),
        w.var("boson_mass"),
        f"{plot_dir}/resonant_background_MC_data_{control_region.value.name}_fit_m_mumugamma.pdf",
        components=[
            (w.pdf("resonant_background_model_res"), "Resonant component: RooDoubleCB"),
            (
                w.pdf("resonant_background_model_non_res"),
                "Non-resonant component: RooJohnson",
            ),
        ],
        nbins=nbins,
        legend=[0.6, 0.6, 0.93, 0.92],  # type: ignore
        is_data=True,
        model_legend_name="Total model: RooAddPdf",
    )

    if cached_norm_para is not None:
        print(
            f"-- > Loaded normalization from cache for {control_region.value.name}: {cached_norm_para}"
        )
        return cached_norm_para

    with open(f"NormParams_{control_region.value.name}.json", "w") as f:
        json.dump(NormPara, f, indent=4)

    return NormPara


def run_resonant_background(args: Namespace):
    os.makedirs(PLOT_DIR, exist_ok=True)

    nbins = getattr(args, "nbins", 60)
    load_from_cache = getattr(args, "use_cache", False)

    build_resonant_background_Higgs_ws(plot_dir=PLOT_DIR, nbins=nbins)
    build_resonant_background_Z_ws(plot_dir=PLOT_DIR, nbins=nbins)

    z_resonant_bkg_parameters = resonant_background_modeling_Z(
        load_from_cache=load_from_cache,
        plot_dir=PLOT_DIR,
        nbins=nbins,
    )

    normalizations_from_CR = [
        get_normalization_from_CR(
            z_resonant_bkg_parameters,
            control_region,
            load_from_cache=load_from_cache,
            plot_dir=PLOT_DIR,
            nbins=nbins,
        )
        for control_region in ControlRegion
    ]

    normalization_extrapolation = fit_and_plot(
        [c.value.midpoint for c in ControlRegion],
        [r["normalization"] for r in normalizations_from_CR],
        [r["normalization_unc"] for r in normalizations_from_CR],
        x0=10.0,
        output_pdf=f"{PLOT_DIR}/normalization_extrapolation.pdf",
        x_lines=[
            ControlRegion.CR1.value.upper,
            ControlRegion.CR2.value.lower,
            ControlRegion.CR2.value.upper,
            ControlRegion.CR3.value.upper,
        ],
        point_labels=["CR1", "CR2", "CR3", "CR4"],
    )
    print(f"Extrapolated normalization: {normalization_extrapolation}")

    return {
        "z_boson_params": z_resonant_bkg_parameters,
        "normalizations_from_cr": normalizations_from_CR,
        "normalization_extrapolation": normalization_extrapolation,
    }
