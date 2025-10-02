from enum import IntEnum, auto
from pprint import pprint

from ROOT import (  # type: ignore
    RooArgSet,  # type: ignore
    RooDataSet,  # type: ignore
    RooFit,  # type: ignore
    RooWorkspace,  # type: ignore
    TFile,  # type: ignore
)

from .fastplot import fastplot
from .ws_helper import *


def resonant_background_modeling_Higgs():
    """
    Resonant Background Modeling
    """

    input_file = "inputs/mass_H_GluGluHToMuMuG_M125_MLL-0To60_Dalitz_012j_Run2.root"

    w = RooWorkspace("resonant_background_ws")

    ### boson model

    w.factory(
        "RooDoubleCB::resonant_background_model_boson("
        "boson_mass[75, 200], "
        "mean_boson[125, 75, 200], "
        "sigma_boson[2, 0.5, 4],"
        "alpha1_upsilon[3, 0, 10],"
        "n1_upsilon[0.5, 0.1, 50],"
        "alpha2_upsilon[3, 0, 10],"
        "n2_upsilon[0.5, 0.1, 50]"
        ")"
    )

    # upsilon model
    w.factory(
        "RooBernstein::resonant_background_model_upsilon(upsilon_mass[8, 12],{1,p1[5, 0, 10]})"
    )

    # 2D model
    w.factory(
        "PROD::resonant_background_model(resonant_background_model_boson,resonant_background_model_upsilon)"
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
    getattr(w, "import")(data)

    # fit to data
    fit_result = w.pdf("resonant_background_model").fitTo(data, RooFit.Save())

    # plot data the and the pdf
    print("\n\n--> Saving plot ")
    nBins = 80
    w.var("boson_mass").SetTitle(r"m_{#mu#mu#gamma}")
    w.var("boson_mass").setUnit(r"GeV")
    w.var("upsilon_mass").SetTitle(r"m_{#mu#mu}")
    w.var("upsilon_mass").setUnit(r"GeV")
    # w.pdf("resonant_background_model_boson_gauss").SetTitle(r"Gaussian Component")
    # w.pdf("resonant_background_model_boson_cb").SetTitle(r"CB Component")
    # w.pdf("signal_model_upsilon").SetTitle(r"CB Component")

    # plot boson mass fit
    fastplot(
        w.pdf("resonant_background_model"),
        w.data("resonant_background_data"),
        w.var("boson_mass"),
        "plots/fit_2d_data/resonant_background_fit_HiggsDalitz_MC_m_mumugamma.pdf",
        nbins=nBins,
        legend=[0.6, 0.6, 0.9, 0.92],  # legend=[0.2, 0.6, 0.5, 0.92], #type: ignore
        is_data=False,
    )

    # plot upsilon mass fit
    fastplot(
        w.pdf("resonant_background_model"),
        w.data("resonant_background_data"),
        w.var("upsilon_mass"),
        "plots/fit_2d_data/resonant_background_fit_HiggsDalitz_MC_m_mumu.pdf",
        # components=[
        #             (upsilon_1S, 10),
        #             (upsilon_2S, 10),
        #             (upsilon_3S, 10),
        #             (background, 20)
        #             ],
        nbins=nBins,
        legend=[0.2, 0.7, 0.53, 0.9],  # legend=[0.6, 0.2, 0.93, 0.4], #type: ignore
        is_data=False,
    )

    # setting all var as constants
    w = set_constant(w)

    print("\n\n--> Fit parameters ")
    boson_parmeters = {}
    boson_parmeters["mean_bosonH"] = (
        fit_result.floatParsFinal().find("mean_boson").getValV()
    )
    boson_parmeters["sigma_bosonH"] = (
        fit_result.floatParsFinal().find("sigma_boson").getValV()
    )
    boson_parmeters["alpha1_upsilonH"] = (
        fit_result.floatParsFinal().find("alpha1_upsilon").getValV()
    )
    boson_parmeters["n1_upsilonH"] = (
        fit_result.floatParsFinal().find("n1_upsilon").getValV()
    )
    boson_parmeters["alpha2_upsilonH"] = (
        fit_result.floatParsFinal().find("alpha2_upsilon").getValV()
    )
    boson_parmeters["n2_upsilonH"] = (
        fit_result.floatParsFinal().find("n2_upsilon").getValV()
    )

    pprint(boson_parmeters)

    print("data.sumEntries(): ", data.sumEntries())

    return w


def resonant_background_modeling_Z():
    ################################
    # Resonant Background Modeling #
    ################################

    input_file = "inputs/mass_Z_ZGTo2MuG_MMuMu-2To15_Run2.root"

    w = RooWorkspace("resonant_background_ws")

    ### boson model

    w.factory(
        "RooDoubleCB::resonant_background_model_boson("
        "boson_mass[75, 200], "
        "mean_boson[91.1876, 70, 150], "
        "sigma_boson[2, 0.5, 4],"
        "alpha1_upsilon[3, 0, 10],"
        "n1_upsilon[0.5, 0.1, 50],"
        "alpha2_upsilon[3, 0, 10],"
        "n2_upsilon[0.5, 0.1, 50]"
        ")"
    )

    # upsilon model
    w.factory(
        "RooBernstein::resonant_background_model_upsilon( upsilon_mass[8, 12], { 1, p1[5, 0, 10] })"
    )

    # 2D model
    w.factory(
        "PROD::resonant_background_model(resonant_background_model_boson,resonant_background_model_upsilon)"
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
    getattr(w, "import")(data)

    # fit to data
    fit_result = w.pdf("resonant_background_model").fitTo(data, RooFit.Save())

    ##################################################################################################################################################################
    # Define the ranges
    left_range = RooFit.Range("LEFT")
    center_range = RooFit.Range("CENTER")
    right_range = RooFit.Range("RIGHT")
    full_range = RooFit.Range("FULL")

    # Get the normalization sets for the ranges
    left_norm_set = RooArgSet(w.var("boson_mass"))
    center_norm_set = RooArgSet(w.var("boson_mass"))
    right_norm_set = RooArgSet(w.var("boson_mass"))
    full_norm_set = RooArgSet(w.var("boson_mass"))

    # Set the ranges for the normalization sets
    w.var("boson_mass").setRange("LEFT", 75, 80)
    w.var("boson_mass").setRange("CENTER", 110, 115)
    w.var("boson_mass").setRange("RIGHT", 135, 200)
    w.var("boson_mass").setRange("FULL", 75, 200)

    # Calculate the integrals for each range
    integral_left = (
        w.pdf("resonant_background_model_boson")
        .createIntegral(left_norm_set, left_range)
        .getVal()
    )
    integral_center = (
        w.pdf("resonant_background_model_boson")
        .createIntegral(center_norm_set, center_range)
        .getVal()
    )
    integral_right = (
        w.pdf("resonant_background_model_boson")
        .createIntegral(right_norm_set, right_range)
        .getVal()
    )
    integral_full = (
        w.pdf("resonant_background_model_boson")
        .createIntegral(full_norm_set, full_range)
        .getVal()
    )

    # Print the results
    print("Integral in LEFT range:", integral_left)
    print("Integral in CENTER range:", integral_center)
    print("Integral in RIGHT range:", integral_right)
    print("Integral in FULL range:", integral_full)

    ##################################################################################################################################################################

    # plot data the and the pdf
    print("\n\n--> Saving plot ")
    nBins = 80
    w.var("boson_mass").SetTitle(r"m_{#mu#mu#gamma}")
    w.var("boson_mass").setUnit(r"GeV")
    w.var("upsilon_mass").SetTitle(r"m_{#mu#mu}")
    w.var("upsilon_mass").setUnit(r"GeV")
    # w.pdf("resonant_background_model_boson_gauss").SetTitle(r"Gaussian Component")
    # w.pdf("resonant_background_model_boson_cb").SetTitle(r"CB Component")
    # w.pdf("signal_model_upsilon").SetTitle(r"CB Component")

    # plot boson mass fit
    fastplot(
        w.pdf("resonant_background_model"),
        w.data("resonant_background_data"),
        w.var("boson_mass"),
        "plots/fit_2d_data/resonant_background_fit_ZG_MC_m_mumugamma.pdf",
        # components=[
        #             (w.pdf("resonant_background_model_boson_cb"), 10),
        #             (w.pdf("resonant_background_model_boson_gauss"), 10),
        #             ],
        nbins=nBins,
        legend=[0.6, 0.6, 0.9, 0.92],  # legend=[0.6, 0.6, 0.93, 0.92], #type:ignore
        is_data=False,
    )

    # plot upsilon mass fit
    fastplot(
        w.pdf("resonant_background_model"),
        w.data("resonant_background_data"),
        w.var("upsilon_mass"),
        "plots/fit_2d_data/resonant_background_fit_ZG_MC_m_mumu.pdf",
        # components=[
        #             (upsilon_1S, 10),
        #             (upsilon_2S, 10),
        #             (upsilon_3S, 10),
        #             (background, 20)
        #             ],
        nbins=nBins,
        legend=[0.6, 0.2, 0.93, 0.4],  # type:ignore
        is_data=False,
    )

    # setting all var as constants
    w = set_constant(w)

    print("\n\n--> Fit parameters ")
    boson_parmeters = {}
    boson_parmeters["mean_boson"] = (
        fit_result.floatParsFinal().find("mean_boson").getValV()
    )  # mean_boson.getValV()
    boson_parmeters["sigma_boson"] = (
        fit_result.floatParsFinal().find("sigma_boson").getValV()
    )  # sigma_boson.getValV()
    boson_parmeters["alpha1_upsilon"] = (
        fit_result.floatParsFinal().find("alpha1_upsilon").getValV()
    )  # alpha1_upsilon.getValV()
    boson_parmeters["n1_upsilon"] = (
        fit_result.floatParsFinal().find("n1_upsilon").getValV()
    )  # n1_upsilon.getValV()
    boson_parmeters["alpha2_upsilon"] = (
        fit_result.floatParsFinal().find("alpha2_upsilon").getValV()
    )  # alpha2_upsilon.getValV()
    boson_parmeters["n2_upsilon"] = (
        fit_result.floatParsFinal().find("n2_upsilon").getValV()
    )  # n2_upsilon.getValV()

    pprint(boson_parmeters)

    print("data.sumEntries(): ", data.sumEntries())

    return w


class ControlRegion(IntEnum):
    CR1 = auto()
    CR2 = auto()
    CR3 = auto()
    CR4 = auto()


def resonant_background_modeling_CR(control_region: ControlRegion):
    print(control_region)
