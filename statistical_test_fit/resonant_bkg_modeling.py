from enum import Enum, auto
from pprint import pprint

from ROOT import (  # type: ignore
    RooArgSet,  # type: ignore
    RooDataSet,  # type: ignore
    RooFit,  # type: ignore
    RooWorkspace,  # type: ignore
    TFile,  # type: ignore
    TMath,  # type: ignore
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
    """
    Resonant Background Modeling
    """

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
    boson_parameters = {}
    boson_parameters["mean_boson"] = (
        fit_result.floatParsFinal().find("mean_boson").getValV()
    )
    boson_parameters["sigma_boson"] = (
        fit_result.floatParsFinal().find("sigma_boson").getValV()
    )
    boson_parameters["alpha1_upsilon"] = (
        fit_result.floatParsFinal().find("alpha1_upsilon").getValV()
    )
    boson_parameters["n1_upsilon"] = (
        fit_result.floatParsFinal().find("n1_upsilon").getValV()
    )
    boson_parameters["alpha2_upsilon"] = (
        fit_result.floatParsFinal().find("alpha2_upsilon").getValV()
    )
    boson_parameters["n2_upsilon"] = (
        fit_result.floatParsFinal().find("n2_upsilon").getValV()
    )

    pprint(boson_parameters)

    print("data.sumEntries(): ", data.sumEntries())

    return w, boson_parameters


class ControlRegion(Enum):
    CR1 = "CR1"
    CR2 = "CR2"
    CR3 = "CR3"
    CR4 = "CR4"


def get_interval(control_region):
    if control_region == ControlRegion.CR1:
        return 5, 8, 3
    if control_region == ControlRegion.CR2:
        return 5, 8, 3
    if control_region == ControlRegion.CR3:
        return 5, 8, 3
    if control_region == ControlRegion.CR4:
        return 5, 8, 3


def get_normalization_from_CR(boson_parameters, control_region: ControlRegion):
    #  """
    # Resonant Background Modeling
    #  """
    lower, upper, width = get_interval(control_region)

    input_file = "inputs/selected_Run2_resonant_background_modeling_MC_data_.root"

    w = RooWorkspace("resonant_background_ws")

    # load upsilon and boson fit parameters from json

    print("###########------------>>>>>>>>>>>  boson_parmeters:  ", boson_parameters)
    ### resonant model
    w.factory(
        "RooDoubleCB::resonant_background_model_res("
        "boson_mass[75, 200], "
        "" + str(boson_parameters["mean_boson"]) + ", "
        "" + str(boson_parameters["sigma_boson"]) + ", "
        "" + str(boson_parameters["alpha1_upsilon"]) + ", "
        "" + str(boson_parameters["n1_upsilon"]) + ", "
        "" + str(boson_parameters["alpha2_upsilon"]) + ", "
        "" + str(boson_parameters["n2_upsilon"]) + ""
        ")"
    )

    # non resonant model
    w.factory(
        "RooBernstein::resonant_background_model_non_res("
        "boson_mass[75, 200],"
        "{"
        "1,"
        "p1[0,-1,1],"
        "p2[0,-1,1]"
        "}"
        ")"
    )

    # w.factory("RooExponential::resonant_background_model_non_res("
    #          "boson_mass[70, 150],"
    #          "alpha[0, -10, 10]"
    #         ")")

    # 1D model
    w.factory(
        "SUM::resonant_background_model(res_frac[0.5, 0., 1.]*resonant_background_model_res,resonant_background_model_non_res)"
    )

    w.factory("weight[-100,100]")

    w.factory("upsilon_mass[4.,35.]")

    # load data
    f = TFile.Open(input_file)
    data1 = RooDataSet(
        "resonant_background_data",
        "resonant_background_data",
        RooArgSet(w.var("boson_mass"), w.var("upsilon_mass"), w.var("weight")),
        RooFit.Import(f.Events),
        RooFit.WeightVar(w.var("weight")),
    )

    data = data1.reduce(
        f"(upsilon_mass < {upper} && upsilon_mass >= {lower} && boson_mass>75. && boson_mass<200.)"
    )
    # data.Print()
    # data = data2.reduce(ROOT.RooArgSet("boson_mass", "weight"))
    getattr(w, "import")(data)

    # fit to data
    fit_result = w.pdf("resonant_background_model").fitTo(data, RooFit.Save())

    # Get normalization
    resonant_background_model_integral = (
        w.pdf("resonant_background_model").createIntegral(data.get()).getVal()
    )
    print("#data: ", data.numEntries())
    print("#model: ", resonant_background_model_integral)
    print("res_frac: ", w.var("res_frac").getVal())
    normalization = w.var("res_frac").getVal() * (data.numEntries() / width)
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
    nBins = 16
    w.var("boson_mass").SetTitle(r"m_{#mu#mu#gamma}")
    w.var("boson_mass").setUnit(r"GeV")

    # w.pdf("resonant_background_model_boson_gauss").SetTitle(r"Gaussian Component")
    # w.pdf("resonant_background_model_boson_cb").SetTitle(r"CB Component")
    # w.pdf("signal_model_upsilon").SetTitle(r"CB Component")

    # plot boson mass fit
    fastplot(
        w.pdf("resonant_background_model"),
        w.data("resonant_background_data"),
        w.var("boson_mass"),
        "plots_2d_data/resonant_background_MC_data_{control_region}_fit_m_mumugamma.pdf",
        components=[
            (w.pdf("resonant_background_model_res"), 10),
            (w.pdf("resonant_background_model_non_res"), 10),
        ],
        nbins=nBins,
        legend=[0.6, 0.6, 0.93, 0.92],  # type: ignore
        is_data=True,
    )

    # setting all var as constants
    w = set_constant(w)

    # Save the workspace into a ROOT file
    print("\n\n--> Parameters for ")
    fit_result.Print()

    return NormPara
