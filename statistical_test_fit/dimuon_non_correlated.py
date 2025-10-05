from pprint import pprint
from ROOT import (  # type: ignore
    RooAddPdf,  # type: ignore
    RooArgList,  # type: ignore
    RooArgSet,  # type: ignore
    RDataFrame,  # type: ignore
    RooCBShape,  # type:ignore
    RooChebychev,  # type: ignore
    RooDataSet,  # type: ignore
    RooFit,  # type: ignore
    RooFormulaVar,  # type: ignore
    RooRealVar,  # type: ignore
    TFile,  # type: ignore
    kTRUE,  # type: ignore
)

from .fastplot import fastplot
from .ws_helper import *


def build_upsilon_model(mass, sufix=None):
    if sufix is not None:
        sufix = f"_{sufix}"
    else:
        sufix = ""

    # set range of observable for MuMu Invariant mass
    M1S = 9.46  # upsilon 1S pgd mass value
    M2S = 10.02  # upsilon 2S pgd mass value
    M3S = 10.35  # upsilon 3S pgd mass value

    # define parameters
    # mean
    m1S = RooRealVar(f"m1S{sufix}", "m1S", M1S, "GeV")
    m2S = RooRealVar(f"m2S{sufix}", "m2S", M2S, "GeV")
    m3S = RooRealVar(f"m3S{sufix}", "m2S", M3S, "GeV")
    # mscale = RooRealVar("mscale", "mass scale factor", 1, 0, 1.5)
    mscale = RooRealVar(f"mscale{sufix}", "mass scale factor", 1, 0, 1.2)

    mean1S = RooFormulaVar(f"mean1S{sufix}", "@0*@1", RooArgList(m1S, mscale))
    mean2S = RooFormulaVar(f"mean2S{sufix}", "@0*@1", RooArgList(m2S, mscale))
    mean3S = RooFormulaVar(f"mean3S{sufix}", "@0*@1", RooArgList(m3S, mscale))

    # Signal + background  Upsilon Model
    sigma1S = RooRealVar(f"sigma1S{sufix}", "sigma1S", 0.1, 0, 1.2)
    sigma2S = RooFormulaVar(
        f"sigma2S{sufix}", "@0*@1/@2", RooArgList(sigma1S, m2S, m1S)
    )
    sigma3S = RooFormulaVar(
        f"sigma3S{sufix}", "@0*@1/@2", RooArgList(sigma1S, m3S, m1S)
    )

    alpha1_upsilon = RooRealVar(f"alpha1_upsilon{sufix}", "alpha1_upsilon", 1, 0, 10)
    alpha2_upsilon = RooRealVar(f"alpha2_upsilon{sufix}", "alpha2_upsilon", 1, 0, 10)
    alpha3_upsilon = RooRealVar(f"alpha3_upsilon{sufix}", "alpha3_upsilon", 1, 0, 10)

    n1_upsilon = RooRealVar(f"n1_upsilon{sufix}", "n1_upsilon", 3, 1, 10)
    n2_upsilon = RooRealVar(f"n2_upsilon{sufix}", "n2_upsilon", 3, 1, 10)
    n3_upsilon = RooRealVar(f"n3_upsilon{sufix}", "n3_upsilon", 3, 1, 10)

    # Upsilon 1S Model
    #  upsilon_1S = RooGaussian("upsilon_1S", "#Upsilon(1S) Model", mass, mean1S, sigma1S)
    upsilon_1S = RooCBShape(
        f"upsilon_1S{sufix}",
        "#Upsilon(1S) Model",
        mass,
        mean1S,
        sigma1S,
        alpha1_upsilon,
        n1_upsilon,
    )

    # Upsilon 2S Model
    #  upsilon_2S = RooGaussian("upsilon_2S", "#Upsilon(2S) Model", mass, mean2S, sigma2S)
    upsilon_2S = RooCBShape(
        f"upsilon_2S{sufix}",
        "#Upsilon(2S) Model",
        mass,
        mean2S,
        sigma2S,
        alpha2_upsilon,
        n2_upsilon,
    )

    # Upsilon 3S Model
    #  upsilon_3S = RooGaussian("upsilon_3S", "#Upsilon(3S) Model", mass, mean3S, sigma3S)
    upsilon_3S = RooCBShape(
        f"upsilon_3S{sufix}",
        "#Upsilon(3S) Model",
        mass,
        mean3S,
        sigma3S,
        alpha3_upsilon,
        n3_upsilon,
    )

    # Build total model
    upsilon_1s_frac = RooRealVar(
        f"upsilon_1s_frac{sufix}", "upsilon_1s_frac", 0.2, 0, 1
    )
    upsilon_2s_frac = RooRealVar(
        f"upsilon_2s_frac{sufix}", "upsilon_2s_frac", 0.2, 0, 1
    )
    upsilon_3s_frac = RooRealVar(
        f"upsilon_3s_frac{sufix}", "upsilon_3s_frac", 0.2, 0, 1
    )

    upsilon_model = RooAddPdf(
        f"upsilon_model{sufix}",
        "upsilon_model",
        RooArgList(upsilon_1S, upsilon_2S, upsilon_3S),
        RooArgList(upsilon_1s_frac, upsilon_2s_frac),
        kTRUE,
    )

    upsilon_model._keepalive = [
        m1S,
        m2S,
        m3S,
        mscale,
        mean1S,
        mean2S,
        mean3S,
        sigma1S,
        sigma2S,
        sigma3S,
        alpha1_upsilon,
        alpha2_upsilon,
        alpha3_upsilon,
        n1_upsilon,
        n2_upsilon,
        n3_upsilon,
        upsilon_1S,
        upsilon_2S,
        upsilon_3S,
        upsilon_1s_frac,
        upsilon_2s_frac,
        upsilon_3s_frac,
    ]

    return upsilon_model, upsilon_1S, upsilon_2S, upsilon_3S


def dimuon_non_correlated(m_mumu_lower=8.0, m_mumu_upper=12.0):
    """Upsilon + Background Model"""

    lowRange = m_mumu_lower
    highRange = m_mumu_upper

    # set mass observble
    #  mass = ROOT.RooRealVar("dimuon_mass_for_upsilon_fit","m_{#mu#mu}",lowRange,highRange, "GeV")
    mass = RooRealVar("upsilon_mass", "m_{#mu#mu}", lowRange, highRange, "GeV")
    mass.setBins(60)

    (
        upsilon_model,
        upsilon_1S,
        upsilon_2S,
        upsilon_3S,
    ) = build_upsilon_model(mass)

    # Background Model
    bkg_a1 = RooRealVar("bkg_a1", "bkg_a1", 0, -1, 1)
    bkg_a2 = RooRealVar("bkg_a2", "bkg_a2", 0, -1, 1)
    # background  = RooChebychev("background","background", mass, RooArgList(bkg_a1, bkg_a2))
    background = RooChebychev(
        "background", "Background Model", mass, RooArgList(bkg_a1)
    )

    upsilon_frac = RooRealVar("upsilon_frac", "upsilon_frac", 0.2, 0, 1)

    total_model = RooAddPdf(
        "total_model",
        "total_model",
        RooArgList(upsilon_model, background),
        RooArgList(upsilon_frac),
        kTRUE,
    )

    in_file = "inputs/selected_Run2_dimuon_non_correlated.root"
    tree_name = "dimuons_masses"
    old_name = "mass"
    new_name = "upsilon_mass"
    out_file = "inputs/selected_Run2_dimuon_non_correlated_renamed_branch.root"

    # Build a dataframe from the original tree
    df = RDataFrame(tree_name, in_file)
    cols = [str(c) for c in df.GetColumnNames()]
    df2 = df.Define(new_name, old_name)
    snapshot_cols = [new_name if c == old_name else c for c in cols if c != new_name]
    df2.Snapshot(tree_name, out_file, snapshot_cols)

    # load data
    f = TFile.Open("inputs/selected_Run2_dimuon_non_correlated_renamed_branch.root")
    # f = TFile.Open("inputs/selected_Run2_dimuon_non_correlated.root")
    f.dimuons_masses.Print()
    # data = RooDataSet("data", "data", f.dimuons_masses, RooArgSet(mass), a)
    data = RooDataSet(
        "data",
        "data",
        RooArgSet(mass),
        RooFit.Import(f.dimuons_masses),
        RooFit.Cut(f"(upsilon_mass < 12.0 && upsilon_mass >= 8.0)"),
    )

    # perform the fit
    _ = total_model.fitTo(data, RooFit.Save())

    # close file
    f.Close()

    fastplot(
        total_model,
        data,
        mass,
        "plots/fit_2d_data/dimuon_mass_fit_Run2.pdf",
        components=[
            (upsilon_1S, 10),
            (upsilon_2S, 10),
            (upsilon_3S, 10),
            (background, 20),
        ],
        nbins=60,
        legend=[0.6, 0.6, 0.93, 0.92],  # type: ignore
    )

    upsilon_model_params = get_pdf_parameters(upsilon_model, RooArgSet(mass))

    return upsilon_model_params
