from ROOT import (  # type: ignore
    RooAddPdf,  # type: ignore
    RooArgList,  # type: ignore
    RooArgSet,  # type: ignore
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


def dimuon_non_correlated(m_mumu_lower=8.0, m_mumu_upper=12.0):
    """Upsilon + Background Model"""

    # set range of observable for MuMu Invariant mass
    M1S = 9.46  # upsilon 1S pgd mass value
    M2S = 10.02  # upsilon 2S pgd mass value
    M3S = 10.35  # upsilon 3S pgd mass value

    lowRange = m_mumu_lower
    highRange = m_mumu_upper

    # set mass observble
    #  mass = ROOT.RooRealVar("dimuon_mass_for_upsilon_fit","m_{#mu#mu}",lowRange,highRange, "GeV")
    mass = RooRealVar("mass", "m_{#mu#mu}", lowRange, highRange, "GeV")
    mass.setBins(60)
    # define parameters
    # mean
    m1S = RooRealVar("m1S", "m1S", M1S, "GeV")
    m2S = RooRealVar("m2S", "m2S", M2S, "GeV")
    m3S = RooRealVar("m3S", "m2S", M3S, "GeV")
    # mscale = RooRealVar("mscale", "mass scale factor", 1, 0, 1.5)
    mscale = RooRealVar("mscale", "mass scale factor", 1, 0, 1.2)

    mean1S = RooFormulaVar("mean1S", "@0*@1", RooArgList(m1S, mscale))
    mean2S = RooFormulaVar("mean2S", "@0*@1", RooArgList(m2S, mscale))
    mean3S = RooFormulaVar("mean3S", "@0*@1", RooArgList(m3S, mscale))

    # Signal + background  Upsilon Model
    sigma1S = RooRealVar("sigma1S", "sigma1S", 0.1, 0, 1.2)
    sigma2S = RooFormulaVar("sigma2S", "@0*@1/@2", RooArgList(sigma1S, m2S, m1S))
    sigma3S = RooFormulaVar("sigma3S", "@0*@1/@2", RooArgList(sigma1S, m3S, m1S))

    alpha1_upsilon = RooRealVar("alpha1_upsilon", "alpha1_upsilon", 1, 0, 10)
    alpha2_upsilon = RooRealVar("alpha2_upsilon", "alpha2_upsilon", 1, 0, 10)
    alpha3_upsilon = RooRealVar("alpha3_upsilon", "alpha3_upsilon", 1, 0, 10)

    n1_upsilon = RooRealVar("n1_upsilon", "n1_upsilon", 3, 1, 10)
    n2_upsilon = RooRealVar("n2_upsilon", "n2_upsilon", 3, 1, 10)
    n3_upsilon = RooRealVar("n3_upsilon", "n3_upsilon", 3, 1, 10)

    # Upsilon 1S Model
    #  upsilon_1S = RooGaussian("upsilon_1S", "#Upsilon(1S) Model", mass, mean1S, sigma1S)
    upsilon_1S = RooCBShape(
        "upsilon_1S",
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
        "upsilon_2S",
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
        "upsilon_3S",
        "#Upsilon(3S) Model",
        mass,
        mean3S,
        sigma3S,
        alpha3_upsilon,
        n3_upsilon,
    )

    # Background Model
    bkg_a1 = RooRealVar("bkg_a1", "bkg_a1", 0, -1, 1)
    bkg_a2 = RooRealVar("bkg_a2", "bkg_a2", 0, -1, 1)
    # background  = RooChebychev("background","background", mass, RooArgList(bkg_a1, bkg_a2))
    background = RooChebychev(
        "background", "Background Model", mass, RooArgList(bkg_a1)
    )

    # Build total model
    upsilon_1s_frac = RooRealVar("upsilon_1s_frac", "upsilon_1s_frac", 0.2, 0, 1)
    upsilon_2s_frac = RooRealVar("upsilon_2s_frac", "upsilon_2s_frac", 0.2, 0, 1)
    upsilon_3s_frac = RooRealVar("upsilon_3s_frac", "upsilon_3s_frac", 0.2, 0, 1)
    upsilon_frac = RooRealVar("upsilon_frac", "upsilon_frac", 0.2, 0, 1)

    upsilon_model = RooAddPdf(
        "upsilon_model",
        "upsilon_model",
        RooArgList(upsilon_1S, upsilon_2S, upsilon_3S),
        RooArgList(upsilon_1s_frac, upsilon_2s_frac),
        kTRUE,
    )

    total_model = RooAddPdf(
        "total_model",
        "total_model",
        RooArgList(upsilon_model, background),
        RooArgList(upsilon_frac),
        kTRUE,
    )

    # load data
    f = TFile.Open("inputs/selected_Run2_dimuon_non_correlated.root")
    data = RooDataSet("data", "data", f.dimuons_masses, RooArgSet(mass))

    # perform the fit
    _ = total_model.fitTo(data, RooFit.Save())

    # close file
    f.Close()

    print("----->>> upsilon_1s_frac.getValV():  ", upsilon_1s_frac.getValV())
    print("----->>> upsilon_2s_frac.getValV():  ", upsilon_2s_frac.getValV())
    print("----->>> upsilon_3s_frac.getValV():  ", upsilon_3s_frac.getValV())
    print("----->>> background.getValV():       ", background.getValV())

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

    # freeze models before return
    observables = upsilon_model.getObservables(data)
    _ = freeze_pdf_params(upsilon_model, observables)

    upsilon_model._keepalive = [
        mass,
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

    return upsilon_model
