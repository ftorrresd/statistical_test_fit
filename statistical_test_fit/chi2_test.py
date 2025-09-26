from dataclasses import dataclass

from ROOT import (
    RooAbsReal,  # type: ignore
    RooArgList,  # type: ignore
    RooArgSet,  # type: ignore
    RooFit,  # type: ignore
    TCanvas,  # type: ignore
    TMath,  # type: ignore
)

from .bkg_model import BkgPdfFamily


@dataclass
class ChiSquareResult:
    chi2: float
    pvalue: float
    ndf: int
    chi2_ndf: float

    @staticmethod
    def compute_chi_square_1d(
        model,
        data,
        x,
        outprefix,
        pdf_family: BkgPdfFamily,
        region_name="left,middle,right",
        nbins=40,
        # norm_range="left,middle,right",
        nfloatpars=None,
    ) -> "ChiSquareResult":
        def _count_float_pars(model, x):
            """Return number of non-constant parameters the model depends on."""
            params_set = model.getParameters(RooArgSet(x))  # RooArgSet
            params_list = RooArgList(params_set)  # make it indexable
            nfree = 0
            for i in range(params_list.getSize()):
                p = params_list[i]
                if hasattr(p, "isConstant") and not p.isConstant():
                    nfree += 1
            return nfree

        # Determine number of floated parameters to subtract in ndf
        if nfloatpars is None:
            nfloatpars = _count_float_pars(model, x)

        # Build a frame and plot only the region of interest
        frame = x.frame(RooFit.Bins(nbins), RooFit.Title(f"GOF in {region_name}"))
        data_name = f"data_{region_name}"
        curve_name = f"pdf_{region_name}"

        # Data: only points inside the region (binned into nbins)
        data.plotOn(
            frame,
            RooFit.CutRange(region_name),
            RooFit.Binning(nbins),
            RooFit.Name(data_name),
        )

        # PDF: evaluated only in the region, but normalized over norm_range
        nData = data.sumEntries("", region_name)
        model.plotOn(
            frame,
            RooFit.Normalization(nData, RooAbsReal.NumEvent),
            RooFit.Name(curve_name),
        )

        # chi2/ndf using the named objects and subtracting nfloatpars
        chi2_ndf = frame.chiSquare(curve_name, data_name, nfloatpars)

        # Recover number of points to get chi2 and ndf
        # (RooHist derives from TGraph --> GetN() is available)
        rhist = frame.findObject(data_name)
        npoints = int(rhist.GetN()) if rhist else 0

        ndf = max(1, npoints - int(nfloatpars))  # guard against <=0

        chi2 = chi2_ndf * ndf

        # p-value for chi2 with ndf dof: upper tail integral
        pval = TMath.Prob(chi2, ndf)

        can = TCanvas("c", "c", 820, 920)
        frame.Draw()

        can.SaveAs(f"plots/fit_1d/control_{outprefix}_{pdf_family}.pdf")

        return ChiSquareResult(chi2=chi2, pvalue=pval, ndf=ndf, chi2_ndf=chi2_ndf)

    @staticmethod
    def compute_chi_square_2d(
        model,
        data,
        x,
        y,
        outprefix,
        pdf_family: BkgPdfFamily,
        region_name="left,middle,right",
        nbins=40,
        # norm_range="left,middle,right",
        nfloatpars=None,
    ) -> "ChiSquareResult":
        def _count_float_pars(model, x):
            """Return number of non-constant parameters the model depends on."""
            params_set = model.getParameters(RooArgSet(x))  # RooArgSet
            params_list = RooArgList(params_set)  # make it indexable
            nfree = 0
            for i in range(params_list.getSize()):
                p = params_list[i]
                if hasattr(p, "isConstant") and not p.isConstant():
                    nfree += 1
            return nfree

        # Determine number of floated parameters to subtract in ndf
        if nfloatpars is None:
            nfloatpars = _count_float_pars(model, x)

        # Build a frame and plot only the region of interest
        frame = y.frame(RooFit.Bins(nbins), RooFit.Title(f"GOF on Y"))
        data_name = f"data_{region_name}"
        curve_name = f"pdf_{region_name}"

        # Data: only points inside the region (binned into nbins)
        data.plotOn(
            frame,
            # RooFit.CutRange(region_name),
            RooFit.Binning(nbins),
            RooFit.Name(data_name),
        )

        # PDF: evaluated only in the region, but normalized over norm_range
        # nData = data.sumEntries("", region_name)
        nData = data.sumEntries("")
        model.plotOn(
            frame,
            # RooFit.ProjWData(RooArgSet(x), data),
            # RooFit.Normalization(nData, RooAbsReal.NumEvent),
            RooFit.NormRange("LEFT,MIDDLE,RIGHT"),
            RooFit.Name(curve_name),
        )

        # chi2/ndf using the named objects and subtracting nfloatpars
        chi2_ndf = frame.chiSquare(curve_name, data_name, nfloatpars)

        # Recover number of points to get chi2 and ndf
        # (RooHist derives from TGraph --> GetN() is available)
        rhist = frame.findObject(data_name)
        npoints = int(rhist.GetN()) if rhist else 0

        ndf = max(1, npoints - int(nfloatpars))  # guard against <=0

        chi2 = chi2_ndf * ndf

        # p-value for chi2 with ndf dof: upper tail integral
        pval = TMath.Prob(chi2, ndf)

        can = TCanvas("c", "c", 820, 920)
        frame.Draw()

        can.SaveAs(f"plots/fit_2d/control_mumugamma_{outprefix}_{pdf_family}.pdf")

        res = ChiSquareResult(chi2=chi2, pvalue=pval, ndf=ndf, chi2_ndf=chi2_ndf)

        # plot on X

        # Build a frame and plot only the region of interest
        frame = x.frame(RooFit.Bins(nbins), RooFit.Title(f"GOF on X"))
        data_name = f"data_{region_name}"
        curve_name = f"pdf_{region_name}"

        # Data: only points inside the region (binned into nbins)
        data.plotOn(
            frame,
            # RooFit.CutRange(region_name),
            RooFit.Binning(nbins),
            RooFit.Name(data_name),
        )

        # PDF: evaluated only in the region, but normalized over norm_range
        # nData = data.sumEntries("", region_name)
        nData = data.sumEntries("")
        model.plotOn(
            frame,
            # RooFit.ProjWData(RooArgSet(x), data),
            RooFit.Normalization(nData, RooAbsReal.NumEvent),
            RooFit.Name(curve_name),
        )

        can = TCanvas("c_x", "c_x", 820, 920)
        frame.Draw()

        can.SaveAs(f"plots/fit_2d/control_mumu_{outprefix}_{pdf_family}.pdf")

        res = ChiSquareResult(chi2=chi2, pvalue=pval, ndf=ndf, chi2_ndf=chi2_ndf)
        return res
