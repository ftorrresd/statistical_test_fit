def configure_root():
    import ROOT  # type: ignore

    ROOT.gROOT.SetBatch(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    ROOT.gSystem.Load("libRooFit")
    ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
    return ROOT
