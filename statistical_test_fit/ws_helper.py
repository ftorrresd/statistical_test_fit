from ROOT import (  # type: ignore
    AddressOf,  # type: ignore
    RooAbsPdf,  # type: ignore
    RooArgSet,  # type: ignore
    RooRealVar,  # type: ignore
    RooWorkspace,  # type: ignore
)


def set_constant(w):
    print("--> Setting all var as constants")
    for v in w.allVars().contentsString().split(","):
        if v != "boson_mass" and v != "upsilon_mass" and v != "evt_weight":
            print("{v} to constant...".format(v=v))
            w.var(v).setConstant()
            w.var(v).Print()
    return w


def freeze_pdf_params(pdf: RooAbsPdf, observables: RooArgSet) -> RooArgSet:
    """Return the (now frozen) parameter set.

    If you have a dataset, you can auto-discover observables.

    Example:
    observables = pdf.getObservables(ata)
    frozen = freeze_pdf_params(model, observables)"""

    params = pdf.getParameters(
        observables
    )  # everything the pdf depends on, except 'observables'

    it = params.createIterator()
    obj = it.Next()
    while obj:
        # Freeze only real-valued fit parameters
        if isinstance(obj, RooRealVar) and not obj.isConstant():
            obj.setConstant(True)
        obj = it.Next()
    return params


def unfreeze_pdf_params(pdf: RooAbsPdf, observables: RooArgSet):
    params = pdf.getParameters(observables)
    it = params.createIterator()
    obj = it.Next()
    while obj:
        if isinstance(obj, RooRealVar) and obj.isConstant():
            obj.setConstant(False)
        obj = it.Next()
