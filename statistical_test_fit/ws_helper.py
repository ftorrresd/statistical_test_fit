from typing import Dict
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


def get_pdf_parameters(
    pdf: RooAbsPdf,
    observables: RooArgSet = RooArgSet(),
    sufix=None,
    skip_constants=False,
) -> Dict[str, float]:
    """
    Build a dictionary mapping parameter names to their current values.

    Parameters
    ----------
    pdf : RooAbsPdf
        The RooFit PDF whose parameters you want to extract.
    observables : RooArgSet
        The set of observables (used to exclude them from parameters).

    Returns
    -------
    Dict[str, float]
        A dictionary of {parameter_name: value}.
    """
    params = pdf.getParameters(observables)
    result = {}
    for p in params:
        # optional: skip observables or fixed constants
        if p.isConstant() and skip_constants:
            continue

        name = p.GetName()
        if sufix is not None:
            name = name.replace(f"_{sufix}", "")
        result[name] = p.getVal()
    return result


def set_pdf_parameters(
    pdf: RooAbsPdf,
    params_dict: Dict[str, float],
    observables: RooArgSet,
    make_constant: bool = False,
    sufix=None,
) -> None:
    """
    Set parameter values of a PDF from a dictionary, optionally making them constant.

    Parameters
    ----------
    pdf : RooAbsPdf
        The RooFit PDF whose parameters will be updated.
    params_dict : Dict[str, float]
        Dictionary of {parameter_name: new_value}.
    observables : RooArgSet
        The set of observables (to exclude them from parameter list).
    make_constant : bool, optional
        If True, parameters are set constant after being updated.
    """
    params = pdf.getParameters(observables)
    for p in params:
        name = p.GetName()
        if sufix is not None:
            name = name.replace(f"_{sufix}", "")

        p.setVal(params_dict[name])
        if make_constant:
            p.setConstant(True)
