from pprint import pprint
import os
import time
from functools import wraps
from argparse import Namespace

from ROOT import (  # type: ignore
    RooArgList,  # type: ignore
    RooArgSet,  # type: ignore
    RooCategory,  # type: ignore
    RooDataSet,  # type: ignore
    RooEffProd,  # type: ignore
    RooFit,  # type: ignore
    RooFormulaVar,  # type: ignore
    RooMultiPdf,  # type: ignore
    RooRealVar,  # type: ignore
    RooWorkspace,  # type: ignore
    TFile,  # type: ignore
)

from .bkg_model import (
    BkgModel,
    compute_winner_and_start_indexes,
)
from .bkg_pdf_families import BkgPdfFamily
from .dimuon_non_correlated import dimuon_non_correlated
from .mass_ranges import (
    BOSON_MASS_LOWER,
    BOSON_MASS_UPPER,
    LEFT_SIDEBAND_LOWER,
    LEFT_SIDEBAND_UPPER,
    MIDDLE_SIDEBAND_LOWER,
    MIDDLE_SIDEBAND_UPPER,
    RIGHT_SIDEBAND_LOWER,
    RIGHT_SIDEBAND_UPPER,
    UPSILON_MASS_LOWER,
    UPSILON_MASS_UPPER,
)
from .make_plots import DataType, ProjDim, make_plots_2d
from .non_resonant_parallel import (
    build_non_resonant_candidate_jobs,
    build_non_resonant_candidate_specs,
    fit_non_resonant_candidate_job,
    reconstruct_non_resonant_candidate,
)
from .parallel_utils import ParallelJob, run_parallel_jobs, run_serial_jobs


def execution_time(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        formatted_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
        print(f"Execution time of {func.__name__}: {formatted_time}")
        return result

    return wrapper


@execution_time
def run_fit_2d_data(args: Namespace):
    w = RooWorkspace("ws")

    # Limits
    upsilon_mass_lower = UPSILON_MASS_LOWER
    upsilon_mass_upper = UPSILON_MASS_UPPER

    left_lower = LEFT_SIDEBAND_LOWER
    left_upper = LEFT_SIDEBAND_UPPER

    middle_lower = MIDDLE_SIDEBAND_LOWER
    middle_upper = MIDDLE_SIDEBAND_UPPER

    right_lower = RIGHT_SIDEBAND_LOWER
    right_upper = RIGHT_SIDEBAND_UPPER

    outprefix = "bkg_only"

    LOAD_FROM_CACHE = args.use_cache
    if not LOAD_FROM_CACHE:
        os.system("rm -rf *.json")

    def _run_dimuon_fit(payload):
        return dimuon_non_correlated(
            payload["upsilon_mass_lower"],
            payload["upsilon_mass_upper"],
            load_from_cache=payload["load_from_cache"],
        )

    dimuon_jobs = [
        ParallelJob(
            key="dimuon_fit",
            label="Dimuon non-correlated fit",
            payload={
                "upsilon_mass_lower": upsilon_mass_lower,
                "upsilon_mass_upper": upsilon_mass_upper,
                "load_from_cache": LOAD_FROM_CACHE,
            },
        )
    ]
    upsilon_params = run_serial_jobs(
        "Non-resonant setup",
        dimuon_jobs,
        _run_dimuon_fit,
    )["dimuon_fit"]
    print("\n\nUpsilon Parameters:")
    pprint(upsilon_params)
    print("\n\n")

    # Observable
    upsilon_mass = RooRealVar(
        "upsilon_mass", "upsilon_mass", upsilon_mass_lower, upsilon_mass_upper
    )
    upsilon_mass.SetTitle("m_{#mu#mu}")  # LaTeX-style title
    upsilon_mass.setUnit("GeV")  # physical unit

    boson_mass = RooRealVar(
        "boson_mass", "boson_mass", BOSON_MASS_LOWER, BOSON_MASS_UPPER
    )
    boson_mass.SetTitle("m_{#mu#mu#gamma}")  # LaTeX-style title
    boson_mass.setUnit("GeV")  # physical unit

    f = TFile.Open("inputs/mass_Run2.root")
    data_full = RooDataSet(
        "data_obs",
        "data_obs",
        RooArgSet(
            upsilon_mass,
            boson_mass,
        ),
        RooFit.Import(f.Events),
        RooFit.Cut(
            f"(upsilon_mass < {upsilon_mass_upper} && upsilon_mass >= {upsilon_mass_lower} && boson_mass>={left_lower} && boson_mass<{right_upper})"
        ),
    )
    getattr(w, "import")(data_full)

    # Named ranges for sidebands
    boson_mass.setRange("LEFT", left_lower, left_upper)
    boson_mass.setRange("MIDDLE", middle_lower, middle_upper)
    boson_mass.setRange("RIGHT", right_lower, right_upper)
    boson_mass.setRange("FULL", left_lower, right_upper)

    # Sideband-only dataset (useful for plotting)
    cut_expr = f"((boson_mass<{left_upper}) || (boson_mass>{right_lower}) || ((boson_mass>{middle_lower}) && (boson_mass<{middle_upper})))"
    # data_sb = data_full.reduce(RooFit.Cut(cut_expr))
    data_sb = data_full.reduce(RooFit.Cut(cut_expr))
    print(f"Sideband entries: {data_sb.numEntries()} (out of {data_full.numEntries()})")

    test_bkg_pdfs: dict[BkgPdfFamily, list[BkgModel]] = {}
    for family in BkgPdfFamily:
        test_bkg_pdfs[family] = []

    candidate_specs = build_non_resonant_candidate_specs(
        upsilon_params=upsilon_params,
        nbins=args.nbins,
        input_file="inputs/mass_Run2.root",
        upsilon_mass_lower=upsilon_mass_lower,
        upsilon_mass_upper=upsilon_mass_upper,
        boson_mass_lower=BOSON_MASS_LOWER,
        boson_mass_upper=BOSON_MASS_UPPER,
        left_lower=left_lower,
        left_upper=left_upper,
        middle_lower=middle_lower,
        middle_upper=middle_upper,
        right_lower=right_lower,
        right_upper=right_upper,
    )
    candidate_results = run_parallel_jobs(
        "Non-resonant candidate fits",
        build_non_resonant_candidate_jobs(candidate_specs),
        fit_non_resonant_candidate_job,
        workers=getattr(args, "workers", None),
    )

    for spec in candidate_specs:
        result = candidate_results[spec.key]
        test_bkg_pdfs[spec.pdf_family].append(
            reconstruct_non_resonant_candidate(result, upsilon_mass, boson_mass)
        )

    winners = {}
    for family in BkgPdfFamily:
        if len(test_bkg_pdfs[family]) > 0:
            print()
            print("##############################################")
            print("##############################################")
            print(f"################## {str(family).upper()} ########################")
            print("##############################################")
            print("##############################################")
            print()
            for test_bkg_pdf in test_bkg_pdfs[family]:
                print(f"Completed worker fit for PDF: {test_bkg_pdf.model.GetName()}")

            print("\n\n=== Test Background-only fit (sidebands) ===")
            for i, test_bkg_pdf in enumerate(test_bkg_pdfs[family]):
                assert test_bkg_pdf.is_complete()
                if i != 0:
                    print("")

                print(test_bkg_pdf)

            # compute winner function
            start, winner = compute_winner_and_start_indexes(
                test_bkg_pdfs[family],
                strict_mode=args.strict_mode,
            )
            winners[family] = winner

            # Plots
            plot_file_name = make_plots_2d(
                ProjDim.Y,
                upsilon_mass,
                boson_mass,
                # data_full,
                None,
                data_sb,
                None,
                test_bkg_pdfs[family],
                family,
                left_lower,
                left_upper,
                middle_lower,
                middle_upper,
                right_lower,
                right_upper,
                f"{outprefix}_{family}",
                nbins=args.nbins,
                data_type=DataType.REAL,
                start=start,
                winner=winner,
                # components=[("Z Resonant BKG", test_bkg_pdfs[family][0].model)],
            )
            print(f"\nPlot saved to: {plot_file_name}")

            plot_file_name = make_plots_2d(
                ProjDim.X,
                upsilon_mass,
                boson_mass,
                # data_full,
                None,
                data_sb,
                None,
                test_bkg_pdfs[family],
                family,
                left_lower,
                left_upper,
                middle_lower,
                middle_upper,
                right_lower,
                right_upper,
                f"{outprefix}_{family}",
                nbins=args.nbins,
                data_type=DataType.REAL,
                start=start,
                winner=winner,
            )
            print(f"\nPlot saved to: {plot_file_name}")

            # Optional: estimate expected background count in the blinded window (from fitted pdf)
            # Compute integral of background over blind window, normalized over full x.
            boson_mass.setRange("BLIND_LEFT", left_upper, middle_lower)
            boson_mass.setRange("BLIND_RIGHT", middle_upper, right_lower)

    pdf_index = RooCategory("pdfIndex", "pdfIndex")
    getattr(w, "import")(pdf_index)
    pdf_list = RooArgList()
    for family in BkgPdfFamily:
        if len(test_bkg_pdfs[family]) > 0:
            w.cat("pdfIndex").defineType(family.value)
            print(f"--> Adding winner for {family}")
            getattr(w, "import")(test_bkg_pdfs[family][winners[family]].model)
            # pdf_list.add(w.pdf(family.value))
            pdf_list.add(w.pdf(test_bkg_pdfs[family][winners[family]].model.GetName()))

    multi_pdf = RooMultiPdf(
        "multiPdf",
        "multiPdf",
        w.cat("pdfIndex"),
        pdf_list,
    )
    getattr(w, "import")(multi_pdf)

    w.Print("v")
    print("pdfIndex: ", w.cat("pdfIndex").states())
