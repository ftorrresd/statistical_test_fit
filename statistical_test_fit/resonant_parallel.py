from dataclasses import dataclass
from typing import Dict, Optional

from .parallel_utils import ParallelJob


@dataclass(frozen=True)
class ResonantJobPayload:
    kind: str
    nbins: int
    plot_dir: str
    load_from_cache: bool = False
    boson_parameters: Optional[Dict[str, float]] = None
    control_region_name: Optional[str] = None


def build_resonant_stage1_jobs(*, nbins: int, plot_dir: str, load_from_cache: bool):
    payloads = [
        ParallelJob(
            key="higgs_workspace",
            label="Higgs resonant workspace",
            payload=ResonantJobPayload(
                kind="higgs_workspace",
                nbins=nbins,
                plot_dir=plot_dir,
            ),
        ),
        ParallelJob(
            key="z_workspace",
            label="Z resonant workspace",
            payload=ResonantJobPayload(
                kind="z_workspace",
                nbins=nbins,
                plot_dir=plot_dir,
            ),
        ),
        ParallelJob(
            key="z_boson_parameters",
            label="Z resonant boson fit",
            payload=ResonantJobPayload(
                kind="z_boson_parameters",
                nbins=nbins,
                plot_dir=plot_dir,
                load_from_cache=load_from_cache,
            ),
        ),
    ]
    return payloads


def build_resonant_cr_jobs(
    *,
    boson_parameters: Dict[str, float],
    nbins: int,
    plot_dir: str,
    load_from_cache: bool,
):
    jobs = []
    for control_region_name in ("CR1", "CR2", "CR3", "CR4"):
        jobs.append(
            ParallelJob(
                key=control_region_name,
                label=f"Normalization fit {control_region_name}",
                payload=ResonantJobPayload(
                    kind="control_region",
                    nbins=nbins,
                    plot_dir=plot_dir,
                    load_from_cache=load_from_cache,
                    boson_parameters=boson_parameters,
                    control_region_name=control_region_name,
                ),
            )
        )
    return jobs


def run_resonant_job(payload: ResonantJobPayload):
    from .root_runtime import configure_root

    configure_root()

    from .resonant_bkg_modeling import (
        ControlRegion,
        build_resonant_background_Higgs_ws,
        build_resonant_background_Z_ws,
        get_normalization_from_CR,
        resonant_background_modeling_Z,
    )

    if payload.kind == "higgs_workspace":
        build_resonant_background_Higgs_ws(
            plot_dir=payload.plot_dir, nbins=payload.nbins
        )
        return {"workspace": "resonant_background_fit_HiggsDalitz.root"}

    if payload.kind == "z_workspace":
        build_resonant_background_Z_ws(plot_dir=payload.plot_dir, nbins=payload.nbins)
        return {"workspace": "resonant_background_fit_ZGamma.root"}

    if payload.kind == "z_boson_parameters":
        return resonant_background_modeling_Z(
            load_from_cache=payload.load_from_cache,
            plot_dir=payload.plot_dir,
            nbins=payload.nbins,
        )

    if payload.kind == "control_region":
        assert payload.boson_parameters is not None
        assert payload.control_region_name is not None
        control_region = ControlRegion[payload.control_region_name]
        return get_normalization_from_CR(
            payload.boson_parameters,
            control_region,
            load_from_cache=payload.load_from_cache,
            plot_dir=payload.plot_dir,
            nbins=payload.nbins,
        )

    raise ValueError(f"Unsupported resonant job kind: {payload.kind}")
