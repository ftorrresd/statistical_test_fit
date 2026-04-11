from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import multiprocessing as mp
import os
import time
import traceback
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple


@dataclass(frozen=True)
class ParallelJob:
    key: str
    label: str
    payload: Any


def resolve_worker_count(requested_workers: Optional[int], n_jobs: int) -> int:
    if n_jobs <= 0:
        return 0

    cpu_count = os.cpu_count() or 1
    if requested_workers is None:
        return max(1, min(cpu_count, n_jobs))

    return max(1, min(int(requested_workers), n_jobs))


def _format_elapsed(seconds: float) -> str:
    return time.strftime("%H:%M:%S", time.gmtime(seconds))


def _render_progress(
    stage_name: str,
    completed: int,
    total: int,
    *,
    last_label: str,
    last_status: str,
    started_at: float,
) -> str:
    bar_width = 28
    fraction = 1.0 if total == 0 else completed / total
    filled = int(round(bar_width * fraction))
    bar = "#" * filled + "-" * (bar_width - filled)
    percent = 100.0 * fraction
    elapsed = _format_elapsed(time.monotonic() - started_at)
    return (
        f"{stage_name} [{bar}] {completed}/{total} {percent:5.1f}% "
        f"last={last_label} {last_status} elapsed={elapsed}"
    )


def print_job_list(stage_name: str, jobs: Iterable[ParallelJob]) -> List[ParallelJob]:
    job_list = list(jobs)
    print(f"\n{stage_name} jobs ({len(job_list)}):")
    for idx, job in enumerate(job_list, start=1):
        print(f"[{idx}/{len(job_list)}] {job.label}")
    return job_list


def run_parallel_jobs(
    stage_name: str,
    jobs: Iterable[ParallelJob],
    worker_fn: Callable[[Any], Any],
    *,
    workers: Optional[int] = None,
) -> Dict[str, Any]:
    job_list = print_job_list(stage_name, jobs)
    if len(job_list) == 0:
        return {}

    max_workers = resolve_worker_count(workers, len(job_list))
    print(f"Using {max_workers} worker(s) for {stage_name}.")

    started_at = time.monotonic()
    results: Dict[str, Any] = {}
    failures: List[Tuple[ParallelJob, str]] = []
    context = mp.get_context("spawn")

    with ProcessPoolExecutor(max_workers=max_workers, mp_context=context) as executor:
        future_to_job = {
            executor.submit(worker_fn, job.payload): job for job in job_list
        }
        for completed, future in enumerate(as_completed(future_to_job), start=1):
            job = future_to_job[future]
            status = "OK"
            try:
                results[job.key] = future.result()
            except Exception:
                status = "FAILED"
                failures.append((job, traceback.format_exc()))

            print(
                "\r"
                + _render_progress(
                    stage_name,
                    completed,
                    len(job_list),
                    last_label=job.label,
                    last_status=status,
                    started_at=started_at,
                ),
                end="",
                flush=True,
            )

    print()

    if failures:
        print(f"\n{stage_name} failures:")
        for job, tb in failures:
            print(f"--- {job.label} ---")
            print(tb)
        raise RuntimeError(f"{stage_name} failed for {len(failures)} job(s).")

    return results


def run_serial_jobs(
    stage_name: str,
    jobs: Iterable[ParallelJob],
    worker_fn: Callable[[Any], Any],
) -> Dict[str, Any]:
    job_list = print_job_list(stage_name, jobs)
    if len(job_list) == 0:
        return {}

    started_at = time.monotonic()
    results: Dict[str, Any] = {}
    for completed, job in enumerate(job_list, start=1):
        results[job.key] = worker_fn(job.payload)
        print(
            "\r"
            + _render_progress(
                stage_name,
                completed,
                len(job_list),
                last_label=job.label,
                last_status="OK",
                started_at=started_at,
            ),
            end="",
            flush=True,
        )

    print()
    return results
