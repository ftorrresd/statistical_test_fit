__all__ = ["run_fit_1d", "run_fit_2d", "run_fit_2d_data"]


def __getattr__(name):
    if name == "run_fit_1d":
        from .fit_1d import run_fit_1d

        return run_fit_1d

    if name == "run_fit_2d":
        from .fit_2d import run_fit_2d

        return run_fit_2d

    if name == "run_fit_2d_data":
        from .fit_2d_data import run_fit_2d_data

        return run_fit_2d_data

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
