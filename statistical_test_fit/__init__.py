__all__ = [
    "run_fit_1d",
    "run_fit_2d",
    "run_fit_2d_data",
    "run_signal_modeling",
    "run_resonant_background",
]


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

    if name == "run_signal_modeling":
        from .signal_modeling import run_signal_modeling

        return run_signal_modeling

    if name == "run_resonant_background":
        from .resonant_bkg_modeling import run_resonant_background

        return run_resonant_background

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
