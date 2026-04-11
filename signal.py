import importlib.util
import os
import sys
import sysconfig


def _load_stdlib_signal_module():
    stdlib_signal = os.path.join(sysconfig.get_path("stdlib"), "signal.py")
    spec = importlib.util.spec_from_file_location("_stdlib_signal", stdlib_signal)
    if spec is None or spec.loader is None:
        raise ImportError("Could not load the standard-library signal module")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


if __name__ != "__main__":
    _stdlib_signal = _load_stdlib_signal_module()
    for _name in dir(_stdlib_signal):
        if _name.startswith("__") and _name not in {"__all__", "__doc__"}:
            continue
        globals()[_name] = getattr(_stdlib_signal, _name)
else:
    import argparse

    def main():
        from statistical_test_fit.root_runtime import configure_root

        parser = argparse.ArgumentParser(
            description="Signal RooFit modeling for H/Z -> Upsilon(nS) + photon."
        )
        parser.add_argument("--nbins", type=int, default=60)
        parser.add_argument("--workers", type=int, default=None)
        args = parser.parse_args()

        os.system("mkdir -p plots")
        os.system("mkdir -p plots/signal_fit")
        os.system(r'find plots/signal_fit -type f ! -name ".gitkeep" -delete')

        configure_root()

        from statistical_test_fit import run_signal_modeling

        run_signal_modeling(args)

    if __name__ == "__main__":
        main()
