UPSILON_MASS_LOWER = 8.0
UPSILON_MASS_UPPER = 12.0

BOSON_MASS_LOWER = 57.0
BOSON_MASS_UPPER = 200.0

LEFT_SIDEBAND_LOWER = BOSON_MASS_LOWER
LEFT_SIDEBAND_UPPER = 75.0
MIDDLE_SIDEBAND_LOWER = 105.0
MIDDLE_SIDEBAND_UPPER = 115.0
RIGHT_SIDEBAND_LOWER = 135.0
RIGHT_SIDEBAND_UPPER = BOSON_MASS_UPPER

RESONANT_CR_UPSILON_MASS_LOWER = 4.0
RESONANT_CR_UPSILON_MASS_UPPER = 35.0

UPSILON_MASS_SEEDS = {
    "1S": 9.46,
    "2S": 10.02,
    "3S": 10.35,
}


def get_signal_boson_plot_range(process: str) -> tuple[float, float]:
    if process == "Z":
        return (75.0, 110.0)
    if process == "H":
        return (110.0, 140.0)
    raise ValueError(f"Unsupported process {process!r}")


def get_signal_upsilon_plot_range(state: str) -> tuple[float, float]:
    center = UPSILON_MASS_SEEDS[state]
    half_width = 0.5
    return (
        max(UPSILON_MASS_LOWER, center - half_width),
        min(UPSILON_MASS_UPPER, center + half_width),
    )
