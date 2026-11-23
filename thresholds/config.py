"""Inputs and choices of the measurement of rho_eq.csv.

The comment above a value names its source when it has one: the GWTC-5.0
injection release and catalogue, or the simulations of this repository. The
chirp-mass bins, MIN_FOUND, NOISE_DRAWS and BOOTSTRAP are choices of the
measurement.
"""

from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent

# GWTC-5.0 cumulative sensitivity estimates, doi:10.5281/zenodo.19500052
# (Essick et al. 2025, arXiv:2508.10638): the injections, and the PSDs the
# release used for their optimal SNRs. Downloaded into data/, see README.md.
DATA = HERE / "data"
INJECTIONS = (
    DATA
    / "mixture-semi_o1_o2-real_o3_o4a_o4b-cartesian_spins_20260410130052UTC-clipped.hdf"
)
INJECTIONS_DOI = "10.5281/zenodo.19500052"
RELEASE_PSDS = DATA / "psds-o1234ab"
SEGMENTS = DATA / "segments"
SNR_SCALE = DATA / "snr_scale.csv"

# Noise curves of the simulations of this repository, downloaded by the Makefile.
NOISE_CURVES = {ifo: REPO / f"o4b_{ifo.lower()}_ref.txt" for ifo in ("H1", "L1", "V1")}

# How each side computes an optimal SNR. The release uses IMRPhenomXPHM from
# 16 Hz with the PSD of the injection's month (Essick et al. 2025, Sec. V and
# App. A). bayestar-inject and bayestar-realize-coincs use the o2-uberbank
# templates (TaylorF2 below 4 Msun total mass, SEOBNRv4_ROM above) from 25 Hz.
RELEASE_WAVEFORM = "IMRPhenomXPHM"
RELEASE_F_LOW = 16.0
SIMULATION_WAVEFORM = "o2-uberbank"
SIMULATION_F_LOW = 25.0
F_HIGH = 2048.0

# Runs of the injection release: GPS boundaries [start, end) and detectors.
# The position of a run in RUNS enters the seed of its noise draws.
RUN_GPS = {
    "O3a": (1238166018, 1253977218),
    "O3b": (1256655618, 1269363618),
    "O4a": (1368975618, 1389456018),
    "O4b": (1396969218, 1422118818),
}
RUN_DETECTORS = {
    "O3a": ("H1", "L1", "V1"),
    "O3b": ("H1", "L1", "V1"),
    "O4a": ("H1", "L1"),
    "O4b": ("H1", "L1", "V1"),
}
RUNS = tuple(RUN_GPS)

# Calendar months of the O4 runs; the release has one PSD per month and
# detector, and one per detector for the whole of O3.
O4_MONTHS = {
    "O4a": [f"2023_{m:02d}" for m in range(5, 13)] + ["2024_01"],
    "O4b": [f"2024_{m:02d}" for m in range(4, 13)] + ["2025_01"],
}

# A detection is a FAR below 1/yr in at least one search, as in GWTC-5.0.
FAR_THRESHOLD = 1.0

# Source classes from the source-frame component masses, with the NS/BH
# boundary of the injections. Injections with a component above BH_MAX_MASS,
# the upper end of the FullPop mass range, are left out.
NS_MAX_MASS = 2.5
BH_MAX_MASS = 100.0
CLASSES = ("BNS", "NSBH", "BBH")

# Detector-frame chirp-mass bins [Msun], about a factor 1.7 wide, halved below
# 2.2 Msun where rho_eq changes faster with mass. They span the lightest BNS of
# the injections (1 + 1 Msun) to 100 + 100 Msun BBH at z ~ 2.4.
CHIRP_MASS_EDGES = (
    0.85,
    1.12,
    1.4,
    1.75,
    2.2,
    4.0,
    7.0,
    12.0,
    20.0,
    35.0,
    60.0,
    120.0,
    300.0,
)

# A bin is measured only with at least this many found injections.
MIN_FOUND = 30

# Observed SNR as bayestar-realize-coincs simulates it (Singer & Price 2016):
# |rho_opt + n| in each observing detector, n complex Gaussian with unit
# variance per quadrature, plus the mean excess of its peak search. A detector
# counts above SINGLE_DETECTOR_MIN (--snr-threshold 1 in the Makefile).
NOISE_DRAWS = 16  # 2048 draws move rho_eq by less than half its error
PEAK_OFFSET = 0.1  # measured on bayestar-realize-coincs, 0.07-0.14
SINGLE_DETECTOR_MIN = 1.0
BOOTSTRAP = 20000  # resamples for the error on rho_eq
