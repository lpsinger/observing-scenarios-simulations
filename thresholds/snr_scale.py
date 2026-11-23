"""Convert the optimal SNRs of the injection release to those of the simulations.

For one detector and one source, the two SNRs differ by the ratio of the
noise-weighted norms of the signal,

    factor = sigma(simulation template, noise curve of the simulations)
             / sigma(release template, release PSD of the injection's month),

which depends on the masses through the frequency content of the signal.
write_table evaluates it on a grid of detector-frame total mass and mass ratio,
without spin and face on, for each detector and month; read_table interpolates
it. The o2-uberbank templates above 4 Msun need SEOBNRv4ROM_v3.0.hdf5 in
LAL_DATA_PATH.
"""

import csv
import gzip
import logging

import numpy as np
from config import (
    F_HIGH,
    O4_MONTHS,
    RELEASE_F_LOW,
    RELEASE_WAVEFORM,
    RUN_DETECTORS,
    SIMULATION_F_LOW,
    SIMULATION_WAVEFORM,
)
from scipy.interpolate import RegularGridInterpolator
from tqdm import tqdm

# Detector-frame total mass [Msun] from the lightest BNS of the release
# (1 + 1 Msun) to 100 + 100 Msun at z ~ 4, and mass ratio down to the lightest
# NSBH (1 + 60 Msun). The o2-uberbank templates switch from TaylorF2 to
# SEOBNRv4_ROM at 4 Msun, where the SNR jumps by about 4%; a grid point on
# each side keeps the jump sharp.
TOTAL_MASS = np.union1d(np.geomspace(1.9, 1000.0, 49), [3.999, 4.001])
MASS_RATIO = np.array([0.015, 0.03, 0.06, 0.12, 0.25, 0.5, 1.0])
MIN_COMPONENT_MASS = 0.9

# Frequency grid of the PSDs written by scripts/pack-psds.py.
PACKED_GRID = np.arange(10.0, 4096.0, 1.0)


def read_release_psd(path):
    """(frequency, PSD) of a gzipped release file, comma- or space-separated."""
    with gzip.open(path, "rt") as f:
        rows = [line.replace(",", " ").split() for line in f if line[0].isdigit()]
    frequency, psd = np.array(rows, float).T
    keep = frequency > 0
    return frequency[keep], psd[keep]


def read_noise_curve(path):
    """(frequency, PSD) of an ASD file, resampled as scripts/pack-psds.py does."""
    frequency, asd = np.loadtxt(path, unpack=True)
    keep = frequency > 0
    log_psd = np.interp(
        np.log(PACKED_GRID), np.log(frequency[keep]), np.log(np.square(asd[keep]))
    )
    return PACKED_GRID, np.exp(log_psd)


def release_psd_files(directory, run):
    """{(ifo, month): path} of the release PSDs of one run; month is 'O3' in O3."""
    if run.startswith("O3"):
        months, pattern = ["O3"], "psd-O3-{d}.*gz"
    else:
        months, pattern = O4_MONTHS[run], "psd-O4-{m}_v1-{d}.*gz"
    files = {}
    for ifo in RUN_DETECTORS[run]:
        for month in months:
            found = sorted(directory.rglob(pattern.format(m=month, d=ifo[0])))
            if found:
                files[ifo, month] = found[0]
    return files


def template(waveform, total_mass, mass_ratio, f_low):
    """(frequency, |h(f)|^2, df) of a non-spinning face-on template, as in bayestar."""
    from ligo.skymap.bayestar.filter import abscissa, sngl_inspiral_psd

    mass1 = total_mass / (1 + mass_ratio)
    h = sngl_inspiral_psd(waveform, mass1, mass_ratio * mass1, f_min=f_low)
    frequency, power = abscissa(h), h.data.data
    keep = (frequency >= f_low) & (frequency < F_HIGH) & (power > 0)
    return frequency[keep], power[keep], h.deltaF


def sigma(signal, psd):
    """Noise-weighted norm sqrt(4 sum |h|^2 / S df), PSD interpolated in log-log."""
    frequency, power, df = signal
    noise = np.exp(np.interp(np.log(frequency), np.log(psd[0]), np.log(psd[1])))
    return np.sqrt(4 * df * np.sum(power / noise))


def norms(waveform, f_low, psds):
    """{key: sigma on the TOTAL_MASS x MASS_RATIO grid} for each PSD of ``psds``.

    A cell whose lighter component would be below MIN_COMPONENT_MASS takes the
    value of the smallest valid mass ratio at the same total mass.
    """
    grids = {key: np.full((len(TOTAL_MASS), len(MASS_RATIO)), np.nan) for key in psds}
    cells = [
        (i, j, total_mass, mass_ratio)
        for i, total_mass in enumerate(TOTAL_MASS)
        for j, mass_ratio in enumerate(MASS_RATIO)
        if mass_ratio * total_mass / (1 + mass_ratio) >= MIN_COMPONENT_MASS
    ]
    for i, j, total_mass, mass_ratio in tqdm(cells, desc=waveform, unit="template"):
        signal = template(waveform, total_mass, mass_ratio, f_low)
        for key, psd in psds.items():
            grids[key][i, j] = sigma(signal, psd)
    for grid in grids.values():
        for row in grid:
            row[: np.argmax(np.isfinite(row))] = row[np.isfinite(row)][0]
    return grids


def write_table(runs, noise_curves, release_psds, path, comments):
    """Write the factors of each run, detector and month to the CSV ``path``."""
    logging.getLogger("BAYESTAR").setLevel(logging.ERROR)
    targets = {ifo: read_noise_curve(p) for ifo, p in noise_curves.items()}
    release = {}
    for run in runs:
        files = release_psd_files(release_psds, run)
        if not files:
            raise FileNotFoundError(f"no release PSD for {run} in {release_psds}")
        release.update(
            {
                (run, ifo, month): read_release_psd(p)
                for (ifo, month), p in files.items()
            }
        )

    simulation = norms(SIMULATION_WAVEFORM, SIMULATION_F_LOW, targets)
    rows = []
    for (run, ifo, month), norm in norms(
        RELEASE_WAVEFORM, RELEASE_F_LOW, release
    ).items():
        # A signal with no power above the release cutoff has no release SNR,
        # so its factor does not matter; it is set to zero.
        factor = np.divide(
            simulation[ifo], norm, out=np.zeros_like(norm), where=norm > 0
        )
        for i, total_mass in enumerate(TOTAL_MASS):
            for j, mass_ratio in enumerate(MASS_RATIO):
                rows.append(
                    {
                        "run": run,
                        "ifo": ifo,
                        "month": month,
                        "total_mass": round(total_mass, 4),
                        "mass_ratio": mass_ratio,
                        "factor": round(factor[i, j], 5),
                    }
                )

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        f.writelines(f"# {line}\n" for line in comments)
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_table(path):
    """{(run, ifo, month): factor(total_mass, mass_ratio)} from write_table.

    The factor is linear in log total mass and log mass ratio, and keeps the
    values of the grid edges outside it.
    """
    with open(path) as f:
        rows = list(csv.DictReader(line for line in f if not line.startswith("#")))
    total_mass = np.unique([float(r["total_mass"]) for r in rows])
    mass_ratio = np.unique([float(r["mass_ratio"]) for r in rows])
    grids = {}
    for r in rows:
        grid = grids.setdefault(
            (r["run"], r["ifo"], r["month"]),
            np.full((len(total_mass), len(mass_ratio)), np.nan),
        )
        i = np.searchsorted(total_mass, float(r["total_mass"]))
        j = np.searchsorted(mass_ratio, float(r["mass_ratio"]))
        grid[i, j] = float(r["factor"])

    axes = np.log(total_mass), np.log(mass_ratio)

    def interpolator(grid):
        function = RegularGridInterpolator(axes, grid)

        def factor(m, q):
            m = np.clip(np.log(m), axes[0][0], axes[0][-1])
            q = np.clip(np.log(q), axes[1][0], axes[1][-1])
            return function(np.column_stack([m, q]))

        return factor

    return {key: interpolator(grid) for key, grid in grids.items()}


def nearest(month, available):
    """The month of ``available`` closest to ``month``, both 'YYYY_MM'."""

    def index(m):
        year, number = map(int, m.split("_"))
        return 12 * year + number

    return min(available, key=lambda m: abs(index(m) - index(month)))
