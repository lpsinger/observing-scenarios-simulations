"""Measure rho_eq.csv, the network SNR threshold equivalent to FAR < 1/yr.

For the injections of one run, source class and detector-frame chirp-mass bin:

1. an injection is found when its smallest FAR over the searches is below
   1/yr;
2. its optimal SNR in each detector is converted from the release to the
   simulations (snr_scale.py), and its network SNR is simulated as
   bayestar-realize-coincs does, with the detectors observing at its GPS time
   (GWOSC segments), over NOISE_DRAWS noise realizations;
3. rho_eq is the SNR above which the fraction of simulated values equals the
   found fraction, so that as many found injections lie below it as missed
   injections above it. Its error is the spread over bootstrap resamples of
   the injections.

    python thresholds/measure.py

README.md lists the inputs to download first. Most of the run time is spent
in the bootstrap.
"""

import argparse
import csv
import shlex
import sys
from itertools import pairwise
from pathlib import Path

import h5py
import numpy as np
import snr_scale
from config import (
    BH_MAX_MASS,
    BOOTSTRAP,
    CHIRP_MASS_EDGES,
    CLASSES,
    FAR_THRESHOLD,
    HERE,
    INJECTIONS,
    INJECTIONS_DOI,
    MIN_FOUND,
    NOISE_CURVES,
    NOISE_DRAWS,
    NS_MAX_MASS,
    PEAK_OFFSET,
    RELEASE_PSDS,
    RUN_DETECTORS,
    RUN_GPS,
    RUNS,
    SEGMENTS,
    SIMULATION_F_LOW,
    SIMULATION_WAVEFORM,
    SINGLE_DETECTOR_MIN,
    SNR_SCALE,
)
from tqdm import tqdm

GPS_EPOCH_UTC = np.datetime64("1980-01-06T00:00:00") - np.timedelta64(18, "s")


def load_injections(path, run, far_threshold):
    """Injections of ``run`` with no source-frame mass above BH_MAX_MASS.

    Returns a dict of arrays with one entry per injection, and the optimal SNR
    of the release in each detector, zero where the detector did not observe.
    """
    prefix = "o3_" if run.startswith("O3") else f"{run.lower()}_"
    with h5py.File(path, "r") as f:
        events = f["events"]
        gps = events["time_geocenter"][:]
        start, end = RUN_GPS[run]
        rows = np.flatnonzero((gps >= start) & (gps < end))

        def column(name):
            return events[name][rows]

        far = [
            column(n)
            for n in events.dtype.names
            if n.startswith(prefix) and n.endswith("_far")
        ]
        if not far:
            raise KeyError(f"no search results for {run} in {path}")
        z = column("redshift")
        m1, m2 = column("mass1_source"), column("mass2_source")
        optimal_snr = {
            ifo: np.nan_to_num(column(f"estimated_optimal_snr_{ifo[0]}"))
            for ifo in RUN_DETECTORS[run]
        }

    heavy, light = np.maximum(m1, m2), np.minimum(m1, m2)
    injections = {
        "gps": gps[rows],
        "found": np.min(far, axis=0) < far_threshold,
        "source_class": np.where(
            heavy < NS_MAX_MASS, "BNS", np.where(light < NS_MAX_MASS, "NSBH", "BBH")
        ),
        "chirp_mass_det": (1 + z) * ((m1 * m2) ** 0.6 / (m1 + m2) ** 0.2),
        "total_mass_det": (1 + z) * (heavy + light),
        "mass_ratio": light / heavy,
    }
    keep = heavy <= BH_MAX_MASS
    injections = {key: value[keep] for key, value in injections.items()}
    optimal_snr = {ifo: snr[keep] for ifo, snr in optimal_snr.items()}
    return injections, optimal_snr


def convert_snr(injections, optimal_snr, run, table):
    """Multiply each optimal SNR by the factor of its detector, month and masses.

    A month without a release PSD for a detector takes the nearest month that
    has one.
    """
    if run.startswith("O3"):
        month = np.full(len(injections["gps"]), "O3")
    else:
        utc = GPS_EPOCH_UTC + injections["gps"].astype("timedelta64[s]")
        month = np.char.replace(utc.astype("datetime64[M]").astype(str), "-", "_")
    for ifo, snr in optimal_snr.items():
        available = sorted(m for r, d, m in table if (r, d) == (run, ifo))
        if not available:
            raise ValueError(f"no SNR conversion for {run} {ifo}")
        factor = np.empty(len(snr))
        for m in np.unique(month):
            use = m if m in available else snr_scale.nearest(m, available)
            selected = month == m
            factor[selected] = table[run, ifo, use](
                injections["total_mass_det"][selected],
                injections["mass_ratio"][selected],
            )
        optimal_snr[ifo] = snr * factor
    return optimal_snr


def observing(gps, run):
    """{ifo: True where ``ifo`` was observing at ``gps``}, from the GWOSC segments."""
    state = {}
    for ifo in RUN_DETECTORS[run]:
        path = SEGMENTS / f"{run}_{ifo}.txt"
        if not path.exists():
            from gwosc.timeline import get_segments

            path.parent.mkdir(parents=True, exist_ok=True)
            segments = np.array(
                get_segments(f"{ifo}_DATA", *RUN_GPS[run]), dtype=np.int64
            )
            np.savetxt(
                path, segments.reshape(-1, 2), fmt="%d", header=f"{ifo}_DATA, {run}"
            )
        table = np.loadtxt(path, dtype=np.int64, ndmin=2)
        start, end = table[np.argsort(table[:, 0])].T
        i = np.searchsorted(start, gps, side="right") - 1
        state[ifo] = (i >= 0) & (gps < end[np.maximum(i, 0)])
    return state


def network_snr(optimal_snr, state, rng):
    """Simulated network SNR, shape (injections, NOISE_DRAWS)."""
    n = len(next(iter(optimal_snr.values())))
    power = np.zeros((n, NOISE_DRAWS))
    for ifo, rho in optimal_snr.items():
        on = np.broadcast_to(state[ifo][:, None], (n, NOISE_DRAWS))
        snr = (
            np.hypot(
                rho[:, None] + rng.standard_normal((n, NOISE_DRAWS)),
                rng.standard_normal((n, NOISE_DRAWS)),
            )
            + PEAK_OFFSET
        )
        power += np.where(on & (snr >= SINGLE_DETECTOR_MIN), snr**2, 0.0)
    return np.sqrt(power)


def equivalent_threshold(snr, found):
    """SNR above which the fraction of simulated values equals the found fraction."""
    return float(np.quantile(snr, 1.0 - found.mean()))


def bootstrap_error(snr, found, rng):
    """Standard deviation of rho_eq over resamples of the injections."""
    estimates = []
    for _ in range(BOOTSTRAP):
        i = rng.integers(0, len(found), len(found))
        estimates.append(equivalent_threshold(snr[i], found[i]))
    return float(np.std(estimates))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--runs", nargs="+", default=["O4a", "O4b"], choices=RUNS)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("-o", "--output", type=Path, default=HERE / "rho_eq.csv")
    args = parser.parse_args()

    curves = ", ".join(f"{ifo} {path.name}" for ifo, path in NOISE_CURVES.items())
    provenance = [
        (
            f"SNRs: optimal SNRs of the simulations ({SIMULATION_WAVEFORM} from "
            f"{SIMULATION_F_LOW:g} Hz) for {curves}"
        ),
        f"command: {shlex.join(['measure.py'] + sys.argv[1:])}",
        f"injections: {INJECTIONS.name}, doi:{INJECTIONS_DOI}",
    ]

    snr_scale.write_table(args.runs, NOISE_CURVES, RELEASE_PSDS, SNR_SCALE, provenance)
    table = snr_scale.read_table(SNR_SCALE)

    results = []
    for run in args.runs:
        injections, optimal_snr = load_injections(INJECTIONS, run, FAR_THRESHOLD)
        optimal_snr = convert_snr(injections, optimal_snr, run, table)
        state = observing(injections["gps"], run)
        inside = np.any(list(state.values()), axis=0)
        injections = {k: v[inside] for k, v in injections.items()}
        optimal_snr = {ifo: snr[inside] for ifo, snr in optimal_snr.items()}
        state = {ifo: on[inside] for ifo, on in state.items()}

        bins = [
            (class_index, cls, bin_index, lo, hi)
            for class_index, cls in enumerate(CLASSES)
            for bin_index, (lo, hi) in enumerate(pairwise(CHIRP_MASS_EDGES))
        ]
        for class_index, cls, bin_index, lo, hi in tqdm(bins, desc=run, unit="bin"):
            selected = (
                (injections["source_class"] == cls)
                & (injections["chirp_mass_det"] >= lo)
                & (injections["chirp_mass_det"] < hi)
            )
            found = injections["found"][selected]
            if found.sum() < MIN_FOUND:
                continue
            # One noise stream per bin, so that a bin does not depend on the
            # ones measured before it.
            rng = np.random.default_rng(
                [args.seed, RUNS.index(run), class_index, bin_index]
            )
            snr = network_snr(
                {ifo: v[selected] for ifo, v in optimal_snr.items()},
                {ifo: v[selected] for ifo, v in state.items()},
                rng,
            )
            rho_eq = equivalent_threshold(snr, found)
            results.append(
                {
                    "run": run,
                    "source_class": cls,
                    "chirp_mass_min": lo,
                    "chirp_mass_max": hi,
                    "n_injections": len(found),
                    "n_found": int(found.sum()),
                    "rho_eq": round(rho_eq, 3),
                    "error": round(bootstrap_error(snr, found, rng), 3),
                }
            )
            tqdm.write(", ".join(f"{k} {v}" for k, v in results[-1].items()))

    with open(args.output, "w", newline="") as f:
        f.writelines(f"# {line}\n" for line in provenance)
        writer = csv.DictWriter(f, fieldnames=list(results[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
