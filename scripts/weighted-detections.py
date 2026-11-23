#!/usr/bin/env python
"""Tabulate the detections of each run that reach the FAR-equivalent threshold."""

import textwrap
from argparse import ArgumentParser, FileType
from pathlib import Path

import numpy as np
import rates
from astropy import units as u
from astropy.table import Table
from ligo.skymap.util import sqlite
from scipy import stats

# Total CBC merger rate density of GWTC-5.0 FullPop [Gpc^-3 yr^-1]: the median
# of its popsummary rate column, and the bounds of the 90% credible interval.
# It is the rate of the whole population rather than of one class, which is
# what lets the fractions of the population cancel in interval() below. A rate
# published per class would not cancel, and would need those fractions, which
# only the population samples of fullpop.h5 give.
FIDUCIAL_RATE = (100.07580795644503, 49.88807658685758, 202.77919891201657)

# Observing time of the runs that have happened [yr]: the calendar lengths of
# the O4a and O4b segments, from the GWOSC detector-status pages. These are
# the runs themselves, not the GPS window of the injection release that
# thresholds/config.py uses. A run absent here, IR1 and O5, has not happened
# and has no observing time of its own, so it is tabulated per year.
SECONDS_PER_YEAR = 86400 * 365.25
DURATIONS = {
    "O4a": (1389456018 - 1368975618) / SECONDS_PER_YEAR,
    "O4b": (1422118818 - 1396796418) / SECONDS_PER_YEAR,
}

# NS/BH boundary of the GWTC-5.0 injections on which the thresholds were
# measured, so that a class means the same thing here as in rho_eq.csv.
NS_MAX_MASS = 2.5

CLASSES = ("BNS", "NSBH", "BBH")

parser = ArgumentParser()
parser.add_argument("runs", nargs="+", type=Path, help="runs/<run>/fullpop directories")
parser.add_argument("-o", "--output", type=FileType("w"), default="-")
parser.add_argument(
    "--population",
    type=Path,
    help="population samples, fullpop.h5, to split the rates by class; the "
    "split cancels in the detection count, so this changes no value",
)
args = parser.parse_args()

# Fraction of the population in each class, as plots-and-tables.ipynb measures
# it. It cancels in interval() below, so it is 1 when the samples are not
# given; the release of the runs does not carry them.
fractions = dict.fromkeys(CLASSES, 1.0)
if args.population is not None:
    samples = Table.read(args.population)
    heavy = np.maximum(samples["mass1"], samples["mass2"])
    light = np.minimum(samples["mass1"], samples["mass2"])
    fractions = {
        name: keep.sum() / len(samples)
        for name, keep in [
            ("BNS", heavy < NS_MAX_MASS),
            ("NSBH", (heavy >= NS_MAX_MASS) & (light < NS_MAX_MASS)),
            ("BBH", light >= NS_MAX_MASS),
        ]
    }

mid, lower, upper = FIDUCIAL_RATE
(width,) = np.diff(stats.norm.interval(0.9))

runs = []
for directory in args.runs:
    name = directory.parent.name
    with sqlite.open(str(directory / "events.sqlite"), "r") as db:
        # Injection rate density of the run, written by bayestar-inject.
        ((comment,),) = db.execute(
            "SELECT comment FROM process WHERE program = ?", ("bayestar-inject",)
        )
        ((ifos,),) = db.execute(
            "SELECT ifos FROM process WHERE program = ?", ("bayestar-realize-coincs",)
        )
    weights = Table.read(directory / "detection-weights.dat", format="ascii.fast_tab")
    passes = np.asarray(weights["passes_rho_eq"], dtype=bool)
    source_class = np.asarray(weights["source_class"])
    runs.append(
        {
            "name": name,
            "network": ifos.replace("1", "").replace(",", ""),
            "years": DURATIONS.get(name, 1.0),
            "simulated": u.Quantity(comment).to_value(u.Gpc**-3 * u.yr**-1),
            "detected": {
                weighted: {
                    name: w[passes & (source_class == name)].sum() for name in CLASSES
                }
                for weighted, w in [
                    (False, np.ones(len(weights))),
                    (True, np.asarray(weights["redshift_weight"])),
                ]
            },
        }
    )


def interval(run, source_class, weighted, years):
    """Detections of one class over ``years``, as median, minus and plus.

    The whole population is drawn from one normalised model at one combined
    rate, so the rate of a class is that combined rate times the fraction of
    the population it holds, and the simulated rate of that class is the
    simulated rate times the same fraction. The fraction cancels between the
    two, as it does in plots-and-tables.ipynb, which is why it need not be
    measured. A rate published per class would not cancel.
    """
    fraction = fractions[source_class]
    detected = run["detected"][weighted][source_class]
    mu = np.log(mid * fraction) + np.log(
        years * detected / (run["simulated"] * fraction)
    )
    sigma = (np.log(upper) - np.log(lower)) / width
    lo, median, hi = rates.poisson_lognormal_rate_quantiles(
        np.array([0.05, 0.5, 0.95]), mu, sigma
    )
    return rates.format_with_errorbars(
        int(np.round(median)), int(np.floor(lo)), int(np.ceil(hi))
    )


# The layout of runs/summary.rst, which this table sits beside.
WIDTHS = [11, 11, 15, 15, 15]
RULE = "+" + "+".join("-" * width for width in WIDTHS) + "+"


def row(cells):
    return "|" + "|".join(f" {c:<{w - 2}} " for c, w in zip(cells, WIDTHS)) + "|"


def section(label, weighted, annual):
    """A full-width label, then one block of three lines per run.

    Only the runs the section has something to say about: a run that has
    happened observed for the time it observed, which is its own block, and a
    run still to come is a forecast, which is the annual one.
    """
    span = sum(WIDTHS) + len(WIDTHS) - 3
    if len(label) > span:
        raise ValueError(f"the label {label!r} is wider than the table")
    if annual:
        listed = [run for run in runs if run["name"] not in DURATIONS]
    else:
        listed = [run for run in runs if run["name"] in DURATIONS]
    if not listed:
        return
    print(f"| {label:<{span}} |", file=args.output)
    print(RULE, file=args.output)
    for run in listed:
        years = 1.0 if annual else run["years"]
        values = [interval(run, c, weighted, years) for c in CLASSES]
        print(
            row([run["name"], run["network"], *(f":math:`{v[0]}" for v in values)]),
            file=args.output,
        )
        # The observing time under the name of the run, where it is not a year.
        print(
            row(
                [
                    "" if annual else f"{years:.2f} yr",
                    "",
                    *(f"^{{+{v[2]}}}" for v in values),
                ]
            ),
            file=args.output,
        )
        print(row(["", "", *(f"_{{-{v[1]}}}`" for v in values)]), file=args.output)
        print(RULE, file=args.output)


# The header of runs/summary.rst, with the three class columns merged.
span = sum(WIDTHS[2:]) + len(WIDTHS[2:]) - 1
print(RULE, file=args.output)
print(
    f"|{' ' * WIDTHS[0]}|{' ' * WIDTHS[1]}| {'Source class':<{span - 2}} |",
    file=args.output,
)
print(
    f"| Observing |{' ' * WIDTHS[1]}+" + "+".join("-" * w for w in WIDTHS[2:]) + "+",
    file=args.output,
)
print(row(["run", "Network", *CLASSES]), file=args.output)
print(RULE.replace("-", "="), file=args.output)

# The runs that have happened first, then the forecasts.
for label, weighted, annual in [
    ("Detections over the observing time of the run, above rho_eq", False, False),
    ("Detections over the observing time of the run, redshift weighted", True, False),
    ("Annual number of detections above rho_eq", False, True),
    ("Annual number of detections above rho_eq, redshift weighted", True, True),
]:
    section(label, weighted, annual)


def paragraph(text):
    """A line of prose under the table, wrapped to the width of the table."""
    print(file=args.output)
    print(textwrap.fill(text, len(RULE)), file=args.output)


if args.population is not None:
    paragraph(
        f"Fraction of {args.population.name} in each class, which cancels "
        "above: " + ", ".join(f"{name} {fractions[name]:.4f}" for name in CLASSES) + "."
    )
