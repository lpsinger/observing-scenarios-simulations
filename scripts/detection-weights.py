#!/usr/bin/env python
"""Tabulate two optional per-event weights for one run and population: the
network SNR threshold that selects as many signals as the FAR < 1/yr cut of the
GWTC catalogues, and the weight of the FullPop redshift evolution. Nothing in
the default products reads this file."""

from argparse import ArgumentParser

import numpy as np
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
from astropy.table import Table, join
from popsummary.popresult import PopulationResult

parser = ArgumentParser()
parser.add_argument("injections", help="injections.dat input file")
parser.add_argument("coincs", help="coincs.dat input file")
parser.add_argument("thresholds", help="rho_eq.csv input file")
parser.add_argument("output", help="ASCII output file")
parser.add_argument(
    "--rho-eq-run",
    default="O4b",
    help="observing run whose thresholds to apply [default: %(default)s]",
)
parser.add_argument(
    "--popsummary", help="PopSummary file of the FullPop model, for the redshift weight"
)
args = parser.parse_args()

# NS/BH boundary of the GWTC-5.0 injections on which the thresholds were measured.
ns_max_mass = 2.5

injections = Table.read(args.injections, format="ascii.fast_tab")
coincs = Table.read(args.coincs, format="ascii.fast_tab")
coincs.rename_column("coinc_event_id", "simulation_id")
table = join(injections, coincs, keys="simulation_id")

# Planck15 is the cosmology of bayestar-inject.
redshift = z_at_value(cosmo.luminosity_distance, table["distance"] * u.Mpc).to_value(
    u.dimensionless_unscaled
)

# The masses in injections.dat are in the detector frame. The thresholds are
# binned in detector-frame chirp mass, the classes are in the source frame.
chirp_mass = (table["mass1"] * table["mass2"]) ** 0.6 / (
    table["mass1"] + table["mass2"]
) ** 0.2
mass1 = table["mass1"] / (1 + redshift)
mass2 = table["mass2"] / (1 + redshift)
source_class = np.where(
    np.maximum(mass1, mass2) < ns_max_mass,
    "BNS",
    np.where(np.minimum(mass1, mass2) < ns_max_mass, "NSBH", "BBH"),
)

thresholds = Table.read(args.thresholds, format="ascii.csv", comment="#")
thresholds = thresholds[thresholds["run"] == args.rho_eq_run]
rho_eq = np.empty(len(table))
for name in np.unique(source_class):
    rows = thresholds[thresholds["source_class"] == name]
    if not len(rows):
        raise ValueError(f"no {name} thresholds for run {args.rho_eq_run}")
    rows.sort("chirp_mass_min")
    keep = source_class == name
    # Events outside the bins take the nearest bin.
    index = np.clip(
        np.searchsorted(rows["chirp_mass_min"], chirp_mass[keep], side="right") - 1,
        0,
        len(rows) - 1,
    )
    rho_eq[keep] = rows["rho_eq"][index]

result = Table(
    {
        "simulation_id": table["simulation_id"],
        "redshift": redshift,
        "rho_eq": rho_eq,
        "passes_rho_eq": (table["snr"] >= rho_eq).astype(int),
    }
)
comments = [f"thresholds: {args.thresholds}, run {args.rho_eq_run}"]

if args.popsummary is not None:
    popresult = PopulationResult(args.popsummary)
    if "lamb" not in popresult.get_metadata("hyperparameters"):
        raise ValueError(f"no lamb hyperparameter in {args.popsummary}")
    # The merger rate of FullPop evolves as (1 + z)**lamb, whereas
    # bayestar-inject draws sources at a constant comoving rate.
    lamb = np.median(popresult.get_hyperparameter_samples(hyperparameters="lamb"))
    result["redshift_weight"] = (1 + redshift) ** lamb
    comments.append(f"redshift_weight = (1 + redshift)**lamb, lamb = {lamb:.4f}")

result.meta["comments"] = comments
result.write(args.output, format="ascii.fast_tab", overwrite=True)
