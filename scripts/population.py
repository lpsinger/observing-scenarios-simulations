#!/usr/bin/env python
"""Convert a FullPop-4.0 distribution to suitable format for bayestar-inject."""

from argparse import ArgumentParser

import numpy as np
from astropy.table import Table
from ligo.skymap.bayestar.ez_emcee import ez_emcee
from popsummary.popresult import PopulationResult
from scipy.interpolate import RegularGridInterpolator

parser = ArgumentParser()
parser.add_argument("input", help="PopSummary input file")
parser.add_argument("output", help="HDF5 output file")
parser.add_argument("--mode", choices=("fullpop", "pixelpop"))
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

np.random.seed(args.seed)
n = 1_000_000

# Read distribution.
popresult = PopulationResult(args.input)

if args.mode == "fullpop":
    (mass_1, mass_2), rates = popresult.get_rates_on_grids(
        "primary_mass_secondary_mass_joint_full_posterior"
    )
    np.testing.assert_array_equal(mass := mass_1, mass_2)
    mass_grid_length = len(mass)
    # FIXME: The dimensions of the FullPop grid in the popsummary file are just wrong.
    rates = rates.reshape(mass_grid_length, mass_grid_length, -1)
    median_rates = np.median(rates, axis=2)
    # Convert to from d^2/(dm1dm2) to d^2p/(dlog(m1)dlog(m2))
    median_rates *= np.prod(np.meshgrid(mass, mass), axis=0)
    log_mass = np.log(mass)
elif args.mode == "pixelpop":
    (log_mass_1, log_mass_2), rates = popresult.get_rates_on_grids(
        "joint_pixelpop_rate"
    )
    median_rates = np.median(rates, axis=0)
    mass_grid_length = int(np.sqrt(median_rates.size))
    median_rates = median_rates.reshape(mass_grid_length, mass_grid_length)
    log_mass_1 = log_mass_1.reshape(mass_grid_length + 1, mass_grid_length + 1)
    log_mass_2 = log_mass_2.reshape(mass_grid_length + 1, mass_grid_length + 1)
    np.testing.assert_array_equal(log_mass := log_mass_1[:, 0], log_mass_2[0, :])
    # Get bin centers from bin edges
    log_mass = 0.5 * (log_mass[1:] + log_mass[:-1])
else:
    raise RuntimeError("This line must not be reached.")

# Fill out entire grid to make it easier to sample correctly near the diagonal.
median_rates = np.where(median_rates == 0, median_rates.T, median_rates)

interp = RegularGridInterpolator((log_mass, log_mass), np.log(median_rates))


def logp(params):
    return interp(params.T)


# Sample masses
m1, m2 = np.exp(
    ez_emcee(
        interp,
        np.full(2, log_mass.min()),
        np.full(2, log_mass.max()),
        vectorize=True,
        nwalkers=20,
        nindep=50000,
    )
).T

# Swap masses as necessary to ensure m1 >= m2
good = m1 >= m2
m1, m2 = np.where(good, m1, m2), np.where(good, m2, m1)

# Keep only events with eta >= 0.01 because that is what SEOBNRv4ROM supports.
eta = m1 * m2 / np.square(m1 + m2)
good = eta >= 0.01
m1 = m1[good]
m2 = m2[good]
n = len(m1)

# From GWTC-4 population paper (https://arxiv.org/abs/2508.18083v2) Section 4.1:
#
# > The BGP analysis fixes...the spin distribution to be uniform in magnitude
# (again truncated over χ∈[0,0.4] for NS masses) and isotropic in orientation.
ns_max_mass = 2.5
x1 = np.random.uniform(0, np.where(m1 <= ns_max_mass, 0.4, 1))
x2 = np.random.uniform(0, np.where(m2 <= ns_max_mass, 0.4, 1))
cos1 = np.random.uniform(0, 1, n)
cos2 = np.random.uniform(0, 1, n)

Table({"mass1": m1, "mass2": m2, "spin1z": x1 * cos2, "spin2z": x2 * cos2}).write(
    args.output, overwrite=True
)
