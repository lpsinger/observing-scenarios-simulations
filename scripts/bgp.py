#!/usr/bin/env python
"""Convert a FullPop-4.0 distribution to suitable format for bayestar-inject."""

from argparse import ArgumentParser

import numpy as np
from astropy.table import Table
from popsummary.popresult import PopulationResult

parser = ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

np.random.seed(args.seed)
n = 1_000_000

# Get maximum a posteriori hyperparameters values.
popresult = PopulationResult(args.input)
(m1, m2), rates = popresult.get_rates_on_grids("ppd_primary_and_secondary_mass")

np.testing.assert_array_equal(m1, m2)
log_m = np.log(m1)

# Fill out entire grid to make it easier to sample correctly near the diagonal.
rates = np.where(rates == 0, rates.T, rates)

# Sample masses
p = rates[:-1, :-1]
p /= p.sum()
i, j = np.unravel_index(np.random.choice(p.size, n, replace=True, p=p.ravel()), p.shape)
m1 = np.exp(np.random.uniform(log_m[i], log_m[i + 1]))
m2 = np.exp(np.random.uniform(log_m[j], log_m[j + 1]))

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
