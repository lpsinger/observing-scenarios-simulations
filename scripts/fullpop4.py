#!/usr/bin/env python
"""Convert a FullPop-4.0 distribution to suitable format for bayestar-inject."""

from argparse import ArgumentParser

import numpy as np
import pandas as pd
from astropy.table import Table
from gwpopulation.utils import truncnorm
from ligo.skymap.bayestar.ez_emcee import ez_emcee
from popsummary.popresult import PopulationResult

parser = ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()

np.random.seed(args.seed)

# Get maximum a posteriori hyperparameters values.
popresult = PopulationResult(args.input)
df = pd.DataFrame(
    popresult.get_hyperparameter_samples(),
    columns=popresult.get_metadata("hyperparameters"),
)
hyper = df.iloc[(df.log_likelihood + df.log_prior).idxmax()]

# FIXME: These parameters are not included in the public release.
m_min = 0.5
m_max = 350


def lopass(m, m_crit, eta):
    return 1 / (1 + (m / m_crit) ** eta)


def hipass(m, m_crit, eta):
    return 1 - lopass(m, m_crit, eta)


def notch(m, m_lo, m_hi, eta_lo, eta_hi, amp):
    return 1 - amp * hipass(m, m_lo, eta_lo) * lopass(m, m_hi, eta_hi)


def mass_dist_1d(m):
    truncnorm1 = truncnorm(m, hyper.mu1, hyper.sig1, m_max, m_min)
    truncnorm2 = truncnorm(m, hyper.mu2, hyper.sig2, m_max, m_min)
    hi = hipass(m, hyper.NSmin, hyper.n0)
    lo = lopass(m, hyper.BHmax, hyper.n5)
    notch1 = notch(
        m,
        hyper.NSmax,
        hyper.BHmin,
        hyper.n1,
        hyper.n2,
        hyper.A,
    )
    notch2 = notch(
        m,
        hyper.UPPERmin,
        hyper.UPPERmax,
        hyper.n3,
        hyper.n4,
        hyper.A2,
    )
    powerlaw = np.piecewise(
        m,
        (
            m < hyper.NSmax,
            (m >= hyper.NSmax) & (m < hyper.BHmin),
        ),
        (
            lambda m: m**hyper.alpha_1,
            lambda m: m**hyper.alpha_dip
            * hyper.NSmax ** (hyper.alpha_1 - hyper.alpha_dip),
            lambda m: m**hyper.alpha_2
            * hyper.NSmax ** (hyper.alpha_1 - hyper.alpha_dip)
            * hyper.BHmin ** (hyper.alpha_dip - hyper.alpha_2),
        ),
    )
    return (
        (1 + hyper.mix1 * truncnorm1 + hyper.mix2 * truncnorm2)
        * notch1
        * notch2
        * hi
        * lo
        * powerlaw
    )


def mass_pairing(m1, m2):
    return (m2 / m1) ** np.where(m2 < 5, hyper.beta_pair_1, hyper.beta_pair_2)


def mass_dist_2d(m1, m2):
    return (m1 >= m2) * mass_dist_1d(m1) * mass_dist_1d(m2) * mass_pairing(m1, m2)


def spin_mag_dist_1d(x, m):
    return truncnorm(x, hyper.mu_chi, hyper.sigma_chi, np.where(m < 2.5, 0.4, 1), 0)


def spin_tilt_dist_1d(cos):
    return truncnorm(cos, hyper.mu_spin, hyper.sigma_spin, 1, -1)


def spin_tilt_dist_2d(cos1, cos2):
    return (
        hyper.xi_spin * spin_tilt_dist_1d(cos1) * spin_tilt_dist_1d(cos2)
        + (1 - hyper.xi_spin) / 4
    )


def logp(params):
    m1, m2, x1, x2, cos1, cos2 = params.T
    return np.log(
        mass_dist_2d(m1, m2)
        * spin_mag_dist_1d(x1, m1)
        * spin_mag_dist_1d(x2, m2)
        * spin_tilt_dist_2d(cos1, cos2)
    )


m1, m2, x1, x2, cos1, cos2 = ez_emcee(
    logp,
    [m_min, m_min, 0, 0, -1, -1],
    [m_max, m_max, 1, 1, 1, 1],
    vectorize=True,
    nwalkers=20,
    nindep=50000,
).T

Table({"mass1": m1, "mass2": m2, "spin1z": x1 * cos1, "spin2z": x2 * cos2}).write(
    args.output, overwrite=True
)
