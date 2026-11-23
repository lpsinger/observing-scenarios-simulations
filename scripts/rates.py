"""Poisson rate with a log-normal prior, and the formatting of its interval.

Reproduced from the notebook `plots-and-tables.ipynb` of this repository, so
that a script can produce the same numbers without running it. The notebook
still carries its own copies.
"""

import numpy as np
from scipy import integrate, optimize, special, stats


@np.vectorize
def poisson_lognormal_rate_cdf(k, mu, sigma):
    """CDF of a Poisson count whose rate has a log-normal prior."""
    lognorm_pdf = stats.lognorm(s=sigma, scale=np.exp(mu)).pdf

    def func(lam):
        return special.gammaincc(k + 1, lam) * lognorm_pdf(lam)

    cdf, _ = integrate.quad(func, 0, np.inf, epsabs=0)
    return cdf


@np.vectorize
def poisson_lognormal_rate_quantiles(p, mu, sigma):
    """Quantiles of a Poisson count whose rate has a log-normal prior.

    The count is treated as a continuous variable so that a root finder can be
    used; ``mu`` and ``sigma`` describe the log of the rate.
    """

    def func(k):
        return poisson_lognormal_rate_cdf(k, mu, sigma) - p

    if func(0) >= 0:
        return 0
    return optimize.root_scalar(func, bracket=[0, 1e6]).root


def format_with_errorbars(mid, lo, hi):
    """``(mid, minus, plus)`` as strings, rounded to the smaller error bar."""
    plus, minus = hi - mid, mid - lo
    smallest = min(max(0, plus), max(0, minus))
    if smallest == 0:
        return str(mid), "0", "0"

    decimals = 1 - int(np.floor(np.log10(smallest)))
    if all(np.issubdtype(type(_), np.integer) for _ in (mid, lo, hi)):
        decimals = min(decimals, 0)

    plus, minus, mid = np.round([plus, minus, mid], decimals)
    fstring = f"%.0{decimals:d}f" if decimals > 0 else "%d"
    return [fstring % _ for _ in [mid, minus, plus]]
