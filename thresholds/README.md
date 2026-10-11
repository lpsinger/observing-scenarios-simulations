# SNR thresholds equivalent to a FAR cut

`rho_eq.csv` holds the network SNR threshold that selects as many signals as the
FAR < 1/yr cut of the GWTC catalogues, per source class and detector-frame
chirp-mass bin. `scripts/detection-weights.py` reads it; nothing else does.

| Column | Content |
|---|---|
| `run` | observing run of the injections, `O4a` or `O4b` |
| `source_class` | `BNS`, `NSBH` or `BBH`, from the source-frame masses with a 2.5 Msun NS/BH boundary |
| `chirp_mass_min`, `chirp_mass_max` | bin edges in detector-frame chirp mass [Msun] |
| `n_injections` | injections in the bin |
| `n_found` | of those, the ones below FAR = 1/yr in at least one search |
| `rho_eq` | the threshold [network SNR] |
| `error` | bootstrap uncertainty on `rho_eq` |

Both runs refer to the sensitivity of the simulations, the `o4b_*_ref.txt`
curves. What the `run` column selects is the injection set, with its searches
and the detectors observing at each injection time.

## How the values were measured

Not in this repository. They come from the GWTC-5.0 cumulative sensitivity
estimates, doi:10.5281/zenodo.19500052, which the simulations never read. For
the injections of one run, class and bin:

1.  take the smallest FAR over the searches of the run, and call an injection
    found below 1/yr;
2.  recover a network SNR as `bayestar-realize-coincs` does, from the optimal
    SNR of the injection in each detector observing at its GPS time, rescaled
    to the reference curves of the simulations, over sixteen noise
    realizations;
3.  `rho_eq` is the SNR above which the fraction of recovered values equals the
    found fraction. Equivalently, as many found injections sit below `rho_eq` as
    missed injections above it.

Both fractions are taken over the same injections, so the population they were
drawn from cancels and no weights enter. LIGO-P2600188 gives the method, the
bins, the bootstrap error and the checks on the noise model and the binning.

The code is `scripts/snr_far_calibration/` of
<https://github.com/weizmannk/cbc-population-distributions>, where `measure.py`
runs the measurement and `threshold.py` holds the noise model and `rho_eq`. This
file is the output of

    $ SNRFAR_OUT=outputs/reference python measure.py --runs O4a O4b \
          --reference ../../data/asd

without the two columns that record the closure of the measurement. Reproducing
it needs the 1.1 GB injection file and the GWOSC segments of each run, so the
table is kept as data here rather than as a `make` target.
