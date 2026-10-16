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

## Which SNRs the thresholds cut

The thresholds cut the network SNR of `coincs.dat`, as `bayestar-realize-coincs`
computes it: o2-uberbank templates from 25 Hz, with the `o4b_*_ref.txt` curves.
Both runs refer to these curves. What the `run` column selects is the injection
set, with its searches and the detectors observing at each injection time.

The injection release computes its optimal SNRs otherwise, with IMRPhenomXPHM
from 16 Hz and, for each month, the 10% quantile in each frequency bin of the
PSDs of that month. The same measurement on those SNRs gives thresholds higher
by 4-12% in O4a and 10-17% in O4b. LIGO-P2600188 tabulates both; only the values
here apply to a simulated SNR.

Within a run the threshold does not follow the detector sensitivity
(LIGO-P2600188, Section 4.3). The O4b rows, the default of
`detection-weights.py --rho-eq-run`, therefore apply unchanged to runs simulated
with other noise curves (IR1, O5), as long as their SNRs are computed with the
same templates and cutoff. This assumes that the searches do not change after
O4.

## How the values were measured

Not in this repository. They come from the GWTC-5.0 cumulative sensitivity
estimates, doi:10.5281/zenodo.19500052, which the simulations never read. For
the injections of one run, class and bin:

1.  take the smallest FAR over the searches of the run, and call an injection
    found below 1/yr;
2.  convert the optimal SNR of the injection in each detector observing at its
    GPS time to the one bayestar computes for the same source, by the ratio of
    the noise-weighted norms of the two templates, then recover a network SNR
    as `bayestar-realize-coincs` does, over sixteen noise realizations;
3.  `rho_eq` is the SNR above which the fraction of recovered values equals the
    found fraction. Equivalently, as many found injections sit below `rho_eq` as
    missed injections above it.

Both fractions are taken over the same injections, so no population weights
enter. Because the found fraction rises gradually with SNR, the threshold
reproduces the number of FAR-selected signals exactly only for the SNR
distribution of the injections; for signals falling as SNR^-4 above 8, the FAR
selection keeps 12% more. LIGO-P2600188 gives the method, the bins, the
bootstrap error and the checks.

The code is `scripts/snr_far_calibration/` of
<https://github.com/weizmannk/cbc-population-distributions>, at the commit named
in the header of `rho_eq.csv`. `snr_scale.py` computes the conversion of step 2
and `measure.py` the thresholds:

    $ python snr_scale.py --runs O4a O4b -o outputs/snr_scale.csv \
          --H1 ../../data/asd/o4b_h1_ref.txt \
          --L1 ../../data/asd/o4b_l1_ref.txt \
          --V1 ../../data/asd/o4b_v1_ref.txt
    $ SNRFAR_OUT=outputs/reference python measure.py --runs O4a O4b \
          --snr-scale outputs/snr_scale.csv

`export_thresholds.py` then writes this file, keeping the provenance lines and
dropping the two columns that record the closure of the measurement. It refuses
a table measured on the SNRs of the release. Reproducing the measurement needs
the 1.1 GB injection file and the GWOSC segments of each run, so the table is
kept as data here rather than as a `make` target.
