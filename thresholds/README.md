# SNR thresholds equivalent to a FAR cut

`rho_eq.csv` holds the network SNR threshold that selects as many signals as the
FAR < 1/yr cut of GWTC-5.0, per source class and detector-frame chirp-mass bin.
`measure.py` writes it and `scripts/detection-weights.py` reads it.

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
The injection release computes its SNRs with other templates and PSDs, which
give thresholds 4-17% higher, so its SNRs are converted before the measurement.
The `run` column selects the injections, with their searches and the detectors
observing at each injection time.

A simulated run takes the thresholds of the observing run it reproduces: the
O4a rows for the HL network of O4, the O4b rows for HLV. Within O4a the
threshold does not change with the detector sensitivity, so the O4b rows, the
default of `detection-weights.py`, also serve IR1 and O5, which assumes that
the searches do not change after O4.

## Using the thresholds

    $ make detection-weights

writes `runs/<run>/fullpop/detection-weights.dat` for every run already
simulated. `passes_rho_eq` is 1 when the network SNR of an event reaches the
threshold of its class and chirp-mass bin, and `redshift_weight` is
`(1 + z)**lamb`, which restores the redshift evolution of the FullPop rate.
Summing `passes_rho_eq`, or its product with `redshift_weight`, over a class
gives the detection count; `scripts/weighted-detections.py` does that sum and
turns it into `runs/weighted-detections.rst`, with the same rate prior and
credible interval as `runs/summary.rst`. `make -B detection-weights` redoes the
weights after a run has been simulated again.

## How the values were measured

For the injections of one run, class and bin of the GWTC-5.0 cumulative
sensitivity estimates (doi:10.5281/zenodo.19500052):

1.  an injection is found when its smallest FAR over the searches is below
    1/yr;
2.  its optimal SNR in each detector observing at its GPS time is converted to
    the one bayestar computes (`snr_scale.py`), and its network SNR is drawn as
    `bayestar-realize-coincs` does, over sixteen noise realizations;
3.  `rho_eq` is the SNR above which the fraction of drawn values equals the
    found fraction, so that as many found injections lie below it as missed
    injections above it.

No population weights enter. Because the found fraction rises gradually with
SNR, the threshold reproduces the number of FAR-selected signals exactly only
for the SNR distribution of the injections; for signals falling as SNR^-4 above
8, the FAR selection keeps 12% more. `config.py` lists every input and choice.

The measurement is not part of `make`. It only needs to be redone for a new
injection set. From the root of the repository, download its inputs into
`thresholds/data/`, which git ignores:

    $ mkdir -p thresholds/data/psds-o1234ab
    $ curl -L -o thresholds/data/mixture-semi_o1_o2-real_o3_o4a_o4b-cartesian_spins_20260410130052UTC-clipped.hdf \
          https://zenodo.org/records/19500052/files/mixture-semi_o1_o2-real_o3_o4a_o4b-cartesian_spins_20260410130052UTC-clipped.hdf
    $ curl -L https://zenodo.org/records/19500052/files/psds-o1234ab.tar.gz | tar -xz -C thresholds/data/psds-o1234ab
    $ make o4b_h1_ref.txt o4b_l1_ref.txt o4b_v1_ref.txt
    $ python thresholds/measure.py

`measure.py` rewrites `rho_eq.csv`; `-o FILE` writes elsewhere. It needs
`SEOBNRv4ROM_v3.0.hdf5` in `LAL_DATA_PATH` and the `gwosc` package, which
downloads the observing segments.
