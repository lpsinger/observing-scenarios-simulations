# Observing scenarios

This repository contains all of the scripts and data to produce the compact
binary coalescence predictions for the LIGO/Virgo/KAGRA Observing Scenarios
Paper.

## Setting up your environment

1.  Set up your cluster environment.

    **[IGWN Grid]**: No special configuration necessary.

    **[NASA Pleiades]**:

    *   Login to a Pleiades front end node.

    *   Install lalsuite-extra by running these commands:

            $ cd /nobackup/$USER
            $ curl -O https://software.igwn.org/lscsoft/source/lalsuite-extra-1.3.0.tar.gz
            $ tar xf lalsuite-extra-1.3.0.tar.gz
            $ cd lalsuite-extra-1.3.0
            $ ./configure --prefix=$HOME/.local
            $ make install

    *   Add the following to your `~/.profile` script:

            export LAL_DATA_PATH=$HOME/.local/share/lalsimulation
            module use -a /nasa/modulefiles/testing
            module load python3-intel/2020.0.014

    *   Log out, and log back in.

2.  Install Poetry.

    I use [Poetry] to manage the installation of the dependencies for this analysis on the cluster. Install Poetry by running the command in the
    [official Poetry installation instructions]:

        $ curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python

    Then log out, and log back in.

3.  Clone this repository (note, on [NASA Pleiades], clone this in your
    `/nobackup/$USER` directory):

        $ git clone https://github.com/lpsinger/observing-scenarios-simulations
        $ cd observing-scenarios-simulations

4.  Install dependencies with Poetry.

    Note, due to [Poetry issue #1651], you need to also make sure that pip is
    up to date _inside_ the poetry environment:

        $ poetry run pip install --upgrade pip

    Finally, install the dependencies and launch a shell inside the poetry env:

        $ poetry install
        $ poetry shell

## To generate source populations

Run this on the idle head node that has the highest core count,
because it uses Python multiprocessing.
This takes about 20-30 minutes per observing scenario.

    make

## To run BAYESTAR

Run this on the head node to submit all of the BAYESTAR jobs to the cluster.

    for subdir in runs/*/*
        do bayestar-localize-coincs $subdir/events.xml.gz -o $subdir/allsky \
            --f-low 11 --cosmology --condor-submit
    done

## To tabulate BAYESTAR localization statistics

Once the above Condor jobs are done, run this on an idle head node with a high
core count to tabulate statistics for all of the sky maps.

    for subdir in runs/*/*
        do ligo-skymap-stats -d $subdir/events.sqlite -o $subdir/allsky.dat \
            $subdir/allsky/\*.fits --cosmology --contour 20 50 90 -j
    done

The `ligo-skymap-stats` script outputs rows in arbitrary order. To avoid large
git diffs, run this command to sort the rows in the data files:

    for f in runs/*/*/allsky.dat
        do ./sort.py $f
    done

## To update human-readable figures and tables

Open the Jupyter notebook `plots-and-tables.ipynb` and run all cells.


[IGWN Grid]: https://computing.docs.ligo.org/guide/grid/
[NASA Pleiades]: https://www.nas.nasa.gov/hecc/
[Poetry]: https://python-poetry.org
[Poetry issue #1651]: https://github.com/python-poetry/poetry/issues/1651
[official Poetry installation instructions]: https://python-poetry.org/docs/#installation
