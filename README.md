# Observing scenarios

This repository contains all of the scripts and data to produce the compact
binary coalescence predictions for the LIGO/Virgo/KAGRA Observing Scenarios
Paper.

## Setting up your environment

I use [Poetry] to manage the installation of the dependencies for this analysis
on the cluster. First, install Poetry by following the
[official Poetry installation instructions].

Next, clone this repository:

    $ git clone https://github.com/lpsinger/observing-scenarios-simulations
    $ cd observing-scenarios-simulations

Next, make sure that pip is up to date, and install poetry:

    $ pip3 install --upgrade --user pip
    $ pip install --user poetry

Note, due to [Poetry issue #1651], you need to also make sure that pip is up to
date _inside_ the poetry environment:

    $ poetry run pip install --upgrade pip

Finally, install the dependencies and launch a shell inside the poetry env:

    $ poetry install
    $ poetry shell

### Cluster-specific environment instructions

#### [NASA Pleiades]

1.  Login to a Pleiades front end node. Add the following to your `~/.profile`
    script:

        module use -a /nasa/modulefiles/testing
        module load python3-intel/2020.0.014

    Log out, and log back in.

2.  Install Poetry following the
    [official Poetry installation instructions for Linux]:

        $ curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python

    Log out, and log back in again.

3.  Clone this repository in your scratch directrory:

        $ cd /nobackup/$USER
        $ git clone https://github.com/lpsinger/observing-scenarios-simulations
        $ cd observing-scenarios-simulations

4.  Update Pip within the Poetry environment:

        $ poetry run pip install --upgrade pip

5.  Finally, install the dependencies and launch a shell inside the poetry env:

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


[Poetry]: https://python-poetry.org
[Poetry issue #1651]: https://github.com/python-poetry/poetry/issues/1651
[official Poetry installation instructions]: https://python-poetry.org/docs/#installation
[official Poetry installation instructions for Linux]: https://python-poetry.org/docs/#osx-linux-bashonwindows-install-instructions
[NASA Pleiades]: https://www.nas.nasa.gov/hecc/
