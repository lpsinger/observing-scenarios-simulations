# Observing scenarios

This repository contains all of the scripts and data to produce the compact
binary coalescence predictions for the LIGO/Virgo/KAGRA Observing Scenarios
Paper.

## Setting up your environment

I use [Poetry](https://python-poetry.org) to manage the installation of the
dependencies for this analysis on the cluster.

First, make sure that the following line is in your shell login script (`.bash_profile`, `.profile`, or equivalent):

    export PATH=$HOME/.local/bin:$PATH

Log out and log back in if necessary. Clone this repository:

    git clone https://github.com/lpsinger/observing-scenarios-simulations
    cd observing-scenarios-simulations

Next, make sure that pip is up to date, and install poetry:

    pip3 install --upgrade --user pip
    pip install --user poetry

Note, due to https://github.com/python-poetry/poetry/issues/1651, you need to also make sure that pip is up to date _inside_ the poetry environment:

    poetry run pip install --upgrade pip

Finally, install the dependencies and launch a shell inside the poetry env:

    poetry install
    poetry shell

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
