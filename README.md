# Observing scenarios

This repository contains all of the scripts and data to produce the compact
binary coalescence predictions for the LIGO/Virgo/KAGRA Observing Scenarios
Paper.

## Setting up your environment

1.  Set up your cluster environment.

    **[IGWN Grid]**: No special configuration necessary.

    **[NASA HECC]**:

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

    I use [Poetry] to manage the installation of the dependencies for this
    analysis on the cluster. Install Poetry by running the command in the
    [official Poetry installation instructions]:

        $ curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python

    Then log out, and log back in.

3.  Clone this repository (note, on [IGWN Grid], clone this in your home
    directory, but on [NASA HECC], clone this in your `/nobackup/$USER`
    directory):

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

**[IGWN Grid]**:

Run this command on the idle head node that has the highest core count,
because it uses Python multiprocessing. This takes about an hour per observing
scenario.

    $ make

**[NASA HECC]**:

Some of the `make` targets are CPU-intensive, and some need to access the
Internet. However, NASA HECC front-end nodes may not run CPU- or
memory-intensive jobs, and worker nodes cannot access the Internet. So at NASA,
the `make` command must be done in two steps.

1.  Download and pack the PSD files by running this command:

        $ make psds

2.  Run the rest of the `make` targets under PBS.

    *   Write this submit script to `make.pbs`:

            #!/bin/sh
            #PBS -Vkoed -l select=1:model=cas_ait -l walltime=10:00:00
            make

    *   Make the submit script executable:

            $ chmod +x make.pbs

    *   Submit the job (note, since it takes more than 8 hours to run, it
        needs to go to the `long` queue):

            $ qsub -q long make.pbs

## To run BAYESTAR

Run this on the head node to submit all of the BAYESTAR jobs to the cluster.

    for eventsfile in runs/*/*/events.xml.gz
        do bayestar-localize-coincs $eventsfile -o $(dirname $eventsfile)/allsky \
            --f-low 11 --cosmology
    done

**[IGWN Grid]**: Add the `--condor-submit` to the command line arguments for
`bayestar-localize-coincs` in order to submit the jobs to the cluster rather
than running them locally.

**[NASA HECC]**: Run under PBS. Use the following submit file:

    #/bin/sh
    #PBS -V -koed -l select=1:model=cas_ait -l walltime=10:00:00
    unset OMP_NUM_THREADS
    ...  # <-- put commands here 

## To tabulate BAYESTAR localization statistics

Once the above Condor jobs are done, run this on an idle head node with a high
core count to tabulate statistics for all of the sky maps.

    for eventsfile in runs/*/*/events.sqlite
        do ligo-skymap-stats -d $eventsfile -o $(dirname $eventsfile)/allsky.dat \
            $(find $(dirname $eventsfile)/allsky -name '*.fits' | sort -V) --cosmology --contour 20 50 90 -j
    done

**[IGWN Grid]**: Run on a high core count head node.

**[NASA HECC]**: Run under PBS. Use the same submit file commands as from the
previous step.

## To update human-readable figures and tables

Open the Jupyter notebook `plots-and-tables.ipynb` and run all cells.


[IGWN Grid]: https://computing.docs.ligo.org/guide/grid/
[NASA HECC]: https://www.nas.nasa.gov/hecc/
[Poetry]: https://python-poetry.org
[Poetry issue #1651]: https://github.com/python-poetry/poetry/issues/1651
[official Poetry installation instructions]: https://python-poetry.org/docs/#installation
