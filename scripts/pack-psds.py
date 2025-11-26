#!/usr/bin/env python
"""Pack ASCII injection files into a psd.xml file."""

import os
from argparse import SUPPRESS, FileType

import glue.ligolw.utils
import lal
import lal.series
import numpy as np
from ligo.skymap.tool import ArgumentParser, register_to_xmldoc


def read_multicolumn_txt_file(filepath, column_name):
    """
    Read .txt file and extract specific column by name.

    Parameters:
        filepath: Path to the .txt file
        column_name: Column name (e.g., 'O5aStrain', 'O5bStrain', 'O5cStrain')

    Returns:
        Tuple: (frequency, asd_data)
    """
    data = np.genfromtxt(filepath, names=True)

    if column_name not in data.dtype.names:
        raise ValueError(
            f"Column '{column_name}' not found. Available: {', '.join(data.dtype.names)}"
        )
    frequency = data[data.dtype.names[0]]
    asd_data = data[column_name]
    return frequency, asd_data


# Command line interface
detector_names = [d.frDetector.prefix for d in lal.CachedDetectors]
detector_long_names = [d.frDetector.name for d in lal.CachedDetectors]
parser = ArgumentParser()
parser.add_argument(
    "-o",
    "--output",
    metavar="OUT.xml[.gz]",
    type=FileType("wb"),
    default="-",
    help="Name of output file [default: stdout]",
)
for name, long_name in zip(detector_names, detector_long_names):
    parser.add_argument(
        "--" + name,
        metavar="PSD.txt",
        type=FileType("r"),
        default=SUPPRESS,
        help="PSD function for {0} detector".format(long_name),
    )
    parser.add_argument(
        f"--{name}-column",
        metavar="COLUMN_NAME",
        type=str,
        default=None,
        help="Column name to read from {0} PSD file (e.g., O5aStrain)".format(
            long_name
        ),
    )

args = parser.parse_args()

psds = {}
for name in detector_names:
    psd_file = getattr(args, name, None)
    if psd_file is None:
        continue

    # Get column argument for this detector
    column = getattr(args, f"{name}_column", None)

    if column is not None:
        # Read specific column from multi-column file
        f, asd = read_multicolumn_txt_file(psd_file.name, column)
    else:
        # Standard mode: read 2-column file
        f, asd = np.loadtxt(psd_file).T

    # Process PSD data
    psd = np.square(asd)

    f0 = 10.0
    fmax = 4096.0
    df = 1.0

    fgrid = np.arange(f0, fmax, df)
    series = lal.CreateREAL8FrequencySeries(
        (psd_file.name).split(".")[0], 0, f0, df, lal.SecondUnit, len(fgrid)
    )
    series.data.data = np.exp(np.interp(np.log(fgrid), np.log(f), np.log(psd)))

    psds[name] = series

xmldoc = lal.series.make_psd_xmldoc(psds)
register_to_xmldoc(xmldoc, parser, args)

with glue.ligolw.utils.SignalsTrap():
    glue.ligolw.utils.write_fileobj(
        xmldoc, args.output, gz=(os.path.splitext(args.output.name)[-1] == ".gz")
    )
