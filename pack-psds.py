#!/usr/bin/env python
"""Pack ASCII injection files into a psd.xml file."""
from ligo.skymap.tool import ArgumentParser, FileType, register_to_xmldoc
from argparse import SUPPRESS
import lal

# Command line interface
detector_names = [d.frDetector.prefix for d in lal.CachedDetectors]
detector_long_names = [d.frDetector.name for d in lal.CachedDetectors]
parser = ArgumentParser()
parser.add_argument(
    '-o', '--output', metavar='OUT.xml[.gz]', type=FileType('wb'),
    default='-', help='Name of output file [default: stdout]')
for name, long_name in zip(detector_names, detector_long_names):
    parser.add_argument(
        '--' + name, metavar='PSD.txt', type=FileType('r'), default=SUPPRESS,
        help='PSD function for {0} detector'.format(long_name))
args = parser.parse_args()

# Late imports
import os
import glue.ligolw.utils
import lal.series
import numpy as np

psds = {}
for name in detector_names:
    psd_file = getattr(args, name, None)
    if psd_file is None:
        continue

    f, asd = np.loadtxt(psd_file).T
    psd = np.square(asd)

    f0 = 10.0
    fmax = 4096.0
    df = 1.0

    fgrid = np.arange(f0, fmax, df)
    series = lal.CreateREAL8FrequencySeries(
        psd_file.name, 0, f0, df, lal.SecondUnit, len(fgrid))
    series.data.data = np.exp(np.interp(np.log(fgrid), np.log(f), np.log(psd)))

    psds[name] = series

xmldoc = lal.series.make_psd_xmldoc(psds)
register_to_xmldoc(xmldoc, parser, args)

with glue.ligolw.utils.SignalsTrap():
    glue.ligolw.utils.write_fileobj(
        xmldoc, args.output,
        gz=(os.path.splitext(args.output.name)[-1] == '.gz'))
