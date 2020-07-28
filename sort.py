#!/usr/bin/env python
from astropy.table import Table
import sys

_, filename = sys.argv
table = Table.read(filename, format='ascii')
table.sort('coinc_event_id')
table.write(filename, format='ascii.tab', overwrite=True)
