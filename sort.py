#!/usr/bin/env python
from astropy.table import Table
import sys

_, filename = sys.argv
table = Table.read(filename, format='ascii')
table['sortindex'] = [int(_.split(':')[-1]) for _ in table['coinc_event_id']]
table.sort('sortindex')
del table['sortindex']
table.write(filename, format='ascii.tab', overwrite=True)
