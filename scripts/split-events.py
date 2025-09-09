#!/usr/bin/env python
"""Split events.xml.gz file into many files with one event per file."""

import logging
from argparse import ArgumentParser
from itertools import groupby
from operator import attrgetter
from pathlib import Path

from igwn_ligolw.ligolw import LIGO_LW, Document, Param
from igwn_ligolw.lsctables import (
    CoincDefTable,
    CoincMapTable,
    CoincTable,
    ProcessParamsTable,
    ProcessTable,
    SnglInspiralTable,
    TimeSlideTable,
)
from igwn_ligolw.utils import load_filename, write_filename
from lalinspiral.thinca import InspiralCoincDef
from tqdm.auto import tqdm

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

table_classes = [
    CoincDefTable,
    CoincMapTable,
    CoincTable,
    ProcessParamsTable,
    ProcessTable,
    SnglInspiralTable,
    TimeSlideTable,
]
table_classes_to_copy = [
    CoincDefTable,
    ProcessParamsTable,
    ProcessTable,
    TimeSlideTable,
]


# Parse command line arguments
parser = ArgumentParser()
parser.add_argument("input")
parser.add_argument("outdir", type=Path)
args = parser.parse_args()

log.info('reading "%s"', args.input)
xmldoc = load_filename(args.input)
tables = {cls: cls.get_table(xmldoc) for cls in table_classes}

# These tables are to be copied without modification
tables_to_copy = [tables[cls] for cls in table_classes_to_copy]

log.info("indexing")

# Find coinc_def_id for sngl_inspiral<->sngl_inspiral coincs
coinc_def_id = tables[CoincDefTable].get_coinc_def_id(
    InspiralCoincDef.search, InspiralCoincDef.search_coinc_type
)

# Create a dictionary mapping coinc_event_id to lists of coinc_event_map rows
keyfunc = attrgetter("coinc_event_id")
coinc_map_dict = {
    key: tuple(items)
    for key, items in groupby(sorted(tables[CoincMapTable], key=keyfunc), key=keyfunc)
}

# Create a dictionary mapping event_id to sngl_inspiral rows
sngl_dict = {row.event_id: row for row in tables[SnglInspiralTable]}

# Create a dictionary mapping event_id to SNR time series
snr_series_dict = {
    param.value: param.parentNode for param in Param.getParamsByName(xmldoc, "event_id")
}

log.info('writing new files to directory "%s"', args.outdir)
args.outdir.mkdir(exist_ok=True)

for coinc in tqdm(tables[CoincTable]):
    if coinc.coinc_def_id != coinc_def_id:
        continue

    new_xmldoc = Document()
    new_xmldoc.appendChild(new_ligolw := LIGO_LW())
    for table in tables_to_copy:
        new_ligolw.appendChild(table)

    new_ligolw.appendChild(new_coinc_table := CoincTable.new())
    new_coinc_table.append(coinc)

    new_ligolw.appendChild(new_coinc_map_table := CoincMapTable.new())
    new_coinc_map_table.extend(coinc_map_dict[coinc.coinc_event_id])

    new_ligolw.appendChild(new_sngl_table := SnglInspiralTable.new())
    new_sngl_table.extend(sngl_dict[row.event_id] for row in new_coinc_map_table)

    for row in new_coinc_map_table:
        new_ligolw.appendChild(snr_series_dict[row.event_id])

    new_filename = str(args.outdir / f"{coinc.coinc_event_id}.xml.gz")
    write_filename(new_xmldoc, new_filename)
