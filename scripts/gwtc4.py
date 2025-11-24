#!/usr/bin/env python
"""Convert Farah's distribution to suitable format for bayestar-inject."""

from argparse import ArgumentParser

from astropy.table import Table

parser = ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
args = parser.parse_args()

data = Table.read(args.input)
Table(
    {
        "mass1": data["mass_1"],
        "mass2": data["mass_2"],
        "spin1z": data["a_1"] * data["cos_tilt_1"],
        "spin2z": data["a_2"] * data["cos_tilt_2"],
    }
).write(args.output, overwrite=True)
