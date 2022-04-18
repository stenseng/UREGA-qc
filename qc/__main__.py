#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lars Stenseng.

@mail: lars@stenseng.net
"""

import argparse
import logging
from signal import SIGINT, SIGTERM, signal
from sys import exit
import georinex as gr

from qc.multipath import Multipath
# from qc.qc import qc
# from qc.slips import slips


def procSigint(signum, frame):
    logging.warning("Received SIGINT. Shutting down, Adjø!")
    exit(3)


def procSigterm(signum, frame):
    logging.warning("Received SIGTERM. Shutting down, Adjø!")
    exit(4)


signal(SIGINT, procSigint)
signal(SIGTERM, procSigterm)


parser = argparse.ArgumentParser()
parser.add_argument("-l",
                    "--logfile",
                    help="Log to file. Default output is terminal.")
parser.add_argument(
    "-v", "--verbosity",
    action="count", default=0, help="Increase verbosity level."
)
args = parser.parse_args()

logLevel = logging.ERROR
if args.verbosity == 1:
    logLevel = logging.WARNING
elif args.verbosity == 2:
    logLevel = logging.INFO
elif args.verbosity > 2:
    logLevel = logging.DEBUG
if args.logfile:
    logging.basicConfig(
        level=logLevel,
        filename=args.logfile,
        format="%(asctime)s;%(levelname)s;%(message)s",
    )
else:
    logging.basicConfig(level=logLevel,
                        format="%(asctime)s;%(levelname)s;%(message)s")

# Load obs file - used for multipath
obs = gr.load(
    'tests/test_data/Rinex3/KLSQ00GRL_R_20213070000_01D_15S_MO.rnx',
    tlim=['2021-11-03T10:32', '2021-11-03T11:32'])

# Load rinex 3 header
hdr = gr.rinexheader(
    'tests/test_data/Rinex3/KLSQ00GRL_R_20213070000_01D_15S_MO.rnx')

# Call multipath object
mptest = Multipath(obs, 'E', hdr)  # 'G'/'R'/'E'/'C' ('M' coming soon)
# Use function get_MP to get MP1 (more coming soon)
MP = mptest.get_MP()
