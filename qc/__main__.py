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
from matplotlib.pyplot import figure, show
import numpy as np
import matplotlib.pyplot as plt

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

# %%
# Load obs file and header (RINEX 3) - used for multipath
obs = gr.load(
    'tests/test_data/Rinex3/KLSQ00GRL_R_20213070000_01D_15S_MO.rnx',
    tlim=['2021-11-03T12:00', '2021-11-03T12:30'])
hdr = gr.rinexheader(
    'tests/test_data/Rinex3/KLSQ00GRL_R_20213070000_01D_15S_MO.rnx')
rnx_version = 3

# %%
# Load obs file and header (RINEX 2)
obs = gr.load("tests/test_data/Rinex2/klsq3070.21o", fast=False,
              tlim=['2021-11-03T12:00:00', '2021-11-03T13:00:30'])
hdr = gr.rinexheader("tests/test_data/Rinex2/klsq3070.21o")
rnx_version = 2

# %% Using Multipath class

MP_eq = 1  # Select MP equation 1/2/5

# Select observation codes: (rnx3/rnx2 examples)
codes = ['C1C', 'C2C', 'C5I', 'L1C', 'L2W']  # GPS test (rnx3)
# codes = ['C1P', 'C2C', 'C3I', 'L1C', 'L2C']  # GLONASS test (rnx3)
# codes = ['C1A', 'C8Q', 'C6A', 'L1C', 'L8I']  # Galileo test (rnx2)
# codes = ['C2X', 'C7I', 'C6I', 'L2I', 'L7I']  # BeiDou test (rnx3)
# codes = ['C1', 'C2', 'C5', 'L1', 'L2']  # GPS/GLONASS test (rnx2)
# codes = ['C1', 'C8', 'C6', 'L1', 'L8']  # Galileo test (rnx2)

# Call multipath object
mp = Multipath(obs,  # Observation file
               hdr,  # Header (from obs)
               'G',  # Constellation G/R/E/C/M
               MP_eq=MP_eq,  # MP equation
               codes=codes,  # Observation codes (don't use for Mixed!)
               rnx_version=rnx_version)  # Rinex version

# Use function get_MP to get MP1
MP = mp.get_MP()
