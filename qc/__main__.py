#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lars Stenseng
@mail: lars@stenseng.net
"""

import argparse
import logging
from signal import SIGINT, SIGTERM, signal
from sys import exit

from qc.multipath import multipath
from qc.qc import qc
from qc.slips import slips


def procSigint(signum, frame):
    logging.warning("Received SIGINT. Shutting down, Adjø!")
    exit(3)


def procSigterm(signum, frame):
    logging.warning("Received SIGTERM. Shutting down, Adjø!")
    exit(4)


signal(SIGINT, procSigint)
signal(SIGTERM, procSigterm)


parser = argparse.ArgumentParser()
parser.add_argument("-l", "--logfile", help="Log to file. Default output is terminal.")
parser.add_argument(
    "-v", "--verbosity", action="count", default=0, help="Increase verbosity level."
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
    logging.basicConfig(level=logLevel, format="%(asctime)s;%(levelname)s;%(message)s")
