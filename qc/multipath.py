#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lars Stenseng
@mail: lars@stenseng.net
"""

import numpy as np
from qc.__version__ import __version__

class Multipath:
    def getmultipath(obs, th, signals, snr, slips):
    ### Inputs:
    # obs: observations
    # th: threshold
    # signals: the used signal combination
    # snr: signal to noise ratio
    # slips: number of slips
    # 
    ### Outputs:
    # sv: satellite value (obj.sv[0])
    # signals: display the signal combination
    # thnr: number of times the multipath exceeded the threshold
    # tstamps: timestamps
    # mp: multipath
    #
    
    ### Testing
    sv = 0
    thnr = 0
    tstamps = 0
    mp = 0
    res = np.array([sv, [signals[0], signals[1]], thnr, [tstamps, mp]])
    
    return res

