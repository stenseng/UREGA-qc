#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lars Stenseng
@mail: lars@stenseng.net
"""

import numpy as np
import georinex as gr
from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
# from qc.__version__ import __version__


class Multipath:

    def __init__(self, obs):  # , signals):
        self.obs = obs
        # self.signals = signals

    def get_MP(self):
        freq = self.load_constellation_freq()
        sv_all = self.sort_sat_types(self.obs)
        self.MP_G = self.append_MP_arrays(sv_all[0], self.obs, freq[0])
        # print(MP_G)
        self.plot_all_MP(sv_all[0], self.MP_G)
        return self.MP_G

    def load_constellation_freq(self):

        '''
        Function that loads frequency variables for each constellation

        Inputs: none
        Outputs: constellation arrays
        '''

        GPS_freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS
        GLONASS_freq = [1603.6875, 1247.3125, 1202.025]  # G1, G2, G3 for GLONASS
        Galileo_freq = [1575.42, 1278.75, 1191.795]  # E1, E6, E5 for Galileo (no L2)
        SBAS_freq = [1575.42, 1176.45]  # L1, L5 for SBAS
        QZSS_freq = [1575.42, 1227.60, 1176.45, 1278.75]  # L1, L2, L5, L6 for QZSS
        BDS_freq = [1561.098, 1207.14, 1268.52]  # B1, B2, B3 for BeiDou
        IRNSS_freq = [1176.45, 2492.028]  # L5, S for IRNSS

        return GPS_freq, GLONASS_freq, Galileo_freq, SBAS_freq, QZSS_freq, BDS_freq, IRNSS_freq

    def sort_sat_types(self, obs):
        '''
        Function that creates arrays for each satellite constellation type
        and then assigns satellites from the RINEX file to the corresponding arrays

        Inputs:
        obs: the loaded observation file

        Outputs:
        7 lists of satellite numbers that were included in the loaded obs file
        G-GPS, R-GLONASS, S-SBAS, E-Galileo, C-BeiDou, J-QZSS, I-IRNSS

        '''

        svG = []
        svR = []
        svS = []
        svE = []
        svC = []
        svJ = []
        svI = []

        # Assign values to distinguish
        for i in range(0, len(obs.sv)):
            if str(obs.sv[i].values)[0] == 'G':
                svG.append(str(obs.sv[i].values))
            elif str(obs.sv[i].values)[0] == 'R':
                svR.append(str(obs.sv[i].values))
            elif str(obs.sv[i].values)[0] == 'S':
                svS.append(str(obs.sv[i].values))
            elif str(obs.sv[i].values)[0] == 'E':
                svE.append(str(obs.sv[i].values))
            elif str(obs.sv[i].values)[0] == 'C':
                svC.append(str(obs.sv[i].values))
            elif str(obs.sv[i].values)[0] == 'J':
                svJ.append(str(obs.sv[i].values))
            elif str(obs.sv[i].values)[0] == 'I':
                svI.append(str(obs.sv[i].values))

        return svG, svR, svS, svE, svC, svJ, svI

    def append_MP_arrays(self, sv, obs, freq):

        '''
        Function for appending results from 'calculate_MP1_G' function.
        Used when we want MP results for many satellites in one variable.

        Inputs:
        sv: a list of GPS satellites that appear in the loaded RINEX file
        obs: the loaded observation file
        freq: list of frequencies for GPS

        Output:
        An array of appended MP values for many satellites which makes it easy to store/plot results

        '''

        MP = []
        for i in range(0, len(sv)):
            MP.append(self.calculate_MP1_G(obs, sv[i], freq))
        return MP

    def calculate_MP1_G(self, obs, sv, freq):

        '''
        Calculate code multipath for first frequency, for GPS

        Inputs:
        obs: the loaded observation file
        sv: chosen satellite vehicle
        freq: list of frequencies for GPS

        Output:
        MP1 value

        '''

        c = 299792458
        f1 = freq[0]*1e6
        f2 = freq[1]*1e6
        obs = obs.sel(sv=sv).dropna(dim='time', how='all')

        P1 = obs.C1C
        L1 = obs.L1C*c/f1
        L2 = obs.L2W*c/f2

        MP = self.MP1(P1, L1, L2, f1, f2)

        # Ambiguities
        navg = np.sum(MP)/len(MP)
        MP = MP - navg

        return MP

    def MP1(self, P, L1, L2, f1, f2):

        '''
        
        Return Multipath equation for the first frequency

        '''
        
        MP = P - (f1**2 + f2**2)/(f1**2 - f2**2)*L1 + (2*(f2**2))/(f1**2 - f2**2)*L2
        return MP

    def plot_all_MP(self, sv, MP):

        '''
        Use the result of 'append_MP_arrays' function to plot multipath.

        Inputs:
        sv: a list of GPS satellites that appear in the loaded RINEX file
        MP: output from the 'append_MP_arrays' function

        Output:
        No output, the function uses the command 'ax.plot()' as many times as there are satellites

        '''

        ax = figure(figsize=(10, 6)).gca()
        for i in range(0, len(sv)):
            ax.plot(MP[i], label=sv[i])
        ax.grid()
        ax.legend()
        plt.xlabel('Time')
        plt.ylabel('MP1 [meters]')
        show()


# %%
# Testing code (comment out rinex file loading (line 151) after first time, as it takes a long time to load)
obs = gr.load(
    '../tests/test_data/Rinex3/KLSQ00GRL_R_20213070000_01D_15S_MO.rnx',
    tlim=['2021-11-03T11:32', '2021-11-03T12:32'])
mptest = Multipath(obs)
MP = mptest.get_MP()
