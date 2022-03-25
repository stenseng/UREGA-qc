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

        # Comment out the GPS or GLONASS lines (3 in total) depending on which is used
        #self.MP_G, sv_legend = self.append_MP_arrays(sv_all[0], self.obs, freq[0])  # GPS
        self.MP_R, sv_legend = self.append_MP_arrays(sv_all[1], self.obs, freq[1])  # GLONASS

        #self.plot_all_MP(sv_legend, self.MP_G)  # GPS
        self.plot_all_MP(sv_legend, self.MP_R)  # GLONASS

        #return self.MP_G  # GPS
        return self.MP_R  # GLONASS

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
        sv_legend = []
        for i in range(0, len(sv)):
            MP1 = self.calculate_MP1_G(obs, sv[i], freq)
            if MP1 is None:
                continue
            else:
                MP.append(MP1)
                sv_legend.append(sv[i])
        return MP, sv_legend


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

        obs_codes = self.select_default_observables_MP1(obs)
        if obs_codes is None:
            return None
        else:
            P1 = obs_codes[0]
            L1 = obs_codes[1]
            L2 = obs_codes[2]

        L1 = L1*c/f1
        L2 = L2*c/f2

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

    def select_default_observables_MP1(self, obs):

        '''

        Set default values for the P1, L1 and L2 observables.
        The default values are C1C, L1C, L2C. 
        If these are not available (they are NaN or not included), 
        use the first available element.
        If P1, L1, or L2 cannot be set at all, don't calculate it and skip sat

        Input:
        obs: the loaded observation file

        Outputs:
        P1: L1 pseudorange
        L1: L1 carrier phase
        L2: L2 carrier phase

        '''

        # Find non-nan data variables
        obsNNan_names  = [var for var in obs.data_vars if not np.isnan(obs.data_vars[var].values).all()]

        # Get all observable types in separate sets
        P1_all = [var for var in obsNNan_names if var[0:2] == 'C1']
        L1_all = [var for var in obsNNan_names if var[0:2] == 'L1']
        L2_all = [var for var in obsNNan_names if var[0:2] == 'L2']

        # Make a check if the '_all' sets above are empty. 
        # If even 1 of them is empty, skip satellite (don't calculate MP for it)
        if P1_all == [] or L1_all == [] or L2_all == []:
            return None

        # Specify default P1, L1 and L2 (will also need to add a condition if the '_all' sets above are empty)
        if 'C1C' in P1_all:
            P1 = obs.C1C  # If C1C is available, use it
        else:
            P1 = obs[P1_all[0]]  # Otherwise select the first element of the pseudorange set

        if 'L1C' in L1_all:
            L1 = obs.L1C  # If L1C is available, use it
        else:
            L1 = obs[L1_all[0]]  # Otherwise select the first element of the L1 set

        if 'L2C' in L2_all:
            L2 = obs.L2C  # If L2C is available, use it
        else:
            L2 = obs[L2_all[0]]  # Otherwise select the first element of the L2 set

        return [P1, L1, L2]

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
