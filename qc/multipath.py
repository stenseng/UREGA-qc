#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lars Stenseng.

@mail: lars@stenseng.net
"""

import numpy as np
from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
# from qc.__version__ import __version__


class Multipath:
    """The class containing functions for computing various MP equations."""

    def __init__(self, obs, const: str, hdr):  # , signals):
        self.obs = obs
        self.const = const  # Constellation type (G/R/E/C). Mixed coming soon
        self.hdr = hdr
        # self.signals = signals

    def get_MP(self):
        """Master function that uses all other functions to return MP."""
        # freq = self.load_constellation_freq()
        sv_all = self.sort_sat_types(self.obs)

        self.MP, sv_legend = self.append_MP_arrays(sv_all,
                                                   self.obs)
        # Comment out the GPS or GLONASS lines depending on which is used
        # self.MP_G, sv_legend =
        # self.append_MP_arrays(sv_all[0], self.obs, freq[0])  # GPS
        # self.MP_R, sv_legend =
        # self.append_MP_arrays(sv_all[1], self.obs, freq[1])  # GLONASS

        self.plot_all_MP(sv_legend, self.MP)
        # self.plot_all_MP(sv_legend, self.MP_G)  # GPS
        # self.plot_all_MP(sv_legend, self.MP_R)  # GLONASS

        return self.MP
        # return self.MP_G  # GPS
        # return self.MP_R  # GLONASS

    def sort_sat_types(self, obs):
        """
        Create arrays for each satellite constellation type.

        Using the chosen satellite constellation (G/R/E/C):
        Assign satellites from the RINEX file to a corresponding array.
        Should also contain an option for a mixed choice. (M)

        Inputs:
        obs: the loaded observation file

        Outputs:
        7 lists of satellite numbers that were included in the loaded obs file
        G-GPS, R-GLONASS, E-Galileo, C-BeiDou, S-SBAS, J-QZSS, I-IRNSS

        """
        # Case: G (GPS)
        if self.const == 'G':
            svG = []
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'G':
                    svG.append(str(obs.sv[i].values))
                else:
                    continue
            return svG
        # Case: R (GLONASS)
        elif self.const == 'R':
            svR = []
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'R':
                    svR.append(str(obs.sv[i].values))
                else:
                    continue
            return svR
        # Case: E (Galileo)
        elif self.const == 'E':
            svE = []
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'E':
                    svE.append(str(obs.sv[i].values))
                else:
                    continue
            return svE
        # Case: C (BeiDou)
        elif self.const == 'C':
            svC = []
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'C':
                    svC.append(str(obs.sv[i].values))
                else:
                    continue
            return svC
        # Case: M (Mixed). Default case, used even if case was not specified
        else:
            svG = []
            svR = []
            svE = []
            svC = []

            # Assign values to distinguish
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'G':
                    svG.append(str(obs.sv[i].values))
                elif str(obs.sv[i].values)[0] == 'R':
                    svR.append(str(obs.sv[i].values))
                elif str(obs.sv[i].values)[0] == 'E':
                    svE.append(str(obs.sv[i].values))
                elif str(obs.sv[i].values)[0] == 'C':
                    svC.append(str(obs.sv[i].values))
                else:
                    continue  # Skip satellites that don't match any category
            return svG, svR, svE, svC

    def append_MP_arrays(self, sv, obs):
        """
        Append results from 'calculate_MP' function.

        Used when we want MP results for many satellites in one variable.

        Inputs:
        sv: a list of GPS satellites that appear in the loaded RINEX file
        obs: the loaded observation file

        Output:
        An array of appended MP values for many satellites
        which makes it easy to store/plot results

        """
        MP = []
        sv_legend = []
        for i in range(0, len(sv)):
            MP1 = self.calculate_MP(obs, sv[i])
            if MP1 is None:
                continue
            else:
                MP.append(MP1)
                sv_legend.append(sv[i])
        return MP, sv_legend

    def calculate_MP(self, obs, sv):
        """
        Calculate code multipath for first frequency, for GPS.

        Inputs:
        obs: the loaded observation file
        sv: chosen satellite vehicle
        const: constellation type for getting proper obs types (G/R/E/C)

        Output:
        MP1 value
        """
        c = 299792458  # Speed of light

        # Get default obs types and frequencies for chosen constellation
        const_def, freq = self.get_const_data_vars(sv)

        f1 = freq[0]*1e6
        f2 = freq[1]*1e6
        # f3 = freq[2]*1e6
        obs = obs.sel(sv=sv).dropna(dim='time', how='all')

        obs_codes = self.select_default_observables_MP1(obs, const_def)
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
        """Return Multipath equation for the first frequency."""
        MP = P - (f1**2 + f2**2)/(f1**2 - f2**2)*L1 + \
            (2*(f2**2))/(f1**2 - f2**2)*L2
        return MP

    def select_default_observables_MP1(self, obs, const_def):
        """
        Set default values for the P1, L1 and L2 observables.

        The default values are C1C, L1C, L2C.
        If these are not available (they are NaN or not included),
        use the first available element.
        If P1, L1, or L2 cannot be set at all, don't calculate it and skip sat

        Input:
        obs: the loaded observation file
        const_def: default obs types for the given constellation

        Outputs:
        P1: L1 pseudorange
        L1: L1 carrier phase
        L2: L2 carrier phase

        """
        # Find non-nan data variables
        obsNNan_names = [var for var in obs.data_vars
                         if not np.isnan(obs.data_vars[var].values).all()]

        # Get all observable types in separate sets
        P1_all = [var for var in obsNNan_names
                  if var[0:2] == const_def[0][0:2]]
        L1_all = [var for var in obsNNan_names
                  if var[0:2] == const_def[1][0:2]]
        L2_all = [var for var in obsNNan_names
                  if var[0:2] == const_def[2][0:2]]

        # Make a check if the '_all' sets above are empty.
        # If even 1 of them is empty, skip satellite, don't calculate MP for it
        if P1_all == [] or L1_all == [] or L2_all == []:
            return None

        # Specify default P1, L1 and L2
        # (will also need to add condition if the '_all' sets above are empty)
        if 'C1C' in P1_all:
            P1 = obs[const_def[0]]  # If C1C is available, use it
        else:
            P1 = obs[P1_all[0]]  # Else select first element of pseudorange set

        if 'L1C' in L1_all:
            L1 = obs[const_def[1]]  # If L1C is available, use it
        else:
            L1 = obs[L1_all[0]]  # Else select first element of L1 set

        if 'L2C' in L2_all:
            L2 = obs[const_def[2]]  # If L2C is available, use it
        else:
            L2 = obs[L2_all[0]]  # Else select first element of L2 set

        return [P1, L1, L2]

    def get_GLONASS_freq_slot(self, sv, hdr):
        """
        Get Glonass frequency slots (k) depending on the satellite.

        Input:
            sv: chosen satellite vehicle
            hdr: Rinex 3 header file
        Output:
            slot_nr: frequency slot (k value) for Glonass frequencies

        """
        # Find the glonass slot properties in the header file,
        # and split the string to make a list. Start reading from idx1
        slot = hdr['GLONASS SLOT / FRQ #'].split()[1:]
        # Separate the data list into satellites ('R01', 'R02', 'R03'...etc)
        # and slot frequencies (k values)
        glonass_sv = slot[::2]
        glonass_slot = slot[1::2]
        # Predetermine return variable
        slot_nr = 0

        for i in range(0, len(glonass_sv)):
            if glonass_sv[i] == sv:
                # Set the return variable to the match, converting it to int
                slot_nr = int(glonass_slot[i])
                return slot_nr

    def get_const_data_vars(self, sv):
        """
        Get data variables depending on the chosen sat constellation.

        (G/R/E/C etc. or M for mixed)

        Input:
            sv: chosen satellite vehicle
        Output:
            const_def: default constellation observation types
            freq: frequencies specific

        """
        const_def = []
        freq = []
        if self.const == 'G':
            const_def = ['C1C', 'L1C', 'L2C']
            freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS
        elif self.const == 'R':
            const_def = ['C1C', 'L1C', 'L2C']
            # k: frequency slot for Glonass, used if self.const = 'R'
            k = self.get_GLONASS_freq_slot(sv, self.hdr)
            f1 = (1602 + (k*9/16))
            f2 = (1246 + (k*7/16))
            f3 = 1202.025
            freq = [f1, f2, f3]
        elif self.const == 'E':
            const_def = ['C1', 'L2', 'L8']
            freq = [1575.42, 1278.75, 1191.795]  # E1,E6,E5 for Galileo(no L2)
        elif self.const == 'C':
            const_def = ['C2', 'L2', 'L7']
            freq = [1561.098, 1207.14, 1268.52]  # B1, B2, B3 for BeiDou
        else:
            # Placeholder: a case for mixed should be added
            const_def = ['C1C', 'L1C', 'L2C']
            freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS

        return const_def, freq

    def plot_all_MP(self, sv, MP):
        """
        Use the result of 'append_MP_arrays' function to plot multipath.

        Inputs:
        sv: a list of GPS satellites that appear in the loaded RINEX file
        MP: output from the 'append_MP_arrays' function

        Output:
        No output, the function uses the command 'ax.plot()',
        as many times as there are satellites

        """
        ax = figure(figsize=(10, 6)).gca()
        for i in range(0, len(sv)):
            ax.plot(MP[i], label=sv[i])
        ax.grid()
        ax.legend()
        plt.xlabel('Time')
        plt.ylabel('MP1 [meters]')
        show()


# Testing code (comment out rinex file loading after first time,
# as it takes a long time to load)
# obs = gr.load(
#    '../tests/test_data/Rinex3/KLSQ00GRL_R_20213070000_01D_15S_MO.rnx',
#    tlim=['2021-11-03T11:32', '2021-11-03T12:32'])
# mptest = Multipath(obs)
# MP = mptest.get_MP()
