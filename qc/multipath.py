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

    def __init__(self, obs, hdr, const: str, MP_eq=1, codes=None, rnx_version = 3):
        self.obs = obs
        self.hdr = hdr
        self.const = const  # Constellation type (G/R/E/C). Mixed coming soon
        self.MP_eq = MP_eq  # MP equation for k frequency (array:1/2/5/all)
        self.codes = codes  # Own obs codes (['C1C','C2C','C5I','L1C','L2C'])
        self.rnx_version = rnx_version  # RINEX 3 or RINEX 2 (3 by default)

    def get_MP(self):
        """Master function that uses all other functions to return MP."""
        #  self.info,
        self.MP, sv_legend = self.append_MP_arrays(self.obs)

        self.plot_all_MP(sv_legend, self.MP)

        return self.MP  # , self.info

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

    def append_MP_arrays(self, obs):
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
        sv = self.sort_sat_types(self.obs)
        # info = []
        MP = []
        sv_legend = []
        print(len(sv))
        for i in range(0, len(sv)):
            # , obs_codes
            MP1 = self.calculate_MP(obs, sv[i])
            if MP1 is None:
                continue
            else:
                MP.append(MP1)
                # info.append([sv[i], obs_codes])
                sv_legend.append(sv[i])
        # info,
        return MP, sv_legend

    def calculate_MP(self, obs, sv):
        """
        Calculate code multipath for chosen constellation and MP frequency

        Inputs:
        obs: the loaded observation file
        sv: chosen satellite vehicle
        const: constellation type for getting proper obs types (G/R/E/C)

        Output:
        MP1 value
        """
        c = 299792458  # Speed of light

        # Get default obs types and frequencies for chosen constellation.
        # Different function based on rnx version
        if self.rnx_version == 3:
            const_def, freq = self.get_const_data_vars(sv)  
            f5 = freq[2]*1e6
        else:
            const_def, freq = self.get_const_data_vars_rnx2(sv)
            # GLONASS doesn't have f5 frequency in rinex 2
            if not self.const == 'R':
                f5 = freq[2]*1e6

        f1 = freq[0]*1e6
        f2 = freq[1]*1e6

        obs = obs.sel(sv=sv).dropna(dim='time', how='all')

        if self.codes is None:
            obs_codes = self.select_default_observables_MP1(obs, const_def)
        #  out_codes = const_def
        else:
            obs_codes = self.select_default_observables_MP1(obs, self.codes)
        #    out_codes = self.codes

        if obs_codes is None:
            return None
        else:
            P1 = obs_codes[0]
            P2 = obs_codes[1]
            if not (self.rnx_version == 2 and self.const == 'R'):
                P5 = obs_codes[2]
            L1 = obs_codes[3]
            L2 = obs_codes[4]

        L1 = L1*c/f1
        L2 = L2*c/f2
        if self.MP_eq == 1:
            MP = self.MP1(P1, L1, L2, f1, f2)
        elif self.MP_eq == 2:
            MP = self.MP2(P2, L1, L2, f1, f2)
        elif self.MP_eq == 5:
            MP = self.MP5(P5, L1, L2, f1, f2, f5)
        else:
            MP = self.MP1(P1, L1, L2, f1, f2)

        # Ambiguities
        navg = np.sum(MP)/len(MP)
        MP = MP - navg

        return MP  # , out_codes

    def MP1(self, P1, L1, L2, f1, f2):
        """Return Multipath equation for the first frequency."""
        MP1 = P1 - L1 - ((2*f2**2)/(f1**2 - f2**2))*(L1-L2)
        return MP1

    def MP2(self, P2, L1, L2, f1, f2):
        """Return Multipath equation for the second frequency."""
        MP2 = P2 - L2 - ((2*f1**2)/(f2**2 - f1**2))*(L2-L1)
        return MP2

    def MP5(self, P5, L1, L2, f1, f2, f5):
        """Return Multipath equation for the frequency band 5."""
        MP5 = P5 - L1 - \
            (f2**2*(f1**2 + f5**2))/(f5**2*(f1**2 - f2**2))*(L1 - L2)
        return MP5

    def select_default_observables_MP1(self, obs, const_def):
        """
        Set default values for the P1, P2, P5, L1 and L2 observables.

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
        if self.rnx_version == '3':
            P1_all = [var for var in obsNNan_names
                      if var[0:2] == const_def[0][0:2]]
            P2_all = [var for var in obsNNan_names
                      if var[0:2] == const_def[1][0:2]]
        else:
        # Rinex2 pseudoranges are described as C1 or P1 etc.
            P1_all = [var for var in obsNNan_names
                      if var[0:2] == const_def[0][0:2] or
                      var[0:2] == 'P1']
            P2_all = [var for var in obsNNan_names
                      if var[0:2] == const_def[1][0:2] or
                      var[0:2] == 'P2']
        P5_all = [var for var in obsNNan_names
                  if var[0:2] == const_def[2][0:2]]
        L1_all = [var for var in obsNNan_names
                  if var[0:2] == const_def[3][0:2]]
        L2_all = [var for var in obsNNan_names
                  if var[0:2] == const_def[4][0:2]]

        # Make a check if the '_all' sets above are empty.
        # If even 1 of them is empty, skip satellite, don't calculate MP for it
        if self.MP_eq == 1:
            if P1_all == [] or L1_all == [] or L2_all == []:
                return None
        elif self.MP_eq == 2:
            if P2_all == [] or L1_all == [] or L2_all == []:
                return None
        elif self.MP_eq == 5:
            if P5_all == [] or L1_all == [] or L2_all == []:
                return None

        # Specify default P1, P2, P5, L1 and L2
        # (will also need to add condition if the '_all' sets above are empty)
        if const_def[0] in P1_all:
            P1 = obs[const_def[0]]  # If C1C is available, use it
        else:
            P1 = obs[P1_all[0]]  # Else select first element of pseudorange set

        if const_def[1] in P2_all:
            P2 = obs[const_def[1]]  # If C1C is available, use it
        else:
            P2 = obs[P1_all[0]]  # Else select first element of pseudorange set

        if self.rnx_version == 3 or (self.rnx_version == 2 and not self.const == 'R'):
            if const_def[2] in P5_all:
                P5 = obs[const_def[2]]  # If C1C is available, use it
            else:
                P5 = obs[P1_all[0]]  # Else select first element of pseudorange set
        else:
        # Can't return without P5 - so set it to P1
            P5 = obs[P1_all[0]]  # (but it will never be used anywhere)

        if const_def[3] in L1_all:
            L1 = obs[const_def[3]]  # If L1C is available, use it
        else:
            L1 = obs[L1_all[0]]  # Else select first element of L1 set

        if const_def[4] in L2_all:
            L2 = obs[const_def[4]]  # If L2C is available, use it
        else:
            L2 = obs[L2_all[0]]  # Else select first element of L2 set

        return [P1, P2, P5, L1, L2]

    def get_GLONASS_freq_slot(self, sv):
        """
        Get Glonass frequency slots (k) depending on the satellite.

        Input:
            sv: chosen satellite vehicle
        Output:
            slot_nr: frequency slot (k value) for Glonass frequencies

        """
        # Find the glonass slot properties in the header file,
        # and split the string to make a list. Start reading from idx1
        slot = self.hdr['GLONASS SLOT / FRQ #'].split()[1:]
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

    def get_GLONASS_freq_slot_rnx2(self, sv):
        """
        Get Glonass frequency slots (k) depending on the satellite.

        RINEX 2 version: header not required, as rnx2 files do not have
        frequency slots in their headers. Manual input is needed.
        This function might also be used when header file in RINEX 3 is not
        available for some reason.
        The frequency slots might change in the future as GLONASS constallation
        development changes in time.
        Note: sat R11(-5) is currently not operational and will not be
        included.

        Input:
            sv: chosen satellite vehicle
        Output:
            slot_nr: frequency slot (k value) for Glonass frequencies
        """

        # Pre-determined frequency slots and their corresponding satellites
        slots = np.array([['R01', 'R02', 'R03', 'R04', 'R05', 'R06', 'R07',
                           'R08', 'R09', 'R10', 'R12', 'R13', 'R14',
                           'R15', 'R16', 'R17', 'R18', 'R19', 'R20', 'R21',
                           'R22', 'R23', 'R24'],
                          [1, -4, 5, 6, 1, -4, 5, 6, -2, -7, -1, -2, -7, 0,
                           -1, 4, -3, 3, 2, 4,  -3, 3, 2]])

        for i in range(0, len(slots[0])):
            if slots[0][i] == sv:
                # Set the return variable to the match, converting it to int
                slot_nr = int(slots[1][i])
                return slot_nr

    def get_const_data_vars(self, sv):
        """
        Get data variables depending on the chosen sat constellation.

        RINEX 3 version.

        P1, P2, P5, L1 and L2
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
            const_def = ['C1C', 'C2C', 'C5I', 'L1C', 'L2C']
            freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS
        elif self.const == 'R':
            const_def = ['C1C', 'C2C', 'C3I', 'L1C', 'L2C']
            # k: frequency slot for Glonass, used if self.const = 'R'
            k = self.get_GLONASS_freq_slot(sv)
            f1 = (1602 + (k*9/16))
            f2 = (1246 + (k*7/16))
            f3 = 1202.025
            freq = [f1, f2, f3]
        elif self.const == 'E':
            const_def = ['C1A', 'C8I', 'C6A', 'L1C', 'L8I']
            freq = [1575.42, 1191.795, 1278.75]  # E1,E5,E6 for Galileo(no L2)
        elif self.const == 'C':
            const_def = ['C2I', 'C7I', 'C6I', 'L2I', 'L7I']
            freq = [1561.098, 1207.14, 1268.52]  # B1, B2, B3 for BeiDou
        else:
            # Placeholder: a case for mixed should be added
            const_def = ['C1C', 'C2C', 'C5I', 'L1C', 'L2C']
            freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS

        return const_def, freq

    def get_const_data_vars_rnx2(self, sv):
        """
        Get data variables depending on the chosen sat constellation.

        RINEX 2 version.

        P1, P2, P5, L1 and L2
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
            const_def = ['C1', 'C2', 'C5', 'L1', 'L2']
            freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS
        elif self.const == 'R':
            # C3 is not valid but is here to keep array size same as others
            const_def = ['C1', 'C2', 'C3', 'L1', 'L2']
            # k: frequency slot for Glonass, used if self.const = 'R'
            k = self.get_GLONASS_freq_slot_rnx2(sv)
            f1 = (1602 + (k*9/16))
            f2 = (1246 + (k*7/16))
            # f3 = 1202.025  # This option does not exist for RINEX 2
            freq = [f1, f2]
        elif self.const == 'E':
            const_def = ['C1', 'C8', 'C6', 'L1', 'L8']
            freq = [1575.42, 1191.795, 1278.75]  # E1,E5,E6 for Galileo(no L2)
        # This option does not exist for RINEX 2
        # elif self.const == 'C':
        #    const_def = ['C2I', 'C7I', 'C6I', 'L2I', 'L7I']
        #    freq = [1561.098, 1207.14, 1268.52]  # B1, B2, B3 for BeiDou
        # if self.const == 'C':
        #    throw exception ? (invalid input)
        else:
            # Placeholder: a case for mixed should be added
            const_def = ['C1', 'C2', 'C5', 'L1', 'L2']
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
        plt.ylabel('MP' + str(self.MP_eq) + ' [meters]')
        show()
