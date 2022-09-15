#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lars Stenseng.

@mail: lars@stenseng.net
"""
# line 853 testing functions
# from qc.__version__ import __version__

import numpy as np
from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
from qc.helper_functions import helper_functions as hf

class Slips:
    """
    Class for cycle slip detection of RINEX files.

    Provides options for different detection algorithms.

    Parameters
    ----------
    L1 : TYPE
        DESCRIPTION.

    Returns
    -------
    L4 : TYPE
        DESCRIPTION.

    """

    def __init__(self, obs, sv, hdr, constellation, codes, rnx_version):

        self.obs = obs
        self.sv = sv
        self.hdr = hdr
        self.constellation = str(constellation)
        self.codes = codes
        self.rnx_version = rnx_version

    def get_slips(self) -> list:

        # self.slips, self.legend = self.__append_slip_arrays(self.obs)
        # self.plot_slips(self.slips, self.legend)

        self.slips = self.slips_MW_single_freq(self.obs, self.sv)

        return self.slips

    '''
    def __append_slip_arrays(self, obs):

        if not self.constellation == 'M':
            sv = self.__sort_sat_types(self.obs)
            # sv = hf().sort_sat_types(self.obs, self.constellation)
            slips = []
            sv_legend = []
            for i in range(0, len(sv)):
                L4_diff, slips_nr, criterion, out_codes = self.slips_MW_single_freq(obs, sv[i])
                if L4_diff is None:
                    continue
                else:
                    # MP.append(xr.DataArray(MP1[:], coords={
                    #    "sv": MP1[:].sv.values, "time": MP1[:].time.values,
                    #    "obs_codes": str(out_codes)}))
                    slips.append(xr.DataArray(L4_diff[:], coords={
                        "sv": L4_diff[:].sv.values, "obs_codes": str(out_codes),
                        "slips_nr": slips_nr, "criterion": criterion}))
                    sv_legend.append(sv[i])
            return slips, sv_legend
        else:
            slips_all = [[], [], [], []]
            sv_legend_all = [[], [], [], []]
            sv_all = self.__sort_sat_types(self.obs)
            # sv_all = hf().sort_sat_types(self.obs, self.constellation)
            for i in range(0, 4):  # Four arrays (svG, svR, svE, svC)
                for j in range(0, len(sv_all[i])):
                    L4_diff, slips_nr, criterion, out_codes = self.slips_MW_single_freq(obs, sv_all[i][j])
                    if L4_diff is None:
                        continue
                    else:
                        slips_all[i].append(xr.DataArray(L4_diff[:], coords={
                            "sv": L4_diff[:].sv.values, "obs_codes": str(out_codes),
                            "slips_nr": slips_nr, "criterion": criterion}))
                        # sv_legend.append(sv[i])
                        # MP_all[i].append(xr.DataArray(MP1[:], coords={
                        #    "sv": MP1[:].sv.values, "time": MP1[:].time.values,
                        #    "obs_codes": str(out_codes)}))
                        sv_legend_all[i].append(sv_all[i][j])
            return slips_all, sv_legend_all
        '''

    def slips_MW_single_freq(self, obs, sv):
        """
        Cycle slip detection algorithm 1.

        Based on Melbourne-Wuebbena,
        but only on carrier phase data (geometry-free combination)
        (from Vaclavovic-Dousa 2016 article)

        Parameters
        ----------
        obs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        '''
        # Select a list of GPS satellites
        svG = []
        for i in range(0, len(obs.sv)):
            if str(obs.sv[i].values)[0] == 'G':
                svG.append(str(obs.sv[i].values))
            else:
                continue
        '''
        # Melbourne-Wuebbena parameters (predetermined)
        I_max = 0.4  # Maximal ionospheric delay [m/h]
        k = 4  # criterion factor

        # For each tracked satellite
        # for i in range(0, len(sv)):

        # Determine rinex version
        if self.rnx_version == 3:
            const_def, freq = hf().get_const_data_vars(
                                            sv,
                                            self.constellation)
        else:
            const_def, freq = hf().get_const_data_vars_rnx2(
                                            sv,
                                            self.constellation)

        # Set obs based on sat
        obs = obs.sel(sv=sv).dropna(dim='time', how='all')

        # Select default obs codes (L1 and L2)
        if self.codes is None:
            obs_codes = self.__select_default_observables(obs, const_def)
            out_codes = const_def
        else:
            obs_codes = self.__select_default_observables(obs, self.codes)
            out_codes = self.codes

        if obs_codes is None:
            return None, None
        else:
            L1 = obs[obs_codes[3]]
            L2 = obs[obs_codes[4]]

        L4 = np.abs(L1 - L2)

        sigma_L4 = np.std(L4)

        criterion = k*sigma_L4 + I_max
        slips_nr = 0

        L4_diff = []
        for j in range(1, len(L4)):
            L4_diff.append(np.abs(L4[j] - L4[j-1]))
            if (np.abs(L4[j] - L4[j-1]) > criterion):
                # If satisfied, raise cycle-slip flag
                slips_nr = slips_nr + 1

        ax = figure(figsize=(10, 6)).gca()
        ax.plot(L2.time[1:], L4_diff, label=sv, linewidth=1.0)
        # labelfull = 'Slip limit: ', criterion.values
        plt.axhline(y=criterion, label='Slip limit', linestyle='-', color='r')
        ax.grid()
        ax.legend()
        plt.xlabel('Time [epochs]')
        plt.ylabel('L4')
        show()

        print('Sat:', sv,
              ', Slips:', slips_nr,
              ', Slip criterion:',  criterion.values)

        return L4_diff, slips_nr, criterion, out_codes

    '''
    def plot_slips(self, slips, legend):
        """
        Plot cycle slips for one satellite vehicle.

        Parameters
        ----------
        obs : TYPE
            DESCRIPTION.
        sat_nr : str
            DESCRIPTION.

        Returns
        -------
        None.

        """
        sat = obs.sel(sv=sat_nr).dropna(dim='time', how='all')

        I_max = 0.4  # Maximal ionospheric delay [m/h]
        k = 4  # criterion factor

        L1 = sat['L1C']
        L2 = sat['L2W']
        L4 = np.abs(L1 - L2)

        sigma_L4 = np.std(L4)

        criterion = k*sigma_L4 + I_max
        slips_nr = 0

        L4_diff = []
        for i in range(1, len(L4)):
            L4_diff.append(np.abs(L4[i] - L4[i-1]))
            if (np.abs(L4[i] - L4[i-1]) > criterion):
                # If satisfied, raise cycle-slip flag
                slips_nr = slips_nr + 1

        ax = figure(figsize=(10, 6)).gca()
        ax.plot(L2.time[1:], L4_diff, label=sat_nr, linewidth=1.0)
        # labelfull = 'Slip limit: ', criterion.values
        plt.axhline(y=criterion, label='Slip limit', linestyle='-', color='r')
        ax.grid()
        ax.legend()
        plt.xlabel('Time [epochs]')
        plt.ylabel('L4')
        show()

        print('Sat:', sat_nr,
              ', Slips:', slips_nr,
              ', Slip criterion:',  criterion.values)
    '''

    def __sort_sat_types(self, obs):
        """
        Create lists for each satellite constellation type.

        Using the chosen satellite constellation (G/R/E/C):
        Assign satellites from the RINEX file to a corresponding list.
        If using Mixed option, it creates four separate lists instead,
        to hold 4 types of constellations. In that case, if a satellite does
        not match any type, it is skipped and not included in any list.

        Parameters
        ----------
        obs : xarray Dataset
            The loaded observation file.

        Returns
        -------
        sv: list
            1 or 4 lists of satellite numbers that were included in obs file.
            G-GPS, R-GLONASS, E-Galileo, C-BeiDou.
            Returns 4 list when using the mixed ('M') option.

        """
        # Case: G (GPS)
        if self.constellation == 'G':
            svG = []
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'G':
                    svG.append(str(obs.sv[i].values))
                else:
                    continue
            return svG
        # Case: R (GLONASS)
        elif self.constellation == 'R':
            svR = []
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'R':
                    svR.append(str(obs.sv[i].values))
                else:
                    continue
            return svR
        # Case: E (Galileo)
        elif self.constellation == 'E':
            svE = []
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'E':
                    svE.append(str(obs.sv[i].values))
                else:
                    continue
            return svE
        # Case: C (BeiDou)
        elif self.constellation == 'C':
            svC = []
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'C':
                    svC.append(str(obs.sv[i].values))
                else:
                    continue
            return svC
        # Case: M (Mixed)
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

    def __select_default_observables(self, obs, const_def):
        """
        Set default values for the P1, P2, P5, L1 and L2 observables.

        If default values are not available (they are NaN or not included),
        use the first available element.
        If P1/P2/P5, L1, or L2 cannot be set at all,
        don't calculate it and skip sat

        Parameters
        ----------
        obs: the loaded observation file
        const_def: default obs types for the given constellation

        Returns
        -------
        P1: L1 pseudorange
        P2: L2 pseudorange
        P5: L5 pseudorange
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
        # if self.MP_eq == 1:
        #    if P1_all == [] or L1_all == [] or L2_all == []:
        #        return None
        # elif self.MP_eq == 2:
        #    if P2_all == [] or L1_all == [] or L2_all == []:
        #        return None
        # elif self.MP_eq == 5:
        #    if P5_all == [] or L1_all == [] or L2_all == []:
        #        return None

        # Specify default P1, P2, P5, L1 and L2
        if const_def[0] in P1_all:
            P1 = obs[const_def[0]]  # If C1C is available, use it
        else:
            P1 = obs[P1_all[0]]  # Else select first element of pseudorange set

        if const_def[1] in P2_all:
            P2 = obs[const_def[1]]  # If C1C is available, use it
        else:
            P2 = obs[P1_all[0]]  # Else select first element of pseudorange set

        if self.rnx_version == 3 or (
                self.rnx_version == 2 and not self.constellation == 'R'):
            if const_def[2] in P5_all:
                P5 = obs[const_def[2]]  # If C1C is available, use it
            else:
                P5 = obs[P1_all[0]]  # Else select first element of pseudorange
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

    def __get_GLONASS_freq_slot(self, sv):
        """
        Get Glonass frequency slots (k) depending on the satellite.

        RINEX 3 version.

        Parameters
        ----------
            sv: chosen satellite vehicle
        Returns
        -------
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

    def __get_GLONASS_freq_slot_rnx2(self, sv):
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

        Parameters
        ----------
            sv: chosen satellite vehicle
        Returns
        -------
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

    def __get_const_data_vars(self, sv):
        """
        Get data variables depending on the chosen sat constellation.

        RINEX 3 version.

        P1, P2, P5, L1 and L2
        (G/R/E/C etc. or M for mixed)

        Parameters
        ----------
            sv: chosen satellite vehicle
        Returns
        -------
            const_def: default constellation observation types
            freq: frequencies specific

        """
        const_def = []
        freq = []
        if self.constellation == 'G':
            const_def = ['C1C', 'C2C', 'C5I', 'L1C', 'L2C']
            freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS
        elif self.constellation == 'R':
            const_def = ['C1C', 'C2C', 'C3I', 'L1C', 'L2C']
            # k: frequency slot for Glonass, used if self.constellation = 'R'
            k = self.__get_GLONASS_freq_slot(sv)
            f1 = (1602 + (k*9/16))
            f2 = (1246 + (k*7/16))
            f3 = 1202.025
            freq = [f1, f2, f3]
        elif self.constellation == 'E':
            const_def = ['C1A', 'C8I', 'C6A', 'L1C', 'L8I']
            freq = [1575.42, 1191.795, 1278.75]  # E1,E5,E6 for Galileo(no L2)
        elif self.constellation == 'C':
            const_def = ['C2I', 'C7I', 'C6I', 'L2I', 'L7I']
            freq = [1561.098, 1207.14, 1268.52]  # B1, B2, B3 for BeiDou
        else:
            if sv[0] == 'G':
                const_def = ['C1C', 'C2C', 'C5I', 'L1C', 'L2C']
                freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS
            elif sv[0] == 'R':
                const_def = ['C1C', 'C2C', 'C3I', 'L1C', 'L2C']
                # k: frequency slot for Glonass, used if self.constel = 'R'
                k = self.__get_GLONASS_freq_slot(sv)
                f1 = (1602 + (k*9/16))
                f2 = (1246 + (k*7/16))
                f3 = 1202.025
                freq = [f1, f2, f3]
            elif sv[0] == 'E':
                const_def = ['C1A', 'C8I', 'C6A', 'L1C', 'L8I']
                freq = [1575.42, 1191.795, 1278.75]  # E1,E5,E6 for Galileo
            elif sv[0] == 'C':
                const_def = ['C2I', 'C7I', 'C6I', 'L2I', 'L7I']
                freq = [1561.098, 1207.14, 1268.52]  # B1, B2, B3 for BeiDou

        return const_def, freq

    def __get_const_data_vars_rnx2(self, sv):
        """
        Get data variables depending on the chosen sat constellation.

        RINEX 2 version.

        P1, P2, P5, L1 and L2
        (G/R/E/C etc. or M for mixed)

        Parameters
        ----------
            sv: chosen satellite vehicle
        Returns
        -------
            const_def: default constellation observation types
            freq: frequencies specific

        """
        const_def = []
        freq = []
        if self.constellation == 'G':
            const_def = ['C1', 'C2', 'C5', 'L1', 'L2']
            freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS
        elif self.constellation == 'R':
            # C3 is not valid but is here to keep array size same as others
            const_def = ['C1', 'C2', 'C3', 'L1', 'L2']
            # k: frequency slot for Glonass, used if self.constellation = 'R'
            k = self.__get_GLONASS_freq_slot_rnx2(sv)
            f1 = (1602 + (k*9/16))
            f2 = (1246 + (k*7/16))
            # f3 = 1202.025  # This option does not exist for RINEX 2
            freq = [f1, f2]
        elif self.constellation == 'E':
            const_def = ['C1', 'C8', 'C6', 'L1', 'L8']
            freq = [1575.42, 1191.795, 1278.75]  # E1,E5,E6 for Galileo(no L2)
        else:
            if sv[0] == 'G':
                const_def = ['C1', 'C2', 'C5', 'L1', 'L2']
                freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS
            elif sv[0] == 'R':
                # C3 is not valid but is here to keep array size same as others
                const_def = ['C1', 'C2', 'C3', 'L1', 'L2']
                # k: frequency slot for Glonass, used if self.constel = 'R'
                k = self.__get_GLONASS_freq_slot_rnx2(sv)
                f1 = (1602 + (k*9/16))
                f2 = (1246 + (k*7/16))
                # f3 = 1202.025  # This option does not exist for RINEX 2
                freq = [f1, f2]
            elif sv[0] == 'E':
                const_def = ['C1', 'C8', 'C6', 'L1', 'L8']
                freq = [1575.42, 1191.795, 1278.75]  # E1,E5,E6 for Galileo

        return const_def, freq
