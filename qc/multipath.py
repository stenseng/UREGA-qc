#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multipath class.

    Using RINEX version 2 or 3 files, calculates multipath equations
    with frequency band 1, 2 or 5 and plots results against time.
    Also allows for specifying constellation type and observation codes.

Author: Magdalena Golofit
02/05/2022
"""

import numpy as np
from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import matplotlib.dates as mdates
from qc.helper_functions import helper_functions as hf

class Multipath:
    """
    The class containing functions for computing MP equations.

    Call Multipath object with the desired settings and then use get_MP()
    function to get results.
    Can return MP1, MP2 or MP5, using available constellations: G(GPS),
    R(GLONASS), E(Galileo), C(BDS) or M(Mixed-return all available ones
    at once).
    Allows using RINEX 2 or 3 versions.
    Allows using chosen observation codes (not possible when using Mixed).

    Parameters
    ----------
    obs : xarray Dataset
        The loaded observation file.
    hdr : dict
        Header file
    constellation: string
        Constellation type G/R/E/C/M
    MP_eq: int
        Multipath equation 1/2/5
    codes: list
        Observation codes (optional)
    rnx_version: int
        Specified RINEX version

    """

    def __init__(self,
                 obs,
                 hdr,
                 constellation: str,  # Constellation type (G/R/E/C/M)
                 MP_eq: int = 1,  # MP equation for k frequency band (k=1/2/5)
                 codes: list = None,  # Own obs codes
                 rnx_version: int = 3):  # RINEX 3 or RINEX 2 (3 by default)

        self.obs = obs
        self.hdr = hdr

        if not isinstance(constellation, str):
            raise TypeError('Constellation parameter must be a string.')
        if constellation not in ['G', 'R', 'E', 'C', 'M']:
            raise TypeError('Invalid constellation parameter: pick G/R/E/C/M.')
        self.constellation = str(constellation)

        if not isinstance(MP_eq, int):
            raise TypeError('MP equation must be an integer.')
        if MP_eq not in [1, 2, 5]:
            raise TypeError('Invalid MP equation: pick band 1, 2 or 5.')
        self.MP_eq = MP_eq

        if constellation == 'M' and codes is not None:
            raise TypeError(
                'When using Mixed case, must set observation codes to None.')
        self.codes = codes

        if not isinstance(rnx_version, int) and rnx_version not in [2, 3]:
            raise TypeError(
                'RINEX version must be an integer equal to 2 or 3.')
        self.rnx_version = rnx_version

    def get_MP(self) -> list:
        """
        Get multipath delay for given Multipath object and plot it.

        Use get_MP() after calling the Multipath object.
        Master function that uses all other functions of the Multipath class
        to return MP and plot the delay versus time.

        Returns
        -------
        MP: list of xarray DataArrays
            1 or 4 lists of xarray DataArrays containing multipath delays.
            Additionally contains information about satellite vehicle for each
            array, time interval, and used observation codes.
            Returns 4 list when using the mixed ('M') option.

        """
        self.MP, sv_legend = self.__append_MP_arrays(self.obs)

        self.plot_all_MP(sv_legend, self.MP)

        return self.MP

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
        '''

    def __append_MP_arrays(self, obs):
        """
        Append results from calculate MP() function.

        Get results from all satellites in one variable.
        If Mixed option was chosen, append results for four variables
        separately, for each constellation list, so that they
        can be plotted separately.

        Parameters
        ----------
        obs : xarray Dataset
            The loaded observation file.

        Returns
        -------
        MP_all: xarray Dateset
            Xarray with all MP values for all satellites.
        sv_legend: list
            A list of included satellites, used for plotting.

        """
        if not self.constellation == 'M':
            # sv = self.__sort_sat_types(self.obs)
            sv = hf().sort_sat_types(self.obs, self.constellation)
            MP = []
            sv_legend = []
            for i in range(0, len(sv)):
                MP1, out_codes = self.__calculate_MP(obs, sv[i])
                if MP1 is None:
                    continue
                else:
                    MP.append(xr.DataArray(MP1[:], coords={
                        "sv": MP1[:].sv.values, "time": MP1[:].time.values,
                        "obs_codes": str(out_codes)}))
                    sv_legend.append(sv[i])
            return MP, sv_legend
        else:
            MP_all = [[], [], [], []]
            sv_legend_all = [[], [], [], []]
            # sv_all = self.__sort_sat_types(self.obs)
            sv_all = hf().sort_sat_types(self.obs, self.constellation)
            for i in range(0, 4):  # Four arrays (svG, svR, svE, svC)
                for j in range(0, len(sv_all[i])):
                    MP1, out_codes = self.__calculate_MP(obs, sv_all[i][j])
                    if MP1 is None:
                        continue
                    else:
                        MP_all[i].append(xr.DataArray(MP1[:], coords={
                            "sv": MP1[:].sv.values, "time": MP1[:].time.values,
                            "obs_codes": str(out_codes)}))
                        sv_legend_all[i].append(sv_all[i][j])
            return MP_all, sv_legend_all

    def __calculate_MP(self, obs, sv):
        """
        Calculate code multipath for chosen constellation and MP frequency.

        Parameters
        ----------
        obs: xarray Dataset
            The loaded observation file
        sv: string
            Chosen satellite vehicle

        Returns
        -------
        MP: xarray DataArray
            Multipath delay

        """
        c = 299792458  # Speed of light

        # Get default obs types and frequencies for chosen constellation.
        # Different function based on rnx version
        if self.rnx_version == 3:
            const_def, freq = hf().get_const_data_vars(sv, self.constellation)
            f5 = freq[2]*1e6
        else:
            const_def, freq = hf().get_const_data_vars_rnx2(sv, self.constellation)
            # GLONASS doesn't have f5 frequency in rinex 2
            if not sv[0] == 'R':
                f5 = freq[2]*1e6

        f1 = freq[0]*1e6
        f2 = freq[1]*1e6

        # One satellite at a time
        obs = obs.sel(sv=sv).dropna(dim='time', how='all')

        if self.codes is None:
            obs_codes = self.__select_default_observables(obs, const_def)
            out_codes = const_def
        else:
            obs_codes = self.__select_default_observables(obs, self.codes)
            out_codes = self.codes

        if obs_codes is None:
            return None, None
        else:
            P1 = obs_codes[0]
            P2 = obs_codes[1]
            if not (self.rnx_version == 2 and self.constellation == 'R'):
                P5 = obs_codes[2]
            L1 = obs_codes[3]
            L2 = obs_codes[4]

        L1 = L1*c/f1
        L2 = L2*c/f2
        if self.MP_eq == 1:
            MP = self.__MP1(P1, L1, L2, f1, f2)
        elif self.MP_eq == 2:
            MP = self.__MP2(P2, L1, L2, f1, f2)
        elif self.MP_eq == 5:
            MP = self.__MP5(P5, L1, L2, f1, f2, f5)
        else:  # Default MP 1
            MP = self.__MP1(P1, L1, L2, f1, f2)

        # Ambiguities
        navg = np.sum(MP)/len(MP)
        MP = MP - navg

        return MP, out_codes

    def __MP1(self, P1, L1, L2, f1, f2):
        """Return Multipath equation for the first frequency."""
        MP1 = P1 - L1 - ((2*f2**2)/(f1**2 - f2**2))*(L1-L2)
        return MP1

    def __MP2(self, P2, L1, L2, f1, f2):
        """Return Multipath equation for the second frequency."""
        MP2 = P2 - L2 - ((2*f1**2)/(f2**2 - f1**2))*(L2-L1)
        return MP2

    def __MP5(self, P5, L1, L2, f1, f2, f5):
        """Return Multipath equation for the frequency band 5."""
        MP5 = P5 - L1 - \
            (f2**2*(f1**2 + f5**2))/(f5**2*(f1**2 - f2**2))*(L1 - L2)
        return MP5

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
    '''
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
    '''

    def plot_all_MP(self, sv, MP):
        """
        Use the result of 'append_MP_arrays' function to plot multipath.

        If only one specifc constellation was selected, return only one plot.
        If the mixed option for constellations was selected, return 4 plots
        for RINEX 3 version and 3 plots for RINEX 2 version
        (as RINEX 2 does not include BDS constellation).

        Parameters
        ----------
        sv: a list of GPS satellites that appear in the loaded RINEX file
        MP: output from the 'append_MP_arrays' function

        Returns
        -------
        No output, the function uses the command 'ax.plot()',
        as many times as there are satellites

        """
        titles = ['GPS', 'GLONASS', 'Galileo', 'BeiDou']
        const_options = ['G', 'R', 'E', 'C']
        if not self.constellation == 'M':
            ax = figure(figsize=(10, 6)).gca()
            for i in range(0, len(sv)):
                ax.plot(MP[i].time, MP[i], label=sv[i])
            ax.grid()
            ax.legend()
            plt.xlabel('Time [epochs]')
            plt.ylabel('MP' + str(self.MP_eq) + ' [meters]')
            title = titles[const_options.index(self.constellation)]
            plt.title(title)
            show()
        else:
            if self.rnx_version == 2:
                loop_range = 3  # Not plotting BeiDou (doesn't exist)
            else:
                loop_range = 4  # Plot all for rnx 3

            for i in range(0, loop_range):
                ax = figure(figsize=(10, 6)).gca()
                for j in range(0, len(sv[i])):
                    ax.plot(MP[i][j].time, MP[i][j], label=sv[i][j])
                ax.grid()
                ax.legend()
                plt.title(titles[i])
                plt.xlabel('Time [epochs]')
                plt.ylabel('MP' + str(self.MP_eq) + ' [meters]')
                show()
