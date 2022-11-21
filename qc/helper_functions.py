#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Helper functions, used in more than one QC classes.

Used in Multipath and CycleSlips

@author: Magdalena
"""
import numpy as np
import xarray as xr

class helper_functions:
    """
    The class containing helper functions for Quality Control.

    Copied directly from the original (synthesis commit) Multipath class.
    The below functions can be used both in Multipath and Cycle Slips classes.
    They were collected in a helper class to avoid repeating functions.

    Parameters
    ----------
        None, provided functions are used in the code.

    """

    def sort_sat_types(self, obs, constellation):
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
        if constellation == 'G':
            svG = []
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'G':
                    svG.append(str(obs.sv[i].values))
                else:
                    continue
            return svG
        # Case: R (GLONASS)
        elif constellation == 'R':
            svR = []
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'R':
                    svR.append(str(obs.sv[i].values))
                else:
                    continue
            return svR
        # Case: E (Galileo)
        elif constellation == 'E':
            svE = []
            for i in range(0, len(obs.sv)):
                if str(obs.sv[i].values)[0] == 'E':
                    svE.append(str(obs.sv[i].values))
                else:
                    continue
            return svE
        # Case: C (BeiDou)
        elif constellation == 'C':
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

    def get_const_data_vars(self, sv, constellation):
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
        if constellation == 'G':
            const_def = ['C1C', 'C2C', 'C5I', 'L1C', 'L2C']
            freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS
        elif constellation == 'R':
            const_def = ['C1C', 'C2C', 'C3I', 'L1C', 'L2C']
            # k: frequency slot for Glonass, used if self.constellation = 'R'
            k = self.__get_GLONASS_freq_slot(sv)
            f1 = (1602 + (k*9/16))
            f2 = (1246 + (k*7/16))
            f3 = 1202.025
            freq = [f1, f2, f3]
        elif constellation == 'E':
            const_def = ['C1A', 'C8I', 'C6A', 'L1C', 'L8I']
            freq = [1575.42, 1191.795, 1278.75]  # E1,E5,E6 for Galileo(no L2)
        elif constellation == 'C':
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

    def get_const_data_vars_rnx2(self, sv, constellation):
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
        if constellation == 'G':
            const_def = ['C1', 'C2', 'C5', 'L1', 'L2']
            freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS
        elif constellation == 'R':
            # C3 is not valid but is here to keep array size same as others
            const_def = ['C1', 'C2', 'C3', 'L1', 'L2']
            # k: frequency slot for Glonass, used if self.constellation = 'R'
            k = self.__get_GLONASS_freq_slot_rnx2(sv)
            f1 = (1602 + (k*9/16))
            f2 = (1246 + (k*7/16))
            # f3 = 1202.025  # This option does not exist for RINEX 2
            freq = [f1, f2]
        elif constellation == 'E':
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

    def get_GLONASS_freq_slot(self, sv):
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

    def data_pre_check(data, const_def):
        """
        Data pre-check function.

        Function that checks if the provided data has adequate size, contains nans,
        or if time differences between measurements are equal in size.
        Gives data check overview (faulty data)
           - 1st row: Satellite going under horizon
           - 2nd row: Measurements not taken at a regular interval
           - 3rd row: Measurements contain NaNs
        Removes nan measurements (if they exist) for correctly computing
        cycle slip algorithms.

        Parameters
        ----------
        data : xarray
            Input rinex file.
        const_def : default obs types for the given constellation

        Returns
        -------
        data2 : xarray
            Xarray clean from nan values or other faulty measurements.
        faulty_idx: array
            Data timestamps showing where bad data is.
            1st row: Satellite going under horizon
            2nd row: Measurements not taken at a regular interval
            3rd row: Measurements contain NaNs

        """
        # Function making basic checks of data quality before cycle slip detection
        out = np.empty(len(data.time))
        out[0] = 0
        # N-consecutive epochs
        # - Long period without measurements (satellite under horizon)
        # - Check for measurement epoch differences: small vs large time difference
        N = 20*data.interval  # (N set to 5 minutes)
        for i in range(1, len(data.time)):  # Moving window method
            # 1st check: differences indicate sat went under horizon
            if int((data.time[i] - data.time[i-1]).values)*1e-9 >= N:
                out[i] = 1  # sat under horizon
            # 2nd check: timestamp differences not a regular interval
            elif round(int((data.time[i] - data.time[i-1]).values)*1e-9, 5) != data.interval:
                out[i] = 2  # faulty measurement intervals
            elif (np.isnan(data.data_vars[const_def[0]][i]) or
                  np.isnan(data.data_vars[const_def[1]][i]) or
                  np.isnan(data.data_vars[const_def[2]][i]) or
                  np.isnan(data.data_vars[const_def[3]][i]) or
                  np.isnan(data.data_vars[const_def[4]][i])):
                out[i] = 3  # nans in measurements
            else:
                out[i] = 0  # ok values

        # take needed variables
        vars_all = [[], [], [], [], []]
        for i in range(0, len(const_def)):
            vars_all[i] = [float(i) for i in data.data_vars[const_def[i]]]

        times = data.time.values

        # Out message, to be displayed after the function is used
        out1 = 0
        out2 = 0
        out3 = 0
        faulty_idx = [[], [], []]
        for i in range(0, len(out)):
            if out[i] == 1:
                out1 = out1 + 1
                faulty_idx[0].append(times[i-1])
            elif out[i] == 2:
                out2 = out2 + 1
                faulty_idx[1].append(times[i])
            elif out[i] == 3:
                out3 = out3 + 1
                faulty_idx[2].append(times[i])
        print('Pre-check complete, LLI:', out1, ', faulty data intervals: ',
              out2, ', NaN measurements:', out3, '. Detecting cycle slips: . . .')

        # Remove NaN measurements by index from [out]

        # this loop is reversed for proper indexing when removing data
        for i in range(len(out)-1, -1, -1):
            if out[i] == 3:
                for j in range(0, 5):
                    del vars_all[j][i]
                times = np.delete(times, i)

        # convert back to xarray format
        data2 = xr.Dataset(
             data_vars=dict(
                 C1=(["time"], vars_all[0]),
                 C2=(["time"], vars_all[1]),
                 C3=(["time"], vars_all[2]),
                 L1=(["time"], vars_all[3]),
                 L2=(["time"], vars_all[4])
                 ),
             coords=dict(
                 sv=str(data.sv.values),
                 time=times),
             attrs=data.attrs)

        return data2, faulty_idx

    def MW_dual_freq(data):
        # define variables
        min_arc = 300  # seconds

        B_w = np.empty(len(data.data_vars['L1']))
        B_w[0] = data.data_vars['L1'][0] - data.data_vars['C1'][0]

        mean_Bw = np.empty(len(data.time))
        mean_Bw[0] = 0

        sigma_Bw = np.empty(len(data.time))
        sigma_Bw[0] = 1

        lambda_w = 0.86  # meter
        k_factor = 5
        slips = []

        # sliding window method
        for i in range(1, len(data.time)):
            arc_length = i*data.interval
            B_w[i] = data.data_vars['L1'][i] - data.data_vars['C1'][i]
            # is arc length long enough? if not, update variables
            if arc_length > min_arc:
                # update difference d
                d = B_w[i] - mean_Bw[i-1]

                # update threshold th
                th = k_factor * np.sqrt(sigma_Bw[i-1])

                # check for cycle slip
                if ((np.abs(d) > th) and (np.sqrt(sigma_Bw[i-1] <= lambda_w))) or (np.isnan(B_w[i])):
                    slips.append(0)  # 0 for slip
                else:
                    slips.append(1)  # 1 for ok values

                mean_Bw[i] = ((i-1)/i)*mean_Bw[i-1] + (1/i)*B_w[i]
                sigma_Bw[i] = ((i-1)/i)*sigma_Bw[i-1] + (1/i)*(B_w[i] - mean_Bw[i-1])**2

            else:
                # update mean, sigma and 300s sliding window mean
                mean_Bw[i] = ((i-1)/i)*mean_Bw[i-1] + (1/i)*B_w[i]
                sigma_Bw[i] = ((i-1)/i)*sigma_Bw[i-1] + (1/i)*(B_w[i] - mean_Bw[i-1])**2

        return slips
