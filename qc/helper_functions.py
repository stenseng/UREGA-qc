#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Helper functions, used in more than one QC classes.

Used in Multipath and CycleSlips

@author: Magdalena
"""
import numpy as np


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

    def sort_sat_types(self, obs):
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

    def get_const_data_vars(self, sv):
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

    def get_const_data_vars_rnx2(self, sv):
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
