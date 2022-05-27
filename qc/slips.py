#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lars Stenseng.

@mail: lars@stenseng.net
"""

# from qc.__version__ import __version__

import georinex as gr
import numpy as np
from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt

obs = gr.load(
    'tests/test_data/Rinex3/KLSQ00GRL_R_20213070000_01D_15S_MO.rnx',
    # tlim=['2021-11-03T12:00', '2021-11-03T12:30'])
    tlim=['2021-11-03T05:30', '2021-11-03T07:30'])
# tlim=['2021-11-03T15:00', '2021-11-03T18:00'])
# hdr = gr.rinexheader(
#    'tests/test_data/Rinex3/KLSQ00GRL_R_20213070000_01D_15S_MO.rnx')
# rnx_version = 3

# %%  Starting test

# Copying helper functions from Multipath class - later on, it could be turned
# into a separate class with helper functions

# Pick GPS satellites
svG = []
for i in range(0, len(obs.sv)):
    if str(obs.sv[i].values)[0] == 'G':
        svG.append(str(obs.sv[i].values))
    else:
        continue

# %%
# 5:30 to 7:30, G08 and G21 give 2 cycle slips # [290:300]

# 'G01','G06','G08','G10','G12','G14','G17','G19','G21','G22','G24','G30','G32'
sat = 'G21'
sattest = obs.sel(sv=sat).dropna(dim='time', how='all')
# G02 data vars with no-nan: C1C, D1C, L1C, S1C, C1W, C2W, D2W, L2W, S1W, S2W

I_max = 0.4  # Maximal ionospheric delay [m/h]
k = 4  # criterion factor

L1 = sattest['L1C']  # GPS
L2 = sattest['L2W']  # GPS

# L1 = sattest['L1C']  # Galileo
# L2 = sattest['L8Q']  # Galileo
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
ax.plot(L2.time[1:], L4_diff, label=sat)
plt.axhline(y=criterion, label='Slip limit', linestyle='-', color='r')
ax.grid()
ax.legend()
plt.xlabel('Time [epochs]')
plt.ylabel('L4')
plt.title('Single-frequency Melbourne-Wuebbena')
show()

print('Slips:', slips_nr, ', Slip criterion:',  criterion.values)

# %%
# Plot all loaded sats, L1 and L2

ax = figure(figsize=(10, 6)).gca()
for i in range(0, len(svG)):
    test = obs.sel(sv=svG[i]).dropna(dim='time', how='all')
    L1test = test['L1C']
    L2test = test['L2W']
    ax.plot(L1test.time, L1test, label=svG[i], linewidth=2.0)
    #ax.plot(L2test.time, L2test, label='L2', linewidth=0.5)
ax.grid()
ax.legend()
plt.xlabel('Time [epochs]')
plt.ylabel('Carrier phases')
show()

# %%
# Plot separate sats, L1 and L2
ax = figure(figsize=(10, 6)).gca()
test = obs.sel(sv='E21').dropna(dim='time', how='all')
L1test = test['L1C']
L2test = test['L2W']
ax.plot(L1test.time, L1test, label='L1', linewidth=2.0)
ax.plot(L2test.time, L2test, label='L2', linewidth=1.0)
ax.grid()
# ax.legend()
plt.xlabel('Time [epochs]')
plt.ylabel('Carrier phases')
show()

# %% Dual-frequency Melbourne-Wuebbena testing

# 'G01','G06','G08','G10','G12','G14','G17','G19','G21','G22','G24','G30','G32'
sat = 'G21'
sattest = obs.sel(sv=sat).dropna(dim='time', how='all')
# G02 data vars with no-nan: C1C, D1C, L1C, S1C, C1W, C2W, D2W, L2W, S1W, S2W

freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS

f1 = freq[0]*1e6
f2 = freq[1]*1e6

P1 = sattest['C1C']
P2 = sattest['C2W']
L1 = sattest['L1C']  # GPS
L2 = sattest['L2W']  # GPS

# L1 = sattest['L1C']  # Galileo
# L2 = sattest['L8Q']  # Galileo
L6 = (1/(f1-f2))*(f1*L1 - f2*L2) - (1/(f1+f2))*(f1*P1 + f2*P2)

sigma_L6 = np.std(L6)
k = 4  # criterion factor

criterion = k*sigma_L6
slips_nr = 0

L6_diff = []
for i in range(1, len(L6)):
    L6_diff.append(np.abs(L6[i] - L6[i-1]))
    if (np.abs(L6[i] - L6[i-1]) > criterion):
        # If satisfied, raise cycle-slip flag
        slips_nr = slips_nr + 1

ax = figure(figsize=(10, 6)).gca()
ax.plot(L2.time[1:], L6_diff, label=sat)
plt.axhline(y=criterion, label='Slip limit', linestyle='-', color='r')
ax.grid()
ax.legend()
plt.xlabel('Time [epochs]')
plt.ylabel('L6')
plt.title('Dual-frequency Melbourne-Wuebbena')
show()

print('Slips:', slips_nr, ', Slip criterion:',  criterion.values)

# %%  Work in Progress


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

    def __init__(self):
        pass

    def slips_MW_single_freq(self, obs):
        """
        Cycle slip detection algorithm 1.

        Based on Melbourne-Wuebbena,
        but only on carrier phase data (single-frequency)
        (from Vaclavovic-Dousa 2016 article)

        Parameters
        ----------
        obs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        # Select a list of GPS satellites
        svG = []
        for i in range(0, len(obs.sv)):
            if str(obs.sv[i].values)[0] == 'G':
                svG.append(str(obs.sv[i].values))
            else:
                continue

        # Melbourne-Wuebbena parameters (predetermined)
        I_max = 0.4  # Maximal ionospheric delay [m/h]
        k = 4  # criterion factor

        # For each tracked satellite
        for i in range(0, len(svG)):
            current_sat = obs.sel(sv=svG[i]).dropna(dim='time', how='all')

            L1 = current_sat['L1C']
            L2 = current_sat['L2W']
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

            print('Sat:', svG[i],
                  ', Slips:', slips_nr,
                  ', Slip criterion:',  criterion.values)

    def plot_slips(self, obs, sat_nr: str):
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


# %% Testing first algorithm
sliptest = Slips().slips_MW_single_freq(obs)

# %% Testing plot function
sliptest = Slips().plot_slips(obs, 'G08')
