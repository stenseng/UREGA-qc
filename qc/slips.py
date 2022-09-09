#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Lars Stenseng.

@mail: lars@stenseng.net
"""
# line 853 testing functions
# from qc.__version__ import __version__

import georinex as gr
import numpy as np
from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt

# %%
obs = gr.load(
    'tests/test_data/Rinex3/KLSQ00GRL_R_20213070000_01D_15S_MO.rnx',
    # tlim=['2021-11-03T12:00', '2021-11-03T12:30'])
    tlim=['2021-11-03T05:30', '2021-11-03T07:30'])
    # tlim=['2021-11-03T01:00', '2021-11-03T02:00'])
# tlim=['2021-11-03T15:00', '2021-11-03T18:00'])
# hdr = gr.rinexheader(
#    'tests/test_data/Rinex3/KLSQ00GRL_R_20213070000_01D_15S_MO.rnx')
# rnx_version = 3

# %% Test with self-collected data which is sure to contain large slips
obs = gr.load('../../survey_data/SEPT0622.22O')

# %% Testing with NASA data
obs = gr.load('../../survey_data/ABMF00GLP_R_20160010000_01D_30S_MO.rnx',
              tlim=['2016-01-01T05:30', '2016-01-01T07:30'])
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
print(svG)
# %%
# 5:30 to 7:30, G08 and G21 give 2 cycle slips # [290:300]
# 01:30 to 01:50, R01 gives cycle slip

# 'G01','G06','G08','G10','G12','G14','G17','G19','G21','G22','G24','G30','G32'
#sat = 'G21'
sat = 'G05'
sattest = obs.sel(sv=sat).dropna(dim='time', how='all')
# G02 data vars with no-nan: C1C, D1C, L1C, S1C, C1W, C2W, D2W, L2W, S1W, S2W

I_max = 0.4  # Maximal ionospheric delay [m/h]
k = 4  # criterion factor

L1 = sattest['L1C']  # GPS
L2 = sattest['L2W']  # GPS

# L1 = sattest['L1C']  # glonass
# L2 = sattest['L2W']  # 

# L1 = L1test
# L2 = L2slip
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
    L2test = test['L2C']
    ax.plot(L1test.time, L1test, label=svG[i], linewidth=2.0)
    #ax.plot(L2test.time, L2test, label='L2', linewidth=0.5)
ax.grid()
ax.legend()
plt.xlabel('Time [epochs]')
plt.ylabel('Carrier phases')
show()

# Glonass info
# const_def = ['C1C', 'C2C', 'C3I', 'L1C', 'L2C']
# f1 = (1602 + (k*9/16))
# f2 = (1246 + (k*7/16))
# f3 = 1202.025
# freq = [f1, f2, f3]
# %%
# Plot separate sats, L1 and L2
# ['C11', 'E22', 'G01', 'G03', 'G06', 'G07', 'G08', 'G09', 'G11', 'G16',
#       'G17', 'G19', 'G23', 'G27', 'G28', 'G30', 'G32', 'R04', 'R05', 'R06',
#       'R09', 'R10', 'R11', 'R15', 'R16', 'R18', 'R19', 'R20', 'R21', 'S20',
#       'S33', 'S35', 'S38']

ax = figure(figsize=(10, 6)).gca()
test = obs.sel(sv='G07').dropna(dim='time', how='all')
L1test = test['L1C']
L2test = test['L2W']
C1test = test['C1C']
C2test = test['C2W']
ax.plot(L1test.time, L1test, label='L1', linewidth=5.0)
ax.plot(L2test.time, L2test, label='L2', linewidth=3.0)
ax.plot(C1test.time, C1test, label='C1', linewidth=5.0)
ax.plot(C2test.time, C2test, label='C2', linewidth=2.0)
ax.grid()
ax.legend()
plt.xlabel('Time [epochs]')
plt.ylabel('Carrier phases')
show()

# %% Detecting and noting location of loss-of-locks
# satnr = 'G13'
satnr = 'R09'
test = obs.sel(sv=satnr).dropna(dim='time', how='all')
L1test = test['L1C']
L2test = test['L2W']
C1test = test['C1C']
C2test = test['C2W']
lls = []
for i in range(0, len(L1test)):
    if np.isnan(L1test[i]):
        lls.append(L1test[i])

print(lls)

# %% Simulating small and large cycle slips, as well as loss of lock
# obs = gr.load(
#    'tests/test_data/Rinex3/KLSQ00GRL_R_20213070000_01D_15S_MO.rnx',
#    tlim=['2021-11-03T05:30', '2021-11-03T07:30'])

satnr = 'G14'
test = obs.sel(sv=satnr).dropna(dim='time', how='all')
ax = figure(figsize=(10, 6)).gca()
L1test = test['L1C']
L2test = test['L2W']
C1test = test['C1C']
C2test = test['C2W']
L2test[190:250] = test['L2W'][190:250] - 2000000
L2test[320:370] = test['L2W'][320:370] - 8000000
L1test[300:] = test['L1C'][300:] - 2000000
C1test[100:150] = float('NaN')
C2test[350:] = test['C2W'][350:] + 4000000
ax.plot(L1test.time, L1test, label='L1', linewidth=5.0)
ax.plot(L2test.time, L2test, label='L2', linewidth=3.0)
ax.plot(C1test.time, C1test, label='C1', linewidth=5.0)
ax.plot(C2test.time, C2test, label='C2', linewidth=2.0)
ax.grid()
ax.legend()
plt.xlabel('Time [epochs]')
plt.ylabel('Carrier phases')
show()

# %% Noisy dataset
obs_noise = obs.sel(sv='G14').dropna(dim='time', how='all')
obs_noise['L2W'][190:250] = obs_noise['L2W'][190:250] - 2000000
obs_noise['L2W'][320:370] = obs_noise['L2W'][320:370] - 8000000
obs_noise['L1C'][300:] = obs_noise['L1C'][300:] - 2000000
obs_noise['L1C'][100:150] = float('NaN')
obs_noise['C2W'][350:] = obs_noise['C2W'][350:] + 4000000
obs_noise['C1C'][100:150] = float('NaN')

ax = figure(figsize=(10, 6)).gca()
L1test = obs_noise['L1C']
L2test = obs_noise['L2W']
C1test = obs_noise['C1C']
C2test = obs_noise['C2W']
ax.plot(L1test.time, L1test, label='L1', linewidth=5.0)
ax.plot(L2test.time, L2test, label='L2', linewidth=3.0)
ax.plot(C1test.time, C1test, label='C1', linewidth=5.0)
ax.plot(C2test.time, C2test, label='C2', linewidth=2.0)
ax.grid()
ax.legend()
plt.xlabel('Time [epochs]')
plt.ylabel('Carrier phases')
plt.title('Noisy dataset')
show()

# %% Testing noisy dataset (geometry-free combination) returning
# small and big cycle slips and loss of lock timestamps
obs_noise = obs.sel(sv='G09').dropna(dim='time', how='all')
L1test = obs_noise['L1C']
L2test = obs_noise['L2W']

I_max = 0.4  # Maximal ionospheric delay [m/h]
k_1 = 3  # criterion factor
k_2 = 0.5  # criterion factor for small cycle slip type

L1 = L1test
L2 = L2test
L4 = np.abs(L1 - L2)

sigma_L4 = np.std(L4)

criterion_1 = k_1*sigma_L4 + I_max
criterion_2 = k_2*sigma_L4 + I_max
slips_nr_1 = 0
slips_nr_2 = 0
lls_time = []

L4_diff = []
for i in range(1, len(L4)):
    L4_diff.append(np.abs(L4[i] - L4[i-1]))
    if (np.abs(L4[i] - L4[i-1]) > criterion_1):
        # If satisfied, raise cycle-slip flag
        slips_nr_1 = slips_nr_1 + 1
    elif (np.abs(L4[i] - L4[i-1]) > criterion_2):
        slips_nr_2 = slips_nr_2 + 1
    elif np.isnan(L4[i].values) and not np.isnan(L4[i-1].values):
        lls_time.append(L4[i].time)
    # elif np.isnan(L4[i].values) and np.isnan(L4[i+1].values)==False:
    #   lls_time.append(L4[i].time)

ax = figure(figsize=(10, 6)).gca()
ax.plot(L2.time[1:], L4_diff, label=satnr)
plt.axhline(y=criterion_1, label='Large slip limit', linestyle='-', color='r')
plt.axhline(y=criterion_2, label='Low slip limit', linestyle='-', color='m')
ax.grid()
ax.legend()
plt.xlabel('Time [epochs]')
plt.ylabel('L4')
plt.title('Single-frequency Melbourne-Wuebbena')
show()

print('Large slips:', slips_nr_1, ', Slip criterion:',  criterion_1.values)
print('Small slips:', slips_nr_2, ', Slip criterion:', criterion_2.values)
print('Loss of locks:', len(lls_time), ', Timestamps:\n',
      [var.time.values for var in lls_time])

# %% Testing slip repair
slip_idx = []
for i in range(1, len(L4)):
    if (np.abs(L4[i] - L4[i-1]) > criterion_1):
        slip_idx = [i, L4[i].values]

# %% Simulating cycle slip in data

ax = figure(figsize=(10, 6)).gca()
test = obs.sel(sv='G13').dropna(dim='time', how='all')

L1test = test['L1C']
L2slip = test['L2W']
L2slip[190:250] = test['L2W'][190:250] - 2000000
L2slip[320:370] = test['L2W'][320:370] - 3000000
#L2slip[19:25] = test['L2W'][19:25] - 2000000
#L2slip[32:37] = test['L2W'][32:37] + 3000000

C1test = test['C1C']
C2slip = test['C2W']
C2slip[190:250] = test['C2W'][190:250] - 2000000
C2slip[320:370] = test['C2W'][320:370] - 3000000
#C2slip[19:25] = test['C2W'][19:25] - 2000000
#C2slip[32:37] = test['C2W'][32:37] + 3000000


ax.plot(L1test.time, L1test, label='L1', linewidth=5.0)
# ax.plot(L2test.time, L2test, label='L2', linewidth=3.0)
ax.plot(L2slip.time, L2slip, label='L2slip', linewidth=3.0)
ax.plot(C1test.time, C1test, label='C1', linewidth=5.0)
ax.plot(C2slip.time, C2slip, label='C2slip', linewidth=2.0)

ax.grid()
ax.legend()
plt.xlabel('Time [epochs]')
plt.ylabel('Carrier phases')
show()

# %% Calculating geometry-free combination

I_max = 0.4  # Maximal ionospheric delay [m/h]
k = 1  # criterion factor

L1 = L1test
L2 = L2slip
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


# %% Wide-lane Melbourne-Wuebbena testing

# 'G01','G06','G08','G10','G12','G14','G17','G19','G21','G22','G24','G30','G32'
sat = 'G13'
sattest = obs.sel(sv=sat).dropna(dim='time', how='all')
# G02 data vars with no-nan: C1C, D1C, L1C, S1C, C1W, C2W, D2W, L2W, S1W, S2W

freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS

f1 = freq[0]*1e6
f2 = freq[1]*1e6

#P1 = sattest['C1C']
#P2 = sattest['C2W']
#L1 = sattest['L1C']  # GPS
#L2 = sattest['L2W']  # GPS
P1 = C1test
P2 = C2slip
L1 = L1test
L2 = L2slip

# L1 = sattest['L1C']  # Galileo
# L2 = sattest['L8Q']  # Galileo
L6 = (1/(f1-f2))*(f1*L1 - f2*L2) - (1/(f1+f2))*(f1*P1 + f2*P2)

sigma_L6 = np.std(L6)
k = 1  # criterion factor

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


# %% Testing first algorithm
sliptest = Slips().slips_MW_single_freq(obs_noise)

# %% Testing plot function
sliptest = Slips().plot_slips(obs_noise, 'G14')

# %% ========================== testing functions =============================

testsats = sort_sat_types(obs, 'G')
data_vars = get_const_data_vars('G14', 'G')
testobs = select_default_observables(obs, data_vars[0], 3, 'G')
# %%

def slips_MW_single_freq(obs):
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

# %%
def plot_slips(obs, sat_nr: str):
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



def select_default_observables(obs, const_def, rnx_version, constellation):
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
    if rnx_version == '3':
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

    if rnx_version == 3 or (
            rnx_version == 2 and not constellation == 'R'):
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

