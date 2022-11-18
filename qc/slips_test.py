# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 10:42:27 2022

@author: Magdalena
"""

# Link to MW algorithm
# https://gssc.esa.int/navipedia/index.php/Detector_based_in_code_and_carrier_phase_data:_The_Melbourne-W%C3%BCbbena_combination

import georinex as gr
import numpy as np
from matplotlib.pyplot import figure, show
import matplotlib.pyplot as plt
import math
import datetime
import xarray as xr
import copy

# %% Load new test data from Sarah

# G01 has 1 cycle slip from 00.00-00.15, 
# 1 cycle slip from 19.00-19.15, 
# 6 cycle slips from 19.15-19.30 
# and 2 cycle slips from 22.15-22.30

obs = gr.load("tests/test_data/Rinex2/klsq0640.22o", fast=False,
              tlim=['2022-03-05T19:15:00', '2022-03-05T19:30:00'])

# hdr = gr.rinexheader("tests/test_data/Rinex2/klsq0640.22o")
# rnx_version = 2

sat = 'G01'
data = obs.sel(sv=sat).dropna(dim='time', how='all')

# %% Plot the data to see cycle slips as suggested

ax = figure(figsize=(10, 6)).gca()
# ax.plot(data.time, data.data_vars['C1'], linestyle='', marker='o', label="C1")
# ax.plot(data.time, data.data_vars['L1'], linestyle='', marker='o', label="L1")
# ax.plot(data.time, data.data_vars['L2'], linestyle='', marker='o', label="L2")
ax.plot(data.time, data.data_vars['P2'], linestyle='', marker='o', label="P2")
# ax.plot(data.time, data.data_vars['P1'], linestyle='', marker='o', label="P1")
# ax.plot(data.time, data.data_vars['C2'], linestyle='', marker='o', label="C2")
# ax.plot(data.time, data.data_vars['C5'], linestyle='', marker='o', label="C5")
# ax.plot(data.time, data.data_vars['L5'], linestyle='', marker='o', label="L5")
ax.grid()
ax.legend()
show()


# %% Slip detection (Melbourne-Wubbena)

# define variables
freq = [1575.42, 1227.60, 1176.45]  # L1, L2, L5 for GPS
f1 = freq[0]*1e6
f2 = freq[1]*1e6

min_arc = 60  # seconds
c = 299792458  # Speed of light
lambda_w = c/(f1-f2)  # meter; wavelength of the wide-lane combination
k_factor = 4

# Defining variables at the initial index (B_w, mean, sigma, th and d)

B_w = np.empty(len(data.data_vars['L1']))
L_w = (1/(f1-f2))*(f1*data.data_vars['L1'][0] - f2*data.data_vars['L2'][0])
C_n = (1/(f1+f2))*(f1*data.data_vars['P1'][0] + f2*data.data_vars['P2'][0])
B_w[0] = L_w - C_n

mean_Bw = np.empty(len(data.time))
mean_Bw[0] = 0

sigma_Bw = np.empty(len(data.time))
sigma_Bw[0] = lambda_w/2

th = np.empty(len(data.time))
th[0] = 0

d = np.empty(len(data.time))
d[0] = 0

slips = []

for i in range(1, len(data.time)):
    arc_length = i*data.interval

    L_w = (1/(f1-f2))*(f1*data.data_vars['L1'][i] -
                       f2*data.data_vars['L2'][i])
    C_n = (1/(f1+f2))*(f1*data.data_vars['P1'][i] +
                       f2*data.data_vars['P2'][i])
    # Wide-lane ambiguity
    B_w[i] = L_w - C_n

    # Is arc length long enough? Are data intervals small enough?
    if int(data.time[i]-data.time[i-1])*1e-9 < min_arc:  # less than 60s:
        # update difference d
        d[i] = B_w[i] - mean_Bw[i-1]

        # update threshold th
        th[i] = k_factor * np.sqrt(sigma_Bw[i-1])

        # check for cycle slip
        if (np.abs(d[i]) > th[i]): #  and (np.sqrt(sigma_Bw[i-1]) <= lambda_w):  # or (np.isnan(B_w[i])):
            slips.append(0)  # 0 for slip
        else:
            slips.append(1)  # 1 for ok values

        mean_Bw[i] = ((i-1)/i)*mean_Bw[i-1] + (1/i)*B_w[i]
        sigma_Bw[i] = ((i-1)/i)*sigma_Bw[i-1] + (1/i)*(d[i]**2)

    else:
        # If not, reset variables
        slips.append(2)  # 2 for LLI

        d[i] = B_w[i] - mean_Bw[i-1]
        mean_Bw[i] = 0
        sigma_Bw[i] = lambda_w/2

# Out message, to be displayed after the function is used
times = data.time.values
out2 = 0  # LLI
out0 = 0  # slip
faulty_out = [[], []]
for i in range(0, len(slips)):
    if slips[i] == 2:
        out2 = out2 + 1
        faulty_out[0].append(times[i])
    elif slips[i] == 0:
        out0 = out0 + 1
        faulty_out[1].append(times[i-1])
print('MW algorithm complete, LLI:', out2, ', slips: ', out0)


# %% Plot different results

# %% Plot slips
xaxis = np.linspace(0, np.size(slips), np.size(slips), endpoint=True)
ax = figure(figsize=(10, 6)).gca()
ax.plot(xaxis, slips, marker='o')
ax.grid()
show()

# %% Plot threshold and differences
xdata = mean_Bw
xaxis = data.time
ax = figure(figsize=(10, 6)).gca()
ax.plot(xaxis, d, label="d", linestyle='', marker='o')
ax.plot(xaxis, th, label="th", linestyle='', marker='o')
ax.grid()
ax.legend()
plt.title("Slip when d > th")
show()

# %% More plotting
# Sigma values must be larger than lambda_w threshold, otherwise a cycle slip would be detected.
xdata = sigma_Bw
lambda_w = 0.86
xaxis = data.time

ax = figure(figsize=(10, 6)).gca()
ax.plot(xaxis, np.sqrt(xdata), linestyle='', marker='o', label="sigma_Bw")
plt.axhline(y=lambda_w, color='r', linestyle='-', label="lambda")
ax.grid()
ax.legend()
plt.title("Slip when sigma_Bw < lambda line")
show()

# %% Plot B_w
xaxis = data.time

ax = figure(figsize=(10, 6)).gca()
ax.plot(xaxis, B_w, linestyle='', marker='o', label="B_w")
ax.grid()
ax.legend()
plt.title("B_w")
show()

