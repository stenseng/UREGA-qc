# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 11:19:25 2022

@author: manue
"""
import georinex as gr
import numpy as np
import math
from ephemeris import Eph
from tec import TEC
import matplotlib.pyplot as plt

# READING FILES
dataPath = "F:/Space/RINEX files/Rinex2\\"
navFile = "klsq3070.21n"
obsFile = "klsq3070.21o"

obs = gr.load(dataPath + obsFile, use='G')
nav = gr.load(dataPath + navFile, use='G')
hdr = gr.rinexheader(dataPath + obsFile)
receiverPos = hdr["position"]
sat = 'G27'
tec = {}
eph = {}

tec[sat] = TEC()
eph[sat] = Eph()
tec[sat].loadData(obs, nav, sat, receiverPos, eph[sat])
tec[sat].getRelativeTEC(eph[sat])
tec[sat].offsetCorrectedTEC(threshold=math.radians(20))
tec[sat].getVerticalTEC()

# dcb = open('CAS0MGXRAP_20213070000_01D_01D_DCB.BSX')

figure = plt.figure()
axes = figure.add_subplot(1, 1, 1)
plt.scatter(tec[sat].obsFile.time, tec[sat].TECpr, marker='.')
plt.title("STEC - Pseudorange")
plt.xlabel('Time') 
plt.ylabel('TECU')
plt.setp(axes.get_xticklabels(), rotation = 35)
plt.show()
for i in range(len(tec[sat].TECr)):
    if tec[sat].sigma[i] < 5:
        plt.plot(tec[sat].obsFile.time[tec[sat].arcStartIndex[i]:tec[sat].arcEndIndex[i]], 
                 tec[sat].TECr[i],'r')
        plt.axvline(tec[sat].obsFile.time.data[tec[sat].arcStartIndex[i]],linewidth=1,color='lime',ls='--')
        plt.axvline(tec[sat].obsFile.time.data[tec[sat].arcEndIndex[i]],linewidth=1,color='lime',ls='--')
        plt.title('STEC (Satellite: %s)' %sat)
        plt.xlabel('Time') 
        plt.ylabel('TECU') 
        plt.setp(axes.get_xticklabels(), rotation = 35)
        plt.show()
# %%
satellites = obs.sv.data
tec = {}
eph = {}

for sat in satellites:
    print(sat)
    tec[sat] = TEC()
    eph[sat] = Eph()
    tec[sat].loadData(obs, nav, sat, receiverPos, eph[sat])
    tec[sat].getRelativeTEC(eph[sat])
    tec[sat].offsetCorrectedTEC(threshold=math.radians(20))

dcb = open('CAS0MGXRAP_20213070000_01D_01D_DCB.BSX')

figure = plt.figure()
axes = figure.add_subplot(1, 1, 1)
for sat in satellites:
    plt.scatter(tec[sat].obsFile.time, tec[sat].TECpr, marker='.')
    plt.title("STEC - Pseudorange")
    plt.xlabel('Time') 
    plt.ylabel('TECU')
    plt.setp(axes.get_xticklabels(), rotation = 35)
    plt.show()
    for i in range(len(tec[sat].TECr)):
        if tec[sat].sigma[i] < 5:
            plt.plot(tec[sat].obsFile.time[tec[sat].arcStartIndex[i]:tec[sat].arcEndIndex[i]], 
                     tec[sat].TECr[i],'r')
            #plt.axvline(tec[sat].obsFile.time.data[tec[sat].arcStartIndex[i]],linewidth=1,color='lime',ls='--')
            #plt.axvline(tec[sat].obsFile.time.data[tec[sat].arcEndIndex[i]],linewidth=1,color='lime',ls='--')
            plt.title("STEC - Carrier Phase")
            plt.xlabel('Time') 
            plt.ylabel('TECU') 
            plt.setp(axes.get_xticklabels(), rotation = 35)
            plt.show()
#%% Sensitivity analysis
x_dif = []
y_dif = []
z_dif = []

for i in range(len(tec[sat].x) - 1):
    x_dif.append(tec[sat].x[i+1] - tec[sat].x[i])
    y_dif.append(tec[sat].y[i+1] - tec[sat].y[i])
    z_dif.append(tec[sat].z[i+1] - tec[sat].z[i])

xMax = max(x_dif)
xMaxI = x_dif.index(xMax)
yMax = max(y_dif)
yMaxI = y_dif.index(yMax)
zMax = max(z_dif)
zMaxI = z_dif.index(zMax)

xMin = min(x_dif)
xMinI = x_dif.index(xMin)
yMin = min(y_dif)
yMinI = y_dif.index(yMin)
zMin = min(z_dif)
zMinI = z_dif.index(zMin)

ax = plt.subplot()
ax.plot(x_dif[0:xMaxI-1])
ax.plot(x_dif[xMaxI+1:xMinI-1])
ax.plot(x_dif[xMinI+1:2329])
plt.show()

ax = plt.subplot()
ax.plot(y_dif)
plt.show()

ax = plt.subplot()
ax.plot(z_dif)
plt.show()

ax = plt.subplot()
ax.scatter(tec[sat].navFile.time, tec[sat].navFile.Toe.data)
plt.setp(ax.get_xticklabels(), rotation=35)
plt.show()
# %% ISMR check
ISMR = open('KLQ2307a00_21.csv')
ismr00 = np.loadtxt(ISMR, delimiter=",")
ISMR = open('KLQ2307a15_21.csv')
ismr15 = np.loadtxt(ISMR, delimiter=",")

x_temp, y_temp, z_temp = eph[sat].getXYZ(259260, 0)
e, a, d = eph[sat].getEAD(x_temp, y_temp, z_temp, tec[sat].receiverPos)
print(np.rad2deg(e))
print(np.rad2deg(a))
# %% TEC Stuff
tec[sat].getRelativeTEC(tec[sat].TECpr, tec[sat].TECcp, tec[sat].elevation)
# aaaa = elevation[2320:]
x, y, z = eph[sat].getXYZ(259260.0, P1=0)
print(eph[sat].getEAD(x, y, z, tec[sat].receiverPos))
# %% PLOTTING
R_e = 6371000  # [m]
ax = plt.subplot(111, projection="3d")
ax.scatter(tec[sat].x, tec[sat].y, tec[sat].z, marker="o", color="b", linewidths=0.1)
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]

x_e = R_e * np.cos(u) * np.sin(v)

y_e = R_e * np.sin(u) * np.sin(v)

z_e = R_e * np.cos(v)
ax.plot_wireframe(x_e, y_e, z_e, color="g")
ax.scatter(tec[sat].receiverPos[0], tec[sat].receiverPos[1], tec[sat].receiverPos[2],
           marker="o", color="r", linewidth=0.1)
ax.set_xlim3d(-5 * R_e, 5 * R_e)
ax.set_ylim3d(-5 * R_e, 5 * R_e)
ax.set_zlim3d(-5 * R_e, 5 * R_e)
plt.show()
ax.set_zlim3d(-5 * R_e, 5 * R_e)
plt.show()
ax.set_zlim3d(-5 * R_e, 5 * R_e)
plt.show()
ax.set_zlim3d(-5 * R_e, 5 * R_e)
plt.show()
plt.show()

ax = plt.subplot()
ax.scatter(tec[sat].obsFile.time, tec[sat].azimuth, linewidths=0.01)
plt.setp(ax.get_xticklabels(), rotation=35)
plt.show()

ax = plt.subplot()
ax.scatter(tec[sat].obsFile.time, tec[sat].distance)
plt.setp(ax.get_xticklabels(), rotation=35)
plt.show()

ax = plt.subplot()
ax.plot(tec[sat].obsFile.time, tec[sat].elevation, "ob")
plt.setp(ax.get_xticklabels(), rotation=35)
plt.show()

ax = plt.subplot()
ax.plot(tec[sat].TECr[0])
plt.show()

"""
figure = plt.figure()
axes = figure.add_subplot(1, 1, 1)

plt.plot(tec[sat].obsFile.time, tec[sat].TECcp,'r')
plt.title("STEC - Carrier Phase")
plt.xlabel('Time') 
plt.ylabel('TECU') 
plt.setp(axes.get_xticklabels(), rotation = 35)
plt.show()
"""