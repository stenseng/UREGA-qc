# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:28:17 2022

@author: manue
"""
import georinex as gr
import numpy as np
import math

k = 40.3082
f1 = 1575.42*10**6  # [1/s]
f2 = 1227.60*10**6  # [1/s]
freqRatio = (f1**2*f2**2)/(k*(f1**2 - f2**2))
tecu = 10**(-16)
c = 299792458.0  # [m/s]
Lambda1 = c/f1  # [m]
Lambda2 = c/f2  # [m]


class TEC:
    def __init__(self) -> None:
        # Data
        self.obsFile = []
        self.navFile = []
        self.receiverPos = []
        # Relevant signals
        self.L1 = []
        self.L2 = []
        self.P1 = []
        self.P2 = []
        # Ephemeris class output
        self.elevation = []
        self.azimuth = []
        self.distance = []
        self.x = []
        self.y = []
        self.z = []
        # TEC values
        self.TECcp = []
        self.TECpr = []
        self.TECr = []
        self.TECs = []
        self.threshold = []
        self.arcStartIndex = []
        self.arcEndIndex = []
        self.offset = []
        self.sigma = []
    
    def loadObsFile(self, obsFile, satellite):
        self.obsFile = obsFile.sel(sv=satellite).dropna(dim='time',how='all')
    
    def loadNavFile(self, navFile, satellite):
        self.navFile = navFile.sel(sv=satellite).dropna(dim='time',how='all')
    
    def loadData(self, obsFile, navFile, satellite, receiverPos, eph):
        self.receiverPos = receiverPos
        self.loadObsFile(obsFile, satellite)
        self.loadNavFile(navFile, satellite)
        eph.loadEph(navFile, satellite)

    def extractObsData(self) -> None:
        # Updating Time Parameters
        self.L1 = self.obsFile.L1.data
        self.L2 = self.obsFile.L2.data
        self.P1 = self.obsFile.P1.data
        self.P2 = self.obsFile.P2.data

    def carrierPhaseTEC(self, L1, L2) -> None:
        if (not np.isnan(L1)) and (not np.isnan(L2)):
            deltaPhi = L1*Lambda1 - L2*Lambda2
            tec = deltaPhi*freqRatio*tecu
        else:
            tec = math.nan
        self.TECcp.append(tec)

    def pseudoRangeTEC(self, P1, P2) -> None:
        if (not np.isnan(P1)) and (not np.isnan(P2)):
            deltaP = P2 - P1
            tec = deltaP*freqRatio*tecu
        else:
            tec = math.nan
        self.TECpr.append(tec)
        
    def getRelativeTEC(self, eph):
        self.extractObsData()
        for index, t in enumerate(self.obsFile.time.data):
            x_temp, y_temp, z_temp = eph.getXYZ(t, self.P1[index])
            e, a, d = eph.getEAD(x_temp, y_temp, z_temp, self.receiverPos)
            self.x.append(x_temp)
            self.y.append(y_temp)
            self.z.append(z_temp)
            self.elevation.append(e)
            self.azimuth.append(a)
            self.distance.append(d)
            self.carrierPhaseTEC(self.L1[index], self.L2[index])
            self.pseudoRangeTEC(self.P1[index], self.P2[index])

    def offsetSingleArc(self, TECpr: float, TECcp: float, elevation: float):
        sumNum = 0
        sumDen = 0
        num1 = 0
        num2 = 0
        num3 = 0
        den1 = 0
        den2 = 0
        idx = 0
        while elevation[idx] >= self.threshold:
            if (not np.isnan(TECpr[idx])) and (not np.isnan(TECcp[idx])):
                diff = TECpr[idx] - TECcp[idx]
                sumNum = sumNum + diff*math.sin(elevation[idx])
                sumDen = sumDen + math.sin(elevation[idx])
                
                # Quality control: computing weighted standard deviation
                num1 = num1 + math.sin(elevation[idx])*diff**2
                num2 = num2 + math.sin(elevation[idx])
                num3 = num3 + math.sin(elevation[idx])*diff
                den1 = den1 + math.sin(elevation[idx])
                den1 = den1 + math.sin(elevation[idx])**2
                
            idx += 1
            try:
                elevation[idx]
            except IndexError:
                break

        offset = sumNum/sumDen
        self.offset.append(offset)
        sigma = (num1*num2 - num3**2)/(den1**2 - den2)
        self.sigma.append(sigma)

        return idx

    def findArc(self, elevation):
        arcFlag = 0
        for idx, elv in enumerate(elevation):
            if elv >= self.threshold:
                arcStart = idx
                arcFlag = 1
                return arcStart, arcFlag
                break
        print("No more arcs were found")
        arcStart = math.nan
        return arcStart, arcFlag

    def computeOffset(self, threshold):
        self.threshold = threshold
        arcStart, arcFLag = self.findArc(self.elevation)
        while arcFLag:
            self.arcStartIndex.append(arcStart)
            lastIndex = self.offsetSingleArc(self.TECpr[arcStart:],
                                             self.TECcp[arcStart:],
                                             self.elevation[arcStart:])
            print(self.offset)
            print(self.sigma)
            self.arcEndIndex.append(arcStart + lastIndex - 1)
            nextArc, arcFLag = self.findArc(self.elevation[arcStart+lastIndex:])
            arcStart = arcStart + lastIndex + nextArc
    
    def offsetCorrectedTEC(self, threshold):
        self.computeOffset(threshold)
        for i, off in enumerate(self.offset):
            self.TECr.append(self.TECcp[self.arcStartIndex[i]: \
                                        self.arcEndIndex[i]] + off)
    def mappingFunction(self):
        pass
                
