# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:28:17 2022

@author: manue
"""
import numpy as np
from numpy import sin, cos
from mpmath import sec
from math import sqrt, nan, asin
from itertools import compress
from itertools import groupby

k = 40.3082
f1 = 1575.42*10**6  # [1/s]
f2 = 1227.60*10**6  # [1/s]
freqRatio = (f1**2*f2**2)/(k*(f1**2 - f2**2))
tecu = 10**(-16)
c = 299792458.0  # [m/s]
Lambda1 = c/f1  # [m]
Lambda2 = c/f2  # [m]
r_E = 6371e3  # Earth's radius [m]


class TEC:
    def __init__(self) -> None:
        # Data
        self.obsFile = []
        self.navFile = []
        self.receiverPos = []
        self.numObs = []
        self.numObsTrim = []
        # Relevant signals
        self.L1 = []
        self.L2 = []
        self.P1 = []
        self.P2 = []
        # Ephemeris class output
        self.elevation = []
        self.elevationMasked = []
        self.mask = []
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
        self.TECv = []
        self.threshold = []
        self.arcStartIndex = []
        self.arcEndIndex = []
        self.offset = []
        self.offsetBool = []
        self.sigma = []
        # Intermediate values
        self.height = []
        self.elevTrim = []
        self.heightTrim = []
        # Tests
        self.key = []
        self.group = []
        self.elevationMask = []
    
    def loadObsFile(self, obsFile, satellite):
        self.obsFile = obsFile.sel(sv=satellite).dropna(dim='time',how='all')
    
    def loadNavFile(self, navFile, satellite):
        self.navFile = navFile.sel(sv=satellite).dropna(dim='time',how='all')
    
    def loadData(self, obsFile, navFile, satellite, receiverPos, eph):
        self.receiverPos = np.array(receiverPos)
        self.loadObsFile(obsFile, satellite)
        self.loadNavFile(navFile, satellite)
        self.numObs = len(self.obsFile.L1.data)
        self.x = np.zeros(self.numObs, dtype = float)
        self.y = np.zeros(self.numObs, dtype = float)
        self.z = np.zeros(self.numObs, dtype = float)
        self.height = np.zeros(self.numObs, dtype = float)
        self.elevation = np.zeros(self.numObs, dtype = float)
        self.distance = np.zeros(self.numObs, dtype = float)
        self.azimuth = np.zeros(self.numObs, dtype = float)
        self.TECcp = np.zeros(self.numObs, dtype = float)
        self.TECpr = np.zeros(self.numObs, dtype = float)
        eph.loadEph(navFile, satellite)

    def extractObsData(self) -> None:
        # Updating Time Parameters
        self.L1 = self.obsFile.L1.data
        self.L2 = self.obsFile.L2.data
        self.P1 = self.obsFile.P1.data
        self.P2 = self.obsFile.P2.data

    def carrierPhaseTEC(self, L1, L2) -> float:
        if (not np.isnan(L1)) and (not np.isnan(L2)):
            deltaPhi = L1*Lambda1 - L2*Lambda2
            tec = deltaPhi*freqRatio*tecu
        else:
            tec = nan
        
        return tec

    def pseudoRangeTEC(self, P1, P2) -> float:
        if (not np.isnan(P1)) and (not np.isnan(P2)):
            deltaP = P2 - P1
            tec = deltaP*freqRatio*tecu
        else:
            tec = nan
        
        return tec
        
    def getRelativeTEC(self, eph):
        """
        Here is where the looping occurs :)

        Parameters
        ----------
        eph : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.extractObsData()
        for index, t in enumerate(self.obsFile.time.data):
            x_temp, y_temp, z_temp = eph.getXYZ(t, self.P1[index])
            e, a, d = eph.getEAD(x_temp, y_temp, z_temp, self.receiverPos)
            h = sqrt(x_temp**2 + y_temp**2 + z_temp**2) - r_E
            self.height[index] = h
            self.x[index] = x_temp
            self.y[index] = y_temp
            self.z[index] = z_temp
            self.elevation[index] = e
            self.azimuth[index] = a
            self.distance[index] = d
            self.TECcp[index] = self.carrierPhaseTEC(self.L1[index], self.L2[index])
            self.TECpr[index] = self.pseudoRangeTEC(self.P1[index], self.P2[index])  
        
    def offsetSingleArc(self, TECpr: float, TECcp: float, elevation: float):
        """
        Loops through the elevation array while the condition that the elevation
        is higher than 20º is met. Computes the offset and standard deviation
        of a given arc.

        Parameters
        ----------
        TECpr : float
            DESCRIPTION.
        TECcp : float
            DESCRIPTION.
        elevation : float
            DESCRIPTION.

        Returns
        -------
        idx : TYPE
            DESCRIPTION.

        """
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
                sumNum = sumNum + diff*sin(elevation[idx])
                sumDen = sumDen + sin(elevation[idx])
                
                # Quality control: computing weighted standard deviation
                num1 = num1 + sin(elevation[idx])*diff**2
                num2 = num2 + sin(elevation[idx])
                num3 = num3 + sin(elevation[idx])*diff
                den1 = den1 + sin(elevation[idx])
                den2 = den2 + sin(elevation[idx])**2
                
            idx += 1
            try:
                elevation[idx]
            except IndexError:
                break

        offset = sumNum/sumDen
        self.offset.append(offset)
        sigma = sqrt((num1*num2 - num3**2)/(den1**2 - den2))
        self.sigma.append(sigma)
        #idx = idx - 1

        return idx

    def findArc(self, elevation):
        """
        This function only loops through an input array of elevations and finds
        the position in the array where a threshold condition is met. All the 
        index management is done in the computeOffset function.

        Parameters
        ----------
        elevation : TYPE
            DESCRIPTION.

        Returns
        -------
        arcStart : TYPE
            DESCRIPTION.
        arcFlag : TYPE
            DESCRIPTION.

        """
        arcFlag = 0
        for idx, elv in enumerate(elevation):
            if elv >= self.threshold:
                arcStart = idx
                arcFlag = 1
                return arcStart, arcFlag
                break
        print("No more arcs were found")
        arcStart = nan
        return arcStart, arcFlag

    def computeOffset(self, threshold):
        """
        This function finds the arcs that meet the elevation condition and
        computes the offset of each arc. It also deals with the index management
        of the start and end position in the complete array.

        Parameters
        ----------
        threshold : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.threshold = threshold
        # Find first arc that meets the elevation condition, passing the whole
        # elevation array first
        arcStart, arcFLag = self.findArc(self.elevation)
        while arcFLag:
            self.arcStartIndex.append(arcStart)
            lastIndex = self.offsetSingleArc(self.TECpr[arcStart:],
                                             self.TECcp[arcStart:],
                                             self.elevation[arcStart:])
            print(self.offset)
            print(self.sigma)
            self.arcEndIndex.append(arcStart + lastIndex - 1)
            # Finding the next arc, now a reduced elevation array is passed
            # starting from the last index of the arc, onwards.
            nextArc, arcFLag = self.findArc(self.elevation[arcStart+lastIndex:])
            arcStart = arcStart + lastIndex + nextArc
    
    def offsetCorrectedTEC(self, threshold):
        """
        Corrects the offset in the arcs where the elevation condition is met.

        Parameters
        ----------
        threshold : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.mask = self.elevation > threshold
        self.mask.astype(np.int)
        self.numObsTrim = np.sum(self.mask)
        
        self.elevationMasked = self.elevation[self.mask]
        # self.elevationMasked = np.array(compress(self.elevation, self.mask))
        self.elevationMask = np.array([self.mask.T, self.elevation.T])
        for key, group in groupby(self.elevationMask, lambda x : x[0]):
            self.key.append(key)
            self.group.append(group)
        
        
        self.computeOffset(threshold)
        for i, off in enumerate(self.offset):
            self.TECr.append(self.TECcp[self.arcStartIndex[i]: \
                                        self.arcEndIndex[i]] + off)
            self.elevTrim.append(self.elevation[self.arcStartIndex[i]: \
                                        self.arcEndIndex[i]])
            self.heightTrim.append(self.height[self.arcStartIndex[i]: \
                                        self.arcEndIndex[i]])
                
    def trimArray(self):
        # for i in range(len(tec[sat].TECr)):
        #    tec[sat].obsFile.time[tec[sat].arcStartIndex[i]:tec[sat].arcEndIndex[i]], 
        #    tec[sat].TECr[i],'r')
        pass

    def getVerticalTEC(self):
        self.offsetBool = sum((self.TECpr[self.mask] - self.TECcp[self.mask])*sin(self.elevationMasked[:]))/sum(sin(self.elevationMasked[:]))
        for i in range((len(self.offset))):
            for j in range(len(self.elevTrim[i][:])):
                # self.TECv[i][j] = self.TECr[i][j] / sec(math.asin(r_E*math.cos(self.elevTrim[i][j])/(r_E + self.heightTrim[i][j])))
                self.TECv.append(self.TECr[i][j] / sec(asin(r_E*cos(self.elevTrim[i][j])/(r_E + self.heightTrim[i][j])))) 
