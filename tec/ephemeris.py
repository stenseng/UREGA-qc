# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 13:09:57 2022

@author: manue
"""
import numpy as np
import math
from datetime import datetime
import pandas


c = 299792458.0  # [m/s]
mu = 3.986005e14  # [m^3/s^2]
Omegadot_Earth = 7.2921151467e-5  # [rad/s]
pi_gps = 3.1415926535898  # Pi in the gps system


class Eph:
    def __init__(self) -> None:
        # Time Parameters
        self.toe = []
        self.toc = []
        self.a0 = []
        self.a1 = []
        self.a2 = []
        # Keplerian Parameters
        self.a = []
        self.e = []
        self.i = []
        self.omega = []
        self.Omega = []
        self.M = []
        self.theta = []
        # Perturbation Parameters
        self.deltaN = []
        self.Omegadot = []
        self.idot = []
        self.Cus = []
        self.Cuc = []
        self.Cis = []
        self.Cic = []
        self.Crs = []
        self.Crc = []
        # Debugging
        self.set = []

    def newton(self, M: float, e: float) -> float:
        x0 = 0.01
        tol = 1e-13
        maxIter = 100
        xn = x0
        for n in range(0, maxIter):
            fxn = self.MtoE(xn, M, e)
            if abs(fxn) < tol:
                return xn
            Dfxn = self.dfMtoE(xn, e)
            if Dfxn == 0:
                print("Zero derivative. No solution found.")
                return np.nan
            xn = xn - fxn / Dfxn
        print("Exceeded maximum iterations. No solution found.")
        return np.nan

    def MtoE(self, E: float, M: float, e: float) -> float:
        f = E - M - e * math.sin(E)
        return f

    def dfMtoE(self, E: float, e: float) -> float:
        df = 1 - e * math.cos(E)
        return df

    def trueAnomaly(self, M: float, e: float) -> float:
        """
        Iteratively solving for the True anomaly knowing the mean anomaly (M).
        First M is converted to Eccentric anomaly (E) and then the Newton's
        method is employed to obtain the True anomaly

        Parameters
        ----------
        M : float
            Mean anomaly [rad].
        e : float
            Eccentricity of the orbit [-].

        Returns
        -------
        theta: float
            True anomaly [rad].

        """
        E = self.newton(M, e)

        sin_theta = math.sqrt(1 - e**2) * math.sin(E) / (1 - e * math.cos(E))
        cos_theta = (math.cos(E) - e) / (1 - e*math.cos(E))
        theta = math.atan2(sin_theta, cos_theta)

        return theta
    
    def loadEph(self, navFile, satellite):
        nav = navFile.sel(sv=satellite).dropna(dim='time',how='all')
        for setOfKepl in range(len(nav.Toe)):
            self.addNav(
                nav.Toe[setOfKepl].data,
                nav.time[setOfKepl].data,
                nav.SVclockBias[setOfKepl].data,
                nav.SVclockDrift[setOfKepl].data,
                nav.SVclockDriftRate[setOfKepl].data,
                nav.sqrtA[setOfKepl].data**2,
                nav.Eccentricity[setOfKepl].data,
                nav.Io[setOfKepl].data,
                nav.omega[setOfKepl].data,
                nav.Omega0[setOfKepl].data,
                nav.M0[setOfKepl].data,
                nav.DeltaN[setOfKepl].data,
                nav.OmegaDot[setOfKepl].data,
                nav.IDOT[setOfKepl].data,
                nav.Cus[setOfKepl].data,
                nav.Cuc[setOfKepl].data,
                nav.Cis[setOfKepl].data,
                nav.Cic[setOfKepl].data,
                nav.Crs[setOfKepl].data,
                nav.Crc[setOfKepl].data,
            )

    def addNav(
        self,
        toe: float,
        toc: float,
        a0: float,
        a1: float,
        a2: float,
        a: float,
        e: float,
        i: float,
        omega: float,
        Omega: float,
        M: float,
        deltaN: float,
        Omegadot: float,
        idot: float,
        Cus: float,
        Cuc: float,
        Cis: float,
        Cic: float,
        Crs: float,
        Crc: float
    ) -> None:
        # Updating Time Parameters
        self.toe.append(toe)
        self.toc.append(toc)
        self.a0.append(a0)
        self.a1.append(a1)
        self.a2.append(a2)
        # Updating Keplerian Parameters
        self.a.append(a)
        self.e.append(e)
        self.i.append(i)
        self.omega.append(omega)
        self.Omega.append(Omega)
        self.M.append(M)
        self.theta.append(self.trueAnomaly(M, e))
        # Updating Perturbation Parameters
        self.deltaN.append(deltaN)
        self.Omegadot.append(Omegadot)
        self.idot.append(idot)
        self.Cus.append(Cus)
        self.Cuc.append(Cuc)
        self.Cis.append(Cis)
        self.Cic.append(Cic)
        self.Crs.append(Crs)
        self.Crc.append(Crc)

    def setOfEphemeris(self, weekSec) -> float:
        """
        Out

        Parameters
        ----------
        tObs : TYPE
            DESCRIPTION.

        Returns
        -------
        arrayPos : float
            Middle position in the array to account for the interpolation.

        """
        for index, t in enumerate(self.toe):
            if (weekSec - t) > 0:
                arrayPos = index + 1/2
            elif weekSec <= self.toe[0]:
                arrayPos = 0 + 1/2
        return arrayPos

    def weekSeconds(self, t) -> float:
        """
        Calculate seconds of the current week from a given epoch. The reference
        epoch is 0 hours (midnight) Sunday 6-Jan-1980. Note that it is being
        dealt with UTC time to be consistent.

        Parameters
        ----------
        t : string
            Epoch of observation.

        Returns
        -------
        weekSeconds : float
            Seconds of the current week.

        """
        secInWeek = 604800
        GPStimestamp = 315964800
        epochTimestamp = self.getTimestamp(t)
        weekSeconds = (epochTimestamp - GPStimestamp) % secInWeek

        return weekSeconds

    def getTimestamp(self, t):
        """
        

        Parameters
        ----------
        t : TYPE
            DESCRIPTION.

        Returns
        -------
        timestamp : TYPE
            DESCRIPTION.

        """
        if type(t) is not str:
            t_str = t.strftime("%Y-%m-%d %H:%M:%S")
        else:
            t_str = t

        t_datetime = datetime.strptime(t_str + " +0000", "%Y-%m-%d %H:%M:%S %z")
        timestamp = datetime.timestamp(t_datetime)

        return timestamp

    def interpolation(self, t, y1, y2, t1, t2):
        y_t = y1 + (t - t1) * ((y2 - y1)/(t2 - t1))
 
        return y_t

    def getXYZ(self, t: float, P1: float):
        """
        Computing the satellite position in cartesian coordinates (i.e. ECEF
        coordinate system) from a given set of ephemeris data at a particular
        time of observation.

        Parameters
        ----------
        t : float
            Time of observation.
        P1 : float
            Pseudo-range, used to apply the pseudo-range correction.

        Returns
        -------
        x : TYPE
            x-coordinate [m].
        y : TYPE
            y-coordinate [m].
        z : TYPE
            z-coordinate [m].
        """
        if type(t) is pandas._libs.tslibs.timestamps.Timestamp:
            weekSec = self.weekSeconds(t)
            ephSet = self.setOfEphemeris(weekSec)
        else:
            weekSec = t
            ephSet = self.setOfEphemeris(weekSec)

        self.set.append(ephSet)
        toe = self.toe[math.floor(ephSet)]
        a = self.a[math.floor(ephSet)]
        e = self.e[math.floor(ephSet)]
        i = self.i[math.floor(ephSet)]
        omega = self.omega[math.floor(ephSet)]
        Omega = self.Omega[math.floor(ephSet)]
        M = self.M[math.floor(ephSet)]
        deltaN = self.deltaN[math.floor(ephSet)]
        Omegadot = self.Omegadot[math.floor(ephSet)]
        idot = self.idot[math.floor(ephSet)]
        Cus = self.Cus[math.floor(ephSet)]
        Cuc = self.Cuc[math.floor(ephSet)]
        Cis = self.Cis[math.floor(ephSet)]
        Cic = self.Cic[math.floor(ephSet)]
        Crs = self.Crs[math.floor(ephSet)]
        Crc = self.Crc[math.floor(ephSet)]

        if np.isnan(P1):
            print("No pseudorange correction available at ", t)
            t_k = weekSec - toe
        else:
            t_k = (weekSec - P1/c) - toe

        # t_k must account for beginning or end of week crossovers
        if t_k > 302400:
            t_k = t_k - 604800
        elif t_k < -302400:
            t_k = t_k + 604800

        n0 = math.sqrt(mu / a**3)
        n = n0 + deltaN
        M_corr = M + n * t_k
        E_corr = self.newton(M_corr, e)
        theta_corr = self.trueAnomaly(M_corr, e)
        phi = theta_corr + omega
        # Second Harmonic Perturbations
        delta_u = Cus * math.sin(2 * phi) + Cuc * math.cos(2 * phi)
        delta_r = Crs * math.sin(2 * phi) + Crc * math.cos(2 * phi)
        delta_i = Cis * math.sin(2 * phi) + Cic * math.cos(2 * phi)
        u = phi + delta_u
        r = a * (1 - e * math.cos(E_corr)) + delta_r
        i_corr = i + t_k * idot + delta_i
        # Position in orbital plane
        x_plane = r * math.cos(u)
        y_plane = r * math.sin(u)

        Omega_corr = Omega + (Omegadot - Omegadot_Earth) * t_k - \
            Omegadot_Earth * toe

        x = x_plane * math.cos(Omega_corr) \
            - y_plane * math.cos(i_corr) * math.sin(Omega_corr)

        y = x_plane * math.sin(Omega_corr) \
            + y_plane * math.cos(i_corr) * math.cos(Omega_corr)

        z = y_plane * math.sin(i_corr)

        return x, y, z

    def ecef2geo(self, x, y, z):
        tol = 1e-13
        # WGS-84 geodetic constants
        a = 6378137.000  # Semi-major axis
        f = 1/298.257223563  # Flattening
        b = a * (1 - f)  # 6356752.314 - semi-minor axis;
        esq = (a**2-b**2)/a**2  # Fist squared eccentricity, e^2
        # Calculate longitude
        lon = math.atan2(y, x)
        # Calculate initial latitude
        p = np.sqrt(x**2 + y**2)
        lat = math.atan2(z * (1 + esq), p)
        N = a/np.sqrt(1 - esq*math.sin(lat)**2)  # Radius of prime vertical

        # Calculate latitude and height
        h = 40  # initial guess of height
        h_old = 0  # initial 'old' height (used for first calc of dh)
        dh = 1  # initial dh value (arbitrary, but must be larger than 0.001)
        while dh > tol:
            if lat == np.pi/2:
                # https://askanydifference.com/difference-between-north-pole-and-south-pole/#:~:text=The%20north%20pole%20lies%20in%20the%20Arctic%20that%20is%20a,2%2C385m%20above%20sea%20level.
                h = z - 6371018
                print("Position at latitude 90º")
                break
            elif lat == - np.pi/2:
                h = abs(z) - 2385
                print("Position at latitude -90º")
                break
            lat = math.atan2(z + N * esq * math.sin(lat), p)
            h = p / math.cos(lat) - N
            # Check "precision" of iteration
            N = a/np.sqrt(1 - esq*math.sin(lat)**2)  # Radius of prime vertical
            dh = np.abs(h_old - h)  # Height difference
            h_old = h  # Store h for next iteration

        # Website to check
        # https://www.oc.nps.edu/oc2902w/coord/llhxyz.htm
        # eph[sat].ecef2geo(1586032.3754, -1932259.093, 5848547.0971)

        """
        Improved values of latitude and height are computed by iterating the
        equations B.6 in Page 198
        https://gssc.esa.int/navipedia/GNSS_Book/ESA_GNSS-Book_TM-23_Vol_I.pdf
        https://archive.psas.pdx.edu/CoordinateSystem/Latitude_to_LocalTangent.pdf
        https://gist.github.com/govert/1b373696c9a27ff4c72a
        https://stackoverflow.com/questions/56945401/converting-xyz-coordinates-to-longitutde-latitude-in-python
        """
        # print(math.degrees(lon), math.degrees(lat), h)
        return lat, lon, h

    def rot2enu(self, latitude, longitude):
        phi = pi_gps/2 - latitude
        lam = pi_gps/2 + longitude

        # Define rotation matices
        R1 = np.array([[1, 0, 0],
                      [0, math.cos(phi), math.sin(phi)],
                      [0, -math.sin(phi), math.cos(phi)]])
        R3 = np.array([[math.cos(lam), math.sin(lam), 0],
                       [-math.sin(lam), math.cos(lam), 0],
                       [0, 0, 1]])

        R = np.dot(R1, R3)

        return R

    def ecef2enu(self, x: float, y: float, z: float, receiverPos: list):
        """
        Conversion to East, North, Up coordinates of a single point given in
        ECEF (i.e. Earth-Centered Earth-Fixed) coordinates.

        Parameters
        ----------
        x : float
            x-coordinate [m].
        y : float
            y-coordinate [m].
        z : float
            z-coordinate [m].
        receiverPos : list
            List containing the position of the receiver [x, y, z] [m].

        Returns
        -------
        enu: array
            Array containing the East, North, Up coordinates.

        """
        satArray = np.array([x, y, z])
        recArray = np.array(receiverPos)
        sat_rec = np.subtract(satArray, recArray)
        lat, lon, h = self.ecef2geo(recArray[0], recArray[1], recArray[2])
        rotMat = self.rot2enu(lat, lon)
        enu = np.dot(rotMat, sat_rec.T)  # Rotate X,Y,Z coordinates to ENU

        return enu.T

    def enu2EAD(self, enu):
        """
        Elevation, azimuth and distance given the East, North, Up coordinates
        of a single point.

        Parameters
        ----------
        enu : TYPE
            DESCRIPTION.

        Returns
        -------
        elevation : float
            units: [rad].
        azimuth : float
            units: [rad].
        distance : float
            units: [m].
        """
        horizontalProj = np.sqrt(enu[0]**2 + enu[1]**2)
        elevation = np.arctan2(enu[2], horizontalProj)
        distance = np.sqrt(enu[0]**2 + enu[1]**2 + enu[2]**2)
        azimuth = np.arctan2(enu[0], enu[1])

        return elevation, azimuth, distance

    def getEAD(self, x: float, y: float, z: float, receiverPos: list):
        """
        Elevation, azimuth and distance given the cartesian coordinates of a 
        single point.

        Parameters
        ----------
        x : float
            x-coordinate [m].
        y : float
            y-coordinate [m].
        z : float
            z-coordinate [m].
        receiverPos : list
            List containing the position of the receiver [x, y, z] [m].

        Returns
        -------
        elevation : float
            units: [rad].
        azimuth : float
            units: [rad].
        distance : float
            units: [m].
        """
        recSat_ENU = self.ecef2enu(x, y, z, receiverPos)
        elevation, azimuth, distance = self.enu2EAD(recSat_ENU)

        return elevation, azimuth, distance
