# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 13:09:57 2022

@author: manue
"""
import numpy as np
import math
import conversion as ccn

c = 299792458.0  # [m/s]
mu = 3.986005e14  # [m^3/s^2]
Omegadot_Earth = 7.2921151467e-5  # [rad/s]


class Eph:
    def __init__(self):
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
        self.x = []
        self.y = []
        self.z = []
        self.elevation = []

    def __newton(self, f, Df, M, e):
        x0 = 0.01
        tol = 1e-13
        maxIter = 100
        xn = x0
        for n in range(0, maxIter):
            fxn = f(xn, M, e)
            if abs(fxn) < tol:
                return xn
            Dfxn = Df(xn, e)
            if Dfxn == 0:
                print("Zero derivative. No solution found.")
                return None
            xn = xn - fxn / Dfxn
        print("Exceeded maximum iterations. No solution found.")
        return None

    def __MtoE(self, E, M, e):
        f = E - M - e * np.sin(E)
        return f

    def __dfMtoE(self, E, e):
        df = 1 - e * np.cos(E)
        return df

    def __trueAnomaly(self, __newton, __MtoE, __dfMtoE, M, e):
        E = self.__newton(self.__MtoE, self.__dfMtoE, M, e)
        theta = 2 * math.atan(math.sqrt((1 + e) /
                                        (1 - e)) * math.tan(E / 2))
        return theta

    def weekSeconds(self, t):
        leapSeconds = 18  # As of 2022
        secInWeek = 604800
        epoch_ref_GPS = np.datetime64('1980-01-06T00:00:00')
        epoch = np.datetime64(t)
        time_diff = epoch - epoch_ref_GPS + np.timedelta64(leapSeconds, 's')
        # Given in [ns]
        weekSeconds = float(time_diff*1e-9) % secInWeek
        return weekSeconds

    def addNav(
        self,
        toe,
        toc,
        a0,
        a1,
        a2,
        a,
        e,
        i,
        omega,
        Omega,
        M,
        deltaN,
        Omegadot,
        idot,
        Cus,
        Cuc,
        Cis,
        Cic,
        Crs,
        Crc
    ):
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
        self.theta.append(self.__trueAnomaly(self.__newton,
                                             self.__MtoE,
                                             self.__dfMtoE,
                                             M,
                                             e)
                          )
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

    def getXYZ(
        self,
        t,
        P1,
        toe,
        toc,
        a0,
        a1,
        a2,
        a,
        e,
        i,
        omega,
        Omega,
        M,
        deltaN,
        Omegadot,
        idot,
        Cus,
        Cuc,
        Cis,
        Cic,
        Crs,
        Crc
    ):
        # for t in self.toe:
        t_week = self.weekSeconds(t)
        t_k = t_week - toe
        # - P1/c
        n0 = math.sqrt(mu / a**3)
        n = n0 + deltaN
        M_corr = M + n * t_k
        E_corr = self.__newton(self.__MtoE,
                               self.__dfMtoE,
                               M_corr,
                               e)
        theta_corr = self.__trueAnomaly(self.__newton,
                                        self.__MtoE,
                                        self.__dfMtoE,
                                        M_corr,
                                        e)
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

        x = (x_plane * math.cos(Omega_corr)
             - y_plane * math.cos(i_corr) * math.sin(Omega_corr))

        y = (x_plane * math.sin(Omega_corr)
             + y_plane * math.cos(i_corr) * math.cos(Omega_corr))

        z = (y_plane * math.sin(i_corr))

        return (x, y, z)

    def getElevation(
        self,
        x,
        y,
        z,
        receiver_pos
    ):
        rec_sat_enu_ccn = ccn.cart2enu(np.array([x, y, z]).T,
                                       np.array(receiver_pos))
        ead = ccn.enu2ead(rec_sat_enu_ccn)
        elevation = ead[:, 0] / np.pi * 180
        return elevation


# (x, y, z) = eph["G15"].getXYZ(34592.2)
