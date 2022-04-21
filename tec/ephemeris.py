# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 13:09:57 2022

@author: manue
"""
import numpy as np
import math


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
                                        (1 - e)) *
                              math.tan(E / 2))
        return theta

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

    def getEphXYZ(
        self,
        t_obs,
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
        x = P1*t_obs

        y = 1.0

        z = 2.0
        return (x, y, z)

# (x, y, z) = eph["G15"].getXYZ(34592.2)
