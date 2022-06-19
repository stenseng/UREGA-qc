# -*- coding: utf-8 -*-
"""
Created on Fri May 20 15:35:20 2022

@author: manue
"""

import unittest
import numpy as np
import math

from ephemeris import Eph
eph = Eph()


class TestEph(unittest.TestCase):
    def test_timeStamp(self):
        """
        https://www.unixtimestamp.com/
        """
        GPS = "1980-01-06 00:00:00"
        manu_sBirthday = "1999-10-12 00:00:00"
        self.assertEqual(eph.getTimestamp(GPS), 315964800)
        self.assertEqual(eph.getTimestamp(manu_sBirthday), 939686400)

    def test_weekSeconds(self):
        """
        Reference epoch test date 2021-11-03 at midnight, correct result known
        from: https://www.labsat.co.uk/index.php/en/gps-time-calculator
        """
        epochRef = "2021-11-03 00:00:00"
        GPS = "1980-01-06 00:00:00"
        lowLim = "2022-05-15 00:00:00"
        upLim = "2022-06-11 23:59:59"
        rightNow = "2022-06-08 11:25:35"
        ISMR = "2021-11-03 00:21:00"
        self.assertEqual(eph.weekSeconds(epochRef), 259200)
        self.assertEqual(eph.weekSeconds(GPS), 0)
        self.assertEqual(eph.weekSeconds(lowLim), 0)
        self.assertEqual(eph.weekSeconds(upLim), 604799)
        self.assertEqual(eph.weekSeconds(rightNow), 300335)
        self.assertEqual(eph.weekSeconds(ISMR), 260460)

    def test_trueAnomaly(self):
        """
        Test that the function returns the correct true anomaly and compare with
        analytic expression.
        """
        # Case of 0 rad
        self.assertEqual(eph.trueAnomaly(0, 0.1), 0)
        # Case of pi rad
        self.assertEqual(eph.trueAnomaly(np.pi, 0.1), np.pi)
        # Case of a circle and arbitrary angle
        self.assertEqual(eph.trueAnomaly(np.pi/4, 0), np.pi/4)
        # Case of a circle and pi/2
        self.assertEqual(eph.trueAnomaly(np.pi/2, 0), np.pi/2)

        M = 2.65
        e = 0.6
        E = eph.newton(M, e)
        theta = 2 * math.atan(math.sqrt((1 + e) /
                                        (1 - e)) * math.tan(E / 2))
        lim = 16
        self.assertEqual(round(eph.trueAnomaly(M, e), lim), round(theta, lim))

    def test_ecef2geo(self):
        """
        
        """
        # print(eph.ecef2geo(0, 0, 6371588))
        # print(eph.ecef2geo(1586032.3754, -1932259.093, 5848547.0971))
    def test_ecef2ead(self):
        northPole = np.array([0, 0, 6371588])
        southPole = np.array([0, 0, -6371588])
        sat1 = np.array([0, 0, 6381588])
        sat1 = np.array([1000, 1000, 6371588])
        sat1 = np.array([0, 0, -6381588])
        recSat_ENU = eph.ecef2enu(sat1[0], sat1[1], sat1[2], southPole)
        print(eph.enu2EAD(recSat_ENU))
    
    def test_getEAD(self):
        """
        x, y, z = eph.getXYZ(259260.0, P1=0)
        receiver_pos = [1586032.3754, -1932259.093, 5848547.0971]
        print(eph.getEAD(x, y, z, receiver_pos))
"""

if __name__ == '__main__':
    unittest.main()
