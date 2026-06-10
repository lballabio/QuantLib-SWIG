# coding=utf-8-unix
"""
 Copyright (C) 2018 Wojciech Ślusarski

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
"""

import QuantLib as ql
import unittest


class IborIndexTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.euribor3m = ql.Euribor3M()

    def setUp(self):
        self.euribor3m.clearFixings()
        # values are not real due to copyrights of the fixing
        self.euribor3m.addFixing(ql.Date(17, 7, 2018), -0.3)
        self.euribor3m.addFixings([ql.Date(12, 7, 2018), ql.Date(13, 7, 2018)], [-0.3, -0.3])

    def testAddFixingFail(self):
        """Testing for RuntimeError while trying to overwrite fixing value"""

        with self.assertRaises(RuntimeError):
            # attempt to overwrite value that is already set at different level
            self.euribor3m.addFixing(ql.Date(17, 7, 2018), -0.4)

        with self.assertRaises(RuntimeError):
            # attempt to overwrite value that is already set at different level
            self.euribor3m.addFixings([ql.Date(12, 7, 2018), ql.Date(13, 7, 2018)], [-0.4, -0.4])

    def testAddFixing(self):
        """Testing for overwriting fixing value"""

        force_overwrite = True
        try:
            # attempt to overwrite value that is already set at different level
            self.euribor3m.addFixing(ql.Date(17, 7, 2018), -0.4, force_overwrite)
            self.euribor3m.addFixings([ql.Date(12, 7, 2018), ql.Date(13, 7, 2018)], [-0.4, -0.4], force_overwrite)
            # try clearFixings and repeat with original levels
            self.euribor3m.clearFixings()
            self.euribor3m.addFixing(ql.Date(17, 7, 2018), -0.3)
            self.euribor3m.addFixings([ql.Date(12, 7, 2018), ql.Date(13, 7, 2018)], [-0.3, -0.3])

        except RuntimeError as err:
            raise AssertionError("Failed to overwrite index fixixng " + "{}".format(err))

    def testTimeSeries(self):
        """Testing for getting time series of the fixing"""

        dates = (ql.Date(12, 7, 2018), ql.Date(13, 7, 2018), ql.Date(17, 7, 2018))
        values = (-0.3, -0.3, -0.3)
        for expected, actual in zip(dates, self.euribor3m.timeSeries().dates()):
            self.assertTrue(expected == actual)
        for expected, actual in zip(values, self.euribor3m.timeSeries().values()):
            self.assertTrue(expected == actual)


class BMAIndexTest(unittest.TestCase):
    def setUp(self):
        self.index = ql.BMAIndex()

    def testInspectors(self):
        """Testing BMA index inspectors"""

        self.assertEqual(self.index.familyName(), "BMA")
        self.assertEqual(self.index.tenor(), ql.Period(1, ql.Weeks))
        self.assertEqual(self.index.fixingDays(), 1)
        self.assertEqual(self.index.currency(), ql.USDCurrency())

    def testFixingDates(self):
        """Testing BMA fixing-date logic"""

        self.assertTrue(self.index.isValidFixingDate(ql.Date(4, 1, 2023)))
        self.assertFalse(self.index.isValidFixingDate(ql.Date(5, 1, 2023)))

    def testMaturityDate(self):
        """Testing BMA maturity-date calculation"""

        self.assertEqual(
            self.index.maturityDate(ql.Date(5, 1, 2023)),
            ql.Date(12, 1, 2023))

    def testFixingSchedule(self):
        """Testing BMA fixing schedule"""

        schedule = self.index.fixingSchedule(
            ql.Date(5, 1, 2023), ql.Date(20, 1, 2023))
        self.assertEqual(schedule[0], ql.Date(4, 1, 2023))
        self.assertEqual(schedule[-1], ql.Date(25, 1, 2023))


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
