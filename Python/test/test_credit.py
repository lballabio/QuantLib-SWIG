"""
 Copyright (C) 2026 Ryan Duguid

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

import unittest

import QuantLib as ql


class CreditTest(unittest.TestCase):
    def setUp(self):
        self.todaysDate = ql.Date(15, ql.May, 2019)
        ql.Settings.instance().evaluationDate = self.todaysDate

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def testUpfrontCdsHelperDefaultSettlementDays(self):
        """Testing default upfront settlement days of UpfrontCdsHelper"""

        dc = ql.Actual365Fixed()
        discountCurve = ql.YieldTermStructureHandle(
            ql.FlatForward(self.todaysDate, 0.05, dc))

        def hazardRate(**kwargs):
            helper = ql.UpfrontCdsHelper(
                0.02, 0.01, ql.Period(5, ql.Years), 0, ql.WeekendsOnly(),
                ql.Quarterly, ql.Following, ql.DateGeneration.CDS,
                ql.Actual360(), 0.4, discountCurve, **kwargs)
            curve = ql.PiecewiseFlatHazardRate(self.todaysDate, [helper], dc)
            return curve.hazardRate(helper.maturityDate())

        self.assertAlmostEqual(
            hazardRate(), hazardRate(upfrontSettlementDays=3), 12)


if __name__ == '__main__':
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
