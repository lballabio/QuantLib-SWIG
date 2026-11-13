"""
 Copyright (C) 2026 Qin Dai

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

import math
import unittest

import QuantLib as ql


class FittedBondCurveTest(unittest.TestCase):
    def setUp(self):
        self.today = ql.Date(16, ql.September, 2015)
        ql.Settings.instance().evaluationDate = self.today
        self.day_counter = ql.Actual365Fixed()

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def testNaturalCubicFittingEvaluation(self):
        """Testing natural cubic fitting method as evaluator"""
        knots = [1.0, 2.0, 5.0, 10.0]
        discounts = [math.exp(-0.03 * t) for t in knots]

        method = ql.NaturalCubicFitting(knots)
        curve = ql.FittedBondDiscountCurve(
            self.today, method, discounts, self.today + ql.Period(12, ql.Years), self.day_counter
        )

        self.assertEqual(curve.fitResults().size(), len(knots))

        # the curve starts at 1, passes through the knots...
        self.assertAlmostEqual(curve.discount(0.0), 1.0, delta=1e-12)
        for t, d in zip(knots, discounts):
            self.assertAlmostEqual(curve.discount(t), d, delta=1e-12)

        # ...interpolates between them with a natural cubic spline...
        spline = ql.CubicNaturalSpline([0.0] + knots, [1.0] + discounts)
        for t in [0.25, 0.5, 1.5, 3.0, 7.5, 9.9]:
            self.assertAlmostEqual(curve.discount(t), spline(t), delta=1e-12)

        # ...and is flat after the last one
        self.assertAlmostEqual(curve.discount(11.0), discounts[-1], delta=1e-12)

    def testNaturalCubicFittingCalibration(self):
        """Testing natural cubic fitting method on bonds priced off a known curve"""
        flat_curve = ql.FlatForward(self.today, 0.03, self.day_counter, ql.Continuous)
        flat_engine = ql.DiscountingBondEngine(ql.YieldTermStructureHandle(flat_curve))

        bonds = []
        helpers = []
        for years in range(1, 11):
            schedule = ql.Schedule(
                self.today,
                self.today + ql.Period(years, ql.Years),
                ql.Period(ql.Annual),
                ql.NullCalendar(),
                ql.Unadjusted,
                ql.Unadjusted,
                ql.DateGeneration.Backward,
                False,
            )
            bond = ql.FixedRateBond(0, 100.0, schedule, [0.03], self.day_counter)
            bond.setPricingEngine(flat_engine)
            helpers.append(ql.BondHelper(ql.QuoteHandle(ql.SimpleQuote(bond.cleanPrice())), bond))
            bonds.append(bond)

        method = ql.NaturalCubicFitting(knotTimes=[1.0, 2.0, 3.0, 5.0, 7.0, 10.0])
        fitted_curve = ql.FittedBondDiscountCurve(self.today, helpers, self.day_counter, method)

        # the fitted curve is close to the one used to price the bonds...
        for i in range(1, 21):
            t = 0.5 * i
            self.assertAlmostEqual(fitted_curve.discount(t), flat_curve.discount(t), delta=1e-3)

        # ...and reprices them
        fitted_engine = ql.DiscountingBondEngine(ql.YieldTermStructureHandle(fitted_curve))
        for bond in bonds:
            quote = bond.cleanPrice()
            bond.setPricingEngine(fitted_engine)
            self.assertAlmostEqual(bond.cleanPrice(), quote, delta=5e-2)


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
