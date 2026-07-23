"""
 Copyright (C) 2026 QuantLib contributors

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


class MultiCurveTest(unittest.TestCase):

    def test_piecewise_and_spreaded_curve(self):
        """Globally bootstrap a dependency cycle of two yield curves."""
        today = ql.Date(23, ql.October, 2025)
        ql.Settings.instance().evaluationDate = today

        internal_ois = ql.RelinkableYieldTermStructureHandle()
        internal_3m = ql.RelinkableYieldTermStructureHandle()
        euribor_3m = ql.Euribor3M(internal_3m)

        quote = ql.SimpleQuote(0.03)
        spread = ql.SimpleQuote(-0.01)
        helpers = []
        for maturity in range(1, 11):
            helpers.append(
                ql.SwapRateHelper(
                    ql.QuoteHandle(quote),
                    ql.Period(maturity, ql.Years),
                    euribor_3m.fixingCalendar(),
                    ql.Annual,
                    ql.Following,
                    ql.Thirty360(ql.Thirty360.BondBasis),
                    euribor_3m,
                    ql.QuoteHandle(),
                    ql.Period(0, ql.Days),
                    internal_ois,
                )
            )

        accuracy = 1.0e-10
        multi_curve = ql.MultiCurve(accuracy)
        curve_3m = ql.GlobalLinearSimpleZeroCurve(
            today, helpers, ql.Actual360(), ql.GlobalBootstrap(accuracy)
        )
        external_3m = multi_curve.addBootstrappedCurve(internal_3m, curve_3m)

        curve_ois = ql.ZeroSpreadedTermStructure(
            internal_3m, ql.QuoteHandle(spread)
        )
        external_ois = multi_curve.addNonBootstrappedCurve(internal_ois, curve_ois)

        # The external handles keep the MultiCurve and all its members alive.
        del multi_curve, curve_3m, curve_ois

        forecast = external_3m.currentLink()
        discount = external_ois.currentLink()
        self.assertAlmostEqual(
            discount.zeroRate(1.0, ql.Continuous).rate()
            - forecast.zeroRate(1.0, ql.Continuous).rate(),
            spread.value(),
            delta=accuracy,
        )
        for helper in helpers:
            self.assertAlmostEqual(helper.quoteError(), 0.0, delta=accuracy)


if __name__ == "__main__":
    unittest.main()
