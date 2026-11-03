"""
 Copyright (C) 2026 Himanshu Gupta

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


class BlackVasicekEngineTest(unittest.TestCase):
    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def testDeterministicRates(self):
        """Testing that the Black-Vasicek engine gives the Black-Scholes value with deterministic rates"""
        today = ql.Date(15, ql.October, 2026)
        ql.Settings.instance().evaluationDate = today
        day_counter = ql.Actual365Fixed()
        r = 0.05

        process = ql.BlackScholesMertonProcess(
            ql.QuoteHandle(ql.SimpleQuote(100.0)),
            ql.YieldTermStructureHandle(ql.FlatForward(today, 0.02, day_counter)),
            ql.YieldTermStructureHandle(ql.FlatForward(today, r, day_counter)),
            ql.BlackVolTermStructureHandle(
                ql.BlackConstantVol(today, ql.NullCalendar(), 0.20, day_counter)),
        )
        # with a vanishing short-rate volatility and a long-term mean equal
        # to the initial rate, the Vasicek rate stays constant
        vasicek = ql.Vasicek(r, 0.3, r, 1.0e-10)

        option = ql.VanillaOption(
            ql.PlainVanillaPayoff(ql.Option.Call, 100.0),
            ql.EuropeanExercise(today + ql.Period(1, ql.Years)),
        )
        option.setPricingEngine(ql.AnalyticEuropeanEngine(process))
        expected = option.NPV()
        option.setPricingEngine(ql.AnalyticBlackVasicekEngine(process, vasicek, 0.0))

        self.assertAlmostEqual(option.NPV(), expected, delta=1.0e-8)


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
