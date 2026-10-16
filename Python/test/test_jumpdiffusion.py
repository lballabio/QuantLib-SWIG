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

import math
import unittest

import QuantLib as ql


class JumpDiffusionTest(unittest.TestCase):
    def setUp(self):
        self.today = ql.Date(15, ql.October, 2026)
        ql.Settings.instance().evaluationDate = self.today
        self.day_counter = ql.Actual360()

        self.spot = ql.SimpleQuote(100.0)
        self.vol = ql.SimpleQuote(0.0)
        self.jump_intensity = ql.SimpleQuote(0.0)
        self.mean_log_jump = ql.SimpleQuote(0.0)
        self.jump_vol = ql.SimpleQuote(0.0)

        self.process = ql.Merton76Process(
            ql.QuoteHandle(self.spot),
            ql.YieldTermStructureHandle(ql.FlatForward(self.today, 0.0, self.day_counter)),
            ql.YieldTermStructureHandle(ql.FlatForward(self.today, 0.08, self.day_counter)),
            ql.BlackVolTermStructureHandle(
                ql.BlackConstantVol(self.today, ql.NullCalendar(),
                                    ql.QuoteHandle(self.vol), self.day_counter)),
            ql.QuoteHandle(self.jump_intensity),
            ql.QuoteHandle(self.mean_log_jump),
            ql.QuoteHandle(self.jump_vol),
        )

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def make_option(self, strike, t):
        exercise = ql.EuropeanExercise(self.today + int(round(t * 360)))
        return ql.VanillaOption(ql.PlainVanillaPayoff(ql.Option.Call, strike), exercise)

    def testMerton76Values(self):
        """Testing the Merton 76 jump-diffusion engine against literature values"""
        # From "Option pricing formulas", E.G. Haug, McGraw-Hill 1998, p. 9,
        # as corrected in the C++ test suite; total volatility 25%, a share gamma
        # of the variance coming from jumps, and zero mean jump.
        values = [
            # strike, t, jump intensity, gamma, value
            (80.0, 0.10, 1.0, 0.25, 20.67),
            (90.0, 0.25, 5.0, 0.25, 12.75),
            (100.0, 0.50, 10.0, 0.25, 9.03),
            (110.0, 0.25, 1.0, 0.25, 2.11),
            (120.0, 0.50, 5.0, 0.25, 2.25),
            (100.0, 0.10, 1.0, 0.50, 3.14),
            (90.0, 0.50, 10.0, 0.50, 15.41),
        ]
        total_vol = 0.25
        engine = ql.JumpDiffusionEngine(self.process)

        for strike, t, intensity, gamma, expected in values:
            jump_vol = total_vol * math.sqrt(gamma / intensity)
            self.jump_intensity.setValue(intensity)
            self.jump_vol.setValue(jump_vol)
            self.mean_log_jump.setValue(-0.5 * jump_vol * jump_vol)
            self.vol.setValue(total_vol * math.sqrt(1.0 - gamma))

            option = self.make_option(strike, t)
            option.setPricingEngine(engine)
            self.assertAlmostEqual(
                option.NPV(), expected, delta=1.0e-2,
                msg=f"strike {strike}, t {t}, intensity {intensity}, gamma {gamma}")

    def testNoJumpsMatchesBlackScholes(self):
        """Testing that the jump-diffusion engine without jumps gives the Black-Scholes value"""
        self.vol.setValue(0.25)
        option = self.make_option(100.0, 0.5)

        option.setPricingEngine(
            ql.JumpDiffusionEngine(self.process, relativeAccuracy=1.0e-8, maxIterations=200))
        calculated = option.NPV()

        black_scholes_process = ql.BlackScholesMertonProcess(
            ql.QuoteHandle(self.spot),
            ql.YieldTermStructureHandle(ql.FlatForward(self.today, 0.0, self.day_counter)),
            ql.YieldTermStructureHandle(ql.FlatForward(self.today, 0.08, self.day_counter)),
            ql.BlackVolTermStructureHandle(
                ql.BlackConstantVol(self.today, ql.NullCalendar(),
                                    ql.QuoteHandle(self.vol), self.day_counter)),
        )
        option.setPricingEngine(ql.AnalyticEuropeanEngine(black_scholes_process))
        expected = option.NPV()

        self.assertAlmostEqual(calculated, expected, delta=1.0e-10)


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
