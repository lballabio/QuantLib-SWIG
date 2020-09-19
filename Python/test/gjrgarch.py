"""
 Copyright (C) 2019 Pedro Coelho

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
"""

import QuantLib as ql
from math import pi, exp, sqrt
import unittest


class GJRGARCHEngineTest(unittest.TestCase):
    def setUp(self):
        self.settle_date = ql.Date.todaysDate()
        ql.Settings.instance().evaluationDate = self.settle_date
        self.dayCounter = ql.ActualActual()
        self.risk_free_handle = ql.YieldTermStructureHandle(ql.FlatForward(self.settle_date, 0.05, self.dayCounter))
        self.dividend_yield_handle = ql.YieldTermStructureHandle(ql.FlatForward(self.settle_date, 0, self.dayCounter))
        self.s0 = 50
        self.omega = 0.000002
        self.alpha = 0.024
        self.beta = 0.93
        self.gamma = 0.059
        self.daysPerYear = 365.0
        self.maturity = [90, 180]
        self.strike = [35, 40, 45, 50, 55, 60]
        self.lambda_values = [0.0, 0.1, 0.2]
        # correct values of analytic approximation
        self.analytic = [
            [[15.4315, 10.5552, 5.9625, 2.3282, 0.5408, 0.0835], [15.8969, 11.2173, 6.9112, 3.4788, 1.3769, 0.4357]],
            [[15.4556, 10.6929, 6.2381, 2.6831, 0.7822, 0.1738], [16.0587, 11.5338, 7.3170, 3.9074, 1.7279, 0.6568]],
            [[15.8000, 11.2734, 7.0376, 3.6767, 1.5871, 0.5934], [16.9286, 12.3170, 8.0405, 4.6348, 2.3429, 1.0590]],
        ]
        # correct values of Monte Carlo
        self.mc_values = [
            [[15.4332, 10.5453, 5.9351, 2.3521, 0.5597, 0.0776], [15.8910, 11.1772, 6.8827, 3.5096, 1.4196, 0.4502]],
            [[15.4580, 10.6433, 6.2019, 2.7513, 0.8374, 0.1706], [15.9884, 11.4139, 7.3103, 4.0497, 1.8862, 0.7322]],
            [[15.6619, 11.1263, 7.0968, 3.9152, 1.8133, 0.7010], [16.5195, 12.3181, 8.6085, 5.5700, 3.3103, 1.8053]],
        ]

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def testOptionPricing(self):

        tolerance = 0.075
        for k in range(3):
            lambda_value = self.lambda_values[k]
            m1 = (
                self.beta
                + (self.alpha + self.gamma * ql.CumulativeNormalDistribution()(lambda_value))
                * (1 + lambda_value * lambda_value)
                + self.gamma * lambda_value * exp(-lambda_value * lambda_value / 2) / sqrt(2 * pi)
            )
            v0 = self.omega / (1 - m1)
            quote = ql.QuoteHandle(ql.SimpleQuote(self.s0))
            garch = ql.GJRGARCHProcess(
                self.risk_free_handle,
                self.dividend_yield_handle,
                quote,
                v0,
                self.omega,
                self.alpha,
                self.beta,
                self.gamma,
                lambda_value,
                self.daysPerYear,
            )
            garch_model = ql.GJRGARCHModel(garch)
            analytic_engine = ql.AnalyticGJRGARCHEngine(garch_model)
            mc_engine = ql.MCEuropeanGJRGARCHEngine(
                process=garch, traits="pseudorandom", timeStepsPerYear=20, requiredTolerance=0.02, seed=1234
            )
            for i in range(2):
                for j in range(6):
                    payoff = ql.PlainVanillaPayoff(ql.Option.Call, self.strike[j])
                    ex_date = self.settle_date + ql.Period(self.maturity[i], ql.Days)
                    exercise = ql.EuropeanExercise(ex_date)
                    option = ql.VanillaOption(payoff, exercise)
                    option.setPricingEngine(analytic_engine)
                    analytic_price = option.NPV()
                    analytic_difference = analytic_price - self.analytic[k][i][j]
                    self.assertTrue(analytic_difference <= 2 * tolerance)
                    option.setPricingEngine(mc_engine)
                    mc_price = option.NPV()
                    mc_difference = mc_price - self.mc_values[k][i][j]
                    self.assertTrue(mc_difference <= 2 * tolerance)


class GJRGARCHCalibrationTest(unittest.TestCase):
    def setUp(self):
        self.settle_date = ql.Date(5, ql.July, 2002)
        ql.Settings.instance().evaluationDate = self.settle_date
        self.dayCounter = ql.Actual365Fixed()
        self.calendar = ql.TARGET()
        self.days = [0, 13, 41, 75, 165, 256, 345, 524, 703]
        self.rates = [0.0357, 0.0357, 0.0349, 0.0341, 0.0355, 0.0359, 0.0368, 0.0386, 0.0401]
        dates = list()
        for day in self.days:
            date = self.settle_date + ql.Period(day, ql.Days)
            dates.append(date)
        self.risk_free_ts = ql.YieldTermStructureHandle(ql.ZeroCurve(dates, self.rates, self.dayCounter))
        self.dividend_yield_handle = ql.YieldTermStructureHandle(ql.FlatForward(self.settle_date, 0, self.dayCounter))
        self.s0 = 4468.17
        self.omega = 0.000002
        self.alpha = 0.024
        self.beta = 0.93
        self.gamma = 0.059
        self.daysPerYear = 365.0
        self.Volatility = [
            0.6625,
            0.4875,
            0.4204,
            0.3667,
            0.3431,
            0.3267,
            0.3121,
            0.3121,
            0.6007,
            0.4543,
            0.3967,
            0.3511,
            0.3279,
            0.3154,
            0.2984,
            0.2921,
            0.5084,
            0.4221,
            0.3718,
            0.3327,
            0.3155,
            0.3027,
            0.2919,
            0.2889,
            0.4541,
            0.3869,
            0.3492,
            0.3149,
            0.2963,
            0.2926,
            0.2819,
            0.2800,
            0.4060,
            0.3607,
            0.3330,
            0.2999,
            0.2887,
            0.2811,
            0.2751,
            0.2775,
            0.3726,
            0.3396,
            0.3108,
            0.2781,
            0.2788,
            0.2722,
            0.2661,
            0.2686,
            0.3550,
            0.3277,
            0.3012,
            0.2781,
            0.2781,
            0.2661,
            0.2661,
            0.2681,
            0.3428,
            0.3209,
            0.2958,
            0.2740,
            0.2688,
            0.2627,
            0.2580,
            0.2620,
            0.3302,
            0.3062,
            0.2799,
            0.2631,
            0.2573,
            0.2533,
            0.2504,
            0.2544,
            0.3343,
            0.2959,
            0.2705,
            0.2540,
            0.2504,
            0.2464,
            0.2448,
            0.2462,
            0.3460,
            0.2845,
            0.2624,
            0.2463,
            0.2425,
            0.2385,
            0.2373,
            0.2422,
            0.3857,
            0.2860,
            0.2578,
            0.2399,
            0.2357,
            0.2327,
            0.2312,
            0.2351,
            0.3976,
            0.2860,
            0.2607,
            0.2356,
            0.2297,
            0.2268,
            0.2241,
            0.2320,
        ]
        self.strike = [3400, 3600, 3800, 4000, 4200, 4400, 4500, 4600, 4800, 5000, 5200, 5400, 5600]
        self.lambda_value = 0.1

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def testCalibration(self):

        m1 = (
            self.beta
            + (self.alpha + self.gamma * ql.CumulativeNormalDistribution()(self.lambda_value))
            * (1 + self.lambda_value * self.lambda_value)
            + self.gamma * self.lambda_value * exp(-self.lambda_value * self.lambda_value / 2) / sqrt(2 * pi)
        )
        v0 = self.omega / (1 - m1)

        helpers = list()
        for s in range(3, 10):
            for m in range(0, 3):
                vol = ql.QuoteHandle(ql.SimpleQuote(self.Volatility[s * 8 + m]))
                maturity = ql.Period(int((self.days[m + 1] + 3) / 7), ql.Weeks)
                heston_helper = ql.HestonModelHelper(
                    maturity,
                    self.calendar,
                    self.s0,
                    self.strike[s],
                    vol,
                    self.risk_free_ts,
                    self.dividend_yield_handle,
                    ql.BlackCalibrationHelper.ImpliedVolError,
                )
                helpers.append(heston_helper)

        new_garch_process = ql.GJRGARCHProcess(
            self.risk_free_ts,
            self.dividend_yield_handle,
            ql.QuoteHandle(ql.SimpleQuote(self.s0)),
            v0,
            self.omega,
            self.alpha,
            self.beta,
            self.gamma,
            self.lambda_value,
            self.daysPerYear,
        )
        new_garch_model = ql.GJRGARCHModel(new_garch_process)
        new_garch_engine = ql.AnalyticGJRGARCHEngine(new_garch_model)
        for helper in helpers:
            helper.setPricingEngine(new_garch_engine)

        om = ql.Simplex(0.05)
        new_garch_model.calibrate(helpers, om, ql.EndCriteria(400, 40, 1.0e-8, 1.0e-8, 1.0e-8))

        sse = 0
        for helper in helpers:
            diff = helper.calibrationError() * 100
            sse += diff * diff

        maxExpected = 15
        self.assertTrue(sse <= maxExpected)


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(GJRGARCHEngineTest, "test"))
    suite.addTest(unittest.makeSuite(GJRGARCHCalibrationTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
