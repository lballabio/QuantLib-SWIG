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


class HestonExpansionTest(unittest.TestCase):
    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def make_model(self, today, r, q, v0, kappa, theta, sigma, rho, spot=100.0):
        day_counter = ql.Actual365Fixed()
        process = ql.HestonProcess(
            ql.YieldTermStructureHandle(ql.FlatForward(today, r, day_counter)),
            ql.YieldTermStructureHandle(ql.FlatForward(today, q, day_counter)),
            ql.QuoteHandle(ql.SimpleQuote(spot)),
            v0, kappa, theta, sigma, rho,
        )
        return ql.HestonModel(process)

    def testExpansionEngineOnAlanLewisReference(self):
        """Testing the Heston expansion engine on Alan Lewis reference prices"""
        # same data and tolerances as testExpansionOnAlanLewisReference
        # in the C++ test suite; Forde is not tested there either, as it
        # doesn't behave well on this example
        today = ql.Date(5, ql.July, 2002)
        ql.Settings.instance().evaluationDate = today
        exercise = ql.EuropeanExercise(ql.Date(5, ql.July, 2003))
        model = self.make_model(today, 0.01, 0.02, v0=0.04, kappa=4.0,
                                theta=0.25, sigma=1.0, rho=-0.5)

        strikes = [80, 90, 100, 110, 120]
        expected = [
            (7.958878113256768, 26.774758743998854),
            (12.017966707346305, 20.933349000596710),
            (17.055270961270109, 16.070154917028834),
            (23.017825898442801, 12.132211516709845),
            (29.811026202682472, 9.024913483457836),
        ]
        formulas = [
            (ql.HestonExpansionEngine.LPP2, 1.003e-2),
            (ql.HestonExpansionEngine.LPP3, 3.645e-3),
        ]

        for formula, tolerance in formulas:
            engine = ql.HestonExpansionEngine(model, formula)
            for strike, (put, call) in zip(strikes, expected):
                for option_type, value in [(ql.Option.Put, put), (ql.Option.Call, call)]:
                    option = ql.VanillaOption(ql.PlainVanillaPayoff(option_type, strike), exercise)
                    option.setPricingEngine(engine)
                    self.assertLessEqual(
                        abs(option.NPV() - value) / value, tolerance,
                        msg=f"formula {formula}, strike {strike}, type {option_type}")

    def testExpansionsOnFordeReference(self):
        """Testing the Heston expansions on Forde reference volatilities"""
        # same data and tolerances as testExpansionOnFordeReference
        # in the C++ test suite
        forward = 100.0
        v0, rho, sigma, kappa, theta = 0.04, -0.4, 0.2, 1.15, 0.04
        terms = [0.1, 1.0, 5.0, 10.0]
        strikes = [60, 80, 90, 100, 110, 120, 140]
        reference_vols = [
            [0.27284673574924445, 0.22360758200372477, 0.21023988547031242, 0.1990674789471587,
             0.19118230678920461, 0.18721342919371017, 0.1899869903378507],
            [0.25200775151345, 0.2127275920953156, 0.20286528150874591, 0.19479398358151515,
             0.18872591728967686, 0.18470857955411824, 0.18204457060905446],
            [0.21637821506229973, 0.20077227130455172, 0.19721753043236154, 0.1942233023784151,
             0.191693211401571, 0.18955229722896752, 0.18491727548069495],
            [0.20672925973965342, 0.198583062164427, 0.19668274423922746, 0.1950420231354201,
             0.193610364344706, 0.1923502827886502, 0.18934360917857015],
        ]
        # rows: LPP2, LPP3, Forde; columns: terms (Forde breaks down for long maturities)
        tolerances = [[0.06, 0.03, 0.03, 0.02], [0.15, 0.08, 0.04, 0.02], [0.06, 0.08, 1.0, 1.0]]
        atm_tolerances = [[4e-6, 7e-4, 2e-3, 9e-4], [7e-6, 4e-4, 9e-4, 4e-4], [4e-4, 3e-2, 0.28, 1.0]]

        for j, term in enumerate(terms):
            expansions = [
                ql.LPP2HestonExpansion(kappa, theta, sigma, v0, rho, term),
                ql.LPP3HestonExpansion(kappa, theta, sigma, v0, rho, term),
                ql.FordeHestonExpansion(kappa, theta, sigma, v0, rho, term),
            ]
            for i, strike in enumerate(strikes):
                for k, expansion in enumerate(expansions):
                    expected = reference_vols[j][i]
                    calculated = expansion.impliedVolatility(strike, forward)
                    tolerance = atm_tolerances[k][j] if strike == forward else tolerances[k][j]
                    self.assertLessEqual(
                        abs(calculated - expected) / expected, tolerance,
                        msg=f"expansion {k}, term {term}, strike {strike}")

    def testEngineMatchesExpansion(self):
        """Testing that the Heston expansion engine prices at the expansion's implied volatility"""
        today = ql.Date(15, ql.October, 2026)
        ql.Settings.instance().evaluationDate = today
        r, q, spot = 0.03, 0.01, 100.0
        v0, kappa, theta, sigma, rho = 0.04, 1.5, 0.05, 0.3, -0.6
        model = self.make_model(today, r, q, v0, kappa, theta, sigma, rho, spot)
        maturity = today + ql.Period(6, ql.Months)
        t = ql.Actual365Fixed().yearFraction(today, maturity)
        exercise = ql.EuropeanExercise(maturity)
        df = math.exp(-r * t)
        forward = spot * math.exp((r - q) * t)

        formulas = [
            (ql.HestonExpansionEngine.LPP2, ql.LPP2HestonExpansion),
            (ql.HestonExpansionEngine.LPP3, ql.LPP3HestonExpansion),
            (ql.HestonExpansionEngine.Forde, ql.FordeHestonExpansion),
        ]
        for formula, expansion_class in formulas:
            expansion = expansion_class(kappa, theta, sigma, v0, rho, t)
            engine = ql.HestonExpansionEngine(model, formula)
            for strike in [80.0, 100.0, 120.0]:
                for option_type in [ql.Option.Put, ql.Option.Call]:
                    option = ql.VanillaOption(ql.PlainVanillaPayoff(option_type, strike), exercise)
                    option.setPricingEngine(engine)
                    vol = expansion.impliedVolatility(strike, forward)
                    expected = ql.blackFormula(option_type, strike, forward, vol * math.sqrt(t), df)
                    self.assertAlmostEqual(
                        option.NPV(), expected, delta=1.0e-10,
                        msg=f"formula {formula}, strike {strike}, type {option_type}")

    def testCloseToAnalyticEngineNearTheMoney(self):
        """Testing the Heston expansion engine against the semi-analytic engine near the money"""
        today = ql.Date(15, ql.October, 2026)
        ql.Settings.instance().evaluationDate = today
        r, q, spot = 0.03, 0.01, 100.0
        model = self.make_model(today, r, q, v0=0.04, kappa=1.5,
                                theta=0.04, sigma=0.3, rho=-0.6, spot=spot)
        maturity = today + ql.Period(3, ql.Months)
        t = ql.Actual365Fixed().yearFraction(today, maturity)
        exercise = ql.EuropeanExercise(maturity)
        df = math.exp(-r * t)
        forward = spot * math.exp((r - q) * t)
        reference_engine = ql.AnalyticHestonEngine(model)

        def implied_vol(option, engine, option_type, strike):
            option.setPricingEngine(engine)
            std_dev = ql.blackFormulaImpliedStdDev(option_type, strike, forward, option.NPV(), df)
            return std_dev / math.sqrt(t)

        # The expansions are approximations, compared in implied volatility.
        # Here LPP2 and LPP3 are within 0.04 vol points of the semi-analytic
        # engine and Forde within 0.25; Forde degrades faster with maturity.
        formulas = [
            (ql.HestonExpansionEngine.LPP2, 1.0e-3),
            (ql.HestonExpansionEngine.LPP3, 1.0e-3),
            (ql.HestonExpansionEngine.Forde, 3.0e-3),
        ]
        for formula, tolerance in formulas:
            engine = ql.HestonExpansionEngine(model, formula)
            for strike in [90.0, 95.0, 100.0, 105.0, 110.0]:
                option_type = ql.Option.Call if strike >= forward else ql.Option.Put
                option = ql.VanillaOption(ql.PlainVanillaPayoff(option_type, strike), exercise)
                expected = implied_vol(option, reference_engine, option_type, strike)
                calculated = implied_vol(option, engine, option_type, strike)
                self.assertAlmostEqual(
                    calculated, expected, delta=tolerance,
                    msg=f"formula {formula}, strike {strike}")


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
