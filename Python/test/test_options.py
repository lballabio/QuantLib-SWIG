"""
 Copyright (C) 2021 Klaus Spanderen

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

class OptionsTest(unittest.TestCase):

    def testFdHestonHullWhite(self):
        """ Testing FDM Heston Hull-White pricing """

        dc = ql.Actual365Fixed()
        todays_date = ql.Settings.instance().evaluationDate
        
        r = ql.YieldTermStructureHandle(ql.FlatForward(todays_date, 0.075, dc))
        d = ql.YieldTermStructureHandle(ql.FlatForward(todays_date, 0.01, dc))

        s0 = 8.0

        v0 = 0.2*0.2
        kappa = 1.0
        theta = v0
        sigma = 0.4
        rho = -0.75

        a = 0.00883
        sig = 0.00631

        underlying = ql.makeQuoteHandle(s0)

        option = ql.VanillaOption(
            ql.PlainVanillaPayoff(ql.Option.Call, s0),
            ql.EuropeanExercise(todays_date + ql.Period(365, ql.Days))
        )

        hull_white_process = ql.HullWhiteProcess(r, a, sig)
        heston_process = ql.HestonProcess(r, d, underlying, v0, kappa, theta, sigma, rho)

        option.setPricingEngine(
            ql.FdHestonHullWhiteVanillaEngine(
                ql.HestonModel(heston_process), hull_white_process, -0.5,
                10, 200, 25, 10, 0, True
            )
        )

        self.assertAlmostEqual(0.87628, option.NPV(), 4)

    def testAnalyticHestonHullWhite(self):
        """ Testing Analytic Heston Hull-White pricing """

        today = ql.Settings.instance().evaluationDate
        dc = ql.Actual365Fixed()

        maturityDate = today + ql.Period(10 * 365, ql.Days)

        v0 = 0.04
        kappa = 0.5
        theta = 0.04
        sigma = 1.0
        sig = 0.09
        rho = -0.9
        a = 0.08

        r = ql.YieldTermStructureHandle(ql.FlatForward(today, 0.05, dc))
        q = ql.YieldTermStructureHandle(ql.FlatForward(today, 0.03, dc))

        option = ql.VanillaOption(
            ql.PlainVanillaPayoff(ql.Option.Call, 100.0),
            ql.EuropeanExercise(maturityDate)
        )

        expected = 40.028973

        s0 = 100
        underlying = ql.makeQuoteHandle(s0)

        hull_white_model = ql.HullWhite(r, a, sig)
        heston_model = ql.HestonModel(
            ql.HestonProcess(r, q, underlying, v0, kappa, theta, sigma, rho)
        )

        option.setPricingEngine(
            ql.AnalyticHestonHullWhiteEngine(heston_model, hull_white_model)
        )
        self.assertAlmostEqual(expected, option.NPV(), 5)

        option.setPricingEngine(
            ql.AnalyticH1HWEngine(heston_model, hull_white_model, 0.0)
        )
        self.assertAlmostEqual(expected, option.NPV(), 5)

    def testCashDividendEuropeanEngine(self):
        """Testing cash dividend European engine"""
        
        today = today = ql.Settings.instance().evaluationDate
        dc = ql.Actual365Fixed()
    
        maturityDate = today + ql.Period(366, ql.Days)
    
        option = ql.VanillaOption(
            ql.PlainVanillaPayoff(ql.Option.Call, 100.0),
            ql.EuropeanExercise(maturityDate)
        )
    
        div_schedule = ql.DividendSchedule()
        div_schedule.append(ql.FixedDividend(5.0, today + ql.Period(92, ql.Days)))
    
        process = ql.BlackScholesMertonProcess(
            ql.makeQuoteHandle(100.0),
            ql.YieldTermStructureHandle(ql.FlatForward(today, 0.075, dc)),
            ql.YieldTermStructureHandle(ql.FlatForward(today, 0.05, dc)),
            ql.BlackVolTermStructureHandle(ql.BlackConstantVol(today, ql.TARGET(), 0.3,  dc))
        )
    
        option.setPricingEngine(
            ql.CashDividendEuropeanEngine(
                process, div_schedule, ql.CashDividendEuropeanEngine.Escrowed
            )
        )
        calculated = option.NPV()
    
        option.setPricingEngine(ql.AnalyticDividendEuropeanEngine(process, div_schedule))
        expected = option.NPV()
        self.assertAlmostEqual(expected, calculated, 8)
    
        option.setPricingEngine(
            ql.CashDividendEuropeanEngine(
                process, div_schedule, ql.CashDividendEuropeanEngine.Spot
            )
        )
        self.assertAlmostEqual(7.9193, option.NPV(), 3)


if __name__ == '__main__':
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
