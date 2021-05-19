"""
 Copyright (C) 2021 Klaus Spanderen

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
 
import unittest

import QuantLib as ql

class OptionsTest(unittest.TestCase):

    def testFdHestonHullWhite(self):
        """ Testing FDM Heston Hull-White pricing """
        
        dc = ql.Actual365Fixed()
        todays_date = ql.Date(19, ql.May, 2021)        
        
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
        
        underlying = ql.QuoteHandle(ql.SimpleQuote(s0))
        
        option = ql.VanillaOption(
            ql.PlainVanillaPayoff(ql.Option.Call, s0),
            ql.EuropeanExercise(todays_date + ql.Period(1, ql.Years))
        )

        hull_white_process = ql.HullWhiteProcess(r, a, sig)
        heston_process = ql.HestonProcess(r, d, underlying, v0, kappa, theta, sigma, rho)        

        option.setPricingEngine(
            ql.FdHestonHullWhiteVanillaEngine(
                ql.HestonModel(heston_process), hull_white_process, -0.5,
                10, 200, 25, 10, 
                controlVariate=True
            )
        )            
        
        self.assertAlmostEqual(0.87628, option.NPV(), 4)


if __name__ == '__main__':
    print('testing QuantLib ' + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(OptionsTest,'test'))
    unittest.TextTestRunner(verbosity=2).run(suite)
 