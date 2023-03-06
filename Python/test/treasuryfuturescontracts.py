"""
 Copyright (C) 2023 Nijaz Kovacevic

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
import unittest

class TreasuryFuturesTest(unittest.TestCase):
    def setUp(self):
        self.calendar = ql.UnitedStates()
        self.businessConvention = ql.ModifiedFollowing
        self.settlementDays = 0
        self.daysCount = ql.ActualActual(ql.ActualActual.Bond)

        self.interestRate = 0.003
        self.calcDate = ql.Date(1,12,2023)
        self.yieldCurve = ql.FlatForward(self.calcDate, self.interestRate, self.daysCount, ql.Compounded, ql.Continuous)

        ql.Settings.instance().evaluationDate = self.calcDate
        self.optionMaturityDate = ql.Date(24,12,2025)
        self.strike = 100
        self.spot = 125 # spot price is the futures price
        self.volatility = 20/100.
        self.optionType = ql.Option.Call

        self.discount = self.yieldCurve.discount(self.optionMaturityDate)
        self.strikepayoff = ql.PlainVanillaPayoff(self.optionType, self.strike)
        self.T = self.yieldCurve.dayCounter().yearFraction(self.calcDate, self.optionMaturityDate)

        self.stddev = self.volatility * self.math.sqrt(self.T)

        self.black = ql.BlackCalculator(self.strikepayoff, self.spot, self.stddev, self.discount)
    
    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()
    
    def testBlackCalculator(self):
        """ Test the Black Calculator Engine """

        self.black = ql.BlackCalculator(self.strikepayoff, self.spot, 
                                    self.stddev, self.discount)
        
        self.assertAlmostEqual()
        self.strike = 100
        self.spot = 125
        self.volatility = 20/100
        self.optionType = ql.Option.Call

        self.discount = self.yieldCurve.discount(self.optionMaturityDate)
        self.strikepayoff = ql.PlainVanillaPayoff(self.optionType, self.strike)
        self.T = self.yieldCurve.dayCounter().yearFraction(self.calcDate, self.optionMaturityDate)

        self.stddev = self.volatility * self.math.sqrt(self.T)

        self.black = ql.BlackCalculator(self.strikepayoff, self.spot, self.stddev, self.discount)

if __name__ == '__main__':
    print('testing QuantLib ' + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TreasuryFuturesTest,'test'))
    unittest.TextTestRunner(verbosity=2).run(suite)