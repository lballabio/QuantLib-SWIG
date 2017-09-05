"""
 Copyright (C) 2017 Wojciech Ślusarski

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
import math

class BlackDeltaCalculatorTest(unittest.TestCase):

    def setUp(self):
        self.todaysDate = ql.Date(5, ql.September, 2017)
        ql.Settings.instance().evaluationDate = self.todaysDate
        self.spotDate = ql.Date(7, ql.September, 2017)
        self.domestic_rate = ql.FlatForward(self.spotDate, 0.017,
                                            ql.Actual365Fixed())
        self.foreign_rate = ql.FlatForward(self.spotDate, 0.013,
                                           ql.Actual365Fixed())

    def test_single_spot_delta(self):
        """Test for a single strike for call spot delta 75"""
        volatility = 0.2
        expiry = 2
        spot_price = 3.6
        domDf = self.domestic_rate.discount(expiry)
        forDf = self.foreign_rate.discount(expiry)
        forward = spot_price * forDf / domDf

        spot_delta_level = 0.75
        stDev = volatility * expiry ** 0.5

        inv_norm_dist = ql.InverseCumulativeNormal()
        expected_strike = inv_norm_dist(spot_delta_level / forDf)
        expected_strike *= stDev
        expected_strike -= 0.5 * stDev ** 2
        expected_strike = math.exp(expected_strike) / forward
        expected_strike = 1 / expected_strike

        option_type = ql.Option.Call
        delta_type = ql.DeltaVolQuote.Spot

        black_calculator = ql.BlackDeltaCalculator(option_type,
                                                   delta_type,
                                                   spot_price,
                                                   domDf,
                                                   forDf,
                                                   stDev)



        strike = black_calculator.strikeFromDelta(spot_delta_level)

        self.assertAlmostEquals(expected_strike, strike, delta=1e-4)

    def test_spot_atm_delta_calculator(self):
        """Test for 0-delta straddle strike"""
        volatility = 0.2
        expiry = 2
        spot_price = 3.6
        domDf = self.domestic_rate.discount(expiry)
        forDf = self.foreign_rate.discount(expiry)
        forward = spot_price * forDf / domDf
        expected_strike = forward * math.exp(-0.5 * volatility ** 2 * expiry)

        option_type = ql.Option.Call
        delta_type = ql.DeltaVolQuote.AtmDeltaNeutral
        stDev = volatility * expiry ** 0.5

        black_calculator = ql.BlackDeltaCalculator(option_type,
                                                   delta_type,
                                                   spot_price,
                                                   domDf,
                                                   forDf,
                                                   stDev)

        strike = black_calculator.atmStrike(ql.DeltaVolQuote.AtmDeltaNeutral)

        self.assertAlmostEquals(expected_strike, strike, delta=1e-4)


if __name__ == '__main__':
    unittest.main()
