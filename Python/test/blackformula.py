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
import math
import QuantLib as ql


class BlackFormulaTest(unittest.TestCase):

    def setUp(self):
        # define the market and option parameters
        self.option_type = ql.Option.Call
        self.spot = 100.0
        self.strike = 100.0
        self.risk_free_rate = 0.05
        self.expiry = 1.0
        self.forward = self.spot * math.exp(self.risk_free_rate * self.expiry)
        self.df = math.exp(-self.risk_free_rate * self.expiry)
        self.vol = 0.2 * math.sqrt(self.expiry)
        self.displacement = 0.0

    def test_blackFormula(self):
        """Testing blackFormula in a simple Black-Scholes World..."""
        #Anyone interested, feel free to provide more accurate number
        expected = 10.4506
        res = ql.blackFormula(self.option_type,
                                 self.strike,
                                 self.forward,
                                 self.vol,
                                 self.df,
                                 self.displacement)
        self.assertAlmostEquals(expected, res, delta=1e-4,
                                msg="Failed to calculate simple  "
                                    "Black-Scholes-Merton price rounded to "
                                    "four decimal places.")

    def test_black_formula_implied_stdev(self):
        """Testing implied volatility calculator"""
        expected = 0.2 * math.sqrt(self.expiry)
        black_price = 10.4506
        res = ql.blackFormulaImpliedStdDev(self.option_type,
                                              self.strike,
                                              self.forward,
                                              black_price,
                                              self.df)
        self.assertAlmostEquals(expected, res, delta=1e-4,
                                msg="Failed to determine Implied Vol rounded " \
                                    "to a single vol bps.")


if __name__ == '__main__':
    unittest.main()
