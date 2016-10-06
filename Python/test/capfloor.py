# coding=utf-8-unix
"""
 Copyright (C) 2016 Wojciech Ślusarski

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

class CapFloorTest(unittest.TestCase):
    def setUp(self):
        self.today_date = ql.Date(9, 9, 2016)
        self.settlement_days = 2
        self.notional = 1e4

        self.calendar = ql.TARGET()
        self.flat_forward_rate = 0.01
        self.rate_day_counter = ql.Actual360()
        self.flat_forward = ql.FlatForward(self.today_date,
                                           self.flat_forward_rate,
                                           self.rate_day_counter,
                                           ql.Continuous, ql.Annual)
        self.term_structure_handle = \
            ql.RelinkableYieldTermStructureHandle(self.flat_forward)

        self.interpolation = ql.Linear()

        self.start_date = ql.Date(13, 9, 2016)
        self.maturity_date = ql.Date(13, 9, 2017)
        self.period = ql.Period(6, ql.Months)
        self.buss_convention = ql.ModifiedFollowing
        self.date_gen_rule = ql.DateGeneration.Forward
        self.eom_rule = False

        self.schedule = ql.Schedule(self.start_date,
                                    self.maturity_date,
                                    self.period,
                                    self.calendar,
                                    self.buss_convention,
                                    self.buss_convention,
                                    self.date_gen_rule,
                                    self.eom_rule)

        ql.Settings.instance().evaluationDate = self.today_date

        self.ibor_index = ql.Euribor(self.period, self.term_structure_handle)

        self.ibor_index.addFixing(ql.Date(9, 9, 2016), 0.01)

        self.ibor_leg = ql.IborLeg([self.notional],
                                   self.schedule,
                                   self.ibor_index)
        self.strike = 0.01
        self.cap = ql.Cap(self.ibor_leg, [self.strike])
        self.cap_npv = 8.54

        self.black_vol = ql.QuoteHandle(ql.SimpleQuote(0.6))

    def testBlackCapFloorEngine(self):
        """ Testing BlackCapFloorEngine """
        black_engine = ql.BlackCapFloorEngine(self.term_structure_handle,
                                        self.black_vol)

        self.cap.setPricingEngine(black_engine)
        npv = self.cap.NPV()
        self.assertAlmostEqual(npv, self.cap_npv,
                               places=1, msg="NPV method is broken")
        vol_guess = 0.5
        imp_vol = self.cap.impliedVolatility(npv,
                                             self.term_structure_handle,
                                             vol_guess)
        self.assertAlmostEqual(self.black_vol.value(),
                               imp_vol, places=4,
                               msg="Implied volatility method is broken")


    def testBachelierCapFloorEngine(self):
        """ Testing BachelierCapFloorEngine """

        bpvol = self.black_vol.value() * self.flat_forward_rate
        bachelier_engine = ql.BachelierCapFloorEngine(self.term_structure_handle,
                                                      ql.QuoteHandle(
                                                         ql.SimpleQuote(bpvol)))
                                                                    
        self.cap.setPricingEngine(bachelier_engine)

        # 50 bps
        vol_guess = 50 / 1e4

        imp_vol = self.cap.impliedVolatility(self.cap_npv, 
                                             self.term_structure_handle,
                                             vol_guess,
                                             type=ql.Normal)

        self.assertAlmostEqual(bpvol, imp_vol, places=4,
                               msg="Normal Implied volatility method is broken")



if __name__ == '__main__':
    print('testing QuantLib ' + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(CapFloorTest,'test'))
    unittest.TextTestRunner(verbosity=2).run(suite)
