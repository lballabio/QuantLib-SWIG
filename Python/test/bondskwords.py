"""
 Copyright (C) 2009 Joseph Malicki

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


class FixedRateBondKwordsTest(unittest.TestCase):
    def setUp(self):
        ql.Settings.instance().setEvaluationDate(ql.Date(2, 1, 2010))
        self.settlement_days = 3
        self.face_amount = 100.0
        self.redemption = 100.0
        self.issue_date = ql.Date(2, 1, 2008)
        self.maturity_date = ql.Date(2, 1, 2018)
        self.calendar = ql.UnitedStates(ql.UnitedStates.GovernmentBond)
        self.day_counter = ql.ActualActual(ql.ActualActual.Bond)
        self.sched = ql.Schedule(
            self.issue_date,
            self.maturity_date,
            ql.Period(ql.Semiannual),
            self.calendar,
            ql.Unadjusted,
            ql.Unadjusted,
            ql.DateGeneration.Backward,
            False,
        )
        self.coupons = [0.05]

    def simple_inspectors(self,bond):
        """ Testing FixedRateBond simple inspectors. """
        self.assertTrue(type(bond) is ql.FixedRateBond)
        self.assertEqual(bond.dayCounter(), self.day_counter)
        self.assertEqual(bond.settlementDays(), self.settlement_days)
        self.assertEqual(bond.notional(), self.face_amount)
        self.assertEqual(bond.issueDate(), self.issue_date)
        self.assertEqual(bond.maturityDate(), self.maturity_date)

    def testFromRates(self):
        """ Testing FixedRateBond from_rates method. """
        bond = ql.FixedRateBond.from_rates(
                                            settlementDays = self.settlement_days,
                                            faceAmount = self.face_amount,
                                            schedule = self.sched,
                                            paymentDayCounter = self.day_counter,
                                            issueDate = self.issue_date,
                                            coupons = self.coupons)
        self.simple_inspectors(bond)
        self.check_notional(bond)
        self.check_redemption(bond)

    def testFromInterestRates(self):
        """ Testing FixedRateBond from_interest_rates method. """
        bond = ql.FixedRateBond.from_interest_rates(
                                            settlementDays = self.settlement_days,
                                            faceAmount = self.face_amount,
                                            schedule = self.sched,
                                            coupons = [ql.InterestRate(0.05,self.day_counter,ql.Continuous,ql.Annual)],
                                            issueDate = self.issue_date)
        self.simple_inspectors(bond)
        self.check_notional(bond)
        self.check_redemption(bond)

    def testFromRatesCalc(self):
        """ Testing FixedRateBond from_interest_rates method. """
        coupons = [ql.InterestRate(0.05,self.day_counter,ql.Continuous,ql.Annual)]
        tenor = ql.Period(3, ql.Months)
        bond = ql.FixedRateBond.from_rates_schedule_calc(
                                            settlementDays = self.settlement_days,
                                            faceAmount = self.face_amount,
                                            coupons = self.coupons,
                                            issueDate = self.issue_date,
                                            couponCalendar = ql.UnitedStates(),
                                            startDate = ql.Date(2, 1, 2010),
                                            maturityDate = self.maturity_date,
                                            tenor = tenor,
                                            accrualDayCounter = self.day_counter
                                            )
        self.simple_inspectors(bond)
        self.check_notional(bond)
        self.check_redemption(bond)
        self.assertEqual(bond.dayCounter(),self.day_counter)

    def check_redemption(self,bond):
        """ Testing FixedRateBond redemption value and date. """
        self.assertEqual(bond.redemption().date(), self.maturity_date)
        self.assertEqual(bond.redemption().amount(), self.redemption)

    def check_notional(self,bond):
        """ Testing FixedRateBond notional values. """
        self.assertEqual(bond.notional(), 100.0)
        self.assertEqual(bond.notionals(), (100.0, 0))

    def tearDown(self):
        ql.Settings.instance().setEvaluationDate(ql.Date())


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(FixedRateBondKwordsTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
