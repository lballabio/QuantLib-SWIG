"""
 Copyright (C) 2009 Joseph Malicki
 Copyright (C) 2019 Prasad Somwanshi

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


class FixedRateBondTest(unittest.TestCase):
    def setUp(self):
        ql.Settings.instance().evaluationDate = ql.Date(2, 1, 2010)
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

        self.bond = ql.FixedRateBond(
            self.settlement_days,
            self.face_amount,
            self.sched,
            self.coupons,
            self.day_counter,
            ql.Following,
            self.redemption,
            self.issue_date,
        )

        self.flat_forward = ql.FlatForward(
            self.issue_date, self.coupons[0], self.day_counter, ql.Compounded, ql.Semiannual
        )
        self.term_structure_handle = ql.RelinkableYieldTermStructureHandle(self.flat_forward)
        bondEngine = ql.DiscountingBondEngine(self.term_structure_handle)
        self.bond.setPricingEngine(bondEngine)

    def testFrequency(self):
        """ Testing FixedRateBond frequency() method. """
        self.assertEqual(self.bond.frequency(), ql.Semiannual)

    def testDayCounter(self):
        """ Testing FixedRateBond dayCounter() method. """
        self.assertEqual(self.bond.dayCounter(), self.day_counter)

    def testSimpleInspectors(self):
        """ Testing FixedRateBond simple inspectors. """
        self.assertEqual(self.bond.settlementDays(), self.settlement_days)
        self.assertEqual(self.bond.notional(), self.face_amount)
        self.assertEqual(self.bond.issueDate(), self.issue_date)
        self.assertEqual(self.bond.maturityDate(), self.maturity_date)

    # def testSettlementValue(self):
    #    """ Testing FixedRateBond settlement value. """
    #    orig_date = ql.Settings.evaluationDate
    #    ql.Settings.evaluationDate = self.issue_date + 1*ql.Months
    #    self.assertEqual(round(self.bond.settlementValue(100.0), 4), 102.3098)
    #    ql.Settings.evaluationDate = orig_date

    def testCashFlows(self):
        """ Testing that the FixedRateBond gives the expected cash flows. """
        self.assertEqual(
            [round(cf.amount(), 4) for cf in self.bond.cashflows()],
            20 * [round(self.face_amount * self.coupons[0] / 2, 4)] + [round(self.redemption, 4)],
        )

    def testRedemption(self):
        """ Testing FixedRateBond redemption value and date. """
        self.assertEqual(self.bond.redemption().date(), self.maturity_date)
        self.assertEqual(self.bond.redemption().amount(), self.redemption)

    def testRedemptions(self):
        """ Testing FixedRateBond redemptions. """
        redemptions = self.bond.redemptions()
        self.assertEqual(len(redemptions), 1)
        self.assertEqual(redemptions[0].date(), self.maturity_date)
        self.assertEqual(redemptions[0].amount(), self.redemption)

    def testNotional(self):
        """ Testing FixedRateBond notional values. """
        self.assertEqual(self.bond.notional(), 100.0)
        self.assertEqual(self.bond.notionals(), (100.0, 0))

    def testNextCoupon(self):
        """ Testing FixedRateBond correct next coupon amount. """
        self.assertEqual(self.bond.nextCouponRate(self.issue_date), 0.05)

    def testPrevCoupon(self):
        """ Testing FixedRateBond correct previous coupon amount. """
        self.assertEqual(self.bond.previousCouponRate(), 0.05)

    def testCleanPrice(self):
        """ Testing FixedRateBond clean price. """
        self.assertEqual(
            round(self.bond.cleanPrice(0.05, self.day_counter, ql.Compounded, ql.Semiannual, self.issue_date), 4),
            99.9964,
        )
        self.assertEqual(
            round(
                self.bond.cleanPrice(
                    0.05, self.day_counter, ql.Compounded, ql.Semiannual, self.issue_date + ql.Period(1, ql.Months)
                ),
                4,
            ),
            99.9921,
        )

        self.assertEqual(
            round(
                self.bond.cleanPrice(
                    0.06, self.day_counter, ql.Compounded, ql.Semiannual, self.issue_date + ql.Period(1, ql.Months)
                ),
                4,
            ),
            92.5985,
        )

    def testDirtyPrice(self):
        """ Testing FixedRateBond dirty price. """
        self.assertEqual(
            round(self.bond.dirtyPrice(0.05, self.day_counter, ql.Compounded, ql.Semiannual, self.issue_date), 4),
            99.9964,
        )
        self.assertEqual(
            round(
                self.bond.dirtyPrice(
                    0.05, self.day_counter, ql.Compounded, ql.Semiannual, self.issue_date + ql.Period(1, ql.Months)
                ),
                4,
            ),
            100.4179,
        )
        self.assertEqual(
            round(
                self.bond.dirtyPrice(
                    0.06, self.day_counter, ql.Compounded, ql.Semiannual, self.issue_date + ql.Period(1, ql.Months)
                ),
                4,
            ),
            93.0244,
        )

    def testCleanPriceFromZSpread(self):
        """ Testing FixedRateBond clean price derived from Z-spread. """
        self.assertEqual(
            round(
                ql.cleanPriceFromZSpread(
                    self.bond,
                    self.flat_forward,
                    0.01,
                    self.day_counter,
                    ql.Compounded,
                    ql.Semiannual,
                    self.issue_date + ql.Period(1, ql.Months),
                ),
                4,
            ),
            92.5926,
        )

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()


class FixedRateBondKwargsTest(unittest.TestCase):
    def setUp(self):
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

    def check_construction(self, bond):
        self.assertTrue(type(bond) is ql.FixedRateBond)
        self.assertEqual(bond.dayCounter(), self.day_counter)
        self.assertEqual(bond.settlementDays(), self.settlement_days)
        self.assertEqual(bond.issueDate(), self.issue_date)
        self.assertEqual(bond.maturityDate(), self.maturity_date)
        self.assertEqual(bond.redemption().date(), self.maturity_date)
        self.assertEqual(bond.redemption().amount(), self.redemption)
        self.assertEqual(bond.notional(self.issue_date), 100.0)
        self.assertEqual(bond.notionals(), (100.0, 0))

    def testFromRates(self):
        """ Testing FixedRateBond from_rates method. """
        bond = ql.FixedRateBond.from_rates(
            settlementDays=self.settlement_days,
            schedule=self.sched,
            paymentDayCounter=self.day_counter,
            issueDate=self.issue_date,
            coupons=self.coupons,
            faceAmount=self.face_amount,
        )
        self.check_construction(bond)

    def testFromInterestRates(self):
        """ Testing FixedRateBond from_interest_rates method. """
        bond = ql.FixedRateBond.from_interest_rates(
            settlementDays=self.settlement_days,
            faceAmount=self.face_amount,
            schedule=self.sched,
            coupons=[ql.InterestRate(0.05, self.day_counter, ql.Continuous, ql.Annual)],
            issueDate=self.issue_date,
        )
        self.check_construction(bond)

    def testFromDateInfo(self):
        """ Testing FixedRateBond from_interest_rates method. """
        bond = ql.FixedRateBond.from_date_info(
            settlementDays=self.settlement_days,
            faceAmount=self.face_amount,
            coupons=self.coupons,
            issueDate=self.issue_date,
            couponCalendar=ql.UnitedStates(),
            startDate=ql.Date(2, 1, 2010),
            maturityDate=self.maturity_date,
            tenor=ql.Period(3, ql.Months),
            accrualDayCounter=self.day_counter,
        )
        self.check_construction(bond)


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(FixedRateBondTest, "test"))
    suite.addTest(unittest.makeSuite(FixedRateBondKwargsTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
