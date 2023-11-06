"""
 Copyright (C) 2009 Joseph Malicki
 Copyright (C) 2019 Prasad Somwanshi
 Copyright (C) 2023 Francois Botha

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


class BondFunctionsTest(unittest.TestCase):
    def setUp(self):
        ql.Settings.instance().evaluationDate = ql.Date(2, 1, 2010)
        self.settlement_days = 3
        self.face_amount = 100.0
        self.redemption = 100.0
        self.issue_date = ql.Date(2, 1, 2008)
        self.maturity_date = ql.Date(2, 1, 2018)
        self.calendar = ql.UnitedStates(ql.UnitedStates.GovernmentBond)
        self.settlement_date = self.calendar.advance(
            ql.Settings.instance().evaluationDate, self.settlement_days, ql.Days)
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
        self.term_structure_handle = ql.RelinkableYieldTermStructureHandle(
            self.flat_forward)
        bondEngine = ql.DiscountingBondEngine(self.term_structure_handle)
        self.bond.setPricingEngine(bondEngine)

    def testStartDate(self):
        """ Testing BondFunctions startDate. """
        self.assertEqual(ql.BondFunctions.startDate(
            self.bond), self.issue_date)

    def testMaturityDate(self):
        """ Testing BondFunctions maturityDate. """
        self.assertEqual(ql.BondFunctions.maturityDate(
            self.bond), self.maturity_date)

    def testIsTradable(self):
        """ Testing BondFunctions isTradable. """
        self.assertTrue(ql.BondFunctions.isTradable(
            self.bond, ql.Date(1, 6, 2010)))

        self.assertFalse(ql.BondFunctions.isTradable(
            self.bond, ql.Date(1, 1, 2028)))

    def testPreviousCashFlowDate(self):
        """ Testing BondFunctions previousCashFlowDate. """
        self.assertEqual(ql.BondFunctions.previousCashFlowDate(self.bond, ql.Date(1, 6, 2010)),
                         ql.Date(4, 1, 2010))

    def testNextCashFlowDate(self):
        """ Testing BondFunctions nextCashFlowDate. """
        self.assertEqual(ql.BondFunctions.nextCashFlowDate(self.bond, ql.Date(1, 6, 2010)),
                         ql.Date(2, 7, 2010))

    def testPreviousCashFlowAmount(self):
        """ Testing BondFunctions previousCashFlowAmount. """
        self.assertEqual(round(ql.BondFunctions.previousCashFlowAmount(
            self.bond, ql.Date(1, 6, 2010)), 4), 2.5)

    def testNextCashFlowAmount(self):
        """ Testing BondFunctions nextCashFlowAmount. """
        self.assertEqual(round(ql.BondFunctions.nextCashFlowAmount(
            self.bond, ql.Date(1, 6, 2010)), 4), 2.5)

    def testPreviousCouponRate(self):
        """ Testing BondFunctions previousCouponRate. """
        self.assertEqual(ql.BondFunctions.previousCouponRate(self.bond), 0.05)

    def testNextCouponRate(self):
        """ Testing BondFunctions nextCouponRate. """
        self.assertEqual(ql.BondFunctions.nextCouponRate(self.bond), 0.05)

    def testAccrualStartDate(self):
        """ Testing BondFunctions accrualStartDate. """
        self.assertEqual(ql.BondFunctions.accrualStartDate(self.bond, ql.Date(1, 6, 2010)),
                         ql.Date(2, 1, 2010))

    def testAccrualEndDate(self):
        """ Testing BondFunctions accrualEndDate. """
        self.assertEqual(ql.BondFunctions.accrualEndDate(self.bond, ql.Date(1, 6, 2010)),
                         ql.Date(2, 7, 2010))

    def testAccrualPeriod(self):
        """ Testing BondFunctions accrualPeriod. """
        self.assertEqual(ql.BondFunctions.accrualPeriod(self.bond, ql.Date(1, 6, 2010)),
                         0.5)

    def testAccrualDays(self):
        """ Testing BondFunctions accrualDays. """
        self.assertEqual(ql.BondFunctions.accrualDays(self.bond, ql.Date(1, 10, 2010)),
                         184)

    def testAccruedPeriod(self):
        """ Testing BondFunctions accruedPeriod. """
        self.assertEqual(round(ql.BondFunctions.accruedPeriod(self.bond, ql.Date(1, 6, 2010)), 8),
                         0.41436464)

    def testAccruedDays(self):
        """ Testing BondFunctions accruedDays. """
        self.assertEqual(ql.BondFunctions.accruedDays(self.bond, ql.Date(1, 6, 2010)),
                         150)

    def testAccruedAmount(self):
        """ Testing BondFunctions accruedAmount. """
        self.assertEqual(round(ql.BondFunctions.accruedAmount(self.bond, ql.Date(1, 6, 2010)), 8),
                         2.0718232)

    def testBps(self):
        """ Testing BondFunctions bps. """
        self.assertEqual(round(ql.BondFunctions.bps(self.bond, self.flat_forward), 8),
                         0.06527501)
        self.assertEqual(round(ql.BondFunctions.bps(self.bond,
                                                    ql.InterestRate(0.03, self.day_counter, ql.Compounded, ql.Annual)), 8),
                         0.07071951)
        self.assertEqual(round(ql.BondFunctions.bps(self.bond,
                                                    0.03, self.day_counter, ql.Compounded, ql.Annual), 8),
                         0.07071951)

    def testCleanPrice(self):
        """ Testing BondFunctions cleanPrice. """
        self.assertEqual(round(ql.BondFunctions.cleanPrice(self.bond, self.flat_forward), 4),
                         99.9448)
        self.assertEqual(round(ql.BondFunctions.cleanPrice(self.bond,
                                                           ql.InterestRate(0.03, self.day_counter, ql.Compounded, ql.Annual)), 4),
                         114.2806)
        self.assertEqual(round(ql.BondFunctions.cleanPrice(self.bond,
                                                           0.03, self.day_counter, ql.Compounded, ql.Annual), 4),
                         114.2806)

    def testAtmRate(self):
        """ Testing BondFunctions atmRate. """
        self.assertEqual(round(ql.BondFunctions.atmRate(self.bond, self.flat_forward,
                                                        self.settlement_date, 99.94475138121548), 4),
                         0.05)

    def testBondYield(self):
        """ Testing BondFunctions bondYield. """
        self.assertEqual(round(ql.BondFunctions.bondYield(self.bond, 110, self.day_counter, ql.Compounded, ql.Annual), 8),
                         0.03582431)

    def testDuration(self):
        """ Testing BondFunctions duration. """
        self.assertEqual(round(ql.BondFunctions.duration(self.bond,
                                                         ql.InterestRate(0.03, self.day_counter, ql.Compounded, ql.Annual)), 4),
                         6.5835)
        self.assertEqual(round(ql.BondFunctions.duration(self.bond,
                                                         0.03, self.day_counter, ql.Compounded, ql.Annual), 4),
                         6.5835)

    def testConvexity(self):
        """ Testing BondFunctions convexity. """
        self.assertEqual(round(ql.BondFunctions.convexity(self.bond,
                                                          ql.InterestRate(0.03, self.day_counter, ql.Compounded, ql.Annual)), 4),
                         54.3498)
        self.assertEqual(round(ql.BondFunctions.convexity(self.bond,
                                                          0.03, self.day_counter, ql.Compounded, ql.Annual), 4),
                         54.3498)

    def testBasisPointValue(self):
        """ Testing BondFunctions basisPointValue. """
        self.assertEqual(round(ql.BondFunctions.basisPointValue(self.bond,
                                                                0.03, self.day_counter, ql.Compounded, ql.Annual,
                                                                self.settlement_date), 8),
                         -0.07527271)
        self.assertEqual(round(ql.BondFunctions.basisPointValue(self.bond,
                                                                ql.InterestRate(
                                                                    0.03, self.day_counter, ql.Compounded, ql.Annual),
                                                                self.settlement_date), 8),
                         -0.07527271)

    def testYieldValueBasisPoint(self):
        """ Testing BondFunctions yieldValueBasisPoint. """
        self.assertEqual(round(ql.BondFunctions.yieldValueBasisPoint(self.bond,
                                                                     ql.InterestRate(
                                                                         0.03, self.day_counter, ql.Compounded, ql.Annual),
                                                                     ql.Date(1, 9, 2010)), 10),
                         -1.44145e-05)
        self.assertEqual(round(ql.BondFunctions.yieldValueBasisPoint(self.bond,
                                                                     0.03, self.day_counter, ql.Compounded, ql.Annual,
                                                                     ql.Date(1, 9, 2010)), 10),
                         -1.44145e-05)

    def testZSpread(self):
        """ Testing BondFunctions zSpread. """
        self.assertEqual(round(ql.BondFunctions.zSpread(self.bond, 87.5, self.flat_forward,
                                                        self.day_counter, ql.Compounded, ql.Annual), 8),
                         0.02125053)


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
