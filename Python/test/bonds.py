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
            couponCalendar=ql.UnitedStates(ql.UnitedStates.GovernmentBond),
            startDate=ql.Date(2, 1, 2010),
            maturityDate=self.maturity_date,
            tenor=ql.Period(3, ql.Months),
            accrualDayCounter=self.day_counter,
        )
        self.check_construction(bond)

class AmortizingFixedRateBondTest(unittest.TestCase):
    def test_interest_rates(self):
        # see AmortizingBondTest::testBrazilianAmortizingFixedRateBond
        # in the C++ test suite

        nominals = [
            1000       	, 983.33300000, 966.66648898, 950.00019204,
            933.33338867, 916.66685434, 900.00001759, 883.33291726,
            866.66619177, 849.99933423, 833.33254728, 816.66589633,
            799.99937871, 783.33299165, 766.66601558, 749.99946306,
            733.33297499, 716.66651646, 699.99971995, 683.33272661,
            666.66624140, 649.99958536, 633.33294599, 616.66615618,
            599.99951997, 583.33273330, 566.66633377, 549.99954356,
            533.33290739, 516.66625403, 499.99963400, 483.33314619,
            466.66636930, 449.99984658, 433.33320226, 416.66634063,
            399.99968700, 383.33290004, 366.66635221, 349.99953317,
            333.33290539, 316.66626012, 299.99948151, 283.33271031,
            266.66594695, 249.99932526, 233.33262024, 216.66590450,
            199.99931312, 183.33277035, 166.66617153, 149.99955437,
            133.33295388, 116.66633464,  99.99973207,  83.33307672,
             66.66646137,  49.99984602,  33.33324734,  16.66662367
        ]

        expected_amortizations = [
            16.66700000, 16.66651102, 16.66629694, 16.66680337,
            16.66653432, 16.66683675, 16.66710033, 16.66672548,
            16.66685753, 16.66678695, 16.66665095, 16.66651761,
            16.66638706, 16.66697606, 16.66655251, 16.66648807,
            16.66645852, 16.66679651, 16.66699333, 16.66648520,
            16.66665604, 16.66663937, 16.66678981, 16.66663620,
            16.66678667, 16.66639952, 16.66679021, 16.66663617,
            16.66665336, 16.66662002, 16.66648780, 16.66677688,
            16.66652271, 16.66664432, 16.66686163, 16.66665363,
            16.66678696, 16.66654783, 16.66681904, 16.66662777,
            16.66664527, 16.66677860, 16.66677119, 16.66676335,
            16.66662168, 16.66670502, 16.66671573, 16.66659137,
            16.66654276, 16.66659882, 16.66661715, 16.66660049,
            16.66661924, 16.66660257, 16.66665534, 16.66661534,
            16.66661534, 16.66659867, 16.66662367, 16.66662367
        ]

        expected_coupons = [
            5.97950399, 4.85474255, 5.27619136, 5.18522454,
            5.33753111, 5.24221882, 4.91231709, 4.59116258,
            4.73037674, 4.63940686, 4.54843737, 3.81920094,
            4.78359948, 3.86733691, 4.38439657, 4.09359456,
            4.00262671, 4.28531030, 3.82068947, 3.55165259,
            3.46502778, 3.71720657, 3.62189368, 2.88388676,
            3.58769952, 2.72800044, 3.38838360, 3.00196900,
            2.91100034, 3.08940793, 2.59877059, 2.63809514,
            2.42551945, 2.45615766, 2.59111761, 1.94857222,
            2.28751141, 1.79268582, 2.19248291, 1.81913832,
            1.90625855, 1.89350716, 1.48110584, 1.62031828,
            1.38600825, 1.23425366, 1.39521333, 1.06968563,
            1.03950542, 1.00065409, 0.90968563, 0.81871706,
            0.79726493, 0.63678002, 0.57187676, 0.49829046,
            0.32913418, 0.27290565, 0.19062560, 0.08662552
        ]

        settlementDays = 0
        issueDate = ql.Date(2, ql.March, 2020)
        maturityDate = ql.Date(2, ql.March, 2025)

        schedule = ql.Schedule(issueDate,
                               maturityDate,
                               ql.Period(ql.Monthly),
                               ql.Brazil(ql.Brazil.Settlement),
                               ql.Unadjusted,
                               ql.Unadjusted,
                               ql.DateGeneration.Backward,
                               False)

        coupons = ql.FixedRateLeg(
            schedule,
            nominals = nominals,
            couponRates = [0.0675],
            dayCount = ql.Business252(ql.Brazil()),
            compounding = ql.Compounded,
            compoundingFrequency = ql.Annual,
            paymentAdjustment = ql.Following,
        )

        bond = ql.Bond(
            settlementDays,
            schedule.calendar(),
            issueDate,
            coupons
        )

        cashflows = bond.cashflows()

        self.assertEqual(len(cashflows), 2 * len(nominals))

        for k in range(len(nominals)):
            self.assertEqual(round(expected_coupons[k], 5),  round(cashflows[2*k].amount(), 5))
            self.assertEqual(round(expected_amortizations[k], 5), round(cashflows[2*k+1].amount(), 5))


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(FixedRateBondTest, "test"))
    suite.addTest(unittest.makeSuite(FixedRateBondKwargsTest, "test"))
    suite.addTest(unittest.makeSuite(AmortizingFixedRateBondTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
