"""
 Copyright (C) 2021 Marcin Rybacki
 Copyright (C) 2023 Marcin Rybacki

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


EPSILON = 1.e-8

CAL = ql.TARGET()

DCT = ql.Actual365Fixed()

IR_FIXINGS = [(ql.Date(3, ql.January, 2023), 0.033),
              (ql.Date(4, ql.January, 2023), 0.033),
              (ql.Date(5, ql.January, 2023), 0.033),
              (ql.Date(6, ql.January, 2023), 0.033),
              (ql.Date(9, ql.January, 2023), 0.03),
              (ql.Date(10, ql.January, 2023), 0.03),
              (ql.Date(11, ql.January, 2023), 0.03),
              (ql.Date(12, ql.January, 2023), 0.03),
              (ql.Date(13, ql.January, 2023), 0.03),
              (ql.Date(17, ql.January, 2023), 0.03),
              (ql.Date(20, ql.January, 2023), 0.03),
              (ql.Date(23, ql.January, 2023), 0.03),
              (ql.Date(24, ql.January, 2023), 0.03),
              (ql.Date(25, ql.January, 2023), 0.03),
              (ql.Date(26, ql.January, 2023), 0.03)]


def flat_rate(rate):
    return ql.FlatForward(
        2, CAL, ql.QuoteHandle(ql.SimpleQuote(rate)), ql.Actual365Fixed())


class ZeroCouponSwapTest(unittest.TestCase):
    def setUp(self):
        valuation_date = CAL.adjust(ql.Date(1, ql.June, 2021))
        ql.Settings.instance().evaluationDate = valuation_date
        self.nominal_ts_handle = ql.YieldTermStructureHandle(flat_rate(0.007))
        self.ibor_idx = ql.Euribor6M(self.nominal_ts_handle)
        self.engine = ql.DiscountingSwapEngine(self.nominal_ts_handle)

    def build_zcs_from_fixed_payment(self, amount):
        return ql.ZeroCouponSwap(ql.Swap.Receiver,
                                 1.0e6,
                                 ql.Date(3, ql.June, 2021),
                                 ql.Date(3, ql.June, 2051),
                                 amount,
                                 self.ibor_idx,
                                 CAL)

    def build_zcs_from_rate(self, rate):
        return ql.ZeroCouponSwap(ql.Swap.Receiver,
                                 1.0e6,
                                 ql.Date(3, ql.June, 2021),
                                 ql.Date(3, ql.June, 2051),
                                 rate,
                                 DCT,
                                 self.ibor_idx,
                                 CAL)

    def test_zero_coupon_swap_inspectors(self):
        """Testing zero coupon swap inspectors"""
        swap = self.build_zcs_from_fixed_payment(1.5e6)
        fail_msg = "Unable to replicate the properties of a ZC swap."

        self.assertEqual(swap.type(), ql.Swap.Receiver,
                         msg=fail_msg)
        self.assertEqual(swap.startDate(), ql.Date(3, ql.June, 2021),
                         msg=fail_msg)
        self.assertEqual(swap.maturityDate(), ql.Date(3, ql.June, 2051),
                         msg=fail_msg)
        self.assertAlmostEqual(swap.baseNominal(), 1.0e6,
                               delta=EPSILON, msg=fail_msg)
        self.assertAlmostEqual(swap.fixedPayment(), 1.5e6,
                               delta=EPSILON, msg=fail_msg)

    def test_npvs_of_par_zero_coupon_swap_with_fixed_payment(self):
        """Testing NPVs of a zero coupon swap with fixed payment"""
        swap = self.build_zcs_from_fixed_payment(1.5e6)
        swap.setPricingEngine(self.engine)
        fair_payment = swap.fairFixedPayment()
        par_swap = self.build_zcs_from_fixed_payment(fair_payment)
        par_swap.setPricingEngine(self.engine)
        npv = par_swap.NPV()
        fail_npv_msg = """ Unable to replicate par zero coupon swap NPV:
                            calculated: {actual}
                            expected: {expected}
                       """.format(actual=npv,
                                  expected=0.0)
        self.assertAlmostEqual(npv, 0.0, delta=EPSILON, msg=fail_npv_msg)

        fxd_leg_npv = par_swap.fixedLegNPV()
        flt_leg_npv = par_swap.floatingLegNPV()
        fail_legs_npv_msg = """ Unable to replicate the NPVs of a par zero coupon swap legs:
                                 fixed leg NPV: {fxd_leg}
                                 floating leg NPV: {flt_leg}
                            """.format(fxd_leg=fxd_leg_npv,
                                       flt_leg=flt_leg_npv)
        self.assertAlmostEqual(abs(fxd_leg_npv), abs(flt_leg_npv),
                               delta=EPSILON,
                               msg=fail_legs_npv_msg)

    def test_npvs_of_par_zero_coupon_swap_with_fixed_rate(self):
        """Testing NPVs of a zero coupon swap with fixed rate"""
        swap = self.build_zcs_from_fixed_payment(1.5e6)
        swap.setPricingEngine(self.engine)
        fair_rate = swap.fairFixedRate(DCT)
        par_swap = self.build_zcs_from_rate(fair_rate)
        par_swap.setPricingEngine(self.engine)
        npv = par_swap.NPV()
        fail_msg = """ Unable to replicate par zero coupon swap NPV:
                        calculated: {actual}
                        expected: {expected}
                   """.format(actual=npv,
                              expected=0.0)
        self.assertAlmostEqual(npv, 0.0, delta=EPSILON, msg=fail_msg)

    def test_zero_coupon_swap_legs(self):
        """Testing zero coupon swap legs"""
        swap = self.build_zcs_from_rate(0.01)
        fxd_leg = swap.fixedLeg()
        fxd_cf = ql.as_fixed_rate_coupon(fxd_leg[0])
        fail_msg_fxd = """Fixed leg cash flow type should be FixedRateCoupon
                          but was {actual}.
                       """.format(actual=type(fxd_cf))
        self.assertTrue(isinstance(fxd_cf, ql.FixedRateCoupon),
                        msg=fail_msg_fxd)

        flt_leg = swap.floatingLeg()
        flt_cf = ql.as_sub_periods_coupon(flt_leg[0])
        fail_msg_flt = """Floating leg cash flow type should be SubPeriodsCoupon
                          but was {actual}.
                       """.format(actual=type(flt_cf))
        self.assertTrue(isinstance(
            flt_cf, ql.SubPeriodsCoupon), msg=fail_msg_flt)


class EquityTotalReturnSwapTest(unittest.TestCase):
    def setUp(self):
        valuation_date = ql.Date(27, ql.January, 2023)
        ql.Settings.instance().evaluationDate = valuation_date

        self.interest_handle = ql.YieldTermStructureHandle(flat_rate(0.03))
        self.dividend_handle = ql.YieldTermStructureHandle(flat_rate(0.0))
        equity_spot = ql.QuoteHandle(ql.SimpleQuote(8690.0))

        self.equity_idx = ql.EquityIndex(
            "eq_idx",
            CAL,
            self.interest_handle,
            self.dividend_handle,
            equity_spot)
        ql.IndexManager.instance().clearHistory(self.equity_idx.name())
        self.equity_idx.addFixing(ql.Date(5, ql.January, 2023), 9010.0)

        self.ibor_idx = ql.USDLibor(
            ql.Period(3, ql.Months), self.interest_handle)
        ql.IndexManager.instance().clearHistory(self.ibor_idx.name())
        self.sofr_idx = ql.Sofr(self.interest_handle)
        ql.IndexManager.instance().clearHistory(self.sofr_idx.name())

        for f_dt, f_val in IR_FIXINGS:
            self.ibor_idx.addFixing(f_dt, f_val)
            self.sofr_idx.addFixing(f_dt, f_val)

    def build_trs(self, interest_idx, start, end, margin=0.025):
        schedule = ql.Schedule(
            start,
            end,
            interest_idx.tenor(),
            interest_idx.fixingCalendar(),
            interest_idx.businessDayConvention(),
            interest_idx.businessDayConvention(),
            ql.DateGeneration.Backward,
            False)
        return ql.EquityTotalReturnSwap(ql.Swap.Receiver,
                                        1.0e6,
                                        schedule,
                                        self.equity_idx,
                                        interest_idx,
                                        DCT,
                                        margin)

    def test_trs_interest_rate_index(self):
        """Testing equity total return swap interest rate index"""
        start = ql.Date(5, ql.January, 2023)
        end = ql.Date(5, ql.April, 2023)

        trs_vs_ibor = self.build_trs(self.ibor_idx, start, end)
        trs_vs_sofr = self.build_trs(self.sofr_idx, start, end)

        fail_msg = "Incorrect interest rate index set to TRS."

        self.assertEqual(trs_vs_ibor.interestRateIndex().name(),
                         "USDLibor3M Actual/360",
                         msg=fail_msg)
        self.assertEqual(trs_vs_sofr.interestRateIndex().name(),
                         "SOFRON Actual/360",
                         msg=fail_msg)

    def test_trs_npv(self):
        """Testing equity total return swap NPV"""
        start = ql.Date(5, ql.January, 2023)
        end = ql.Date(5, ql.April, 2023)

        pricer = ql.DiscountingSwapEngine(self.interest_handle)

        trs_vs_ibor = self.build_trs(self.ibor_idx, start, end)
        trs_vs_ibor.setPricingEngine(pricer)

        trs_vs_sofr = self.build_trs(self.sofr_idx, start, end)
        trs_vs_sofr.setPricingEngine(pricer)

        par_trs_vs_ibor = self.build_trs(
            self.ibor_idx, start, end, trs_vs_ibor.fairMargin())
        par_trs_vs_ibor.setPricingEngine(pricer)
        par_trs_vs_sofr = self.build_trs(
            self.sofr_idx, start, end, trs_vs_sofr.fairMargin())
        par_trs_vs_sofr.setPricingEngine(pricer)

        fail_msg = "Par TRS expected to have NPV equal to zero."

        self.assertAlmostEqual(
            par_trs_vs_ibor.NPV(), 0.0, delta=EPSILON, msg=fail_msg)
        self.assertAlmostEqual(
            par_trs_vs_sofr.NPV(), 0.0, delta=EPSILON, msg=fail_msg)


if __name__ == '__main__':
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
