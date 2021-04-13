"""
 Copyright (C) 2021 Marcin Rybacki

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


EPSILON = 1.e-9

CAL = ql.TARGET()

VALUATION_DATE = CAL.adjust(ql.Date(15, ql.March, 2021))

RATE_AVERAGING_MAP = {ql.RateAveraging.Compound: 'Compounded',
                      ql.RateAveraging.Simple: 'Simple'}


def flat_rate(rate):
    return ql.FlatForward(
        2, CAL, ql.QuoteHandle(ql.SimpleQuote(rate)), ql.Actual365Fixed())


def create_ibor_leg(ibor_idx, start, end, payment_lag=0):
    day_count = ibor_idx.dayCounter()
    bdc = ibor_idx.businessDayConvention()
    cal = ibor_idx.fixingCalendar()
    sch = ql.MakeSchedule(effectiveDate=start,
                          terminationDate=end,
                          tenor=ibor_idx.tenor(),
                          calendar=cal,
                          convention=bdc,
                          backwards=True)
    return ql.IborLeg([1.0], sch, ibor_idx, day_count, bdc, cal, payment_lag)


def create_overnight_leg(overnight_idx, start, end, payment_lag=0):
    day_count = ql.Actual365Fixed()
    sch = ql.MakeSchedule(effectiveDate=start,
                          terminationDate=end,
                          tenor=ql.Period(1, ql.Years),
                          calendar=CAL,
                          convention=ql.Following,
                          backwards=True)
    return ql.OvernightLeg(
        [1.0], sch, overnight_idx, day_count, ql.Following, CAL, payment_lag)


def create_fixed_rate_leg(start, end, payment_lag=0):
    day_count = ql.Actual365Fixed()
    sch = ql.MakeSchedule(effectiveDate=start,
                          terminationDate=end,
                          tenor=ql.Period(1, ql.Years),
                          calendar=CAL,
                          convention=ql.Following,
                          backwards=True)
    return ql.FixedRateLeg(
        sch, day_count, [1.0], [0.005], ql.Following, CAL, payment_lag)


class IborCouponTest(unittest.TestCase):
    def setUp(self):
        ql.Settings.instance().evaluationDate = VALUATION_DATE
        self.nominal_ts_handle = ql.YieldTermStructureHandle(flat_rate(0.007))
        self.ibor_idx = ql.Euribor6M(self.nominal_ts_handle)

    def test_payment_lag(self):
        """Testing payment lag of an Ibor leg"""
        start = ql.Date(17, ql.March, 2021)
        end = ql.Date(17, ql.March, 2031)
        pay_lag = 2
        leg_without_lag = create_ibor_leg(self.ibor_idx, start, end)
        leg_with_lag = create_ibor_leg(self.ibor_idx, start, end, pay_lag)
        for c_f_without_lag, c_f_with_lag in zip(leg_without_lag, leg_with_lag):
            actual_payment_date = c_f_with_lag.date()
            expected_payment_date = CAL.advance(
                c_f_without_lag.date(),
                pay_lag,
                ql.Days,
                self.ibor_idx.businessDayConvention())

            fail_msg = """ Unable to replicate Ibor coupon payment date:
                            calculated: {actual}
                            expected: {expected}
                       """.format(actual=actual_payment_date,
                                  expected=expected_payment_date)
            self.assertTrue(actual_payment_date == expected_payment_date,
                            msg=fail_msg)


class OvernightCouponTest(unittest.TestCase):
    def setUp(self):
        ql.Settings.instance().evaluationDate = VALUATION_DATE
        self.nominal_ts_handle = ql.YieldTermStructureHandle(flat_rate(0.007))
        self.overnight_idx = ql.Eonia(self.nominal_ts_handle)

    def test_payment_lag(self):
        """Testing payment lag of an overnight leg"""
        start = ql.Date(17, ql.March, 2021)
        end = ql.Date(17, ql.March, 2031)
        pay_lag = 2
        leg_without_lag = create_overnight_leg(self.overnight_idx, start, end)
        leg_with_lag = create_overnight_leg(
            self.overnight_idx, start, end, pay_lag)
        for c_f_without_lag, c_f_with_lag in zip(leg_without_lag, leg_with_lag):
            actual_payment_date = c_f_with_lag.date()
            expected_payment_date = CAL.advance(
                c_f_without_lag.date(), pay_lag, ql.Days, ql.Following)

            fail_msg = """ Unable to replicate overnight coupon payment date:
                            calculated: {actual}
                            expected: {expected}
                       """.format(actual=actual_payment_date,
                                  expected=expected_payment_date)
            self.assertTrue(actual_payment_date == expected_payment_date,
                            msg=fail_msg)


class FixedRateCouponTest(unittest.TestCase):
    def setUp(self):
        ql.Settings.instance().evaluationDate = VALUATION_DATE
        self.nominal_ts_handle = ql.YieldTermStructureHandle(flat_rate(0.007))

    def test_payment_lag(self):
        """Testing payment lag of a fixed rate leg"""
        start = ql.Date(17, ql.March, 2021)
        end = ql.Date(17, ql.March, 2031)
        pay_lag = 2
        leg_without_lag = create_fixed_rate_leg(start, end)
        leg_with_lag = create_fixed_rate_leg(start, end, pay_lag)
        for c_f_without_lag, c_f_with_lag in zip(leg_without_lag, leg_with_lag):
            actual_payment_date = c_f_with_lag.date()
            expected_payment_date = CAL.advance(
                c_f_without_lag.date(), pay_lag, ql.Days, ql.Following)

            fail_msg = """ Unable to replicate fixed rate coupon payment date:
                            calculated: {actual}
                            expected: {expected}
                       """.format(actual=actual_payment_date,
                                  expected=expected_payment_date)
            self.assertTrue(actual_payment_date == expected_payment_date,
                            msg=fail_msg)


def create_sub_periods_coupon(ibor_idx, start, end, averaging_method):
    payment_calendar = ibor_idx.fixingCalendar()
    payment_bdc = ibor_idx.businessDayConvention()
    payment_date = payment_calendar.adjust(end, payment_bdc)
    fixing_delay = ibor_idx.fixingDays()
    cpn = ql.SubPeriodsCoupon(
        payment_date, 1.0, start, end, fixing_delay, ibor_idx)
    use_compounded_rate = (averaging_method == ql.RateAveraging.Compound)
    if use_compounded_rate:
        cpn.setPricer(ql.CompoundingRatePricer())
    else:
        cpn.setPricer(ql.AveragingRatePricer())
    return cpn


def create_sub_periods_leg(
        ibor_idx, start, end, cpn_frequency, averaging_method):
    sch = ql.MakeSchedule(effectiveDate=start,
                          terminationDate=end,
                          tenor=cpn_frequency,
                          calendar=ibor_idx.fixingCalendar(),
                          convention=ibor_idx.businessDayConvention(),
                          backwards=True)
    day_count = ibor_idx.dayCounter()
    return ql.SubPeriodsLeg(
        [1.0],
        sch,
        ibor_idx,
        day_count,
        ql.Following,
        CAL,
        0,
        [2],
        [1.0],
        [0.0],
        [0.0],
        ql.Period(),
        CAL,
        ql.Unadjusted,
        False,
        averaging_method)


def sum_leg_payments(leg):
    return sum([cf.amount() for cf in leg])


def compounded_leg_payment(leg):
    compound = 1.0
    for cf in leg:
        floating_cf = ql.as_floating_rate_coupon(cf)
        year_fraction = floating_cf.accrualPeriod()
        fixing = floating_cf.indexFixing()
        compound *= (1.0 + year_fraction * fixing)
    return compound - 1.0


def averaged_leg_payment(leg):
    acc = 0.0
    for cf in leg:
        floating_cf = ql.as_floating_rate_coupon(cf)
        year_fraction = floating_cf.accrualPeriod()
        fixing = floating_cf.indexFixing()
        acc += year_fraction * fixing
    return acc


class SubPeriodsCouponTest(unittest.TestCase):
    def setUp(self):
        ql.Settings.instance().evaluationDate = VALUATION_DATE
        self.nominal_ts_handle = ql.YieldTermStructureHandle(flat_rate(0.007))
        self.ibor_idx = ql.Euribor6M(self.nominal_ts_handle)
        self.ibor_idx.addFixing(ql.Date(10, ql.February, 2021), 0.0085)

    def check_single_period_coupon_replication(self, start, end, averaging):
        ibor_leg = create_ibor_leg(self.ibor_idx, start, end)
        sub_periods_cpn = create_sub_periods_coupon(
            self.ibor_idx, start, end, averaging)

        actual_payment = sub_periods_cpn.amount()
        expected_payment = sum_leg_payments(ibor_leg)

        fail_msg = """ Unable to replicate single period coupon payment:
                            calculated: {actual}
                            expected: {expected}
                            start: {start}
                            end: {end}
                   """.format(actual=actual_payment,
                              expected=expected_payment,
                              start=start,
                              end=end)
        self.assertTrue(
            abs(actual_payment - expected_payment) < EPSILON,
            msg=fail_msg)

    def check_multiple_compounded_sub_periods_coupon_replication(
            self, start, end):
        ibor_leg = create_ibor_leg(self.ibor_idx, start, end)
        sub_periods_cpn = create_sub_periods_coupon(
            self.ibor_idx, start, end, ql.RateAveraging.Compound)

        actual_payment = sub_periods_cpn.amount()
        expected_payment = compounded_leg_payment(ibor_leg)

        fail_msg = """ Unable to replicate compounded multiple sub-period coupon payment:
                            calculated: {actual}
                            expected: {expected}
                            start: {start}
                            end: {end}
                   """.format(actual=actual_payment,
                              expected=expected_payment,
                              start=start,
                              end=end)
        self.assertTrue(
            abs(actual_payment - expected_payment) < EPSILON,
            msg=fail_msg)

    def check_multiple_averaged_sub_periods_coupon_replication(
            self, start, end):
        ibor_leg = create_ibor_leg(self.ibor_idx, start, end)
        sub_periods_cpn = create_sub_periods_coupon(
            self.ibor_idx, start, end, ql.RateAveraging.Simple)

        actual_payment = sub_periods_cpn.amount()
        expected_payment = averaged_leg_payment(ibor_leg)

        fail_msg = """ Unable to replicate averaged multiple sub-period coupon payment:
                            calculated: {actual}
                            expected: {expected}
                            start: {start}
                            end: {end}
                   """.format(actual=actual_payment,
                              expected=expected_payment,
                              start=start,
                              end=end)
        self.assertTrue(
            abs(actual_payment - expected_payment) < EPSILON,
            msg=fail_msg)

    def check_sub_periods_leg_replication(self, averaging_method):
        start = ql.Date(18, ql.March, 2021)
        end = ql.Date(18, ql.March, 2022)

        sub_periods_cpn = create_sub_periods_coupon(
            self.ibor_idx, start, end, averaging_method)
        sub_periods_leg = create_sub_periods_leg(
            self.ibor_idx, start, end, ql.Period(1, ql.Years), averaging_method)

        actual_payment = sum_leg_payments(sub_periods_leg)
        expected_payment = sub_periods_cpn.amount()

        fail_msg = """ Unable to replicate sub-period leg payments:
                            calculated: {actual}
                            expected: {expected}
                            averaging: {averaging}
                   """.format(actual=actual_payment,
                              expected=expected_payment,
                              averaging=RATE_AVERAGING_MAP[averaging_method])
        self.assertTrue(
            abs(actual_payment - expected_payment) < EPSILON,
            msg=fail_msg)

    def test_regular_single_period_forward_starting_coupon(self):
        """Testing regular single period forward starting coupon"""
        start = ql.Date(15, ql.April, 2021)
        end = ql.Date(15, ql.October, 2021)

        self.check_single_period_coupon_replication(
            start, end, ql.RateAveraging.Simple)
        self.check_single_period_coupon_replication(
            start, end, ql.RateAveraging.Compound)

    def test_regular_single_period_coupon_after_fixing(self):
        """Testing regular single period coupon after fixing"""
        start = ql.Date(12, ql.February, 2021)
        end = ql.Date(12, ql.August, 2021)

        self.check_single_period_coupon_replication(
            start, end, ql.RateAveraging.Simple)
        self.check_single_period_coupon_replication(
            start, end, ql.RateAveraging.Compound)

    def test_irregular_single_period_coupon_after_fixing(self):
        """Testing irregular single period coupon after fixing"""
        start = ql.Date(12, ql.February, 2021)
        end = ql.Date(12, ql.June, 2021)

        self.check_single_period_coupon_replication(
            start, end, ql.RateAveraging.Simple)
        self.check_single_period_coupon_replication(
            start, end, ql.RateAveraging.Compound)

    def test_regular_compounded_forward_starting_coupon_with_multiple_sub_periods(self):
        """Testing regular forward starting coupon with multiple compounded sub-periods"""
        start = ql.Date(15, ql.April, 2021)
        end = ql.Date(15, ql.April, 2022)

        self.check_multiple_compounded_sub_periods_coupon_replication(
            start, end)

    def test_regular_averaged_forward_starting_coupon_with_multiple_sub_periods(self):
        """Testing regular forward starting coupon with multiple averaged sub-periods"""
        start = ql.Date(15, ql.April, 2021)
        end = ql.Date(15, ql.April, 2022)

        self.check_multiple_averaged_sub_periods_coupon_replication(start, end)

    def test_sub_periods_leg_cash_flows(self):
        """Testing sub-periods leg replication"""
        self.check_sub_periods_leg_replication(ql.RateAveraging.Compound)
        self.check_sub_periods_leg_replication(ql.RateAveraging.Simple)

    def test_casting(self):
        """Testing casting to sub periods coupon"""
        start = ql.Date(18, ql.March, 2021)
        end = ql.Date(18, ql.March, 2022)
        sub_periods_leg = create_sub_periods_leg(
            self.ibor_idx, start, end, ql.Period(1, ql.Years), ql.RateAveraging.Compound)
        cf = sub_periods_leg[0]
        self.assertTrue(not isinstance(cf, ql.SubPeriodsCoupon))
        self.assertTrue(isinstance(
            ql.as_sub_periods_coupon(cf), ql.SubPeriodsCoupon))


if __name__ == '__main__':
    print('testing QuantLib ' + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(SubPeriodsCouponTest, 'test'))
    suite.addTest(unittest.makeSuite(IborCouponTest, 'test'))
    suite.addTest(unittest.makeSuite(OvernightCouponTest, 'test'))
    suite.addTest(unittest.makeSuite(FixedRateCouponTest, 'test'))
    unittest.TextTestRunner(verbosity=2).run(suite)
