"""
 Copyright (C) 2020 Marcin Rybacki

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

from typing import List, Tuple


UK_NOMINAL_DATA = [
    (ql.Date(26, ql.November, 2009), 0.475),
    (ql.Date(2, ql.December, 2009), 0.47498),
    (ql.Date(29, ql.December, 2009), 0.49988),
    (ql.Date(25, ql.February, 2010), 0.59955),
    (ql.Date(18, ql.March, 2010), 0.65361),
    (ql.Date(25, ql.May, 2010), 0.82830),
    (ql.Date(16, ql.September, 2010), 0.78960),
    (ql.Date(16, ql.December, 2010), 0.93762),
    (ql.Date(17, ql.March, 2011), 1.12037),
    (ql.Date(16, ql.June, 2011), 1.31308),
    (ql.Date(22, ql.September, 2011), 1.52011),
    (ql.Date(25, ql.November, 2011), 1.78399),
    (ql.Date(26, ql.November, 2012), 2.41170),
    (ql.Date(25, ql.November, 2013), 2.83935),
    (ql.Date(25, ql.November, 2014), 3.12888),
    (ql.Date(25, ql.November, 2015), 3.34298),
    (ql.Date(25, ql.November, 2016), 3.50632),
    (ql.Date(27, ql.November, 2017), 3.63666),
    (ql.Date(26, ql.November, 2018), 3.74723),
    (ql.Date(25, ql.November, 2019), 3.83988),
    (ql.Date(25, ql.November, 2021), 4.00508),
    (ql.Date(25, ql.November, 2024), 4.16042),
    (ql.Date(26, ql.November, 2029), 4.15577),
    (ql.Date(27, ql.November, 2034), 4.04933),
    (ql.Date(25, ql.November, 2039), 3.95217),
    (ql.Date(25, ql.November, 2049), 3.80932),
    (ql.Date(25, ql.November, 2059), 3.80849),
    (ql.Date(25, ql.November, 2069), 3.72677),
    (ql.Date(27, ql.November, 2079), 3.63082)]

EU_NOMINAL_DATA = [(ql.Date(28, ql.October, 2018), 0.475),
                   (ql.Date(28, ql.October, 2019), 0.47498),
                   (ql.Date(28, ql.October, 2020), 0.49988)]

UK_FIXING_DATA = [
    (ql.Date(20, ql.July, 2007), 206.1),
    (ql.Date(20, ql.August, 2007), 207.3),
    (ql.Date(20, ql.September, 2007), 208.0),
    (ql.Date(22, ql.October, 2007), 208.9),
    (ql.Date(20, ql.November, 2007), 209.7),
    (ql.Date(20, ql.December, 2007), 210.9),
    (ql.Date(21, ql.January, 2008), 209.8),
    (ql.Date(20, ql.February, 2008), 211.4),
    (ql.Date(20, ql.March, 2008), 212.1),
    (ql.Date(21, ql.April, 2008), 214.0),
    (ql.Date(20, ql.May, 2008), 215.1),
    (ql.Date(21, ql.June, 2008), 216.8),
    (ql.Date(20, ql.August, 2008), 216.5),
    (ql.Date(22, ql.September, 2008), 217.2),
    (ql.Date(20, ql.October, 2008), 218.4),
    (ql.Date(20, ql.November, 2008), 217.7),
    (ql.Date(22, ql.December, 2008), 216.0),
    (ql.Date(20, ql.January, 2009), 212.9),
    (ql.Date(20, ql.February, 2009), 210.1),
    (ql.Date(20, ql.March, 2009), 211.4),
    (ql.Date(20, ql.April, 2009), 211.3),
    (ql.Date(20, ql.May, 2009), 211.5),
    (ql.Date(22, ql.June, 2009), 212.8),
    (ql.Date(20, ql.July, 2009), 213.4),
    (ql.Date(20, ql.August, 2009), 213.4),
    (ql.Date(21, ql.September, 2009), 213.4)]

EU_FIXING_DATA = [(ql.Date(1, ql.June, 2018), 103.76),
                  (ql.Date(1, ql.July, 2018), 103.41),
                  (ql.Date(1, ql.August, 2018), 103.58)]

EU_ZERO_COUPON_DATA = [(ql.Date(28, ql.September, 2019), 3.0495),
                       (ql.Date(28, ql.September, 2020), 2.93),
                       (ql.Date(28, ql.September, 2021), 2.9795)]

UK_ZERO_COUPON_DATA = [
    (ql.Date(25, ql.November, 2010), 3.0495),
    (ql.Date(25, ql.November, 2011), 2.93),
    (ql.Date(26, ql.November, 2012), 2.9795),
    (ql.Date(25, ql.November, 2013), 3.029),
    (ql.Date(25, ql.November, 2014), 3.1425),
    (ql.Date(25, ql.November, 2015), 3.211),
    (ql.Date(25, ql.November, 2016), 3.2675),
    (ql.Date(25, ql.November, 2017), 3.3625),
    (ql.Date(25, ql.November, 2018), 3.405),
    (ql.Date(25, ql.November, 2019), 3.48),
    (ql.Date(25, ql.November, 2021), 3.576),
    (ql.Date(25, ql.November, 2024), 3.649),
    (ql.Date(26, ql.November, 2029), 3.751),
    (ql.Date(27, ql.November, 2034), 3.77225),
    (ql.Date(25, ql.November, 2039), 3.77),
    (ql.Date(25, ql.November, 2049), 3.734),
    (ql.Date(25, ql.November, 2059), 3.714)]


def create_inflation_helper(
    inflation_data: Tuple[ql.Date, float],
    inflation_index: ql.ZeroInflationIndex,
    observation_lag: ql.Period,
    calendar: ql.Calendar,
    business_day_convention,
    day_counter: ql.DayCounter,
    discount_curve_handle: ql.YieldTermStructureHandle
) -> ql.ZeroCouponInflationSwapHelper:
    maturity = inflation_data[0]
    quote = ql.QuoteHandle(ql.SimpleQuote(inflation_data[1] / 100.0))
    return ql.ZeroCouponInflationSwapHelper(
        quote,
        observation_lag,
        maturity,
        calendar,
        business_day_convention,
        day_counter,
        inflation_index,
        discount_curve_handle)


def build_nominal_term_structure(
        nominal_data: List[Tuple[ql.Date, float]]) -> ql.YieldTermStructure:
    nominal_dc = ql.Actual365Fixed()
    dates = [x[0] for x in nominal_data]
    rates = [x[1] for x in nominal_data]
    return ql.ZeroCurve(dates, rates, nominal_dc)


def build_ukrpi_index(
        fixing_data: List[Tuple[ql.Date, float]],
        inflation_crv_handle: ql.ZeroInflationTermStructureHandle,
        interpolated: bool = False) -> ql.ZeroInflationIndex:
    # uk rpi index fixing data
    index = ql.UKRPI(interpolated, inflation_crv_handle)
    for x in fixing_data:
        # force override in case of multiple use
        index.addFixing(x[0], x[1], True)

    return index


def build_hicp_index(
        fixing_data: List[Tuple[ql.Date, float]],
        inflation_crv_handle: ql.ZeroInflationTermStructureHandle,
        interpolated: bool = False) -> ql.ZeroInflationIndex:
    # uk rpi index fixing data
    index = ql.EUHICP(interpolated, inflation_crv_handle)
    for x in fixing_data:
        index.addFixing(x[0], x[1], True)

    return index


SEASONAL = {ql.January: 1.0, ql.February: 1.01, ql.March: 1.011,
            ql.April: 1.009, ql.May: 1.008, ql.June: 1.012,
            ql.July: 1.0078, ql.August: 1.006,
            ql.September: 1.0085, ql.October: 1.0096,
            ql.November: 1.0067, ql.December: 1.0055}


def construct_seasonality(evaluation_date: ql.Date) -> ql.Seasonality:
    frequency = ql.Monthly
    seasonality_base_date = ql.Date(1, ql.January, evaluation_date.year())
    factors = list(SEASONAL.values())
    return ql.MultiplicativePriceSeasonality(
        seasonality_base_date, frequency, factors)


def get_seasonality_factor(d: ql.Date):
    return SEASONAL[d.month()]


def build_inflation_term_structure(
        reference_date: ql.Date,
        zero_coupon_data: List[Tuple[ql.Date, float]],
        inflation_index: ql.ZeroInflationIndex,
        nominal_term_structure_handle: ql.YieldTermStructureHandle,
        observation_lag: ql.Period,
        include_seasonality: bool = False) -> ql.ZeroInflationTermStructure:
    calendar = ql.TARGET()
    payment_convention = ql.ModifiedFollowing
    day_counter = ql.ActualActual()
    helpers = [create_inflation_helper(x,
                                       inflation_index,
                                       observation_lag,
                                       calendar,
                                       payment_convention,
                                       day_counter,
                                       nominal_term_structure_handle)
               for x in zero_coupon_data]
    base_zero_rate = zero_coupon_data[0][1] / 100.0
    cpi_term_structure = ql.PiecewiseZeroInflation(
        reference_date,
        calendar,
        day_counter,
        observation_lag,
        inflation_index.frequency(),
        inflation_index.interpolated(),
        base_zero_rate,
        nominal_term_structure_handle,
        helpers)
    if include_seasonality:
        seasonality = construct_seasonality(reference_date)
        cpi_term_structure.setSeasonality(seasonality)
    return cpi_term_structure


def create_inflation_swap(
        index: ql.ZeroInflationIndex,
        start_date: ql.Date,
        end_date: ql.Date,
        rate: float,
        observation_lag: ql.Period) -> ql.ZeroCouponInflationSwap:
    payer = ql.ZeroCouponInflationSwap.Payer
    calendar = ql.UnitedKingdom()
    payment_convention = ql.ModifiedFollowing
    day_counter = ql.ActualActual()
    nominal = 1.e6
    return ql.ZeroCouponInflationSwap(payer, nominal, start_date,
                                      end_date, calendar, payment_convention,
                                      day_counter, rate, index,
                                      observation_lag)


class InflationTest(unittest.TestCase):
    def test_fom_indexation_without_seasonality(self):
        """Testing first of month indexation without seasonality"""

        today = ql.Date(25, ql.November, 2009)
        ql.Settings.instance().setEvaluationDate(today)

        # Nominal term structure handle
        nominal_ts_handle = ql.RelinkableYieldTermStructureHandle()
        nominal_ts = build_nominal_term_structure(UK_NOMINAL_DATA)
        nominal_ts_handle.linkTo(nominal_ts)
        discount_engine = ql.DiscountingSwapEngine(nominal_ts_handle)

        # Inflation curve handle
        inflation_ts_handle = ql.RelinkableZeroInflationTermStructureHandle()
        inflation_index = build_ukrpi_index(
            UK_FIXING_DATA, inflation_ts_handle)
        observation_lag = ql.Period(2, ql.Months)
        inflation_ts = build_inflation_term_structure(
            today,
            UK_ZERO_COUPON_DATA,
            inflation_index,
            nominal_ts_handle,
            observation_lag)
        inflation_ts_handle.linkTo(inflation_ts)

        # Create par inflation swap
        zciis = create_inflation_swap(
            inflation_index, today, ql.Date(25, ql.November, 2059),
            0.03714, observation_lag)
        zciis.setPricingEngine(discount_engine)

        # Check whether swap prices to par
        self.assertTrue(
            abs(zciis.NPV() < 1.e-5),
            msg="Failed to price zero coupon inflation swap to par.")

        inflation_cf = ql.as_indexed_cashflow(
            zciis.inflationLeg()[0])

        # Obtaining base index for the inflation swap
        swap_base_d = inflation_cf.baseDate()
        swap_base_index = inflation_index.fixing(swap_base_d)

        # Replicate fixing projection
        fixing_d = inflation_cf.fixingDate()
        ts_base_d = inflation_ts.baseDate()
        ts_base_index = inflation_index.fixing(ts_base_d)

        # Apply FOM indexation rule
        effective_fixing_d = ql.Date(
            1, fixing_d.month(), fixing_d.year())
        fraction = inflation_ts.dayCounter().yearFraction(
            ts_base_d, effective_fixing_d)
        t = inflation_ts.timeFromReference(effective_fixing_d)
        zero_rate = inflation_ts.zeroRate(t)

        expected_fixing = ts_base_index * (
            1.0 + zero_rate)**fraction

        expected_inf_leg_payment = (
            expected_fixing / swap_base_index - 1.0) * inflation_cf.notional()
        actual_inf_leg_payment = inflation_cf.amount()
        self.assertAlmostEquals(
            first=actual_inf_leg_payment,
            second=expected_inf_leg_payment,
            delta=1.e-10,
            msg="Failed to replicate the inflation leg expected flow.")

    def test_linear_indexation_without_seasonality(self):
        """Testing linear indexation without seasonality"""

        today = ql.Date(25, ql.November, 2009)
        ql.Settings.instance().setEvaluationDate(today)

        # Nominal term structure handle
        nominal_ts_handle = ql.RelinkableYieldTermStructureHandle()
        nominal_ts = build_nominal_term_structure(UK_NOMINAL_DATA)
        nominal_ts_handle.linkTo(nominal_ts)
        discount_engine = ql.DiscountingSwapEngine(nominal_ts_handle)

        # Inflation curve handle
        inflation_ts_handle = ql.RelinkableZeroInflationTermStructureHandle()
        index = build_ukrpi_index(
            UK_FIXING_DATA, inflation_ts_handle, interpolated=True)
        observation_lag = ql.Period(3, ql.Months)
        inflation_ts = build_inflation_term_structure(
            today, UK_ZERO_COUPON_DATA, index, nominal_ts_handle, observation_lag)
        inflation_ts_handle.linkTo(inflation_ts)

        # Create par inflation swap
        zciis = create_inflation_swap(
            index,
            ql.Date(25, ql.April, 2009),
            ql.Date(25, ql.November, 2059),
            0.038,
            observation_lag)
        zciis.setPricingEngine(discount_engine)

        inflation_cf = ql.as_indexed_cashflow(
            zciis.inflationLeg()[0])

        # Replicate base index for the inflation swap
        def interpolate_historic_index(fixing_date: ql.Date):
            f_d = ql.Date(1, fixing_date.month(), fixing_date.year())
            s_d = ql.Date.endOfMonth(fixing_date) + 1
            slope = (fixing_date - f_d) / (
                (s_d + observation_lag) - (f_d + observation_lag))
            return index.fixing(f_d) + slope * (
                index.fixing(s_d) - index.fixing(f_d))

        swap_base_d = inflation_cf.baseDate()
        swap_base_index = index.fixing(swap_base_d)
        expected_swap_base_index = interpolate_historic_index(swap_base_d)
        self.assertAlmostEquals(
            first=swap_base_index,
            second=expected_swap_base_index,
            msg="Failed to replicate the base index on the inflation swap.")

        # Replicate projected swap fixing
        fixing_d = inflation_cf.fixingDate()

        ts_base_d = inflation_ts.baseDate()
        ts_base_index = index.fixing(ts_base_d)
        expected_ts_base_index = interpolate_historic_index(ts_base_d)
        self.assertAlmostEquals(
            first=ts_base_index,
            second=expected_ts_base_index,
            msg="Failed to replicate the base index on the inflation curve.")

        # Apply linear indexation rule
        fraction = inflation_ts.dayCounter().yearFraction(
            ts_base_d, fixing_d)
        t = inflation_ts.timeFromReference(fixing_d)
        zero_rate = inflation_ts.zeroRate(t)

        expected_fixing = ts_base_index * (
            1.0 + zero_rate)**fraction

        # Assert inflation leg projected amount
        expected_inf_leg_payment = (
            expected_fixing / swap_base_index - 1.0) * inflation_cf.notional()
        actual_inf_leg_payment = inflation_cf.amount()
        self.assertAlmostEquals(
            first=actual_inf_leg_payment,
            second=expected_inf_leg_payment,
            delta=1.e-10,
            msg="Failed to replicate the inflation leg expected flow.")

    def test_linear_indexation_with_seasonality(self):
        """Testing linear indexation with seasonality"""

        today = ql.Date(25, ql.November, 2009)
        ql.Settings.instance().setEvaluationDate(today)

        # Nominal term structure handle
        nominal_ts_handle = ql.RelinkableYieldTermStructureHandle()
        nominal_ts = build_nominal_term_structure(UK_NOMINAL_DATA)
        nominal_ts_handle.linkTo(nominal_ts)
        discount_engine = ql.DiscountingSwapEngine(nominal_ts_handle)

        # Inflation curve handle
        inflation_ts_handle = ql.RelinkableZeroInflationTermStructureHandle()
        inflation_index = build_ukrpi_index(
            UK_FIXING_DATA, inflation_ts_handle, interpolated=True)
        observation_lag = ql.Period(3, ql.Months)
        inflation_ts = build_inflation_term_structure(
            today,
            UK_ZERO_COUPON_DATA,
            inflation_index,
            nominal_ts_handle,
            observation_lag,
            include_seasonality=True)
        inflation_ts_handle.linkTo(inflation_ts)

        # Create par inflation swap
        zciis = create_inflation_swap(
            inflation_index,
            ql.Date(25, ql.April, 2009),
            ql.Date(25, ql.November, 2059),
            0.037,
            observation_lag)
        zciis.setPricingEngine(discount_engine)

        inflation_cf = ql.as_indexed_cashflow(
            zciis.inflationLeg()[0])

        # Obtaining base index for the inflation swap
        swap_base_d = inflation_cf.baseDate()
        swap_base_index = inflation_index.fixing(swap_base_d)

        # Replicate fixing projection
        fixing_d = inflation_cf.fixingDate()
        ts_base_d = inflation_ts.baseDate()
        ts_base_index = inflation_index.fixing(ts_base_d)

        # Apply linear indexation rule
        fraction = inflation_ts.dayCounter().yearFraction(
            ts_base_d, fixing_d)
        t = inflation_ts.timeFromReference(fixing_d)
        zero_rate = inflation_ts.zeroRate(t)

        # Calculate seasonality adjustment
        seasonality_b = get_seasonality_factor(ts_base_d)
        seasonality_f = get_seasonality_factor(fixing_d)

        expected_fixing = ts_base_index * (
            seasonality_f / seasonality_b) * (1.0 + zero_rate)**fraction

        expected_inf_leg_payment = (
            expected_fixing / swap_base_index - 1.0) * inflation_cf.notional()
        actual_inf_leg_payment = inflation_cf.amount()
        self.assertAlmostEquals(
            first=actual_inf_leg_payment,
            second=expected_inf_leg_payment,
            delta=1.e-10,
            msg="Failed to replicate the inflation leg expected flow.")


if __name__ == '__main__':
    print('testing QuantLib ' + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(InflationTest, 'test'))
    unittest.TextTestRunner(verbosity=2).run(suite)
