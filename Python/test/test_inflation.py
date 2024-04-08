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


EPSILON = 1.e-9

# Hypothetical market data
EUR_ZERO_RATES = [(ql.Period(1, ql.Days), 0.0048),
                  (ql.Period(1, ql.Years), 0.0048),
                  (ql.Period(2, ql.Years), 0.00475),
                  (ql.Period(3, ql.Years), 0.005),
                  (ql.Period(5, ql.Years), 0.0055),
                  (ql.Period(10, ql.Years), 0.007)]

EUR_BEI_SWAP_RATES = [(ql.Period(1, ql.Years), 0.0301),
                      (ql.Period(2, ql.Years), 0.0299),
                      (ql.Period(3, ql.Years), 0.0305),
                      (ql.Period(5, ql.Years), 0.0315),
                      (ql.Period(10, ql.Years), 0.0355)]

# Source:
# https://ec.europa.eu/eurostat/web/products-datasets/-/teicp240.
EU_FIXING_DATA = [(ql.Date(1, ql.April, 2018), 103.11),
                  (ql.Date(1, ql.May, 2018), 103.64),
                  (ql.Date(1, ql.June, 2018), 103.76),
                  (ql.Date(1, ql.July, 2018), 103.41),
                  (ql.Date(1, ql.August, 2018), 103.58)]

CAL = ql.TARGET()

DAY_COUNTER = ql.ActualActual(ql.ActualActual.ISDA)

BDC = ql.ModifiedFollowing

VALUATION_DATE = CAL.adjust(ql.Date(10, ql.September, 2018))

OBSERVATION_LAG = ql.Period(3, ql.Months)


def create_inflation_swap_helper(
        reference_date,
        inflation_data,
        inflation_index,
        interpolation,
        discount_curve_handle,
        observation_lag=OBSERVATION_LAG,
        calendar=CAL,
        business_day_convention=BDC,
        day_counter=DAY_COUNTER):
    maturity = CAL.advance(reference_date, inflation_data[0])
    quote = ql.makeQuoteHandle(inflation_data[1])
    return ql.ZeroCouponInflationSwapHelper(
        quote,
        observation_lag,
        maturity,
        calendar,
        business_day_convention,
        day_counter,
        inflation_index,
        interpolation,
        discount_curve_handle)


def build_nominal_term_structure(
        reference_date,
        nominal_data):
    nominal_dc = ql.Actual365Fixed()
    dates = [CAL.advance(reference_date, x[0]) for x in nominal_data]
    rates = [x[1] for x in nominal_data]
    return ql.ZeroCurve(dates, rates, nominal_dc)


def build_hicp_index(
        fixing_data,
        inflation_crv_handle):
    index = ql.EUHICP(inflation_crv_handle)
    for x in fixing_data:
        # force override in case of multiple use
        index.addFixing(x[0], x[1], True)
    return index


SEASONAL = {ql.January: 1.0, ql.February: 1.01, ql.March: 1.011,
            ql.April: 1.009, ql.May: 1.008, ql.June: 1.012,
            ql.July: 1.0078, ql.August: 1.006,
            ql.September: 1.0085, ql.October: 1.0096,
            ql.November: 1.0067, ql.December: 1.0055}


def construct_seasonality(reference_date):
    frequency = ql.Monthly
    seasonality_base_date = ql.Date(1, ql.January, reference_date.year())
    factors = list(SEASONAL.values())
    return ql.MultiplicativePriceSeasonality(
        seasonality_base_date, frequency, factors)


def build_inflation_term_structure(
        reference_date,
        zero_coupon_swaps_data,
        inflation_index,
        interpolation,
        nominal_term_structure_handle,
        observation_lag=OBSERVATION_LAG,
        include_seasonality=False):
    helpers = [create_inflation_swap_helper(reference_date,
                                            x,
                                            inflation_index,
                                            interpolation,
                                            nominal_term_structure_handle)
               for x in zero_coupon_swaps_data]
    base_zero_rate = zero_coupon_swaps_data[0][1]
    cpi_term_structure = ql.PiecewiseZeroInflation(
        reference_date,
        inflation_index.lastFixingDate(),
        inflation_index.frequency(),
        DAY_COUNTER,
        helpers)
    if include_seasonality:
        seasonality = construct_seasonality(reference_date)
        cpi_term_structure.setSeasonality(seasonality)
    return cpi_term_structure


def create_inflation_swap(
        inflation_idx,
        start_date,
        end_date,
        rate,
        interpolation,
        observation_lag=OBSERVATION_LAG,
        nominal=1.e6,
        payer=ql.Swap.Payer):
    return ql.ZeroCouponInflationSwap(
        payer,
        nominal,
        start_date,
        end_date,
        CAL,
        BDC,
        DAY_COUNTER,
        rate,
        inflation_idx,
        observation_lag,
        interpolation)


def interpolate_historic_index(
        inflation_idx, fixing_date, observation_lag=OBSERVATION_LAG):
    first_dt = ql.Date(1, fixing_date.month(), fixing_date.year())
    second_dt = ql.Date.endOfMonth(fixing_date) + 1
    slope_numerator = fixing_date - first_dt
    slope_denominator = (
        (second_dt + observation_lag) - (first_dt + observation_lag))
    slope = float(slope_numerator) / float(slope_denominator)
    return inflation_idx.fixing(first_dt) + slope * (
        inflation_idx.fixing(second_dt) - inflation_idx.fixing(first_dt))


class InflationTest(unittest.TestCase):
    def setUp(self):
        ql.Settings.instance().evaluationDate = VALUATION_DATE
        self.inflation_ts_handle = ql.RelinkableZeroInflationTermStructureHandle()
        self.nominal_ts_handle = ql.RelinkableYieldTermStructureHandle()
        self.nominal_ts_handle.linkTo(
            build_nominal_term_structure(VALUATION_DATE, EUR_ZERO_RATES))
        self.discount_engine = ql.DiscountingSwapEngine(self.nominal_ts_handle)

    def test_par_swap_pricing_fom_indexation_without_seasonality(self):
        """Testing pricing of par inflation swap for First-Of-Month indexation"""

        inflation_idx = build_hicp_index(
            EU_FIXING_DATA, self.inflation_ts_handle)
        inflation_ts = build_inflation_term_structure(
            VALUATION_DATE,
            EUR_BEI_SWAP_RATES,
            inflation_idx,
            ql.CPI.Flat,
            self.nominal_ts_handle)
        self.inflation_ts_handle.linkTo(inflation_ts)

        zciis = create_inflation_swap(
            inflation_idx,
            VALUATION_DATE,
            CAL.advance(VALUATION_DATE, ql.Period(10, ql.Years)),
            0.0355,
            ql.CPI.Flat)
        zciis.setPricingEngine(self.discount_engine)
        npv = zciis.NPV()
        # Check whether swap prices to par
        fail_msg = """ Failed to price zero coupon inflation swap to par:
                            index: {inflation_idx}
                            end date: {end_date}
                            observation lag: {observation_lag}
                            npv: {npv}
                            expected npv: {expected_npv}
                            tolerance: {tolerance}
                   """.format(inflation_idx=inflation_idx.familyName(),
                              end_date=zciis.maturityDate(),
                              observation_lag=OBSERVATION_LAG,
                              npv=npv,
                              expected_npv=0.0,
                              tolerance=EPSILON)
        self.assertTrue(
            abs(npv) < EPSILON,
            msg=fail_msg)

    def test_inflation_leg_payment_fom_indexation_without_seasonality(self):
        """Testing inflation leg payment for First-Of-Month indexation"""

        inflation_idx = build_hicp_index(
            EU_FIXING_DATA, self.inflation_ts_handle)
        inflation_ts = build_inflation_term_structure(
            VALUATION_DATE,
            EUR_BEI_SWAP_RATES,
            inflation_idx,
            ql.CPI.Flat,
            self.nominal_ts_handle)
        self.inflation_ts_handle.linkTo(inflation_ts)

        zciis = create_inflation_swap(
            inflation_idx,
            VALUATION_DATE,
            CAL.advance(VALUATION_DATE, ql.Period(10, ql.Years)),
            0.0355,
            ql.CPI.Flat)
        zciis.setPricingEngine(self.discount_engine)

        inflation_cf = ql.as_indexed_cashflow(
            zciis.inflationLeg()[0])
        # Obtaining base index for the inflation swap
        swap_base_dt = inflation_cf.baseDate()
        swap_base_fixing = inflation_idx.fixing(swap_base_dt)
        # Replicate fixing projection
        fixing_dt = inflation_cf.fixingDate()
        ts_base_dt = inflation_ts.baseDate()
        ts_base_fixing = inflation_idx.fixing(ts_base_dt)
        # Apply FOM indexation rule
        effective_fixing_dt = ql.Date(
            1, fixing_dt.month(), fixing_dt.year())
        fraction = inflation_ts.dayCounter().yearFraction(
            ts_base_dt, effective_fixing_dt)
        t = inflation_ts.timeFromReference(effective_fixing_dt)
        zero_rate = inflation_ts.zeroRate(t)
        expected_fixing = ts_base_fixing * (
            1.0 + zero_rate)**fraction

        expected_inflation_leg_payment = (
            expected_fixing / swap_base_fixing - 1.0) * inflation_cf.notional()
        actual_inflation_leg_payment = inflation_cf.amount()

        fail_msg = """ Failed to replicate inflation leg payment
                       for First-Of-Month indexation:
                            index: {inflation_idx}
                            end date: {end_date}
                            observation lag: {observation_lag}
                            inflation leg payment: {actual_payment}
                            replicated payment: {expected_payment}
                            tolerance: {tolerance}
                   """.format(inflation_idx=inflation_idx.familyName(),
                              end_date=zciis.maturityDate(),
                              observation_lag=OBSERVATION_LAG,
                              actual_payment=actual_inflation_leg_payment,
                              expected_payment=expected_inflation_leg_payment,
                              tolerance=EPSILON)
        self.assertAlmostEqual(
            first=actual_inflation_leg_payment,
            second=expected_inflation_leg_payment,
            delta=EPSILON,
            msg=fail_msg)

    def test_lagged_fixing_method(self):
        """Testing lagged fixing method"""

        inflation_idx = build_hicp_index(
            EU_FIXING_DATA, self.inflation_ts_handle)
        inflation_ts = build_inflation_term_structure(
            VALUATION_DATE,
            EUR_BEI_SWAP_RATES,
            inflation_idx,
            ql.CPI.Flat,
            self.nominal_ts_handle)
        self.inflation_ts_handle.linkTo(inflation_ts)

        maturity_date = ql.Date(25, ql.October, 2027)
        lag = ql.Period(3, ql.Months)
        indexation = ql.CPI.Flat

        actual_fixing = ql.CPI.laggedFixing(inflation_idx, maturity_date, lag, indexation)
        expected_fixing = inflation_idx.fixing(ql.Date(1, ql.July, 2027))

        fail_msg = """ Failed to replicate lagged fixing:
                            index: {inflation_idx}
                            actual fixing: {actual_fixing}
                            expected fixing: {expected_fixing}
                            tolerance: {tolerance}
                   """.format(inflation_idx=inflation_idx.familyName(),
                              actual_fixing=actual_fixing,
                              expected_fixing=expected_fixing,
                              tolerance=EPSILON)
        self.assertAlmostEqual(
            first=actual_fixing,
            second=expected_fixing,
            msg=fail_msg,
            delta=EPSILON)


if __name__ == '__main__':
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
