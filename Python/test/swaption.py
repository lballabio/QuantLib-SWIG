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


BASIS_POINT = 1e-4

EPSILON = 1.e-10

OPTION_TYPE_MAP = {ql.VanillaSwap.Receiver: 'Receiver',
                   ql.VanillaSwap.Payer: 'Payer'}


SETTLEMENT_TYPE_MAP = {ql.Settlement.Physical: 'Physical',
                       ql.Settlement.Cash: 'Cash'}


SETTLEMENT_METHOD_MAP = {ql.Settlement.PhysicalOTC: 'Physical OTC',
                         ql.Settlement.CollateralizedCashPrice: (
                             'Collateralized Cash Price'),
                         ql.Settlement.ParYieldCurve: 'Par Yield Curve'}


def compounded_annual_constant_rate_discount(
        rate: float,
        day_counter: ql.DayCounter):
    def _calc(start: ql.Date, end: ql.Date):
        time = day_counter.yearFraction(start, end)
        return (1.0 + rate) ** (-time)
    return _calc


def par_yield_bps(underlying: ql.VanillaSwap,
                  discount_handle: ql.YieldTermStructureHandle):
    fixed_leg = underlying.fixedLeg()
    first_coupon = ql.as_fixed_rate_coupon(fixed_leg[0])
    discount_date = first_coupon.accrualStartDate()
    discount = discount_handle.discount(discount_date)
    fixed_rate = underlying.fixedRate()
    fixed_dct = underlying.fixedDayCount()
    fair_rate = underlying.fairRate()
    ir_func = compounded_annual_constant_rate_discount(fair_rate, fixed_dct)
    bps = sum([ir_func(discount_date, c_f.date()) * c_f.amount() / fixed_rate
               for c_f in fixed_leg
               if c_f.date() > discount_date])
    return abs(bps) * discount


def swap_pv01(underlying: ql.VanillaSwap):
    return abs(underlying.fixedLegBPS()) / BASIS_POINT


def make_const_black_vol_engine(discount_handle, volatility):
    h = ql.QuoteHandle(ql.SimpleQuote(volatility))
    return ql.BlackSwaptionEngine(discount_handle, h)


def make_const_bachelier_vol_engine(discount_handle, volatility):
    h = ql.QuoteHandle(ql.SimpleQuote(volatility))
    return ql.BachelierSwaptionEngine(discount_handle, h)


class SwaptionTest(unittest.TestCase):
    def setUp(self):
        self.calendar = ql.TARGET()
        self.today = self.calendar.adjust(ql.Date.todaysDate())
        ql.Settings.instance().setEvaluationDate(self.today)

        projection_curve_handle = ql.RelinkableYieldTermStructureHandle()
        self.projection_rate = 0.01
        self.projection_quote_handle = ql.RelinkableQuoteHandle()
        projection_curve = ql.FlatForward(
            self.today, self.projection_quote_handle, ql.Actual365Fixed())
        projection_curve_handle.linkTo(projection_curve)

        self.discount_handle = ql.YieldTermStructureHandle(ql.FlatForward(
            self.today, ql.QuoteHandle(ql.SimpleQuote(0.0085)), ql.Actual365Fixed()))
        self.swap_engine = ql.DiscountingSwapEngine(self.discount_handle)

        self.idx = ql.Euribor6M(projection_curve_handle)

        self.exercises = [ql.Period(1, ql.Years), ql.Period(2, ql.Years),
                          ql.Period(3, ql.Years), ql.Period(5, ql.Years),
                          ql.Period(7, ql.Years), ql.Period(10, ql.Years)]
        self.lengths = [ql.Period(1, ql.Years), ql.Period(2, ql.Years),
                        ql.Period(3, ql.Years), ql.Period(5, ql.Years),
                        ql.Period(7, ql.Years), ql.Period(10, ql.Years),
                        ql.Period(15, ql.Years), ql.Period(20, ql.Years)]
        self.swap_type = [ql.VanillaSwap.Receiver, ql.VanillaSwap.Payer]

    def _assert_swaption_delta(self,
                               swaption_pricer_func,
                               use_bachelier_vol: bool):
        strikes = [0.03, 0.04, 0.05, 0.06, 0.07]
        vols = [0.01, 0.10, 0.20, 0.30, 0.70, 0.90]

        settle_map = {ql.Settlement.PhysicalOTC: ql.Settlement.Physical,
                      ql.Settlement.CollateralizedCashPrice: ql.Settlement.Cash}

        for v in vols:
            for e in self.exercises:
                for l in self.lengths:
                    for s in strikes:
                        for t in self.swap_type:
                            for s_m, s_t in settle_map.items():
                                volatility = v / 100.0 if use_bachelier_vol else v
                                swaption_engine = swaption_pricer_func(
                                    self.discount_handle, volatility)
                                exercise_date = self.calendar.advance(
                                    self.today, e)
                                start_date = self.calendar.advance(
                                    exercise_date, ql.Period(2, ql.Days))

                                self.projection_quote_handle.linkTo(
                                    ql.SimpleQuote(self.projection_rate))

                                underlying = ql.MakeVanillaSwap(
                                    l, self.idx, s, ql.Period(0, ql.Days),
                                    effectiveDate=start_date,
                                    fixedLegTenor=ql.Period(1, ql.Years),
                                    fixedLegDayCount=ql.Thirty360(),
                                    floatingLegSpread=0.0,
                                    swapType=t)
                                underlying.setPricingEngine(self.swap_engine)

                                fair_rate = underlying.fairRate()

                                swaption = ql.Swaption(underlying,
                                                       ql.EuropeanExercise(
                                                           exercise_date),
                                                       s_t,
                                                       s_m)
                                swaption.setPricingEngine(swaption_engine)

                                value = swaption.NPV()
                                delta = swaption.delta() * BASIS_POINT

                                self.projection_quote_handle.linkTo(
                                    ql.SimpleQuote(self.projection_rate + BASIS_POINT))

                                bumped_fair_rate = underlying.fairRate()
                                bumped_value = swaption.NPV()
                                bumped_delta = swaption.delta() * BASIS_POINT

                                delta_bump = bumped_fair_rate - fair_rate
                                approx_delta = (bumped_value - value) / \
                                    delta_bump * BASIS_POINT

                                lower_bound = min(
                                    delta, bumped_delta) - EPSILON
                                upper_bound = max(
                                    delta, bumped_delta) + EPSILON

                                # Based on the Mean Value Theorem, the below inequality
                                # should hold for any function that is monotonic in the
                                # area of the bump.
                                check_is_correct = (lower_bound < approx_delta) and (
                                    approx_delta < upper_bound)

                                fail_msg = f""" Swaption delta test failed for:
                                                    option tenor: {e}
                                                    volatility : {volatility}
                                                    option type: {OPTION_TYPE_MAP[t]}
                                                    swap tenor: {l}
                                                    strike: {s}
                                                    settlement: {SETTLEMENT_TYPE_MAP[s_t]}
                                                    method: {SETTLEMENT_METHOD_MAP[s_m]}
                                                    delta: {delta}
                                                    approx delta: {approx_delta}
                                            """
                                self.assertTrue(check_is_correct, msg=fail_msg)

    def _assert_swaption_annuity(self,
                                 swaption_pricer_func,
                                 use_bachelier_vol: bool):
        self.projection_quote_handle.linkTo(
            ql.SimpleQuote(self.projection_rate))

        settle_type = ql.Settlement.Cash
        methods = [ql.Settlement.ParYieldCurve,
                   ql.Settlement.CollateralizedCashPrice]

        for e in self.exercises:
            for l in self.lengths:
                for t in self.swap_type:
                    for m in methods:
                        volatility = 0.003 if use_bachelier_vol else 0.3
                        strike = 0.03
                        swaption_engine = swaption_pricer_func(
                            self.discount_handle, volatility)
                        exercise_date = self.calendar.advance(
                            self.today, e)
                        start_date = self.calendar.advance(
                            exercise_date, ql.Period(2, ql.Days))

                        underlying = ql.MakeVanillaSwap(
                            l, self.idx, strike, ql.Period(0, ql.Days),
                            effectiveDate=start_date,
                            fixedLegTenor=ql.Period(1, ql.Years),
                            fixedLegDayCount=ql.Thirty360(),
                            floatingLegSpread=0.0,
                            swapType=t)
                        underlying.setPricingEngine(self.swap_engine)

                        swaption = ql.Swaption(underlying,
                                               ql.EuropeanExercise(
                                                   exercise_date),
                                               settle_type,
                                               m)
                        swaption.setPricingEngine(swaption_engine)

                        annuity = swaption.annuity()
                        expected_annuity = 0.0
                        if (m == ql.Settlement.CollateralizedCashPrice):
                            expected_annuity = swap_pv01(underlying)
                        if (m == ql.Settlement.ParYieldCurve):
                            expected_annuity = par_yield_bps(
                                underlying, self.discount_handle)

                        fail_msg = f""" Swaption annuity test failed for:
                                            option tenor: {e}
                                            volatility : {volatility}
                                            option type: {OPTION_TYPE_MAP[t]}
                                            swap tenor: {l}
                                            strike: {strike}
                                            settlement: {SETTLEMENT_TYPE_MAP[settle_type]}
                                            method: {SETTLEMENT_METHOD_MAP[m]}
                                            annuity: {annuity}
                                            replicated annuity: {expected_annuity}
                                        """
                        self.assertAlmostEquals(
                            first=annuity,
                            second=expected_annuity,
                            delta=EPSILON,
                            msg=fail_msg)

    def test_swaption_delta_black_volatility(self):
        """Testing swaption delta in Black model"""
        self._assert_swaption_delta(
            swaption_pricer_func=make_const_black_vol_engine,
            use_bachelier_vol=False)

    def test_swaption_delta_bachelier_volatility(self):
        """Testing swaption delta in Bachelier model"""
        self._assert_swaption_delta(
            swaption_pricer_func=make_const_bachelier_vol_engine,
            use_bachelier_vol=True)

    def test_swaption_annuity_black_model(self):
        """Testing swaption annuity in Black model"""
        self._assert_swaption_annuity(
            swaption_pricer_func=make_const_black_vol_engine,
            use_bachelier_vol=False)

    def test_swaption_annuity_bachelier_model(self):
        """Testing swaption annuity in Bachelier model"""
        self._assert_swaption_annuity(
            swaption_pricer_func=make_const_bachelier_vol_engine,
            use_bachelier_vol=True)


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(SwaptionTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
