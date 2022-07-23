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

OPTION_TYPE_MAP = {ql.Swap.Receiver: 'Receiver',
                   ql.Swap.Payer: 'Payer'}

SETTLEMENT_TYPE_MAP = {ql.Settlement.Physical: 'Physical',
                       ql.Settlement.Cash: 'Cash'}

SETTLEMENT_METHOD_MAP = {ql.Settlement.PhysicalOTC: 'Physical OTC',
                         ql.Settlement.CollateralizedCashPrice: (
                             'Collateralized Cash Price'),
                         ql.Settlement.ParYieldCurve: 'Par Yield Curve'}


def compounded_annual_constant_rate_discount(
        rate,
        day_counter):
    def _calc(start, end):
        time = day_counter.yearFraction(start, end)
        return (1.0 + rate) ** (-time)
    return _calc


def par_yield_bps(underlying,
                  discount_handle):
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


def swap_pv01(underlying):
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
        ql.Settings.instance().evaluationDate = self.today

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
        self.swap_type = [ql.Swap.Receiver, ql.Swap.Payer]

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def _assert_swaption_annuity(self,
                                 swaption_pricer_func,
                                 use_bachelier_vol):
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
                            fixedLegDayCount=ql.Thirty360(ql.Thirty360.BondBasis),
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

                        fail_msg = """ Swaption annuity test failed for:
                                            option tenor: {option_tenor}
                                            volatility : {volatility}
                                            option type: {option_type}
                                            swap tenor: {swap_tenor}
                                            strike: {strike}
                                            settlement: {settle_type}
                                            method: {method}
                                            annuity: {annuity}
                                            replicated annuity: {expected_annuity}
                                   """.format(option_tenor=e,
                                              volatility=volatility,
                                              option_type=OPTION_TYPE_MAP[t],
                                              swap_tenor=l,
                                              strike=strike,
                                              settle_type=SETTLEMENT_TYPE_MAP[settle_type],
                                              method=SETTLEMENT_METHOD_MAP[m],
                                              annuity=annuity,
                                              expected_annuity=expected_annuity)
                        self.assertAlmostEquals(
                            first=annuity,
                            second=expected_annuity,
                            delta=EPSILON,
                            msg=fail_msg)

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
