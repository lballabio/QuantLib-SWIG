"""
 Copyright (C) 2020 tbd

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


def make_const_black_vol_engine(discount_handle, volatility):
    h = ql.QuoteHandle(ql.SimpleQuote(volatility))
    return ql.BlackSwaptionEngine(discount_handle, h)


def make_const_bachelier_vol_engine(discount_handle, volatility):
    h = ql.QuoteHandle(ql.SimpleQuote(volatility))
    return ql.BachelierSwaptionEngine(discount_handle, h)


class SwaptionTest(unittest.TestCase):

    def _assert_swaption_delta(self,
                               swaption_pricer_func,
                               use_bachelier_vol: bool):
        calendar = ql.TARGET()
        today = calendar.adjust(ql.Date.todaysDate())
        ql.Settings.instance().setEvaluationDate(today)

        bump = 1.e-4
        epsilon = 1.e-10

        projection_curve_handle = ql.RelinkableYieldTermStructureHandle()
        projection_rate = 0.01
        projection_quote_handle = ql.RelinkableQuoteHandle()
        projection_curve = ql.FlatForward(
            today, projection_quote_handle, ql.Actual365Fixed())
        projection_curve_handle.linkTo(projection_curve)

        discount_handle = ql.YieldTermStructureHandle(ql.FlatForward(
            today, ql.QuoteHandle(ql.SimpleQuote(0.0085)), ql.Actual365Fixed()))
        swap_engine = ql.DiscountingSwapEngine(discount_handle)

        idx = ql.Euribor6M(projection_curve_handle)

        exercises = [ql.Period(1, ql.Years), ql.Period(2, ql.Years),
                     ql.Period(3, ql.Years), ql.Period(5, ql.Years),
                     ql.Period(7, ql.Years), ql.Period(10, ql.Years)]
        lengths = [ql.Period(1, ql.Years), ql.Period(2, ql.Years),
                   ql.Period(3, ql.Years), ql.Period(5, ql.Years),
                   ql.Period(7, ql.Years), ql.Period(10, ql.Years),
                   ql.Period(15, ql.Years), ql.Period(20, ql.Years)]
        swap_type = [ql.VanillaSwap.Receiver, ql.VanillaSwap.Payer]

        types = [ql.Settlement.Physical, ql.Settlement.Cash]
        methods = [ql.Settlement.PhysicalOTC,
                   ql.Settlement.CollateralizedCashPrice]

        strikes = [0.03, 0.04, 0.05, 0.06, 0.07]
        vols = [0.01, 0.10, 0.20, 0.30, 0.70, 0.90]

        for v in vols:
            for e in exercises:
                for l in lengths:
                    for s in strikes:
                        for h in range(2):
                            volatility = v / 100.0 if use_bachelier_vol else v
                            swaption_engine = swaption_pricer_func(
                                discount_handle, volatility)
                            exercise_date = calendar.advance(today, e)
                            start_date = calendar.advance(
                                exercise_date, ql.Period(2, ql.Days))

                            projection_quote_handle.linkTo(
                                ql.SimpleQuote(projection_rate))

                            underlying = ql.MakeVanillaSwap(
                                l, idx, s, ql.Period(0, ql.Days),
                                effectiveDate=start_date,
                                fixedLegTenor=ql.Period(1, ql.Years),
                                fixedLegDayCount=ql.Thirty360(),
                                floatingLegSpread=0.0,
                                swapType=swap_type[h])
                            underlying.setPricingEngine(swap_engine)

                            fair_rate = underlying.fairRate()

                            swaption = ql.Swaption(underlying,
                                                   ql.EuropeanExercise(
                                                       exercise_date),
                                                   types[h],
                                                   methods[h])
                            swaption.setPricingEngine(swaption_engine)

                            value = swaption.NPV()
                            delta = swaption.delta() * bump

                            projection_quote_handle.linkTo(
                                ql.SimpleQuote(projection_rate + bump))

                            bumped_fair_rate = underlying.fairRate()
                            bumped_value = swaption.NPV()
                            bumped_delta = swaption.delta() * bump

                            delta_bump = bumped_fair_rate - fair_rate
                            approx_delta = (bumped_value - value) / \
                                delta_bump * bump

                            lower_bound = min(delta, bumped_delta) - epsilon
                            upper_bound = max(delta, bumped_delta) + epsilon

                            # Based on the Mean Value Theorem, the below equality
                            # should hold for any function that is monotonic in the
                            # area of the bump.
                            check_is_correct = (lower_bound <= approx_delta) and (
                                approx_delta <= upper_bound)
                            self.assertTrue(check_is_correct)

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


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(SwaptionTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
