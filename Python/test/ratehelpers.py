# coding=utf-8-unix
"""
 Copyright (C) 2009 Joseph Malicki
 Copyright (C) 2016 Wojciech Ślusarski


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


class FixedRateBondHelperTest(unittest.TestCase):
    def setUp(self):
        ql.Settings.instance().setEvaluationDate(ql.Date(2, 1, 2010))
        self.settlement_days = 3
        self.face_amount = 100.0
        self.redemption = 100.0
        self.quote_handle = ql.QuoteHandle(ql.SimpleQuote(100.0))

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

        self.bond_helper = ql.FixedRateBondHelper(
            self.quote_handle,
            self.settlement_days,
            self.face_amount,
            self.sched,
            self.coupons,
            self.day_counter,
            ql.Following,
            self.redemption,
            self.issue_date,
        )

    def testBond(self):
        """ Testing FixedRateBondHelper bond() method. """
        bond = self.bond_helper.bond()
        self.assertEqual(bond.issueDate(), self.issue_date)
        self.assertEqual(bond.nextCouponRate(), self.coupons[0])


class FxSwapRateHelperTest(unittest.TestCase):
    def setUp(self):

        # Market rates are artificial, just close to real ones.
        self.default_quote_date = ql.Date(26, 8, 2016)

        self.fx_swap_quotes = {
            (1, ql.Months): 20e-4,
            (3, ql.Months): 60e-4,
            (6, ql.Months): 120e-4,
            (1, ql.Years): 240e-4,
        }

        # Valid only for the quote date of ql.Date(26, 8, 2016)
        self.maturities = [ql.Date(30, 9, 2016), ql.Date(30, 11, 2016), ql.Date(28, 2, 2017), ql.Date(30, 8, 2017)]

        self.fx_spot_quote_EURPLN = 4.3
        self.fx_spot_quote_EURUSD = 1.1

    def build_eur_curve(self, quotes_date):
        """
        Builds the EUR OIS curve as the collateral currency discount curve
        :param quotes_date: date fro which it is assumed all market data are
            valid
        :return: tuple consisting of objects related to EUR OIS discounting
            curve: ql.PiecewiseFlatForward,
                   ql.YieldTermStructureHandle
                   ql.RelinkableYieldTermStructureHandle
        """
        calendar = ql.TARGET()
        settlementDays = 2

        todaysDate = quotes_date
        ql.Settings.instance().evaluationDate = todaysDate

        todays_Eonia_quote = -0.00341

        # market quotes
        # deposits, key structure as (settlement_days_number, number_of_units_
        # for_maturity, unit)
        deposits = {(0, 1, ql.Days): todays_Eonia_quote}

        discounting_yts_handle = ql.RelinkableYieldTermStructureHandle()
        on_index = ql.Eonia(discounting_yts_handle)
        on_index.addFixing(todaysDate, todays_Eonia_quote / 100.0)

        ois = {
            (1, ql.Weeks): -0.342,
            (1, ql.Months): -0.344,
            (3, ql.Months): -0.349,
            (6, ql.Months): -0.363,
            (1, ql.Years): -0.389,
        }

        # convert them to Quote objects
        for sett_num, n, unit in deposits.keys():
            deposits[(sett_num, n, unit)] = ql.SimpleQuote(deposits[(sett_num, n, unit)] / 100.0)

        for n, unit in ois.keys():
            ois[(n, unit)] = ql.SimpleQuote(ois[(n, unit)] / 100.0)

        # build rate helpers
        dayCounter = ql.Actual360()
        # looping left if somone wants two add more deposits to tests, e.g. T/N

        depositHelpers = [
            ql.DepositRateHelper(
                ql.QuoteHandle(deposits[(sett_num, n, unit)]),
                ql.Period(n, unit),
                sett_num,
                calendar,
                ql.ModifiedFollowing,
                True,
                dayCounter,
            )
            for sett_num, n, unit in deposits.keys()
        ]

        oisHelpers = [
            ql.OISRateHelper(
                settlementDays, ql.Period(n, unit), ql.QuoteHandle(ois[(n, unit)]), on_index, discounting_yts_handle
            )
            for n, unit in ois.keys()
        ]

        rateHelpers = depositHelpers + oisHelpers

        # term-structure construction
        oisSwapCurve = ql.PiecewiseFlatForward(todaysDate, rateHelpers, ql.Actual360())
        oisSwapCurve.enableExtrapolation()
        return (
            oisSwapCurve,
            ql.YieldTermStructureHandle(oisSwapCurve),
            ql.RelinkableYieldTermStructureHandle(oisSwapCurve),
        )

    def build_pln_fx_swap_curve(self, base_ccy_yts, fx_swaps, fx_spot):
        """
        Build curve implied from fx swap curve.
        :param base_ccy_yts:
            Relinkable yield term structure handle to curve in base currency.
        :param fx_swaps:
            Dictionary with swap points, already divided by 10,000
        :param fx_spot:
            Float value of fx spot exchange rate.
        :return: tuple consisting of objects related to fx swap implied curve:
                ql.PiecewiseFlatForward,
                ql.YieldTermStructureHandle
                ql.RelinkableYieldTermStructureHandle
                list of ql.FxSwapRateHelper
        """
        todaysDate = base_ccy_yts.referenceDate()
        # I am not sure if that is required, but I guss it is worth setting
        # up just in case somewhere another thread updates this setting.
        ql.Settings.instance().evaluationDate = todaysDate

        calendar = ql.JointCalendar(ql.TARGET(), ql.Poland())
        spot_date_lag = 2
        trading_calendar = ql.UnitedStates()

        # build rate helpers

        spotFx = ql.SimpleQuote(fx_spot)

        fxSwapHelpers = [
            ql.FxSwapRateHelper(
                ql.QuoteHandle(ql.SimpleQuote(fx_swaps[(n, unit)])),
                ql.QuoteHandle(spotFx),
                ql.Period(n, unit),
                spot_date_lag,
                calendar,
                ql.ModifiedFollowing,
                True,
                True,
                base_ccy_yts,
                trading_calendar,
            )
            for n, unit in fx_swaps.keys()
        ]

        # term-structure construction
        fxSwapCurve = ql.PiecewiseFlatForward(todaysDate, fxSwapHelpers, ql.Actual365Fixed())
        fxSwapCurve.enableExtrapolation()
        return (
            fxSwapCurve,
            ql.YieldTermStructureHandle(fxSwapCurve),
            ql.RelinkableYieldTermStructureHandle(fxSwapCurve),
            fxSwapHelpers,
        )

    def build_curves(self, quote_date):
        """
        Build all the curves in one call for a specified quote date

        :param quote_date: date for which quotes are valid,
            e.g. ql.Date(26, 8, 2016)
        """
        self.today = quote_date
        self.eur_ois_curve, self.eur_ois_handle, self.eur_ois_rel_handle = self.build_eur_curve(self.today)

        self.pln_eur_implied_curve, self.pln_eur_implied_curve_handle, self.pln_eur_implied_curve_relinkable_handle, self.eur_pln_fx_swap_helpers = self.build_pln_fx_swap_curve(
            self.eur_ois_rel_handle, self.fx_swap_quotes, self.fx_spot_quote_EURPLN
        )

    def testQuote(self):
        """ Testing FxSwapRateHelper.quote()  method. """
        self.build_curves(self.default_quote_date)
        # Not sure if all Python versions and machine will guarantee that the
        #  lists are not messed, probably some ordered maps should be used
        # here while retrieving values from fx_swap_quotes dictionary
        original_quotes = list(self.fx_swap_quotes.values())
        for n in range(len(original_quotes)):
            original_quote = original_quotes[n]
            rate_helper_quote = self.eur_pln_fx_swap_helpers[n].quote().value()
            self.assertEquals(original_quote, rate_helper_quote)

    def testLatestDate(self):
        """ Testing FxSwapRateHelper.latestDate()  method. """
        self.build_curves(self.default_quote_date)
        # Check if still the test date is unchanged, otherwise all other
        # tests here make no sense.
        self.assertEquals(self.today, ql.Date(26, 8, 2016))

        # Hard coded expected maturities of fx swaps
        for n in range(len(self.maturities)):
            self.assertEquals(self.maturities[n], self.eur_pln_fx_swap_helpers[n].latestDate())

    def testImpliedRates(self):
        """
        Testing if rates implied from the curve are returning fx forwards
        very close to those used for bootstrapping
        """
        self.build_curves(self.default_quote_date)
        # Not sure if all Python versions and machine will guarantee that the
        #  lists are not messed, probably some ordered maps should be used
        # here while retrieving values from fx_swap_quotes dictionary
        original_quotes = list(self.fx_swap_quotes.values())
        spot_date = ql.Date(30, 8, 2016)
        spot_df = self.eur_ois_curve.discount(spot_date) / self.pln_eur_implied_curve.discount(spot_date)

        for n in range(len(original_quotes)):
            original_quote = original_quotes[n]
            maturity = self.maturities[n]
            original_forward = self.fx_spot_quote_EURPLN + original_quote
            curve_impl_forward = (
                self.fx_spot_quote_EURPLN
                * self.eur_ois_curve.discount(maturity)
                / self.pln_eur_implied_curve.discount(maturity)
                / spot_df
            )

            self.assertAlmostEqual(original_forward, curve_impl_forward, places=6)

    def testFxMarketConventionsForCrossRate(self):
        """
        Testing if ql.FxSwapRateHelper obeys the fx spot market
        conventions for cross rates.
        """
        today = ql.Date(1, 7, 2016)
        spot_date = ql.Date(5, 7, 2016)
        self.build_curves(today)

        us_calendar = ql.UnitedStates()

        joint_calendar = ql.JointCalendar(ql.TARGET(), ql.Poland())

        settlement_calendar = ql.JointCalendar(joint_calendar, us_calendar)

        # Settlement should be on a day where all three centers are operating
        #  and follow EndOfMonth rule
        maturities = [
            settlement_calendar.advance(spot_date, n, unit, ql.ModifiedFollowing, True)
            for n, unit in self.fx_swap_quotes.keys()
        ]

        for n in range(len(maturities)):
            self.assertEqual(maturities[n], self.eur_pln_fx_swap_helpers[n].latestDate())

    def testFxMarketConventionsForCrossRateONPeriod(self):
        """
        Testing if ql.FxSwapRateHelper obeys the fx spot market
        conventions for cross rates' ON Period.
        """
        today = ql.Date(1, 7, 2016)
        ql.Settings.instance().evaluationDate = today

        spot_date = ql.Date(5, 7, 2016)
        fwd_points = 4.0
        # critical for ON rate helper
        on_period = ql.Period("1d")
        fixing_days = 0

        # empty RelinkableYieldTermStructureHandle is sufficient for testing
        # dates
        base_ccy_yts = ql.RelinkableYieldTermStructureHandle()

        us_calendar = ql.UnitedStates()

        joint_calendar = ql.JointCalendar(ql.TARGET(), ql.Poland())

        # Settlement should be on a day where all three centers are operating
        #  and follow EndOfMonth rule
        on_rate_helper = ql.FxSwapRateHelper(
            ql.QuoteHandle(ql.SimpleQuote(fwd_points)),
            ql.QuoteHandle(ql.SimpleQuote(self.fx_spot_quote_EURPLN)),
            on_period,
            fixing_days,
            joint_calendar,
            ql.ModifiedFollowing,
            False,
            True,
            base_ccy_yts,
            us_calendar,
        )

        self.assertEqual(spot_date, on_rate_helper.latestDate())

    def testFxMarketConventionsForCrossRateAdjustedSpotDate(self):
        """
        Testing if ql.FxSwapRateHelper obeys the fx spot market
        conventions
        """
        today = ql.Date(30, 6, 2016)
        spot_date = ql.Date(5, 7, 2016)
        self.build_curves(today)
        us_calendar = ql.UnitedStates()
        joint_calendar = ql.JointCalendar(ql.TARGET(), ql.Poland())

        settlement_calendar = ql.JointCalendar(joint_calendar, us_calendar)
        # Settlement should be on a day where all three centers are operating
        #  and follow EndOfMonth rule
        maturities = [
            joint_calendar.advance(spot_date, n, unit, ql.ModifiedFollowing, True)
            for n, unit in self.fx_swap_quotes.keys()
        ]

        maturities = [settlement_calendar.adjust(date) for date in maturities]

        for n in range(len(maturities)):
            self.assertEquals(maturities[n], self.eur_pln_fx_swap_helpers[n].latestDate())

    def testFxMarketConventionsForDatesInEURUSD_ON_Period(self):
        """
        Testing if ql.FxSwapRateHelper obeys the fx spot market
        conventions for EURUSD settlement dates on the ON Period.
        """
        today = ql.Date(1, 7, 2016)
        ql.Settings.instance().evaluationDate = today

        spot_date = ql.Date(5, 7, 2016)
        fwd_points = 4.0
        # critical for ON rate helper
        on_period = ql.Period("1d")
        fixing_days = 0

        # empty RelinkableYieldTermStructureHandle is sufficient for testing
        # dates
        base_ccy_yts = ql.RelinkableYieldTermStructureHandle()

        # In EURUSD, there must be two days to spot date in Target calendar
        # and one day in US, therefore it is sufficient to pass only Target
        # as a base calendar
        calendar = ql.TARGET()
        trading_calendar = ql.UnitedStates()

        on_rate_helper = ql.FxSwapRateHelper(
            ql.QuoteHandle(ql.SimpleQuote(fwd_points)),
            ql.QuoteHandle(ql.SimpleQuote(self.fx_spot_quote_EURUSD)),
            on_period,
            fixing_days,
            calendar,
            ql.ModifiedFollowing,
            False,
            True,
            base_ccy_yts,
            trading_calendar,
        )

        self.assertEqual(spot_date, on_rate_helper.latestDate())

    def testFxMarketConventionsForDatesInEURUSD_ShortEnd(self):
        """
        Testing if ql.FxSwapRateHelper obeys the fx spot market
        conventions for EURUSD settlement dates on the 3M tenor.
        """
        today = ql.Date(1, 7, 2016)
        ql.Settings.instance().evaluationDate = today

        expected_3M_date = ql.Date(5, 10, 2016)
        fwd_points = 4.0
        # critical for ON rate helper
        period = ql.Period("3M")
        fixing_days = 2

        # empty RelinkableYieldTermStructureHandle is sufficient for testing
        # dates
        base_ccy_yts = ql.RelinkableYieldTermStructureHandle()

        # In EURUSD, there must be two days to spot date in Target calendar
        # and one day in US, therefore it is sufficient to pass only Target
        # as a base calendar. Passing joint calendar would result in wrong
        # spot date of the trade
        calendar = ql.TARGET()
        trading_calendar = ql.UnitedStates()

        rate_helper = ql.FxSwapRateHelper(
            ql.QuoteHandle(ql.SimpleQuote(fwd_points)),
            ql.QuoteHandle(ql.SimpleQuote(self.fx_spot_quote_EURUSD)),
            period,
            fixing_days,
            calendar,
            ql.ModifiedFollowing,
            True,
            True,
            base_ccy_yts,
            trading_calendar,
        )

        self.assertEqual(expected_3M_date, rate_helper.latestDate())

    def tearDown(self):
        ql.Settings.instance().setEvaluationDate(ql.Date())


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(FixedRateBondHelperTest, "test"))
    suite.addTest(unittest.makeSuite(FxSwapRateHelperTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
