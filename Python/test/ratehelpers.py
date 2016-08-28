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

import QuantLib
import unittest


class FixedRateBondHelperTest(unittest.TestCase):
    def setUp(self):
        self.settlement_days = 3
        self.face_amount = 100.0
        self.redemption = 100.0
        self.quote_handle = QuantLib.QuoteHandle(QuantLib.SimpleQuote(100.0))

        self.issue_date = QuantLib.Date(2, 1, 2008)
        self.maturity_date = QuantLib.Date(2, 1, 2018)
        self.calendar = QuantLib.UnitedStates(
            QuantLib.UnitedStates.GovernmentBond)
        self.day_counter = QuantLib.ActualActual(QuantLib.ActualActual.Bond)
        self.sched = QuantLib.Schedule(self.issue_date, self.maturity_date,
                                       QuantLib.Period(QuantLib.Semiannual),
                                       self.calendar,
                                       QuantLib.Unadjusted, QuantLib.Unadjusted,
                                       QuantLib.DateGeneration.Backward, False)
        self.coupons = [0.05]

        self.bond_helper = QuantLib.FixedRateBondHelper(self.quote_handle,
                                                        self.settlement_days,
                                                        self.face_amount,
                                                        self.sched,
                                                        self.coupons,
                                                        self.day_counter,
                                                        QuantLib.Following,
                                                        self.redemption,
                                                        self.issue_date)

    def testBond(self):
        """ Testing FixedRateBondHelper bond() method. """
        bond = self.bond_helper.bond()
        self.assertEqual(bond.issueDate(), self.issue_date)
        self.assertEqual(bond.nextCouponRate(), self.coupons[0])


class FxSwapRateHelperTest(unittest.TestCase):
    def build_eur_curve(self, quotes_date):
        """
        Builds the EUR OIS curve as the collateral currency discount curve
        :param quotes_date: date fro which it is assumed all market data are
            valid
        :return: tuple consisting of objects related to EUR OIS discounting
            curve: QuantLib.PiecewiseFlatForward,
                   QuantLib.YieldTermStructureHandle
                   QuantLib.RelinkableYieldTermStructureHandle
        """
        calendar = QuantLib.TARGET()
        settlementDays = 2

        todaysDate = quotes_date
        QuantLib.Settings.instance().evaluationDate = todaysDate

        todays_Eonia_quote = -0.00341

        # market quotes
        # deposits, key structure as (settlement_days_number, number_of_units_
        # for_maturity, unit)
        deposits = {(0, 1, QuantLib.Days): todays_Eonia_quote}

        discounting_yts_handle = QuantLib.RelinkableYieldTermStructureHandle()
        on_index = QuantLib.Eonia(discounting_yts_handle)
        on_index.addFixing(todaysDate, todays_Eonia_quote / 100.0)

        ois = {(1, QuantLib.Weeks): -0.342,
               (1, QuantLib.Months): -0.344,
               (3, QuantLib.Months): -0.349,
               (6, QuantLib.Months): -0.363,
               (1, QuantLib.Years): -0.389}

        # convert them to Quote objects
        for sett_num, n, unit in deposits.keys():
            deposits[(sett_num, n, unit)] = \
                QuantLib.SimpleQuote(deposits[(sett_num, n, unit)] / 100.0)

        for n, unit in ois.keys():
            ois[(n, unit)] = QuantLib.SimpleQuote(ois[(n, unit)] / 100.0)

        # build rate helpers
        dayCounter = QuantLib.Actual360()
        # looping left if somone wants two add more deposits to tests, e.g. T/N

        depositHelpers = [QuantLib.DepositRateHelper(
            QuantLib.QuoteHandle(deposits[(sett_num, n,
                                           unit)]),
            QuantLib.Period(n, unit), sett_num,
            calendar, QuantLib.ModifiedFollowing,
            True, dayCounter)
                          for sett_num, n, unit in deposits.keys()]

        oisHelpers = [QuantLib.OISRateHelper(settlementDays,
                                             QuantLib.Period(n, unit),
                                             QuantLib.QuoteHandle(
                                                 ois[(n, unit)]),
                                             on_index,
                                             discounting_yts_handle)
                      for n, unit in ois.keys()]

        rateHelpers = depositHelpers + oisHelpers

        # term-structure construction
        oisSwapCurve = QuantLib.PiecewiseFlatForward(todaysDate, rateHelpers,
                                                     QuantLib.Actual360())
        oisSwapCurve.enableExtrapolation()
        return oisSwapCurve, QuantLib.YieldTermStructureHandle(oisSwapCurve), \
               QuantLib.RelinkableYieldTermStructureHandle(oisSwapCurve)

    def build_pln_fx_swap_curve(self, base_ccy_yts, fx_swaps, fx_spot):

        todaysDate = base_ccy_yts.referenceDate()
        # I am not sure if that is required, but I guss it is worth setting
        # up just in case somewhere another thread updates this setting.
        QuantLib.Settings.instance().evaluationDate = todaysDate

        calendar = QuantLib.JointCalendar(QuantLib.TARGET(), QuantLib.Poland())
        spot_date_lag = 2
        spot_date = calendar.advance(todaysDate, spot_date_lag, QuantLib.Days)

        # market quotes
        # Update deposit Rates ( usual source will be Euribor Fixings on the Curve Date


        # build rate helpers

        spotFx = QuantLib.SimpleQuote(fx_spot)

        fxSwapHelpers = [QuantLib.FxSwapRateHelper(
            QuantLib.QuoteHandle(
                QuantLib.SimpleQuote(fx_swaps[(n, unit)])),
            QuantLib.QuoteHandle(spotFx),
            QuantLib.Period(n, unit),
            spot_date_lag,
            calendar,
            QuantLib.ModifiedFollowing,
            True, True,
            base_ccy_yts)
                         for n, unit in fx_swaps.keys()]

        # term structure handles

        discountTermStructure = QuantLib.RelinkableYieldTermStructureHandle()

        # term-structure construction

        fxSwapCurve = QuantLib.PiecewiseFlatForward(todaysDate, fxSwapHelpers,
                                                    QuantLib.Actual365Fixed())
        fxSwapCurve.enableExtrapolation()
        return fxSwapCurve, QuantLib.YieldTermStructureHandle(fxSwapCurve), \
               QuantLib.RelinkableYieldTermStructureHandle(fxSwapCurve), \
               fxSwapHelpers

    def setUp(self):
        self.today = QuantLib.Date(26, 8, 2016)
        self.eur_ois_curve, self.eur_ois_handle, self.eur_ois_rel_handle = \
            self.build_eur_curve(self.today)

        self.fx_swap_quotes = {(1, QuantLib.Months): 20e-4,
                               (3, QuantLib.Months): 60e-4,
                               (6, QuantLib.Months): 120e-4,
                               (1, QuantLib.Years): 240e-4}
        self.maturities = [QuantLib.Date(30, 9, 2016),
                           QuantLib.Date(30, 11, 2016),
                           QuantLib.Date(28, 2, 2017),
                           QuantLib.Date(30, 8, 2017)]

        self.fx_spot_quote = 4.3

        self.pln_eur_implied_curve, self.pln_eur_implied_curve_handle, \
        self.pln_eur_implied_curve_relinkable_handle, \
        self.eur_pln_fx_swap_helpers = self.build_pln_fx_swap_curve(
                                                        self.eur_ois_rel_handle,
                                                        self.fx_swap_quotes,
                                                        self.fx_spot_quote)

    def testQuote(self):

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

        # Check if still the test date is unchanged, otherwise all other
        # tests here make no sense.
        self.assertEquals(self.today, QuantLib.Date(26, 8, 2016))

        # Hard coded expected maturities of fx swaps


        for n in range(len(self.maturities)):
            self.assertEquals(self.maturities[n],
                              self.eur_pln_fx_swap_helpers[n].latestDate())

    def testImpliedRates(self):

        # Not sure if all Python versions and machine will guarantee that the
        #  lists are not messed, probably some ordered maps should be used
        # here while retrieving values from fx_swap_quotes dictionary
        original_quotes = list(self.fx_swap_quotes.values())
        spot_date = QuantLib.Date(30, 8, 2016)
        spot_df = self.eur_ois_curve.discount(spot_date) / \
                  self.pln_eur_implied_curve.discount(spot_date)

        for n in range(len(original_quotes)):
            original_quote = original_quotes[n]
            maturity = self.maturities[n]
            original_forward = self.fx_spot_quote + original_quote
            curve_impl_forward = self.fx_spot_quote * \
                                 self.eur_ois_curve.discount(maturity) / \
                                 self.pln_eur_implied_curve.discount(
                                     maturity) / spot_df

            self.assertAlmostEqual(original_forward, curve_impl_forward,
                                   places=6)


if __name__ == '__main__':
    print('testing QuantLib ' + QuantLib.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(FixedRateBondHelperTest, 'test'))
    unittest.TextTestRunner(verbosity=2).run(suite)
