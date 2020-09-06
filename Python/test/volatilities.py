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


EPSILON = 1.e-10

CAL = ql.TARGET()

ATM_NORM_VOLS = (
    (0.0086, 0.00128, 0.00195, 0.00269, 0.00327, 0.00361, 0.00387,
     0.00409, 0.00427, 0.00443, 0.00488, 0.00504, 0.00508, 0.00504),
    (0.0092, 0.00134, 0.00197, 0.00264, 0.00319, 0.00352, 0.00383,
     0.00402, 0.00419, 0.00431, 0.00478, 0.00499, 0.00507, 0.00503),
    (0.00112, 0.00153, 0.00210, 0.00276, 0.00327, 0.00353, 0.00384,
     0.00408, 0.00426, 0.00445, 0.00486, 0.00505, 0.00509, 0.00510),
    (0.00129, 0.00171, 0.00226, 0.00288, 0.00335, 0.00360, 0.00388,
     0.00410, 0.00430, 0.00446, 0.00487, 0.00506, 0.00511, 0.00510),
    (0.00146, 0.00187, 0.00246, 0.00301, 0.00342, 0.00369, 0.00393,
     0.00413, 0.00432, 0.00449, 0.00489, 0.00510, 0.00513, 0.00515),
    (0.00165, 0.00209, 0.00263, 0.00313, 0.00350, 0.00376, 0.00400,
     0.00420, 0.00437, 0.00453, 0.00488, 0.00509, 0.00514, 0.00517),
    (0.00209, 0.00253, 0.00300, 0.00340, 0.00370, 0.00395, 0.00419,
     0.00434, 0.00450, 0.00464, 0.00493, 0.00510, 0.00513, 0.00519),
    (0.00251, 0.00289, 0.00332, 0.00362, 0.00392, 0.00412, 0.00432,
     0.00447, 0.00460, 0.00473, 0.00496, 0.00510, 0.00513, 0.00516),
    (0.00340, 0.00366, 0.00392, 0.00411, 0.00432, 0.00445, 0.00461,
     0.00472, 0.00480, 0.00490, 0.00503, 0.00513, 0.00513, 0.00512),
    (0.00403, 0.00418, 0.00436, 0.00449, 0.00461, 0.00471, 0.00482,
     0.00492, 0.00499, 0.00505, 0.00512, 0.00513, 0.00509, 0.00507),
    (0.00440, 0.00448, 0.00460, 0.00471, 0.00484, 0.00491, 0.00499,
     0.00507, 0.00514, 0.00519, 0.00516, 0.00514, 0.00506, 0.00502),
    (0.00496, 0.00497, 0.00504, 0.00512, 0.00518, 0.00522, 0.00526,
     0.00529, 0.00533, 0.00538, 0.00526, 0.00517, 0.00504, 0.00496),
    (0.00539, 0.00537, 0.00540, 0.00542, 0.00544, 0.00545, 0.00545,
     0.00544, 0.00544, 0.00549, 0.00531, 0.00518, 0.00501, 0.00491),
    (0.00540, 0.00537, 0.00538, 0.00537, 0.00535, 0.00536, 0.00535,
     0.00533, 0.00535, 0.00537, 0.00514, 0.00498, 0.00479, 0.00466),
    (0.00528, 0.00524, 0.00526, 0.00523, 0.00522, 0.00523, 0.00520,
     0.00519, 0.00518, 0.00518, 0.00495, 0.00474, 0.00454, 0.00438),
    (0.00514, 0.00512, 0.00513, 0.00510, 0.00508, 0.00507, 0.00503,
     0.00499, 0.00498, 0.00497, 0.00476, 0.00453, 0.00431, 0.00414),
    (0.00496, 0.00496, 0.00497, 0.00495, 0.00495, 0.00492, 0.00486,
     0.00479, 0.00474, 0.00471, 0.00451, 0.00429, 0.00408, 0.00392))

ATM_NORM_VOL_OPT_TENORS = (ql.Period(1, ql.Months), ql.Period(2, ql.Months),
                           ql.Period(3, ql.Months), ql.Period(6, ql.Months),
                           ql.Period(9, ql.Months), ql.Period(1, ql.Years),
                           ql.Period(18, ql.Months), ql.Period(2, ql.Years),
                           ql.Period(3, ql.Years), ql.Period(4, ql.Years),
                           ql.Period(5, ql.Years), ql.Period(7, ql.Years),
                           ql.Period(10, ql.Years), ql.Period(15, ql.Years),
                           ql.Period(20, ql.Years), ql.Period(25, ql.Years),
                           ql.Period(30, ql.Years))

ATM_NORM_VOL_SWAP_TENORS = (ql.Period(1, ql.Years), ql.Period(2, ql.Years),
                            ql.Period(3, ql.Years), ql.Period(4, ql.Years),
                            ql.Period(5, ql.Years), ql.Period(6, ql.Years),
                            ql.Period(7, ql.Years), ql.Period(8, ql.Years),
                            ql.Period(9, ql.Years), ql.Period(10, ql.Years),
                            ql.Period(15, ql.Years), ql.Period(20, ql.Years),
                            ql.Period(25, ql.Years), ql.Period(30, ql.Years))

ATM_LOGNORM_VOLS = (
    (0.1300, 0.1560, 0.1390, 0.1220),
    (0.1440, 0.1580, 0.1460, 0.1260),
    (0.1600, 0.1590, 0.1470, 0.1290),
    (0.1640, 0.1470, 0.1370, 0.1220),
    (0.1400, 0.1300, 0.1250, 0.1100),
    (0.1130, 0.1090, 0.1070, 0.0930))

ATM_LOGNORM_VOL_OPT_TENORS = (ql.Period(1, ql.Months),
                              ql.Period(6, ql.Months),
                              ql.Period(3, ql.Years),
                              ql.Period(5, ql.Years),
                              ql.Period(10, ql.Years),
                              ql.Period(25, ql.Years))

ATM_LOGNORM_VOL_SWAP_TENORS = (ql.Period(1, ql.Years),
                               ql.Period(3, ql.Years),
                               ql.Period(10, ql.Years),
                               ql.Period(25, ql.Years))

SMILE_OPT_TENORS = (ql.Period(1, ql.Years),
                    ql.Period(10, ql.Years),
                    ql.Period(30, ql.Years))

SMILE_SWAP_TENORS = (ql.Period(2, ql.Years),
                     ql.Period(10, ql.Years),
                     ql.Period(30, ql.Years))

STRIKE_SPREADS = (-0.01, -0.005, 0.0, 0.005, 0.01)

NORM_VOL_SPREADS = (
    (-0.0006, 0.0005, 0.0, 0.0006, 0.0006),
    (-0.0006, 0.0005, 0.0, 0.00065, 0.00066),
    (-0.0006, 0.0001, 0.0, 0.0006, 0.0006),
    (-0.0006, 0.0005, 0.0, 0.0006, 0.0006),
    (-0.0006, 0.0005, 0.0, 0.0006, 0.0006),
    (-0.0003, 0.0005, 0.0, 0.0003, 0.0006),
    (-0.0006, 0.0005, 0.0, 0.0006, 0.0006),
    (-0.0006, 0.0005, 0.0, 0.0006, 0.0006),
    (-0.0003, 0.0005, 0.0, 0.0003, 0.0006))

LOGNORM_VOL_SPREADS = (
    (0.00599, 0.0049, 0.0000, -0.0001, -0.0049),
    (0.00729, 0.0086, 0.0000, -0.0024, -0.0049),
    (0.00738, 0.0102, 0.0000, -0.0039, -0.0049),
    (0.00465, 0.0063, 0.0000, -0.0032, -0.0049),
    (0.00558, 0.0084, 0.0000, -0.0050, -0.0049),
    (0.00576, 0.0083, 0.0000, -0.0043, -0.0049),
    (0.0049, 0.0059, 0.0000, -0.0030, -0.0049),
    (0.0049, 0.0078, 0.0000, -0.0045, -0.0049),
    (0.0049, 0.0079, 0.0000, -0.0042, -0.0049))

ZERO_COUPON_DATA = (
    (ql.Period(1, ql.Years), 0.003),
    (ql.Period(2, ql.Years), 0.005),
    (ql.Period(3, ql.Years), 0.006),
    (ql.Period(4, ql.Years), 0.007),
    (ql.Period(5, ql.Years), 0.009),
    (ql.Period(10, ql.Years), 0.011),
    (ql.Period(15, ql.Years), 0.014),
    (ql.Period(20, ql.Years), 0.016),
    (ql.Period(30, ql.Years), 0.019))

NORM_VOL_MATRIX = ql.SwaptionVolatilityMatrix(
    CAL,
    ql.ModifiedFollowing,
    ATM_NORM_VOL_OPT_TENORS,
    ATM_NORM_VOL_SWAP_TENORS,
    ql.Matrix(ATM_NORM_VOLS),
    ql.Actual365Fixed(),
    False,
    ql.Normal)

LOGNORM_VOL_MATRIX = ql.SwaptionVolatilityMatrix(
    CAL,
    ql.ModifiedFollowing,
    ATM_LOGNORM_VOL_OPT_TENORS,
    ATM_LOGNORM_VOL_SWAP_TENORS,
    ql.Matrix(ATM_LOGNORM_VOLS),
    ql.Actual365Fixed(),
    False,
    ql.ShiftedLognormal)


def build_euribor_swap_idx(
        projectionCurveHandle):
    return ql.EuriborSwapIsdaFixA(ql.Period(1, ql.Years),
                                  projectionCurveHandle)


def build_nominal_term_structure(valuation_date, nominal_quotes):
    dates, rates = zip(*[(CAL.advance(valuation_date, x[0]), x[1])
                         for x in nominal_quotes])
    crv = ql.ZeroCurve(dates, rates, ql.Actual365Fixed())
    crv.enableExtrapolation()
    return crv


def build_linear_swaption_cube(
        volatility_matrix,
        spread_opt_tenors,
        spread_swap_tenors,
        strike_spreads,
        vol_spreads,
        swap_index_base,
        short_swap_index_base=None,
        vega_weighted_smile_fit=False):
    vol_spreads = [[ql.QuoteHandle(ql.SimpleQuote(v)) for v in row]
                   for row in vol_spreads]
    cube = ql.SwaptionVolCube2(
        ql.SwaptionVolatilityStructureHandle(volatility_matrix),
        spread_opt_tenors,
        spread_swap_tenors,
        strike_spreads,
        vol_spreads,
        swap_index_base,
        short_swap_index_base if short_swap_index_base else swap_index_base,
        vega_weighted_smile_fit)
    cube.enableExtrapolation()
    return cube


class SwaptionVolatilityCubeTest(unittest.TestCase):
    def setUp(self):
        self.today = CAL.adjust(ql.Date.todaysDate())
        ql.Settings.instance().setEvaluationDate(self.today)

        curve_handle = ql.RelinkableYieldTermStructureHandle()
        curve = build_nominal_term_structure(self.today, ZERO_COUPON_DATA)
        curve_handle.linkTo(curve)

        self.idx = ql.Euribor6M(curve_handle)
        self.swap_idx = build_euribor_swap_idx(curve_handle)
        self.swap_engine = ql.DiscountingSwapEngine(curve_handle)

    def _get_fair_rate(self, option_tenor, swap_tenor):
        exercise_date = CAL.advance(self.today, option_tenor)
        start_date = CAL.advance(exercise_date, ql.Period(2, ql.Days))
        underlying = ql.MakeVanillaSwap(
            swap_tenor, self.idx, 0.0, ql.Period(0, ql.Days),
            effectiveDate=start_date,
            fixedLegTenor=ql.Period(1, ql.Years),
            fixedLegDayCount=ql.Thirty360(),
            floatingLegSpread=0.0,
            swapType=ql.VanillaSwap.Receiver)
        underlying.setPricingEngine(self.swap_engine)
        return underlying.fairRate()

    def _assert_atm_strike(
            self, cube, interpolation, vol_type):
        opt_tenor = ql.Period(1, ql.Years)
        swap_tenor = ql.Period(10, ql.Years)
        expected_atm_strike = self._get_fair_rate(opt_tenor, swap_tenor)
        actual_atm_strike = cube.atmStrike(
            cube.optionDateFromTenor(opt_tenor), swap_tenor)
        fail_msg = """ 
                    ATM strike test failed for:
                        cube interpolation: {interpolation}
                        volatility_type: {vol_type}
                        option tenor: {option_tenor}
                        swap tenor: {swap_tenor}
                        strike: {strike}
                        replicated strike: {replicated_strike}
                   """.format(interpolation=interpolation,
                              vol_type=vol_type,
                              option_tenor=opt_tenor,
                              swap_tenor=swap_tenor,
                              strike=actual_atm_strike,
                              replicated_strike=expected_atm_strike)
        self.assertAlmostEquals(
            first=actual_atm_strike,
            second=expected_atm_strike,
            delta=EPSILON,
            msg=fail_msg)

    def _assert_atm_vol(
            self,
            cube,
            opt_tenor,
            swap_tenor,
            expected_vol,
            interpolation,
            vol_type):
        option_date = cube.optionDateFromTenor(opt_tenor)
        strike = cube.atmStrike(option_date, swap_tenor)
        actual_vol = cube.volatility(option_date, swap_tenor, strike)
        fail_msg = """ 
                    ATM vol test failed for:
                        cube interpolation: {interpolation}
                        volatility_type: {vol_type}
                        option tenor: {option_tenor}
                        swap tenor: {swap_tenor}
                        strike: {strike}
                        volatility: {vol}
                        replicated volatility: {replicated_vol}
                   """.format(interpolation=interpolation,
                              vol_type=vol_type,
                              option_tenor=opt_tenor,
                              swap_tenor=swap_tenor,
                              strike=strike,
                              vol=actual_vol,
                              replicated_vol=expected_vol)
        self.assertAlmostEquals(
            first=actual_vol,
            second=expected_vol,
            delta=EPSILON,
            msg=fail_msg)

    def _assert_vol_spread(
            self,
            cube,
            opt_tenor,
            swap_tenor,
            strike_spread,
            expected_vol,
            interpolation,
            vol_type):
        option_date = cube.optionDateFromTenor(opt_tenor)
        strike = cube.atmStrike(option_date, swap_tenor) + strike_spread
        actual_vol = cube.volatility(option_date, swap_tenor, strike)
        fail_msg = """ 
                    Vol spread test failed for:
                        cube interpolation: {interpolation}
                        volatility_type: {vol_type}
                        option tenor: {option_tenor}
                        swap tenor: {swap_tenor}
                        strike: {strike}
                        volatility: {vol}
                        replicated volatility: {replicated_vol}
                   """.format(interpolation=interpolation,
                              vol_type=vol_type,
                              option_tenor=opt_tenor,
                              swap_tenor=swap_tenor,
                              strike=strike,
                              vol=actual_vol,
                              replicated_vol=expected_vol)
        self.assertAlmostEquals(
            first=actual_vol,
            second=expected_vol,
            delta=EPSILON,
            msg=fail_msg)

    def test_linear_normal_cube_at_the_money_strike(self):
        """Testing ATM strike for linearly interpolated normal vol cube"""
        linear_cube = build_linear_swaption_cube(
            NORM_VOL_MATRIX,
            SMILE_OPT_TENORS,
            SMILE_SWAP_TENORS,
            STRIKE_SPREADS,
            NORM_VOL_SPREADS,
            self.swap_idx)
        self._assert_atm_strike(
            cube=linear_cube,
            interpolation='linear',
            vol_type='normal')

    def test_linear_lognormal_cube_at_the_money_strike(self):
        """Testing ATM strike for linearly interpolated log-normal vol cube"""
        linear_cube = build_linear_swaption_cube(
            LOGNORM_VOL_MATRIX,
            SMILE_OPT_TENORS,
            SMILE_SWAP_TENORS,
            STRIKE_SPREADS,
            LOGNORM_VOL_SPREADS,
            self.swap_idx)
        self._assert_atm_strike(
            cube=linear_cube,
            interpolation='linear',
            vol_type='log-normal')

    def test_linear_normal_cube_at_the_money_vol(self):
        """Testing ATM volatility for linearly interpolated normal vol cube"""
        linear_cube = build_linear_swaption_cube(
            NORM_VOL_MATRIX,
            SMILE_OPT_TENORS,
            SMILE_SWAP_TENORS,
            STRIKE_SPREADS,
            NORM_VOL_SPREADS,
            self.swap_idx)
        self._assert_atm_vol(
            cube=linear_cube,
            opt_tenor=ql.Period(1, ql.Years),
            swap_tenor=ql.Period(10, ql.Years),
            expected_vol=0.00453,
            interpolation='linear',
            vol_type='normal')

    def test_linear_lognormal_cube_at_the_money_vol(self):
        """Testing ATM volatility for linearly interpolated log-normal vol cube"""
        linear_cube = build_linear_swaption_cube(
            LOGNORM_VOL_MATRIX,
            SMILE_OPT_TENORS,
            SMILE_SWAP_TENORS,
            STRIKE_SPREADS,
            LOGNORM_VOL_SPREADS,
            self.swap_idx)
        self._assert_atm_vol(
            cube=linear_cube,
            opt_tenor=ql.Period(10, ql.Years),
            swap_tenor=ql.Period(10, ql.Years),
            expected_vol=0.1250,
            interpolation='linear',
            vol_type='log-normal')

    def test_linear_normal_cube_spread_vol(self):
        """Testing spread volatility for linearly interpolated normal cube"""
        linear_cube = build_linear_swaption_cube(
            NORM_VOL_MATRIX,
            SMILE_OPT_TENORS,
            SMILE_SWAP_TENORS,
            STRIKE_SPREADS,
            NORM_VOL_SPREADS,
            self.swap_idx)
        self._assert_vol_spread(
            cube=linear_cube,
            opt_tenor=ql.Period(1, ql.Years),
            swap_tenor=ql.Period(10, ql.Years),
            strike_spread=0.005,
            expected_vol=0.00453 + 0.00065,
            interpolation='linear',
            vol_type='normal')

    def test_linear_lognormal_cube_spread_vol(self):
        """Testing spread volatility for linearly interpolated log-normal cube"""
        linear_cube = build_linear_swaption_cube(
            LOGNORM_VOL_MATRIX,
            SMILE_OPT_TENORS,
            SMILE_SWAP_TENORS,
            STRIKE_SPREADS,
            LOGNORM_VOL_SPREADS,
            self.swap_idx)
        self._assert_vol_spread(
            cube=linear_cube,
            opt_tenor=ql.Period(10, ql.Years),
            swap_tenor=ql.Period(10, ql.Years),
            strike_spread=0.005,
            expected_vol=0.125 - 0.005,
            interpolation='linear',
            vol_type='log-normal')


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(SwaptionVolatilityCubeTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
