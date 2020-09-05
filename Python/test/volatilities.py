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


EURIBOR_IDX = ql.Euribor6M()

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

SMILE_OPT_TENORS = (ql.Period(1, ql.Years),
                    ql.Period(10, ql.Years),
                    ql.Period(30, ql.Years))

SMILE_SWAP_TENORS = (ql.Period(2, ql.Years),
                     ql.Period(10, ql.Years),
                     ql.Period(30, ql.Years))

NORM_VOL_SPREADS = (
    (-0.0006, 0.0005, 0.0, 0.0006, 0.0006),
    (-0.0006, 0.0005, 0.0, 0.00065, 0.0006),
    (-0.0006, 0.0001, 0.0, 0.0006, 0.0006),
    (-0.0006, 0.0005, 0.0, 0.0006, 0.0006),
    (-0.0006, 0.0005, 0.0, 0.0006, 0.0006),
    (-0.0003, 0.0005, 0.0, 0.0003, 0.0003),
    (-0.0006, 0.0005, 0.0, 0.0006, 0.0006),
    (-0.0006, 0.0005, 0.0, 0.0006, 0.0006),
    (-0.0003, 0.0005, 0.0, 0.0003, 0.0003))

ZERO_COUPON_DATA = (
    (ql.Period(1, ql.Years), -0.003),
    (ql.Period(2, ql.Years), 0.0001),
    (ql.Period(3, ql.Years), 0.001),
    (ql.Period(4, ql.Years), 0.0013),
    (ql.Period(5, ql.Years), 0.0015),
    (ql.Period(6, ql.Years), 0.0019),
    (ql.Period(7, ql.Years), 0.0021),
    (ql.Period(8, ql.Years), 0.0025),
    (ql.Period(9, ql.Years), 0.003),
    (ql.Period(10, ql.Years), 0.005))

NORM_VOL_MATRIX = ql.SwaptionVolatilityMatrix(
    CAL,
    ql.ModifiedFollowing,
    ATM_NORM_VOL_OPT_TENORS,
    ATM_NORM_VOL_SWAP_TENORS,
    ql.Matrix(ATM_NORM_VOLS),
    ql.Actual365Fixed(),
    False,
    ql.Normal)


def create_euribor_swap_idx(
        projectionCurveHandle):
    return ql.EuriborSwapIsdaFixA(ql.Period(1, ql.Years),
                                  projectionCurveHandle)


def build_nominal_term_structure(valuation_date, nominal_quotes):
    dates, rates = zip(*[(CAL.advance(valuation_date, x[0]), x[1])
                         for x in nominal_quotes])
    return ql.ZeroCurve(dates, rates, ql.Actual365Fixed())


def create_linear_swaption_cube(
        volatility_matrix,
        spread_opt_tenors,
        spread_swap_tenors,
        vol_spreads):
    return
