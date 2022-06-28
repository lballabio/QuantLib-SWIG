"""
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2009 Joseph Malicki

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

import sys
import unittest

from date import DateTest
from daycounters import DayCountersTest
from instruments import InstrumentTest
from marketelements import MarketElementTest
from integrals import IntegralTest
from solvers1d import Solver1DTest
from termstructures import TermStructureTest
from bonds import FixedRateBondTest, FixedRateBondKwargsTest
from ratehelpers import (
    FixedRateBondHelperTest,
    FxSwapRateHelperTest,
    OISRateHelperTest,
    CrossCurrencyBasisSwapRateHelperTest)
from cms import CmsTest
from assetswap import AssetSwapTest
from capfloor import CapFloorTest
from blackformula import BlackFormulaTest
from blackformula import BlackDeltaCalculatorTest
from iborindex import IborIndexTest
from sabr import SabrTest
from slv import SlvTest
from ode import OdeTest
from americanquantooption import AmericanQuantoOptionTest
from extrapolation import ExtrapolationTest
from fdm import FdmTest
from swaption import SwaptionTest
from volatilities import SviSmileSectionTest, SwaptionVolatilityCubeTest
from inflation import InflationTest
from coupons import (
    CashFlowsTest,
    SubPeriodsCouponTest,
    IborCouponTest,
    OvernightCouponTest,
    FixedRateCouponTest)
from options import OptionsTest
from swap import ZeroCouponSwapTest
from currencies import CurrencyTest


def test():
    import QuantLib
    print('testing QuantLib ' + QuantLib.__version__)

    suite = unittest.TestSuite()

    suite.addTest(unittest.makeSuite(DateTest, 'test'))
    suite.addTest(DayCountersTest())
    suite.addTest(unittest.makeSuite(InstrumentTest, 'test'))
    suite.addTest(unittest.makeSuite(MarketElementTest, 'test'))
    suite.addTest(unittest.makeSuite(IntegralTest, 'test'))
    suite.addTest(Solver1DTest())
    suite.addTest(unittest.makeSuite(TermStructureTest, 'test'))
    suite.addTest(unittest.makeSuite(FixedRateBondTest, 'test'))
    suite.addTest(unittest.makeSuite(FixedRateBondKwargsTest, 'test'))
    suite.addTest(unittest.makeSuite(FixedRateBondHelperTest, 'test'))
    suite.addTest(unittest.makeSuite(CmsTest, 'test'))
    suite.addTest(unittest.makeSuite(AssetSwapTest, 'test'))
    suite.addTest(unittest.makeSuite(OISRateHelperTest, "test"))
    suite.addTest(unittest.makeSuite(FxSwapRateHelperTest, 'test'))
    suite.addTest(unittest.makeSuite(CapFloorTest, 'test'))
    suite.addTest(unittest.makeSuite(BlackFormulaTest, 'test'))
    suite.addTest(unittest.makeSuite(BlackDeltaCalculatorTest, 'test'))
    suite.addTest(unittest.makeSuite(IborIndexTest, 'test'))
    suite.addTest(unittest.makeSuite(SabrTest, 'test'))
    suite.addTest(unittest.makeSuite(SlvTest, 'test'))
    suite.addTest(unittest.makeSuite(OdeTest, 'test'))
    suite.addTest(unittest.makeSuite(AmericanQuantoOptionTest, 'test'))
    suite.addTest(unittest.makeSuite(ExtrapolationTest, 'test'))
    suite.addTest(unittest.makeSuite(FdmTest, 'test'))
    suite.addTest(unittest.makeSuite(SwaptionTest, "test"))
    suite.addTest(unittest.makeSuite(SwaptionVolatilityCubeTest, 'test'))
    suite.addTest(unittest.makeSuite(InflationTest, "test"))
    suite.addTest(unittest.makeSuite(CrossCurrencyBasisSwapRateHelperTest, "test"))
    suite.addTest(unittest.makeSuite(CashFlowsTest, "test"))
    suite.addTest(unittest.makeSuite(SubPeriodsCouponTest, "test"))
    suite.addTest(unittest.makeSuite(IborCouponTest, "test"))
    suite.addTest(unittest.makeSuite(OvernightCouponTest, "test"))
    suite.addTest(unittest.makeSuite(FixedRateCouponTest, "test"))
    suite.addTest(unittest.makeSuite(OptionsTest, "test"))
    suite.addTest(unittest.makeSuite(ZeroCouponSwapTest, "test"))
    suite.addTest(unittest.makeSuite(CurrencyTest, "test"))
    suite.addTest(unittest.makeSuite(SviSmileSectionTest, "test"))

    result = unittest.TextTestRunner(verbosity=2).run(suite)

    if not result.wasSuccessful():
        sys.exit(1)


if __name__ == '__main__':
    test()
