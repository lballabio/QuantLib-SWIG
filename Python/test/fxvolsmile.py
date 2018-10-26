# -*- coding: utf-8 -*-
"""
 Copyright (C) 2018 Wojciech Ślusarski

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


class FxVolSmileTest(unittest.TestCase):

    def setUp(self):
        # Castagna, Antonio and Mercurio, Fabio, Consistent Pricing of FX Options.
        # Available at SSRN: https://ssrn.com/abstract=873788
        # or http://dx.doi.org/10.2139/ssrn.873788
        self.today = ql.Date(1, 7, 2005)
        ql.Settings.instance().evaluationDate = self.today
        self.t = 94 / 365.0
        self.spot = ql.QuoteHandle(ql.SimpleQuote(1.205))
        self.expiry = [ql.Date(3, 10, 2005),
                       ql.Date(3, 7, 2006)]
        """
        Quotes of the 25-delta Risk Reversal and Vega Weighted ButterFly are 
        missing in the paper. Values typed here were selected to approximately 
        match results in Table 2.
        This test was designed to automate check if the C++ implementation 
        exposed to Python works.
        For comprehensive tests check the unittest fxvolsmile.cpp in QuantLib, 
        which was developed by Quaternion Risk Management Ltd.
        In order to understand how FxBlackVannaVolgaVolatilitySurface works, 
        have a look into examples/fxoptionvolsurfaceplot.py
        """
        self.sigma_atm = [0.0905, 0.0940]
        self.sigma_rr = [-0.005, -0.0022]
        self.sigma_vwb = [0.0013, 0.0014]

        self.df_usd = [0.9902752, 0.9585801]
        self.df_eur = [0.9945049, 0.9785056]
        self.day_count = ql.Actual365Fixed()
        self.crv_usd = ql.DiscountCurve([self.today] + self.expiry,
                                        [1.0] + self.df_usd,
                                        self.day_count)
        self.crv_eur = ql.DiscountCurve([self.today] + self.expiry,
                                        [1.0] + self.df_eur, self.day_count)
        self.domcrv_handle = ql.RelinkableYieldTermStructureHandle()
        self.forcrv_handle = ql.RelinkableYieldTermStructureHandle()
        """
        The paper does not mention calendar, TARGET does not disturb in 
        calculations required by this test. Proper calendar handling for 
        fx options is a separate issue, that should be carefully implemented 
        (e.g. see implementation for FxSwapRateHelper class).
        """
        self.calendar = ql.TARGET()


    def test_recover_standard_delta_strikes_volatility(self):
        """Testing fx volatility at volatility surface nodes"""
        self.domcrv_handle.linkTo(self.crv_usd)
        self.forcrv_handle.linkTo(self.crv_eur)
        vol_surface = ql.FxBlackVannaVolgaVolatilitySurface(self.today,
                                                            self.expiry,
                                                            self.sigma_atm,
                                                            self.sigma_rr,
                                                            self.sigma_vwb,
                                                            self.day_count,
                                                            self.calendar,
                                                            self.spot,
                                                            self.domcrv_handle,
                                                            self.forcrv_handle)
        # Table 2, page 10 in the article
        strikes_3M = [1.1733, 1.2114, 1.2487]
        strikes_1Y = [1.1597, 1.2355, 1.3148]
        expected_vols_3M = [0.0943, 0.0905,  0.0893]
        expected_vols_1Y = [0.0965, 0.0940, 0.0943]
        expiry_3M =  self.expiry[0]
        expiry_1Y = self.expiry[1]
        diff_vols_3M = [vol_surface.blackVol(expiry_3M, strike) - expected_vol
                        for strike, expected_vol in zip(strikes_3M,
                                                        expected_vols_3M)]
        diff_vols_1Y = [vol_surface.blackVol(expiry_1Y, strike) - expected_vol
                        for strike, expected_vol in zip(strikes_1Y,
                                                        expected_vols_1Y)]
        for difference, strike in zip(diff_vols_3M, strikes_3M):
            self.assertAlmostEqual(0.0, difference, delta=0.00001,
                                   msg="Failed to determine volatility at 3M "
                                       "expiry, strike: {},\n"
                                       "observed difference: "
                                       "{}".format(strike, difference))
        for difference, strike in zip(diff_vols_1Y, strikes_1Y):
            self.assertAlmostEqual(0.0, difference, delta=0.00001,
                                   msg="Failed to determine volatility at 1Y "
                                       "expiry, strike: {},\n"
                                       "observed difference: "
                                       "{}".format(strike, difference))

if __name__ == '__main__':
    print('testing QuantLib ' + ql.QuantLib.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(FxVolSmileTest, 'test'))
    unittest.TextTestRunner(verbosity=2).run(suite)