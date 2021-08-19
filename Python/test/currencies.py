"""
 Copyright (C) 2021 Marcin Rybacki

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


class CurrencyTest(unittest.TestCase):

    def test_default_currency_constructor(self):
        """Testing default currency constructor"""
        fail_msg = "Failed to create default currency."
        default_ccy = ql.Currency()
        self.assertTrue(default_ccy.empty(), fail_msg)

    def test_eur_constructor(self):
        """Testing EUR constructor"""
        fail_msg = "Failed to create EUR currency."
        eur = ql.EURCurrency()
        self.assertFalse(eur.empty(), fail_msg)

    def test_bespoke_currency_constructor(self):
        """Testing bespoke currency constructor"""
        fail_msg = "Failed to create bespoke currency."
        custom_ccy = ql.Currency(
            "CCY", "CCY", 100, "#", "", 100, ql.Rounding(), "")
        self.assertFalse(custom_ccy.empty(), fail_msg)


if __name__ == '__main__':
    print('testing QuantLib ' + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(CurrencyTest, 'test'))
    unittest.TextTestRunner(verbosity=2).run(suite)
