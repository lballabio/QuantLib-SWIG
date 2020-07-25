"""
 Copyright (C) 2019 Klaus Spanderen

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

import math
import unittest

import QuantLib as ql


class ExtrapolationTest(unittest.TestCase):
    def testKnownExpExtrapolation(self):
        """Testing Richardson extrapolation of e^x at x->1 with known order of convergence"""
        f = lambda x: math.exp(1+x)
        x = ql.RichardsonExtrapolation(f, 0.01, 1.0)(4.0)

        self.assertAlmostEqual(x, math.exp(1), 4,
            msg="Unable to extrapolate exp(x) at x->1")

    def testUnknownExpExtrapolation(self):
        """Testing Richardson extrapolation of e^x at x->1 with unknown order of convergence"""
        f = lambda x: math.exp(1+x)
        x = ql.RichardsonExtrapolation(f, 0.01)(4.0, 2.0)

        self.assertAlmostEqual(x, math.exp(1), 4,
            msg="Unable to extrapolate exp(x) at x->1")


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(ExtrapolationTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
