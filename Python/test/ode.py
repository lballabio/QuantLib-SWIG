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

class OdeTest(unittest.TestCase):

    def test1dODE(self):
        """ Testing one dimesnional ODE """

        yEnd = ql.RungeKutta(1e-8)(lambda x, y : y, 1, 0, 1)

        self.assertAlmostEqual(yEnd, math.exp(1), 5,
            msg="Unable to reproduce one dimensional ODE solution.")


    def test2dODE(self):
        """ Testing multi-dimesnional ODE """

        yEnd = ql.RungeKutta(1e-8)(lambda x, y : [y[1], -y[0]],
                                   [0, 1], 0, 0.5*math.pi)[0]

        self.assertAlmostEqual(yEnd, 1.0, 5,
            msg="Unable to reproduce multi-dimensional ODE solution.")


if __name__ == '__main__':
    print('testing QuantLib ' + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(OdeTest,'test'))
    unittest.TextTestRunner(verbosity=2).run(suite)
