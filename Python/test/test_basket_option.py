"""
 Copyright (C) 2024 Klaus Spanderen

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


class BasketOptionTest(unittest.TestCase):
    def setUp(self):
        self.todaysDate = ql.Date(26, ql.October, 2024)
        ql.Settings.instance().evaluationDate = self.todaysDate

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def testThreeAssetSpreadOption(self):
        """Testing three asset spread option"""

        def build_process(s: float, q: float, v: float) -> ql.BlackScholesMertonProcess:
            return ql.BlackScholesMertonProcess(
                ql.QuoteHandle(ql.SimpleQuote(s)),
                ql.YieldTermStructureHandle(
                    ql.FlatForward(self.todaysDate, q, ql.Actual365Fixed())
                ),
                ql.YieldTermStructureHandle(
                    ql.FlatForward(self.todaysDate, 0.05, ql.Actual365Fixed())
                ),
                ql.BlackVolTermStructureHandle(
                    ql.BlackConstantVol(
                        self.todaysDate, ql.TARGET(), v, ql.Actual365Fixed()
                    )
                ),
            )

        processes = [
            build_process(100, 0.05, 0.3),
            build_process(50, 0.07, 0.45),
            build_process(50, 0.025, 0.2),
        ]

        processes_vector = ql.GeneralizedBlackScholesProcessVector(processes)

        rho = ql.Matrix([[1.0, 0.2, -0.1], [0.2, 1.0, -0.3], [-0.1, -0.3, 1.0]])

        exercise = ql.EuropeanExercise(self.todaysDate + ql.Period(1, ql.Years))
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, 2.0)

        basket_option = ql.BasketOption(
            ql.AverageBasketPayoff(payoff, ql.Array([1, -1, -1])), exercise
        )

        expected = 11.932739641

        basket_option.setPricingEngine(ql.ChoiBasketEngine(processes_vector, rho, 10))
        self.assertAlmostEqual(basket_option.NPV(), expected)

        basket_option.setPricingEngine(ql.DengLiZhouBasketEngine(processes_vector, rho))
        self.assertAlmostEqual(basket_option.NPV(), expected, 1)

        basket_option.setPricingEngine(
            ql.MCEuropeanBasketEngine(
                ql.StochasticProcessArray(processes, rho),
                "lowdiscrepancy",
                timeSteps=1,
                requiredTolerance=0.1,
            )
        )
        self.assertAlmostEqual(basket_option.NPV(), expected, 1)

        basket_option.setPricingEngine(
            ql.FdndimBlackScholesVanillaEngine(
                processes_vector, rho, ql.UnsignedIntVector([25, 15, 15]), 15
            )
        )
        self.assertAlmostEqual(basket_option.NPV(), expected, 1)

    def testTwoAssetSpreadOption(self):
        """Testing two asset spread option"""

        def build_process(s: float, v: float) -> ql.BlackProcess:
            return ql.BlackProcess(
                ql.QuoteHandle(ql.SimpleQuote(s)),
                ql.YieldTermStructureHandle(
                    ql.FlatForward(self.todaysDate, 0.05, ql.Actual365Fixed())
                ),
                ql.BlackVolTermStructureHandle(
                    ql.BlackConstantVol(
                        self.todaysDate, ql.TARGET(), v, ql.Actual365Fixed()
                    )
                ),
            )

        p1 = build_process(100, 0.3)
        p2 = build_process(90, 0.45)
        rho = -0.75
        rho_m = ql.Matrix([[1, rho], [rho, 1]])

        processes_vector = ql.GeneralizedBlackScholesProcessVector([p1, p2])

        exercise = ql.EuropeanExercise(self.todaysDate + ql.Period(6, ql.Months))
        payoff = ql.PlainVanillaPayoff(ql.Option.Put, 10.0)
        basket_option = ql.BasketOption(ql.SpreadBasketPayoff(payoff), exercise)

        expected = 17.96241322097977

        basket_option.setPricingEngine(ql.ChoiBasketEngine(processes_vector, rho_m, 15))
        self.assertAlmostEqual(basket_option.NPV(), expected, 10)

        basket_option.setPricingEngine(
            ql.DengLiZhouBasketEngine(processes_vector, rho_m)
        )
        self.assertAlmostEqual(basket_option.NPV(), expected, 4)

        basket_option.setPricingEngine(ql.KirkEngine(p1, p2, rho))
        self.assertAlmostEqual(basket_option.NPV(), expected, 1)

        basket_option.setPricingEngine(ql.BjerksundStenslandSpreadEngine(p1, p2, rho))
        self.assertAlmostEqual(basket_option.NPV(), expected, 2)

        basket_option.setPricingEngine(
            ql.OperatorSplittingSpreadEngine(
                p1, p2, rho, ql.OperatorSplittingSpreadEngine.First
            )
        )
        self.assertAlmostEqual(basket_option.NPV(), expected, 1)

        basket_option.setPricingEngine(
            ql.OperatorSplittingSpreadEngine(
                p1, p2, rho, ql.OperatorSplittingSpreadEngine.Second
            )
        )
        self.assertAlmostEqual(basket_option.NPV(), expected, 2)

        basket_option.setPricingEngine(
            ql.FdndimBlackScholesVanillaEngine(
                processes_vector, rho_m, ql.UnsignedIntVector([25, 25]), 15
            )
        )
        self.assertAlmostEqual(basket_option.NPV(), expected, 1)

        basket_option.setPricingEngine(
            ql.Fd2dBlackScholesVanillaEngine(p1, p2, rho, xGrid=25, yGrid=25, tGrid=15)
        )
        self.assertAlmostEqual(basket_option.NPV(), expected, 1)

        basket_option.setPricingEngine(
            ql.MCEuropeanBasketEngine(
                ql.StochasticProcessArray([p1, p2], rho_m),
                "lowdiscrepancy",
                timeSteps=1,
                requiredTolerance=0.1,
            )
        )
        self.assertAlmostEqual(basket_option.NPV(), expected, 1)


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
