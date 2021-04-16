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

import unittest

import QuantLib as ql


class SlvTest(unittest.TestCase):
    def setUp(self):
        self.todaysDate = ql.Date(15, ql.May, 2019)
        ql.Settings.instance().evaluationDate = self.todaysDate
        self.settlementDate = self.todaysDate + ql.Period(2, ql.Days)
        self.dc = ql.Actual365Fixed()
        self.riskFreeRate = ql.YieldTermStructureHandle(ql.FlatForward(self.settlementDate, 0.05, self.dc))
        self.dividendYield = ql.YieldTermStructureHandle(ql.FlatForward(self.settlementDate, 0.025, self.dc))
        self.underlying = ql.QuoteHandle(ql.SimpleQuote(100.0))

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def constVol(self, vol):
        return ql.BlackVolTermStructureHandle(ql.BlackConstantVol(self.settlementDate, ql.TARGET(), vol, self.dc))

    def testSlvProcess(self):
        """ Testing HestonSLVProcess generation """

        hestonProcess = ql.HestonProcess(
            self.riskFreeRate, self.riskFreeRate, self.underlying, 0.1 * 0.1, 1.0, 0.25 * 0.25, 0.15, -0.75
        )

        localVol = ql.LocalVolSurface(
            ql.BlackVolTermStructureHandle(ql.BlackConstantVol(self.settlementDate, ql.TARGET(), 0.10, self.dc)),
            self.riskFreeRate,
            self.riskFreeRate,
            self.underlying,
        )

        ql.HestonSLVProcess(hestonProcess, localVol)

    def testSlvProcessAsBlackScholes(self):
        """ Testing HestonSLVProcess equal to Black-Scholes process """

        hestonProcess = ql.HestonProcess(
            self.riskFreeRate, self.dividendYield, self.underlying, 0.01, 1.0, 0.01, 1e-4, 0.0
        )

        exercise = ql.EuropeanExercise(self.todaysDate + ql.Period(1, ql.Years))
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, self.underlying.value())

        option = ql.VanillaOption(payoff, exercise)

        hestonModel = ql.HestonModel(hestonProcess)
        option.setPricingEngine(ql.FdHestonVanillaEngine(hestonModel, 20, 100, 3))

        hestonNPV = option.NPV()

        option.setPricingEngine(
            ql.AnalyticEuropeanEngine(
                ql.BlackScholesMertonProcess(self.underlying, self.dividendYield, self.riskFreeRate, self.constVol(0.1))
            )
        )

        bsNPV = option.NPV()

        self.assertAlmostEqual(
            hestonNPV, bsNPV, 2, msg="Unable to reproduce Heston vanilla option price with Black-Scholes process"
        )

        leverageFct = ql.LocalVolSurface(self.constVol(2.0), self.riskFreeRate, self.dividendYield, self.underlying)

        option.setPricingEngine(
            ql.FdHestonVanillaEngine(
                hestonModel,
                20,
                100,
                3,
                1,
                ql.FdmSchemeDesc.Hundsdorfer(),
                leverageFct,
            )
        )

        slvNPV = option.NPV()

        bsmProcess = ql.BlackScholesMertonProcess(self.underlying, self.dividendYield, self.riskFreeRate, self.constVol(0.2))

        option.setPricingEngine(ql.AnalyticEuropeanEngine(bsmProcess))

        bsNPV = option.NPV()

        self.assertAlmostEqual(
            slvNPV,
            bsNPV,
            2,
            msg="Unable to reproduce Heston plus constant local vol option price with Black-Scholes formula",
        )

        barrier_lo = 70.0
        barrier_hi = 130.0

        barrierOption = ql.DoubleBarrierOption(
            ql.DoubleBarrier.KnockOut,
            barrier_lo,
            barrier_hi,
            0.0,
            ql.CashOrNothingPayoff(ql.Option.Call, 0.0, 1.0),
            exercise);

        barrierOption.setPricingEngine(
            ql.FdHestonDoubleBarrierEngine(
                hestonModel,
                400,
                100,
                2,
                1,
                ql.FdmSchemeDesc.Hundsdorfer(),
                leverageFct,
            )
        )

        slvBarrierNPV = barrierOption.NPV()

        barrierOption.setPricingEngine(ql.AnalyticDoubleBarrierBinaryEngine(bsmProcess))

        bsmBarrierNPV = barrierOption.NPV()

        self.assertAlmostEqual(
            slvBarrierNPV,
            bsmBarrierNPV,
            2,
            msg="Unable to reproduce Heston plus constant local vol "
            "double barrier option price with Black-Scholes Double Barrier Binary Engine",
        )



if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(SlvTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
