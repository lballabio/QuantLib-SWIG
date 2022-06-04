"""
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2007 StatPro Italia srl
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

import QuantLib as ql
import unittest
import math

flag = None


def raiseFlag():
    global flag
    flag = 1


def binaryFunction(x, y):
    return 2.0 * x + y


def extrapolatedForwardRate(
        firstSmoothingPoint,
        lastLiquidForwardRate,
        ultimateForwardRate,
        alpha):

    def calculate(t):
        deltaT = t - firstSmoothingPoint
        beta = (1.0 - math.exp(-alpha * deltaT)) / (alpha * deltaT)
        return ultimateForwardRate + (
            lastLiquidForwardRate - ultimateForwardRate) * beta

    return calculate


class TermStructureTest(unittest.TestCase):
    def setUp(self):
        self.calendar = ql.TARGET()
        today = self.calendar.adjust(ql.Date.todaysDate())
        self.settlementDays = 2
        self.dayCounter = ql.Actual360()
        settlement = self.calendar.advance(today, self.settlementDays, ql.Days)
        deposits = [
            ql.DepositRateHelper(
                ql.QuoteHandle(ql.SimpleQuote(rate / 100)),
                ql.Period(n, units),
                self.settlementDays,
                self.calendar,
                ql.ModifiedFollowing,
                False,
                self.dayCounter,
            )
            for (n, units, rate) in [
                (1, ql.Months, 4.581),
                (2, ql.Months, 4.573),
                (3, ql.Months, 4.557),
                (6, ql.Months, 4.496),
                (9, ql.Months, 4.490),
            ]
        ]
        swaps = [
            ql.SwapRateHelper(
                ql.QuoteHandle(ql.SimpleQuote(rate / 100)),
                ql.Period(years, ql.Years),
                self.calendar,
                ql.Annual,
                ql.Unadjusted,
                ql.Thirty360(),
                ql.Euribor6M(),
            )
            for (years, rate) in [(1, 4.54), (5, 4.99), (10, 5.47), (20, 5.89), (30, 5.96)]
        ]

        self.termStructure = ql.PiecewiseFlatForward(
            settlement, deposits + swaps, self.dayCounter)

    def testImpliedObs(self):
        "Testing observability of implied term structure"
        global flag
        flag = None
        h = ql.RelinkableYieldTermStructureHandle()
        settlement = self.termStructure.referenceDate()
        new_settlement = self.calendar.advance(settlement, 3, ql.Years)
        implied = ql.ImpliedTermStructure(h, new_settlement)
        obs = ql.Observer(raiseFlag)
        obs.registerWith(implied)
        h.linkTo(self.termStructure)
        if not flag:
            self.fail("Observer was not notified of term structure change")

    def testFSpreadedObs(self):
        "Testing observability of forward-spreaded term structure"
        global flag
        flag = None
        me = ql.SimpleQuote(0.01)
        mh = ql.QuoteHandle(me)
        h = ql.RelinkableYieldTermStructureHandle()
        spreaded = ql.ForwardSpreadedTermStructure(h, mh)
        obs = ql.Observer(raiseFlag)
        obs.registerWith(spreaded)
        h.linkTo(self.termStructure)
        if not flag:
            self.fail("Observer was not notified of term structure change")
        flag = None
        me.setValue(0.005)
        if not flag:
            self.fail("Observer was not notified of spread change")

    def testZSpreadedObs(self):
        "Testing observability of zero-spreaded term structure"
        global flag
        flag = None
        me = ql.SimpleQuote(0.01)
        mh = ql.QuoteHandle(me)
        h = ql.RelinkableYieldTermStructureHandle()
        spreaded = ql.ZeroSpreadedTermStructure(h, mh)
        obs = ql.Observer(raiseFlag)
        obs.registerWith(spreaded)
        h.linkTo(self.termStructure)
        if not flag:
            self.fail("Observer was not notified of term structure change")
        flag = None
        me.setValue(0.005)
        if not flag:
            self.fail("Observer was not notified of spread change")

    def testCompositeZeroYieldStructure(self):
        """Testing composite zero yield structure"""
        settlement = self.termStructure.referenceDate()
        compounding = ql.Compounded
        freq = ql.Semiannual
        flatTs = ql.FlatForward(
            settlement,
            ql.QuoteHandle(ql.SimpleQuote(0.0085)),
            self.dayCounter)
        firstHandle = ql.YieldTermStructureHandle(flatTs)
        secondHandle = ql.YieldTermStructureHandle(self.termStructure)
        compositeTs = ql.CompositeZeroYieldStructure(
            firstHandle, secondHandle, binaryFunction, compounding, freq)
        maturity = settlement + ql.Period(20, ql.Years)
        expectedZeroRate = binaryFunction(
            firstHandle.zeroRate(
                maturity, self.dayCounter, compounding, freq).rate(),
            secondHandle.zeroRate(
                maturity, self.dayCounter, compounding, freq).rate())
        actualZeroRate = compositeTs.zeroRate(
            maturity, self.dayCounter, compounding, freq).rate()
        failMsg = """ Composite zero yield structure rate replication failed:
                        expected zero rate: {expected}
                        actual zero rate: {actual}
                  """.format(expected=expectedZeroRate,
                             actual=actualZeroRate)
        self.assertAlmostEquals(
            first=expectedZeroRate,
            second=actualZeroRate,
            delta=1.0e-12,
            msg=failMsg)

    def testUltimateForwardTermStructure(self):
        """Testing ultimate forward term structure"""
        settlement = self.termStructure.referenceDate()
        ufr = ql.QuoteHandle(ql.SimpleQuote(0.06))
        llfr = ql.QuoteHandle(ql.SimpleQuote(0.05))
        fsp = ql.Period(20, ql.Years)
        alpha = 0.05
        baseCrvHandle = ql.YieldTermStructureHandle(self.termStructure)
        ufrCrv = ql.UltimateForwardTermStructure(
            baseCrvHandle, llfr, ufr, fsp, alpha)
        cutOff = ufrCrv.timeFromReference(settlement + fsp)
        forwardCalculator = extrapolatedForwardRate(
            cutOff, llfr.value(), ufr.value(), alpha)
        times = [ufrCrv.timeFromReference(settlement + ql.Period(x, ql.Years))
                 for x in [21, 30, 40, 50, 60, 70, 80, 90, 100]]
        for t in times:
            actualForward = ufrCrv.forwardRate(
                cutOff, t, ql.Continuous, ql.NoFrequency, True).rate()
            expectedForward = forwardCalculator(t)
            failMsg = """ UFR term structure forward replication failed for:
                            time to maturity: {timeToMaturity}
                            expected forward rate: {expected}
                            actual forward rate: {actual}
                      """.format(timeToMaturity=t,
                                 expected=expectedForward,
                                 actual=actualForward)
            self.assertAlmostEquals(
                first=expectedForward,
                second=actualForward,
                delta=1.0e-12,
                msg=failMsg)

    def testQuantoTermStructure(self):
        """Testing quanto term structure"""
        today = ql.Date.todaysDate()

        dividend_ts = ql.YieldTermStructureHandle(
            ql.FlatForward(
                today,
                ql.QuoteHandle(ql.SimpleQuote(0.055)),
                self.dayCounter
            )
        )
        r_domestic_ts = ql.YieldTermStructureHandle(
            ql.FlatForward(
                today,
                ql.QuoteHandle(ql.SimpleQuote(-0.01)),
                self.dayCounter
            )
        )
        r_foreign_ts = ql.YieldTermStructureHandle(
            ql.FlatForward(
                today,
                ql.QuoteHandle(ql.SimpleQuote(0.02)),
                self.dayCounter
            )
        )
        sigma_s = ql.BlackVolTermStructureHandle(
            ql.BlackConstantVol(
                today,
                self.calendar,
                ql.QuoteHandle(ql.SimpleQuote(0.25)),
                self.dayCounter
            )
        )
        sigma_fx = ql.BlackVolTermStructureHandle(
            ql.BlackConstantVol(
                today,
                self.calendar,
                ql.QuoteHandle(ql.SimpleQuote(0.05)),
                self.dayCounter
            )
        )
        rho = ql.QuoteHandle(ql.SimpleQuote(0.3))
        s_0 = ql.QuoteHandle(ql.SimpleQuote(100.0))

        exercise = ql.EuropeanExercise(self.calendar.advance(today, 6, ql.Months))
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, 95.0)

        vanilla_option = ql.VanillaOption(payoff, exercise)
        quanto_ts = ql.YieldTermStructureHandle(
            ql.QuantoTermStructure(
                dividend_ts,
                r_domestic_ts,
                r_foreign_ts,
                sigma_s,
                ql.nullDouble(),
                sigma_fx,
                ql.nullDouble(),
                rho.value()
            )
        )
        gbm_quanto = ql.BlackScholesMertonProcess(s_0, quanto_ts, r_domestic_ts, sigma_s)
        vanilla_engine = ql.AnalyticEuropeanEngine(gbm_quanto)
        vanilla_option.setPricingEngine(vanilla_engine)

        quanto_option = ql.QuantoVanillaOption(payoff, exercise)
        gbm_vanilla = ql.BlackScholesMertonProcess(s_0, dividend_ts, r_domestic_ts, sigma_s)
        quanto_engine = ql.QuantoEuropeanEngine(gbm_vanilla, r_foreign_ts, sigma_fx, rho)
        quanto_option.setPricingEngine(quanto_engine)

        quanto_option_pv = quanto_option.NPV()
        vanilla_option_pv = vanilla_option.NPV()

        message = """Failed to reproduce QuantoOption / EuropeanQuantoEngine NPV:
                      {quanto_pv}
                      by using the QuantoTermStructure as the dividend together with
                      VanillaOption / AnalyticEuropeanEngine:
                      {vanilla_pv}
                  """.format(
            quanto_pv=quanto_option_pv,
            vanilla_pv=vanilla_option_pv
        )

        self.assertAlmostEquals(
            quanto_option_pv,
            vanilla_option_pv,
            delta=1e-12,
            msg=message
        )


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TermStructureTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
