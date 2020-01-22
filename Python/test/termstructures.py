"""
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2007 StatPro Italia srl

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

flag = None


def raiseFlag():
    global flag
    flag = 1


class TermStructureTest(unittest.TestCase):
    def setUp(self):
        self.calendar = ql.TARGET()
        today = self.calendar.adjust(ql.Date.todaysDate())
        self.settlementDays = 2
        settlement = self.calendar.advance(today, self.settlementDays, ql.Days)
        deposits = [
            ql.DepositRateHelper(
                ql.QuoteHandle(ql.SimpleQuote(rate / 100)),
                ql.Period(n, units),
                self.settlementDays,
                self.calendar,
                ql.ModifiedFollowing,
                False,
                ql.Actual360(),
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

        self.termStructure = ql.PiecewiseFlatForward(settlement, deposits + swaps, ql.Actual360())

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


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TermStructureTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
