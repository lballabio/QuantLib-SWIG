"""
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl

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


class MarketElementTest(unittest.TestCase):
    def testObservable(self):
        "Testing observability of market elements"
        global flag
        flag = None
        me = ql.SimpleQuote(0.0)
        obs = ql.Observer(raiseFlag)
        obs.registerWith(me)
        me.setValue(3.14)
        if not flag:
            self.fail("Observer was not notified of market element change")
        flag = None
        obs.unregisterWith(me)
        me.setValue(2.71)
        if flag:
            self.fail("Observer was notified after unregistering")

    def testObservableHandle(self):
        "Testing observability of market element handles"
        global flag
        flag = None
        me1 = ql.SimpleQuote(0.0)
        h = ql.RelinkableQuoteHandle(me1)
        obs = ql.Observer(raiseFlag)
        obs.registerWith(h)
        me1.setValue(3.14)
        if not flag:
            self.fail("Observer was not notified of market element change")
        flag = None
        me2 = ql.SimpleQuote(0.0)
        h.linkTo(me2)
        if not flag:
            self.fail("Observer was not notified of market element change")
        flag = None
        obs.unregisterWith(h)
        me2.setValue(2.71)
        if flag:
            self.fail("Observer was notified after unregistering")

    def testObservableErrors(self):
        class Handle:
            def __init__(self, x):
                self.x = x

            def asObservable(self):
                return 10 / self.x

        obs = ql.Observer(raiseFlag)
        self.assertRaises(TypeError, obs.registerWith, 123)
        self.assertRaises(TypeError, obs.registerWith, obs)
        self.assertRaises(TypeError, obs.registerWith, Handle)
        self.assertRaises(ZeroDivisionError, obs.registerWith, Handle(0))
        self.assertRaises(TypeError, obs.registerWith, Handle(1))

    def test_SimpleQuote(self):
        for value in (100, 2.71, 10**100):
            self.assertAlmostEqual(ql.SimpleQuote(value).value(), value)


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
