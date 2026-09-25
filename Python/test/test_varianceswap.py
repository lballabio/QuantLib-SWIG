"""
 Copyright (C) 2026 Mahimn Patel

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
"""

import unittest

import QuantLib as ql


class VarianceSwapTest(unittest.TestCase):
    def setUp(self):
        # data from "A Guide to Volatility and Variance Swaps", Derman, Kamal & Zou, 1999,
        # with maturity corrected from 0.25 to 0.246575 (Jan 1, 1999 to Apr 1, 1999)
        self.today = ql.Date(1, ql.January, 1999)
        ql.Settings.instance().evaluationDate = self.today
        self.maturity = self.today + 90
        self.day_counter = ql.Actual365Fixed()

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def make_process(self, vol_ts):
        return ql.BlackScholesMertonProcess(
            ql.QuoteHandle(ql.SimpleQuote(100.0)),
            ql.YieldTermStructureHandle(ql.FlatForward(self.today, 0.0, self.day_counter)),
            ql.YieldTermStructureHandle(ql.FlatForward(self.today, 0.05, self.day_counter)),
            ql.BlackVolTermStructureHandle(vol_ts),
        )

    def testReplicatingVarianceSwap(self):
        """Testing variance swap with replicating cost engine"""
        put_strikes = [50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0]
        put_vols = [0.30, 0.29, 0.28, 0.27, 0.26, 0.25, 0.24, 0.23, 0.22, 0.21, 0.20]
        call_strikes = [100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0]
        call_vols = [0.20, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13]

        # the put and call strikes share the at-the-money point
        vol_ts = ql.BlackVarianceSurface(
            self.today,
            ql.NullCalendar(),
            [self.maturity],
            put_strikes + call_strikes[1:],
            ql.Matrix([[v] for v in put_vols + call_vols[1:]]),
            self.day_counter,
        )
        engine = ql.ReplicatingVarianceSwapEngine(
            self.make_process(vol_ts), 5.0, call_strikes, put_strikes
        )
        swap = ql.VarianceSwap(ql.Position.Long, 0.04, 50000.0, self.today, self.maturity)
        swap.setPricingEngine(engine)

        self.assertAlmostEqual(swap.variance(), 0.04189, delta=1.0e-4)

    def testMCVarianceSwap(self):
        """Testing variance swap with Monte Carlo engine"""
        # with a variance curve, the fair variance is v*v for any intermediate
        # (t1, v1) with 0 <= t1 < t and 0 <= v1 < v
        vol_ts = ql.BlackVarianceCurve(
            self.today, [self.today + 37, self.maturity], [0.1, 0.2], self.day_counter, True
        )
        engine = ql.MCVarianceSwapEngine(
            self.make_process(vol_ts),
            "pseudorandom",
            timeStepsPerYear=250,
            requiredSamples=1023,
            seed=42,
        )
        swap = ql.VarianceSwap(ql.Position.Long, 0.04, 50000.0, self.today, self.maturity)
        swap.setPricingEngine(engine)

        self.assertAlmostEqual(swap.variance(), 0.04, delta=3.0e-4)


if __name__ == "__main__":
    unittest.main()
