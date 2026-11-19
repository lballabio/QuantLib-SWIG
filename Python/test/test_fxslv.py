"""
 Copyright (C) 2026 Raghav Subbarao

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

# Tests of the wrappers of ql/experimental/fxslv (smile quotes, targets,
# sections, factory and sticky-delta solvers). Market data, limits and
# reference values are those of test-suite/fxsmilesections.cpp, so the
# wrapped classes are checked against the C++ tests.

import math
import unittest

import QuantLib as ql


class MarketData:
    """MarketData of test-suite/fxsmilesections.cpp."""

    def __init__(self):
        self.today = ql.Date(2, ql.January, 2024)
        self.expiry = ql.Date(2, ql.January, 2025)
        ql.Settings.instance().evaluationDate = self.today
        self.forDiscount = ql.YieldTermStructureHandle(ql.FlatForward(self.today, 0.05, ql.Actual365Fixed()))
        self.domDiscount = ql.YieldTermStructureHandle(ql.FlatForward(self.today, 0.03, ql.Actual365Fixed()))
        self.spotQuote = ql.SimpleQuote(1.7554)
        self.spot = ql.QuoteHandle(self.spotQuote)
        self.atm = ql.SimpleQuote(0.14483)
        self.rr = [ql.SimpleQuote(0.02770), ql.SimpleQuote(0.048752)]
        self.bf = [ql.SimpleQuote(0.007425), ql.SimpleQuote(0.02376)]
        self.deltas = [0.25, 0.1]
        self.deltaType = ql.DeltaVolQuote.PaSpot
        self.atmType = ql.DeltaVolQuote.AtmFwd
        self.settlement = ql.FxSettleConvention(ql.WeekendsOnly(), 2)

    def quotes(self, flyType=ql.FxRrBfQuotes.SmileStrangle):
        return ql.FxRrBfQuotes(ql.QuoteHandle(self.atm),
                               [ql.QuoteHandle(q) for q in self.rr],
                               [ql.QuoteHandle(q) for q in self.bf],
                               self.deltas, flyType)

    def section(self, cls, deltaType=None, quotes=None):
        return cls(self.expiry, self.spot, quotes or self.quotes(), self.forDiscount, self.domDiscount,
                   deltaType if deltaType is not None else self.deltaType, self.atmType,
                   ql.Actual365Fixed(), self.settlement)


# model -> (class, calibration error limit of the C++ test)
MODELS = {
    ql.FxSmileModel.Polynomial: (ql.FxPolynomialSmileSection, 15.0e-4),
    ql.FxSmileModel.Sabr: (ql.FxSabrSmileSection, 25.0e-4),
    ql.FxSmileModel.Svi: (ql.FxSviSmileSection, 1.0e-6),
    ql.FxSmileModel.Quadratic: (ql.FxQuadraticSmileSection, 30.0e-4),
    ql.FxSmileModel.CostFlatDynamics: (ql.FxCostSmileSectionFlatDynamics, 17.0e-4),
    ql.FxSmileModel.CostScaledDynamics: (ql.FxCostSmileSectionScaledDynamics, 35.0e-4),
    ql.FxSmileModel.CostStickyDelta: (ql.FxCostSmileSectionStickyDelta, 27.0e-4),
    ql.FxSmileModel.CostStickyDeltaConvex: (ql.FxCostSmileSectionStickyDeltaConvex, 27.0e-4),
}

STICKY_COSTS = [1.0, -15.3166, -0.06278, 0.28742]
XS = [-0.08, -0.03, 0.0, 0.03, 0.08]


class FxSmileSectionTest(unittest.TestCase):

    def setUp(self):
        self.md = MarketData()

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()

    def testEveryModelCalibrates(self):
        """Every model fits the five quotes within the C++ test's limit."""
        for model, (cls, limit) in MODELS.items():
            with self.subTest(model=cls.__name__):
                s = self.md.section(cls)
                self.assertEqual(len(s.calibrationResiduals()), 5)
                self.assertLess(s.calibrationError(), limit)
                self.assertAlmostEqual(s.atmLevel(), s.forward(), places=14)

    def testStrikesAndDeltasInvert(self):
        for model, (cls, _) in MODELS.items():
            s = self.md.section(cls)
            for d in self.md.deltas:
                for parity, delta in ((ql.Option.Call, d), (ql.Option.Put, -d)):
                    with self.subTest(model=cls.__name__, delta=delta):
                        k = s.strikeByDelta(delta, parity)
                        self.assertAlmostEqual(s.deltaByStrike(k, parity), delta, delta=1e-5)
                        self.assertAlmostEqual(s.volByDelta(delta, parity), s.volByStrike(k), delta=1e-10)

    def testFactoryBuildsEachModel(self):
        """makeFxSmileSection gives the named class, and downcasts reach it."""
        casts = {
            ql.FxSmileModel.Polynomial: ql.as_fx_polynomial_smile_section,
            ql.FxSmileModel.Sabr: ql.as_fx_sabr_smile_section,
            ql.FxSmileModel.Svi: ql.as_fx_svi_smile_section,
            ql.FxSmileModel.Quadratic: ql.as_fx_quadratic_smile_section,
            ql.FxSmileModel.CostFlatDynamics: ql.as_fx_cost_smile_section_flat_dynamics,
            ql.FxSmileModel.CostScaledDynamics: ql.as_fx_cost_smile_section_scaled_dynamics,
            ql.FxSmileModel.CostStickyDelta: ql.as_fx_cost_smile_section_sticky_delta,
            ql.FxSmileModel.CostStickyDeltaConvex: ql.as_fx_cost_smile_section_sticky_delta_convex,
        }
        md = self.md
        for model, (cls, _) in MODELS.items():
            with self.subTest(model=cls.__name__):
                built = ql.makeFxSmileSection(model, md.expiry, md.spot, md.quotes(), md.forDiscount,
                                              md.domDiscount, md.deltaType, md.atmType, md.settlement)
                self.assertIsNotNone(casts[model](built))
                direct = md.section(cls)
                for m in (0.8, 0.9, 1.0, 1.1, 1.2):
                    k = m * built.forward()
                    self.assertAlmostEqual(built.volByStrike(k), direct.volByStrike(k), delta=1e-12)
        sabr = md.section(ql.FxSabrSmileSection)
        self.assertIsNone(ql.as_fx_svi_smile_section(sabr))

    def testCostSectionAccessors(self):
        s = self.md.section(ql.FxCostSmileSectionStickyDelta)
        self.assertEqual(len(s.costs()), 4)
        self.assertTrue(s.weightedCalibration())
        self.assertEqual(ql.FxCostSmileSectionStickyDelta.smileRange(), 8.0)

    def testStickySlopeIsTheSmileSlope(self):
        """Mirrors testCostSmileSectionStickyDelta."""
        for deltaType in (ql.DeltaVolQuote.PaSpot, ql.DeltaVolQuote.Spot):
            s = self.md.section(ql.FxCostSmileSectionStickyDelta, deltaType)
            tau = s.exerciseTime()
            for m in (0.8, 0.9, 1.0, 1.1, 1.2):
                with self.subTest(deltaType=deltaType, m=m):
                    k, h = m * s.forward(), 1.0e-5
                    fd = (s.volByStrike(k * math.exp(h)) - s.volByStrike(k * math.exp(-h))) * math.sqrt(tau) / (2 * h)
                    self.assertAlmostEqual(s.slope(k), fd, delta=1e-6)
                    dk = 1.0e-5 * k
                    fdk = (s.volByStrike(k + dk) - s.volByStrike(k - dk)) / (2 * dk)
                    self.assertAlmostEqual(s.volDerivative(k), fdk, delta=1e-6)

    def testSectionReactsToQuotes(self):
        """Python-side SimpleQuotes drive recalibration through the wrapped handles."""
        s = self.md.section(ql.FxPolynomialSmileSection)
        f, v = s.forward(), s.volByStrike(0.9 * s.forward())
        self.md.spotQuote.setValue(1.80)
        self.assertNotAlmostEqual(s.forward(), f, places=6)
        k = 0.9 * s.forward()
        before = s.volByStrike(k)
        self.md.bf[0].setValue(self.md.bf[0].value() + 0.002)
        self.assertNotAlmostEqual(s.volByStrike(k), before, places=8)
        self.assertNotAlmostEqual(v, 0.0)

    def testMarketStrangleQuotes(self):
        quotes = self.md.quotes(ql.FxRrBfQuotes.MarketStrangle)
        self.assertEqual(quotes.flyType(), ql.FxRrBfQuotes.MarketStrangle)
        s = self.md.section(ql.FxSviSmileSection, quotes=quotes)
        self.assertLess(s.calibrationError(), 1e-6)

    def testDeltaVolQuotes(self):
        """FxDeltaVolQuotes read off a calibrated smile give the same smile back."""
        md = self.md
        base = md.section(ql.FxSviSmileSection)
        t, dt = base.exerciseTime(), base.deltaType()
        q = [ql.DeltaVolQuoteHandle(ql.DeltaVolQuote(ql.QuoteHandle(ql.SimpleQuote(base.atmVol())), dt, t,
                                                     base.atmType()))]
        for d in (0.25, 0.10):
            for parity, delta in ((ql.Option.Call, d), (ql.Option.Put, -d)):
                vol = ql.QuoteHandle(ql.SimpleQuote(base.volByDelta(delta, parity)))
                q.append(ql.DeltaVolQuoteHandle(ql.DeltaVolQuote(delta, vol, t, dt)))
        quotes = ql.FxDeltaVolQuotes(q)
        self.assertEqual(len(quotes.quotes()), 5)
        again = md.section(ql.FxSviSmileSection, quotes=quotes)
        for m in (0.9, 1.0, 1.1):
            k = m * base.forward()
            self.assertAlmostEqual(again.volByStrike(k), base.volByStrike(k), delta=1e-6)

    def testTargets(self):
        md = self.md
        s = md.section(ql.FxSviSmileSection)
        dvq = ql.DeltaVolQuoteHandle(ql.DeltaVolQuote(0.25, ql.QuoteHandle(ql.SimpleQuote(0.16)), s.exerciseTime(),
                                                      s.deltaType()))
        target = ql.FxDeltaVolTarget(dvq)
        self.assertTrue(target.hasPoint(s))
        strike, vol = target.point(s)
        self.assertEqual(vol, 0.16)
        self.assertAlmostEqual(target.residual(s), s.volResidual(strike, vol), places=14)
        self.assertGreater(target.weight(s), 0.0)
        rr = ql.FxRiskReversalTarget(0.25, 0.0277)
        self.assertFalse(rr.hasPoint(s))
        with self.assertRaises(RuntimeError):
            rr.point(s)
        targets = ql.FxSmileTargets([target, rr, ql.FxBrokerStrangleTarget(0.14483, 0.007425, 0.25)])
        self.assertEqual(len(targets), 3)

    def testDateAndTimeModes(self):
        md = self.md
        s = md.section(ql.FxSabrSmileSection)
        self.assertTrue(s.hasSettleConvention())
        self.assertEqual(s.settleConvention().spotLag(), 2)
        self.assertEqual(s.spotDate(), md.settlement.valueDate(md.today))
        self.assertEqual(s.deliveryDate(), md.settlement.valueDate(md.expiry))
        timed = ql.FxSabrSmileSection(1.0, md.spot, md.quotes(), md.forDiscount, md.domDiscount, md.deltaType,
                                      md.atmType)
        self.assertFalse(timed.hasSettleConvention())
        with self.assertRaises(RuntimeError):
            timed.settleConvention()
        expected = md.spot.value() * md.forDiscount.discount(1.0) / md.domDiscount.discount(1.0)
        self.assertAlmostEqual(timed.forward(), expected, places=12)

    def testDeltaConverter(self):
        s = self.md.section(ql.FxSabrSmileSection)
        c = s.deltaConverter()
        self.assertTrue(c.premiumAdjusted())
        self.assertAlmostEqual(c.forward(), s.forward(), places=12)
        w = 0.15
        for parity, delta in ((ql.Option.Call, 0.25), (ql.Option.Put, -0.25)):
            k = c.strike(parity, delta, w)
            self.assertAlmostEqual(c.delta(parity, k, w), delta, places=10)
        k = 1.7
        self.assertAlmostEqual(c.delta(ql.Option.Call, k, w) - c.delta(ql.Option.Put, k, w),
                               c.callMinusPutDelta(k), places=12)
        self.assertFalse(c.withType(ql.DeltaVolQuote.Fwd).premiumAdjusted())


class FxStickyDeltaSolverTest(unittest.TestCase):

    def testFlatSolverReferenceValues(self):
        """Reference values of testStickyDeltaSmileSolver."""
        ref = {False: [0.040779300983, 0.032947238825, 0.032620876353, 0.035379624984, 0.045066283516],
               True: [0.040868871525, 0.033074416691, 0.032755032154, 0.035508341595, 0.045174432271]}
        for pa in (False, True):
            solver = ql.FxStickyDeltaSmileSolver(STICKY_COSTS, pa, -0.1, 0.1)
            ws = solver.singularTotalVol()
            expected = 1.0 / (2.0 * (STICKY_COSTS[2] - STICKY_COSTS[1])) if pa else -1.0 / (2.0 * STICKY_COSTS[1])
            self.assertAlmostEqual(ws, expected, places=14)
            self.assertAlmostEqual(solver.singularLogMoneyness(), ws * ws / 2, places=14)
            for x, w in zip(XS, ref[pa]):
                self.assertAlmostEqual(solver.totalVol(x), w, delta=1e-9)
            self.assertEqual(solver.convexity(0.0), 0.0)
            self.assertTrue(solver.hasCosts([0.37 * c for c in STICKY_COSTS]))
            self.assertFalse(solver.hasCosts([-c for c in STICKY_COSTS]))
        with self.assertRaises(RuntimeError):
            ql.FxStickyDeltaSmileSolver([1.0, 0.5, -0.06, 0.29], False, -0.1, 0.1)

    def testConvexSolverReferenceValues(self):
        """Reference values of testStickyDeltaConvexSmileSolver (scipy solve_bvp)."""
        refW = {False: [0.043161310158, 0.036456627963, 0.036126201696, 0.038485343052, 0.046943186387],
                True: [0.043273312194, 0.036599102919, 0.036274554535, 0.038628792457, 0.047067508358]}
        for pa in (False, True):
            flat = ql.FxStickyDeltaSmileSolver(STICKY_COSTS, pa, -0.1, 0.1)
            lo = flat.singularLogMoneyness() - 8.0 * flat.singularTotalVol()
            hi = flat.singularLogMoneyness() + 8.0 * flat.singularTotalVol()
            smile = ql.FxStickyDeltaConvexSmileSolver(STICKY_COSTS, pa, lo, hi)
            self.assertAlmostEqual(smile.minLogMoneyness(), lo, places=14)
            for x, w in zip(XS, refW[pa]):
                self.assertAlmostEqual(smile.totalVol(x), w, delta=1e-8)
            for x in XS:  # kappa = -F / (2 w C_gamma), F the cost equation without convexity
                w, beta = smile.totalVol(x), smile.slope(x)
                self.assertAlmostEqual(smile.convexity(x), -smile.costEquation(x, w, beta) / (2 * w), delta=1e-8)


if __name__ == "__main__":
    print("testing QuantLib", ql.__version__)
    unittest.main(verbosity=2)
