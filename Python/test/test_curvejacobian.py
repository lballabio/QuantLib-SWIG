"""
 Copyright (C) 2026 Kyrylo Protsenko

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


BUMP = 1.0e-6


def node_derivatives(curves, quote, bump=BUMP):
    """Finite differences of every curve's nodes to one quote."""
    q0 = quote.value()
    quote.setValue(q0 + bump)
    up = [list(c.data()) for c in curves]
    quote.setValue(q0 - bump)
    down = [list(c.data()) for c in curves]
    quote.setValue(q0)
    # the data()[0] entry is the value at the reference date, which is
    # not a free variable and is dropped from the Jacobians
    return [
        [(u[j] - d[j]) / (2.0 * bump) for j in range(1, len(u))]
        for u, d in zip(up, down)
    ]


class CurveJacobianTest(unittest.TestCase):
    def setUp(self):
        self.today = ql.Date(23, ql.October, 2025)
        ql.Settings.instance().evaluationDate = self.today

    def build_ois_curve(self, curve_type=ql.PiecewiseLogLinearDiscount):
        """ESTR OIS curve, discounting itself."""
        quotes, helpers = [], []
        estr = ql.Estr()
        for months, rate in [(3, 0.0195), (12, 0.0210), (24, 0.0225),
                             (60, 0.0250), (120, 0.0260)]:
            q = ql.SimpleQuote(rate)
            quotes.append(q)
            helpers.append(
                ql.OISRateHelper(2, ql.Period(months, ql.Months),
                                 ql.QuoteHandle(q), estr)
            )
        curve = curve_type(
            0, ql.TARGET(), helpers, ql.Actual360()
        )
        return curve, quotes, helpers

    def build_projection_curve(self, discount_curve):
        """6M curve bootstrapped with exogenous OIS discounting."""
        euribor6m = ql.Euribor6M()
        discount_handle = ql.YieldTermStructureHandle(discount_curve)
        quotes, helpers = [], []

        depo = ql.SimpleQuote(0.0230)
        quotes.append(depo)
        helpers.append(ql.DepositRateHelper(ql.QuoteHandle(depo), euribor6m))

        for years, rate in [(1, 0.0235), (2, 0.0248), (5, 0.0270), (10, 0.0280)]:
            q = ql.SimpleQuote(rate)
            quotes.append(q)
            helpers.append(
                ql.SwapRateHelper(
                    ql.QuoteHandle(q),
                    ql.Period(years, ql.Years),
                    ql.TARGET(),
                    ql.Annual,
                    ql.Unadjusted,
                    ql.Thirty360(ql.Thirty360.BondBasis),
                    euribor6m,
                    ql.QuoteHandle(),
                    ql.Period(0, ql.Days),
                    discount_handle,
                )
            )
        curve = ql.PiecewiseLogLinearDiscount(
            0, ql.TARGET(), helpers, ql.Actual360()
        )
        return curve, quotes, helpers

    def test_jacobian_is_analytical_and_square(self):
        """jacobian() covers OIS helpers analytically."""
        curve, _, helpers = self.build_ois_curve()

        analytic = ql.BoolVector()
        jacobian = curve.jacobian(analytic)

        self.assertEqual(jacobian.rows(), len(helpers))
        self.assertEqual(jacobian.columns(), len(curve.data()) - 1)
        self.assertEqual(len(analytic), len(helpers))
        for i, flag in enumerate(analytic):
            self.assertTrue(flag, "row %d was not computed analytically" % i)

        # jacobian() and inverseJacobian() must be mutual inverses
        inverse = curve.inverseJacobian()
        product = jacobian * inverse
        for i in range(product.rows()):
            for j in range(product.columns()):
                self.assertAlmostEqual(
                    product[i][j], 1.0 if i == j else 0.0, delta=1.0e-8,
                    msg="element (%d,%d) of J * inverse(J)" % (i, j),
                )

    def test_inverse_jacobian_matches_bumped_quotes(self):
        """inverseJacobian() reproduces node moves under quote bumps."""
        curve, quotes, _ = self.build_ois_curve()
        inverse = curve.inverseJacobian()

        self.assertEqual(inverse.rows(), len(curve.data()) - 1)
        self.assertEqual(inverse.columns(), len(quotes))

        for k, quote in enumerate(quotes):
            expected = node_derivatives([curve], quote)[0]
            for j, fd in enumerate(expected):
                self.assertAlmostEqual(
                    inverse[j][k], fd, delta=1.0e-5,
                    msg="node %d vs quote %d" % (j, k),
                )

    def test_cross_jacobian(self):
        """Projection helpers respond analytically to the discount nodes."""
        ois, _, _ = self.build_ois_curve()
        projection, _, projection_helpers = self.build_projection_curve(ois)

        graph = ql.CurveJacobianGraph()
        graph.add(ois)
        graph.add(projection)

        analytic = ql.BoolVector()
        cross = graph.crossJacobian(projection, ois, analytic)

        self.assertEqual(cross.rows(), len(projection_helpers))
        self.assertEqual(cross.columns(), len(ois.data()) - 1)
        for i, flag in enumerate(analytic):
            self.assertTrue(flag, "row %d was not computed analytically" % i)

        non_zero = any(
            cross[i][j] != 0.0
            for i in range(cross.rows())
            for j in range(cross.columns())
        )
        self.assertTrue(
            non_zero,
            "the cross-Jacobian of the projection helpers with respect to "
            "the discount nodes is identically zero",
        )

        # the discount curve does not depend on the projection curve
        reverse = graph.crossJacobian(ois, projection)
        for i in range(reverse.rows()):
            for j in range(reverse.columns()):
                self.assertAlmostEqual(reverse[i][j], 0.0, delta=1.0e-12)

    def test_inverse_jacobian_matches_bumped_quotes(self):
        """Inverse-Jacobian blocks match a re-bootstrap."""
        ois, ois_quotes, _ = self.build_ois_curve()
        projection, projection_quotes, _ = self.build_projection_curve(ois)

        graph = ql.CurveJacobianGraph()
        graph.add(ois)
        graph.add(projection)
        self.assertTrue(graph.isComplete())

        partial_graph = ql.CurveJacobianGraph()
        partial_graph.add(projection)
        self.assertFalse(partial_graph.isComplete())
        analytic = ql.BoolVector()
        partial = partial_graph.inverseJacobian(
            projection, projection, analytic
        )
        self.assertTrue(all(analytic))

        strict_graph = ql.CurveJacobianGraph(True)
        strict_graph.add(projection)
        self.assertFalse(strict_graph.isComplete())
        with self.assertRaises(RuntimeError):
            strict_graph.inverseJacobian(projection, projection)
        strict_graph.add(ois)
        self.assertTrue(strict_graph.isComplete())
        strict = strict_graph.inverseJacobian(projection, projection)
        self.assertEqual(partial.rows(), strict.rows())
        self.assertEqual(partial.columns(), strict.columns())
        for i in range(partial.rows()):
            for j in range(partial.columns()):
                self.assertAlmostEqual(partial[i][j], strict[i][j], delta=1.0e-10)

        curves = [ois, projection]
        for name, quotes, wrt in (
            ("discount", ois_quotes, ois),
            ("projection", projection_quotes, projection),
        ):
            sensitivities = []
            for c in curves:
                analytic = ql.BoolVector()
                sensitivities.append(graph.inverseJacobian(c, wrt, analytic))
                self.assertEqual(len(analytic), len(quotes))
                self.assertTrue(all(analytic))
            for k, quote in enumerate(quotes):
                expected = node_derivatives(curves, quote)
                for c, fd_nodes in enumerate(expected):
                    for j, fd in enumerate(fd_nodes):
                        self.assertAlmostEqual(
                            sensitivities[c][j][k], fd, delta=1.0e-5,
                            msg="node %d of curve %d vs quote %d of the %s "
                                "curve" % (j, c, k, name),
                        )

    def test_zero_node_transformations(self):
        """Continuous-zero and stored-node Jacobians are mutual inverses."""
        for curve_type, expected_analytic in (
            (ql.PiecewiseLogLinearDiscount, True),
            (ql.PiecewiseLinearZero, True),
            (ql.PiecewiseLinearForward, False),
        ):
            curve, quotes, _ = self.build_ois_curve(curve_type)
            curve.enableExtrapolation()
            graph = ql.CurveJacobianGraph(True)
            graph.add(curve)
            self.assertTrue(graph.isComplete())

            analytic = ql.BoolVector()
            zero_node = graph.zeroNodeJacobian(curve, analytic)
            node_zero = graph.nodeZeroJacobian(curve)
            nodes = len(curve.data()) - 1
            self.assertEqual(zero_node.rows(), nodes)
            self.assertEqual(zero_node.columns(), nodes)
            self.assertEqual(len(analytic), nodes)
            self.assertTrue(all(flag == expected_analytic for flag in analytic))

            identity = zero_node * node_zero
            for i in range(nodes):
                for j in range(nodes):
                    self.assertAlmostEqual(
                        identity[i][j], 1.0 if i == j else 0.0,
                        delta=1.0e-8,
                    )

            inverse = graph.inverseJacobian(curve, curve)
            dates = list(curve.dates())[1:]
            for k, quote in enumerate(quotes):
                value = quote.value()
                quote.setValue(value + BUMP)
                up = [
                    curve.zeroRate(
                        d, ql.Actual360(), ql.Continuous,
                        ql.NoFrequency, True,
                    ).rate()
                    for d in dates
                ]
                quote.setValue(value - BUMP)
                down = [
                    curve.zeroRate(
                        d, ql.Actual360(), ql.Continuous,
                        ql.NoFrequency, True,
                    ).rate()
                    for d in dates
                ]
                quote.setValue(value)
                curve.data()

                for i in range(nodes):
                    expected = (up[i] - down[i]) / (2.0 * BUMP)
                    calculated = sum(
                        zero_node[i][j] * inverse[j][k]
                        for j in range(nodes)
                    )
                    self.assertAlmostEqual(calculated, expected, delta=2.0e-5)

    def test_par_risk_propagates_across_curves(self):
        """Node risk on one curve becomes par risk on every curve."""
        ois, ois_quotes, _ = self.build_ois_curve()
        projection, projection_quotes, _ = self.build_projection_curve(ois)

        graph = ql.CurveJacobianGraph()
        graph.add(ois)
        graph.add(projection)

        # arbitrary node risk on the projection curve only
        nodes = len(projection.data()) - 1
        node_risk = ql.Array(nodes)
        for j in range(nodes):
            node_risk[j] = 1.0 + 0.5 * j

        analytic = ql.BoolVector()
        risk = graph.parRisk([projection], [node_risk], analytic)
        self.assertEqual(len(risk), 2)
        self.assertEqual(len(analytic), len(ois_quotes) + len(projection_quotes))
        self.assertTrue(all(analytic))
        self.assertEqual(len(risk[0]), len(ois_quotes))
        self.assertEqual(len(risk[1]), len(projection_quotes))

        # the same result, curve by curve
        on_ois = graph.parRisk([projection], [node_risk], ois)
        for k in range(len(on_ois)):
            self.assertAlmostEqual(on_ois[k], risk[0][k], delta=1.0e-12)

        # it must equal transpose(inverseJacobian) * nodeRisk
        for c, curve in enumerate([ois, projection]):
            jacobian = graph.inverseJacobian(projection, curve)
            for k in range(jacobian.columns()):
                expected = sum(
                    jacobian[j][k] * node_risk[j] for j in range(jacobian.rows())
                )
                self.assertAlmostEqual(
                    risk[c][k], expected, delta=1.0e-10,
                    msg="quote %d of curve %d" % (k, c),
                )

        # the discount quotes must carry a non-trivial part of the risk
        self.assertTrue(
            any(abs(risk[0][k]) > 1.0e-10 for k in range(len(risk[0]))),
            "no risk was propagated to the discount curve",
        )

        # Multiple trade-risk contributions on the same curve are summed.
        doubled = graph.parRisk(
            [projection, projection], [node_risk, node_risk]
        )
        for c in range(len(risk)):
            for k in range(len(risk[c])):
                self.assertAlmostEqual(doubled[c][k], 2.0 * risk[c][k],
                                       delta=1.0e-12)

    def test_supported_derived_curve_is_inspected_by_add(self):
        """add() records a derived curve exposing baseCurve()."""
        ois, _, _ = self.build_ois_curve()
        wrapper = ql.SpreadDiscountCurve(
            ql.YieldTermStructureHandle(ois),
            [ois.referenceDate(), ois.maxDate()],
            [1.0, 0.99],
        )
        projection, _, _ = self.build_projection_curve(wrapper)

        graph = ql.CurveJacobianGraph()
        graph.add(ois)
        graph.add(projection)

        open_analytic = ql.BoolVector()
        open_cross = graph.crossJacobian(projection, ois, open_analytic)
        self.assertEqual(open_cross.rows(), len(open_analytic))
        self.assertTrue(open_analytic[0])
        self.assertFalse(any(open_analytic[1:]))

        node_risk = ql.Array(len(projection.data()) - 1)
        for j in range(len(node_risk)):
            node_risk[j] = 1.0 + j
        open_risk = graph.parRisk([projection], [node_risk])

        graph.add(wrapper)
        self.assertEqual(len(graph.curves()), 2)
        analytic = ql.BoolVector()
        cross = graph.crossJacobian(projection, ois, analytic)
        self.assertEqual(cross.rows(), len(analytic))
        # The deposit row has no discount-curve dependence.
        # Swap rows follow the wrapper numerically.
        self.assertTrue(analytic[0])
        self.assertFalse(any(analytic[1:]))

        closed_risk = graph.parRisk([projection], [node_risk])
        for c in range(len(open_risk)):
            for k in range(len(open_risk[c])):
                self.assertAlmostEqual(open_risk[c][k], closed_risk[c][k],
                                       delta=1.0e-10)

    def test_opaque_derived_curve_is_not_supported(self):
        """Derived curves without an underlying-curve getter are rejected."""
        ois, _, _ = self.build_ois_curve()
        wrapper = ql.ZeroSpreadedTermStructure(
            ql.YieldTermStructureHandle(ois),
            ql.QuoteHandle(ql.SimpleQuote(0.0010)),
        )

        graph = ql.CurveJacobianGraph()
        graph.add(ois)
        with self.assertRaisesRegex(RuntimeError, "not supported"):
            graph.add(wrapper)

    def test_piecewise_spread_curve_can_join_graph(self):
        """Spread curves report their base dependency to the graph."""
        base, base_quotes, _ = self.build_ois_curve()
        estr = ql.Estr()
        helpers, spread_quotes = [], []
        for months, rate in [(3, 0.0205), (12, 0.0220), (24, 0.0235)]:
            quote = ql.SimpleQuote(rate)
            spread_quotes.append(quote)
            helpers.append(
                ql.OISRateHelper(
                    2, ql.Period(months, ql.Months),
                    ql.QuoteHandle(quote), estr,
                )
            )
        spread = ql.PiecewiseLogLinearSpreadDiscount(
            ql.YieldTermStructureHandle(base), helpers
        )
        spread.enableExtrapolation()

        graph = ql.CurveJacobianGraph()
        graph.add(base)
        graph.add(spread)

        analytic = ql.BoolVector()
        cross = graph.crossJacobian(spread, base, analytic)
        self.assertEqual(cross.rows(), len(helpers))
        self.assertTrue(any(abs(cross[i][j]) > 1.0e-12
                            for i in range(cross.rows())
                            for j in range(cross.columns())))
        self.assertTrue(all(analytic))

        # The local spread-node scale includes the base discount factor.
        local_analytic = ql.BoolVector()
        inverse = spread.inverseJacobian(local_analytic)
        self.assertTrue(all(local_analytic))
        for k, quote in enumerate(spread_quotes):
            expected = node_derivatives([spread], quote)[0]
            for j, value in enumerate(expected):
                self.assertAlmostEqual(inverse[j][k], value, delta=1.0e-5)

        # The graph-level derivative includes rebootstrap of the spread curve
        # when a quote on its base curve moves.
        propagated = graph.inverseJacobian(spread, base)
        for k, quote in enumerate(base_quotes):
            expected = node_derivatives([spread], quote)[0]
            for j, value in enumerate(expected):
                self.assertAlmostEqual(propagated[j][k], value, delta=1.0e-5)

        # A further curve using the spread curve must preserve the complete
        # base -> spread -> projection dependency chain.
        projection, _, _ = self.build_projection_curve(spread)
        graph.add(projection)
        layered_analytic = ql.BoolVector()
        layered = graph.inverseJacobian(projection, base, layered_analytic)
        self.assertTrue(all(layered_analytic))
        for k, quote in enumerate(base_quotes):
            expected = node_derivatives([projection], quote)[0]
            for j, value in enumerate(expected):
                self.assertAlmostEqual(layered[j][k], value, delta=2.0e-5)

    def test_graph_rejects_unregistered_and_unsupported_curves(self):
        ois, _, _ = self.build_ois_curve()
        projection, _, _ = self.build_projection_curve(ois)

        graph = ql.CurveJacobianGraph()
        graph.add(ois)

        self.assertEqual(len(graph.curves()), 1)
        with self.assertRaises(RuntimeError):
            graph.crossJacobian(projection, ois)

        flat = ql.FlatForward(
            self.today, ql.QuoteHandle(ql.SimpleQuote(0.02)), ql.Actual360()
        )
        with self.assertRaises(RuntimeError):
            graph.add(flat)

    def test_helper_quote_sensitivities(self):
        ois, _, helpers = self.build_ois_curve()
        ois.discount(1.0)  # trigger the bootstrap

        helper = helpers[-1]
        self.assertTrue(helper.hasAnalyticQuoteSensitivities())
        self.assertTrue(helper.hasCompleteQuoteSensitivities(ois))

        own = helper.impliedQuoteSensitivities()
        self.assertTrue(len(own) > 0)
        by_curve = helper.impliedQuoteSensitivities(ois)
        self.assertEqual(len(own), len(by_curve))
        for (time, value), (date, other) in zip(own, by_curve):
            self.assertAlmostEqual(time, ois.timeFromReference(date), delta=1.0e-12)
            self.assertAlmostEqual(value, other, delta=1.0e-12)


class InterpolationNodeWeightsTest(unittest.TestCase):
    """nodeWeights() is the analytical primitive the Jacobians build on."""

    def test_linear_interpolation_weights(self):
        x = ql.Array(4)
        y = ql.Array(4)
        for j, (xj, yj) in enumerate([(0.0, 1.0), (1.0, 2.0),
                                      (2.0, 1.5), (3.0, 3.0)]):
            x[j], y[j] = xj, yj
        f = ql.LinearInterpolation(x, y)

        # halfway between nodes 1 and 2: equal weights, no other node
        weights = dict(f.nodeWeights(1.5))
        self.assertEqual(sorted(weights), [1, 2])
        self.assertAlmostEqual(weights[1], 0.5, delta=1.0e-12)
        self.assertAlmostEqual(weights[2], 0.5, delta=1.0e-12)

        # the weights must reproduce a bump of a node value
        bump = 1.0e-7
        for j in range(len(y)):
            bumped = ql.Array(len(y))
            for k in range(len(y)):
                bumped[k] = y[k]
            bumped[j] = y[j] + bump
            up = ql.LinearInterpolation(x, bumped)(1.5)
            expected = (up - f(1.5)) / bump
            self.assertAlmostEqual(
                weights.get(j, 0.0), expected, delta=1.0e-6,
                msg="node %d" % j,
            )

    def test_unsupported_interpolation_returns_nothing(self):
        x = ql.Array(4)
        y = ql.Array(4)
        for j in range(4):
            x[j], y[j] = float(j), 1.0 + 0.1 * j
        # node weights are not implemented for every interpolation
        self.assertEqual(len(ql.LagrangeInterpolation(x, y).nodeWeights(1.5)), 0)


class CoDependentCurves:
    """3M and 6M curves depending on each other, bootstrapped jointly.

    3M curve: FRAs, then 3s6s basis swaps; 6M curve: 3s6s basis swaps,
    then 6M swaps.  Mirrors the C++ test-suite fixture.
    """

    ACCURACY = 1.0e-11

    def __init__(self, today, analytic_jacobian):
        self.int3m = ql.RelinkableYieldTermStructureHandle()
        self.int6m = ql.RelinkableYieldTermStructureHandle()
        euribor3m = ql.Euribor3M(self.int3m)
        euribor6m = ql.Euribor6M(self.int6m)
        discount = ql.YieldTermStructureHandle(
            ql.FlatForward(today, 0.02, ql.Actual360())
        )

        self.quotes3m, self.quotes6m = [], []
        helpers3m, helpers6m = [], []

        for i in (1, 4, 7, 10, 13, 16, 19):
            q = ql.SimpleQuote(0.024 + 0.0001 * i)
            self.quotes3m.append(q)
            helpers3m.append(
                ql.FraRateHelper(
                    ql.QuoteHandle(q), i, i + 3, euribor3m.fixingDays(),
                    euribor3m.fixingCalendar(),
                    euribor3m.businessDayConvention(),
                    euribor3m.endOfMonth(), euribor3m.dayCounter(),
                    ql.Pillar.LastRelevantDate,
                )
            )
        for y in (2, 3, 4):
            q = ql.SimpleQuote(0.0015 + 0.0001 * y)
            self.quotes3m.append(q)
            helpers3m.append(
                ql.IborIborBasisSwapRateHelper(
                    ql.QuoteHandle(q), ql.Period(y, ql.Years),
                    euribor3m.fixingDays(), euribor3m.fixingCalendar(),
                    euribor3m.businessDayConvention(), euribor3m.endOfMonth(),
                    euribor3m, euribor6m, discount, True,
                )
            )
        for m in (6, 12, 18):
            q = ql.SimpleQuote(0.0012 + 0.0001 * (m // 6))
            self.quotes6m.append(q)
            helpers6m.append(
                ql.IborIborBasisSwapRateHelper(
                    ql.QuoteHandle(q), ql.Period(m, ql.Months),
                    euribor3m.fixingDays(), euribor3m.fixingCalendar(),
                    euribor3m.businessDayConvention(), euribor3m.endOfMonth(),
                    euribor3m, euribor6m, discount, False,
                )
            )
        for y in (5, 7, 10):
            q = ql.SimpleQuote(0.026 + 0.0002 * y)
            self.quotes6m.append(q)
            helpers6m.append(
                ql.SwapRateHelper(
                    ql.QuoteHandle(q), ql.Period(y, ql.Years),
                    euribor6m.fixingCalendar(), ql.Annual, ql.Following,
                    ql.Thirty360(ql.Thirty360.BondBasis), euribor6m,
                    ql.QuoteHandle(), ql.Period(0, ql.Days), discount,
                )
            )

        self.curve3m = ql.GlobalPiecewiseLogLinearDiscount(
            today, helpers3m, ql.Actual360(),
            ql.GlobalBootstrap(self.ACCURACY),
        )
        self.curve6m = ql.GlobalPiecewiseLogLinearDiscount(
            today, helpers6m, ql.Actual360(),
            ql.GlobalBootstrap(self.ACCURACY),
        )

        # keep every piece alive: the handles outlive the curves, and the
        # MultiCurve owns the joint bootstrap
        self.discount = discount
        self.index3m, self.index6m = euribor3m, euribor6m
        self.helpers3m, self.helpers6m = helpers3m, helpers6m
        self.multi_curve = ql.MultiCurve(self.ACCURACY, analytic_jacobian)
        self.external3m = self.multi_curve.addBootstrappedCurve(
            self.int3m, self.curve3m
        )
        self.external6m = self.multi_curve.addBootstrappedCurve(
            self.int6m, self.curve6m
        )

    @property
    def curves(self):
        return [self.curve3m, self.curve6m]

    @property
    def quotes(self):
        return self.quotes3m + self.quotes6m


class AnalyticJacobianBootstrapTest(unittest.TestCase):
    def setUp(self):
        self.today = ql.Date(23, ql.October, 2025)
        ql.Settings.instance().evaluationDate = self.today

    def build_ois_helpers(self):
        estr = ql.Estr()
        helpers = []
        for months, rate in [(3, 0.0195), (12, 0.0210), (24, 0.0225),
                             (60, 0.0250), (120, 0.0260)]:
            helpers.append(
                ql.OISRateHelper(2, ql.Period(months, ql.Months),
                                 ql.QuoteHandle(ql.SimpleQuote(rate)), estr)
            )
        return helpers

    def test_global_bootstrap_with_analytic_jacobian(self):
        """The analytical Jacobian converges to the same curve."""
        accuracy = 1.0e-12

        default = ql.GlobalPiecewiseLogLinearDiscount(
            0, ql.TARGET(), self.build_ois_helpers(), ql.Actual360(),
            ql.GlobalBootstrap(accuracy),
        )
        analytic = ql.GlobalPiecewiseLogLinearDiscount(
            0, ql.TARGET(), self.build_ois_helpers(), ql.Actual360(),
            ql.GlobalBootstrap(accuracy, None, None, True),
        )

        expected = list(default.data())
        obtained = list(analytic.data())
        self.assertEqual(len(expected), len(obtained))
        for j, (e, o) in enumerate(zip(expected, obtained)):
            self.assertAlmostEqual(o, e, delta=1.0e-8, msg="node %d" % j)

    def test_supplied_optimizer_must_consume_analytic_jacobian(self):
        """The analytic flag cannot silently use an FD optimizer."""
        accuracy = 1.0e-12
        ignored = ql.GlobalPiecewiseLogLinearDiscount(
            0, ql.TARGET(), self.build_ois_helpers(), ql.Actual360(),
            ql.GlobalBootstrap(
                accuracy, ql.LevenbergMarquardt(), None, True
            ),
        )
        with self.assertRaises(RuntimeError):
            list(ignored.data())

        optimizer = ql.LevenbergMarquardt(
            accuracy, accuracy, accuracy, True
        )
        self.assertTrue(optimizer.usesCostFunctionJacobian())
        analytic = ql.GlobalPiecewiseLogLinearDiscount(
            0, ql.TARGET(), self.build_ois_helpers(), ql.Actual360(),
            ql.GlobalBootstrap(accuracy, optimizer, None, True),
        )
        self.assertTrue(len(analytic.data()) > 1)

    def test_multi_curve_with_analytic_jacobian(self):
        """MultiCurve accepts and honours the analyticJacobian flag."""
        default = CoDependentCurves(self.today, False)
        expected = [list(c.data()) for c in default.curves]

        analytic = CoDependentCurves(self.today, True)
        obtained = [list(c.data()) for c in analytic.curves]

        for c, (a, b) in enumerate(zip(expected, obtained)):
            self.assertEqual(len(a), len(b))
            for j, (x, y) in enumerate(zip(a, b)):
                self.assertAlmostEqual(
                    x, y, delta=1.0e-8, msg="node %d of curve %d" % (j, c)
                )

    def test_group_inverse_jacobian_spans_all_quotes(self):
        """inverseJacobian() of a group member covers the whole group."""
        pair = CoDependentCurves(self.today, False)
        curves, quotes = pair.curves, pair.quotes

        # own quotes first, then the other member's
        ordered = {
            id(pair.curve3m): pair.quotes3m + pair.quotes6m,
            id(pair.curve6m): pair.quotes6m + pair.quotes3m,
        }

        for c, curve in enumerate(curves):
            inverse = curve.inverseJacobian()
            self.assertEqual(inverse.rows(), len(curve.data()) - 1)
            self.assertEqual(inverse.columns(), len(quotes))

            for k, quote in enumerate(ordered[id(curve)]):
                expected = node_derivatives([curve], quote)[0]
                for j, fd in enumerate(expected):
                    self.assertAlmostEqual(
                        inverse[j][k], fd, delta=1.0e-4,
                        msg="node %d of curve %d vs group quote %d"
                            % (j, c, k),
                    )


if __name__ == "__main__":
    unittest.main()
