"""
 Copyright (C) 2020 Klaus Spanderen

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


class FdmTest(unittest.TestCase):
    def setUp(self):
        self.todaysDate = ql.Date(15, ql.May, 2019)
        ql.Settings.instance().evaluationDate = self.todaysDate

    def tearDown(self):
        ql.Settings.instance().evaluationDate = ql.Date()


    def test1dMesher(self):
        """Testing one dimensional mesher"""

        m = ql.Concentrating1dMesher(0, 1, 10)
        self.assertEqual(m.size(), 10)
        for i in range(0,10):
            self.assertAlmostEqual(m.location(i), i/9.0, 14)

        m = ql.Concentrating1dMesher(0, 1, 10,
            [ql.Concentrating1dMesherPoint(0.75, 0.01,False),
             ql.Concentrating1dMesherPoint(0.5, 0.01, True)])

        self.assertEqual(m.size(), 10)
        self.assertAlmostEqual(m.location(0), 0.0, 14)
        self.assertAlmostEqual(m.location(9), 1.0, 14)

        p = list(x for x in m.locations() if ql.close_enough(x, 0.5))
        self.assertEqual(len(p), 1)
        p = list(x for x in m.locations() if ql.close_enough(x, 0.75))
        self.assertEqual(len(p), 0)

        m = ql.Predefined1dMesher([0,2,4])
        self.assertEqual(m.size(), 3)
        self.assertEqual(m.location(0), 0)
        self.assertEqual(m.location(1), 2)
        self.assertEqual(m.location(2), 4)

    def testFdmLinearOpIterator(self):
        """Testing iterators for linear operators"""

        dim = [2,2,3]
        pos = [0,0,0]
        idx = 0
        opIter = ql.FdmLinearOpIterator(dim, pos, idx)

        self.assertEqual(opIter.index(), 0)

        opIter.increment()
        self.assertEqual(opIter.index(), 1)
        self.assertEqual(opIter.coordinates(), (1, 0, 0))
        opIter.increment()
        self.assertEqual(opIter.coordinates(), (0, 1, 0))

        opIter2 = ql.FdmLinearOpIterator(dim, pos, idx)
        self.assertEqual(opIter.notEqual(opIter2), True)
        self.assertEqual(opIter.notEqual(opIter), False)


    def testFdmLinearOpLayout(self):
        """Testing memory layout for linear operators"""

        dim = [2,2,3]

        m = ql.FdmLinearOpLayout(dim)

        self.assertEqual(m.size(), 2*2*3)
        self.assertEqual(m.dim(), (2, 2, 3))
        self.assertEqual(m.spacing(), (1, 2, 4))
        self.assertEqual(m.index((0,1,2)), 10)
        self.assertEqual(m.neighbourhood(m.begin(), 0, 1), 1)
        self.assertEqual(m.neighbourhood(m.begin(), 2, 2), 8)
        self.assertEqual(m.neighbourhood(m.begin(), 0, 1, 2, 2), 9)

        n = m.iter_neighbourhood(m.begin(), 0, 1)
        opIter = m.begin()
        opIter.increment()

        self.assertEqual(opIter.notEqual(n), False)

    def testFdmMesherComposite(self):
        """Testing mesher composites"""

        m1 = ql.Concentrating1dMesher(0, 1, 2)
        m2 = ql.Uniform1dMesher(0, 2, 3)

        m = ql.FdmMesherComposite(m1, m2)
        self.assertEqual(len(m.getFdm1dMeshers()), 2)

        locations = m.locations(0)
        self.assertEqual(len(locations), 6)

        self.assertEqual(list(map(lambda x: int(x+0.5), locations)), [0, 1, 0, 1, 0, 1])

        locations = m.locations(1)
        self.assertEqual(list(map(lambda x: int(x+0.5), locations)), [0, 0, 1, 1, 2, 2])

    def testFdmLinearOpComposite(self):
        """Testing linear operator composites"""

        class Foo:
            t1 = 0.0
            t2 = 0.0

            @classmethod
            def size(self):
                return 42

            def setTime(self, t1, t2):
                self.t1 = t1
                self.t2 = t2

            @classmethod
            def apply(self, r):
                return 2*r

            @classmethod
            def apply_mixed(self, r):
                return 3*r

            @classmethod
            def apply_direction(self, direction , r):
                return direction*r

            @classmethod
            def solve_splitting(self, direction , r, s):
                return direction*s*r

            @classmethod
            def preconditioner(self, r, s):
                return s*r


        foo = Foo()

        c = ql.FdmLinearOpCompositeProxy(foo)

        self.assertEqual(c.size(), foo.size())

        c.setTime(1.0, 2.0)
        self.assertAlmostEqual(foo.t1, 1.0, 14)
        self.assertAlmostEqual(foo.t2, 2.0, 14)

        r = ql.Array([1,2,3,4])
        self.assertEqual(list(c.apply(r)), list(2*r))
        self.assertEqual(list(c.apply_mixed(r)), list(3*r))
        self.assertEqual(list(c.apply_direction(7, r)), list(7*r))

        s = list(c.solve_splitting(7, r, 0.5))
        self.assertEqual(len(s), len(r))
        for i, x in enumerate(s):
            self.assertAlmostEqual(x, 3.5*r[i], 14)

        self.assertEqual(list(c.preconditioner(r, 4)), list(4*r))

        class Bar:
            @classmethod
            def apply(self, r):
                return 1

            def apply_mixed(self, r):
                pass

        with self.assertRaises(RuntimeError):
            ql.FdmLinearOpCompositeProxy(Bar()).apply(r)

        with self.assertRaises(RuntimeError):
            ql.FdmLinearOpCompositeProxy(Bar()).apply_mixed(r)


    def testFdmBlackScholesOp(self):
        """Testing linear Black-Scholes operator"""

        todaysDate = ql.Date(1, ql.January, 2020)
        ql.Settings.instance().evaluationDate = todaysDate
        dc = ql.Actual365Fixed()

        settlementDate = todaysDate + 2
        riskFreeRate = ql.FlatForward(settlementDate, 0.05, dc)

        exercise = ql.EuropeanExercise(ql.Date(27, ql.December, 2020))
        maturity = dc.yearFraction(todaysDate, exercise.lastDate())

        strike = 110.0
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, strike)

        underlying = ql.SimpleQuote(100.0)
        volatility = ql.BlackConstantVol(settlementDate, ql.TARGET(), 0.10, dc)
        dividendYield = ql.FlatForward(settlementDate, 0.05, dc)

        process = ql.BlackScholesMertonProcess(
            ql.QuoteHandle(underlying),
            ql.YieldTermStructureHandle(dividendYield),
            ql.YieldTermStructureHandle(riskFreeRate),
            ql.BlackVolTermStructureHandle(volatility)
        )

        mesher = ql.FdmMesherComposite(
            ql.FdmBlackScholesMesher(10, process, maturity, strike))

        op = ql.FdmBlackScholesOp(mesher, process, strike)
        self.assertEqual(op.size(), 1)

        op.setTime(0, 0.1)

        c = list(map(lambda x: payoff(math.exp(x)), mesher.locations(0)))
        p = op.apply(c)

        e = [ 0.0, 0.0, 0.0, 0.0, 0.0,
              3.18353, 0.755402, -1.30583, -2.19881, -4.0271 ]

        for i, x in enumerate(e):
            self.assertAlmostEqual(x, p[i], 5)


    def testFdmFirstOrderOperator(self):
        """Testing first order operator"""

        mesher = ql.Uniform1dMesher(0.0, math.pi, 1000)

        op = ql.FirstDerivativeOp(0, ql.FdmMesherComposite(mesher))

        l = mesher.locations()

        x = list(map(math.sin, l))

        y = op.apply(x)

        for u, v in zip(l, y):
            self.assertAlmostEqual(v, math.cos(u), 4)


    def testFdmSecondOrderOperator(self):
        """Testing second order operator"""

        mesher = ql.Uniform1dMesher(0.0, math.pi, 1000)

        op = ql.SecondDerivativeOp(0, ql.FdmMesherComposite(mesher))

        x = list(map(math.sin, mesher.locations()))

        y = op.apply(x)

        for u, v in zip(x, y):
            self.assertAlmostEqual(v, -u, 4)

    def testFdmBoundaryCondition(self):
        """Testing Dirichlet Boundary conditions"""

        m = ql.FdmMesherComposite(
            ql.Uniform1dMesher(0.0, 1.0, 5))

        b = ql.FdmDirichletBoundary(
            m, math.pi, 0, ql.FdmBoundaryCondition.Upper)

        x = ql.Array(len(m.locations(0)), 0.0)

        b.applyAfterApplying(x)

        self.assertEqual(list(x), [0,0,0,0, math.pi])

        s = ql.FdmBoundaryConditionSet()
        s.push_back(b)

        self.assertEqual(len(s), 1)

    def testFdmStepConditionCallBack(self):
        """Testing step condition call back function"""

        class Foo:
            @classmethod
            def applyTo(self, a, t):
                for i in range(5):
                    a[i] = t+1.0

        m = ql.FdmStepConditionProxy(Foo())

        x = ql.Array(5)

        m.applyTo(x, 2.0)

        self.assertEqual(len(x), 5)
        self.assertEqual(list(x), [3.0, 3.0, 3.0, 3.0, 3.0])

    def testFdmInnerValueCalculatorCallBack(self):
        """Testing inner value call back function"""

        class Foo:
            @classmethod
            def innerValue(self, opIter, t):
                return opIter.index() + t

            @classmethod
            def avgInnerValue(self, opIter, t):
                return opIter.index() + 2*t

        m = ql.FdmInnerValueCalculatorProxy(Foo())

        dim = [2,2,3]
        pos = [0,0,0]

        opIter = ql.FdmLinearOpIterator(dim, pos, 0)

        while (opIter.index() < 2*2*3):
            idx = opIter.index()

            self.assertEqual(m.innerValue(opIter, 2.0), idx + 2.0)
            self.assertEqual(m.avgInnerValue(opIter, 2.0), idx + 4.0)

            opIter.increment()


    def testFdmLogInnerValueCalculator(self):
        """Testing log inner value calculator"""

        m = ql.FdmMesherComposite(
            ql.Uniform1dMesher(math.log(50), math.log(150), 11))

        p = ql.PlainVanillaPayoff(ql.Option.Call, 100)

        v = ql.FdmLogInnerValue(p, m, 0)

        opIter = m.layout().begin()
        while opIter.notEqual(m.layout().end()):
            x = math.exp(m.location(opIter, 0));
            self.assertAlmostEqual(p(x), v.innerValue(opIter, 1.0), 14)
            opIter.increment()


    def testAmericanOptionPricing(self):
        """Testing Black-Scholes and Heston American Option pricing"""

        xSteps = 100
        tSteps = 25
        dampingSteps = 0

        todaysDate = ql.Date(15, ql.January, 2020)
        ql.Settings.instance().evaluationDate = todaysDate

        dc = ql.Actual365Fixed()

        riskFreeRate = ql.YieldTermStructureHandle(
            ql.FlatForward(todaysDate, 0.06, dc))
        dividendYield = ql.YieldTermStructureHandle(
            ql.FlatForward(todaysDate, 0.02, dc))

        strike = 110.0
        payoff = ql.PlainVanillaPayoff(ql.Option.Put, strike)

        maturityDate = todaysDate + ql.Period(1, ql.Years)
        maturity = dc.yearFraction(todaysDate, maturityDate)

        exercise = ql.AmericanExercise(todaysDate, maturityDate)

        spot = ql.QuoteHandle(ql.SimpleQuote(100.0))
        volatility = ql.BlackConstantVol(todaysDate, ql.TARGET(), 0.20, dc)

        process = ql.BlackScholesMertonProcess(
            spot, dividendYield, riskFreeRate,
            ql.BlackVolTermStructureHandle(volatility)
        )

        option = ql.VanillaOption(payoff, exercise)
        option.setPricingEngine(ql.FdBlackScholesVanillaEngine.make(
            process, xGrid = xSteps, tGrid = tSteps,
            dampingSteps = dampingSteps)
        )

        expected = option.NPV()

        equityMesher = ql.FdmBlackScholesMesher(
            xSteps, process, maturity,
            strike, cPoint = (strike, 0.1)
        )

        mesher = ql.FdmMesherComposite(equityMesher)

        op = ql.FdmBlackScholesOp(mesher, process, strike)

        innerValueCalculator = ql.FdmLogInnerValue(payoff, mesher, 0)

        x = []
        rhs = []
        layout = mesher.layout()
        opIter = layout.begin()
        while (opIter.notEqual(layout.end())):
            x.append(mesher.location(opIter, 0))
            rhs.append(innerValueCalculator.avgInnerValue(opIter, maturity))
            opIter.increment()

        rhs = ql.Array(rhs)

        bcSet = ql.FdmBoundaryConditionSet()
        stepCondition = ql.FdmStepConditionComposite.vanillaComposite(
            ql.DividendSchedule(), exercise, mesher,
            innerValueCalculator, todaysDate, dc
        )

        # only to test an Operator defined in python
        class OperatorProxy:
            def __init__(self, op):
                self.op = op

            def size(self):
                return self.op.size()

            def setTime(self, t1, t2):
                return self.op.setTime(t1, t2)

            def apply(self, r):
                return self.op.apply(r)

            def apply_direction(self, i, r):
                return self.op.apply_direction(i, r)

            def solve_splitting(self, i, r, s):
                return self.op.solve_splitting(i, r, s)


        proxyOp = ql.FdmLinearOpCompositeProxy(OperatorProxy(op))

        solver = ql.FdmBackwardSolver(
            proxyOp, bcSet, stepCondition, ql.FdmSchemeDesc.Douglas()
        )

        solver.rollback(rhs, maturity, 0.0, tSteps, dampingSteps)

        spline = ql.CubicNaturalSpline(x, rhs);

        logS = math.log(spot.value())

        calculated = spline(logS)

        self.assertAlmostEqual(calculated, expected, 1)

        solverDesc = ql.FdmSolverDesc(
            mesher, bcSet, stepCondition, innerValueCalculator,
            maturity, tSteps, dampingSteps)

        calculated = ql.Fdm1DimSolver(
            solverDesc, ql.FdmSchemeDesc.Douglas(), op).interpolateAt(logS)

        self.assertAlmostEqual(calculated, expected, 2)

        v0 = 0.4*0.4
        kappa = 1.0
        theta = v0
        sigma = 1e-4
        rho = 0.0

        hestonProcess = ql.HestonProcess(
            riskFreeRate, dividendYield,
            spot, v0, kappa, theta, sigma, rho)

        leverageFct = ql.LocalVolSurface(
            ql.BlackVolTermStructureHandle(
                ql.BlackConstantVol(todaysDate, ql.TARGET(), 0.50, dc)),
            riskFreeRate,
            dividendYield,
            spot.value()
        )

        vSteps = 3

        vMesher = ql.FdmHestonLocalVolatilityVarianceMesher(
            vSteps, hestonProcess, leverageFct, maturity)

        avgVolaEstimate = vMesher.volaEstimate()

        self.assertAlmostEqual(avgVolaEstimate, 0.2, 5)

        mesher = ql.FdmMesherComposite(equityMesher, vMesher)

        innerValueCalculator = ql.FdmLogInnerValue(payoff, mesher, 0)

        stepCondition = ql.FdmStepConditionComposite.vanillaComposite(
            ql.DividendSchedule(), exercise, mesher,
            innerValueCalculator, todaysDate, dc
        )

        solverDesc = ql.FdmSolverDesc(
            mesher, bcSet, stepCondition, innerValueCalculator,
            maturity, tSteps, dampingSteps)

        calculated = ql.FdmHestonSolver(
            hestonProcess, solverDesc, leverageFct = leverageFct).valueAt(
                spot.value(), 0.16)

        self.assertAlmostEqual(calculated, expected, 1)


    def testBSMRNDCalculator(self):
        """Testing Black-Scholes risk neutral density calculator"""

        dc = ql.Actual365Fixed()
        todaysDate = ql.Date(15, ql.January, 2020)

        r   = 0.0
        q   = 0.0
        vol = 0.2
        s0  = 100

        process = ql.BlackScholesMertonProcess(
            ql.QuoteHandle(ql.SimpleQuote(s0)),
            ql.YieldTermStructureHandle(
                ql.FlatForward(todaysDate, q, dc)),
            ql.YieldTermStructureHandle(
                ql.FlatForward(todaysDate, r, dc)),
            ql.BlackVolTermStructureHandle(
                ql.BlackConstantVol(todaysDate, ql.TARGET(), vol, dc))
        )

        rnd = ql.BSMRNDCalculator(process)

        t = 1.2
        x = math.log(80.0)

        mu = math.log(s0) + (r-q-0.5*vol*vol)*t

        calculated = rnd.pdf(x, t)

        stdev = vol * math.sqrt(t)

        expected = (1.0/(math.sqrt(2*math.pi)*stdev) *
            math.exp( -0.5*math.pow((x-mu)/stdev, 2.0) ))

        self.assertAlmostEqual(calculated, expected, 8)


    def testOrnsteinUhlenbeckVsBachelier(self):
        """Testing Fdm Ornstein-Uhlenbeck pricing"""

        todaysDate = ql.Date(15, ql.January, 2020)
        ql.Settings.instance().evaluationDate = todaysDate

        dc = ql.Actual365Fixed()

        rTS = ql.FlatForward(todaysDate, 0.06, dc)

        strike = 110.0
        payoff = ql.PlainVanillaPayoff(ql.Option.Put, strike)

        maturityDate = todaysDate + ql.Period(2, ql.Years)

        exercise = ql.EuropeanExercise(maturityDate)

        option = ql.VanillaOption(payoff, exercise)

        x0 = 100
        sigma = 20.0
        speed = 5

        pdeEngine = ql.FdOrnsteinUhlenbeckVanillaEngine(
            ql.OrnsteinUhlenbeckProcess(speed, sigma, x0, x0), rTS, 50
        )

        option.setPricingEngine(pdeEngine)
        calculated = option.NPV()

        stdev = math.sqrt(sigma*sigma/(2*speed))

        expected = ql.bachelierBlackFormula(
            ql.Option.Put,
            strike, x0, stdev,
            rTS.discount(maturityDate)
        )

        self.assertAlmostEqual(calculated, expected, 2)


    def testSparseLinearMatrixSolver(self):
        """Testing sparse linear matrix solver"""

        A = ql.Matrix([
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.5],
            [1.0, 0.5, 1.0]
        ])

        b = ql.Array([ 1.0, 0.2, 0.5 ])

        expected = ql.inverse(A)*b

        def foo(x):
            return A*x

        calculated = ql.BiCGstab(
            ql.MatrixMultiplicationProxy(foo), 100, 1e-6).solve(b)

        for i in range(3):
            self.assertAlmostEqual(expected[i], calculated[i], 4)

        calculated = ql.GMRES(
            ql.MatrixMultiplicationProxy(foo), 100, 1e-6).solve(b)

        for i in range(3):
            self.assertAlmostEqual(expected[i], calculated[i], 4)

        def preconditioner(x):
            return ql.inverse(A)*x

        calculated = ql.BiCGstab(
            ql.MatrixMultiplicationProxy(foo), 100, 1e-6,
            ql.MatrixMultiplicationProxy(preconditioner)).solve(b)

        for i in range(3):
            self.assertAlmostEqual(expected[i], calculated[i], 4)

    def testGlued1dMesher(self):
        """Testing sparse linear matrix solver"""

        m1 = ql.Uniform1dMesher(0, 2, 3)
        m2 = ql.Uniform1dMesher(2, 4, 3)

        m3 = ql.Glued1dMesher(m1, m2)

        self.assertEqual(m3.locations(), (0,1,2,3,4))

    def testFdmZeroInnerValue(self):
        """Testing FdmZeroInnerValue"""
        opIter = ql.FdmLinearOpIterator([1], [0], 0)

        self.assertEqual(ql.FdmZeroInnerValue().innerValue(opIter, 1.0), 0.0)


if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(FdmTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)
