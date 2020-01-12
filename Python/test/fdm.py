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
        ql.Settings.instance().setEvaluationDate(ql.Date())
        
    def testDisposableArray(self):
        """Testing disposable array"""
        a = ql.Array([1,2,3,4])
        b = ql.DisposableArray(a)
        self.assertEqual([1,2,3,4], list(b))
        self.assertEqual(len(a), 0)
        
        c = ql.Array(b)
        self.assertEqual([1,2,3,4], list(c))
        self.assertEqual(len(b), 0)

        with self.assertRaises(TypeError) as context:
            a + 1
        
        
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
        
        dim = ql.UnsignedIntVector([2,2,3])
        pos = ql.UnsignedIntVector([0,0,0])
        idx = 0
        iter = ql.FdmLinearOpIterator(dim, pos, idx)
        
        self.assertEqual(iter.index(), 0)
        
        iter.increment()
        self.assertEqual(iter.index(), 1)        
        self.assertEqual(iter.coordinates(), (1, 0, 0))
        iter.increment()
        self.assertEqual(iter.coordinates(), (0, 1, 0))

        iter2 = ql.FdmLinearOpIterator(dim, pos, idx)
        self.assertEqual(iter.notEqual(iter2), True)
        self.assertEqual(iter.notEqual(iter), False)


    def testFdmLinearOpLayout(self):
        """Testing memory layout for linear operators"""
        
        dim = ql.UnsignedIntVector([2,2,3])
        
        m = ql.FdmLinearOpLayout(dim)
        
        self.assertEqual(m.size(), 2*2*3)
        self.assertEqual(m.dim(), (2, 2, 3))
        self.assertEqual(m.spacing(), (1, 2, 4))
        self.assertEqual(m.index((0,1,2)), 10)
        self.assertEqual(m.neighbourhood(m.begin(), 0, 1), 1)
        self.assertEqual(m.neighbourhood(m.begin(), 2, 2), 8)
        self.assertEqual(m.neighbourhood(m.begin(), 0, 1, 2, 2), 9)
        
        n = m.iter_neighbourhood(m.begin(), 0, 1)        
        iter = m.begin()
        iter.increment()
        
        self.assertEqual(iter.notEqual(n), False)

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
            def size(self):
                return 42

            def setTime(self, t1, t2):
                self.t1 = t1
                self.t2 = t2
                
            def apply(self, r):
                return 2*r

            def apply_mixed(self, r):            
                return 3*r

            def apply_direction(self, direction , r):
                return direction*r

            def solve_splitting(self, direction , r, s):
                return direction*s*r

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
            def apply(self, r):
                return 1
            
            def apply_mixed(self, r):
                pass
           
        with self.assertRaises(RuntimeError) as context:
            ql.FdmLinearOpCompositeProxy(Bar()).apply(r)

        with self.assertRaises(RuntimeError) as context:            
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
        
        x = list(map(lambda x: math.sin(x), l))
        
        y = op.apply(x)
        
        for u, v in zip(l, y):
            self.assertAlmostEqual(v, math.cos(u), 4)
            
            
    def testFdmSecondOrderOperator(self):
        """Testing second order operator"""
                        
        mesher = ql.Uniform1dMesher(0.0, math.pi, 1000)
            
        op = ql.SecondDerivativeOp(0, ql.FdmMesherComposite(mesher))
        
        x = list(map(lambda x: math.sin(x), mesher.locations()))
        
        y = op.apply(x)
        
        for u, v in zip(x, y):
            self.assertAlmostEqual(v, -u, 4)
            
    def testFdmBoundaryCondition(self):
        """Testing Diriclet Boundary conditions"""

        m = ql.FdmMesherComposite(
            ql.Uniform1dMesher(0.0, 1.0, 5))
        
        b = ql.FdmDirichletBoundary(
            m, math.pi, 0, ql.BoundaryConditionFdmLinearOp.Upper)
        
        x = ql.Array(len(m.locations(0)), 0.0)
        
        b.applyAfterApplying(x)
        
        self.assertEqual(list(x), [0,0,0,0, math.pi])
        
        s = ql.FdmBoundaryConditionSet()
        s.push_back(b)
        
        self.assertEqual(len(s), 1)
        
        ql.BoundaryConditionFdmLinearOp.NoSide # None is replaced by _None
        
    def testFdmStepConditionCallBack(self):
        """Testing step condition call back function"""

        class Foo:
            def applyTo(self, a, t):                
                for i,x in enumerate(a):
                    a[i] = t+1.0
            
        m = ql.FdmStepConditionProxy(Foo())
        
        x = ql.Array(5)
        
        m.applyTo(x, 2.0)
        
        self.assertEqual(len(x), 5)
        self.assertEqual(list(x), [3.0, 3.0, 3.0, 3.0, 3.0])       
        
    def testFdmInnerValueCalculatorCallBack(self):
        """Testing inner value call back function"""

        class Foo:
            def innerValue(self, iter, t):                
                return iter.index() + t
              
            def avgInnerValue(self, iter, t):                
                return iter.index() + 2*t

        m = ql.FdmInnerValueCalculatorProxy(Foo())
        
        dim = ql.UnsignedIntVector([2,2,3])
        pos = ql.UnsignedIntVector([0,0,0])
        
        iter = ql.FdmLinearOpIterator(dim, pos, 0)
        
        for i in range(2*2*3):
            idx = iter.index()
            
            self.assertEqual(m.innerValue(iter, 2.0), idx + 2.0)
            self.assertEqual(m.avgInnerValue(iter, 2.0), idx + 4.0)
            
            iter.increment()
    
    
    def testFdmLogInnerValueCalculator(self):
        """Testing log inner value calculator"""
        
        m = ql.FdmMesherComposite(
            ql.Uniform1dMesher(math.log(50), math.log(150), 11))
        
        p = ql.PlainVanillaPayoff(ql.Option.Call, 100)
        
        v = ql.FdmLogInnerValue(p, m, 0)
        
        iter = m.layout().begin()
        for i in range(m.layout().size()):
            x = math.exp(m.location(iter, 0));
            self.assertAlmostEqual(p(x), v.innerValue(iter, 1.0), 14)
            iter.increment()

            
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
        iter = layout.begin()
        while (iter.notEqual(layout.end())):
            x.append(mesher.location(iter, 0))
            rhs.append(innerValueCalculator.avgInnerValue(iter, maturity))
            iter.increment()
            
        rhs = ql.Array(rhs)    
        
        bcSet = ql.FdmBoundaryConditionSet()
        stepCondition = ql.FdmStepConditionComposite.vanillaComposite(
            ql.DividendSchedule(), exercise, mesher,
            innerValueCalculator, todaysDate, dc
        )
            
        solver = ql.FdmBackwardSolver(
            op, bcSet, stepCondition, ql.FdmSchemeDesc.Douglas()
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
        
        
                    
if __name__ == "__main__":
    print("testing QuantLib " + ql.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(FdmTest, "test"))
    unittest.TextTestRunner(verbosity=2).run(suite)