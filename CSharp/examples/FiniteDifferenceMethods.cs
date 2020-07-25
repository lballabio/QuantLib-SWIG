/*
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
*/

using System;

using QuantLib;
using static QuantLib.NQuantLibc;

namespace FiniteDifferenceMethods
{
    class FdmDemo
    {
        // test class for the operator delegate and proxy
        public class FdmBSDelegate : FdmLinearOpCompositeDelegate 
        {
            private readonly FdmLinearOpComposite op;

            public FdmBSDelegate(FdmLinearOpComposite op)
            {
                this.op = op;
            }

            override public uint size() 
            {
                return op.size();
            }

            override public void setTime(double t1, double t2)
            {
                op.setTime(t1, t2);
            }

            override public QlArray apply(QlArray r)
            {
                return op.apply(r);
            }

            override public QlArray apply_direction(uint i, QlArray r)
            {
                return op.apply_direction(i, r);
            }

            override public QlArray solve_splitting(uint i, QlArray r, double s)
            {
                return op.solve_splitting(i, r, s);
            }
        };

        static void Main(string[] args) 
        {
            const int xSteps = 100;
            const int tSteps = 25;
            const int dampingSteps = 0;

            Date today = new Date(15, Month.January, 2020);
            Settings.instance().setEvaluationDate(today);

            DayCounter dc = new Actual365Fixed();

            YieldTermStructureHandle rTS = new YieldTermStructureHandle(
                new FlatForward(today, 0.06, dc));
            YieldTermStructureHandle qTS = new YieldTermStructureHandle(
                new FlatForward(today, 0.02, dc));

            const double strike = 110.0;
            StrikedTypePayoff payoff = new PlainVanillaPayoff(Option.Type.Put, strike);

            Date maturityDate = today.Add(new Period(1, TimeUnit.Years));
            double maturity = dc.yearFraction(today, maturityDate);

            Exercise exercise = new AmericanExercise(today, maturityDate);

            Instrument vanillaOption = new VanillaOption(payoff, exercise);

            QuoteHandle spot = new QuoteHandle(new SimpleQuote(100.0));
            BlackVolTermStructureHandle volatility = new BlackVolTermStructureHandle(
                new BlackConstantVol(today, new TARGET(), 0.20, dc));

            BlackScholesMertonProcess process = 
                new BlackScholesMertonProcess(spot, qTS, rTS, volatility);

            vanillaOption.setPricingEngine(new FdBlackScholesVanillaEngine(
                process, tSteps, xSteps, dampingSteps));

            double expected = vanillaOption.NPV();

            // build an PDE engine from scratch
            Fdm1dMesher equityMesher = new FdmBlackScholesMesher(
                xSteps, process, maturity, strike, 
                nullDouble(), nullDouble(), 0.0001, 1.5, 
                new DoublePair(strike, 0.1));

            FdmMesherComposite mesher = new FdmMesherComposite(equityMesher);

            FdmLinearOpComposite op = new FdmBlackScholesOp(mesher, process, strike);

            FdmInnerValueCalculator calc = new FdmLogInnerValue(payoff, mesher, 0);

            QlArray x = new QlArray(equityMesher.size());
            QlArray rhs = new QlArray(equityMesher.size());

            FdmLinearOpIterator iter = mesher.layout().begin();            
            for (uint i = 0; i < rhs.size(); ++i, iter.increment()) 
            {
                x.set(i, mesher.location(iter, 0));
                rhs.set(i, calc.avgInnerValue(iter, maturity));
            }

            FdmBoundaryConditionSet bcSet = new FdmBoundaryConditionSet();

            FdmStepConditionComposite stepCondition = 
                FdmStepConditionComposite.vanillaComposite(
                    new DividendSchedule(), exercise, mesher, calc, today, dc);


            FdmLinearOpComposite proxyOp = new FdmLinearOpCompositeProxy(
                new FdmBSDelegate(op));

            FdmBackwardSolver solver = new FdmBackwardSolver(
                proxyOp, bcSet, stepCondition, FdmSchemeDesc.Douglas());

            solver.rollback(rhs, maturity, 0.0, tSteps, dampingSteps);

            double logS = Math.Log(spot.value());

            double calculated = new CubicNaturalSpline(x, rhs).call(logS);

            Console.WriteLine("Homebrew PDE engine        : {0:0.0000}", calculated);
            Console.WriteLine("FdBlackScholesVanillaEngine: {0:0.0000}", expected);
        }
    }
}