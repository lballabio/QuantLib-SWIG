/*
 Copyright (C) 2018 Klaus Spanderen


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

package examples;

import org.quantlib.{Array => QArray, _}

object Swing {
    def main(args: Array[String]) : Unit = {
        val todaysDate = new Date(30, Month.September, 2018)
        val dc = new Actual365Fixed
        
        val riskFreeRate = new FlatForward(todaysDate, 0.0, dc)
        val dividendYield = new FlatForward(todaysDate, 0.0, dc)
        
        val underlying = new SimpleQuote(30.0)
        val vol = new BlackConstantVol(todaysDate, new TARGET, 0.20, dc)
        
        val exerciseDates = new DateVector
        (0 until 31).foreach(i => exerciseDates.add(
            new Date(1, Month.January, 2019).add(new Period(i, TimeUnit.Days))))
            
        val swingOption = new VanillaSwingOption(
            new VanillaForwardPayoff(Option.Type.Call, underlying.value), 
            new SwingExercise(exerciseDates), 0, exerciseDates.size)
            
        val bsProcess = new BlackScholesMertonProcess(
            new QuoteHandle(underlying),
            new YieldTermStructureHandle(dividendYield),
            new YieldTermStructureHandle(riskFreeRate),
            new BlackVolTermStructureHandle(vol))
            
        swingOption.setPricingEngine(new FdSimpleBSSwingEngine(bsProcess))

        println(f"Black Scholes Price: ${swingOption.NPV}%2.4f")
        
                
        val x0 = 0d
        val x1 = 0d

        val beta = 4d
        val eta = 4d
        val jumpIntensity = 1d
        val speed = 1d
        val volatility = 0.1
        
        val curveShape = new DoublePairVector()
        (0 until exerciseDates.size.toInt).foreach{i =>  
            val t = dc.yearFraction(todaysDate, exerciseDates.get(i))
            val gs = math.log(underlying.value) -                                  
                volatility*volatility/(4*speed)*(1-math.exp(-2*speed*t)) -
                jumpIntensity/beta*math.log((eta-math.exp(-beta*t))/(eta-1.0))
                
            curveShape.add(new DoublePair(t, gs))
        }
        
        val ouProcess = new ExtendedOrnsteinUhlenbeckProcess(
            speed, volatility, x0, new UnaryFunctionDelegate() {
                override def value(x: Double): Double = x0
            })
            
        val jProcess = new ExtOUWithJumpsProcess(
            ouProcess, x1, beta, jumpIntensity, eta)
            
        swingOption.setPricingEngine(new FdSimpleExtOUJumpSwingEngine(
            jProcess, riskFreeRate, 25, 25, 200, curveShape))

        println(f"Kluge Model Price  : ${swingOption.NPV}%2.4f")    
    }
}
