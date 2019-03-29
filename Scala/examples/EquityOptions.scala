/*
 Copyright (C) 2011 Klaus Spanderen


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

import scala.actors.Actor
import org.quantlib.{Array => QArray, _}


/**
 * EquityOption Test app - simple multithreading Scala version of 
 * QuantLib/Examples/EquityOption
 * to illustrate use of Quantlib through supplied SWIG interfaces.
 *
 * You need to run this using the Java Jar file and JNI library
  */


class VanillaPricingService(payoff: PlainVanillaPayoff, 
                            exercise: Exercise) extends Actor {
    start()

    def act() {
        react {
            case (engine: PricingEngine) => {

                // copy instrument data to ensure thread safe execution

                val t = exercise exerciseType match {
                    case Exercise.Type.European => new EuropeanExercise(
                                                    exercise.dates().get(0))
                    case Exercise.Type.Bermudan => new BermudanExercise(
                                                    exercise.dates());
                    case Exercise.Type.American =>
                        new AmericanExercise(
                            exercise.dates get 0, exercise.dates().get(
                            (exercise.dates.size() - 1).toInt))
                }
                val instrument = new VanillaOption(
                    new PlainVanillaPayoff(payoff.optionType, payoff.strike),
                                           exercise)
                instrument.setPricingEngine(engine)
                reply(instrument.NPV())
            }
        }
    }
}

object SimpleFactory {
    val optionType = Option.Type.Put
    val strike     = 40.0
    val underlying = 36.0
    val riskFreeRate = 0.06
    val dividendYield = 0.00
    val volatility = 0.2

    val calendar = new TARGET()
    val dayCounter = new Actual365Fixed()
    val settlementDate = new Date(17, Month.May, 1998)

    def bsProcess() : BlackScholesMertonProcess = {
        val flatVolatility = new BlackVolTermStructureHandle(
                    new BlackConstantVol(settlementDate, calendar, 
                                         volatility, dayCounter))
        new BlackScholesMertonProcess(spot, divYield,
                                      rTS, flatVolatility)
    }

    def hestonProcess() : HestonProcess = {
        new HestonProcess(rTS, divYield, spot, volatility*volatility,
                          1.0, volatility*volatility, 0.0001, 0.0)
    }

    def batesProcess() : BatesProcess = {
        new BatesProcess(rTS, divYield, spot, volatility*volatility,
                         1.0, volatility*volatility, 0.0001, 0.0,
                         1e-14, 1e-14, 1e-14)
    }
    private def rTS : YieldTermStructureHandle = {
        new YieldTermStructureHandle(
                    new FlatForward(settlementDate, riskFreeRate, dayCounter))
    }
    private def divYield : YieldTermStructureHandle = {
        new YieldTermStructureHandle(
                    new FlatForward(settlementDate, dividendYield, dayCounter))
    }
    private def spot : QuoteHandle = {
        new QuoteHandle(new SimpleQuote(underlying))
    }
}

object EquityOptions {
    def main(args: Array[String]) : Unit = {

        try {
            System.loadLibrary("QuantLibJNI");
        } 
        catch {
            case ex: UnsatisfiedLinkError => {
                  println("please check your LD_LIBRARY_PATH variable")
                throw ex
            }
        }

        val beginTime  = System.currentTimeMillis()
        
        val optionType = Option.Type.Put
        val strike     = 40.0
        val todaysDate = new Date(15, Month.May, 1998)
        val settlementDate = SimpleFactory.settlementDate
        Settings.instance setEvaluationDate todaysDate

        val maturity = new Date(17, Month.May, 1999)
        val dayCounter = new Actual365Fixed()
        val calendar = new TARGET()

        // define European, Bermudan, and American exercises
        val exerciseDates = new DateVector()
        (1 to 4).foreach(i => exerciseDates add settlementDate.
                                     add(new Period(3*i, TimeUnit.Months)))
        
        val europeanExercise = new EuropeanExercise(maturity)
        val bermudanExercise = new BermudanExercise(exerciseDates)
        val americanExercise = new AmericanExercise(settlementDate, maturity)

        val payoff = new PlainVanillaPayoff(optionType, strike)

        // Black-Scholes for European
        val analyticEuropeanNpv =
                 new VanillaPricingService(payoff, europeanExercise) !! 
                        new AnalyticEuropeanEngine(SimpleFactory.bsProcess())

        val hestonModel = new HestonModel(SimpleFactory.hestonProcess())
             
        // Heston for European                            
        val analyticHestonNpv = 
            new VanillaPricingService(payoff, europeanExercise) !!
                new AnalyticHestonEngine(new HestonModel(
                    SimpleFactory.hestonProcess()))
                
        val fdEuropeanHestonNpv =
            new VanillaPricingService(payoff, europeanExercise) !!
                new FdHestonVanillaEngine(new HestonModel(
                    SimpleFactory.hestonProcess()), 50, 150) 

        val fdAmericanHestonNpv =
            new VanillaPricingService(payoff, americanExercise) !!
                new FdHestonVanillaEngine(new HestonModel(
                    SimpleFactory.hestonProcess()), 100, 150) 

        val fdBermudanHestonNpv =
            new VanillaPricingService(payoff, bermudanExercise) !!
                new FdHestonVanillaEngine(new HestonModel(
                    SimpleFactory.hestonProcess()), 100, 150) 

        val batesModel = new BatesModel(SimpleFactory.batesProcess())
        
        val cosHestonNpv = 
            new VanillaPricingService(payoff, europeanExercise) !!
                new COSHestonEngine(new HestonModel(
                                    SimpleFactory.hestonProcess()))

        // Bates for European
        val analyticBatesNpv =
            new VanillaPricingService(payoff, europeanExercise) !!
                new BatesEngine(new BatesModel(
                	SimpleFactory.batesProcess()))

        val fdEuropeanBatesNpv =
            new VanillaPricingService(payoff, europeanExercise) !!
                new FdBatesVanillaEngine(new BatesModel(
                    SimpleFactory.batesProcess()), 50, 150)

        val fdAmericanBatesNpv =
            new VanillaPricingService(payoff, americanExercise) !!
                new FdBatesVanillaEngine(new BatesModel(
                    SimpleFactory.batesProcess()))

        val fdBermudanBatesNpv =
            new VanillaPricingService(payoff, bermudanExercise) !!
                new FdBatesVanillaEngine(new BatesModel(
                    SimpleFactory.batesProcess()))

        // Barone-Adesi and Whaley approximation for American
        val baroneAdesiWhaleyNpv = 
            new VanillaPricingService(payoff, americanExercise) !!
                new BaroneAdesiWhaleyEngine(SimpleFactory.bsProcess())

        // Bjerksund and Stensland approximation for American
        val bjerksundStenslandNpv = 
            new VanillaPricingService(payoff, americanExercise) !!
                new BjerksundStenslandEngine(SimpleFactory.bsProcess())

        // Integral
        val integralNpv = new VanillaPricingService(payoff,europeanExercise) !!
                                new IntegralEngine(SimpleFactory.bsProcess())

        // Finite Difference
        var timeSteps : Int = 801;
        val fdEuropeanNpv =   
            new VanillaPricingService(payoff, europeanExercise) !! 
                new FDEuropeanEngine(SimpleFactory.bsProcess(),     
                                     timeSteps, timeSteps-1)
        val fdBermudanNpv = 
            new VanillaPricingService(payoff, bermudanExercise) !!
                  new FDBermudanEngine(SimpleFactory.bsProcess(), 
                                     timeSteps, timeSteps-1)
        val fdAmericanNpv = 
            new VanillaPricingService(payoff, americanExercise) !!
                  new FDAmericanEngine(SimpleFactory.bsProcess(), 
                                     timeSteps, timeSteps-1)

        // Binomial method      
        val jarrowRuddEuropeanNpv = 
            new VanillaPricingService(payoff, europeanExercise) !!
                       new BinomialJRVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val jarrowRuddBermudanNpv = 
            new VanillaPricingService(payoff, bermudanExercise) !! 
                    new BinomialJRVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val jarrowRuddAmericanNpv = 
            new VanillaPricingService(payoff, americanExercise) !!
                    new BinomialJRVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val coxRossRubinsteinEuropeanNpv = 
            new VanillaPricingService(payoff, europeanExercise) !!
                   new BinomialCRRVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val coxRossRubinsteinBermudanNpv = 
            new VanillaPricingService(payoff, bermudanExercise) !!
                   new BinomialCRRVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val coxRossRubinsteinAmericanNpv = 
            new VanillaPricingService(payoff, americanExercise) !!
                   new BinomialCRRVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val additiveEqpEuropeanNpv =                                          
            new VanillaPricingService(payoff, europeanExercise) !!
             new BinomialEQPVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val additiveEqpBermudanNpv = 
            new VanillaPricingService(payoff, bermudanExercise) !!
                 new BinomialEQPVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val additiveEqpAmericanNpv = 
            new VanillaPricingService(payoff, americanExercise) !!
                 new BinomialEQPVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val trigeirgisEuropeanNpv = 
            new VanillaPricingService(payoff, europeanExercise) !!
                        new BinomialTrigeorgisVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val trigeirgisBermudanNpv = 
            new VanillaPricingService(payoff, bermudanExercise) !!
                          new BinomialTrigeorgisVanillaEngine(SimpleFactory.bsProcess(), timeSteps) 
        val trigeirgisAmericanNpv =      
            new VanillaPricingService(payoff, americanExercise) !!
                          new BinomialTrigeorgisVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val tianEuropeanNpv = 
            new VanillaPricingService(payoff, europeanExercise) !!
                          new BinomialTianVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val tianBermudanNpv = 
            new VanillaPricingService(payoff, bermudanExercise) !!
                          new BinomialTianVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val tianAmericanNpv = 
            new VanillaPricingService(payoff, americanExercise) !!
                          new BinomialTianVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val leisenReimerEuropeanNpv = 
            new VanillaPricingService(payoff, europeanExercise) !!
                          new BinomialTianVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val leisenReimerBermudanNpv = 
            new VanillaPricingService(payoff, bermudanExercise) !!
                          new BinomialLRVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val leisenReimerAmericanNpv = 
            new VanillaPricingService(payoff, americanExercise) !!
                          new BinomialLRVanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val joshiEuropeanNpv = 
            new VanillaPricingService(payoff, europeanExercise) !!
                          new BinomialJ4VanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val joshiBermudanNpv = 
            new VanillaPricingService(payoff, bermudanExercise) !!
                          new BinomialJ4VanillaEngine(SimpleFactory.bsProcess(), timeSteps)
        val joshiAmericanNpv = 
            new VanillaPricingService(payoff, americanExercise) !!
                          new BinomialJ4VanillaEngine(SimpleFactory.bsProcess(), timeSteps)

        // Monte-Carlo methods
        timeSteps = 1;
        val americanTimeSteps = 25
        val mcSeed = 42;
        val nSamples = 32768; // 2^15
        val maxSamples = 1048576; // 2^20

        val pseudoMcEuropeanNpv = 
            new VanillaPricingService(payoff, europeanExercise) !!
                new MCPREuropeanEngine(SimpleFactory.bsProcess(),
                                       timeSteps,
                                       QuantLib.nullInt(),
                                       true, false,
                                       nSamples, 0.02, maxSamples, mcSeed)

        val pseudoMcAmericanNpv = 
            new VanillaPricingService(payoff, americanExercise) !!
                new MCPRAmericanEngine(SimpleFactory.bsProcess(),
                                       americanTimeSteps,
                                       QuantLib.nullInt(),
                                       true, false, 
                                       nSamples, 0.02, maxSamples, mcSeed)

        val quasiMcEuropeanNpv = 
            new VanillaPricingService(payoff, europeanExercise) !!
                new MCLDEuropeanEngine(SimpleFactory.bsProcess(),
                                       timeSteps,
                                       QuantLib.nullInt(),
                                       false, false,
                                       nSamples, 0.02, maxSamples, mcSeed)

        val quasiMcAmericanNpv = 
            new VanillaPricingService(payoff, americanExercise) !!
                new MCLDAmericanEngine(SimpleFactory.bsProcess(),
                                       americanTimeSteps,
                                       QuantLib.nullInt(),
                                       true, false, 
                                       nSamples, 0.02, maxSamples, mcSeed)


        // write column headings
        printf("\n%-35s %-14s %-14s %-14s\n" + "="*76+ "\n", 
               "Method", "European", "Bermudan", "American")

        val fmt = "%34s %13.9f %13.9f %13.9f\n";        
        printf(fmt, "Black-Scholes", analyticEuropeanNpv(), 
                                     Double.NaN, Double.NaN)
        printf(fmt, "Heston Semi-Analytic", analyticHestonNpv(), 
                                            Double.NaN, Double.NaN)
        printf(fmt, "Heston Finite-Difference", 
        	fdEuropeanHestonNpv(), fdBermudanHestonNpv(), fdAmericanHestonNpv())
        printf(fmt, "COS Heston Method", cosHestonNpv(), 
                                            Double.NaN, Double.NaN)
        printf(fmt, "Bates Semi-Analytic", analyticBatesNpv(), 
                                            Double.NaN, Double.NaN)
        printf(fmt, "Bates Finite-Difference", 
        	fdEuropeanBatesNpv(), fdBermudanBatesNpv(), fdAmericanBatesNpv())
        printf(fmt, "Barone-Adesi/Whaley", Double.NaN, Double.NaN,
                                           baroneAdesiWhaleyNpv());
        printf(fmt, "Bjerksund/Stensland", Double.NaN, Double.NaN,
                                           bjerksundStenslandNpv())
        printf(fmt, "Integral", integralNpv(),
                                           Double.NaN, Double.NaN)
        printf(fmt, "Finite differences", fdEuropeanNpv(),fdBermudanNpv(),
                                          fdAmericanNpv())
        printf(fmt, "Binomial Jarrow-Rudd",jarrowRuddEuropeanNpv(),
                            jarrowRuddBermudanNpv(),jarrowRuddAmericanNpv()) 
        printf(fmt, "Binomial Cox-Ross-Rubinstein",
                                        coxRossRubinsteinEuropeanNpv(),
                                        coxRossRubinsteinBermudanNpv(),    
                                        coxRossRubinsteinAmericanNpv())
        printf(fmt, "Additive equiprobabilities",additiveEqpEuropeanNpv(),
                        additiveEqpBermudanNpv(),additiveEqpAmericanNpv())
        printf(fmt, "Binomial Trigeorgis", trigeirgisEuropeanNpv(),
                            trigeirgisBermudanNpv(),trigeirgisAmericanNpv())
        printf(fmt, "Binomial Tian", tianEuropeanNpv(),
                                        tianBermudanNpv(),tianAmericanNpv())
        printf(fmt, "Binomial Leisen-Reimer", leisenReimerEuropeanNpv(),
                        leisenReimerBermudanNpv(), leisenReimerAmericanNpv())
        printf(fmt, "Binomial Joshi", joshiEuropeanNpv(), 
                                      joshiBermudanNpv(), joshiAmericanNpv())
        printf(fmt, "MC (crude)", pseudoMcEuropeanNpv(), 
                                  Double.NaN, pseudoMcAmericanNpv())
        printf(fmt, "MC (Sobol)", quasiMcEuropeanNpv(), 
                                  Double.NaN, quasiMcAmericanNpv())

        val msecs = (System.currentTimeMillis()-beginTime)
        println("Run completed in "+msecs+" ms.")
    }
}

