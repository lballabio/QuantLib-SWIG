/*
  Copyright (C) 2008 Allen Kuo
  Copyright (C) 2014 Wondersys Srl
  Copyright (C) 2017 BN Algorithms Ltd

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

import org.quantlib.CallableFixedRateBond;
import org.quantlib.Date;
import org.quantlib.Period;
import org.quantlib.Settings;
import org.quantlib.DayCounter;
import org.quantlib.ActualActual;
import org.quantlib.Compounding;
import org.quantlib.Frequency;
import org.quantlib.FlatForward;
import org.quantlib.QuoteHandle;
import org.quantlib.Quote;
import org.quantlib.SimpleQuote;
import org.quantlib.YieldTermStructure;
import org.quantlib.YieldTermStructureHandle;
import org.quantlib.InterestRate;
import org.quantlib.Month;
import org.quantlib.Date;
import org.quantlib.CallabilitySchedule;
import org.quantlib.Schedule;
import org.quantlib.NullCalendar;
import org.quantlib.UnitedStates;
import org.quantlib.Calendar;
import org.quantlib.Callability;
import org.quantlib.CallabilityPrice;
import org.quantlib.TimeUnit;
import org.quantlib.BusinessDayConvention;
import org.quantlib.DateGeneration;
import org.quantlib.ShortRateModel;
import org.quantlib.HullWhite;
import org.quantlib.PricingEngine;
import org.quantlib.TreeCallableFixedRateBondEngine;
import org.quantlib.DoubleVector;

public class CallableBonds {

    public static YieldTermStructure flatRate(Date today,
                                              Quote forward,
                                              DayCounter dc,
                                              Compounding compounding,
                                              Frequency frequency)
    {
        return new FlatForward(today,
                               new QuoteHandle(forward),
                               dc,
                               compounding,
                               frequency);
    }

    public static void priceCallable(double sigma,
                                     YieldTermStructureHandle termStructure,
                                     CallableFixedRateBond callableBond,
                                     Frequency frequency,
                                     DayCounter bondDayCounter)
    {
        long maxIterations = 1000;
        double accuracy = 1e-8;
        long gridIntervals = 40;
        double reversionParameter = .03;

        ShortRateModel hw0 =
            new HullWhite(termStructure,
                          reversionParameter,
                          sigma);

        PricingEngine engine0 =
            new TreeCallableFixedRateBondEngine(hw0,
                                                gridIntervals);

        callableBond.setPricingEngine(engine0);

        System.out.printf("sigma/vol (%%) =  %6.2f \n",
                          100.*sigma);

        System.out.printf("QuantLib price/yld (%%)  %10.2f  /  %10.2f \n",
                          callableBond.cleanPrice(),
                          100.0 * callableBond.yield(bondDayCounter,
                                                     Compounding.Compounded,
                                                     frequency,
                                                     accuracy,
                                                     maxIterations));

	System.out.printf("QuantLib OAS from model clean price (bp)  %10.2f \n",
			  10000.0 * callableBond.OAS(callableBond.cleanPrice(),
						     termStructure,
						     bondDayCounter,
						     Compounding.Compounded,
						     frequency));


	double cp=callableBond.cleanPriceOAS(10*1e-4,
					   termStructure,
					   bondDayCounter,
					   Compounding.Compounded,
					   frequency);
	System.out.printf("QuantLib spreaded clean price with 10bp OAS %f \n",
			  cp);

	System.out.printf("QuantLib OAS from spreaded clean price (bp)  %10.2f \n",
			  10000.0 * callableBond.OAS(cp,
						     termStructure,
						     bondDayCounter,
						     Compounding.Compounded,
						     frequency));

	System.out.printf("QuantLib effectiveDuration / convexity for 10bp OAS %f / %f \n",
			  callableBond.effectiveDuration(10*1e-4,
							 termStructure,
							 bondDayCounter,
							 Compounding.Compounded,
							 frequency),
			  callableBond.effectiveConvexity(10*1e-4,
							  termStructure,
							  bondDayCounter,
							  Compounding.Compounded,
							  frequency));
    }


    public static void main(String[] args) throws Exception {
        System.out.println("Callable Bonds example:");

        Date today = new Date(16, Month.October, 2007);
        Settings.instance().setEvaluationDate(today);

        double bbCurveRate = 0.055;
        DayCounter bbDayCounter =
            new ActualActual(ActualActual.Convention.Bond);

        InterestRate bbIR =
            new InterestRate(bbCurveRate,
                             bbDayCounter,
                             Compounding.Compounded,
                             Frequency.Semiannual);

        YieldTermStructureHandle termStructure =
            new YieldTermStructureHandle(flatRate(today,
                                                  new SimpleQuote(bbIR.rate()),
                                                  bbIR.dayCounter(),
                                                  bbIR.compounding(),
                                                  bbIR.frequency()));


        // set up the call schedule
        CallabilitySchedule callSchedule =
            new CallabilitySchedule();
        double callPrice = 100.;
        long numberOfCallDates = 24;
        Date callDate =
            new Date(15, Month.September, 2006);

        Calendar nullCalendar = new NullCalendar();

        for (long i=0; i< numberOfCallDates; ++i) {

            CallabilityPrice myPrice=
                new CallabilityPrice(callPrice,
                                     CallabilityPrice.Type.Clean);
            callSchedule.add(new Callability(myPrice,
                                             Callability.Call,
                                             callDate));
            callDate = nullCalendar.advance(callDate, 3, TimeUnit.Months);
        }

        // set up the callable bond
        Date dated =
            new Date(16, Month.September, 2004);
        Date issue = dated;
        Date maturity = new Date(15, Month.September, 2012);
        int settlementDays = 3;  // Bloomberg OAS1 settle is Oct 19, 2007
        Calendar bondCalendar = new
            UnitedStates(UnitedStates.Market.GovernmentBond);
        double coupon = .0465;
        Frequency frequency =  Frequency.Quarterly;
        double redemption = 100.0;
        double faceAmount = 100.0;

        /* The 30/360 day counter Bloomberg uses for this bond cannot
           reproduce the US Bond/ISMA (constant) cashflows used in PFC1.
           Therefore use ActAct(Bond)
        */
        DayCounter bondDayCounter =
            new ActualActual(ActualActual.Convention.Bond);

        // PFC1 shows no indication dates are being adjusted
        // for weekends/holidays for vanilla bonds
        BusinessDayConvention accrualConvention = BusinessDayConvention.Unadjusted;
        BusinessDayConvention paymentConvention = BusinessDayConvention.Unadjusted;

        Schedule sch =
            new Schedule(dated,
                         maturity,
                         new Period(frequency),
                         bondCalendar,
                         accrualConvention,
                         accrualConvention,
                         DateGeneration.Rule.Backward,
                         false);

        DoubleVector couponsVector = new DoubleVector();
        couponsVector.add(coupon);

        CallableFixedRateBond callableBond =
            new CallableFixedRateBond(settlementDays,
                                      faceAmount,
                                      sch,
                                      couponsVector,
                                      bondDayCounter,
                                      paymentConvention,
                                      redemption,
                                      issue,
                                      callSchedule);


        priceCallable(1e-20,
                      termStructure,
                      callableBond,
                      frequency,
                      bondDayCounter);

        System.out.printf("Bloomberg price/yld (%%) 96.50 / 5.47 \n");

        priceCallable(0.01,
                      termStructure,
                      callableBond,
                      frequency,
                      bondDayCounter);

        System.out.printf("Bloomberg price/yld (%%) 95.68 / 5.66 \n");

        priceCallable(0.03,
                      termStructure,
                      callableBond,
                      frequency,
                      bondDayCounter);

        System.out.printf("Bloomberg price/yld (%%) 92.34 / 6.49 \n");

        priceCallable(0.06,
                      termStructure,
                      callableBond,
                      frequency,
                      bondDayCounter);

        System.out.printf("Bloomberg price/yld (%%) 87.16 / 7.83 \n");

        priceCallable(0.12,
                      termStructure,
                      callableBond,
                      frequency,
                      bondDayCounter);

        System.out.printf("Bloomberg price/yld (%%) 77.31 / 10.65 \n");


        System.out.println("Done");
    }

}
