/*

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

import java.util.stream.IntStream;

import org.quantlib.Settings;

import org.quantlib.CallableFixedRateBond;
import org.quantlib.CallabilitySchedule;
import org.quantlib.CallabilityPrice;
import org.quantlib.Callability;
import org.quantlib.TreeCallableFixedRateBondEngine;
import org.quantlib.HullWhite;

import org.quantlib.YieldTermStructure;
import org.quantlib.YieldTermStructureHandle;
import org.quantlib.RelinkableYieldTermStructureHandle;

import org.quantlib.PiecewiseFlatForward;
import org.quantlib.TimeUnit;
import org.quantlib.RateHelperVector;
import org.quantlib.Calendar;
import org.quantlib.Compounding;
import org.quantlib.UnitedStates;
import org.quantlib.NullCalendar;
import org.quantlib.BusinessDayConvention;
import org.quantlib.DateGeneration;
import org.quantlib.Period;
import org.quantlib.Schedule;
import org.quantlib.DoubleVector;
    
import org.quantlib.ZeroCouponBond;
import org.quantlib.FlatForward;
import org.quantlib.Actual360;
import org.quantlib.Thirty360;
import org.quantlib.ActualActual;
import org.quantlib.Month;
import org.quantlib.Frequency;
import org.quantlib.Date;
import org.quantlib.DiscountingBondEngine;
import org.quantlib.BondHelper;
import org.quantlib.QuoteHandle;
import org.quantlib.SimpleQuote;

/** Examples focusing on the OAS and related functions of callable
 * bonds.
 */
public class CallableBondsOAS {

    static public final Date today = new Date(28, Month.March, 2017);
	


    static public YieldTermStructure mkYC()
    {
	Calendar cal = new UnitedStates(UnitedStates.Market.GovernmentBond);
	    
	final int np = 12;
	double ycrate[] = { 0.772,                
			    0.917,              
			    1.011,               
			    1.302,
			    1.549,               
			    1.822,               
			    1.964,               
			    2.242,               
			    2.418,
			    2.763,              
			    3.025,              
			    3.025};
	int ycInt[]    =  { 3, 6, 1, 2, 3, 4, 5, 7, 10, 20, 30, 50 };
	TimeUnit ycUnit[]   = { TimeUnit.Months, TimeUnit.Months,
				TimeUnit.Years, TimeUnit.Years, TimeUnit.Years, TimeUnit.Years, TimeUnit.Years,
				TimeUnit.Years, TimeUnit.Years, TimeUnit.Years, TimeUnit.Years, TimeUnit.Years};

	RateHelperVector helpers = new  RateHelperVector();
	for (int i=0; i<np; ++i) {
	    final double yield = ycrate[i] * 0.01;
	    FlatForward r = new FlatForward(1, cal, yield,
					    new Actual360(),
					    Compounding.Compounded,
					    Frequency.Annual);
	    ZeroCouponBond z = new ZeroCouponBond( 1,
						   cal,
						   100.0,
						   cal.advance(today, ycInt[i], ycUnit[i]));
	    DiscountingBondEngine pe = new DiscountingBondEngine(new YieldTermStructureHandle(r));

	    z.setPricingEngine(pe);
	    double cleanprice = z.cleanPrice();
	    helpers.add(new BondHelper( new QuoteHandle( new SimpleQuote(cleanprice)), z));
        }

	return new PiecewiseFlatForward(today,
					helpers,
					// Check this convention
					new ActualActual( ActualActual.Convention.Bond));
    }

    static public CallableFixedRateBond mkUST1()
    {
	double coupon = 0.07125;
	Date issue = new Date(15, Month.February, 1993);
	Date matur = new Date(15, Month.February, 2023);

	DoubleVector coupons = new DoubleVector();
        coupons.add(coupon);

	return new CallableFixedRateBond(1,
					 100.0,
					 new Schedule(issue,
						      matur,
						      new Period(Frequency.Semiannual),
						      new UnitedStates(UnitedStates.Market.GovernmentBond),
						      BusinessDayConvention.Unadjusted,
						      BusinessDayConvention.Unadjusted,						      
						      DateGeneration.Rule.Backward,
						      false),
					 coupons,
					 new ActualActual(ActualActual.Convention.Bond),
					 BusinessDayConvention.ModifiedFollowing,
					 100.0,
					 issue,
					 new CallabilitySchedule());
    }

    static public CallableFixedRateBond mkBLRDG()
    {
	double coupon = 0.0115;
	Date issue = new Date(29, Month.June, 2012);
	Date matur = new Date(29, Month.June, 2017);

	DoubleVector coupons = new DoubleVector();
        coupons.add(coupon);

        CallabilitySchedule callSchedule =
	    new CallabilitySchedule();

        double callPrice = 100.;
        long numberOfCallDates = 60;
        Date callDate =issue;
	Calendar nullCalendar = new NullCalendar();	

        for (long i=0; i< numberOfCallDates; ++i) {

            CallabilityPrice myPrice=
		new CallabilityPrice(callPrice,
				      CallabilityPrice.Type.Clean);
            callSchedule.add(new Callability(myPrice,
					     Callability.Call,
					     callDate));
            callDate = nullCalendar.advance(callDate, 1, TimeUnit.Months);
        }	

	return new CallableFixedRateBond(1,
					 100.0,
					 new Schedule(issue,
						      matur,
						      new Period(Frequency.Monthly),
						      new UnitedStates(UnitedStates.Market.GovernmentBond),
						      BusinessDayConvention.Unadjusted,
						      BusinessDayConvention.Unadjusted,						      
						      DateGeneration.Rule.Backward,
						      false),
					 coupons,
					 new ActualActual(ActualActual.Convention.Bond),
					 BusinessDayConvention.ModifiedFollowing,
					 100.0,
					 issue,
					 callSchedule);
    }

    static public CallableFixedRateBond mkFHMLC()
    {
	double coupon = 0.0275;
	Date issue = new Date(26, Month.September, 2016);
	Date matur = new Date(26, Month.September, 2036);

	DoubleVector coupons = new DoubleVector();
        coupons.add(coupon);

        CallabilitySchedule callSchedule =
	    new CallabilitySchedule();

        double callPrice = 100.; 
	long numberOfCallDates = 76;
        Date callDate =new Date(26, Month.September, 2017);
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

	return new CallableFixedRateBond(1,
					 100.0,
					 new Schedule(issue,
						      matur,
						      new Period(Frequency.Semiannual),
						      new UnitedStates(UnitedStates.Market.GovernmentBond),
						      BusinessDayConvention.Unadjusted,
						      BusinessDayConvention.Unadjusted,						      
						      DateGeneration.Rule.Backward,
						      false),
					 coupons,
					 new Thirty360(Thirty360.Convention.USA),
					 BusinessDayConvention.ModifiedFollowing,
					 100.0,
					 issue,
					 callSchedule);
    }

    static public CallableFixedRateBond mkWRI()
    {
	double coupon = 0.0646;
	Date issue = new Date(11, Month.August, 1998);
	Date matur = new Date(11, Month.August, 2028);

	DoubleVector coupons = new DoubleVector();
        coupons.add(coupon);

        CallabilitySchedule callSchedule =
	    new CallabilitySchedule();

        double callPrice = 100.; 
        Date callDate =new Date(11, Month.August, 2018);

	CallabilityPrice myPrice=
	    new CallabilityPrice(callPrice,
				 CallabilityPrice.Type.Clean);
	callSchedule.add(new Callability(myPrice,
					 Callability.Put,
					 callDate));

	return new CallableFixedRateBond(1,
					 100.0,
					 new Schedule(issue,
						      matur,
						      new Period(Frequency.Semiannual),
						      new UnitedStates(UnitedStates.Market.GovernmentBond),
						      BusinessDayConvention.Unadjusted,
						      BusinessDayConvention.Unadjusted,						      
						      DateGeneration.Rule.Backward,
						      false),
					 coupons,
					 new Thirty360(Thirty360.Convention.USA),
					 BusinessDayConvention.ModifiedFollowing,
					 100.0,
					 issue,
					 callSchedule);
    }        

    public static TreeCallableFixedRateBondEngine mkEngine(RelinkableYieldTermStructureHandle yc)
    {
	double reversionParameter=0.03;
	double sigma=0.012;
	int gridIntervals=40;

	HullWhite hw0 = new HullWhite(yc,
				      reversionParameter,
				      sigma);
	return new TreeCallableFixedRateBondEngine(hw0,
						   gridIntervals);
    }

    public static void printPricing(CallableFixedRateBond b,
				    double OAS,
				    RelinkableYieldTermStructureHandle ych,
				    Date settlementDate)
    {
	double cleanOAS=b.cleanPriceOAS(OAS * 1e-4, ych, new ActualActual(ActualActual.Convention.Bond), Compounding.Compounded, Frequency.Semiannual, settlementDate);
        System.out.printf("OAS (round-trip): %6.4f \t NPV: %8.4f \t Clean: %8.4f \t Clean w/OAS:%8.4f \t EffDur: %6.4f \t EffCvx: %6.4f  \n",
			  1e4* b.OAS(cleanOAS, ych, new ActualActual(ActualActual.Convention.Bond), Compounding.Compounded, Frequency.Semiannual), 
			  b.NPV(),
			  b.cleanPrice(),
			  cleanOAS,
			  b.effectiveDuration(OAS * 1e-4, ych, new ActualActual(ActualActual.Convention.Bond), Compounding.Compounded, Frequency.Semiannual),
			  b.effectiveConvexity(OAS * 1e-4, ych, new ActualActual(ActualActual.Convention.Bond), Compounding.Compounded, Frequency.Semiannual, 1e-2)			  
			  );
    }

    /** Example OAS calculation function that is thread safe

	Note that separate bonds, engines and a relinkable handle are
	made for each valuation which is a requirement for thread
	safety: i.e., any bond, execution engine or relinkable handle
	should called by a single thread at any one time
    */
    public static double valuePar(double OAS,
				  YieldTermStructure yc)
    {
	CallableFixedRateBond b2 =mkWRI();
	RelinkableYieldTermStructureHandle ych= new RelinkableYieldTermStructureHandle(yc);
	TreeCallableFixedRateBondEngine eng =mkEngine(ych);
	b2.setPricingEngine(eng);
	double cleanOAS=b2.cleanPriceOAS(OAS * 1e-4, ych, new ActualActual(ActualActual.Convention.Bond), Compounding.Compounded, Frequency.Semiannual);
	return 1e4* b2.OAS(cleanOAS, ych, new ActualActual(ActualActual.Convention.Bond), Compounding.Compounded, Frequency.Semiannual);
    }
			   
			   

    public static void main(String[] args) throws Exception {

	Settings.instance().setEvaluationDate(today);

	YieldTermStructure yc=mkYC();

	CallableFixedRateBond ust1= mkUST1();

	RelinkableYieldTermStructureHandle ych= new RelinkableYieldTermStructureHandle(yc);

	TreeCallableFixedRateBondEngine eng =mkEngine(ych);

	ust1.setPricingEngine(eng);

	System.out.println("* UST1 *:");	
	printPricing(ust1, 0.3, ych, today);

	Calendar nullCalendar = new NullCalendar();
	CallableFixedRateBond blrdg= mkBLRDG();
	blrdg.setPricingEngine(eng);

	System.out.println("* BLRDG *:");	
	printPricing(blrdg, 28.01, ych, nullCalendar.advance(today, 3, TimeUnit.Days));

	CallableFixedRateBond fhmlc = mkFHMLC();
	fhmlc.setPricingEngine(eng);
	System.out.println("* FHMLC *:");	
	printPricing(fhmlc, -42.0, ych, today);

	CallableFixedRateBond wri = mkWRI();
	wri.setPricingEngine(eng);
	System.out.println("* WRI *:");	
	printPricing(wri, 558.8, ych, today);

	// Example multi-threaded parallel execution
	IntStream.range(0,100).parallel()
	    .mapToDouble( e -> { return valuePar(558.8-0.0000001*e, yc); })
	    .forEachOrdered( e ->  { System.out.print( e +" ");} );
	
    }

    
}
