/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2010, 2011 Lluis Pujol Bajador
 Copyright (C) 2017, 2018, 2019, 2020 Matthias Lungwitz

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


#ifndef quantlib_cash_flows_i
#define quantlib_cash_flows_i

%include date.i
%include types.i
%include calendars.i
%include daycounters.i
%include indexes.i
%include termstructures.i
%include scheduler.i
%include vectors.i
%include volatilities.i

%{
using QuantLib::CashFlow;
%}

%shared_ptr(CashFlow)
class CashFlow : public Observable {
  private:
    CashFlow();
  public:
    Real amount() const;
    Date date() const;
    bool hasOccurred(const Date& refDate = Date()) const;
};

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<CashFlow> )
#endif
%template(Leg) std::vector<ext::shared_ptr<CashFlow> >;
typedef std::vector<ext::shared_ptr<CashFlow> > Leg;


// implementations

%{
using QuantLib::SimpleCashFlow;
using QuantLib::IndexedCashFlow;
using QuantLib::Redemption;
using QuantLib::AmortizingPayment;
using QuantLib::Coupon;
using QuantLib::FixedRateCoupon;
using QuantLib::Leg;
using QuantLib::FloatingRateCoupon;
using QuantLib::OvernightIndexedCoupon;
%}

%shared_ptr(SimpleCashFlow)
class SimpleCashFlow : public CashFlow {
  public:
    SimpleCashFlow(Real amount, const Date& date);
};

%shared_ptr(Redemption)
class Redemption : public CashFlow {
  public:
    Redemption(Real amount, const Date& date);
};

%shared_ptr(AmortizingPayment)
class AmortizingPayment : public CashFlow {
  public:
    AmortizingPayment(Real amount, const Date& date);
};


%shared_ptr(IndexedCashFlow)
class IndexedCashFlow : public CashFlow {
  public:
    IndexedCashFlow(Real notional,
                    const ext::shared_ptr<Index>& index,
                    const Date& baseDate,
                    const Date& fixingDate,
                    const Date& paymentDate,
                    bool growthOnly = false);
    Real notional() const;
    Date baseDate() const;
    Date fixingDate() const;
    ext::shared_ptr<Index> index() const;
    bool growthOnly() const;
};

%inline %{
    ext::shared_ptr<IndexedCashFlow> as_indexed_cashflow(const ext::shared_ptr<CashFlow>& cf) {
        return ext::dynamic_pointer_cast<IndexedCashFlow>(cf);
    }
%}


%shared_ptr(Coupon)
class Coupon : public CashFlow {
  private:
    Coupon();
  public:
    Real nominal() const;
    Date accrualStartDate() const;
    Date accrualEndDate() const;
    Date referencePeriodStart() const;
    Date referencePeriodEnd() const;
    Date exCouponDate() const;
    Real rate() const;
    Time accrualPeriod() const;
    BigInteger accrualDays() const;
    DayCounter dayCounter() const;
    Real accruedAmount(const Date& date) const;
};

%inline %{
    ext::shared_ptr<Coupon> as_coupon(const ext::shared_ptr<CashFlow>& cf) {
        return ext::dynamic_pointer_cast<Coupon>(cf);
    }
%}


%shared_ptr(FixedRateCoupon)
class FixedRateCoupon : public Coupon {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") FixedRateCoupon;
    #endif
  public:
    FixedRateCoupon(const Date& paymentDate, Real nominal,
                    Rate rate, const DayCounter& dayCounter,
                    const Date& startDate, const Date& endDate,
                    const Date& refPeriodStart = Date(),
                    const Date& refPeriodEnd = Date(),
                    const Date& exCouponDate = Date());
    InterestRate interestRate() const;
};

%inline %{
    ext::shared_ptr<FixedRateCoupon> as_fixed_rate_coupon(
                                      const ext::shared_ptr<CashFlow>& cf) {
        return ext::dynamic_pointer_cast<FixedRateCoupon>(cf);
    }
%}


%{
using QuantLib::FloatingRateCouponPricer;
%}

%shared_ptr(FloatingRateCouponPricer)
class FloatingRateCouponPricer {
  private:
    FloatingRateCouponPricer();
  public:
    virtual Real swapletPrice() const;
    virtual Rate swapletRate() const;
    virtual Real capletPrice(Rate effectiveCap) const;
    virtual Rate capletRate(Rate effectiveCap) const;
    virtual Real floorletPrice(Rate effectiveFloor) const;
    virtual Rate floorletRate(Rate effectiveFloor) const;
};

void setCouponPricer(const Leg&,
                     const ext::shared_ptr<FloatingRateCouponPricer>&);

%shared_ptr(FloatingRateCoupon)
class FloatingRateCoupon : public Coupon {
  private:
    FloatingRateCoupon();
  public:
    Date fixingDate() const;
    Integer fixingDays() const;
    bool isInArrears() const;
    Real gearing() const;
    Rate spread() const;
    Rate indexFixing() const;
    Rate adjustedFixing() const;
    Rate convexityAdjustment() const;
    Real price(const Handle<YieldTermStructure>& discountCurve) const;
    ext::shared_ptr<InterestRateIndex> index() const;
    void setPricer(const ext::shared_ptr<FloatingRateCouponPricer>& p);
};

%inline %{
    ext::shared_ptr<FloatingRateCoupon> as_floating_rate_coupon(
                                      const ext::shared_ptr<CashFlow>& cf) {
        return ext::dynamic_pointer_cast<FloatingRateCoupon>(cf);
    }
%}

%shared_ptr(OvernightIndexedCoupon)
class OvernightIndexedCoupon : public FloatingRateCoupon {
  public:
    OvernightIndexedCoupon(
                const Date& paymentDate,
                Real nominal,
                const Date& startDate,
                const Date& endDate,
                const ext::shared_ptr<OvernightIndex>& overnightIndex,
                Real gearing = 1.0,
                Spread spread = 0.0,
                const Date& refPeriodStart = Date(),
                const Date& refPeriodEnd = Date(),
                const DayCounter& dayCounter = DayCounter(),
                bool telescopicValueDates = false);
    const std::vector<Date>& fixingDates() const;
    const std::vector<Time>& dt() const;
    const std::vector<Rate>& indexFixings() const;
    const std::vector<Date>& valueDates() const;
};

%{
using QuantLib::CappedFlooredCoupon;
%}

%shared_ptr(CappedFlooredCoupon)
class CappedFlooredCoupon : public FloatingRateCoupon {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CappedFlooredCoupon;
    #endif
  public:
    CappedFlooredCoupon(const ext::shared_ptr<FloatingRateCoupon>& underlying,
                        Rate cap = Null<Rate>(),
                        Rate floor = Null<Rate>());
    Rate cap() const;
    Rate floor() const;
    Rate effectiveCap() const;
    Rate effectiveFloor() const;
    bool isCapped() const;
    bool isFloored() const;
    void setPricer(const ext::shared_ptr<FloatingRateCouponPricer>& p);
};


// specialized floating-rate coupons

%{
using QuantLib::IborCoupon;
using QuantLib::CappedFlooredIborCoupon;
%}

%shared_ptr(IborCoupon)
class IborCoupon : public FloatingRateCoupon {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") IborCoupon;
    #endif
  public:
    IborCoupon(const Date& paymentDate, Real nominal,
               const Date& startDate, const Date& endDate,
               Integer fixingDays,
               ext::shared_ptr<IborIndex>& index,
               Real gearing = 1.0, Spread spread = 0.0,
               const Date& refPeriodStart = Date(),
               const Date& refPeriodEnd = Date(),
               const DayCounter& dayCounter = DayCounter(),
               bool isInArrears = false,
               const Date& exCouponDate = Date());
    static void createAtParCoupons();
    static void createIndexedCoupons();
    static bool usingAtParCoupons();
};

%shared_ptr(CappedFlooredIborCoupon)
class CappedFlooredIborCoupon : public CappedFlooredCoupon {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CappedFlooredIborCoupon;
    #endif
  public:
    CappedFlooredIborCoupon(const Date& paymentDate, Real nominal,
                            const Date& startDate, const Date& endDate,
                            Integer fixingDays,
                            ext::shared_ptr<IborIndex>& index,
                            Real gearing = 1.0, Spread spread = 0.0,
                            const Rate cap = Null<Rate>(),
                            const Rate floor = Null<Rate>(),
                            const Date& refPeriodStart = Date(),
                            const Date& refPeriodEnd = Date(),
                            const DayCounter& dayCounter = DayCounter(),
                            bool isInArrears = false,
                            const Date& exCouponDate = Date());
};

%{
using QuantLib::IborCouponPricer;
using QuantLib::BlackIborCouponPricer;
%}

%shared_ptr(IborCouponPricer)
class IborCouponPricer : public FloatingRateCouponPricer {
  private:
    IborCouponPricer();
  public:
    Handle<OptionletVolatilityStructure> capletVolatility() const;
    void setCapletVolatility(const Handle<OptionletVolatilityStructure>& v =
                                     Handle<OptionletVolatilityStructure>());
};

%shared_ptr(BlackIborCouponPricer)
class BlackIborCouponPricer : public IborCouponPricer {
  public:
    enum TimingAdjustment { Black76, BivariateLognormal };
    BlackIborCouponPricer(const Handle<OptionletVolatilityStructure>& v =
                                    Handle<OptionletVolatilityStructure>(),
                                    const TimingAdjustment timingAdjustment = Black76,
                                    const Handle<Quote> correlation =
                                    Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(1.0))));
};

%{
using QuantLib::CmsCoupon;
using QuantLib::CappedFlooredCmsCoupon;
using QuantLib::CmsSpreadCoupon;
using QuantLib::CappedFlooredCmsSpreadCoupon;
%}

%shared_ptr(CmsCoupon)
class CmsCoupon : public FloatingRateCoupon {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CmsCoupon;
    #endif
  public:
    CmsCoupon(const Date& paymentDate, Real nominal,
              const Date& startDate, const Date& endDate,
              Integer fixingDays, const ext::shared_ptr<SwapIndex>& index,
              Real gearing = 1.0, Spread spread = 0.0,
              const Date& refPeriodStart = Date(),
              const Date& refPeriodEnd = Date(),
              const DayCounter& dayCounter = DayCounter(),
              bool isInArrears = false,
              const Date& exCouponDate = Date());
};

%shared_ptr(CmsSpreadCoupon)
class CmsSpreadCoupon : public FloatingRateCoupon {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CmsSpreadCoupon;
    #endif
  public:
    CmsSpreadCoupon(const Date& paymentDate,
                    Real nominal,
                    const Date& startDate,
                    const Date& endDate,
                    Natural fixingDays,
                    const ext::shared_ptr<SwapSpreadIndex>& index,
                    Real gearing = 1.0,
                    Spread spread = 0.0,
                    const Date& refPeriodStart = Date(),
                    const Date& refPeriodEnd = Date(),
                    const DayCounter& dayCounter = DayCounter(),
                    bool isInArrears = false,
                    const Date& exCouponDate = Date());
};

%{
using QuantLib::CmsCouponPricer;
using QuantLib::AnalyticHaganPricer;
using QuantLib::NumericHaganPricer;
using QuantLib::GFunctionFactory;
using QuantLib::LinearTsrPricer;
using QuantLib::CmsSpreadCouponPricer;
using QuantLib::LognormalCmsSpreadPricer;
%}

%shared_ptr(CmsCouponPricer)
class CmsCouponPricer : public FloatingRateCouponPricer {
  private:
    CmsCouponPricer();
  public:
    Handle<SwaptionVolatilityStructure> swaptionVolatility() const;
    void setSwaptionVolatility(const Handle<SwaptionVolatilityStructure>& v =
                                      Handle<SwaptionVolatilityStructure>());
};

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<CmsCouponPricer> )
#endif
namespace std {
    %template(CmsCouponPricerVector)
        vector<ext::shared_ptr<CmsCouponPricer> >;
}

class GFunctionFactory {
  private:
    GFunctionFactory();
  public:
    enum YieldCurveModel { Standard,
                           ExactYield,
                           ParallelShifts,
                           NonParallelShifts };
};

%shared_ptr(AnalyticHaganPricer)
class AnalyticHaganPricer : public CmsCouponPricer {
  public:
    AnalyticHaganPricer(const Handle<SwaptionVolatilityStructure>& v,
                        GFunctionFactory::YieldCurveModel model,
                        const Handle<Quote>& meanReversion);
};

%shared_ptr(NumericHaganPricer)
class NumericHaganPricer : public CmsCouponPricer {
  public:
    NumericHaganPricer(const Handle<SwaptionVolatilityStructure>& v,
                       GFunctionFactory::YieldCurveModel model,
                       const Handle<Quote>& meanReversion,
                       Rate lowerLimit = 0.0,
                       Rate upperLimit = 1.0,
                       Real precision = 1.0e-6);
};

%shared_ptr(CappedFlooredCmsCoupon)
class CappedFlooredCmsCoupon: public CappedFlooredCoupon {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CappedFlooredCoupon;
    #endif
  public:
    CappedFlooredCmsCoupon(
                  const Date& paymentDate, Real nominal,
                  const Date& startDate, const Date& endDate,
                  Natural fixingDays, const ext::shared_ptr<SwapIndex>& index,
                  Real gearing = 1.0, Spread spread = 0.0,
                  const Rate cap = Null<Rate>(),
                  const Rate floor = Null<Rate>(),
                  const Date& refPeriodStart = Date(),
                  const Date& refPeriodEnd = Date(),
                  const DayCounter& dayCounter = DayCounter(),
                  bool isInArrears = false,
                  const Date& exCouponDate = Date());
};

%shared_ptr(CappedFlooredCmsSpreadCoupon)
class CappedFlooredCmsSpreadCoupon: public CappedFlooredCoupon {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CappedFlooredCoupon;
    #endif
  public:
    CappedFlooredCmsSpreadCoupon(
                  const Date& paymentDate, Real nominal,
                  const Date& startDate, const Date& endDate,
                  Natural fixingDays,
                  const ext::shared_ptr<SwapSpreadIndex>& index,
                  Real gearing = 1.0, Spread spread = 0.0,
                  const Rate cap = Null<Rate>(),
                  const Rate floor = Null<Rate>(),
                  const Date& refPeriodStart = Date(),
                  const Date& refPeriodEnd = Date(),
                  const DayCounter& dayCounter = DayCounter(),
                  bool isInArrears = false,
                  const Date& exCouponDate = Date());
};

%rename (LinearTsrPricerSettings) LinearTsrPricer::Settings;
%feature ("flatnested") Settings;

%shared_ptr(LinearTsrPricer)
class LinearTsrPricer : public CmsCouponPricer {
  public:
    struct Settings {

        Settings();
        Settings &withRateBound(const Real lowerRateBound = 0.0001,
                                const Real upperRateBound = 2.0000);
        Settings &withVegaRatio(const Real vegaRatio = 0.01);                                
        Settings &withVegaRatio(const Real vegaRatio,
                                const Real lowerRateBound,
                                const Real upperRateBound);
        Settings &withPriceThreshold(const Real priceThreshold = 1.0E-8);
        Settings &withPriceThreshold(const Real priceThreshold,
                                     const Real lowerRateBound,
                                     const Real upperRateBound);
        Settings &withBSStdDevs(const Real stdDevs = 3.0);
        Settings &withBSStdDevs(const Real stdDevs,
                        const Real lowerRateBound,
                        const Real upperRateBound);
        enum Strategy {
            RateBound,
            VegaRatio,
            PriceThreshold,
            BSStdDevs
        };
        
    };
  
    LinearTsrPricer(
            const Handle<SwaptionVolatilityStructure> &swaptionVol,
            const Handle<Quote> &meanReversion,
            const Handle<YieldTermStructure> &couponDiscountCurve =
                                                 Handle<YieldTermStructure>(),
            const LinearTsrPricer::Settings &settings =
                                                LinearTsrPricer::Settings());
};

%shared_ptr(CmsSpreadCouponPricer)
class CmsSpreadCouponPricer : public FloatingRateCouponPricer {
  private:
    CmsSpreadCouponPricer();
  public:
    Handle<Quote> correlation() const;
    void setCorrelation(const Handle<Quote> &correlation = Handle<Quote>());
};

%shared_ptr(LognormalCmsSpreadPricer)
class LognormalCmsSpreadPricer : public CmsSpreadCouponPricer {
  public:
    LognormalCmsSpreadPricer(
            const ext::shared_ptr<CmsCouponPricer>& cmsPricer,
            const Handle<Quote> &correlation,
            const Handle<YieldTermStructure> &couponDiscountCurve =
                Handle<YieldTermStructure>(),
            const Size IntegrationPoints = 16,
            const boost::optional<VolatilityType> volatilityType = boost::none,
            const Real shift1 = Null<Real>(), const Real shift2 = Null<Real>());
    Real swapletPrice() const;
    Rate swapletRate() const;
    Real capletPrice(Rate effectiveCap) const;
    Rate capletRate(Rate effectiveCap) const;
    Real floorletPrice(Rate effectiveFloor) const;
    Rate floorletRate(Rate effectiveFloor) const;
};

// cash flow vector builders

%{
Leg _FixedRateLeg(const Schedule& schedule,
                  const DayCounter& dayCount,
                  const std::vector<Real>& nominals,
                  const std::vector<Rate>& couponRates,
                  BusinessDayConvention paymentAdjustment = Following,
                  const DayCounter& firstPeriodDayCount = DayCounter(),
                  const Period& exCouponPeriod = Period(),
                  const Calendar& exCouponCalendar = Calendar(),
                  BusinessDayConvention exCouponConvention = Unadjusted,
                  bool exCouponEndOfMonth = false) {
    return QuantLib::FixedRateLeg(schedule)
        .withNotionals(nominals)
        .withCouponRates(couponRates,dayCount)
        .withPaymentAdjustment(paymentAdjustment)
        .withFirstPeriodDayCounter(firstPeriodDayCount)
        .withExCouponPeriod(exCouponPeriod,
                            exCouponCalendar,
                            exCouponConvention,
                            exCouponEndOfMonth);
}
%}
#if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
%feature("kwargs") _FixedRateLeg;
#endif
%rename(FixedRateLeg) _FixedRateLeg;
Leg _FixedRateLeg(const Schedule& schedule,
                  const DayCounter& dayCount,
                  const std::vector<Real>& nominals,
                  const std::vector<Rate>& couponRates,
                  BusinessDayConvention paymentAdjustment = Following,
                  const DayCounter& firstPeriodDayCount = DayCounter(),
                  const Period& exCouponPeriod = Period(),
                  const Calendar& exCouponCalendar = Calendar(),
                  BusinessDayConvention exCouponConvention = Unadjusted,
                  bool exCouponEndOfMonth = false);

%{
Leg _IborLeg(const std::vector<Real>& nominals,
             const Schedule& schedule,
             const ext::shared_ptr<IborIndex>& index,
             const DayCounter& paymentDayCounter = DayCounter(),
             const BusinessDayConvention paymentConvention = Following,
             const std::vector<Natural>& fixingDays = std::vector<Natural>(),
             const std::vector<Real>& gearings = std::vector<Real>(),
             const std::vector<Spread>& spreads = std::vector<Spread>(),
             const std::vector<Rate>& caps = std::vector<Rate>(),
             const std::vector<Rate>& floors = std::vector<Rate>(),
             bool isInArrears = false,
             const Period& exCouponPeriod = Period(),
             const Calendar& exCouponCalendar = Calendar(),
             BusinessDayConvention exCouponConvention = Unadjusted,
             bool exCouponEndOfMonth = false) {
    return QuantLib::IborLeg(schedule, index)
        .withNotionals(nominals)
        .withPaymentDayCounter(paymentDayCounter)
        .withPaymentAdjustment(paymentConvention)
        .withFixingDays(fixingDays)
        .withGearings(gearings)
        .withSpreads(spreads)
        .withCaps(caps)
        .withFloors(floors)
        .inArrears(isInArrears)
        .withExCouponPeriod(exCouponPeriod,
                            exCouponCalendar,
                            exCouponConvention,
                            exCouponEndOfMonth);
}
%}
#if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
%feature("kwargs") _IborLeg;
#endif
%rename(IborLeg) _IborLeg;
Leg _IborLeg(const std::vector<Real>& nominals,
             const Schedule& schedule,
             const ext::shared_ptr<IborIndex>& index,
             const DayCounter& paymentDayCounter = DayCounter(),
             const BusinessDayConvention paymentConvention = Following,
             const std::vector<Natural>& fixingDays = std::vector<Natural>(),
             const std::vector<Real>& gearings = std::vector<Real>(),
             const std::vector<Spread>& spreads = std::vector<Spread>(),
             const std::vector<Rate>& caps = std::vector<Rate>(),
             const std::vector<Rate>& floors = std::vector<Rate>(),
             bool isInArrears = false,
             const Period& exCouponPeriod = Period(),
             const Calendar& exCouponCalendar = Calendar(),
             BusinessDayConvention exCouponConvention = Unadjusted,
             bool exCouponEndOfMonth = false);

%{
Leg _OvernightLeg(const std::vector<Real>& nominals,
             const Schedule& schedule,
             const ext::shared_ptr<OvernightIndex>& index,
             const DayCounter& paymentDayCounter = DayCounter(),
             const BusinessDayConvention paymentConvention = Following,
             const std::vector<Real>& gearings = std::vector<Real>(),
             const std::vector<Spread>& spreads = std::vector<Spread>(),
             bool telescopicValueDates = false) {
    return QuantLib::OvernightLeg(schedule, index)
        .withNotionals(nominals)
        .withPaymentDayCounter(paymentDayCounter)
        .withPaymentAdjustment(paymentConvention)
        .withGearings(gearings)
        .withSpreads(spreads)
        .withTelescopicValueDates(telescopicValueDates);
}
%}
#if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
%feature("kwargs") _OvernightLeg;
#endif
%rename(OvernightLeg) _OvernightLeg;
Leg _OvernightLeg(const std::vector<Real>& nominals,
             const Schedule& schedule,
             const ext::shared_ptr<OvernightIndex>& index,
             const DayCounter& paymentDayCounter = DayCounter(),
             const BusinessDayConvention paymentConvention = Following,
             const std::vector<Real>& gearings = std::vector<Real>(),
             const std::vector<Spread>& spreads = std::vector<Spread>(),
             bool telescopicValueDates = false);

%{
Leg _CmsLeg(const std::vector<Real>& nominals,
            const Schedule& schedule,
            const ext::shared_ptr<SwapIndex>& index,
            const DayCounter& paymentDayCounter = DayCounter(),
            const BusinessDayConvention paymentConvention = Following,
            const std::vector<Natural>& fixingDays = std::vector<Natural>(),
            const std::vector<Real>& gearings = std::vector<Real>(),
            const std::vector<Spread>& spreads = std::vector<Spread>(),
            const std::vector<Rate>& caps = std::vector<Rate>(),
            const std::vector<Rate>& floors = std::vector<Rate>(),
            bool isInArrears = false,
            const Period& exCouponPeriod = Period(),
            const Calendar& exCouponCalendar = Calendar(),
            const BusinessDayConvention exCouponConvention = Unadjusted,
            bool exCouponEndOfMonth = false) {
    return QuantLib::CmsLeg(schedule, index)
        .withNotionals(nominals)
        .withPaymentDayCounter(paymentDayCounter)
        .withPaymentAdjustment(paymentConvention)
        .withFixingDays(fixingDays)
        .withGearings(gearings)
        .withSpreads(spreads)
        .withCaps(caps)
        .withFloors(floors)
        .withExCouponPeriod(exCouponPeriod, exCouponCalendar,
                            exCouponConvention, exCouponEndOfMonth)
        .inArrears(isInArrears);
}
%}
#if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
%feature("kwargs") _CmsLeg;
#endif
%rename(CmsLeg) _CmsLeg;
Leg _CmsLeg(const std::vector<Real>& nominals,
            const Schedule& schedule,
            const ext::shared_ptr<SwapIndex>& index,
            const DayCounter& paymentDayCounter = DayCounter(),
            const BusinessDayConvention paymentConvention = Following,
            const std::vector<Natural>& fixingDays = std::vector<Natural>(),
            const std::vector<Real>& gearings = std::vector<Real>(),
            const std::vector<Spread>& spreads = std::vector<Spread>(),
            const std::vector<Rate>& caps = std::vector<Rate>(),
            const std::vector<Rate>& floors = std::vector<Rate>(),
            bool isInArrears = false,
            const Period& exCouponPeriod = Period(),
            const Calendar& exCouponCalendar = Calendar(),
            const BusinessDayConvention exCouponConvention = Unadjusted,
            bool exCouponEndOfMonth = false);

%{
Leg _CmsZeroLeg(const std::vector<Real>& nominals,
                const Schedule& schedule,
                const ext::shared_ptr<SwapIndex>& index,
                const DayCounter& paymentDayCounter = DayCounter(),
                const BusinessDayConvention paymentConvention = Following,
                const std::vector<Natural>& fixingDays = std::vector<Natural>(),
                const std::vector<Real>& gearings = std::vector<Real>(),
                const std::vector<Spread>& spreads = std::vector<Spread>(),
                const std::vector<Rate>& caps = std::vector<Rate>(),
                const std::vector<Rate>& floors = std::vector<Rate>(),
                const Period& exCouponPeriod = Period(),
                const Calendar& exCouponCalendar = Calendar(),
                const BusinessDayConvention exCouponConvention = Unadjusted,
                bool exCouponEndOfMonth = false) {
    return QuantLib::CmsLeg(schedule, index)
        .withNotionals(nominals)
        .withPaymentDayCounter(paymentDayCounter)
        .withPaymentAdjustment(paymentConvention)
        .withFixingDays(fixingDays)
        .withGearings(gearings)
        .withSpreads(spreads)
        .withCaps(caps)
        .withFloors(floors)
        .withExCouponPeriod(exCouponPeriod, exCouponCalendar,
                            exCouponConvention, exCouponEndOfMonth)
        .withZeroPayments();
}
%}
#if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
%feature("kwargs") _CmsZeroLeg;
#endif
%rename(CmsZeroLeg) _CmsZeroLeg;
Leg _CmsZeroLeg(const std::vector<Real>& nominals,
                const Schedule& schedule,
                const ext::shared_ptr<SwapIndex>& index,
                const DayCounter& paymentDayCounter = DayCounter(),
                const BusinessDayConvention paymentConvention = Following,
                const std::vector<Natural>& fixingDays = std::vector<Natural>(),
                const std::vector<Real>& gearings = std::vector<Real>(),
                const std::vector<Spread>& spreads = std::vector<Spread>(),
                const std::vector<Rate>& caps = std::vector<Rate>(),
                const std::vector<Rate>& floors = std::vector<Rate>(),
                const Period& exCouponPeriod = Period(),
                const Calendar& exCouponCalendar = Calendar(),
                const BusinessDayConvention exCouponConvention = Unadjusted,
                bool exCouponEndOfMonth = false);

%{
Leg _CmsSpreadLeg(const std::vector<Real>& nominals,
            const Schedule& schedule,
            const ext::shared_ptr<SwapSpreadIndex>& index,
            const DayCounter& paymentDayCounter = DayCounter(),
            const BusinessDayConvention paymentConvention = Following,
            const std::vector<Natural>& fixingDays = std::vector<Natural>(),
            const std::vector<Real>& gearings = std::vector<Real>(),
            const std::vector<Spread>& spreads = std::vector<Spread>(),
            const std::vector<Rate>& caps = std::vector<Rate>(),
            const std::vector<Rate>& floors = std::vector<Rate>(),
            bool isInArrears = false) {
    return QuantLib::CmsSpreadLeg(schedule, index)
        .withNotionals(nominals)
        .withPaymentDayCounter(paymentDayCounter)
        .withPaymentAdjustment(paymentConvention)
        .withFixingDays(fixingDays)
        .withGearings(gearings)
        .withSpreads(spreads)
        .withCaps(caps)
        .withFloors(floors)
        .inArrears(isInArrears);
}
%}
#if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
%feature("kwargs") _CmsSpreadLeg;
#endif
%rename(CmsSpreadLeg) _CmsSpreadLeg;
Leg _CmsSpreadLeg(const std::vector<Real>& nominals,
            const Schedule& schedule,
            const ext::shared_ptr<SwapSpreadIndex>& index,
            const DayCounter& paymentDayCounter = DayCounter(),
            const BusinessDayConvention paymentConvention = Following,
            const std::vector<Natural>& fixingDays = std::vector<Natural>(),
            const std::vector<Real>& gearings = std::vector<Real>(),
            const std::vector<Spread>& spreads = std::vector<Spread>(),
            const std::vector<Rate>& caps = std::vector<Rate>(),
            const std::vector<Rate>& floors = std::vector<Rate>(),
            bool isInArrears = false);
                
// cash-flow analysis

%{
using QuantLib::CashFlows;
using QuantLib::Duration;
%}

struct Duration {
    enum Type { Simple, Macaulay, Modified };
};

class CashFlows {
    #if defined(SWIGPYTHON)
    %rename("yieldRate")   yield;
    #endif
  private:
    CashFlows();
    CashFlows(const CashFlows&);
  public:
    static Date startDate(const Leg &);
    static Date maturityDate(const Leg &);
    static Date
        previousCashFlowDate(const Leg& leg,
                             bool includeSettlementDateFlows,
                             Date settlementDate = Date());
    static Date
        nextCashFlowDate(const Leg& leg,
                         bool includeSettlementDateFlows,
                         Date settlementDate = Date());

    %extend {

        static ext::shared_ptr<CashFlow>
        previousCashFlow(const Leg& leg,
                         bool includeSettlementDateFlows,
                         Date settlementDate = Date()) {
            Leg::const_reverse_iterator i =
                QuantLib::CashFlows::previousCashFlow(
                    leg, includeSettlementDateFlows, settlementDate);

            if (i == leg.rend())
                return ext::shared_ptr<CashFlow>();
            else
                return *i;
        }

        static ext::shared_ptr<CashFlow>
        nextCashFlow(const Leg& leg,
                     bool includeSettlementDateFlows,
                     Date settlementDate = Date()) {
            Leg::const_iterator i =
                QuantLib::CashFlows::nextCashFlow(
                    leg, includeSettlementDateFlows, settlementDate);

            if (i == leg.end())
                return ext::shared_ptr<CashFlow>();
            else
                return *i;
        }

    }

    %extend {
        static Real npv(
                   const Leg& leg,
                   const ext::shared_ptr<YieldTermStructure>& discountCurve,
                   Spread zSpread,
                   const DayCounter &dayCounter,
                   Compounding compounding,
                   Frequency frequency,
                   bool includeSettlementDateFlows,
                   const Date& settlementDate = Date(),
                   const Date& npvDate = Date()) {
            return QuantLib::CashFlows::npv(leg, discountCurve,
                                            zSpread,
                                            dayCounter,
                                            compounding,
                                            frequency,
                                            includeSettlementDateFlows,
                                            settlementDate,
                                            npvDate);
        }
        static Real npv(
                   const Leg& leg,
                   const Handle<YieldTermStructure>& discountCurve,
                   bool includeSettlementDateFlows,
                   const Date& settlementDate = Date(),
                   const Date& npvDate = Date()) {
            return QuantLib::CashFlows::npv(leg, **discountCurve,
                                            includeSettlementDateFlows,
                                            settlementDate, npvDate);
        }
    }
    static Real npv(const Leg&,
                    const InterestRate&,
                    bool includeSettlementDateFlows,
                    Date settlementDate = Date(),
                    Date npvDate = Date());

    static Real npv(const Leg&,
                    Rate yield,
                    const DayCounter&dayCounter,
                    Compounding compounding,
                    Frequency frequency,
                    bool includeSettlementDateFlows,
                    Date settlementDate = Date(),
                    Date npvDate = Date());
    %extend {
        static Real bps(
                   const Leg& leg,
                   const ext::shared_ptr<YieldTermStructure>& discountCurve,
                   bool includeSettlementDateFlows,
                   const Date& settlementDate = Date(),
                   const Date& npvDate = Date()) {
            return QuantLib::CashFlows::bps(leg, *discountCurve,
                                            includeSettlementDateFlows,
                                            settlementDate, npvDate);
        }
        static Real bps(
                   const Leg& leg,
                   const Handle<YieldTermStructure>& discountCurve,
                   bool includeSettlementDateFlows,
                   const Date& settlementDate = Date(),
                   const Date& npvDate = Date()) {
            return QuantLib::CashFlows::bps(leg, **discountCurve,
                                            includeSettlementDateFlows,
                                            settlementDate, npvDate);
        }
    }
    static Real bps(const Leg&,
                    const InterestRate &,
                    bool includeSettlementDateFlows,
                    Date settlementDate = Date(),
                    Date npvDate = Date());
    static Real bps(const Leg&,
                    Rate yield,
                    const DayCounter&dayCounter,
                    Compounding compounding,
                    Frequency frequency,
                    bool includeSettlementDateFlows,
                    Date settlementDate = Date(),
                    Date npvDate = Date());

    %extend {
        static Rate atmRate(
                   const Leg& leg,
                   const ext::shared_ptr<YieldTermStructure>& discountCurve,
                   bool includeSettlementDateFlows,
                   const Date& settlementDate = Date(),
                   const Date& npvDate = Date(),
                   Real npv = Null<Real>()) {
            return QuantLib::CashFlows::atmRate(leg, *discountCurve,
                                                includeSettlementDateFlows,
                                                settlementDate, npvDate,
                                                npv);
        }
    }
    static Rate yield(const Leg&,
                      Real npv,
                      const DayCounter& dayCounter,
                      Compounding compounding,
                      Frequency frequency,
                      bool includeSettlementDateFlows,
                      Date settlementDate = Date(),
                      Date npvDate = Date(),
                      Real accuracy = 1.0e-10,
                      Size maxIterations = 10000,
                      Rate guess = 0.05);
    static Time duration(const Leg&,
                         const InterestRate&,
                         Duration::Type type,
                         bool includeSettlementDateFlows,
                         Date settlementDate = Date());

    static Time duration(const Leg&,
             Rate yield,
             const DayCounter& dayCounter,
             Compounding compounding,
             Frequency frequency,
             Duration::Type type,
             bool includeSettlementDateFlows,
             Date settlementDate = Date(),
             Date npvDate = Date());

    static Real convexity(const Leg&,
                          const InterestRate&,
                          bool includeSettlementDateFlows,
                          Date settlementDate = Date(),
                          Date npvDate = Date());

    static Real convexity(const Leg&,
             Rate yield,
             const DayCounter& dayCounter,
             Compounding compounding,
             Frequency frequency,
             bool includeSettlementDateFlows,
             Date settlementDate = Date(),
             Date npvDate = Date());

    static Real basisPointValue(const Leg& leg,
             const InterestRate& yield,
             bool includeSettlementDateFlows,
             Date settlementDate = Date(),
             Date npvDate = Date());

    static Real basisPointValue(const Leg& leg,
             Rate yield,
             const DayCounter& dayCounter,
             Compounding compounding,
             Frequency frequency,
             bool includeSettlementDateFlows,
             Date settlementDate = Date(),
             Date npvDate = Date());

    static Spread zSpread(const Leg& leg,
             Real npv,
             const ext::shared_ptr<YieldTermStructure>&,
             const DayCounter& dayCounter,
             Compounding compounding,
             Frequency frequency,
             bool includeSettlementDateFlows,
             Date settlementDate = Date(),
             Date npvDate = Date(),
             Real accuracy = 1.0e-10,
             Size maxIterations = 100,
             Rate guess = 0.0);

};


#endif
