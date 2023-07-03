/*
 Copyright (C) 2010 Joseph Wang
 Copyright (C) 2010, 2011, 2014 StatPro Italia srl
 Copyright (C) 2018, 2019, 2020 Matthias Lungwitz
 
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

#ifndef quantlib_inflation_i
#define quantlib_inflation_i

%include termstructures.i
%include cashflows.i
%include swap.i
%include interpolation.i

%{
  using QuantLib::InflationTermStructure;
  using QuantLib::Seasonality;
  using QuantLib::MultiplicativePriceSeasonality;
  using QuantLib::KerkhofSeasonality;
%}

%shared_ptr(Seasonality);
class Seasonality {
  private:
    Seasonality();
  public:
    virtual Rate correctZeroRate(const Date &d, const Rate r,
                                 const InflationTermStructure& iTS) const;
    virtual Rate correctYoYRate(const Date &d, const Rate r,
                                const InflationTermStructure& iTS) const;
    virtual bool isConsistent(const InflationTermStructure& iTS);
};

%shared_ptr(MultiplicativePriceSeasonality)
class MultiplicativePriceSeasonality : public Seasonality {
  public:
    MultiplicativePriceSeasonality(const Date& seasonalityBaseDate,
                                   Frequency frequency,
                                   const std::vector<Rate>& seasonalityFactors);
};

%shared_ptr(KerkhofSeasonality)
class KerkhofSeasonality : public Seasonality {
  public:
    KerkhofSeasonality(const Date& seasonalityBaseDate,
                       const std::vector<Rate>& seasonalityFactors);
};


%{
  using QuantLib::InflationTermStructure;
  using QuantLib::YoYInflationTermStructure;
  using QuantLib::ZeroInflationTermStructure;
%}

%shared_ptr(InflationTermStructure);
class InflationTermStructure : public TermStructure {
  private:
    InflationTermStructure();
  public:
    virtual Period observationLag() const;
    virtual Frequency frequency() const;
    virtual Rate baseRate() const;
    virtual Date baseDate() const;
    void setSeasonality(const ext::shared_ptr<Seasonality>& seasonality =
                                            ext::shared_ptr<Seasonality>());
    ext::shared_ptr<Seasonality> seasonality() const;
    bool hasSeasonality() const;
};

%shared_ptr(YoYInflationTermStructure);
class YoYInflationTermStructure : public InflationTermStructure {
  private:
    YoYInflationTermStructure();
  public:
    Rate yoyRate(const Date &d, const Period& instObsLag = Period(-1,Days),
                 bool forceLinearInterpolation = false,
                 bool extrapolate = false) const;
    Rate yoyRate(Time t,
                 bool extrapolate = false) const;
    bool indexIsInterpolated() const;
};

%template(YoYInflationTermStructureHandle) Handle<YoYInflationTermStructure>;
%template(RelinkableYoYInflationTermStructureHandle)
    RelinkableHandle<YoYInflationTermStructure>;


%shared_ptr(ZeroInflationTermStructure);
class ZeroInflationTermStructure : public InflationTermStructure {
  private:
    ZeroInflationTermStructure();
  public:
    Rate zeroRate(const Date &d, const Period& instObsLag = Period(-1,Days),
                  bool forceLinearInterpolation = false,
                  bool extrapolate = false) const;
    Rate zeroRate(Time t,
                  bool extrapolate = false) const;
};

%template(ZeroInflationTermStructureHandle) Handle<ZeroInflationTermStructure>;
%template(RelinkableZeroInflationTermStructureHandle)
    RelinkableHandle<ZeroInflationTermStructure>;



// inflation indexes

%fragment("zeroinflationindex", "header") {
using QuantLib::Region;
using QuantLib::CustomRegion;
using QuantLib::InflationIndex;
using QuantLib::ZeroInflationIndex;
using QuantLib::YoYInflationIndex;
}
%fragment("zeroinflationindex");

class Region {
  public:
    std::string name() const;
    std::string code() const;
  protected:
    Region();
};

class CustomRegion : public Region {
  public:
    CustomRegion(const std::string& name,
                 const std::string& code);
};

%shared_ptr(InflationIndex)
class InflationIndex : public Index {
  protected:
    InflationIndex();
  public:
    std::string familyName() const;
    Region region() const;
    bool revised() const;
    bool interpolated() const;
    Frequency frequency() const;
    Period availabilityLag() const;
    Currency currency() const;
};

%shared_ptr(ZeroInflationIndex)
class ZeroInflationIndex : public InflationIndex {
  public:
      ZeroInflationIndex(const std::string& familyName,
                         const Region& region,
                         bool revised,
                         Frequency frequency,
                         const Period& availabilityLag,
                         const Currency& currency,
                         const Handle<ZeroInflationTermStructure>& h = {});
      ZeroInflationIndex(const std::string& familyName,
                         const Region& region,
                         bool revised,
                         bool interpolated,
                         Frequency frequency,
                         const Period& availabilityLag,
                         const Currency& currency,
                         const Handle<ZeroInflationTermStructure>& h = {});
      Handle<ZeroInflationTermStructure> zeroInflationTermStructure() const;
      ext::shared_ptr<ZeroInflationIndex> clone(const Handle<ZeroInflationTermStructure>& h) const;
};

%inline %{
    ext::shared_ptr<ZeroInflationIndex> as_zero_inflation_index(const ext::shared_ptr<Index>& i) {
        return ext::dynamic_pointer_cast<ZeroInflationIndex>(i);
    }
%}

%shared_ptr(YoYInflationIndex)
class YoYInflationIndex : public InflationIndex {
  public:
    // Constructor for year-on-year indices defined as a ratio.
    YoYInflationIndex(
            const ext::shared_ptr<ZeroInflationIndex>& underlyingIndex,
            bool interpolated,
            Handle<YoYInflationTermStructure> ts = {});

    // Constructor for quoted year-on-year indices.
    YoYInflationIndex(
            const std::string& familyName,
            const Region& region,
            bool revised,
            bool interpolated,
            Frequency frequency,
            const Period& availabilityLag,
            const Currency& currency,
            Handle<YoYInflationTermStructure> ts = {});

    // Deprecated constructor.
    YoYInflationIndex(const std::string& familyName,
                      const Region& region,
                      bool revised,
                      bool interpolated,
                      bool ratio,
                      Frequency frequency,
                      const Period& availabilityLag,
                      const Currency& currency,
                      const Handle<YoYInflationTermStructure>& ts = {});
    bool ratio() const;
    ext::shared_ptr<ZeroInflationIndex> underlyingIndex() const;
    Handle<YoYInflationTermStructure> yoyInflationTermStructure() const;
    ext::shared_ptr<YoYInflationIndex> clone(const Handle<YoYInflationTermStructure>& h) const;
};

%define export_zii_instance(Name)
%{
using QuantLib::Name;
%}
%shared_ptr(Name)
class Name : public ZeroInflationIndex {
  public:
    Name(const Handle<ZeroInflationTermStructure>& h = {});
    Name(bool interpolated,
         const Handle<ZeroInflationTermStructure>& h = {});
};
%enddef

%define export_yii_instance(Name)
%fragment("Name","header") {
using QuantLib::Name;
}
%fragment("Name");
%shared_ptr(Name)
class Name : public YoYInflationIndex {
  public:
    Name(bool interpolated,
         const Handle<YoYInflationTermStructure>& h = {});
};
%enddef

export_zii_instance(EUHICP);
export_zii_instance(EUHICPXT);
export_zii_instance(FRHICP);
export_zii_instance(UKRPI);
export_zii_instance(USCPI);
export_zii_instance(ZACPI);

%{
using QuantLib::UKHICP;
%}
%shared_ptr(UKHICP)
class UKHICP : public ZeroInflationIndex {
  public:
    UKHICP(const Handle<ZeroInflationTermStructure>& h = {});
};


export_yii_instance(YYEUHICP);
export_yii_instance(YYEUHICPXT);
export_yii_instance(YYEUHICPr);
export_yii_instance(YYFRHICP);
export_yii_instance(YYFRHICPr);
export_yii_instance(YYUKRPI);
export_yii_instance(YYUKRPIr);
export_yii_instance(YYUSCPI);
export_yii_instance(YYUSCPIr);
export_yii_instance(YYZACPI);
export_yii_instance(YYZACPIr);

%{
using QuantLib::AUCPI;
%}
%shared_ptr(AUCPI)
class AUCPI : public ZeroInflationIndex {
  public:
    AUCPI(Frequency frequency,
          bool revised,
          const Handle<ZeroInflationTermStructure>& h = {});
    AUCPI(Frequency frequency,
          bool revised,
          bool interpolated,
          const Handle<ZeroInflationTermStructure>& h = {});
};


// utilities

%{
using QuantLib::CPI;
%}

struct CPI {
    enum InterpolationType { AsIndex, Flat, Linear };

    static Real laggedFixing(const ext::shared_ptr<ZeroInflationIndex>& index,
                             const Date& date,
                             const Period& observationLag,
                             InterpolationType interpolationType);
};

// cashflows

%{
using QuantLib::InflationCoupon;
using QuantLib::CPICoupon;
using QuantLib::CPICouponPricer;
using QuantLib::CPICashFlow;
using QuantLib::ZeroInflationCashFlow;
%}

%shared_ptr(InflationCoupon)
class InflationCoupon : public Coupon {
  private:
    InflationCoupon();
  public:
    Date fixingDate() const;
    Integer fixingDays() const;
    Period observationLag() const;
    Rate indexFixing() const;
    ext::shared_ptr<InflationIndex> index() const;
};

%inline %{
    ext::shared_ptr<InflationCoupon> as_inflation_coupon(
                                      const ext::shared_ptr<CashFlow>& cf) {
        return ext::dynamic_pointer_cast<InflationCoupon>(cf);
    }
%}

%shared_ptr(CPICouponPricer)
class CPICouponPricer {
  public:
    CPICouponPricer();
};

%shared_ptr(CPICoupon)
class CPICoupon : public InflationCoupon {
  public:
    CPICoupon(Real baseCPI,
              const Date& paymentDate,
              Real nominal,
              const Date& startDate,
              const Date& endDate,
              const ext::shared_ptr<ZeroInflationIndex>& index,
              const Period& observationLag,
              CPI::InterpolationType observationInterpolation,
              const DayCounter& dayCounter,
              Real fixedRate,
              const Date& refPeriodStart = Date(),
              const Date& refPeriodEnd = Date(),
              const Date& exCouponDate = Date());
    CPICoupon(const Date& baseDate,
              const Date& paymentDate,
              Real nominal,
              const Date& startDate,
              const Date& endDate,
              const ext::shared_ptr<ZeroInflationIndex>& index,
              const Period& observationLag,
              CPI::InterpolationType observationInterpolation,
              const DayCounter& dayCounter,
              Real fixedRate,
              const Date& refPeriodStart = Date(),
              const Date& refPeriodEnd = Date(),
              const Date& exCouponDate = Date());
    CPICoupon(Real baseCPI,
              const Date& baseDate,
              const Date& paymentDate,
              Real nominal,
              const Date& startDate,
              const Date& endDate,
              const ext::shared_ptr<ZeroInflationIndex>& index,
              const Period& observationLag,
              CPI::InterpolationType observationInterpolation,
              const DayCounter& dayCounter,
              Real fixedRate,
              const Date& refPeriodStart = Date(),
              const Date& refPeriodEnd = Date(),
              const Date& exCouponDate = Date());
    CPICoupon(Real baseCPI,
              const Date& paymentDate,
              Real nominal,
              const Date& startDate,
              const Date& endDate,
              const ext::shared_ptr<ZeroInflationIndex>& index,
              const Period& observationLag,
              CPI::InterpolationType observationInterpolation,
              const DayCounter& dayCounter,
              Real fixedRate,
              Spread spread,
              const Date& refPeriodStart = Date(),
              const Date& refPeriodEnd = Date(),
              const Date& exCouponDate = Date());
    CPICoupon(const Date& baseDate,
              const Date& paymentDate,
              Real nominal,
              const Date& startDate,
              const Date& endDate,
              const ext::shared_ptr<ZeroInflationIndex>& index,
              const Period& observationLag,
              CPI::InterpolationType observationInterpolation,
              const DayCounter& dayCounter,
              Real fixedRate,
              Spread spread,
              const Date& refPeriodStart = Date(),
              const Date& refPeriodEnd = Date(),
              const Date& exCouponDate = Date());
    CPICoupon(Real baseCPI,
              const Date& baseDate,
              const Date& paymentDate,
              Real nominal,
              const Date& startDate,
              const Date& endDate,
              const ext::shared_ptr<ZeroInflationIndex>& index,
              const Period& observationLag,
              CPI::InterpolationType observationInterpolation,
              const DayCounter& dayCounter,
              Real fixedRate,
              Spread spread,
              const Date& refPeriodStart = Date(),
              const Date& refPeriodEnd = Date(),
              const Date& exCouponDate = Date());
    Rate fixedRate() const;
    Spread spread() const;
    Rate adjustedIndexGrowth() const;
    Rate indexFixing() const;
    Real indexRatio(Date d) const;
    Rate baseCPI() const;
    Date baseDate() const;
    CPI::InterpolationType observationInterpolation() const;
    ext::shared_ptr<ZeroInflationIndex> cpiIndex() const;
    void setPricer(const ext::shared_ptr<CPICouponPricer>&);
};

%inline %{
    ext::shared_ptr<CPICoupon> as_cpi_coupon(
                                      const ext::shared_ptr<CashFlow>& cf) {
        return ext::dynamic_pointer_cast<CPICoupon>(cf);
    }
%}

%shared_ptr(CPICashFlow)
class CPICashFlow : public IndexedCashFlow {
  public:
    CPICashFlow(Real notional,
                const ext::shared_ptr<ZeroInflationIndex>& index,
                const Date& baseDate,
                Real baseFixing,
                const Date& observationDate,
                const Period& observationLag,
                CPI::InterpolationType interpolation,
                const Date& paymentDate,
                bool growthOnly = false);
    CPI::InterpolationType interpolation() const;
    Frequency frequency() const;
};

%inline %{
    ext::shared_ptr<CPICashFlow> as_cpi_cashflow(
                                      const ext::shared_ptr<CashFlow>& cf) {
        return ext::dynamic_pointer_cast<CPICashFlow>(cf);
    }
%}

%{
Leg _CPILeg(const std::vector<Real>& nominals,
            const Schedule& schedule,
            const ext::shared_ptr<ZeroInflationIndex>& index,
            Real baseCPI,
            const Period& observationLag,
            const DayCounter& paymentDayCounter = DayCounter(),
            const BusinessDayConvention paymentConvention = Following,
            const std::vector<Real>& fixedRates = std::vector<Real>(),
            const std::vector<Spread>& spreads = std::vector<Spread>(),
            const std::vector<Rate>& caps = std::vector<Rate>(),
            const std::vector<Rate>& floors = std::vector<Rate>(),
            const Period& exCouponPeriod = Period(),
            const Calendar& exCouponCalendar = Calendar(),
            BusinessDayConvention exCouponConvention = Unadjusted,
            bool exCouponEndOfMonth = false,
            const Calendar& paymentCalendar = Calendar(),
            bool growthOnly = true,
            CPI::InterpolationType observationInterpolation = CPI::AsIndex) {
    return QuantLib::CPILeg(schedule, index, baseCPI, observationLag)
        .withNotionals(nominals)
        .withPaymentDayCounter(paymentDayCounter)
        .withPaymentAdjustment(paymentConvention)
        .withPaymentCalendar(paymentCalendar.empty() ? schedule.calendar() : paymentCalendar)
        .withFixedRates(fixedRates)
        .withSpreads(spreads)
        .withCaps(caps)
        .withFloors(floors)
        .withExCouponPeriod(exCouponPeriod,
                            exCouponCalendar,
                            exCouponConvention,
                            exCouponEndOfMonth)
        .withSubtractInflationNominal(growthOnly)
        .withObservationInterpolation(observationInterpolation);
}
%}
#if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
%feature("kwargs") _CPILeg;
#endif
%rename(CPILeg) _CPILeg;
Leg _CPILeg(const std::vector<Real>& nominals,
            const Schedule& schedule,
            const ext::shared_ptr<ZeroInflationIndex>& index,
            Real baseCPI,
            const Period& observationLag,
            const DayCounter& paymentDayCounter = DayCounter(),
            const BusinessDayConvention paymentConvention = Following,
            const std::vector<Real>& fixedRates = std::vector<Real>(),
            const std::vector<Spread>& spreads = std::vector<Spread>(),
            const std::vector<Rate>& caps = std::vector<Rate>(),
            const std::vector<Rate>& floors = std::vector<Rate>(),
            const Period& exCouponPeriod = Period(),
            const Calendar& exCouponCalendar = Calendar(),
            BusinessDayConvention exCouponConvention = Unadjusted,
            bool exCouponEndOfMonth = false,
            const Calendar& paymentCalendar = Calendar(),
            bool growthOnly = true,
            CPI::InterpolationType observationInterpolation = CPI::AsIndex);

%shared_ptr(ZeroInflationCashFlow)
class ZeroInflationCashFlow : public CashFlow {
  public:
    ZeroInflationCashFlow(Real notional,
                          const ext::shared_ptr<ZeroInflationIndex>& index,
                          CPI::InterpolationType observationInterpolation,
                          const Date& startDate,
                          const Date& endDate,
                          const Period& observationLag,
                          const Date& paymentDate,
                          bool growthOnly = false);
    Real notional() const;
    Date baseDate() const;
    Date fixingDate() const;
    bool growthOnly() const;
    CPI::InterpolationType observationInterpolation() const;
    ext::shared_ptr<ZeroInflationIndex> zeroInflationIndex() const;
};

%inline %{
    ext::shared_ptr<ZeroInflationCashFlow> as_zero_inflation_cash_flow(
                                      const ext::shared_ptr<CashFlow>& cf) {
        return ext::dynamic_pointer_cast<ZeroInflationCashFlow>(cf);
    }
%}

// bootstrapped curves
%{
using QuantLib::BootstrapHelper;
using QuantLib::ZeroCouponInflationSwapHelper;
using QuantLib::YearOnYearInflationSwapHelper;
using QuantLib::YoYOptionletHelper;
%}

%shared_ptr(BootstrapHelper<ZeroInflationTermStructure>)
%shared_ptr(BootstrapHelper<YoYInflationTermStructure>)
%shared_ptr(BootstrapHelper<YoYOptionletVolatilitySurface>)

template <class TS>
class BootstrapHelper : public Observable {
  public:
    Handle<Quote> quote() const;
    Date latestDate() const;
	Date earliestDate() const;
	Date maturityDate() const;
	Date latestRelevantDate() const;
	Date pillarDate() const;
	Real impliedQuote() const;
	Real quoteError() const;
  private:
    BootstrapHelper();
};

%template(ZeroHelper) BootstrapHelper<ZeroInflationTermStructure>;
%template(YoYHelper) BootstrapHelper<YoYInflationTermStructure>;
%template(YoYOptionHelper) BootstrapHelper<YoYOptionletVolatilitySurface>;

 
#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<BootstrapHelper<ZeroInflationTermStructure> > )
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<BootstrapHelper<YoYInflationTermStructure> > )
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<BootstrapHelper<YoYOptionletVolatilitySurface> > )
#endif
namespace std {
    %template(ZeroHelperVector) vector<ext::shared_ptr<BootstrapHelper<ZeroInflationTermStructure> > >;
    %template(YoYHelperVector) vector<ext::shared_ptr<BootstrapHelper<YoYInflationTermStructure> > >;
    %template(YoYOptionHelperVector) vector<ext::shared_ptr<BootstrapHelper<YoYOptionletVolatilitySurface> > >;
}

%shared_ptr(ZeroCouponInflationSwapHelper)
class ZeroCouponInflationSwapHelper : public BootstrapHelper<ZeroInflationTermStructure> {
  public:
    ZeroCouponInflationSwapHelper(
            const Handle<Quote>& quote,
            const Period& lag,   // lag on swap observation of index
            const Date& maturity,
            const Calendar& calendar,
            BusinessDayConvention bcd,
            const DayCounter& dayCounter,
            const ext::shared_ptr<ZeroInflationIndex>& index,
            CPI::InterpolationType observationInterpolation,
            const Handle<YieldTermStructure>& nominalTS);
};

%shared_ptr(YearOnYearInflationSwapHelper)
class YearOnYearInflationSwapHelper : public BootstrapHelper<YoYInflationTermStructure> {
  public:
    YearOnYearInflationSwapHelper(const Handle<Quote>& quote,
                                  const Period& lag,
                                  const Date& maturity,
                                  const Calendar& calendar,
                                  BusinessDayConvention bdc,
                                  const DayCounter& dayCounter,
                                  const ext::shared_ptr<YoYInflationIndex>& index,
                                  const Handle<YieldTermStructure>& nominalTS);
};


%{
using QuantLib::PiecewiseZeroInflationCurve;
using QuantLib::PiecewiseYoYInflationCurve;
%}

%shared_ptr(PiecewiseZeroInflationCurve<Linear>);

template <class Interpolator>
class PiecewiseZeroInflationCurve : public ZeroInflationTermStructure {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    //%feature("kwargs") PiecewiseZeroInflationCurve;
    #endif
  public:
    PiecewiseZeroInflationCurve(
              const Date& referenceDate,
              const Calendar& calendar,
              const DayCounter& dayCounter,
              const Period& lag,
              Frequency frequency,
              Rate baseRate,
              const std::vector<ext::shared_ptr<BootstrapHelper<ZeroInflationTermStructure> > >& instruments,
              Real accuracy = 1.0e-12,
              const Interpolator& i = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};

%template(PiecewiseZeroInflation) PiecewiseZeroInflationCurve<Linear>;


%shared_ptr(PiecewiseYoYInflationCurve<Linear>);

template <class Interpolator>
class PiecewiseYoYInflationCurve : public YoYInflationTermStructure {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") PiecewiseYoYInflationCurve;
    #endif
  public:
    PiecewiseYoYInflationCurve(
              const Date& referenceDate,
              const Calendar& calendar,
              const DayCounter& dayCounter,
              const Period& lag,
              Frequency frequency,
              bool indexIsInterpolated,
              Rate baseRate,
              const std::vector<ext::shared_ptr<BootstrapHelper<YoYInflationTermStructure> > >& instruments,
              Real accuracy = 1.0e-12,
              const Interpolator& i = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};

%template(PiecewiseYoYInflation) PiecewiseYoYInflationCurve<Linear>;


// utilities

%{
using QuantLib::inflationPeriod;
using QuantLib::inflationYearFraction;
%}

%template(DatePair) std::pair<Date,Date>;

std::pair<Date,Date> inflationPeriod(const Date & d,
                                     Frequency f);

Time inflationYearFraction(Frequency f,
                           bool indexIsInterpolated,
                           const DayCounter& dayCount,
                           const Date& d1, const Date& d2);


%inline %{

    Date inflationBaseDate(const Date& referenceDate,
                           const Period& observationLag,
                           Frequency frequency,
                           bool indexIsInterpolated) {
        if (indexIsInterpolated) {
            return referenceDate - observationLag;
        } else {
            return QuantLib::inflationPeriod(referenceDate - observationLag,
                                             frequency).first;
        }
    }

%}

// year-on-year inflation coupons

%{
using QuantLib::YoYInflationCoupon;
using QuantLib::CappedFlooredYoYInflationCoupon;
using QuantLib::YoYInflationCouponPricer;
using QuantLib::setCouponPricer;
using QuantLib::BlackYoYInflationCouponPricer;
using QuantLib::UnitDisplacedBlackYoYInflationCouponPricer;
using QuantLib::BachelierYoYInflationCouponPricer;
%}

%shared_ptr(YoYInflationCouponPricer)
class YoYInflationCouponPricer {
  private:
    YoYInflationCouponPricer();
};

void setCouponPricer(const Leg&, const ext::shared_ptr<YoYInflationCouponPricer>&);

%shared_ptr(YoYInflationCoupon)
class YoYInflationCoupon : public InflationCoupon {
  public:
    YoYInflationCoupon(const Date& paymentDate,
                       Real nominal,
                       const Date& startDate,
                       const Date& endDate,
                       Natural fixingDays,
                       const ext::shared_ptr<YoYInflationIndex>& index,
                       const Period& observationLag,
                       const DayCounter& dayCounter,
                       Real gearing = 1.0,
                       Spread spread = 0.0,
                       const Date& refPeriodStart = Date(),
                       const Date& refPeriodEnd = Date());
    Real gearing() const;
    Spread spread() const;
    Rate adjustedFixing() const;
    ext::shared_ptr<YoYInflationIndex> yoyIndex() const;
};

%inline %{
    ext::shared_ptr<YoYInflationCoupon> as_yoy_inflation_coupon(
                                      const ext::shared_ptr<CashFlow>& cf) {
        return ext::dynamic_pointer_cast<YoYInflationCoupon>(cf);
    }
%}

%shared_ptr(CappedFlooredYoYInflationCoupon)
class CappedFlooredYoYInflationCoupon : public YoYInflationCoupon {
  public:
    CappedFlooredYoYInflationCoupon(const Date& paymentDate,
                                    Real nominal,
                                    const Date& startDate,
                                    const Date& endDate,
                                    Natural fixingDays,
                                    const ext::shared_ptr<YoYInflationIndex>& index,
                                    const Period& observationLag,
                                    const DayCounter& dayCounter,
                                    Real gearing = 1.0,
                                    Spread spread = 0.0,
                                    const Rate cap = Null<Rate>(),
                                    const Rate floor = Null<Rate>(),
                                    const Date& refPeriodStart = Date(),
                                    const Date& refPeriodEnd = Date());
    Rate rate() const;
    Rate cap() const;
    Rate floor() const;
    Rate effectiveCap() const;
    Rate effectiveFloor() const;
    Rate underlyingRate() const;
    bool isCapped() const;
    bool isFloored() const;
};

%inline %{
    ext::shared_ptr<CappedFlooredYoYInflationCoupon> as_capped_floored_yoy_inflation_coupon(
                                      const ext::shared_ptr<CashFlow>& cf) {
        return ext::dynamic_pointer_cast<CappedFlooredYoYInflationCoupon>(cf);
    }
%}

%{
Leg _yoyInflationLeg(const Schedule& schedule,
                     const Calendar& calendar,
                     const ext::shared_ptr<YoYInflationIndex>& index,
                     const Period& observationLag,
                     const std::vector<Real>& notionals,
                     const DayCounter& paymentDayCounter,
                     BusinessDayConvention paymentAdjustment = Following,
                     Natural fixingDays = 0,
                     const std::vector<Real>& gearings = std::vector<Real>(),
                     const std::vector<Spread>& spreads = std::vector<Spread>(),
                     const std::vector<Rate>& caps = std::vector<Rate>(),
                     const std::vector<Rate>& floors = std::vector<Rate>()) {
        return QuantLib::yoyInflationLeg(schedule, calendar, index, observationLag)
            .withNotionals(notionals)
            .withPaymentDayCounter(paymentDayCounter)
            .withPaymentAdjustment(paymentAdjustment)
            .withFixingDays(fixingDays)
            .withGearings(gearings)
            .withSpreads(spreads)
            .withCaps(caps)
            .withFloors(floors);
}
%}
%feature("kwargs") _yoyInflationLeg;
%rename(yoyInflationLeg) _yoyInflationLeg;
Leg _yoyInflationLeg(const Schedule& schedule,
                     const Calendar& calendar,
                     const ext::shared_ptr<YoYInflationIndex>& index,
                     const Period& observationLag,
                     const std::vector<Real>& notionals,
                     const DayCounter& paymentDayCounter,
                     BusinessDayConvention paymentAdjustment = Following,
                     Natural fixingDays = 0,
                     const std::vector<Real>& gearings = std::vector<Real>(),
                     const std::vector<Spread>& spreads = std::vector<Spread>(),
                     const std::vector<Rate>& caps = std::vector<Rate>(),
                     const std::vector<Rate>& floors = std::vector<Rate>());

%shared_ptr(BlackYoYInflationCouponPricer)
class BlackYoYInflationCouponPricer : public YoYInflationCouponPricer {
  public:
    BlackYoYInflationCouponPricer(
        const Handle<YoYOptionletVolatilitySurface>& capletVol,
        const Handle<YieldTermStructure>& nominalTermStructure);
};

%shared_ptr(UnitDisplacedBlackYoYInflationCouponPricer)
class UnitDisplacedBlackYoYInflationCouponPricer : public YoYInflationCouponPricer {
  public:
    UnitDisplacedBlackYoYInflationCouponPricer(
        const Handle<YoYOptionletVolatilitySurface>& capletVol,
        const Handle<YieldTermStructure>& nominalTermStructure);
};

%shared_ptr(BachelierYoYInflationCouponPricer)
class BachelierYoYInflationCouponPricer : public YoYInflationCouponPricer {
  public:
    BachelierYoYInflationCouponPricer(
        const Handle<YoYOptionletVolatilitySurface>& capletVol,
        const Handle<YieldTermStructure>& nominalTermStructure);
};


// inflation instruments

%{
using QuantLib::ZeroCouponInflationSwap;
using QuantLib::YearOnYearInflationSwap;
using QuantLib::CPISwap;
%}

%shared_ptr(ZeroCouponInflationSwap)
class ZeroCouponInflationSwap : public Swap {
  public:
    ZeroCouponInflationSwap(
                   Type type,
                   Real nominal,
                   const Date& start,
                   const Date& maturity,
                   const Calendar& calendar,
                   BusinessDayConvention convention,
                   const DayCounter& dayCounter,
                   Rate fixedRate,
                   const ext::shared_ptr<ZeroInflationIndex>& index,
                   const Period& lag,
                   CPI::InterpolationType observationInterpolation,
                   bool adjustInfObsDates = false,
                   Calendar infCalendar = Calendar(),
                   BusinessDayConvention infConvention = BusinessDayConvention());
    Rate fairRate();
    Real fixedLegNPV();
    Real inflationLegNPV();
    std::vector<ext::shared_ptr<CashFlow> > fixedLeg();
    std::vector<ext::shared_ptr<CashFlow> > inflationLeg();
    Type type();
};

%shared_ptr(YearOnYearInflationSwap)
class YearOnYearInflationSwap : public Swap {
  public:
    YearOnYearInflationSwap(
               Type type,
               Real nominal,
               const Schedule& fixedSchedule,
               Rate fixedRate,
               const DayCounter& fixedDayCounter,
               const Schedule& yoySchedule,
               const ext::shared_ptr<YoYInflationIndex>& index,
               const Period& lag,
               Spread spread,
               const DayCounter& yoyDayCounter,
               const Calendar& paymentCalendar,
               BusinessDayConvention paymentConvention = Following);
    Rate fairRate();
    Real fixedLegNPV();
    Real yoyLegNPV();
    Spread fairSpread();
    const Leg& fixedLeg();
    const Leg& yoyLeg();
};

%shared_ptr(CPISwap)
class CPISwap : public Swap {
  public:
    CPISwap(
            Type type,
            Real nominal,
            bool subtractInflationNominal,
            Spread spread,
            const DayCounter& floatDayCount,
            const Schedule& floatSchedule,
            const BusinessDayConvention& floatRoll,
            Natural fixingDays,
            const ext::shared_ptr<IborIndex>& floatIndex,
            Rate fixedRate,
            Real baseCPI,
            const DayCounter& fixedDayCount,
            const Schedule& fixedSchedule,
            const BusinessDayConvention& fixedRoll,
            const Period& observationLag,
            const ext::shared_ptr<ZeroInflationIndex>& fixedIndex,
            CPI::InterpolationType observationInterpolation = CPI::AsIndex,
            Real inflationNominal = Null<Real>() );
    Rate fairRate();
    Real floatLegNPV();
    Spread fairSpread();
    Real fixedLegNPV();
    const Leg& cpiLeg();
    const Leg& floatLeg();
};


%{
using QuantLib::YoYInflationCapFloor;
using QuantLib::YoYInflationCap;
using QuantLib::YoYInflationFloor;
using QuantLib::YoYInflationCollar;
%}

%shared_ptr(YoYInflationCapFloor)
class YoYInflationCapFloor : public Instrument {
  public:
     enum Type { Cap, Floor, Collar };
     YoYInflationCapFloor(YoYInflationCapFloor::Type type,
                          const Leg& yoyLeg,
                          const std::vector<Rate>& strikes);
     Volatility impliedVolatility(
                           Real price,
                           const Handle<YoYInflationTermStructure>& curve,
                           Volatility guess,
                           Real accuracy = 1.0e-4,
                           Size maxEvaluations = 100,
                           Volatility minVol = 1.0e-7,
                           Volatility maxVol = 4.0) const;
     %extend {
         std::vector<Real> optionletPrices() {
             return self->result<std::vector<Real> >("optionletsPrice");
         }
     }
};

%shared_ptr(YoYInflationCap)
class YoYInflationCap : public YoYInflationCapFloor {
  public:
    YoYInflationCap(
            const std::vector<ext::shared_ptr<CashFlow> >& leg,
            const std::vector<Rate>& capRates);
};

%shared_ptr(YoYInflationFloor)
class YoYInflationFloor : public YoYInflationCapFloor {
  public:
    YoYInflationFloor(
            const std::vector<ext::shared_ptr<CashFlow> >& leg,
            const std::vector<Rate>& floorRates);
};

%shared_ptr(YoYInflationCollar)
class YoYInflationCollar : public YoYInflationCapFloor {
  public:
    YoYInflationCollar(
            const std::vector<ext::shared_ptr<CashFlow> >& leg,
            const std::vector<Rate>& capRates,
            const std::vector<Rate>& floorRates);
};

%{
using QuantLib::InterpolatedZeroInflationCurve;
using QuantLib::InterpolatedYoYInflationCurve;
%}

%shared_ptr(InterpolatedZeroInflationCurve<Linear>);
%shared_ptr(InterpolatedYoYInflationCurve<Linear>);

template <class Interpolator>
class InterpolatedZeroInflationCurve : public ZeroInflationTermStructure {
    //%feature("kwargs") InterpolatedZeroInflationCurve;
  public:
    InterpolatedZeroInflationCurve(const Date& referenceDate,
                                   const Calendar& calendar,
                                   const DayCounter& dayCounter,
                                   const Period& lag,
                                   Frequency frequency,
                                   const std::vector<Date>& dates,
                                   const std::vector<Rate>& rates,
                                   const Interpolator &interpolator = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    const std::vector<Real>& data() const;
    const std::vector<Rate>& rates() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Rate> > nodes() const;
    #endif
};

template <class Interpolator>
class InterpolatedYoYInflationCurve : public YoYInflationTermStructure {
    %feature("kwargs") InterpolatedYoYInflationCurve;
  public:
    InterpolatedYoYInflationCurve(const Date& referenceDate,
                                   const Calendar& calendar,
                                   const DayCounter& dayCounter,
                                   const Period& lag,
                                   Frequency frequency,
                                   bool indexIsInterpolated,
                                   const std::vector<Date>& dates,
                                   const std::vector<Rate>& rates,
                                   const Interpolator &interpolator
                                                        = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    const std::vector<Real>& data() const;
    const std::vector<Rate>& rates() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Rate> > nodes() const;
    #endif
};

%template(ZeroInflationCurve) InterpolatedZeroInflationCurve<Linear>;
%template(YoYInflationCurve) InterpolatedYoYInflationCurve<Linear>;

%{
using QuantLib::YoYCapFloorTermPriceSurface;
using QuantLib::InterpolatedYoYCapFloorTermPriceSurface;
%}

%shared_ptr(YoYCapFloorTermPriceSurface)
class YoYCapFloorTermPriceSurface : public InflationTermStructure {
  private:
    YoYCapFloorTermPriceSurface();
  public:
    virtual std::pair<std::vector<Time>, std::vector<Rate> > atmYoYSwapTimeRates() const;
    virtual std::pair<std::vector<Date>, std::vector<Rate> > atmYoYSwapDateRates() const;
    virtual ext::shared_ptr<YoYInflationTermStructure> YoYTS() const;
    ext::shared_ptr<YoYInflationIndex> yoyIndex();
    virtual BusinessDayConvention businessDayConvention() const;
    virtual Natural fixingDays() const;
    virtual Real price(const Date& d, Rate k);
    virtual Real capPrice(const Date& d, Rate k);
    virtual Real floorPrice(const Date& d, Rate k);
    virtual Rate atmYoYSwapRate(const Date &d,
                                bool extrapolate = true);
    virtual Rate atmYoYRate(const Date &d,
                            const Period &obsLag = Period(-1,Days),
                            bool extrapolate = true);

    virtual Real price(const Period& d, Rate k) const;
    virtual Real capPrice(const Period& d, Rate k) const;
    virtual Real floorPrice(const Period& d, Rate k) const;
    virtual Rate atmYoYSwapRate(const Period &d,
                                bool extrapolate = true) const;
    virtual Rate atmYoYRate(const Period &d,
                            const Period &obsLag = Period(-1,Days),
                            bool extrapolate = true) const;

    virtual std::vector<Rate> strikes();
    virtual std::vector<Rate> capStrikes();
    virtual std::vector<Rate> floorStrikes();
    virtual std::vector<Period> maturities();
    virtual Rate minStrike() const;
    virtual Rate maxStrike() const;
    virtual Date minMaturity() const;
    virtual Date maxMaturity() const;

    virtual Date yoyOptionDateFromTenor(const Period& p) const;
};

%define export_yoy_capfloor_termpricesurface(Name,Interpolator2D, Interpolator1D)

%{
typedef InterpolatedYoYCapFloorTermPriceSurface<Interpolator2D, Interpolator1D> Name;
%}

%shared_ptr(Name);
class Name : public YoYCapFloorTermPriceSurface {
  public:
    %extend {
        Name(Natural fixingDays,
          const Period &yyLag,  // observation lag
          const ext::shared_ptr<YoYInflationIndex>& yii,
          Rate baseRate,
          const Handle<YieldTermStructure> &nominal,
          const DayCounter &dc,
          const Calendar &cal,
          const BusinessDayConvention &bdc,
          const std::vector<Rate> &cStrikes,
          const std::vector<Rate> &fStrikes,
          const std::vector<Period> &cfMaturities,
          const Matrix &cPrice,
          const Matrix &fPrice,
          const Interpolator2D &interpolator2d = Interpolator2D(),
          const Interpolator1D &interpolator1d = Interpolator1D()) {
            return new Name(fixingDays, yyLag, yii, baseRate, nominal,
                            dc, cal, bdc, cStrikes, fStrikes, cfMaturities,
                            cPrice, fPrice);
        }
    }
};
%enddef

export_yoy_capfloor_termpricesurface(YoYInflationCapFloorTermPriceSurface,Bicubic,Cubic);


%{
using QuantLib::YoYInflationCapFloorEngine;
using QuantLib::YoYInflationBlackCapFloorEngine;
using QuantLib::YoYInflationUnitDisplacedBlackCapFloorEngine;
using QuantLib::YoYInflationBachelierCapFloorEngine;
%}

%shared_ptr(YoYInflationBlackCapFloorEngine)
class YoYInflationBlackCapFloorEngine : public PricingEngine {
  public:
    YoYInflationBlackCapFloorEngine(const ext::shared_ptr<YoYInflationIndex>&,
                                    const Handle<YoYOptionletVolatilitySurface>& vol,
                                    const Handle<YieldTermStructure>& nominalTermStructure);
};

%shared_ptr(YoYInflationUnitDisplacedBlackCapFloorEngine)
class YoYInflationUnitDisplacedBlackCapFloorEngine : public PricingEngine {
  public:
    YoYInflationUnitDisplacedBlackCapFloorEngine(const ext::shared_ptr<YoYInflationIndex>&,
                                    const Handle<YoYOptionletVolatilitySurface>& vol,
                                    const Handle<YieldTermStructure>& nominalTermStructure);
};

%shared_ptr(YoYInflationBachelierCapFloorEngine)
class YoYInflationBachelierCapFloorEngine : public PricingEngine {
  public:
    YoYInflationBachelierCapFloorEngine(const ext::shared_ptr<YoYInflationIndex>&,
                                    const Handle<YoYOptionletVolatilitySurface>& vol,
                                    const Handle<YieldTermStructure>& nominalTermStructure);
};

%shared_ptr(YoYOptionletHelper)
class YoYOptionletHelper : public BootstrapHelper<YoYOptionletVolatilitySurface> {
  public:
      %extend {
        YoYOptionletHelper(
         const Handle<Quote>& price,
         Real notional,
         YoYInflationCapFloor::Type capFloorType,
         Period &lag,
         const DayCounter& yoyDayCounter,
         const Calendar& paymentCalendar,
         Natural fixingDays,
         const ext::shared_ptr<YoYInflationIndex>& index,
         Rate strike, Size n,
         const ext::shared_ptr<PricingEngine> &pricer) {
            ext::shared_ptr<QuantLib::YoYInflationCapFloorEngine> engine = ext::dynamic_pointer_cast<YoYInflationCapFloorEngine>(pricer);
            return new YoYOptionletHelper(price, notional, capFloorType, lag, yoyDayCounter, paymentCalendar, fixingDays, index, strike, n, engine);
         }
     }
};

%{
using QuantLib::YoYOptionletStripper;
using QuantLib::InterpolatedYoYOptionletStripper;
%}

%shared_ptr(YoYOptionletStripper)
class YoYOptionletStripper {
  private:
    YoYOptionletStripper();
  public:
      %extend {
        virtual void initialize(const ext::shared_ptr<YoYCapFloorTermPriceSurface>& surf,
                                const ext::shared_ptr<PricingEngine>& pricer,
                                Real slope) const {
            ext::shared_ptr<QuantLib::YoYInflationCapFloorEngine> engine = ext::dynamic_pointer_cast<YoYInflationCapFloorEngine>(pricer);
            return (self)->initialize(surf, engine, slope);
        }
    }
    virtual Rate maxStrike() const;
    virtual std::vector<Rate> strikes() const;
    virtual std::pair<std::vector<Rate>, std::vector<Volatility> > slice(const Date &d) const;
};

%shared_ptr(InterpolatedYoYOptionletStripper<Linear>)
template <class Interpolator1D>
class InterpolatedYoYOptionletStripper: public YoYOptionletStripper {
  public:
    InterpolatedYoYOptionletStripper();
};

%template(InterpolatedYoYInflationOptionletStripper) InterpolatedYoYOptionletStripper<Linear>;

%{
using QuantLib::InterpolatedYoYOptionletVolatilityCurve;
using QuantLib::KInterpolatedYoYOptionletVolatilitySurface;
using QuantLib::YoYInflationCapFloorEngine;
%}

%shared_ptr(InterpolatedYoYOptionletVolatilityCurve<Linear>);

template <class Interpolator1D>
class InterpolatedYoYOptionletVolatilityCurve : public YoYOptionletVolatilitySurface {
  public:
    InterpolatedYoYOptionletVolatilityCurve(Natural settlementDays,
                                            const Calendar&,
                                            BusinessDayConvention bdc,
                                            const DayCounter& dc,
                                            const Period &lag,
                                            Frequency frequency,
                                            bool indexIsInterpolated,
                                            const std::vector<Date> &d,
                                            const std::vector<Volatility> &v,
                                            Rate minStrike,
                                            Rate maxStrike,
                                            const Interpolator1D &i =
                                                        Interpolator1D());
};

%template(InterpolatedYoYInflationOptionletVolatilityCurve) InterpolatedYoYOptionletVolatilityCurve<Linear>;

%shared_ptr(KInterpolatedYoYOptionletVolatilitySurface<Linear>);
template <class Interpolator1D>
class KInterpolatedYoYOptionletVolatilitySurface : public YoYOptionletVolatilitySurface {
  public:
      %extend {
        KInterpolatedYoYOptionletVolatilitySurface(
            Natural settlementDays,
            const Calendar& calendar,
            BusinessDayConvention bdc,
            const DayCounter& dc,
            const Period& lag,
            const ext::shared_ptr<YoYCapFloorTermPriceSurface>& capFloorPrices,
            const ext::shared_ptr<PricingEngine>& pricer,
            const ext::shared_ptr<YoYOptionletStripper>& yoyOptionletStripper,
            Real slope,
            const Interpolator1D& interpolator = Interpolator1D()) {
                 ext::shared_ptr<QuantLib::YoYInflationCapFloorEngine> engine = ext::dynamic_pointer_cast<YoYInflationCapFloorEngine>(pricer);
                     return new KInterpolatedYoYOptionletVolatilitySurface<Interpolator1D>(settlementDays,
                                 calendar, bdc, dc, lag, capFloorPrices, engine, yoyOptionletStripper,
                                 slope, interpolator);
        }
     }
     std::pair<std::vector<Rate>, std::vector<Volatility> > Dslice(
                                                 const Date &d) const;
};

%template(KInterpolatedYoYInflationOptionletVolatilitySurface) KInterpolatedYoYOptionletVolatilitySurface<Linear>;

#endif
