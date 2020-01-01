/*
 Copyright (C) 2010 Joseph Wang
 Copyright (C) 2010, 2011, 2014 StatPro Italia srl
 Copyright (C) 2018, 2019 Matthias Lungwitz
 
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
%include swap.i

%{
  using QuantLib::Seasonality;
  using QuantLib::MultiplicativePriceSeasonality;
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
    virtual bool indexIsInterpolated() const;
    virtual Rate baseRate() const;
    virtual Handle<YieldTermStructure> nominalTermStructure() const;
    virtual Date baseDate() const;
    void setSeasonality(const boost::shared_ptr<Seasonality>& seasonality =
                                            boost::shared_ptr<Seasonality>());
    boost::shared_ptr<Seasonality> seasonality() const;
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
%#include <ql/indexes/inflationindex.hpp>
%#include <ql/indexes/inflation/all.hpp>
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
                         bool interpolated,
                         Frequency frequency,
                         const Period& availabilityLag,
                         const Currency& currency,
                         const Handle<ZeroInflationTermStructure>& h =
                                       Handle<ZeroInflationTermStructure>());
};

%shared_ptr(YoYInflationIndex)

class YoYInflationIndex : public InflationIndex {
  protected:
    YoYInflationIndex();
};

%define export_zii_instance(Name)
%{
using QuantLib::Name;
%}
%shared_ptr(Name)
class Name : public ZeroInflationIndex {
  public:
    Name(bool interpolated,
         const Handle<ZeroInflationTermStructure>& h =
                                    Handle<ZeroInflationTermStructure>());
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
         const Handle<YoYInflationTermStructure>& h =
                                    Handle<YoYInflationTermStructure>());
};
%enddef

export_zii_instance(EUHICP);
export_zii_instance(EUHICPXT);
export_zii_instance(FRHICP);
export_zii_instance(UKRPI);
export_zii_instance(USCPI);
export_zii_instance(ZACPI);

export_yii_instance(YYEUHICP);
export_yii_instance(YYEUHICPXT);
export_yii_instance(YYFRHICP);
export_yii_instance(YYUKRPI);
export_yii_instance(YYUSCPI);
export_yii_instance(YYZACPI);

// utilities

%{
    using QuantLib::CPI;
%}

struct CPI {
    enum InterpolationType { AsIndex, Flat, Linear };
};

// bootstrapped curves


%{

using QuantLib::BootstrapHelper;
using QuantLib::ZeroCouponInflationSwapHelper;
using QuantLib::YearOnYearInflationSwapHelper;
%}

%shared_ptr(BootstrapHelper<ZeroInflationTermStructure>)
%shared_ptr(BootstrapHelper<YoYInflationTermStructure>)

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

 
#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( boost::shared_ptr<BootstrapHelper<ZeroInflationTermStructure> > )
SWIG_STD_VECTOR_ENHANCED( boost::shared_ptr<BootstrapHelper<YoYInflationTermStructure> > )
#endif
namespace std {
    %template(ZeroHelperVector) vector<boost::shared_ptr<BootstrapHelper<ZeroInflationTermStructure> > >;
    %template(YoYHelperVector) vector<boost::shared_ptr<BootstrapHelper<YoYInflationTermStructure> > >;
}

%shared_ptr(ZeroCouponInflationSwapHelper)
class ZeroCouponInflationSwapHelper : public BootstrapHelper<ZeroInflationTermStructure> {
  public:
    // using extend to prevent deprecation warning
    %extend {
        ZeroCouponInflationSwapHelper(
            const Handle<Quote>& quote,
            const Period& lag,   // lag on swap observation of index
            const Date& maturity,
            const Calendar& calendar,
            BusinessDayConvention bcd,
            const DayCounter& dayCounter,
            const boost::shared_ptr<ZeroInflationIndex>& index,
            const Handle<YieldTermStructure>& nominalTS = Handle<YieldTermStructure>()) {

            return new ZeroCouponInflationSwapHelper(quote,lag,maturity,
                                                     calendar,bcd,
                                                     dayCounter,index,
                                                     nominalTS);
        }
    }
};

%shared_ptr(YearOnYearInflationSwapHelper)
class YearOnYearInflationSwapHelper : public BootstrapHelper<YoYInflationTermStructure> {
  public:
    // using extend to prevent deprecation warning
    %extend {
        YearOnYearInflationSwapHelper(const Handle<Quote>& quote,
                                      const Period& lag,
                                      const Date& maturity,
                                      const Calendar& calendar,
                                      BusinessDayConvention bdc,
                                      const DayCounter& dayCounter,
                                      const boost::shared_ptr<YoYInflationIndex>& index,
                                      const Handle<YieldTermStructure>& nominalTS =
                                                        Handle<YieldTermStructure>()) {
            return new YearOnYearInflationSwapHelper(quote,lag,maturity,
                                                     calendar,bdc,
                                                     dayCounter,index,
                                                     nominalTS);
        }
    }
};


%{
using QuantLib::PiecewiseZeroInflationCurve;
using QuantLib::PiecewiseYoYInflationCurve;
%}

%shared_ptr(PiecewiseZeroInflationCurve<Linear>);

template <class Interpolator>
class PiecewiseZeroInflationCurve : public ZeroInflationTermStructure {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") PiecewiseZeroInflationCurve;
    #endif
  public:
    PiecewiseZeroInflationCurve(
              const Date& referenceDate,
              const Calendar& calendar,
              const DayCounter& dayCounter,
              const Period& lag,
              Frequency frequency,
              bool indexIsInterpolated,
              Rate baseRate,
              const Handle<YieldTermStructure>& nominalTS,
              const std::vector<boost::shared_ptr<BootstrapHelper<ZeroInflationTermStructure> > >& instruments,
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
              const Handle<YieldTermStructure>& nominalTS,
              const std::vector<boost::shared_ptr<BootstrapHelper<YoYInflationTermStructure> > >& instruments,
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


// inflation instruments

%{
using QuantLib::ZeroCouponInflationSwap;
using QuantLib::YearOnYearInflationSwap;
using QuantLib::CPISwap;
%}

%shared_ptr(ZeroCouponInflationSwap)
class ZeroCouponInflationSwap : public Swap {
  public:
    enum Type { Receiver = -1, Payer = 1 };
    ZeroCouponInflationSwap(
                   ZeroCouponInflationSwap::Type type,
                   Real nominal,
                   const Date& start,
                   const Date& maturity,
                   const Calendar& calendar,
                   BusinessDayConvention convention,
                   const DayCounter& dayCounter,
                   Rate fixedRate,
                   const boost::shared_ptr<ZeroInflationIndex>& index,
                   const Period& lag,
                   bool adjustInfObsDates = false,
                   Calendar infCalendar = Calendar(),
                   BusinessDayConvention infConvention = Following);
    Rate fairRate();
    Real fixedLegNPV();
    Real inflationLegNPV();
    std::vector<boost::shared_ptr<CashFlow> > fixedLeg();
    std::vector<boost::shared_ptr<CashFlow> > inflationLeg();
    ZeroCouponInflationSwap::Type type();
};

%shared_ptr(YearOnYearInflationSwap)
class YearOnYearInflationSwap : public Swap {
  public:
    enum Type { Receiver = -1, Payer = 1 };
    YearOnYearInflationSwap(
               YearOnYearInflationSwap::Type type,
               Real nominal,
               const Schedule& fixedSchedule,
               Rate fixedRate,
               const DayCounter& fixedDayCounter,
               const Schedule& yoySchedule,
               const boost::shared_ptr<YoYInflationIndex>& index,
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
    enum Type { Receiver = -1, Payer = 1 };
    CPISwap(
            CPISwap::Type type,
            Real nominal,
            bool subtractInflationNominal,
            Spread spread,
            const DayCounter& floatDayCount,
            const Schedule& floatSchedule,
            const BusinessDayConvention& floatRoll,
            Natural fixingDays,
            const boost::shared_ptr<IborIndex>& floatIndex,
            Rate fixedRate,
            Real baseCPI,
            const DayCounter& fixedDayCount,
            const Schedule& fixedSchedule,
            const BusinessDayConvention& fixedRoll,
            const Period& observationLag,
            const boost::shared_ptr<ZeroInflationIndex>& fixedIndex,
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
     Volatility impliedVolatility(
                           Real price,
                           const Handle<YoYInflationTermStructure>& curve,
                           Volatility guess,
                           Real accuracy = 1.0e-4,
                           Size maxEvaluations = 100,
                           Volatility minVol = 1.0e-7,
                           Volatility maxVol = 4.0) const;
  private:
    YoYInflationCapFloor();
};

%shared_ptr(YoYInflationCap)
class YoYInflationCap : public YoYInflationCapFloor {
  public:
    YoYInflationCap(
            const std::vector<boost::shared_ptr<CashFlow> >& leg,
            const std::vector<Rate>& capRates);
};

%shared_ptr(YoYInflationFloor)
class YoYInflationFloor : public YoYInflationCapFloor {
  public:
    YoYInflationFloor(
            const std::vector<boost::shared_ptr<CashFlow> >& leg,
            const std::vector<Rate>& floorRates);
};

%shared_ptr(YoYInflationCollar)
class YoYInflationCollar : public YoYInflationCapFloor {
  public:
    YoYInflationCollar(
            const std::vector<boost::shared_ptr<CashFlow> >& leg,
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
  public:
    InterpolatedZeroInflationCurve(const Date& referenceDate,
                                   const Calendar& calendar,
                                   const DayCounter& dayCounter,
                                   const Period& lag,
                                   Frequency frequency,
                                   bool indexIsInterpolated,
                                   const Handle<YieldTermStructure>& yTS,
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

template <class Interpolator>
class InterpolatedYoYInflationCurve : public YoYInflationTermStructure {
  public:
    InterpolatedYoYInflationCurve(const Date& referenceDate,
                                   const Calendar& calendar,
                                   const DayCounter& dayCounter,
                                   const Period& lag,
                                   Frequency frequency,
                                   bool indexIsInterpolated,
                                   const Handle<YieldTermStructure>& yTS,
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

#endif
