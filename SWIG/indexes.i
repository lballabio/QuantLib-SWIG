
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007 StatPro Italia srl
 Copyright (C) 2015, 2018 Matthias Groncki
 Copyright (C) 2016 Peter Caspers
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

#ifndef quantlib_indexes_i
#define quantlib_indexes_i

%include date.i
%include calendars.i
%include daycounters.i
%include currencies.i
%include types.i
%include termstructures.i
%include timeseries.i
%include vectors.i
%include boost_shared_ptr.i

%{
using QuantLib::IndexManager;
%}

class IndexManager {
    #if defined(SWIGRUBY)
    %rename("hasHistory?")  hasHistory;
    #endif
  private:
    IndexManager();
  public:
    static IndexManager& instance();
    void setHistory(const std::string& name, const TimeSeries<Real>& fixings);
    const TimeSeries<Real>& getHistory(const std::string& name) const;
    bool hasHistory(const std::string& name) const;
    std::vector<std::string> histories() const;
    void clearHistory(const std::string& name);
    void clearHistories();
};


// base index class

%{
using QuantLib::Index;
%}

%shared_ptr(Index)

class Index : public Observable {
    #if defined(SWIGRUBY)
    %rename("isValidFixingDate?") isValidFixingDate;
    %rename("addFixing!") addFixing;
    #endif
  private:
    Index();
  public:
    std::string name() const;
    Calendar fixingCalendar() const;
    bool isValidFixingDate(const Date& fixingDate) const;
    Real fixing(const Date& fixingDate,
                bool forecastTodaysFixing = false) const;
    void addFixing(const Date& fixingDate, Rate fixing,
                   bool forceOverwrite = false);
    const TimeSeries<Real>& timeSeries() const;
    void clearFixings();
    %extend {
        #if defined(SWIGRUBY)
        %rename("addFixings!") addFixings;
        #endif
        void addFixings(const std::vector<Date>& fixingDates,
                        const std::vector<Rate>& fixings,
                        bool forceOverwrite = false) {
            self->addFixings(fixingDates.begin(),fixingDates.end(),
                             fixings.begin(),
                             forceOverwrite);
        }
        std::string __str__() {
            return self->name()+" index";
        }
    }
};


// interest-rate indexes
%{
using QuantLib::InterestRateIndex;
%}

%shared_ptr(InterestRateIndex)

class InterestRateIndex : public Index {
  protected:
    InterestRateIndex();
  public:
    std::string familyName() const;
    Period tenor() const;
    Natural fixingDays() const;
    Date fixingDate(const Date& valueDate) const;
    Currency currency() const;
    DayCounter dayCounter() const;
    Date maturityDate(const Date& valueDate) const;
    Date valueDate(const Date& fixingDate) const;
};


// IborIndex indexes
%{
using QuantLib::IborIndex;
using QuantLib::OvernightIndex;
%}

%shared_ptr(IborIndex)

class IborIndex : public InterestRateIndex {
    #if defined(SWIGRUBY)
    %rename("isAdjusted?") isAdjusted;
    #endif
  public:
    IborIndex(const std::string& familyName,
              const Period& tenor,
              Integer settlementDays,
              const Currency& currency,
              const Calendar& calendar,
              BusinessDayConvention convention,
              bool endOfMonth,
              const DayCounter& dayCounter,
              const Handle<YieldTermStructure>& h =
                                    Handle<YieldTermStructure>());
    BusinessDayConvention businessDayConvention() const;
    bool endOfMonth() const;
    Handle<YieldTermStructure> forwardingTermStructure() const;
    boost::shared_ptr<IborIndex> clone(
                                   const Handle<YieldTermStructure>&) const;
};

%inline %{
    boost::shared_ptr<IborIndex> as_iborindex(
                          const boost::shared_ptr<InterestRateIndex>& index) {
        return boost::dynamic_pointer_cast<IborIndex>(index);
    }
%}

%shared_ptr(OvernightIndex)

class OvernightIndex : public IborIndex {
  public:
    OvernightIndex(const std::string& familyName,
                   Integer settlementDays,
                   const Currency& currency,
                   const Calendar& calendar,
                   const DayCounter& dayCounter,
                   const Handle<YieldTermStructure>& h =
                                    Handle<YieldTermStructure>());
};

%{
using QuantLib::Libor;
%}

%shared_ptr(Libor)

class Libor : public IborIndex {
  public:
    Libor(const std::string& familyName,
          const Period& tenor,
          Natural settlementDays,
          const Currency& currency,
          const Calendar& financialCenterCalendar,
          const DayCounter& dayCounter,
          const Handle<YieldTermStructure>& h =
                                     Handle<YieldTermStructure>());
};


%define export_xibor_instance(Name)
%{
using QuantLib::Name;
%}
%shared_ptr(Name)

class Name : public IborIndex {
  public:
      Name(const Period& tenor,
           const Handle<YieldTermStructure>& h =
                                    Handle<YieldTermStructure>());
};
%enddef

%define export_quoted_xibor_instance(Name,Base)
%{
using QuantLib::Name;
%}
%shared_ptr(Name)

class Name : public Base {
  public:
      Name(const Handle<YieldTermStructure>& h =
                                    Handle<YieldTermStructure>());
};
%enddef

%define export_overnight_instance(Name)
%{
using QuantLib::Name;
%}
%shared_ptr(Name)

class Name : public OvernightIndex {
  public:
      Name(const Handle<YieldTermStructure>& h =
                                    Handle<YieldTermStructure>());
};
%enddef


%{
using QuantLib::SwapIndex;
%}

%shared_ptr(SwapIndex)

class SwapIndex : public InterestRateIndex {
    #if defined(SWIGRUBY)
    %rename("isAdjusted?") isAdjusted;
    #endif
  public:
    SwapIndex(const std::string& familyName,
              const Period& tenor,
              Integer settlementDays,
              const Currency& currency,
              const Calendar& calendar,
              const Period& fixedLegTenor,
              BusinessDayConvention fixedLegConvention,
              const DayCounter& fixedLegDayCounter,
              const boost::shared_ptr<IborIndex>& iborIndex);
    SwapIndex(const std::string& familyName,
              const Period& tenor,
              Integer settlementDays,
              const Currency& currency,
              const Calendar& calendar,
              const Period& fixedLegTenor,
              BusinessDayConvention fixedLegConvention,
              const DayCounter& fixedLegDayCounter,
              const boost::shared_ptr<IborIndex>& iborIndex,
              const Handle<YieldTermStructure>& discountCurve);
    Period fixedLegTenor() const;
    BusinessDayConvention fixedLegConvention() const;
    boost::shared_ptr<IborIndex> iborIndex() const;
    Handle<YieldTermStructure> forwardingTermStructure() const;
    Handle<YieldTermStructure> discountingTermStructure() const;
    boost::shared_ptr<SwapIndex> clone(const Handle<YieldTermStructure>& h) const;
    boost::shared_ptr<SwapIndex> clone(const Handle<YieldTermStructure>& forwarding,
                                       const Handle<YieldTermStructure>& discounting) const;
    boost::shared_ptr<SwapIndex> clone(const Period& tenor) const;
};

%define export_swap_instance(Name)
%{
using QuantLib::Name;
%}
%shared_ptr(Name)
class Name : public SwapIndex {
  public:
    Name(const Period &tenor,
         const Handle<YieldTermStructure>& h =
                                    Handle<YieldTermStructure>());
    Name(const Period &tenor,
         const Handle<YieldTermStructure>& h1,
         const Handle<YieldTermStructure>& h2);
};
%enddef

%define export_quoted_swap_instance(Name,Base)
%{
using QuantLib::Name;
%}
%shared_ptr(Name)
class Name : public Base {
  public:
    Name(const Handle<YieldTermStructure>& h =
                                    Handle<YieldTermStructure>());
    Name(const Handle<YieldTermStructure>& h1,
         const Handle<YieldTermStructure>& h2);
};
%enddef

%inline %{
    boost::shared_ptr<SwapIndex> as_swap_index(
                          const boost::shared_ptr<InterestRateIndex>& index) {
        return boost::dynamic_pointer_cast<SwapIndex>(index);
    }
%}

%{
using QuantLib::SwapSpreadIndex;
%}

%shared_ptr(SwapSpreadIndex)

class SwapSpreadIndex : public InterestRateIndex {
  public:
    SwapSpreadIndex(const std::string& familyName,
                    const boost::shared_ptr<SwapIndex>& swapIndex1,
                    const boost::shared_ptr<SwapIndex>& swapIndex2,
                    const Real gearing1 = 1.0,
                    const Real gearing2 = -1.0);
    Rate forecastFixing(const Date& fixingDate) const;
    Rate pastFixing(const Date& fixingDate) const;
    boost::shared_ptr<SwapIndex> swapIndex1();
    boost::shared_ptr<SwapIndex> swapIndex2();
    Real gearing1();
    Real gearing2();
};

export_xibor_instance(AUDLibor);
export_xibor_instance(CADLibor);
export_xibor_instance(Cdor);
export_xibor_instance(CHFLibor);
export_xibor_instance(DKKLibor);

export_xibor_instance(Bbsw);
export_quoted_xibor_instance(Bbsw1M,Bbsw);
export_quoted_xibor_instance(Bbsw2M,Bbsw);
export_quoted_xibor_instance(Bbsw3M,Bbsw);
export_quoted_xibor_instance(Bbsw4M,Bbsw);
export_quoted_xibor_instance(Bbsw5M,Bbsw);
export_quoted_xibor_instance(Bbsw6M,Bbsw);

export_xibor_instance(Bkbm);
export_quoted_xibor_instance(Bkbm1M,Bkbm);
export_quoted_xibor_instance(Bkbm2M,Bkbm);
export_quoted_xibor_instance(Bkbm3M,Bkbm);
export_quoted_xibor_instance(Bkbm4M,Bkbm);
export_quoted_xibor_instance(Bkbm5M,Bkbm);
export_quoted_xibor_instance(Bkbm6M,Bkbm);

export_xibor_instance(Euribor);
export_quoted_xibor_instance(EuriborSW,Euribor);
export_quoted_xibor_instance(Euribor2W,Euribor);
export_quoted_xibor_instance(Euribor3W,Euribor);
export_quoted_xibor_instance(Euribor1M,Euribor);
export_quoted_xibor_instance(Euribor2M,Euribor);
export_quoted_xibor_instance(Euribor3M,Euribor);
export_quoted_xibor_instance(Euribor4M,Euribor);
export_quoted_xibor_instance(Euribor5M,Euribor);
export_quoted_xibor_instance(Euribor6M,Euribor);
export_quoted_xibor_instance(Euribor7M,Euribor);
export_quoted_xibor_instance(Euribor8M,Euribor);
export_quoted_xibor_instance(Euribor9M,Euribor);
export_quoted_xibor_instance(Euribor10M,Euribor);
export_quoted_xibor_instance(Euribor11M,Euribor);
export_quoted_xibor_instance(Euribor1Y,Euribor);

export_xibor_instance(Euribor365);
export_quoted_xibor_instance(Euribor365_SW,Euribor365);
export_quoted_xibor_instance(Euribor365_2W,Euribor365);
export_quoted_xibor_instance(Euribor365_3W,Euribor365);
export_quoted_xibor_instance(Euribor365_1M,Euribor365);
export_quoted_xibor_instance(Euribor365_2M,Euribor365);
export_quoted_xibor_instance(Euribor365_3M,Euribor365);
export_quoted_xibor_instance(Euribor365_4M,Euribor365);
export_quoted_xibor_instance(Euribor365_5M,Euribor365);
export_quoted_xibor_instance(Euribor365_6M,Euribor365);
export_quoted_xibor_instance(Euribor365_7M,Euribor365);
export_quoted_xibor_instance(Euribor365_8M,Euribor365);
export_quoted_xibor_instance(Euribor365_9M,Euribor365);
export_quoted_xibor_instance(Euribor365_10M,Euribor365);
export_quoted_xibor_instance(Euribor365_11M,Euribor365);
export_quoted_xibor_instance(Euribor365_1Y,Euribor365);

export_xibor_instance(EURLibor);
export_quoted_xibor_instance(EURLiborSW,EURLibor);
export_quoted_xibor_instance(EURLibor2W,EURLibor);
export_quoted_xibor_instance(EURLibor1M,EURLibor);
export_quoted_xibor_instance(EURLibor2M,EURLibor);
export_quoted_xibor_instance(EURLibor3M,EURLibor);
export_quoted_xibor_instance(EURLibor4M,EURLibor);
export_quoted_xibor_instance(EURLibor5M,EURLibor);
export_quoted_xibor_instance(EURLibor6M,EURLibor);
export_quoted_xibor_instance(EURLibor7M,EURLibor);
export_quoted_xibor_instance(EURLibor8M,EURLibor);
export_quoted_xibor_instance(EURLibor9M,EURLibor);
export_quoted_xibor_instance(EURLibor10M,EURLibor);
export_quoted_xibor_instance(EURLibor11M,EURLibor);
export_quoted_xibor_instance(EURLibor1Y,EURLibor);

export_xibor_instance(GBPLibor);
export_xibor_instance(Jibar);
export_xibor_instance(JPYLibor);
export_xibor_instance(Mosprime);
export_xibor_instance(NZDLibor);
export_xibor_instance(Pribor);
export_xibor_instance(Robor);
export_xibor_instance(SEKLibor);
export_xibor_instance(Shibor);
export_xibor_instance(Tibor);
export_xibor_instance(THBFIX);
export_xibor_instance(TRLibor);
export_xibor_instance(USDLibor);
export_xibor_instance(Wibor);
export_xibor_instance(Zibor);

export_overnight_instance(Aonia);
export_overnight_instance(Eonia);
export_overnight_instance(Sonia);
export_overnight_instance(FedFunds);
export_overnight_instance(Nzocr);
export_overnight_instance(Sofr);

export_swap_instance(EuriborSwapIsdaFixA);
export_swap_instance(EuriborSwapIsdaFixB);
export_swap_instance(EuriborSwapIfrFix);

export_swap_instance(EurLiborSwapIsdaFixA);
export_swap_instance(EurLiborSwapIsdaFixB);
export_swap_instance(EurLiborSwapIfrFix);

export_swap_instance(ChfLiborSwapIsdaFix);
export_swap_instance(GbpLiborSwapIsdaFix);
export_swap_instance(JpyLiborSwapIsdaFixAm);
export_swap_instance(JpyLiborSwapIsdaFixPm);
export_swap_instance(UsdLiborSwapIsdaFixAm);
export_swap_instance(UsdLiborSwapIsdaFixPm);

export_xibor_instance(Bibor);
export_quoted_xibor_instance(BiborSW,Bibor);
export_quoted_xibor_instance(Bibor1M,Bibor);
export_quoted_xibor_instance(Bibor2M,Bibor);
export_quoted_xibor_instance(Bibor3M,Bibor);
export_quoted_xibor_instance(Bibor6M,Bibor);
export_quoted_xibor_instance(Bibor9M,Bibor);
export_quoted_xibor_instance(Bibor1Y,Bibor);

#endif
