
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2014 StatPro Italia srl
 Copyright (C) 2018, 2019 Matthias Lungwitz
 Copyright (C) 2020 Marcin Rybacki
 
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

#ifndef quantlib_term_structures_i
#define quantlib_term_structures_i

%include common.i
%include types.i
%include interestrate.i
%include date.i
%include calendars.i
%include daycounters.i
%include currencies.i
%include observer.i
%include marketelements.i
%include interpolation.i
%include functions.i


%{
using QuantLib::TermStructure;
%}

%shared_ptr(TermStructure);
class TermStructure : public Observable {
  private:
    TermStructure();
  public:
    DayCounter dayCounter() const;
    Time timeFromReference(const Date& date) const;
    Calendar calendar() const;
    Date referenceDate() const;
    Date maxDate() const;
    Time maxTime() const;
    // from Extrapolator, since we can't use multiple inheritance
    // and we're already inheriting from Observable
    void enableExtrapolation();
    void disableExtrapolation();
    bool allowsExtrapolation();
};


%{
using QuantLib::YieldTermStructure;
%}

%shared_ptr(YieldTermStructure);
class YieldTermStructure : public TermStructure {
  private:
    YieldTermStructure();
  public:
    DiscountFactor discount(const Date&, bool extrapolate = false);
    DiscountFactor discount(Time, bool extrapolate = false);
    InterestRate zeroRate(const Date& d,
                          const DayCounter&, Compounding, Frequency f = Annual,
                          bool extrapolate = false) const;
    InterestRate zeroRate(Time t,
                          Compounding, Frequency f = Annual,
                          bool extrapolate = false) const;
    InterestRate forwardRate(const Date& d1, const Date& d2,
                             const DayCounter&, Compounding,
                             Frequency f = Annual,
                             bool extrapolate = false) const;
    InterestRate forwardRate(Time t1, Time t2,
                             Compounding, Frequency f = Annual,
                             bool extrapolate = false) const;
};

%template(YieldTermStructureHandle) Handle<YieldTermStructure>;
%template(RelinkableYieldTermStructureHandle) RelinkableHandle<YieldTermStructure>;


// implied term structure

%{
using QuantLib::ImpliedTermStructure;
%}

%shared_ptr(ImpliedTermStructure);
class ImpliedTermStructure : public YieldTermStructure {
  public:
    ImpliedTermStructure(const Handle<YieldTermStructure>& curveHandle,
                         const Date& referenceDate);
};

// spreaded term structures

%{
using QuantLib::ZeroSpreadedTermStructure;
using QuantLib::ForwardSpreadedTermStructure;
%}

%shared_ptr(ZeroSpreadedTermStructure);
class ZeroSpreadedTermStructure : public YieldTermStructure {
  public:
    ZeroSpreadedTermStructure(const Handle<YieldTermStructure>& curveHandle,
                              const Handle<Quote>& spreadHandle,
                              Compounding comp = QuantLib::Continuous,
                              Frequency freq = QuantLib::NoFrequency,
                              const DayCounter& dc = DayCounter());
};

%shared_ptr(ForwardSpreadedTermStructure);
class ForwardSpreadedTermStructure : public YieldTermStructure {
  public:
    ForwardSpreadedTermStructure(const Handle<YieldTermStructure>& curveHandle,
                                 const Handle<Quote>& spreadHandle);
};

%{
using QuantLib::InterpolatedPiecewiseZeroSpreadedTermStructure;
%}

%shared_ptr(InterpolatedPiecewiseZeroSpreadedTermStructure<Linear>);
%shared_ptr(InterpolatedPiecewiseZeroSpreadedTermStructure<BackwardFlat>);

template <class Interpolator>
class InterpolatedPiecewiseZeroSpreadedTermStructure : public YieldTermStructure {
  public:
    InterpolatedPiecewiseZeroSpreadedTermStructure(
                const Handle<YieldTermStructure>& curveHandle,
                const std::vector< Handle<Quote> >& spreadHandles,
                const std::vector<Date>& dates,
                Compounding comp = QuantLib::Continuous,
                Frequency freq = QuantLib::NoFrequency,
                const DayCounter& dc = DayCounter(),
                const Interpolator& factory = Interpolator());
};

%template(SpreadedLinearZeroInterpolatedTermStructure) InterpolatedPiecewiseZeroSpreadedTermStructure<Linear>;
%template(SpreadedBackwardFlatZeroInterpolatedTermStructure) InterpolatedPiecewiseZeroSpreadedTermStructure<BackwardFlat>;


// flat forward curve

%{
using QuantLib::FlatForward;
%}

%shared_ptr(FlatForward);
class FlatForward : public YieldTermStructure {
  public:
    FlatForward(const Date& referenceDate,
                const Handle<Quote>& forward,
                const DayCounter& dayCounter,
                Compounding compounding = QuantLib::Continuous,
                Frequency frequency = QuantLib::Annual);
    FlatForward(const Date& referenceDate,
                Rate forward,
                const DayCounter& dayCounter,
                Compounding compounding = QuantLib::Continuous,
                Frequency frequency = QuantLib::Annual);
    FlatForward(Integer settlementDays, const Calendar& calendar,
                const Handle<Quote>& forward,
                const DayCounter& dayCounter,
                Compounding compounding = QuantLib::Continuous,
                Frequency frequency = QuantLib::Annual);
    FlatForward(Integer settlementDays, const Calendar& calendar,
                Rate forward,
                const DayCounter& dayCounter,
                Compounding compounding = QuantLib::Continuous,
                Frequency frequency = QuantLib::Annual);
};


%{
using QuantLib::UltimateForwardTermStructure;
%}

%shared_ptr(UltimateForwardTermStructure);
class UltimateForwardTermStructure : public YieldTermStructure {
  public:
    UltimateForwardTermStructure(const Handle<YieldTermStructure>& curveHandle,
                                 const Handle<Quote>& lastLiquidForwardRate,
                                 const Handle<Quote>& ultimateForwardRate,
                                 const Period& firstSmoothingPoint,
                                 Real alpha);
};


#if defined(SWIGPYTHON)
%{
using QuantLib::CompositeZeroYieldStructure;
%}

%shared_ptr(CompositeZeroYieldStructure<BinaryFunction>);

template <class F>
class CompositeZeroYieldStructure : public YieldTermStructure {
  public:
    %extend {
        CompositeZeroYieldStructure(
                const Handle<YieldTermStructure>& h1,
                const Handle<YieldTermStructure>& h2,
                PyObject* function,
                Compounding comp = QuantLib::Continuous,
                Frequency freq = QuantLib::NoFrequency) {
            return new CompositeZeroYieldStructure<F>(h1, h2, F(function), comp, freq);
        }
    }
};

%template(CompositeZeroYieldStructure) CompositeZeroYieldStructure<BinaryFunction>;
#endif


%{
using QuantLib::QuantoTermStructure;
%}

%shared_ptr(QuantoTermStructure);
class QuantoTermStructure : public YieldTermStructure {
  public:
    QuantoTermStructure(const Handle<YieldTermStructure>& underlyingDividendTS,
                        Handle<YieldTermStructure> riskFreeTS,
                        Handle<YieldTermStructure> foreignRiskFreeTS,
                        Handle<BlackVolTermStructure> underlyingBlackVolTS,
                        Real strike,
                        Handle<BlackVolTermStructure> exchRateBlackVolTS,
                        Real exchRateATMlevel,
                        Real underlyingExchRateCorrelation);
};


#endif
