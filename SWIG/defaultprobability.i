/*
 Copyright (C) 2008, 2009 StatPro Italia srl
 Copyright (C) 2018 Matthias Lungwitz

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

#ifndef quantlib_default_probability_structures_i
#define quantlib_default_probability_structures_i

%include common.i
%include types.i
%include date.i
%include calendars.i
%include daycounters.i
%include scheduler.i
%include observer.i
%include marketelements.i
%include interpolation.i
%include termstructures.i
%include piecewiseyieldcurve.i
%include bonds.i

%{
using QuantLib::DefaultProbabilityTermStructure;
%}

%shared_ptr(DefaultProbabilityTermStructure);
class DefaultProbabilityTermStructure : public TermStructure {
  private:
    DefaultProbabilityTermStructure();
  public:
    Probability defaultProbability(const Date&, bool extrapolate = false);
    Probability defaultProbability(Time, bool extrapolate = false);
    Probability defaultProbability(const Date&, const Date&,
                                   bool extrapolate = false);
    Probability defaultProbability(Time, Time, bool extrapolate = false);

    Probability survivalProbability(const Date&, bool extrapolate = false);
    Probability survivalProbability(Time, bool extrapolate = false);

    Real defaultDensity(const Date&, bool extrapolate = false);
    Real defaultDensity(Time, bool extrapolate = false);

    Real hazardRate(const Date&, bool extrapolate = false);
    Real hazardRate(Time, bool extrapolate = false);
};


%template(DefaultProbabilityTermStructureHandle)
Handle<DefaultProbabilityTermStructure>;
%template(RelinkableDefaultProbabilityTermStructureHandle)
RelinkableHandle<DefaultProbabilityTermStructure>;


// concrete curves


// flat forward curve

%{
using QuantLib::FlatHazardRate;
%}

%shared_ptr(FlatHazardRate);
class FlatHazardRate : public DefaultProbabilityTermStructure {
  public:
    FlatHazardRate(Integer settlementDays,
                   const Calendar& calendar,
                   const Handle<Quote>& hazardRate,
                   const DayCounter& dayCounter);
    FlatHazardRate(const Date& todaysDate,
                   const Handle<Quote>& hazardRate,
                   const DayCounter& dayCounter);
};


%{
using QuantLib::InterpolatedHazardRateCurve;
%}

// add other instantiations both here and below the class
%shared_ptr(InterpolatedHazardRateCurve<BackwardFlat>);

template <class Interpolator>
class InterpolatedHazardRateCurve : public DefaultProbabilityTermStructure {
  public:
    InterpolatedHazardRateCurve(const std::vector<Date>& dates,
                                const std::vector<Real>& hazardRates,
                                const DayCounter& dayCounter,
                                const Calendar& calendar = Calendar(),
                                const Interpolator& i = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Real>& hazardRates() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};

%template(HazardRateCurve) InterpolatedHazardRateCurve<BackwardFlat>;


%{
using QuantLib::InterpolatedDefaultDensityCurve;
%}

// add other instantiations both here and below the class
%shared_ptr(InterpolatedDefaultDensityCurve<Linear>);

template <class Interpolator>
class InterpolatedDefaultDensityCurve : public DefaultProbabilityTermStructure {
  public:
    InterpolatedDefaultDensityCurve(const std::vector<Date>& dates,
                                    const std::vector<Real>& densities,
                                    const DayCounter& dayCounter,
                                    const Calendar& calendar = Calendar(),
                                    const Interpolator& i = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Real>& defaultDensities() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};

%template(DefaultDensityCurve) InterpolatedDefaultDensityCurve<Linear>;


%{
using QuantLib::InterpolatedSurvivalProbabilityCurve;
%}

// add other instantiations both here and below the class
%shared_ptr(InterpolatedSurvivalProbabilityCurve<Linear>);

template <class Interpolator>
class InterpolatedSurvivalProbabilityCurve : public DefaultProbabilityTermStructure {
  public:
    InterpolatedSurvivalProbabilityCurve(const std::vector<Date>& dates,
                                         const std::vector<Probability>& probabilities,
                                         const DayCounter& dayCounter,
                                         const Calendar& calendar = Calendar(),
                                         const Interpolator& i = Interpolator());
    const std::vector<Date>& dates() const;
    const std::vector<Probability>& survivalProbabilities() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};

%template(SurvivalProbabilityCurve) InterpolatedSurvivalProbabilityCurve<Linear>;


%{
using QuantLib::DefaultProbabilityHelper;
using QuantLib::SpreadCdsHelper;
using QuantLib::UpfrontCdsHelper;
%}

// rate helpers for curve bootstrapping

%shared_ptr(DefaultProbabilityHelper)
class DefaultProbabilityHelper : public Observable {
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
    DefaultProbabilityHelper();
};

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<DefaultProbabilityHelper> )
#endif
namespace std {
    %template(DefaultProbabilityHelperVector)
    vector<ext::shared_ptr<DefaultProbabilityHelper> >;
}


%shared_ptr(SpreadCdsHelper)
class SpreadCdsHelper : public DefaultProbabilityHelper {
  public:
    SpreadCdsHelper(
            const Handle<Quote>& spread,
            const Period& tenor,
            Integer settlementDays,
            const Calendar& calendar,
            Frequency frequency,
            BusinessDayConvention convention,
            DateGeneration::Rule rule,
            const DayCounter& dayCounter,
            Real recoveryRate,
            const Handle<YieldTermStructure>& discountCurve,
            bool settlesAccrual = true,
            bool paysAtDefaultTime = true,
            const Date& startDate = Date(),
            const DayCounter& lastPeriodDayCounter = DayCounter(),
            bool rebatesAccrual = true,
            CreditDefaultSwap::PricingModel model = CreditDefaultSwap::Midpoint);
    SpreadCdsHelper(
            Rate spread,
            const Period& tenor,
            Integer settlementDays,
            const Calendar& calendar,
            Frequency frequency,
            BusinessDayConvention convention,
            DateGeneration::Rule rule,
            const DayCounter& dayCounter,
            Real recoveryRate,
            const Handle<YieldTermStructure>& discountCurve,
            bool settlesAccrual = true,
            bool paysAtDefaultTime = true,
            const Date& startDate = Date(),
            const DayCounter& lastPeriodDayCounter = DayCounter(),
            bool rebatesAccrual = true,
            CreditDefaultSwap::PricingModel model = CreditDefaultSwap::Midpoint);
};


%shared_ptr(UpfrontCdsHelper)
class UpfrontCdsHelper : public DefaultProbabilityHelper {
  public:
    UpfrontCdsHelper(
            const Handle<Quote>& upfront,
            Rate spread,
            const Period& tenor,
            Integer settlementDays,
            const Calendar& calendar,
            Frequency frequency,
            BusinessDayConvention convention,
            DateGeneration::Rule rule,
            const DayCounter& dayCounter,
            Real recoveryRate,
            const Handle<YieldTermStructure>& discountCurve,
            Natural upfrontSettlementDays=0,
            bool settlesAccrual = true,
            bool paysAtDefaultTime = true,
            const Date& startDate = Date(),
            const DayCounter& lastPeriodDayCounter = DayCounter(),
            bool rebatesAccrual = true,
            CreditDefaultSwap::PricingModel model = CreditDefaultSwap::Midpoint);
    UpfrontCdsHelper(
            Rate upfront,
            Rate spread,
            const Period& tenor,
            Integer settlementDays,
            const Calendar& calendar,
            Frequency frequency,
            BusinessDayConvention convention,
            DateGeneration::Rule rule,
            const DayCounter& dayCounter,
            Real recoveryRate,
            const Handle<YieldTermStructure>& discountCurve,
            Natural upfrontSettlementDays=0,
            bool settlesAccrual = true,
            bool paysAtDefaultTime = true,
            const Date& startDate = Date(),
            const DayCounter& lastPeriodDayCounter = DayCounter(),
            bool rebatesAccrual = true,
            CreditDefaultSwap::PricingModel model = CreditDefaultSwap::Midpoint);
};



// bootstrap traits

%{
using QuantLib::HazardRate;
using QuantLib::DefaultDensity;
%}

struct HazardRate {};
struct DefaultDensity {};

// curve

%{
using QuantLib::PiecewiseDefaultCurve;
%}

/* We have to resort to a macro, because the R implementation of shared_ptr
   can't take class templates with two or more template arguments. */

%define export_piecewise_default_curve(Name,Traits,Interpolator)

%{
typedef PiecewiseDefaultCurve<Traits, Interpolator> Name;
%}

%shared_ptr(Name);
class Name : public DefaultProbabilityTermStructure {
  public:
    %extend {
        Name(const Date& referenceDate,
             const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
             const DayCounter& dayCounter,
             const Interpolator& i = Interpolator(),
             const _IterativeBootstrap& b = _IterativeBootstrap()) {
            return new Name(referenceDate, instruments, dayCounter,
                            i, Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
        Name(Integer settlementDays, const Calendar& calendar,
             const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
             const DayCounter& dayCounter,
             const Interpolator& i = Interpolator(),
             const _IterativeBootstrap& b = _IterativeBootstrap()) {
            return new Name(settlementDays, calendar, instruments, dayCounter,
                            i, Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
        Name(const Date& referenceDate,
             const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
             const DayCounter& dayCounter,
             const _IterativeBootstrap& b) {
            return new Name(referenceDate, instruments, dayCounter,
                            Interpolator(), Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
        Name(Integer settlementDays, const Calendar& calendar,
             const std::vector<ext::shared_ptr<DefaultProbabilityHelper> >& instruments,
             const DayCounter& dayCounter,
             const _IterativeBootstrap& b) {
            return new Name(settlementDays, calendar, instruments, dayCounter,
                            Interpolator(), Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
    }
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};

%enddef

// add other instantiations if you need them
export_piecewise_default_curve(PiecewiseFlatHazardRate,HazardRate,BackwardFlat);


// bond engine based on default probability

%{
using QuantLib::RiskyBondEngine;
%}

%shared_ptr(RiskyBondEngine)
class RiskyBondEngine : public PricingEngine {
  public:
    RiskyBondEngine(const Handle<DefaultProbabilityTermStructure>& defaultCurve,
                    Real recoveryRate,
                    const Handle<YieldTermStructure>& riskFreeCurve);
};


#endif
