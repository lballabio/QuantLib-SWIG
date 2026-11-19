/*
 Copyright (C) 2026 Raghav Subbarao

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/* FX smile sections of the custom QuantLib build (ql/experimental/fxslv/).

   First slice: settlement convention, delta conversions, smile quotes and
   calibration targets, smile sections (by strike, by delta and the cost
   models), the model factory and the sticky-delta solvers.
*/

#ifndef quantlib_fxslv_i
#define quantlib_fxslv_i

%include common.i
%include types.i
%include date.i
%include calendars.i
%include daycounters.i
%include marketelements.i
%include options.i
%include termstructures.i
%include volatilities.i
%include linearalgebra.i
%include vectors.i

%{
#include <ql/experimental/fxslv/conventions/fxsettlementconvention.hpp>
#include <ql/experimental/fxslv/smilesections/fxdeltaconverter.hpp>
#include <ql/experimental/fxslv/smilesections/fxsmilequotes.hpp>
#include <ql/experimental/fxslv/smilesections/fxsmilesection.hpp>
#include <ql/experimental/fxslv/smilesections/fxsmilesectionbystrike.hpp>
#include <ql/experimental/fxslv/smilesections/fxsmilesectionbydelta.hpp>
#include <ql/experimental/fxslv/smilesections/fxcostsmilesection.hpp>
#include <ql/experimental/fxslv/smilesections/fxsmilesectionfactory.hpp>
#include <ql/experimental/fxslv/smilesections/fxstickydeltasolver.hpp>

using QuantLib::FxSettleConvention;
using QuantLib::FxDeltaConverter;
using QuantLib::FxSmileTarget;
using QuantLib::FxDeltaVolTarget;
using QuantLib::FxRiskReversalTarget;
using QuantLib::FxBrokerStrangleTarget;
using QuantLib::FxSmileQuotes;
using QuantLib::FxRrBfQuotes;
using QuantLib::FxDeltaVolQuotes;
using QuantLib::FxSmileSection;
using QuantLib::FxSmileSectionByStrike;
using QuantLib::FxPolynomialSmileSection;
using QuantLib::FxSabrSmileSection;
using QuantLib::FxSviSmileSection;
using QuantLib::FxSmileSectionByDelta;
using QuantLib::FxQuadraticSmileSection;
using QuantLib::FxCostSmileSection;
using QuantLib::FxCostSmileSectionFlatDynamics;
using QuantLib::FxCostSmileSectionScaledDynamics;
using QuantLib::FxCostSmileSectionStickyDelta;
using QuantLib::FxCostSmileSectionStickyDeltaConvex;
using QuantLib::FxSmileModel;
using QuantLib::makeFxSmileSection;
using QuantLib::FxStickyDeltaSolver;
using QuantLib::FxStickyDeltaSmileSolver;
using QuantLib::FxStickyDeltaConvexSmileSolver;
%}


// ---------------------------------------------------------------------------
// Conventions and delta conversions
// ---------------------------------------------------------------------------

class FxSettleConvention {
  public:
    explicit FxSettleConvention(Calendar calendar, Natural spotLag = 2);
    const Calendar& calendar() const;
    Natural spotLag() const;
    Date valueDate(const Date& tradeDate) const;
};

class FxDeltaConverter {
  public:
    FxDeltaConverter(DeltaVolQuote::DeltaType type,
                     Real spot,
                     DiscountFactor domesticDiscount,
                     DiscountFactor foreignDiscount);
    DeltaVolQuote::DeltaType type() const;
    Real spot() const;
    Real forward() const;
    DiscountFactor domesticDiscount() const;
    DiscountFactor foreignDiscount() const;
    bool premiumAdjusted() const;
    FxDeltaConverter withType(DeltaVolQuote::DeltaType type) const;
    Real delta(Option::Type parity, Rate strike, Real stdDev) const;
    Rate atmStrike(DeltaVolQuote::AtmType type, Real stdDev) const;
    Rate peakCallStrike(Real stdDev) const;
    Rate strike(Option::Type parity, Real delta, Real stdDev) const;
    Real callMinusPutDelta(Rate strike) const;
    Real putDeltaLimit(Rate strike) const;
    Real callDeltaLimit() const;
};


// ---------------------------------------------------------------------------
// Smile quotes and calibration targets
// ---------------------------------------------------------------------------

%template(DeltaVolQuoteHandleVector) std::vector<Handle<DeltaVolQuote> >;

class FxSmileSection;

%shared_ptr(FxSmileTarget)
class FxSmileTarget {
  private:
    FxSmileTarget();
  public:
    Real residual(const FxSmileSection& section) const;
    Real weight(const FxSmileSection& section) const;
    %extend {
        //! Whether the target fixes a (strike, vol) point on the given section's smile.
        bool hasPoint(const FxSmileSection& section) const {
            return self->point(section).has_value();
        }
        //! The (strike, vol) the target fixes; raises if it fixes none.
        std::pair<Real, Real> point(const FxSmileSection& section) const {
            auto p = self->point(section);
            QL_REQUIRE(p.has_value(), "this target fixes no point on the smile");
            return *p;
        }
    }
};

%shared_ptr(FxDeltaVolTarget)
class FxDeltaVolTarget : public FxSmileTarget {
  public:
    explicit FxDeltaVolTarget(Handle<DeltaVolQuote> quote);
};

%shared_ptr(FxRiskReversalTarget)
class FxRiskReversalTarget : public FxSmileTarget {
  public:
    FxRiskReversalTarget(Real delta, Volatility riskReversal);
};

%shared_ptr(FxBrokerStrangleTarget)
class FxBrokerStrangleTarget : public FxSmileTarget {
  public:
    FxBrokerStrangleTarget(Volatility atmVol, Volatility brokerFly, Real delta);
};

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<FxSmileTarget> )
#endif
%template(FxSmileTargets) std::vector<ext::shared_ptr<FxSmileTarget> >;

%shared_ptr(FxSmileQuotes)
class FxSmileQuotes : public Observable {
  private:
    FxSmileQuotes();
  public:
    Volatility referenceVol() const;
};

%shared_ptr(FxRrBfQuotes)
class FxRrBfQuotes : public FxSmileQuotes {
  public:
    enum FlyType { SmileStrangle, MarketStrangle };
    FxRrBfQuotes(Handle<Quote> atm,
                 std::vector<Handle<Quote> > riskReversals,
                 std::vector<Handle<Quote> > butterflies,
                 std::vector<Real> deltas,
                 FxRrBfQuotes::FlyType flyType);
    const Handle<Quote>& atm() const;
    const std::vector<Handle<Quote> >& riskReversals() const;
    const std::vector<Handle<Quote> >& butterflies() const;
    const std::vector<Real>& deltas() const;
    FxRrBfQuotes::FlyType flyType() const;
};

%shared_ptr(FxDeltaVolQuotes)
class FxDeltaVolQuotes : public FxSmileQuotes {
  public:
    explicit FxDeltaVolQuotes(std::vector<Handle<DeltaVolQuote> > quotes);
    const std::vector<Handle<DeltaVolQuote> >& quotes() const;
};


// ---------------------------------------------------------------------------
// Smile sections
// ---------------------------------------------------------------------------

/* Every section has a date-mode constructor (expiry date, day counter and
   settlement convention, optional fixed reference date) and a time-mode
   constructor (expiry time; no settlement lag). */

%define FX_SECTION_CONSTRUCTORS(Name)
    Name(const Date& exerciseDate,
         const Handle<Quote>& spot,
         const ext::shared_ptr<FxSmileQuotes>& quotes,
         const Handle<YieldTermStructure>& foreignDiscount,
         const Handle<YieldTermStructure>& domesticDiscount,
         DeltaVolQuote::DeltaType deltaType,
         DeltaVolQuote::AtmType atmType,
         const DayCounter& dayCounter,
         const FxSettleConvention& settlement,
         const Date& referenceDate = Date());
    Name(Time exerciseTime,
         const Handle<Quote>& spot,
         const ext::shared_ptr<FxSmileQuotes>& quotes,
         const Handle<YieldTermStructure>& foreignDiscount,
         const Handle<YieldTermStructure>& domesticDiscount,
         DeltaVolQuote::DeltaType deltaType,
         DeltaVolQuote::AtmType atmType,
         const DayCounter& dayCounter = DayCounter());
%enddef

%define FX_COST_SECTION_CONSTRUCTORS(Name)
    Name(const Date& exerciseDate,
         const Handle<Quote>& spot,
         const ext::shared_ptr<FxSmileQuotes>& quotes,
         const Handle<YieldTermStructure>& foreignDiscount,
         const Handle<YieldTermStructure>& domesticDiscount,
         DeltaVolQuote::DeltaType deltaType,
         DeltaVolQuote::AtmType atmType,
         const DayCounter& dayCounter,
         const FxSettleConvention& settlement,
         const Date& referenceDate = Date(),
         bool weightedCalibrationFlag = true);
    Name(Time exerciseTime,
         const Handle<Quote>& spot,
         const ext::shared_ptr<FxSmileQuotes>& quotes,
         const Handle<YieldTermStructure>& foreignDiscount,
         const Handle<YieldTermStructure>& domesticDiscount,
         DeltaVolQuote::DeltaType deltaType,
         DeltaVolQuote::AtmType atmType,
         const DayCounter& dayCounter = DayCounter(),
         bool weightedCalibrationFlag = true);
%enddef

/* FxSmileSection is also a LazyObject in C++; only SmileSection is declared
   as a base here (both derive from Observable), and the LazyObject methods
   are added below. */
%shared_ptr(FxSmileSection)
class FxSmileSection : public SmileSection {
  private:
    FxSmileSection();
  public:
    DeltaVolQuote::DeltaType deltaType() const;
    DeltaVolQuote::AtmType atmType() const;
    bool premiumAdjust() const;

    Handle<Quote> spot() const;
    Volatility atmVol() const;
    Rate atmStrike() const;
    Real forward() const;
    Handle<YieldTermStructure> foreignDiscount() const;
    Handle<YieldTermStructure> domesticDiscount() const;
    Date spotDate() const;
    Date deliveryDate() const;
    DiscountFactor domesticDiscountFactor() const;
    DiscountFactor foreignDiscountFactor() const;
    FxDeltaConverter deltaConverter() const;
    Real maxCallDelta() const;
    const ext::shared_ptr<FxSmileQuotes>& smileQuotes() const;

    const Array& calibrationResiduals() const;
    Real calibrationError() const;

    Volatility volByStrike(Rate strike) const;
    Volatility volByDelta(Real delta, Option::Type parity) const;
    Real deltaByStrike(Rate strike, Option::Type parity) const;
    Rate strikeByDelta(Real delta, Option::Type parity) const;
    Real volResidual(Rate strike, Volatility vol) const;
    Real volDerivative(Rate strike) const;
    Real atmVolForwardSensitivity() const;

    Real normedCallPrice(Real moneyness) const;
    Probability exerciseProbability(Real moneyness) const;
    Real moneynessFromProbability(Probability p) const;

    %extend {
        //! Whether the section was built in date mode (with a settlement convention).
        bool hasSettleConvention() const {
            return self->settleConvention().has_value();
        }
        FxSettleConvention settleConvention() const {
            QL_REQUIRE(self->settleConvention().has_value(),
                       "no settlement convention: the section was built in time mode");
            return *(self->settleConvention());
        }
        // LazyObject interface
        void recalculate() { self->recalculate(); }
        void freeze() { self->freeze(); }
        void unfreeze() { self->unfreeze(); }
    }
};

%shared_ptr(FxSmileSectionByStrike)
class FxSmileSectionByStrike : public FxSmileSection {
  private:
    FxSmileSectionByStrike();
};

%shared_ptr(FxPolynomialSmileSection)
class FxPolynomialSmileSection : public FxSmileSectionByStrike {
  public:
    FX_SECTION_CONSTRUCTORS(FxPolynomialSmileSection)
    Real a() const;
    Real b() const;
    Real c() const;
};

%shared_ptr(FxSabrSmileSection)
class FxSabrSmileSection : public FxSmileSectionByStrike {
  public:
    FX_SECTION_CONSTRUCTORS(FxSabrSmileSection)
    Real alpha() const;
    Real beta() const;
    Real nu() const;
    Real rho() const;
};

%shared_ptr(FxSviSmileSection)
class FxSviSmileSection : public FxSmileSectionByStrike {
  public:
    FX_SECTION_CONSTRUCTORS(FxSviSmileSection)
    Real a() const;
    Real b() const;
    Real rho() const;
    Real m() const;
    Real sigma() const;
    Real minimumVariance() const;
};

%shared_ptr(FxSmileSectionByDelta)
class FxSmileSectionByDelta : public FxSmileSection {
  private:
    FxSmileSectionByDelta();
};

%shared_ptr(FxQuadraticSmileSection)
class FxQuadraticSmileSection : public FxSmileSectionByDelta {
  public:
    FX_SECTION_CONSTRUCTORS(FxQuadraticSmileSection)
    Real a() const;
    Real b() const;
    Real c() const;
};

%shared_ptr(FxCostSmileSection)
class FxCostSmileSection : public FxSmileSectionByStrike {
  public:
    FxCostSmileSection(const Date& exerciseDate,
                       const Handle<Quote>& spot,
                       const ext::shared_ptr<FxSmileQuotes>& quotes,
                       const Handle<YieldTermStructure>& foreignDiscount,
                       const Handle<YieldTermStructure>& domesticDiscount,
                       DeltaVolQuote::DeltaType deltaType,
                       DeltaVolQuote::AtmType atmType,
                       Real alpha,
                       const DayCounter& dayCounter,
                       const FxSettleConvention& settlement,
                       const Date& referenceDate = Date(),
                       bool weightedCalibrationFlag = true);
    FxCostSmileSection(Time exerciseTime,
                       const Handle<Quote>& spot,
                       const ext::shared_ptr<FxSmileQuotes>& quotes,
                       const Handle<YieldTermStructure>& foreignDiscount,
                       const Handle<YieldTermStructure>& domesticDiscount,
                       DeltaVolQuote::DeltaType deltaType,
                       DeltaVolQuote::AtmType atmType,
                       Real alpha = 0.5,
                       const DayCounter& dayCounter = DayCounter(),
                       bool weightedCalibrationFlag = true);
    bool weightedCalibration() const;
    Real alpha() const;
    const std::vector<Real>& costs() const;
};

%shared_ptr(FxCostSmileSectionFlatDynamics)
class FxCostSmileSectionFlatDynamics : public FxCostSmileSection {
  public:
    FX_COST_SECTION_CONSTRUCTORS(FxCostSmileSectionFlatDynamics)
};

%shared_ptr(FxCostSmileSectionScaledDynamics)
class FxCostSmileSectionScaledDynamics : public FxCostSmileSection {
  public:
    FX_COST_SECTION_CONSTRUCTORS(FxCostSmileSectionScaledDynamics)
};

%shared_ptr(FxCostSmileSectionStickyDelta)
class FxCostSmileSectionStickyDelta : public FxCostSmileSection {
  public:
    FX_COST_SECTION_CONSTRUCTORS(FxCostSmileSectionStickyDelta)
    Real slope(Rate strike) const;
    Real convexity(Rate strike) const;
    %extend {
        //! Total vol range, in units of the reference total vol, over which the smile is solved.
        static Real smileRange() { return FxCostSmileSectionStickyDelta::smileRange; }
    }
};

%shared_ptr(FxCostSmileSectionStickyDeltaConvex)
class FxCostSmileSectionStickyDeltaConvex : public FxCostSmileSectionStickyDelta {
  public:
    FX_COST_SECTION_CONSTRUCTORS(FxCostSmileSectionStickyDeltaConvex)
};

enum class FxSmileModel {
    Sabr,
    Svi,
    Polynomial,
    Quadratic,
    CostFlatDynamics,
    CostScaledDynamics,
    CostStickyDelta,
    CostStickyDeltaConvex
};

ext::shared_ptr<FxSmileSection>
makeFxSmileSection(FxSmileModel model,
                   const Date& expiry,
                   const Handle<Quote>& spot,
                   const ext::shared_ptr<FxSmileQuotes>& quotes,
                   const Handle<YieldTermStructure>& foreignCurve,
                   const Handle<YieldTermStructure>& domesticCurve,
                   DeltaVolQuote::DeltaType deltaType,
                   DeltaVolQuote::AtmType atmType,
                   const FxSettleConvention& settlement,
                   const Date& referenceDate = Date());


/* Downcasts: the factory and FxVarianceSurface hold sections as FxSmileSection.
   Each returns None if the section is not of that class. */
%define FX_SECTION_DOWNCAST(Name, function)
%inline %{
    ext::shared_ptr<Name> function(const ext::shared_ptr<FxSmileSection>& s) {
        return ext::dynamic_pointer_cast<Name>(s);
    }
%}
%enddef

FX_SECTION_DOWNCAST(FxSmileSectionByStrike, as_fx_smile_section_by_strike)
FX_SECTION_DOWNCAST(FxPolynomialSmileSection, as_fx_polynomial_smile_section)
FX_SECTION_DOWNCAST(FxSabrSmileSection, as_fx_sabr_smile_section)
FX_SECTION_DOWNCAST(FxSviSmileSection, as_fx_svi_smile_section)
FX_SECTION_DOWNCAST(FxSmileSectionByDelta, as_fx_smile_section_by_delta)
FX_SECTION_DOWNCAST(FxQuadraticSmileSection, as_fx_quadratic_smile_section)
FX_SECTION_DOWNCAST(FxCostSmileSection, as_fx_cost_smile_section)
FX_SECTION_DOWNCAST(FxCostSmileSectionFlatDynamics, as_fx_cost_smile_section_flat_dynamics)
FX_SECTION_DOWNCAST(FxCostSmileSectionScaledDynamics, as_fx_cost_smile_section_scaled_dynamics)
FX_SECTION_DOWNCAST(FxCostSmileSectionStickyDelta, as_fx_cost_smile_section_sticky_delta)
FX_SECTION_DOWNCAST(FxCostSmileSectionStickyDeltaConvex, as_fx_cost_smile_section_sticky_delta_convex)


// ---------------------------------------------------------------------------
// Sticky-delta smile solvers
// ---------------------------------------------------------------------------

%shared_ptr(FxStickyDeltaSolver)
class FxStickyDeltaSolver {
  private:
    FxStickyDeltaSolver();
  public:
    Real totalVol(Real x) const;
    Real slope(Real x) const;
    Real convexity(Real x) const;
    const std::vector<Real>& costs() const;
    bool hasCosts(const std::vector<Real>& costs) const;
    bool premiumAdjust() const;
    Real minLogMoneyness() const;
    Real maxLogMoneyness() const;
    Real costEquation(Real x, Real w, Real beta) const;
};

%shared_ptr(FxStickyDeltaSmileSolver)
class FxStickyDeltaSmileSolver : public FxStickyDeltaSolver {
  public:
    FxStickyDeltaSmileSolver(std::vector<Real> costs,
                             bool premiumAdjust,
                             Real xMin,
                             Real xMax,
                             Real stepFraction = 0.1,
                             Real accuracy = 1.0e-8);
    Real singularLogMoneyness() const;
    Real singularTotalVol() const;
    Real singularSlope() const;
    Real quadraticRoot(Real x, Real w, Real betaRef) const;
};

%shared_ptr(FxStickyDeltaConvexSmileSolver)
class FxStickyDeltaConvexSmileSolver : public FxStickyDeltaSolver {
  public:
    FxStickyDeltaConvexSmileSolver(std::vector<Real> costs,
                                   bool premiumAdjust,
                                   Real xMin,
                                   Real xMax,
                                   Size nodes = 401,
                                   Real tolerance = 1.0e-13);
};

#endif
