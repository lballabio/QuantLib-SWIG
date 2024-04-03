/*
 Copyright (C) 2004, 2005, 2006, 2007 StatPro Italia srl
 Copyright (C) 2009 Joseph Malicki
 Copyright (C) 2011 Lluis Pujol Bajador
 Copyright (C) 2014 Simon Mazzucca
 Copyright (C) 2016 Gouthaman Balaraman
 Copyright (C) 2017 BN Algorithms Ltd
 Copyright (C) 2018 Matthias Groncki
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

#ifndef quantlib_bonds_i
#define quantlib_bonds_i

%include instruments.i
%include calendars.i
%include daycounters.i
%include cashflows.i
%include interestrate.i
%include indexes.i
%include inflation.i
%include shortratemodels.i

%{
using QuantLib::Bond;
typedef Bond::Price BondPrice;
using QuantLib::ZeroCouponBond;
using QuantLib::FixedRateBond;
using QuantLib::AmortizingFixedRateBond;
using QuantLib::FloatingRateBond;
using QuantLib::AmortizingFloatingRateBond;
using QuantLib::DiscountingBondEngine;
using QuantLib::simplifyNotificationGraph;
%}

class BondPrice {
  public:
    enum Type { Dirty, Clean };
    BondPrice(Real amount, Type type);
    Real amount() const;
    Type type() const;
    bool isValid() const;
};

%shared_ptr(Bond)
class Bond : public Instrument {
    #if defined(SWIGPYTHON)
    %rename(bondYield) yield;
    #endif
  public:
    Bond(Natural settlementDays,
            const Calendar& calendar,
            Real faceAmount,
            const Date& maturityDate,
            const Date& issueDate = Date(),
            const Leg& cashflows = Leg());
    Bond(Natural settlementDays,
            const Calendar& calendar,
            const Date& issueDate = Date(),
            const Leg& coupons = Leg());
    // public functions
    Rate nextCouponRate(const Date& d = Date());
    Rate previousCouponRate(const Date& d = Date());
    // inspectors
    Natural settlementDays() const;
    Date settlementDate(Date d = Date());
    Date startDate() const;
    Date maturityDate() const;
    Date issueDate() const;
    std::vector<ext::shared_ptr<CashFlow> > cashflows() const;
    std::vector<ext::shared_ptr<CashFlow> > redemptions() const;
    ext::shared_ptr<CashFlow> redemption() const;
    Calendar calendar() const;
    std::vector<Real> notionals() const;
    Real notional(Date d = Date()) const;
    // calculations
    Real cleanPrice();
    Real cleanPrice(Rate yield,
                    const DayCounter &dc,
                    Compounding compounding,
                    Frequency frequency,
                    const Date& settlement = Date());
    Real dirtyPrice();
    Real dirtyPrice(Rate yield,
                    const DayCounter &dc,
                    Compounding compounding,
                    Frequency frequency,
                    const Date& settlement = Date());
    Real yield(const DayCounter& dc,
               Compounding compounding,
               Frequency freq,
               Real accuracy = 1.0e-8,
               Size maxEvaluations = 100);
    Real yield(BondPrice price,
               const DayCounter& dc,
               Compounding compounding,
               Frequency freq,
               const Date& settlement = Date(),
               Real accuracy = 1.0e-8,
               Size maxEvaluations = 100,
               Real guess = 0.05);
    Real yield(Real cleanPrice,
               const DayCounter& dc,
               Compounding compounding,
               Frequency freq,
               const Date& settlement = Date(),
               Real accuracy = 1.0e-8,
               Size maxEvaluations = 100);
    Real accruedAmount(const Date& settlement = Date());
    Real settlementValue() const;
    Real settlementValue(Real cleanPrice) const;
};

void simplifyNotificationGraph(Bond& bond, bool unregisterCoupons = false);


%inline %{
    Real cleanPriceFromZSpread(
                   const Bond& bond,
                   const ext::shared_ptr<YieldTermStructure>& discountCurve,
                   Spread zSpread,
                   const DayCounter& dc,
                   Compounding compounding,
                   Frequency freq,
                   const Date& settlementDate = Date()) {
        return QuantLib::BondFunctions::cleanPrice(
                                  bond,
                                  discountCurve,
                                  zSpread, dc, compounding,
                                  freq, settlementDate);
    }

%}

namespace QuantLib {

    Schedule sinkingSchedule(const Date& startDate,
                             const Period& bondLength,
                             const Frequency& frequency,
                             const Calendar& paymentCalendar);

    std::vector<Real> sinkingNotionals(const Period& bondLength,
                                       const Frequency& frequency,
                                       Rate couponRate,
                                       Real initialNotional);

}

%shared_ptr(ZeroCouponBond)
class ZeroCouponBond : public Bond {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") ZeroCouponBond;
    #endif
  public:
    ZeroCouponBond(
            Natural settlementDays,
            const Calendar &calendar,
            Real faceAmount,
            const Date & maturityDate,
            BusinessDayConvention paymentConvention = QuantLib::Following,
            Real redemption = 100.0,
            const Date& issueDate = Date());
};

%shared_ptr(FixedRateBond)
class FixedRateBond : public Bond {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") FixedRateBond;
    #endif
  public:
    FixedRateBond(
            Integer settlementDays,
            Real faceAmount,
            const Schedule &schedule,
            const std::vector<Rate>& coupons,
            const DayCounter& paymentDayCounter,
            BusinessDayConvention paymentConvention = QuantLib::Following,
            Real redemption = 100.0,
            Date issueDate = Date(),
            const Calendar& paymentCalendar = Calendar(),
            const Period& exCouponPeriod = Period(),
            const Calendar& exCouponCalendar = Calendar(),
            BusinessDayConvention exCouponConvention = Unadjusted,
            bool exCouponEndOfMonth = false);
    Frequency frequency() const;
    DayCounter dayCounter() const;
};


%shared_ptr(AmortizingFixedRateBond)
class AmortizingFixedRateBond : public Bond {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") AmortizingFixedRateBond;
    #endif
  public:
    AmortizingFixedRateBond(
            Integer settlementDays,
            const std::vector<Real>& notionals,
            const Schedule& schedule,
            const std::vector<Rate>& coupons,
            const DayCounter& accrualDayCounter,
            BusinessDayConvention paymentConvention = QuantLib::Following,
            Date issueDate = Date(),
            const Period& exCouponPeriod = Period(),
            const Calendar& exCouponCalendar = Calendar(),
            const BusinessDayConvention exCouponConvention = Unadjusted,
            bool exCouponEndOfMonth = false,
            const std::vector<Real>& redemptions = { 100.0 },
            Integer paymentLag = 0);
    Frequency frequency() const;
    DayCounter dayCounter() const;
};


%shared_ptr(AmortizingFloatingRateBond)
class AmortizingFloatingRateBond : public Bond {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") AmortizingFloatingRateBond;
    #endif
  public:
    AmortizingFloatingRateBond(
        Size settlementDays,
        const std::vector<Real>& notional,
        const Schedule& schedule,
        const ext::shared_ptr<IborIndex>& index,
        const DayCounter& accrualDayCounter,
        BusinessDayConvention paymentConvention = Following,
        Size fixingDays = Null<Size>(),
        const std::vector<Real>& gearings = std::vector<Real>(1, 1.0),
        const std::vector<Spread>& spreads = std::vector<Spread>(1, 0.0),
        const std::vector<Rate>& caps = std::vector<Rate>(),
        const std::vector<Rate>& floors = std::vector<Rate>(),
        bool inArrears = false,
        const Date& issueDate = Date(),
        const Period& exCouponPeriod = Period(),
        const Calendar& exCouponCalendar = Calendar(),
        const BusinessDayConvention exCouponConvention = Unadjusted,
        bool exCouponEndOfMonth = false,
        const std::vector<Real>& redemptions = { 100.0 },
        Integer paymentLag = 0);
};


%shared_ptr(FloatingRateBond)
class FloatingRateBond : public Bond {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") FloatingRateBond;
    #endif
  public:
    FloatingRateBond(
        Size settlementDays,
        Real faceAmount,
        const Schedule& schedule,
        const ext::shared_ptr<IborIndex>& index,
        const DayCounter& paymentDayCounter,
        BusinessDayConvention paymentConvention = Following,
        Size fixingDays = Null<Size>(),
        const std::vector<Real>& gearings = std::vector<Real>(),
        const std::vector<Spread>& spreads = std::vector<Spread>(),
        const std::vector<Rate>& caps = std::vector<Rate>(),
        const std::vector<Rate>& floors = std::vector<Rate>(),
        bool inArrears = false,
        Real redemption = 100.0,
        const Date& issueDate = Date(),
        const Period& exCouponPeriod = Period(),
        const Calendar& exCouponCalendar = Calendar(),
        BusinessDayConvention exCouponConvention = Unadjusted,
        bool exCouponEndOfMonth = false);
};


%{
using QuantLib::CmsRateBond;
using QuantLib::AmortizingCmsRateBond;
%}

%shared_ptr(CmsRateBond)
class CmsRateBond : public Bond {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CmsRateBond;
    #endif
  public:
    CmsRateBond(Size settlementDays,
                   Real faceAmount,
                   const Schedule& schedule,
                   const ext::shared_ptr<SwapIndex>& index,
                   const DayCounter& paymentDayCounter,
                   BusinessDayConvention paymentConvention,
                   Natural fixingDays,
                   const std::vector<Real>& gearings,
                   const std::vector<Spread>& spreads,
                   const std::vector<Rate>& caps,
                   const std::vector<Rate>& floors,
                   bool inArrears = false,
                   Real redemption = 100.0,
                   const Date& issueDate = Date());
};


%shared_ptr(AmortizingCmsRateBond)
class AmortizingCmsRateBond : public Bond {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") AmortizingCmsRateBond;
    #endif
  public:
    AmortizingCmsRateBond(
                Natural settlementDays,
                const std::vector<Real>& notionals,
                const Schedule& schedule,
                const ext::shared_ptr<SwapIndex>& index,
                const DayCounter& paymentDayCounter,
                BusinessDayConvention paymentConvention = Following,
                Natural fixingDays = Null<Natural>(),
                const std::vector<Real>& gearings = { 1.0 },
                const std::vector<Spread>& spreads = { 0.0 },
                const std::vector<Rate>& caps = {},
                const std::vector<Rate>& floors = {},
                bool inArrears = false,
                const Date& issueDate = Date());
};


%shared_ptr(DiscountingBondEngine)
class DiscountingBondEngine : public PricingEngine {
  public:
    DiscountingBondEngine(const Handle<YieldTermStructure>& discountCurve);
};


%{
using QuantLib::CallableBond;
using QuantLib::Callability;
using QuantLib::SoftCallability;
using QuantLib::CallabilitySchedule;

using QuantLib::CallableFixedRateBond;
using QuantLib::CallableZeroCouponBond;
using QuantLib::TreeCallableFixedRateBondEngine;
using QuantLib::BlackCallableFixedRateBondEngine;
%}

%shared_ptr(Callability)
class Callability {
  public:
    enum Type { Call, Put };
    Callability(const BondPrice& price,
                Type type,
                const Date& date);
    const BondPrice& price() const;
    Type type() const;
    Date date() const;
};

%shared_ptr(SoftCallability)
class SoftCallability : public Callability {
  public:
    SoftCallability(const BondPrice& price,
                    const Date& date,
                    Real trigger);
};

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<Callability> )
#endif
namespace std {
    %template(CallabilitySchedule) vector<ext::shared_ptr<Callability> >;
}


%shared_ptr(CallableBond)
class CallableBond : public Bond {
  private:
    CallableBond();
  public:
    const std::vector<ext::shared_ptr<Callability> >& callability() const;

    Volatility impliedVolatility(const BondPrice& targetPrice,
                                 const Handle<YieldTermStructure>& discountCurve,
                                 Real accuracy,
                                 Size maxEvaluations,
                                 Volatility minVol,
                                 Volatility maxVol) const;

    Real OAS(Real cleanPrice,
             const Handle<YieldTermStructure>& engineTS,
             const DayCounter& dc,
             Compounding compounding,
             Frequency freq,
             const Date& settlementDate = Date(),
             Real accuracy =1e-10,
             Size maxIterations = 100,
             Spread guess = 0.0);

    Real cleanPriceOAS(Real oas,
                       const Handle<YieldTermStructure>& engineTS,
                       const DayCounter& dayCounter,
                       Compounding compounding,
                       Frequency frequency,
                       Date settlementDate = Date());

    Real effectiveDuration(Real oas,
                           const Handle<YieldTermStructure>& engineTS,
                           const DayCounter& dayCounter,
                           Compounding compounding,
                           Frequency frequency,
                           Real bump=2e-4);

    Real effectiveConvexity(Real oas,
                            const Handle<YieldTermStructure>& engineTS,
                            const DayCounter& dayCounter,
                            Compounding compounding,
                            Frequency frequency,
                            Real bump=2e-4);
};


%shared_ptr(CallableFixedRateBond)
class CallableFixedRateBond : public CallableBond {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CallableFixedRateBond;
    #endif
  public:
    CallableFixedRateBond(
            Integer settlementDays,
            Real faceAmount,
            const Schedule &schedule,
            const std::vector<Rate>& coupons,
            const DayCounter& accrualDayCounter,
            BusinessDayConvention paymentConvention,
            Real redemption,
            Date issueDate,
            const std::vector<ext::shared_ptr<Callability> >& putCallSchedule,
            const Period& exCouponPeriod = Period(),
            const Calendar& exCouponCalendar = Calendar(),
            BusinessDayConvention exCouponConvention = Unadjusted,
            bool exCouponEndOfMonth = false);
};


%shared_ptr(CallableZeroCouponBond)
class CallableZeroCouponBond : public CallableBond {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CallableZeroCouponBond;
    #endif
  public:
    CallableZeroCouponBond(
            Integer settlementDays,
            Real faceAmount,
            const Calendar& calendar,
            const Date& maturityDate,
            const DayCounter& dayCounter,
            BusinessDayConvention paymentConvention = Following,
            Real redemption = 100.0,
            const Date& issueDate = Date(),
            const std::vector<ext::shared_ptr<Callability> >& putCallSchedule
                           = std::vector<ext::shared_ptr<Callability> >());
};


%shared_ptr(TreeCallableFixedRateBondEngine)
class TreeCallableFixedRateBondEngine : public PricingEngine {
  public:
    TreeCallableFixedRateBondEngine(
                         const ext::shared_ptr<ShortRateModel>& model,
                         Size timeSteps,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>());
    TreeCallableFixedRateBondEngine(
                         const ext::shared_ptr<ShortRateModel>& model,
                         const TimeGrid& grid,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>());
};

%shared_ptr(BlackCallableFixedRateBondEngine)
class BlackCallableFixedRateBondEngine : public PricingEngine {
  public:
    BlackCallableFixedRateBondEngine(
                const Handle<Quote>& fwdYieldVol,
                const Handle<YieldTermStructure>& discountCurve);
};

%{
using QuantLib::CPIBond;
%}

%shared_ptr(CPIBond)
class CPIBond : public Bond {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CPIBond;
    #endif
  public:
    CPIBond(
            Natural settlementDays,
            Real faceAmount,
            bool growthOnly,
            Real baseCPI,
            const Period& observationLag,
            const ext::shared_ptr<ZeroInflationIndex>& cpiIndex,
            CPI::InterpolationType observationInterpolation,
            const Schedule& schedule,
            const std::vector<Rate>& coupons,
            const DayCounter& accrualDayCounter,
            BusinessDayConvention paymentConvention = ModifiedFollowing,
            const Date& issueDate = Date(),
            const Calendar& paymentCalendar = Calendar(),
            const Period& exCouponPeriod = Period(),
            const Calendar& exCouponCalendar = Calendar(),
            BusinessDayConvention exCouponConvention = Unadjusted,
            bool exCouponEndOfMonth = false);
};


#endif
