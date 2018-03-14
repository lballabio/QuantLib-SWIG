/*
 Copyright (C) 2008, 2009 StatPro Italia srl

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

#ifndef quantlib_credit_default_swap_i
#define quantlib_credit_default_swap_i

%include instruments.i
%include credit.i
%include termstructures.i
%include null.i

%{
using QuantLib::CreditDefaultSwap;
using QuantLib::MidPointCdsEngine;
using QuantLib::IntegralCdsEngine;
using QuantLib::IsdaCdsEngine;
using QuantLib::Claim;

typedef boost::shared_ptr<Instrument> CreditDefaultSwapPtr;
typedef boost::shared_ptr<PricingEngine> MidPointCdsEnginePtr;
typedef boost::shared_ptr<PricingEngine> IntegralCdsEnginePtr;
typedef boost::shared_ptr<PricingEngine> IsdaCdsEnginePtr;
%}

#if defined(SWIGJAVA) || defined(SWIGCSHARP)
%rename(_CreditDefaultSwap) CreditDefaultSwap;
#else
%ignore CreditDefaultSwap;
#endif
class CreditDefaultSwap {
  public:
    enum PricingModel { Midpoint, ISDA };
#if defined(SWIGJAVA) || defined(SWIGCSHARP)
  private:
    CreditDefaultSwap();
#endif
};

%rename(CreditDefaultSwap) CreditDefaultSwapPtr;
class CreditDefaultSwapPtr : public boost::shared_ptr<Instrument> {
  public:
    %extend {

        static const CreditDefaultSwap::PricingModel Midpoint = CreditDefaultSwap::Midpoint;
        static const CreditDefaultSwap::PricingModel ISDA = CreditDefaultSwap::ISDA;

        CreditDefaultSwapPtr(Protection::Side side,
                             Real notional,
                             Rate spread,
                             const Schedule& schedule,
                             BusinessDayConvention paymentConvention,
                             const DayCounter& dayCounter,
                             bool settlesAccrual = true,
                             bool paysAtDefaultTime = true,
                             const Date& protectionStart = Date(),
                             const DayCounter& lastPeriodDayCounter = DayCounter(),
                             const bool rebatesAccrual = true) {
            return new CreditDefaultSwapPtr(
                    new CreditDefaultSwap(side, notional, spread, schedule,
                                          paymentConvention, dayCounter,
                                          settlesAccrual, paysAtDefaultTime,
                                          protectionStart, boost::shared_ptr<Claim>(),
                                          lastPeriodDayCounter, rebatesAccrual));
        }
        CreditDefaultSwapPtr(Protection::Side side,
                             Real notional,
                             Rate upfront,
                             Rate spread,
                             const Schedule& schedule,
                             BusinessDayConvention paymentConvention,
                             const DayCounter& dayCounter,
                             bool settlesAccrual = true,
                             bool paysAtDefaultTime = true,
                             const Date& protectionStart = Date(),
                             const Date& upfrontDate = Date(),
                             const DayCounter& lastPeriodDayCounter = DayCounter(),
                             const bool rebatesAccrual = true) {
            return new CreditDefaultSwapPtr(
                    new CreditDefaultSwap(side, notional, upfront, spread,
                                          schedule, paymentConvention,
                                          dayCounter, settlesAccrual,
                                          paysAtDefaultTime, protectionStart,
                                          upfrontDate, boost::shared_ptr<Claim>(),
                                          lastPeriodDayCounter, rebatesAccrual));
        }
        Protection::Side side() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->side();
        }
        Real notional() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->notional();
        }
        Rate runningSpread() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->runningSpread();
        }
        doubleOrNull upfront() const {
            boost::optional<Rate> result =
                boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->upfront();
            if (result)
                return *result;
            else
                return Null<double>();
        }
        bool settlesAccrual() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->settlesAccrual();
        }
        bool paysAtDefaultTime() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->paysAtDefaultTime();
        }
        Rate fairSpread() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->fairSpread();
        }
        Rate fairUpfront() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->fairUpfront();
        }
        Real couponLegBPS() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->couponLegBPS();
        }
        Real couponLegNPV() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->couponLegNPV();
        }
        Real defaultLegNPV() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->defaultLegNPV();
        }
        Real upfrontBPS() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->upfrontBPS();
        }
        Real upfrontNPV() const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->upfrontNPV();
        }
        Rate impliedHazardRate(Real targetNPV,
                               const Handle<YieldTermStructure>& discountCurve,
                               const DayCounter& dayCounter,
                               Real recoveryRate = 0.4,
                               Real accuracy = 1.0e-6,
                               CreditDefaultSwap::PricingModel model = CreditDefaultSwap::ISDA) const {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->impliedHazardRate(targetNPV, discountCurve, dayCounter,
                                    recoveryRate, accuracy, model);
        }
        std::vector<boost::shared_ptr<CashFlow> > coupons() {
            return boost::dynamic_pointer_cast<CreditDefaultSwap>(*self)
                ->coupons();
        }
    }
};


%rename(MidPointCdsEngine) MidPointCdsEnginePtr;
class MidPointCdsEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        MidPointCdsEnginePtr(
                   const Handle<DefaultProbabilityTermStructure>& probability,
                   Real recoveryRate,
                   const Handle<YieldTermStructure>& discountCurve) {
            return new MidPointCdsEnginePtr(
                              new MidPointCdsEngine(probability, recoveryRate,
                                                    discountCurve));
        }
    }
};

%rename(IntegralCdsEngine) IntegralCdsEnginePtr;
class IntegralCdsEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        IntegralCdsEnginePtr(
				   const Period &integrationStep,
                   const Handle<DefaultProbabilityTermStructure>& probability,
                   Real recoveryRate,
                   const Handle<YieldTermStructure>& discountCurve,
				   bool includeSettlementDateFlows = false) {
            return new IntegralCdsEnginePtr(
                              new IntegralCdsEngine(integrationStep, probability,
                                                    recoveryRate, discountCurve,
													includeSettlementDateFlows));
        }
    }
};

#if defined(SWIGJAVA) || defined(SWIGCSHARP)
%rename(_IsdaCdsEngine) IsdaCdsEngine;
#else
%ignore IsdaCdsEngine;
#endif
class IsdaCdsEngine {
  public:
    enum NumericalFix { None, Taylor };
    enum AccrualBias { HalfDayBias, NoBias };
    enum ForwardsInCouponPeriod { Flat, Piecewise };
#if defined(SWIGJAVA) || defined(SWIGCSHARP)
  private:
    IsdaCdsEngine();
#endif
};

%rename(IsdaCdsEngine) IsdaCdsEnginePtr;
class IsdaCdsEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
      
      static const IsdaCdsEngine::NumericalFix None = IsdaCdsEngine::None;
      static const IsdaCdsEngine::NumericalFix Taylor = IsdaCdsEngine::Taylor;

      static const IsdaCdsEngine::AccrualBias HalfDayBias = IsdaCdsEngine::HalfDayBias;
      static const IsdaCdsEngine::AccrualBias NoBias = IsdaCdsEngine::NoBias;

      static const IsdaCdsEngine::ForwardsInCouponPeriod Flat = IsdaCdsEngine::Flat;
      static const IsdaCdsEngine::ForwardsInCouponPeriod Piecewise = IsdaCdsEngine::Piecewise;
      
      IsdaCdsEnginePtr(
                   const Handle<DefaultProbabilityTermStructure>& probability,
                   Real recoveryRate,
                   const Handle<YieldTermStructure>& discountCurve,
                   bool includeSettlementDateFlows = false,
                   const IsdaCdsEngine::NumericalFix numericalFix = IsdaCdsEngine::Taylor,
                   const IsdaCdsEngine::AccrualBias accrualBias = IsdaCdsEngine::HalfDayBias,
                   const IsdaCdsEngine::ForwardsInCouponPeriod forwardsInCouponPeriod = IsdaCdsEngine::Piecewise) {
          
          return new IsdaCdsEnginePtr(
            new IsdaCdsEngine(probability, recoveryRate,
                              discountCurve, includeSettlementDateFlows,
                              numericalFix, accrualBias, forwardsInCouponPeriod)
          );
        }
    }
};


#endif
