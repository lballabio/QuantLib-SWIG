/*
 Copyright (C) 2008, 2009 StatPro Italia srl
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

#ifndef quantlib_credit_default_swap_i
#define quantlib_credit_default_swap_i

%include instruments.i
%include credit.i
%include termstructures.i
%include bonds.i
%include null.i
%include options.i
%include exercise.i

%{
using QuantLib::CreditDefaultSwap;
using QuantLib::MidPointCdsEngine;
using QuantLib::IntegralCdsEngine;
using QuantLib::IsdaCdsEngine;
using QuantLib::Claim;
using QuantLib::FaceValueClaim;
using QuantLib::FaceValueAccrualClaim;
using QuantLib::CdsOption;
using QuantLib::BlackCdsOptionEngine;
%}

%shared_ptr(Claim);
class Claim {
  private:
    Claim();
  public:
    Real amount(const Date& defaultDate,
                Real notional,
                Real recoveryRate) const;
};

%shared_ptr(FaceValueClaim)
class FaceValueClaim : public Claim {
  public:
    FaceValueClaim();
};

%shared_ptr(FaceValueAccrualClaim)
class FaceValueAccrualClaim : public Claim {
  public:
    FaceValueAccrualClaim(const ext::shared_ptr<Bond>& bond);
};


%shared_ptr(CreditDefaultSwap)
class CreditDefaultSwap : public Instrument {
  public:
    enum PricingModel {
        Midpoint,
        ISDA
    };

    CreditDefaultSwap(Protection::Side side,
                         Real notional,
                         Rate spread,
                         const Schedule& schedule,
                         BusinessDayConvention paymentConvention,
                         const DayCounter& dayCounter,
                         bool settlesAccrual = true,
                         bool paysAtDefaultTime = true,
                         const Date& protectionStart = Date());
    CreditDefaultSwap(Protection::Side side,
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
                         const ext::shared_ptr<Claim>& claim =
                                                    ext::shared_ptr<Claim>(),
                         const DayCounter& lastPeriodDayCounter = DayCounter(),
                         const bool rebatesAccrual = true);
    Protection::Side side() const;
    Real notional() const;
    Rate runningSpread() const;
    %extend {
    doubleOrNull upfront() const {
            boost::optional<Rate> result =
                self->upfront();
            if (result)
                return *result;
            else
                return Null<double>();
        }
    }
    bool settlesAccrual() const;
    bool paysAtDefaultTime() const;
    std::vector<ext::shared_ptr<CashFlow> > coupons();
    const Date& protectionStartDate() const;
    const Date& protectionEndDate() const;
    bool rebatesAccrual() const;
    ext::shared_ptr<CashFlow> upfrontPayment() const;
    ext::shared_ptr<CashFlow> accrualRebate() const;
    const Date& tradeDate() const;
    Natural cashSettlementDays() const;
    Rate fairUpfront() const;
    Rate fairSpread() const;
    Real couponLegBPS() const;
    Real upfrontBPS() const;
    Real couponLegNPV() const;
    Real defaultLegNPV() const;
    Real upfrontNPV() const;
    Real accrualRebateNPV() const;
    Rate impliedHazardRate(Real targetNPV,
                           const Handle<YieldTermStructure>& discountCurve,
                           const DayCounter& dayCounter,
                           Real recoveryRate = 0.4,
                           Real accuracy = 1.0e-6,
                           CreditDefaultSwap::PricingModel model = CreditDefaultSwap::Midpoint) const;
    Rate conventionalSpread(Real conventionalRecovery,
                            const Handle<YieldTermStructure>& discountCurve,
                            const DayCounter& dayCounter,
                            CreditDefaultSwap::PricingModel model = CreditDefaultSwap::Midpoint) const;
};


%shared_ptr(MidPointCdsEngine)
class MidPointCdsEngine : public PricingEngine {
  public:
    MidPointCdsEngine(const Handle<DefaultProbabilityTermStructure>& probability,
                      Real recoveryRate,
                      const Handle<YieldTermStructure>& discountCurve);
};

%shared_ptr(IntegralCdsEngine)
class IntegralCdsEngine : public PricingEngine {
  public:
    IntegralCdsEngine(const Period &integrationStep,
                      const Handle<DefaultProbabilityTermStructure>& probability,
                      Real recoveryRate,
                      const Handle<YieldTermStructure>& discountCurve,
                      bool includeSettlementDateFlows = false);
};

%shared_ptr(IsdaCdsEngine)
class IsdaCdsEngine : public PricingEngine {
    #if defined(SWIGPYTHON)
    %rename(NoFix) None;
    #endif
  public:
    enum NumericalFix {None, Taylor};
    enum AccrualBias {HalfDayBias, NoBias};
    enum ForwardsInCouponPeriod {Flat, Piecewise};
    IsdaCdsEngine(
            const Handle<DefaultProbabilityTermStructure> &probability,
            Real recoveryRate,
            const Handle<YieldTermStructure> &discountCurve,
            bool includeSettlementDateFlows = false,
            const IsdaCdsEngine::NumericalFix numericalFix = IsdaCdsEngine::Taylor,
            const IsdaCdsEngine::AccrualBias accrualBias = IsdaCdsEngine::HalfDayBias,
            const IsdaCdsEngine::ForwardsInCouponPeriod forwardsInCouponPeriod = IsdaCdsEngine::Piecewise);
};

%shared_ptr(CdsOption)
class CdsOption : public Option {
    public:
        CdsOption(const ext::shared_ptr<CreditDefaultSwap>& swap,
                  const ext::shared_ptr<Exercise>& exercise,
                  bool knocksOut = true);
        Rate atmRate() const;
        Real riskyAnnuity() const;
        Volatility impliedVolatility(
                              Real price,
                              const Handle<YieldTermStructure>& termStructure,
                              const Handle<DefaultProbabilityTermStructure>&,
                              Real recoveryRate,
                              Real accuracy = 1.e-4,
                              Size maxEvaluations = 100,
                              Volatility minVol = 1.0e-7,
                              Volatility maxVol = 4.0) const;
};

%shared_ptr(BlackCdsOptionEngine)
class BlackCdsOptionEngine : public PricingEngine {
    public:
        BlackCdsOptionEngine(const Handle<DefaultProbabilityTermStructure>&,
                             Real recoveryRate,
                             const Handle<YieldTermStructure>& termStructure,
                             const Handle<Quote>& vol);
        Handle<YieldTermStructure> termStructure();
        Handle<Quote> volatility();
};

#endif
