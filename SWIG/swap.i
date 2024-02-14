/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2007 StatPro Italia srl
 Copyright (C) 2011 Lluis Pujol Bajador
 Copyright (C) 2015 Gouthaman Balaraman
 Copyright (C) 2016 Peter Caspers
 Copyright (C) 2017, 2018, 2019 Matthias Lungwitz
 Copyright (C) 2018 Matthias Groncki
 Copyright (C) 2023 Marcin Rybacki

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

#ifndef quantlib_swap_i
#define quantlib_swap_i

%include instruments.i
%include termstructures.i
%include cashflows.i
%include timebasket.i
%include indexes.i
%include bonds.i

%{
using QuantLib::Swap;
using QuantLib::VanillaSwap;
using QuantLib::MakeVanillaSwap;
using QuantLib::NonstandardSwap;
using QuantLib::DiscountingSwapEngine;
using QuantLib::FloatFloatSwap;
using QuantLib::OvernightIndexedSwap;
using QuantLib::MakeOIS;
using QuantLib::ZeroCouponSwap;
using QuantLib::EquityTotalReturnSwap;
using QuantLib::ArithmeticAverageOIS;
using QuantLib::simplifyNotificationGraph;
%}

%shared_ptr(Swap)
class Swap : public Instrument {
    %warnfilter(509) Swap;
  public:
    enum Type { Receiver = -1, Payer = 1 };
    Swap(const std::vector<ext::shared_ptr<CashFlow> >& firstLeg,
         const std::vector<ext::shared_ptr<CashFlow> >& secondLeg);
    Swap(const std::vector<Leg>& legs,
         const std::vector<bool>& payer);
    Size numberOfLegs() const;
    Date startDate() const;
    Date maturityDate() const;
    const Leg & leg(Size i);
    Real legNPV(Size j) const;
    Real legBPS(Size k) const;
    DiscountFactor startDiscounts(Size j) const;
    DiscountFactor endDiscounts(Size j) const;
    DiscountFactor npvDateDiscount() const;
    bool payer(Size j) const;
};

void simplifyNotificationGraph(Swap& swap, bool unregisterCoupons = false);


%shared_ptr(VanillaSwap)
class VanillaSwap : public Swap {
  public:
    %extend {
        VanillaSwap(Type type, Real nominal,
                    const Schedule& fixedSchedule, Rate fixedRate,
                    const DayCounter& fixedDayCount,
                    const Schedule& floatSchedule,
                    const ext::shared_ptr<IborIndex>& index,
                    Spread spread,
                    const DayCounter& floatingDayCount,
                    ext::optional<bool> withIndexedCoupons = ext::nullopt) {
            // work around the lack of typemap for this argument
            ext::optional<BusinessDayConvention> paymentConvention = ext::nullopt;

            return new VanillaSwap(type, nominal, fixedSchedule, fixedRate, fixedDayCount,
                                   floatSchedule, index, spread, floatingDayCount,
                                   paymentConvention, withIndexedCoupons);
        }
    }
    Type type() const;
    Rate fairRate();
    Spread fairSpread();
    Real fixedLegBPS();
    Real floatingLegBPS();
    Real fixedLegNPV();
    Real floatingLegNPV();
    // Inspectors 
    const Leg& fixedLeg();
    const Leg& floatingLeg();
    Real nominal();
    const Schedule& fixedSchedule();
    const Schedule& floatingSchedule();
    Rate fixedRate();
    Spread spread();
    const DayCounter& floatingDayCount();
    const DayCounter& fixedDayCount();
};

#if defined(SWIGPYTHON)
%rename (_MakeVanillaSwap) MakeVanillaSwap;
#endif
class MakeVanillaSwap {
      public:
        MakeVanillaSwap& receiveFixed(bool flag = true);
        MakeVanillaSwap& withType(Swap::Type type);
        MakeVanillaSwap& withNominal(Real n);

        MakeVanillaSwap& withSettlementDays(Natural settlementDays);
        MakeVanillaSwap& withEffectiveDate(const Date&);
        MakeVanillaSwap& withTerminationDate(const Date&);
        MakeVanillaSwap& withRule(DateGeneration::Rule r);
        MakeVanillaSwap& withPaymentConvention(BusinessDayConvention bdc);

        MakeVanillaSwap& withFixedLegTenor(const Period& t);
        MakeVanillaSwap& withFixedLegCalendar(const Calendar& cal);
        MakeVanillaSwap& withFixedLegConvention(BusinessDayConvention bdc);
        MakeVanillaSwap& withFixedLegTerminationDateConvention(
                                                   BusinessDayConvention bdc);
        MakeVanillaSwap& withFixedLegRule(DateGeneration::Rule r);
        MakeVanillaSwap& withFixedLegEndOfMonth(bool flag = true);
        MakeVanillaSwap& withFixedLegFirstDate(const Date& d);
        MakeVanillaSwap& withFixedLegNextToLastDate(const Date& d);
        MakeVanillaSwap& withFixedLegDayCount(const DayCounter& dc);

        MakeVanillaSwap& withFloatingLegTenor(const Period& t);
        MakeVanillaSwap& withFloatingLegCalendar(const Calendar& cal);
        MakeVanillaSwap& withFloatingLegConvention(BusinessDayConvention bdc);
        MakeVanillaSwap& withFloatingLegTerminationDateConvention(
                                                   BusinessDayConvention bdc);
        MakeVanillaSwap& withFloatingLegRule(DateGeneration::Rule r);
        MakeVanillaSwap& withFloatingLegEndOfMonth(bool flag = true);
        MakeVanillaSwap& withFloatingLegFirstDate(const Date& d);
        MakeVanillaSwap& withFloatingLegNextToLastDate(const Date& d);
        MakeVanillaSwap& withFloatingLegDayCount(const DayCounter& dc);
        MakeVanillaSwap& withFloatingLegSpread(Spread sp);

        MakeVanillaSwap& withDiscountingTermStructure(
                              const Handle<YieldTermStructure>& discountCurve);
        MakeVanillaSwap& withPricingEngine(
                              const ext::shared_ptr<PricingEngine>& engine);

        MakeVanillaSwap& withIndexedCoupons(bool flag = true);
        MakeVanillaSwap& withAtParCoupons(bool flag = true);

        MakeVanillaSwap(const Period& swapTenor,
                        const ext::shared_ptr<IborIndex>& index,
                        Rate fixedRate,
                        const Period& forwardStart);
        
        %extend {
            ext::shared_ptr<VanillaSwap> makeVanillaSwap() {
                return (ext::shared_ptr<VanillaSwap>)(* $self);
            }
        }
};

#if defined(SWIGPYTHON)
%pythoncode {
_MAKEVANILLA_METHODS = {
    "receiveFixed": "receiveFixed",
    "swapType": "withType",
    "Nominal": "withNominal",
    "settlementDays": "withSettlementDays",
    "effectiveDate": "withEffectiveDate",
    "terminationDate": "withTerminationDate",
    "dateGenerationRule": "withRule",
    "paymentConvention": "withPaymentConvention",
    "fixedLegTenor": "withFixedLegTenor",
    "fixedLegCalendar": "withFixedLegCalendar",
    "fixedLegConvention": "withFixedLegConvention",
    "fixedLegTerminationDateConvention": "withFixedLegTerminationDateConvention",
    "fixedLegDateGenRule": "withFixedLegRule",
    "fixedLegEndOfMonth": "withFixedLegEndOfMonth",
    "fixedLegFirstDate": "withFixedLegFirstDate",
    "fixedLegNextToLastDate": "withFixedLegNextToLastDate",
    "fixedLegDayCount": "withFixedLegDayCount",
    "floatingLegTenor": "withFloatingLegTenor",
    "floatingLegCalendar": "withFloatingLegCalendar",
    "floatingLegConvention": "withFloatingLegConvention",
    "floatingLegTerminationDateConvention": "withFloatingLegTerminationDateConvention",
    "floatingLegDateGenRule": "withFloatingLegRule",
    "floatingLegEndOfMonth": "withFloatingLegEndOfMonth",
    "floatingLegFirstDate": "withFloatingLegFirstDate",
    "floatingLegNextToLastDate": "withFloatingLegNextToLastDate",
    "floatingLegDayCount": "withFloatingLegDayCount",
    "floatingLegSpread": "withFloatingLegSpread",
    "discountingTermStructure": "withDiscountingTermStructure",
    "pricingEngine": "withPricingEngine",
    "withIndexedCoupons": "withIndexedCoupons",
    "atParCoupons": "withAtParCoupons",
}

def MakeVanillaSwap(swapTenor, iborIndex, fixedRate, forwardStart, **kwargs):
    mv = _MakeVanillaSwap(swapTenor, iborIndex, fixedRate, forwardStart)
    _SetSwapAttrs("MakeVanillaSwap", _MAKEVANILLA_METHODS, mv, kwargs)
    return mv.makeVanillaSwap()

def _SetSwapAttrs(func_name, method_map, mv, attrs):
    for name, value in attrs.items():
        try:
            method = method_map[name]
        except KeyError:
            raise TypeError(f"{func_name}() got an unexpected keyword argument {name!r}") from None
        if value is not None:
            getattr(mv, method)(value)
}
#endif

%shared_ptr(NonstandardSwap)
class NonstandardSwap : public Swap {
  public:
    NonstandardSwap(Type type,
                    const std::vector<Real> &fixedNominal,
                    const std::vector<Real> &floatingNominal,
                    const Schedule &fixedSchedule,
                    const std::vector<Real> &fixedRate,
                    const DayCounter &fixedDayCount,
                    const Schedule &floatSchedule,
                    const ext::shared_ptr<IborIndex> &index,
                    const std::vector<Real> &gearing,
                    const std::vector<Spread> &spread,
                    const DayCounter &floatDayCount,
                    const bool intermediateCapitalExchange = false,
                    const bool finalCapitalExchange = false,
                    BusinessDayConvention paymentConvention = Following);
    // Inspectors
    Type type() const;
    const std::vector<Real> &fixedNominal() const;
    const std::vector<Real> &floatingNominal() const;

    const Schedule &fixedSchedule() const;
    const std::vector<Real> &fixedRate() const;
    const DayCounter &fixedDayCount() const;

    const Schedule &floatingSchedule() const;
    const ext::shared_ptr<IborIndex> &iborIndex() const;
    Spread spread() const;
    Real gearing() const;
    const std::vector<Spread>& spreads() const;
    const std::vector<Real>& gearings() const;
    const DayCounter &floatingDayCount() const;

    BusinessDayConvention paymentConvention() const;

    const Leg &fixedLeg() const;
    const Leg &floatingLeg() const;
};

%shared_ptr(DiscountingSwapEngine)
class DiscountingSwapEngine : public PricingEngine {
  public:
    DiscountingSwapEngine(const Handle<YieldTermStructure>& discountCurve,
                          bool includeSettlementDateFlows,
                          const Date& settlementDate = Date(),
                          const Date& npvDate = Date());
    %extend {
        DiscountingSwapEngine(const Handle<YieldTermStructure>& discountCurve,
                              const Date& settlementDate = Date(),
                              const Date& npvDate = Date()) {
            return new DiscountingSwapEngine(discountCurve,
                                             ext::nullopt,
                                             settlementDate,
                                             npvDate);
        }
    }
};


%{
using QuantLib::AssetSwap;
using QuantLib::OvernightIndexedSwapIndex;
%}

%shared_ptr(AssetSwap)
class AssetSwap : public Swap {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") AssetSwap;
    #endif
  public:
    AssetSwap(bool payFixedRate,
                 const ext::shared_ptr<Bond>& bond,
                 Real bondCleanPrice,
                 const ext::shared_ptr<IborIndex>& index,
                 Spread spread,
                 const Schedule& floatSchedule = Schedule(),
                 const DayCounter& floatingDayCount = DayCounter(),
                 bool parAssetSwap = true);
    Real fairCleanPrice();
    Spread fairSpread();
};

%shared_ptr(FloatFloatSwap)
class FloatFloatSwap : public Swap {
  public:
    FloatFloatSwap(Type type, const std::vector<Real> &nominal1,
        const std::vector<Real> &nominal2, const Schedule &schedule1,
        const ext::shared_ptr<InterestRateIndex> &index1,
        const DayCounter &dayCount1, const Schedule &schedule2,
        const ext::shared_ptr<InterestRateIndex> &index2,
        const DayCounter &dayCount2,
        const bool intermediateCapitalExchange = false,
        const bool finalCapitalExchange = false,
        const std::vector<Real> &gearing1 = std::vector<Real>(),
        const std::vector<Real> &spread1 = std::vector<Real>(),
        const std::vector<Real> &cappedRate1 = std::vector<Real>(),
        const std::vector<Real> &flooredRate1 = std::vector<Real>(),
        const std::vector<Real> &gearing2 = std::vector<Real>(),
        const std::vector<Real> &spread2 = std::vector<Real>(),
        const std::vector<Real> &cappedRate2 = std::vector<Real>(),
        const std::vector<Real> &flooredRate2 = std::vector<Real>(),
        BusinessDayConvention paymentConvention1 = Following,
        BusinessDayConvention paymentConvention2 = Following);
};

%shared_ptr(OvernightIndexedSwap)
class OvernightIndexedSwap : public Swap {
  public:
    OvernightIndexedSwap(
            Type type,
            Real nominal,
            const Schedule& schedule,
            Rate fixedRate,
            const DayCounter& fixedDC,
            const ext::shared_ptr<OvernightIndex>& index,
            Spread spread = 0.0,
            Integer paymentLag = 0,
            BusinessDayConvention paymentAdjustment = Following,
            Calendar paymentCalendar = Calendar(),
            bool telescopicValueDates = false,
            RateAveraging::Type averagingMethod = RateAveraging::Compound);
    
    OvernightIndexedSwap(
            Type type,
            std::vector<Real> nominals,
            const Schedule& schedule,
            Rate fixedRate,
            const DayCounter& fixedDC,
            const ext::shared_ptr<OvernightIndex>& index,
            Spread spread = 0.0,
            Integer paymentLag = 0,
            BusinessDayConvention paymentAdjustment = Following,
            Calendar paymentCalendar = Calendar(),
            bool telescopicValueDates = false,
            RateAveraging::Type averagingMethod = RateAveraging::Compound);

    OvernightIndexedSwap(Type type,
                         const std::vector<Real>& fixedNominals,
                         const Schedule& fixedSchedule,
                         Rate fixedRate,
                         const DayCounter& fixedDC,
                         const std::vector<Real>& overnightNominals,
                         const Schedule& overnightSchedule,
                         const ext::shared_ptr<OvernightIndex>& overnightIndex,
                         Spread spread = 0.0,
                         Integer paymentLag = 0,
                         BusinessDayConvention paymentAdjustment = Following,
                         const Calendar& paymentCalendar = Calendar(),
                         bool telescopicValueDates = false,
                         RateAveraging::Type averagingMethod = RateAveraging::Compound);

    Rate fixedLegBPS();
    Real fixedLegNPV();
    Real fairRate();
    Real overnightLegBPS();
    Real overnightLegNPV();
    Spread fairSpread();
    // Inspectors
    Type type();
    Real nominal();
    std::vector<Real> nominals();
    Frequency paymentFrequency();
    Rate fixedRate();
    const DayCounter& fixedDayCount();
    ext::shared_ptr<OvernightIndex> overnightIndex() const;
    Spread spread();
    const Leg& fixedLeg();
    const Leg& overnightLeg();
    RateAveraging::Type averagingMethod();
};

#if defined(SWIGPYTHON)
%rename (_MakeOIS) MakeOIS;
#endif
class MakeOIS {
      public:
        MakeOIS(const Period& swapTenor,
                const ext::shared_ptr<OvernightIndex>& overnightIndex,
                Rate fixedRate = Null<Rate>(),
                const Period& fwdStart = 0*Days);

        %extend {
            ext::shared_ptr<OvernightIndexedSwap> makeOIS() {
                return (ext::shared_ptr<OvernightIndexedSwap>)(* $self);
            }
        }

        MakeOIS& receiveFixed(bool flag = true);
        MakeOIS& withType(Swap::Type type);
        MakeOIS& withNominal(Real n);
        MakeOIS& withSettlementDays(Natural settlementDays);
        MakeOIS& withEffectiveDate(const Date&);
        MakeOIS& withTerminationDate(const Date&);
        MakeOIS& withRule(DateGeneration::Rule r);
        MakeOIS& withFixedLegRule(DateGeneration::Rule r);
        MakeOIS& withOvernightLegRule(DateGeneration::Rule r);
        MakeOIS& withPaymentFrequency(Frequency f);
        MakeOIS& withFixedLegPaymentFrequency(Frequency f);
        MakeOIS& withOvernightLegPaymentFrequency(Frequency f);
        MakeOIS& withPaymentAdjustment(BusinessDayConvention convention);
        MakeOIS& withPaymentLag(Integer lag);
        MakeOIS& withPaymentCalendar(const Calendar& cal);
        MakeOIS& withCalendar(const Calendar& cal);
        MakeOIS& withFixedLegCalendar(const Calendar& cal);
        MakeOIS& withOvernightLegCalendar(const Calendar& cal);
        MakeOIS& withConvention(BusinessDayConvention bdc);
        MakeOIS& withFixedLegConvention(BusinessDayConvention bdc);
        MakeOIS& withOvernightLegConvention(BusinessDayConvention bdc);
        MakeOIS& withTerminationDateConvention(BusinessDayConvention bdc);
        MakeOIS& withFixedLegTerminationDateConvention(BusinessDayConvention bdc);
        MakeOIS& withOvernightLegTerminationDateConvention(BusinessDayConvention bdc);
        MakeOIS& withEndOfMonth(bool flag = true);
        MakeOIS& withFixedLegEndOfMonth(bool flag = true);
        MakeOIS& withOvernightLegEndOfMonth(bool flag = true);
        MakeOIS& withFixedLegDayCount(const DayCounter& dc);
        MakeOIS& withOvernightLegSpread(Spread sp);
        MakeOIS& withDiscountingTermStructure(
                  const Handle<YieldTermStructure>& discountingTermStructure);
        MakeOIS& withTelescopicValueDates(bool telescopicValueDates);
        MakeOIS& withAveragingMethod(RateAveraging::Type averagingMethod);
        MakeOIS& withPricingEngine(
                              const ext::shared_ptr<PricingEngine>& engine);
};

#if defined(SWIGPYTHON)
%pythoncode {
_MAKEOIS_METHODS = {
    "receiveFixed": "receiveFixed",
    "swapType": "withType",
    "nominal": "withNominal",
    "settlementDays": "withSettlementDays",
    "effectiveDate": "withEffectiveDate",
    "terminationDate": "withTerminationDate",
    "dateGenerationRule": "withRule",
    "fixedLegRule": "withFixedLegRule",
    "overnightLegRule": "withOvernightLegRule",
    "paymentFrequency": "withPaymentFrequency",
    "fixedLegPaymentFrequency": "withFixedLegPaymentFrequency",
    "overnightLegPaymentFrequency": "withOvernightLegPaymentFrequency",
    "paymentAdjustmentConvention": "withPaymentAdjustment",
    "paymentLag": "withPaymentLag",
    "paymentCalendar": "withPaymentCalendar",
    "calendar": "withCalendar",
    "fixedLegCalendar": "withFixedLegCalendar",
    "overnightLegCalendar": "withOvernightLegCalendar",
    "convention": "withConvention",
    "fixedLegConvention": "withFixedLegConvention",
    "overnightLegConvention": "withOvernightLegConvention",
    "terminationDateConvention": "withTerminationDateConvention",
    "fixedLegTerminationDateConvention": "withFixedLegTerminationDateConvention",
    "overnightLegTerminationDateConvention": "withOvernightLegTerminationDateConvention",
    "endOfMonth": "withEndOfMonth",
    "fixedLegEndOfMonth": "withFixedLegEndOfMonth",
    "overnightLegEndOfMonth": "withOvernightLegEndOfMonth",
    "fixedLegDayCount": "withFixedLegDayCount",
    "overnightLegSpread": "withOvernightLegSpread",
    "discountingTermStructure": "withDiscountingTermStructure",
    "telescopicValueDates": "withTelescopicValueDates",
    "averagingMethod": "withAveragingMethod",
    "pricingEngine": "withPricingEngine",
}

def MakeOIS(swapTenor, overnightIndex, fixedRate, fwdStart=Period(0, Days), **kwargs):
    mv = _MakeOIS(swapTenor, overnightIndex, fixedRate, fwdStart)
    _SetSwapAttrs("MakeOIS", _MAKEOIS_METHODS, mv, kwargs)
    return mv.makeOIS()
}
#endif


%shared_ptr(OvernightIndexedSwapIndex)
class OvernightIndexedSwapIndex : public SwapIndex {
  public:
    OvernightIndexedSwapIndex(
              const std::string& familyName,
              const Period& tenor,
              Natural settlementDays,
              Currency currency,
              const ext::shared_ptr<OvernightIndex>& overnightIndex,
              bool telescopicValueDates = false,
              RateAveraging::Type averagingMethod = RateAveraging::Compound);
    //! \name Inspectors
    //@{
    ext::shared_ptr<OvernightIndex> overnightIndex() const;
    /*! \warning Relinking the term structure underlying the index will
                 not have effect on the returned swap.
    */
    ext::shared_ptr<OvernightIndexedSwap> underlyingSwap(
                                            const Date& fixingDate) const;
};

%inline %{
    ext::shared_ptr<OvernightIndexedSwap> as_overnight_swap_index(
                          const ext::shared_ptr<InterestRateIndex>& index) {
        return ext::dynamic_pointer_cast<OvernightIndexedSwap>(index);
    }
%}

%shared_ptr(ZeroCouponSwap)
class ZeroCouponSwap : public Swap {
  public:
    ZeroCouponSwap(Type type,
                   Real baseNominal,
                   const Date& startDate,
                   const Date& maturityDate, 
                   Real fixedPayment,
                   ext::shared_ptr<IborIndex> iborIndex,
                   const Calendar& paymentCalendar,
                   BusinessDayConvention paymentConvention = Following,
                   Natural paymentDelay = 0);

    ZeroCouponSwap(Type type,
                   Real baseNominal,
                   const Date& startDate,
                   const Date& maturityDate,
                   Rate fixedRate,
                   const DayCounter& fixedDayCounter,
                   ext::shared_ptr<IborIndex> iborIndex,
                   const Calendar& paymentCalendar,
                   BusinessDayConvention paymentConvention = Following,
                   Natural paymentDelay = 0);

    // Inspectors
    Type type() const;
    Real baseNominal() const;
    const ext::shared_ptr<IborIndex>& iborIndex() const;
    
    const Leg& fixedLeg() const;
    const Leg& floatingLeg() const;
    Real fixedPayment() const;

    Real fixedLegNPV() const;
    Real floatingLegNPV() const;
    Real fairFixedPayment() const;
    Rate fairFixedRate(const DayCounter& dayCounter) const;
};

%shared_ptr(EquityTotalReturnSwap)
class EquityTotalReturnSwap : public Swap {
  public:
    EquityTotalReturnSwap(Type type,
                          Real nominal,
                          Schedule schedule,
                          ext::shared_ptr<EquityIndex> equityIndex,
                          const ext::shared_ptr<IborIndex>& interestRateIndex,
                          DayCounter dayCounter,
                          Rate margin,
                          Real gearing = 1.0,
                          Calendar paymentCalendar = Calendar(),
                          BusinessDayConvention paymentConvention = Unadjusted,
                          Natural paymentDelay = 0);

    EquityTotalReturnSwap(Type type,
                          Real nominal,
                          Schedule schedule,
                          ext::shared_ptr<EquityIndex> equityIndex,
                          const ext::shared_ptr<OvernightIndex>& interestRateIndex,
                          DayCounter dayCounter,
                          Rate margin,
                          Real gearing = 1.0,
                          Calendar paymentCalendar = Calendar(),
                          BusinessDayConvention paymentConvention = Unadjusted,
                          Natural paymentDelay = 0);

    // Inspectors
    Type type() const;
    Real nominal() const;
    
    const ext::shared_ptr<EquityIndex>& equityIndex() const;
    const ext::shared_ptr<InterestRateIndex>& interestRateIndex() const;
    
    const Schedule& schedule() const;
    const DayCounter& dayCounter() const;
    Rate margin() const;
    Real gearing() const;
    const Calendar& paymentCalendar() const;
    BusinessDayConvention paymentConvention() const;
    Natural paymentDelay() const;

    const Leg& equityLeg() const;
    const Leg& interestRateLeg() const;

    Real equityLegNPV() const;
    Real interestRateLegNPV() const;
    Real fairMargin() const;
};

%shared_ptr(ArithmeticAverageOIS)
class ArithmeticAverageOIS : public Swap {
  public:
    ArithmeticAverageOIS(Type type,
                         Real nominal,
                         const Schedule& fixedLegSchedule,
                         Rate fixedRate,
                         DayCounter fixedDC,
                         ext::shared_ptr<OvernightIndex> overnightIndex,
                         const Schedule& overnightLegSchedule,
                         Spread spread = 0.0,
                         Real meanReversionSpeed = 0.03,
                         Real volatility = 0.00, // NO convexity adjustment by default
                         bool byApprox = false); // TRUE to use Katsumi Takada approximation
    ArithmeticAverageOIS(Type type,
                         std::vector<Real> nominals,
                         const Schedule& fixedLegSchedule,
                         Rate fixedRate,
                         DayCounter fixedDC,
                         ext::shared_ptr<OvernightIndex> overnightIndex,
                         const Schedule& overnightLegSchedule,
                         Spread spread = 0.0,
                         Real meanReversionSpeed = 0.03,
                         Real volatility = 0.00, // NO convexity adjustment by default
                         bool byApprox = false); // TRUE to use Katsumi Takada approximation

    Type type() const;
    Real nominal() const;
    std::vector<Real> nominals() const;

    Frequency fixedLegPaymentFrequency();
    Frequency overnightLegPaymentFrequency();

    Rate fixedRate() const;
    const DayCounter& fixedDayCount();

    ext::shared_ptr<OvernightIndex> overnightIndex();
    Spread spread() const;

    const Leg& fixedLeg() const;
    const Leg& overnightLeg() const;

    Real fixedLegBPS() const;
    Real fixedLegNPV() const;
    Real fairRate() const;

    Real overnightLegBPS() const;
    Real overnightLegNPV() const;
    Spread fairSpread() const;
};
#endif
