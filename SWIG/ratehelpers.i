
/*
 Copyright (C) 2005, 2006, 2007, 2008 StatPro Italia srl
 Copyright (C) 2009 Joseph Malicki
 Copyright (C) 2018 Matthias Lungwitz
 Copyright (C) 2021 Marcin Rybacki

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

#ifndef quantlib_rate_helpers_i
#define quantlib_rate_helpers_i

%include bonds.i
%include date.i
%include calendars.i
%include daycounters.i
%include futures.i
%include marketelements.i
%include types.i
%include vectors.i
%include swap.i

%{
using QuantLib::Pillar;
using QuantLib::RateHelper;
using QuantLib::DepositRateHelper;
using QuantLib::FraRateHelper;
using QuantLib::FuturesRateHelper;
using QuantLib::SwapRateHelper;
using QuantLib::BondHelper;
using QuantLib::FixedRateBondHelper;
using QuantLib::OISRateHelper;
using QuantLib::DatedOISRateHelper;
using QuantLib::FxSwapRateHelper;
using QuantLib::OvernightIndexFutureRateHelper;
using QuantLib::SofrFutureRateHelper;
using QuantLib::ConstNotionalCrossCurrencyBasisSwapRateHelper;
using QuantLib::MtMCrossCurrencyBasisSwapRateHelper;
using QuantLib::IborIborBasisSwapRateHelper;
using QuantLib::OvernightIborBasisSwapRateHelper;
using QuantLib::ArithmeticOISRateHelper;
%}

struct Pillar {
    enum Choice { MaturityDate, LastRelevantDate, CustomDate};
};

%shared_ptr(RateHelper)
class RateHelper : public Observable {
  public:
    const Handle<Quote>& quote() const;
    Date latestDate() const;
	Date earliestDate() const;
	Date maturityDate() const;
	Date latestRelevantDate() const;
	Date pillarDate() const;
	Real impliedQuote() const;
	Real quoteError() const;
  private:
    RateHelper();
};

%shared_ptr(DepositRateHelper)
class DepositRateHelper : public RateHelper {
  public:
    DepositRateHelper(
            const Handle<Quote>& rate,
            const Period& tenor,
            Natural fixingDays,
            const Calendar& calendar,
            BusinessDayConvention convention,
            bool endOfMonth,
            const DayCounter& dayCounter);
    DepositRateHelper(
            Rate rate,
            const Period& tenor,
            Natural fixingDays,
            const Calendar& calendar,
            BusinessDayConvention convention,
            bool endOfMonth,
            const DayCounter& dayCounter);
    DepositRateHelper(const Handle<Quote>& rate,
                         const ext::shared_ptr<IborIndex>& index);
    DepositRateHelper(Rate rate,
                         const ext::shared_ptr<IborIndex>& index);
};

%shared_ptr(FraRateHelper)
class FraRateHelper : public RateHelper {
  public:
    FraRateHelper(
            const Handle<Quote>& rate,
            Natural monthsToStart,
            Natural monthsToEnd,
            Natural fixingDays,
            const Calendar& calendar,
            BusinessDayConvention convention,
            bool endOfMonth,
            const DayCounter& dayCounter,
            Pillar::Choice pillar = Pillar::LastRelevantDate,
            Date customPillarDate = Date(),
            bool useIndexedCoupon = true);
    FraRateHelper(
            Rate rate,
            Natural monthsToStart,
            Natural monthsToEnd,
            Natural fixingDays,
            const Calendar& calendar,
            BusinessDayConvention convention,
            bool endOfMonth,
            const DayCounter& dayCounter,
            Pillar::Choice pillar = Pillar::LastRelevantDate,
            Date customPillarDate = Date(),
            bool useIndexedCoupon = true);
    FraRateHelper(const Handle<Quote>& rate,
                  Natural monthsToStart,
                  const ext::shared_ptr<IborIndex>& index,
                  Pillar::Choice pillar = Pillar::LastRelevantDate,
                  Date customPillarDate = Date(),
                  bool useIndexedCoupon = true);
    FraRateHelper(Rate rate,
                  Natural monthsToStart,
                  const ext::shared_ptr<IborIndex>& index,
                  Pillar::Choice pillar = Pillar::LastRelevantDate,
                  Date customPillarDate = Date(),
                  bool useIndexedCoupon = true);
    FraRateHelper(const Handle<Quote>& rate,
                  Natural immOffsetStart,
                  Natural immOffsetEnd,
                  const ext::shared_ptr<IborIndex>& iborIndex,
                  Pillar::Choice pillar = Pillar::LastRelevantDate,
                  Date customPillarDate = Date(),
                  bool useIndexedCoupon = true);
    FraRateHelper(Rate rate,
                  Natural immOffsetStart,
                  Natural immOffsetEnd,
                  const ext::shared_ptr<IborIndex>& iborIndex,
                  Pillar::Choice pillar = Pillar::LastRelevantDate,
                  Date customPillarDate = Date(),
                  bool useIndexedCoupon = true);
    FraRateHelper(Rate rate,
                  Period periodToStart,
                  const ext::shared_ptr<IborIndex>& iborIndex,
                  Pillar::Choice pillar = Pillar::LastRelevantDate,
                  Date customPillarDate = Date(),
                  bool useIndexedCoupon = true);
    FraRateHelper(const Handle<Quote>& rate,
                  Period periodToStart,
                  const ext::shared_ptr<IborIndex>& iborIndex,
                  Pillar::Choice pillar = Pillar::LastRelevantDate,
                  Date customPillarDate = Date(),
                  bool useIndexedCoupon = true);
};

%shared_ptr(FuturesRateHelper)
class FuturesRateHelper : public RateHelper {
  public:
    FuturesRateHelper(
            const Handle<Quote>& price,
            const Date& iborStartDate,
            Natural nMonths,
            const Calendar& calendar,
            BusinessDayConvention convention,
            bool endOfMonth,
            const DayCounter& dayCounter,
            const Handle<Quote>& convexityAdjustment = Handle<Quote>(),
            Futures::Type type = Futures::IMM);
    FuturesRateHelper(
            Real price,
            const Date& iborStartDate,
            Natural nMonths,
            const Calendar& calendar,
            BusinessDayConvention convention,
            bool endOfMonth,
            const DayCounter& dayCounter,
            Rate convexityAdjustment = 0.0,
            Futures::Type type = Futures::IMM);
    FuturesRateHelper(
            const Handle<Quote>& price,
            const Date& iborStartDate,
            const Date& iborEndDate,
            const DayCounter& dayCounter,
            const Handle<Quote>& convexityAdjustment = Handle<Quote>(),
            Futures::Type type = Futures::IMM);
    FuturesRateHelper(
            Real price,
            const Date& iborStartDate,
            const Date& iborEndDate,
            const DayCounter& dayCounter,
            Rate convexityAdjustment = 0.0,
            Futures::Type type = Futures::IMM);
    FuturesRateHelper(
            const Handle<Quote>& price,
            const Date& iborStartDate,
            const ext::shared_ptr<IborIndex>& index,
            const Handle<Quote>& convexityAdjustment = Handle<Quote>(),
            Futures::Type type = Futures::IMM);
    FuturesRateHelper(
            Real price,
            const Date& iborStartDate,
            const ext::shared_ptr<IborIndex>& index,
            Real convexityAdjustment = 0.0,
            Futures::Type type = Futures::IMM);
};

%shared_ptr(SwapRateHelper)
class SwapRateHelper : public RateHelper {
  public:
    SwapRateHelper(
            const Handle<Quote>& rate,
            const Period& tenor,
            const Calendar& calendar,
            Frequency fixedFrequency,
            BusinessDayConvention fixedConvention,
            const DayCounter& fixedDayCount,
            const ext::shared_ptr<IborIndex>& index,
            const Handle<Quote>& spread = Handle<Quote>(),
            const Period& fwdStart = 0*Days,
            const Handle<YieldTermStructure>& discountingCurve = {},
            Natural settlementDays = Null<Natural>(),
            Pillar::Choice pillar = Pillar::LastRelevantDate,
            Date customPillarDate = Date(),
            bool endOfMonth = false,
            ext::optional<bool> withIndexedCoupons = ext::nullopt);
    SwapRateHelper(
            Rate rate,
            const Period& tenor,
            const Calendar& calendar,
            Frequency fixedFrequency,
            BusinessDayConvention fixedConvention,
            const DayCounter& fixedDayCount,
            const ext::shared_ptr<IborIndex>& index,
            const Handle<Quote>& spread = Handle<Quote>(),
            const Period& fwdStart = 0*Days,
            const Handle<YieldTermStructure>& discountingCurve = {},
            Natural settlementDays = Null<Natural>(),
            Pillar::Choice pillar = Pillar::LastRelevantDate,
            Date customPillarDate = Date(),
            bool endOfMonth = false,
            ext::optional<bool> withIndexedCoupons = ext::nullopt);
    SwapRateHelper(
            const Handle<Quote>& rate,
            const ext::shared_ptr<SwapIndex>& index,
            const Handle<Quote>& spread = Handle<Quote>(),
            const Period& fwdStart = 0*Days,
            const Handle<YieldTermStructure>& discountingCurve = {},
            Pillar::Choice pillar = Pillar::LastRelevantDate,
            Date customPillarDate = Date(),
            bool endOfMonth = false,
            ext::optional<bool> withIndexedCoupons = ext::nullopt);
    SwapRateHelper(
            Rate rate,
            const ext::shared_ptr<SwapIndex>& index,
            const Handle<Quote>& spread = Handle<Quote>(),
            const Period& fwdStart = 0*Days,
            const Handle<YieldTermStructure>& discountingCurve = {},
            Pillar::Choice pillar = Pillar::LastRelevantDate,
            Date customPillarDate = Date(),
            bool endOfMonth = false,
            ext::optional<bool> withIndexedCoupons = ext::nullopt);
    Spread spread();
    ext::shared_ptr<VanillaSwap> swap();
};

%shared_ptr(BondHelper)
class BondHelper : public RateHelper {
  public:
    BondHelper(const Handle<Quote>& cleanPrice,
               const ext::shared_ptr<Bond>& bond,
               BondPrice::Type priceType = BondPrice::Clean);

    ext::shared_ptr<Bond> bond();
};

%shared_ptr(FixedRateBondHelper)
class FixedRateBondHelper : public BondHelper {
  public:
    FixedRateBondHelper(
                  const Handle<Quote>& cleanPrice,
                  Size settlementDays,
                  Real faceAmount,
                  const Schedule& schedule,
                  const std::vector<Rate>& coupons,
                  const DayCounter& paymentDayCounter,
                  BusinessDayConvention paymentConvention = Following,
                  Real redemption = 100.0,
                  const Date& issueDate = Date(),
                  const Calendar& paymentCalendar = Calendar(),
                  const Period& exCouponPeriod = Period(),
                  const Calendar& exCouponCalendar = Calendar(),
                  BusinessDayConvention exCouponConvention = Unadjusted,
                  bool exCouponEndOfMonth = false,
                  BondPrice::Type priceType = BondPrice::Clean);

    ext::shared_ptr<FixedRateBond> fixedRateBond();
};


%shared_ptr(OISRateHelper)
class OISRateHelper : public RateHelper {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") OISRateHelper;
    #endif
  public:
    OISRateHelper(
            Natural settlementDays,
            const Period& tenor,
            const Handle<Quote>& rate,
            const ext::shared_ptr<OvernightIndex>& index,
            const Handle<YieldTermStructure>& discountingCurve = {},
            bool telescopicValueDates = false,
            Integer paymentLag = 0,
            BusinessDayConvention paymentConvention = Following,
            Frequency paymentFrequency = Annual,
            const Calendar& paymentCalendar = Calendar(),
            const Period& forwardStart = 0 * Days, 
            const Spread overnightSpread = 0.0,
            Pillar::Choice pillar = Pillar::LastRelevantDate,
            Date customPillarDate = Date(), 
            RateAveraging::Type averagingMethod = RateAveraging::Compound,
            ext::optional<bool> endOfMonth = ext::nullopt,
            ext::optional<Frequency> fixedPaymentFrequency = ext::nullopt,
            const Calendar& fixedCalendar = Calendar());
    ext::shared_ptr<OvernightIndexedSwap> swap();
};

%shared_ptr(DatedOISRateHelper)
class DatedOISRateHelper : public RateHelper {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") DatedOISRateHelper;
    #endif
  public:
    DatedOISRateHelper(
            const Date& startDate,
            const Date& endDate,
            const Handle<Quote>& rate,
            const ext::shared_ptr<OvernightIndex>& index,
            const Handle<YieldTermStructure>& discountingCurve = {},
            bool telescopicValueDates = false, 
            RateAveraging::Type averagingMethod = RateAveraging::Compound,
            Integer paymentLag = 0,
            BusinessDayConvention paymentConvention = Following,
            Frequency paymentFrequency = Annual,
            const Calendar& paymentCalendar = Calendar(),
            const Period& forwardStart = 0 * Days,
            Spread overnightSpread = 0.0,
            ext::optional<bool> endOfMonth = ext::nullopt,
            ext::optional<Frequency> fixedPaymentFrequency = ext::nullopt,
            const Calendar& fixedCalendar = Calendar());
};

%shared_ptr(FxSwapRateHelper)
class FxSwapRateHelper : public RateHelper {
  public:
    FxSwapRateHelper(
            const Handle<Quote>& fwdPoint,
            const Handle<Quote>& spotFx,
            const Period& tenor,
            Natural fixingDays,
            const Calendar& calendar,
            BusinessDayConvention convention,
            bool endOfMonth,
            bool isFxBaseCurrencyCollateralCurrency,
            const Handle<YieldTermStructure>& collateralCurve,
            const Calendar& tradingCalendar = Calendar());
};

%shared_ptr(OvernightIndexFutureRateHelper)
class OvernightIndexFutureRateHelper : public RateHelper {
  public:
    OvernightIndexFutureRateHelper(
            const Handle<Quote>& price,
            const Date& valueDate,
            const Date& maturityDate,
            const ext::shared_ptr<OvernightIndex>& index,
            const Handle<Quote>& convexityAdjustment = Handle<Quote>(), 
            RateAveraging::Type averagingMethod = RateAveraging::Compound);
};

%shared_ptr(SofrFutureRateHelper)
class SofrFutureRateHelper : public OvernightIndexFutureRateHelper {
  public:
    SofrFutureRateHelper(const Handle<Quote>& price,
                         Month referenceMonth,
                         Year referenceYear,
                         Frequency referenceFreq,
                         const Handle<Quote>& convexityAdjustment = Handle<Quote>());

    SofrFutureRateHelper(Real price,
                         Month referenceMonth,
                         Year referenceYear,
                         Frequency referenceFreq,
                         Real convexityAdjustment = 0.0);
};

%shared_ptr(ConstNotionalCrossCurrencyBasisSwapRateHelper)
class ConstNotionalCrossCurrencyBasisSwapRateHelper : public RateHelper {
  public:
    ConstNotionalCrossCurrencyBasisSwapRateHelper(const Handle<Quote>& basis,
                                                  const Period& tenor,
                                                  Natural fixingDays,
                                                  Calendar calendar,
                                                  BusinessDayConvention convention,
                                                  bool endOfMonth,
                                                  ext::shared_ptr<IborIndex> baseCurrencyIndex,
                                                  ext::shared_ptr<IborIndex> quoteCurrencyIndex,
                                                  Handle<YieldTermStructure> collateralCurve,
                                                  bool isFxBaseCurrencyCollateralCurrency,
                                                  bool isBasisOnFxBaseCurrencyLeg);
};

%shared_ptr(MtMCrossCurrencyBasisSwapRateHelper)
class MtMCrossCurrencyBasisSwapRateHelper : public RateHelper {
  public:
    MtMCrossCurrencyBasisSwapRateHelper(const Handle<Quote>& basis,
                                        const Period& tenor,
                                        Natural fixingDays,
                                        Calendar calendar,
                                        BusinessDayConvention convention,
                                        bool endOfMonth,
                                        ext::shared_ptr<IborIndex> baseCurrencyIndex,
                                        ext::shared_ptr<IborIndex> quoteCurrencyIndex,
                                        Handle<YieldTermStructure> collateralCurve,
                                        bool isFxBaseCurrencyCollateralCurrency,
                                        bool isBasisOnFxBaseCurrencyLeg,
                                        bool isFxBaseCurrencyLegResettable);
};

%shared_ptr(IborIborBasisSwapRateHelper)
class IborIborBasisSwapRateHelper : public RateHelper {
  public:
    IborIborBasisSwapRateHelper(const Handle<Quote>& basis,
                                const Period& tenor,
                                Natural settlementDays,
                                Calendar calendar,
                                BusinessDayConvention convention,
                                bool endOfMonth,
                                const ext::shared_ptr<IborIndex>& baseIndex,
                                const ext::shared_ptr<IborIndex>& otherIndex,
                                Handle<YieldTermStructure> discountHandle,
                                bool bootstrapBaseCurve);
};

%shared_ptr(OvernightIborBasisSwapRateHelper)
class OvernightIborBasisSwapRateHelper : public RateHelper {
  public:
    OvernightIborBasisSwapRateHelper(const Handle<Quote>& basis,
                                     const Period& tenor,
                                     Natural settlementDays,
                                     Calendar calendar,
                                     BusinessDayConvention convention,
                                     bool endOfMonth,
                                     const ext::shared_ptr<OvernightIndex>& baseIndex,
                                     const ext::shared_ptr<IborIndex>& otherIndex,
                                     Handle<YieldTermStructure> discountHandle = Handle<YieldTermStructure>());
};

%shared_ptr(ArithmeticOISRateHelper)
class ArithmeticOISRateHelper : public RateHelper {
  public:
    ArithmeticOISRateHelper(
        Natural settlementDays,
        const Period& tenor, // swap maturity
        Frequency fixedLegPaymentFrequency,
        const Handle<Quote>& fixedRate,
        ext::shared_ptr<OvernightIndex> overnightIndex,
        Frequency overnightLegPaymentFrequency,
        Handle<Quote> spread,
        Real meanReversionSpeed = 0.03,
        Real volatility = 0.00, // NO convexity adjustment by default
        bool byApprox = false,  // TRUE to use Katsumi Takada approximation
        // exogenous discounting curve
        Handle<YieldTermStructure> discountingCurve = Handle<YieldTermStructure>());

    ext::shared_ptr<ArithmeticAverageOIS> swap() const;
};

// allow use of RateHelper vectors
#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<RateHelper> )
#endif
namespace std {
    %template(RateHelperVector) vector<ext::shared_ptr<RateHelper> >;
}

// allow use of RateHelper vectors
#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<BondHelper> )
#endif
namespace std {
    %template(BondHelperVector) vector<ext::shared_ptr<BondHelper> >;
}

%inline %{
    const ext::shared_ptr<DepositRateHelper> as_depositratehelper(const ext::shared_ptr<RateHelper> helper) {
        return ext::dynamic_pointer_cast<DepositRateHelper>(helper);
    }
	const ext::shared_ptr<FraRateHelper> as_fraratehelper(const ext::shared_ptr<RateHelper> helper) {
        return ext::dynamic_pointer_cast<FraRateHelper>(helper);
    }
    const ext::shared_ptr<SwapRateHelper> as_swapratehelper(const ext::shared_ptr<RateHelper> helper) {
        return ext::dynamic_pointer_cast<SwapRateHelper>(helper);
    }
    const ext::shared_ptr<OISRateHelper> as_oisratehelper(const ext::shared_ptr<RateHelper> helper) {
        return ext::dynamic_pointer_cast<OISRateHelper>(helper);
    }
    const ext::shared_ptr<ConstNotionalCrossCurrencyBasisSwapRateHelper> as_constnotionalcrosscurrencybasisswapratehelper(
            const ext::shared_ptr<RateHelper> helper) {
        return ext::dynamic_pointer_cast<ConstNotionalCrossCurrencyBasisSwapRateHelper>(helper);
    }
    const ext::shared_ptr<MtMCrossCurrencyBasisSwapRateHelper> as_mtmcrosscurrencybasisswapratehelper(
            const ext::shared_ptr<RateHelper> helper) {
        return ext::dynamic_pointer_cast<MtMCrossCurrencyBasisSwapRateHelper>(helper);
    }
%}

#endif
