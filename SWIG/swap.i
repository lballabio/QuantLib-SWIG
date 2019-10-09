/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2007 StatPro Italia srl
 Copyright (C) 2011 Lluis Pujol Bajador
 Copyright (C) 2015 Gouthaman Balaraman
 Copyright (C) 2016 Peter Caspers
 Copyright (C) 2017, 2018, 2019 Matthias Lungwitz
 Copyright (C) 2018 Matthias Groncki

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
%}

%shared_ptr(Swap)
class Swap : public Instrument {
  public:
    Swap(const std::vector<boost::shared_ptr<CashFlow> >& firstLeg,
         const std::vector<boost::shared_ptr<CashFlow> >& secondLeg);
    Date startDate();
    Date maturityDate();
    const Leg & leg(Size i);
    Real legNPV(Size j) const;
};

%shared_ptr(VanillaSwap)
class VanillaSwap : public Swap {
  public:
    enum Type { Receiver = -1, Payer = 1 };
    VanillaSwap(VanillaSwap::Type type, Real nominal,
                   const Schedule& fixedSchedule, Rate fixedRate,
                   const DayCounter& fixedDayCount,
                   const Schedule& floatSchedule,
                   const boost::shared_ptr<IborIndex>& index,
                   Spread spread,
                   const DayCounter& floatingDayCount);
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
        MakeVanillaSwap& withType(VanillaSwap::Type type);
        MakeVanillaSwap& withNominal(Real n);

        MakeVanillaSwap& withSettlementDays(Natural settlementDays);
        MakeVanillaSwap& withEffectiveDate(const Date&);
        MakeVanillaSwap& withTerminationDate(const Date&);
        MakeVanillaSwap& withRule(DateGeneration::Rule r);

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
                              const boost::shared_ptr<PricingEngine>& engine);

        MakeVanillaSwap(const Period& swapTenor,
                        const boost::shared_ptr<IborIndex>& index,
                        Rate fixedRate,
                        const Period& forwardStart);
        
        %extend{
            boost::shared_ptr<VanillaSwap> makeVanillaSwap(){
                return (boost::shared_ptr<VanillaSwap>)(* $self);
            }
        }
};

#if defined(SWIGPYTHON)
%pythoncode{
def MakeVanillaSwap(swapTenor,iborIndex,fixedRate,forwardStart,
    receiveFixed=None, swapType=None, Nominal=None, settlementDays=None,
    effectiveDate=None, terminationDate=None, dateGenerationRule=None,
    fixedLegTenor=None, fixedLegCalendar=None, fixedLegConvention=None,
    fixedLegDayCount=None, floatingLegTenor=None, floatingLegCalendar=None,
    floatingLegConvention=None, floatingLegDayCount=None, floatingLegSpread=None,
    discountingTermStructure=None, pricingEngine=None,
    fixedLegTerminationDateConvention=None,  fixedLegDateGenRule=None,
    fixedLegEndOfMonth=None, fixedLegFirstDate=None, fixedLegNextToLastDate=None,
    floatingLegTerminationDateConvention=None,  floatingLegDateGenRule=None,
    floatingLegEndOfMonth=None, floatingLegFirstDate=None, floatingLegNextToLastDate=None):
    mv = _MakeVanillaSwap(swapTenor, iborIndex, fixedRate, forwardStart)
    if receiveFixed is not None:
        mv.receiveFixed(receiveFixed)
    if swapType is not None:
        mv.withType(swapType)
    if Nominal is not None:
        mv.withNominal(Nominal)
    if settlementDays is not None:
        mv.withSettlementDays(settlementDays)
    if effectiveDate is not None:
        mv.withEffectiveDate(effectiveDate)
    if terminationDate is not None:
        mv.withTerminationDate(terminationDate)
    if dateGenerationRule is not None:
        mv.withRule(dateGenerationRule)
    if fixedLegTenor is not None:
        mv.withFixedLegTenor(fixedLegTenor)
    if fixedLegCalendar is not None:
        mv.withFixedLegCalendar(fixedLegCalendar)
    if fixedLegConvention is not None:
        mv.withFixedLegConvention(fixedLegConvention)
    if fixedLegDayCount is not None:
        mv.withFixedLegDayCount(fixedLegDayCount)
    if floatingLegTenor is not None:
        mv.withFloatingLegTenor(floatingLegTenor)
    if floatingLegCalendar is not None:
        mv.withFloatingLegCalendar(floatingLegCalendar)
    if floatingLegConvention is not None:
        mv.withFloatingLegConvention(floatingLegConvention)
    if floatingLegDayCount is not None:
        mv.withFloatingLegDayCount(floatingLegDayCount)
    if floatingLegSpread is not None:
        mv.withFloatingLegSpread(floatingLegSpread)
    if discountingTermStructure is not None:
        mv.withDiscountingTermStructure(discountingTermStructure)
    if pricingEngine is not None:
        mv.withPricingEngine(pricingEngine)
    if fixedLegTerminationDateConvention is not None:
        mv.withFixedLegTerminationDateConvention(fixedLegTerminationDateConvention)
    if fixedLegDateGenRule is not None:
        mv.withFixedLegRule(fixedLegDateGenRule)
    if fixedLegEndOfMonth is not None:
        mv.withFixedLegEndOfMonth(fixedLegEndOfMonth)
    if fixedLegFirstDate is not None:
        mv.withFixedLegFirstDate(fixedLegFirstDate)
    if fixedLegNextToLastDate is not None:
        mv.withFixedLegNextToLastDate(fixedLegNextToLastDate)
    if floatingLegTerminationDateConvention is not None:
        mv.withFloatingLegTerminationDateConvention(floatingLegTerminationDateConvention)
    if floatingLegDateGenRule is not None:
        mv.withFloatingLegRule(floatingLegDateGenRule)
    if floatingLegEndOfMonth is not None:
        mv.withFloatingLegEndOfMonth(floatingLegEndOfMonth)
    if floatingLegFirstDate is not None:
        mv.withFloatingLegFirstDate(floatingLegFirstDate)
    if floatingLegNextToLastDate is not None:
        mv.withFloatingLegNextToLastDate(floatingLegNextToLastDate)
    return mv.makeVanillaSwap()
}
#endif

%shared_ptr(NonstandardSwap)
class NonstandardSwap : public Swap {
  public:
    NonstandardSwap(VanillaSwap::Type type,
                       const std::vector<Real> &fixedNominal,
                       const std::vector<Real> &floatingNominal,
                       const Schedule &fixedSchedule,
                       const std::vector<Real> &fixedRate,
                       const DayCounter &fixedDayCount,
                       const Schedule &floatSchedule,
                       const boost::shared_ptr<IborIndex> &index,
                       const std::vector<Real> &gearing,
                       const std::vector<Spread> &spread,
                       const DayCounter &floatDayCount,
                       const bool intermediateCapitalExchange = false,
                       const bool finalCapitalExchange = false,
                       BusinessDayConvention paymentConvention = Following);
    // Inspectors
    VanillaSwap::Type type() const;
    const std::vector<Real> &fixedNominal() const;
    const std::vector<Real> &floatingNominal() const;

    const Schedule &fixedSchedule() const;
    const std::vector<Real> &fixedRate() const;
    const DayCounter &fixedDayCount() const;

    const Schedule &floatingSchedule() const;
    const boost::shared_ptr<IborIndex> &iborIndex() const;
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
                                             boost::none,
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
                 const boost::shared_ptr<Bond>& bond,
                 Real bondCleanPrice,
                 const boost::shared_ptr<IborIndex>& index,
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
    FloatFloatSwap(VanillaSwap::Type type, const std::vector<Real> &nominal1,
        const std::vector<Real> &nominal2, const Schedule &schedule1,
        const boost::shared_ptr<InterestRateIndex> &index1,
        const DayCounter &dayCount1, const Schedule &schedule2,
        const boost::shared_ptr<InterestRateIndex> &index2,
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
    enum Type { Receiver = -1, Payer = 1 };
    
    OvernightIndexedSwap(
            OvernightIndexedSwap::Type type,
            Real nominal,
            const Schedule& schedule,
            Rate fixedRate,
            const DayCounter& fixedDC,
            const boost::shared_ptr<OvernightIndex>& index,
            Spread spread = 0.0,
            Natural paymentLag = 0,
            BusinessDayConvention paymentAdjustment = Following,
            Calendar paymentCalendar = Calendar(),
            bool telescopicValueDates = false);
    
    OvernightIndexedSwap(
            OvernightIndexedSwap::Type type,
            std::vector<Real> nominals,
            const Schedule& schedule,
            Rate fixedRate,
            const DayCounter& fixedDC,
            const boost::shared_ptr<OvernightIndex>& index,
            Spread spread = 0.0,
            Natural paymentLag = 0,
            BusinessDayConvention paymentAdjustment = Following,
            Calendar paymentCalendar = Calendar(),
            bool telescopicValueDates = false);


    Rate fixedLegBPS();
    Real fixedLegNPV();
    Real fairRate();
    Real overnightLegBPS();
    Real overnightLegNPV();
    Spread fairSpread();
    // Inspectors
    OvernightIndexedSwap::Type type();
    Real nominal();
    std::vector<Real> nominals();
    Frequency paymentFrequency();
    Rate fixedRate();
    const DayCounter& fixedDayCount();
    Spread spread();
    const Leg& fixedLeg();
    const Leg& overnightLeg();
};


%shared_ptr(OvernightIndexedSwapIndex)
class OvernightIndexedSwapIndex : public SwapIndex {
  public:
    OvernightIndexedSwapIndex(
              const std::string& familyName,
              const Period& tenor,
              Natural settlementDays,
              Currency currency,
              const boost::shared_ptr<OvernightIndex>& overnightIndex,
              bool telescopicValueDates = false);
    //! \name Inspectors
    //@{
    boost::shared_ptr<OvernightIndex> overnightIndex() const;
    /*! \warning Relinking the term structure underlying the index will
                 not have effect on the returned swap.
    */
    boost::shared_ptr<OvernightIndexedSwap> underlyingSwap(
                                            const Date& fixingDate) const;
};

%inline %{
    boost::shared_ptr<OvernightIndexedSwap> as_overnight_swap_index(
                          const boost::shared_ptr<InterestRateIndex>& index) {
        return boost::dynamic_pointer_cast<OvernightIndexedSwap>(index);
    }
%}

#endif
