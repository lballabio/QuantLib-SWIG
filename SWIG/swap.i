/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2007 StatPro Italia srl
 Copyright (C) 2011 Lluis Pujol Bajador
 Copyright (C) 2015 Gouthaman Balaraman
 Copyright (C) 2016 Peter Caspers
 Copyright (C) 2017, 2018 Matthias Lungwitz

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

%{
using QuantLib::Swap;
using QuantLib::VanillaSwap;
using QuantLib::NonstandardSwap;
using QuantLib::DiscountingSwapEngine;
using QuantLib::FloatFloatSwap;
using QuantLib::OvernightIndexedSwap;

typedef boost::shared_ptr<Instrument> SwapPtr;
typedef boost::shared_ptr<Instrument> VanillaSwapPtr;
typedef boost::shared_ptr<Instrument> NonstandardSwapPtr;
typedef boost::shared_ptr<PricingEngine> DiscountingSwapEnginePtr;
typedef boost::shared_ptr<Instrument> FloatFloatSwapPtr;
typedef boost::shared_ptr<Instrument> OvernightIndexedSwapPtr;
%}

%rename(Swap) SwapPtr;
class SwapPtr : public boost::shared_ptr<Instrument> {
  public:
    %extend {
        SwapPtr(const std::vector<boost::shared_ptr<CashFlow> >& firstLeg,
                const std::vector<boost::shared_ptr<CashFlow> >& secondLeg) {
            return new SwapPtr(new Swap(firstLeg, secondLeg));
        }
        Date startDate() {
            return boost::dynamic_pointer_cast<Swap>(*self)->startDate();
        }
        Date maturityDate() {
            return boost::dynamic_pointer_cast<Swap>(*self)->maturityDate();
        }
        const Leg & leg(Size i){
            return boost::dynamic_pointer_cast<Swap>(*self)->leg(i);
        }
        Real legNPV(Size j) const {
            return boost::dynamic_pointer_cast<Swap>(*self)->legNPV(j);
        }
    }
};


#if defined(SWIGJAVA) || defined(SWIGCSHARP)
%rename(_VanillaSwap) VanillaSwap;
#else
%ignore VanillaSwap;
#endif
class VanillaSwap {
  public:
    enum Type { Receiver = -1, Payer = 1 };
#if defined(SWIGJAVA) || defined(SWIGCSHARP)
  private:
    VanillaSwap();
#endif
};

%rename(VanillaSwap) VanillaSwapPtr;
class VanillaSwapPtr : public SwapPtr {
  public:
    %extend {
        static const VanillaSwap::Type Receiver = VanillaSwap::Receiver;
        static const VanillaSwap::Type Payer = VanillaSwap::Payer;
        VanillaSwapPtr(VanillaSwap::Type type, Real nominal,
                       const Schedule& fixedSchedule, Rate fixedRate,
                       const DayCounter& fixedDayCount,
                       const Schedule& floatSchedule,
                       const IborIndexPtr& index,
                       Spread spread,
                       const DayCounter& floatingDayCount) {
            boost::shared_ptr<IborIndex> libor =
                boost::dynamic_pointer_cast<IborIndex>(index);
            return new VanillaSwapPtr(
                    new VanillaSwap(type, nominal,fixedSchedule,fixedRate,
                                    fixedDayCount,floatSchedule,libor,
                                    spread, floatingDayCount));
        }
        Rate fairRate() {
            return boost::dynamic_pointer_cast<VanillaSwap>(*self)->fairRate();
        }
        Spread fairSpread() {
            return boost::dynamic_pointer_cast<VanillaSwap>(*self)
                 ->fairSpread();
        }
        Real fixedLegBPS() {
            return boost::dynamic_pointer_cast<VanillaSwap>(*self)
                 ->fixedLegBPS();
        }
        Real floatingLegBPS() {
            return boost::dynamic_pointer_cast<VanillaSwap>(*self)
                 ->floatingLegBPS();
        }   
        Real fixedLegNPV() {
            return boost::dynamic_pointer_cast<VanillaSwap> (*self)
                ->fixedLegNPV();
        }
        Real floatingLegNPV() {
            return boost::dynamic_pointer_cast<VanillaSwap> (*self)
                ->floatingLegNPV();
        }
        // Inspectors 
        const Leg& fixedLeg() {
            return boost::dynamic_pointer_cast<VanillaSwap> (*self)
                ->fixedLeg();
        }
        const Leg& floatingLeg() {
            return boost::dynamic_pointer_cast<VanillaSwap> (*self)
                ->floatingLeg();
        }
        Real nominal() {
            return boost::dynamic_pointer_cast<VanillaSwap> (*self)
                ->nominal();
        }
        const Schedule& fixedSchedule() {
            return boost::dynamic_pointer_cast<VanillaSwap> (*self)
                ->fixedSchedule();
        }
        const Schedule& floatingSchedule() {
            return boost::dynamic_pointer_cast<VanillaSwap> (*self)
                ->floatingSchedule();
        }
        Rate fixedRate() {
            return boost::dynamic_pointer_cast<VanillaSwap> (*self)
                ->fixedRate();
        }
        Spread spread() {
            return boost::dynamic_pointer_cast<VanillaSwap> (*self)
                ->spread();
        }
        const DayCounter& floatingDayCount() {
            return boost::dynamic_pointer_cast<VanillaSwap> (*self)
                ->floatingDayCount();
        }
        const DayCounter& fixedDayCount() {
            return boost::dynamic_pointer_cast<VanillaSwap> (*self)
                ->fixedDayCount();
        }
    }
};

#if defined(SWIGJAVA) || defined(SWIGCSHARP)
%rename(_NonstandardSwap) NonstandardSwap;
#else
%ignore NonstandardSwap;
#endif
class NonstandardSwap {
  public:
#if defined(SWIGJAVA) || defined(SWIGCSHARP)
  private:
    NonstandardSwap();
#endif
};

%rename(NonstandardSwap) NonstandardSwapPtr;
class NonstandardSwapPtr : public SwapPtr {
  public:
    %extend {
        NonstandardSwapPtr(VanillaSwap::Type type,
                           const std::vector<Real> &fixedNominal,
                           const std::vector<Real> &floatingNominal,
                           const Schedule &fixedSchedule,
                           const std::vector<Real> &fixedRate,
                           const DayCounter &fixedDayCount,
                           const Schedule &floatSchedule,
                           const IborIndexPtr &index,
                           const std::vector<Real> &gearing,
                           const std::vector<Spread> &spread,
                           const DayCounter &floatDayCount,
                           const bool intermediateCapitalExchange = false,
                           const bool finalCapitalExchange = false,
                           BusinessDayConvention paymentConvention = Following) {
            boost::shared_ptr<IborIndex> libor =
                boost::dynamic_pointer_cast<IborIndex>(index);
            return new NonstandardSwapPtr(
                new NonstandardSwap(type,fixedNominal,floatingNominal,fixedSchedule,fixedRate,
                                    fixedDayCount,floatSchedule,libor,gearing,
                                    spread,floatDayCount,intermediateCapitalExchange,
                                    finalCapitalExchange, paymentConvention));
        }
        // Inspectors
        const Leg& fixedLeg() {
            return boost::dynamic_pointer_cast<NonstandardSwap> (*self)
                ->fixedLeg();
        }
        const Leg& floatingLeg() {
            return boost::dynamic_pointer_cast<NonstandardSwap> (*self)
                ->floatingLeg();
        }
        std::vector<Real> fixedNominals() {
            return boost::dynamic_pointer_cast<NonstandardSwap> (*self)
                ->fixedNominal();
        }
        std::vector<Real> floatingNominals() {
            return boost::dynamic_pointer_cast<NonstandardSwap> (*self)
                ->floatingNominal();
        }
        const Schedule& fixedSchedule() {
            return boost::dynamic_pointer_cast<NonstandardSwap> (*self)
                ->fixedSchedule();
        }
        const Schedule& floatingSchedule() {
            return boost::dynamic_pointer_cast<NonstandardSwap> (*self)
                ->floatingSchedule();
        }
        std::vector<Rate> fixedRate() {
            return boost::dynamic_pointer_cast<NonstandardSwap> (*self)
                ->fixedRate();
        }
        std::vector<Spread> spreads() {
            return boost::dynamic_pointer_cast<NonstandardSwap> (*self)
                ->spreads();
        }
        std::vector<Spread> gearings() {
            return boost::dynamic_pointer_cast<NonstandardSwap> (*self)
                ->gearings();
        }
        const DayCounter& floatingDayCount() {
            return boost::dynamic_pointer_cast<NonstandardSwap> (*self)
                ->floatingDayCount();
        }
        const DayCounter& fixedDayCount() {
            return boost::dynamic_pointer_cast<NonstandardSwap> (*self)
                ->fixedDayCount();
        }
    }
};

%rename(DiscountingSwapEngine) DiscountingSwapEnginePtr;
class DiscountingSwapEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        DiscountingSwapEnginePtr(
                            const Handle<YieldTermStructure>& discountCurve,
                            const Date& settlementDate = Date(),
                            const Date& npvDate = Date()) {
            return new DiscountingSwapEnginePtr(
                                    new DiscountingSwapEngine(discountCurve,
                                                              boost::none,
                                                              settlementDate,
                                                              npvDate));
        }
        DiscountingSwapEnginePtr(
                            const Handle<YieldTermStructure>& discountCurve,
                            bool includeSettlementDateFlows,
                            const Date& settlementDate = Date(),
                            const Date& npvDate = Date()) {
            return new DiscountingSwapEnginePtr(
                         new DiscountingSwapEngine(discountCurve,
                                                   includeSettlementDateFlows,
                                                   settlementDate,
                                                   npvDate));
        }
    }
};


%{
using QuantLib::AssetSwap;
typedef boost::shared_ptr<Instrument> AssetSwapPtr;
%}

%rename(AssetSwap) AssetSwapPtr;
class AssetSwapPtr : public SwapPtr {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") AssetSwapPtr;
    #endif
  public:
    %extend {
        AssetSwapPtr(bool payFixedRate,
                     const BondPtr& bond,
                     Real bondCleanPrice,
                     const InterestRateIndexPtr& index,
                     Spread spread,
                     const Schedule& floatSchedule = Schedule(),
                     const DayCounter& floatingDayCount = DayCounter(),
                     bool parAssetSwap = true) {
            const boost::shared_ptr<Bond> b =
                boost::dynamic_pointer_cast<Bond>(bond);
            const boost::shared_ptr<IborIndex> i =
                boost::dynamic_pointer_cast<IborIndex>(index);
            return new AssetSwapPtr(
                new AssetSwap(payFixedRate,b,bondCleanPrice,i,spread,
                              floatSchedule,floatingDayCount,parAssetSwap));
        }
        Real fairCleanPrice() {
            return boost::dynamic_pointer_cast<AssetSwap>(*self)
                ->fairCleanPrice();
        }
        Spread fairSpread() {
            return boost::dynamic_pointer_cast<AssetSwap>(*self)
                ->fairSpread();
        }
    }
};

%rename(FloatFloatSwap) FloatFloatSwapPtr;
class FloatFloatSwapPtr : public SwapPtr {

  public:
    %extend {
        FloatFloatSwapPtr(VanillaSwap::Type type, const std::vector<Real> &nominal1,
            const std::vector<Real> &nominal2, const Schedule &schedule1,
            const InterestRateIndexPtr &indexPtr1,
            const DayCounter &dayCount1, const Schedule &schedule2,
            const InterestRateIndexPtr &indexPtr2,
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
            BusinessDayConvention paymentConvention2 = Following) {
            boost::shared_ptr<InterestRateIndex> index1 =
                boost::dynamic_pointer_cast<InterestRateIndex>(indexPtr1);
            boost::shared_ptr<InterestRateIndex> index2 =
                boost::dynamic_pointer_cast<InterestRateIndex>(indexPtr2);
            return new FloatFloatSwapPtr(
                    new FloatFloatSwap(type, nominal1,nominal2,schedule1,
                                    index1,dayCount1,schedule2,
                                    index2, dayCount2,
                                    intermediateCapitalExchange, finalCapitalExchange,
                                    gearing1, spread1, cappedRate1,
                                    flooredRate1, gearing2, spread2,
                                    cappedRate2, flooredRate2,
                                    paymentConvention1, paymentConvention2));
        }

    }
};

#if defined(SWIGJAVA) || defined(SWIGCSHARP)
%rename(_OvernightIndexedSwap) OvernightIndexedSwap;
#else
%ignore OvernightIndexedSwap;
#endif
class OvernightIndexedSwap {
  public:
    enum Type { Receiver = -1, Payer = 1 };
#if defined(SWIGJAVA) || defined(SWIGCSHARP)
  private:
    OvernightIndexedSwap();
#endif
};

%rename(OvernightIndexedSwap) OvernightIndexedSwapPtr;
class OvernightIndexedSwapPtr : public SwapPtr {
  public:
    %extend {
        static const OvernightIndexedSwap::Type Receiver = OvernightIndexedSwap::Receiver;
        static const OvernightIndexedSwap::Type Payer = OvernightIndexedSwap::Payer;
		
		OvernightIndexedSwapPtr(
				OvernightIndexedSwap::Type type,
				Real nominal,
				const Schedule& schedule,
				Rate fixedRate,
				const DayCounter& fixedDC,
				const OvernightIndexPtr& overnightIndex,
				Spread spread = 0.0,
				Natural paymentLag = 0,
				BusinessDayConvention paymentAdjustment = Following,
				Calendar paymentCalendar = Calendar(),
				bool telescopicValueDates = false) {
				boost::shared_ptr<OvernightIndex> index =
					boost::dynamic_pointer_cast<OvernightIndex>(overnightIndex);
				return new OvernightIndexedSwapPtr(
				 new OvernightIndexedSwap(type, nominal, schedule, fixedRate, fixedDC,
				index, spread, paymentLag, paymentAdjustment, paymentCalendar, telescopicValueDates));
		}
		
		OvernightIndexedSwapPtr(
				OvernightIndexedSwap::Type type,
				std::vector<Real> nominals,
				const Schedule& schedule,
				Rate fixedRate,
				const DayCounter& fixedDC,
				const OvernightIndexPtr& overnightIndex,
				Spread spread = 0.0,
				Natural paymentLag = 0,
				BusinessDayConvention paymentAdjustment = Following,
				Calendar paymentCalendar = Calendar(),
				bool telescopicValueDates = false) {
				boost::shared_ptr<OvernightIndex> index =
					boost::dynamic_pointer_cast<OvernightIndex>(overnightIndex);
				return new OvernightIndexedSwapPtr(
				 new OvernightIndexedSwap(type, nominals, schedule, fixedRate, fixedDC,
				index, spread, paymentLag, paymentAdjustment, paymentCalendar, telescopicValueDates));
		}


        Rate fixedLegBPS() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap>(*self)->fixedLegBPS();
        }
        Real fixedLegNPV() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap>(*self)
                 ->fixedLegNPV();
        }
        Real fairRate() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap>(*self)
                 ->fairRate();
        }
        Real overnightLegBPS() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap>(*self)
                 ->overnightLegBPS();
        }
        Real overnightLegNPV() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap> (*self)
                ->overnightLegNPV();
        }
        Spread fairSpread() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap> (*self)
                ->fairSpread();
        }
        // Inspectors
		OvernightIndexedSwap::Type type() {
			return boost::dynamic_pointer_cast<OvernightIndexedSwap>(*self)->type();
		}
        Real nominal() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap> (*self)
                ->nominal();
        }
		std::vector<Real> nominals() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap> (*self)
                ->nominals();
        }
		Frequency paymentFrequency() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap> (*self)
                ->paymentFrequency();
        }
		Rate fixedRate() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap> (*self)
                ->fixedRate();
        }
		const DayCounter& fixedDayCount() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap> (*self)
                ->fixedDayCount();
        }
		Spread spread() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap> (*self)
                ->spread();
        }
        const Leg& fixedLeg() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap> (*self)
                ->fixedLeg();
        }
        const Leg& overnightLeg() {
            return boost::dynamic_pointer_cast<OvernightIndexedSwap> (*self)
                ->overnightLeg();
        }
    }
};

#endif
