
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2016 Peter Caspers
 Copyright (C) 2017, 2018, 2019 Matthias Lungwitz
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

#ifndef quantlib_swaption_i
#define quantlib_swaption_i

%include options.i
%include marketelements.i
%include termstructures.i
%include volatilities.i
%include swap.i
%include old_volatility.i
%include daycounters.i

%{
using QuantLib::Actual365Fixed;
using QuantLib::Swaption;
using QuantLib::NonstandardSwaption;
using QuantLib::Settlement;
using QuantLib::FloatFloatSwaption;
%}

struct Settlement {
    enum Type { Physical, Cash };
    enum Method { PhysicalOTC, PhysicalCleared, CollateralizedCashPrice, ParYieldCurve };
};
    
%shared_ptr(Swaption)
class Swaption : public Option {
  public:
    Swaption(const ext::shared_ptr<FixedVsFloatingSwap>& swap,
             const ext::shared_ptr<Exercise>& exercise,
             Settlement::Type type = Settlement::Physical,
             Settlement::Method settlementMethod = Settlement::PhysicalOTC);
    
    Settlement::Type settlementType() const;       
    Settlement::Method settlementMethod() const;
    VanillaSwap::Type type() const;
    const ext::shared_ptr<FixedVsFloatingSwap>& underlying() const;
    const ext::shared_ptr<VanillaSwap>& underlyingSwap() const;
    
    //! implied volatility
    Volatility impliedVolatility(
                          Real price,
                          const Handle<YieldTermStructure>& discountCurve,
                          Volatility guess,
                          Real accuracy = 1.0e-4,
                          Natural maxEvaluations = 100,
                          Volatility minVol = 1.0e-7,
                          Volatility maxVol = 4.0,
                          VolatilityType type = ShiftedLognormal,
                          Real displacement = 0.0) const;
    %extend {
        Real vega() {
            return self->result<Real>("vega");
        }

        Real delta() {
            return self->result<Real>("delta");
        }

        Real annuity() {
            return self->result<Real>("annuity");
        }
    }
};

%{
using QuantLib::BasketGeneratingEngine;
%}

%shared_ptr(NonstandardSwaption)
class NonstandardSwaption : public Instrument {
  public:
    NonstandardSwaption(const ext::shared_ptr<NonstandardSwap>& swap,
                const ext::shared_ptr<Exercise>& exercise,
                Settlement::Type type = Settlement::Physical,
                Settlement::Method settlementMethod = Settlement::PhysicalOTC);
                
    const ext::shared_ptr<NonstandardSwap> &underlyingSwap() const;

    %extend {                
        std::vector<ext::shared_ptr<BlackCalibrationHelper> > calibrationBasket(
            ext::shared_ptr<SwapIndex> swapIndex,
            ext::shared_ptr<SwaptionVolatilityStructure> swaptionVolatility,
            std::string typeStr) {

            BasketGeneratingEngine::CalibrationBasketType type;
            if(typeStr == "Naive")
                type = BasketGeneratingEngine::Naive;
            else if(typeStr == "MaturityStrikeByDeltaGamma")
                type = BasketGeneratingEngine::MaturityStrikeByDeltaGamma;
            else
                QL_FAIL("type " << typeStr << "unknown.");

            std::vector<ext::shared_ptr<BlackCalibrationHelper> > hs =
                self->calibrationBasket(swapIndex, swaptionVolatility, type);
            std::vector<ext::shared_ptr<BlackCalibrationHelper> > helpers(hs.size());
            for (Size i=0; i<hs.size(); ++i)
                helpers[i] = hs[i];
            return helpers;
        }


        std::vector<Real> probabilities() {
            return self->result<std::vector<Real> >("probabilities");
        }
    }
};

%shared_ptr(FloatFloatSwaption)
class FloatFloatSwaption : public Instrument {
public:
    FloatFloatSwaption(const ext::shared_ptr<FloatFloatSwap>& swap,
                const ext::shared_ptr<Exercise>& exercise,
                Settlement::Type delivery = Settlement::Physical,
                Settlement::Method settlementMethod = Settlement::PhysicalOTC);

    const ext::shared_ptr<FloatFloatSwap> &underlyingSwap();
    
    %extend {

        std::vector<ext::shared_ptr<BlackCalibrationHelper> > calibrationBasket(
        ext::shared_ptr<SwapIndex> swapIndex,
        ext::shared_ptr<SwaptionVolatilityStructure> swaptionVolatility,
        std::string typeStr) {

        BasketGeneratingEngine::CalibrationBasketType type;
        if(typeStr == "Naive")
            type = BasketGeneratingEngine::Naive;
        else if(typeStr == "MaturityStrikeByDeltaGamma")
            type = BasketGeneratingEngine::MaturityStrikeByDeltaGamma;
        else
            QL_FAIL("type " << typeStr << "unknown.");

        std::vector<ext::shared_ptr<BlackCalibrationHelper> > hs =
            self->calibrationBasket(swapIndex, swaptionVolatility, type);
        std::vector<ext::shared_ptr<BlackCalibrationHelper> > helpers(hs.size());
        for (Size i=0; i<hs.size(); ++i)
            helpers[i] = hs[i];
        return helpers;
        }

        Real underlyingValue() {
            return self->result<Real>("underlyingValue");
        }

        std::vector<Real> probabilities() {
            return self->result<std::vector<Real> >("probabilities");
        }
    }
};

// pricing engines

%{
using QuantLib::BlackSwaptionEngine;
using QuantLib::BachelierSwaptionEngine;
%}

%shared_ptr(BlackSwaptionEngine)
class BlackSwaptionEngine : public PricingEngine {
  public:
    enum CashAnnuityModel { SwapRate, DiscountCurve };
    BlackSwaptionEngine(const Handle<YieldTermStructure> & discountCurve,
                        const Handle<Quote>& vol,
                        const DayCounter& dc = Actual365Fixed(),
                        Real displacement = 0.0,
                        CashAnnuityModel model = DiscountCurve);
    BlackSwaptionEngine(const Handle<YieldTermStructure> & discountCurve,
                        const Handle<SwaptionVolatilityStructure>& v,
                        CashAnnuityModel model = DiscountCurve);
};

%shared_ptr(BachelierSwaptionEngine)
class BachelierSwaptionEngine : public PricingEngine {
  public:
    enum CashAnnuityModel { SwapRate, DiscountCurve };
    BachelierSwaptionEngine(const Handle<YieldTermStructure> & discountCurve,
                            const Handle<Quote>& vol,
                            const DayCounter& dc = Actual365Fixed(),
                            CashAnnuityModel model = DiscountCurve);
    BachelierSwaptionEngine(const Handle<YieldTermStructure> & discountCurve,
                            const Handle<SwaptionVolatilityStructure>& v,
                            CashAnnuityModel model = DiscountCurve);
};

#endif

