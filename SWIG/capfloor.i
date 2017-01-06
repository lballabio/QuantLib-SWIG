/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl

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

#ifndef quantlib_cap_floor_i
#define quantlib_cap_floor_i

%include options.i
%include marketelements.i
%include termstructures.i
%include cashflows.i
%include volatilities.i

%{
using QuantLib::CapFloor;
using QuantLib::Cap;
using QuantLib::Floor;
using QuantLib::Collar;

typedef boost::shared_ptr<Instrument> CapFloorPtr;
typedef boost::shared_ptr<Instrument> CapPtr;
typedef boost::shared_ptr<Instrument> FloorPtr;
typedef boost::shared_ptr<Instrument> CollarPtr;
%}

%rename(CapFloor) CapFloorPtr;
class CapFloorPtr : public boost::shared_ptr<Instrument> {
    #if defined(SWIGMZSCHEME) || defined(SWIGGUILE)
    %rename("implied-volatility") impliedVolatility;
    %rename("floating-leg") floatingLeg;
    %rename("cap-rates") capRates;
    %rename("floor-rates") floorRates;
    %rename("start-date") startDate;
    %rename("maturity-date") maturityDate;
    %rename("atm-rate") atmRate;
    #endif
  public:
     %extend {
        Volatility impliedVolatility(Real price,
                                     const Handle<YieldTermStructure>& disc,
                                     Volatility guess,
                                     Real accuracy = 1.0e-4,
                                     Natural maxEvaluations = 100,
                                     Volatility minVol = 1.0e-7,
                                     Volatility maxVol = 4.0,
                                     VolatilityType type = ShiftedLognormal,
                                     Real displacement = 0.0) const {
            return boost::dynamic_pointer_cast<CapFloor>(*self)->
                impliedVolatility(price, disc, guess, accuracy,
                                  maxEvaluations, minVol, maxVol,
                                  type, displacement);
        }
        const Leg& floatingLeg() const {
            return boost::dynamic_pointer_cast<CapFloor>(*self)->floatingLeg();
        }

        const std::vector<Rate>& capRates() const {
            return boost::dynamic_pointer_cast<CapFloor>(*self)->capRates();
        }
        const std::vector<Rate>& floorRates() const {
            return boost::dynamic_pointer_cast<CapFloor>(*self)->floorRates();
        }
        Date startDate() const {
            return boost::dynamic_pointer_cast<CapFloor>(*self)->startDate();
        }
        Date maturityDate() const {
            return boost::dynamic_pointer_cast<CapFloor>(*self)->maturityDate();
        }

        Rate atmRate(const boost::shared_ptr<YieldTermStructure>& discountCurve) const {
            return boost::dynamic_pointer_cast<CapFloor>(*self)->atmRate(*discountCurve);
        }
    }
};



%rename(Cap) CapPtr;
class CapPtr : public CapFloorPtr {
  public:
    %extend {
        CapPtr(const std::vector<boost::shared_ptr<CashFlow> >& leg,
               const std::vector<Rate>& capRates) {
            return new CapPtr(new Cap(leg,capRates));
        }
    }
};

%rename(Floor) FloorPtr;
class FloorPtr : public CapFloorPtr {
  public:
    %extend {
        FloorPtr(const std::vector<boost::shared_ptr<CashFlow> >& leg,
                 const std::vector<Rate>& floorRates) {
            return new FloorPtr(new Floor(leg,floorRates));
        }
    }
};

%rename(Collar) CollarPtr;
class CollarPtr : public CapFloorPtr {
  public:
    %extend {
        CollarPtr(const std::vector<boost::shared_ptr<CashFlow> >& leg,
                  const std::vector<Rate>& capRates,
                  const std::vector<Rate>& floorRates) {
            return new CollarPtr(new Collar(leg,capRates,floorRates));
        }
    }
};

%{
using QuantLib::BlackCapFloorEngine;
typedef boost::shared_ptr<PricingEngine> BlackCapFloorEnginePtr;
%}

%rename(BlackCapFloorEngine) BlackCapFloorEnginePtr;
class BlackCapFloorEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        BlackCapFloorEnginePtr(
                           const Handle<YieldTermStructure>& termStructure,
                           const Handle<Quote>& vol) {
            return new BlackCapFloorEnginePtr(
                                  new BlackCapFloorEngine(termStructure,vol));
        }
        BlackCapFloorEnginePtr(
                           const Handle<YieldTermStructure>& termStructure,
                           const Handle<OptionletVolatilityStructure>& vol) {
            return new BlackCapFloorEnginePtr(
                                  new BlackCapFloorEngine(termStructure,vol));
        }
    }
};

%{
using QuantLib::BachelierCapFloorEngine;
typedef boost::shared_ptr<PricingEngine> BachelierCapFloorEnginePtr;
%}

%rename(BachelierCapFloorEngine) BachelierCapFloorEnginePtr;
class BachelierCapFloorEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        BachelierCapFloorEnginePtr(
                           const Handle<YieldTermStructure>& termStructure,
                           const Handle<Quote>& vol) {
            return new BachelierCapFloorEnginePtr(
                                new BachelierCapFloorEngine(termStructure,vol));
        }
        BachelierCapFloorEnginePtr(
                           const Handle<YieldTermStructure>& termStructure,
                           const Handle<OptionletVolatilityStructure>& vol) {
            return new BachelierCapFloorEnginePtr(
                                new BachelierCapFloorEngine(termStructure,vol));
        }
    }
};

#endif
