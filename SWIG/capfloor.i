/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
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
%}

%shared_ptr(CapFloor)
class CapFloor : public Instrument {
  public:
    Volatility impliedVolatility(Real price,
                                 const Handle<YieldTermStructure>& disc,
                                 Volatility guess,
                                 Real accuracy = 1.0e-4,
                                 Natural maxEvaluations = 100,
                                 Volatility minVol = 1.0e-7,
                                 Volatility maxVol = 4.0,
                                 VolatilityType type = ShiftedLognormal,
                                 Real displacement = 0.0) const;
    const Leg& floatingLeg() const;

    const std::vector<Rate>& capRates();
    const std::vector<Rate>& floorRates();
    Date startDate() const;
    Date maturityDate() const;

    Rate atmRate(const boost::shared_ptr<YieldTermStructure>& discountCurve);
};


%shared_ptr(Cap)
class Cap : public CapFloor {
  public:
    Cap(const std::vector<boost::shared_ptr<CashFlow> >& leg,
           const std::vector<Rate>& capRates);
};

%shared_ptr(Floor)
class Floor : public CapFloor {
  public:
    Floor(const std::vector<boost::shared_ptr<CashFlow> >& leg,
             const std::vector<Rate>& floorRates);
};

%shared_ptr(Collar)
class Collar : public CapFloor {
  public:
    Collar(const std::vector<boost::shared_ptr<CashFlow> >& leg,
              const std::vector<Rate>& capRates,
              const std::vector<Rate>& floorRates);
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
