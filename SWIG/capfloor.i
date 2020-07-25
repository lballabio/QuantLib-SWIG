/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2018 Matthias Lungwitz
 Copyright (C) 2019 Wojciech Ślusarski

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
    enum Type { Cap, Floor, Collar };
    const Leg& floatingLeg() const;

    const std::vector<Rate>& capRates();
    const std::vector<Rate>& floorRates();
    Date startDate() const;
    Date maturityDate() const;
    Type type() const;
    
    Rate atmRate(const YieldTermStructure& discountCurve) const;

    %extend {
      const Real vega() {
        return self->result<Real>("vega");
      }
      const std::vector<Real> optionletsPrice() {
        return self->result<std::vector<Real> >("optionletsPrice");
      }
      const std::vector<Real> optionletsVega() {
        return self->result<std::vector<Real> >("optionletsVega");
      }
      const std::vector<Real> optionletsDelta() {
        return self->result<std::vector<Real> >("optionletsDelta");
      }
      const std::vector<DiscountFactor> optionletsDiscountFactor() {
        return self->result<std::vector<DiscountFactor> >("optionletsDiscountFactor");
      }
      const std::vector<Rate> optionletsAtmForward(){
        return self->result<std::vector<Real> >("optionletsAtmForward");
      }
      const std::vector<Rate> optionletsStdDev(){
        return self->result<std::vector<Real> >("optionletsStdDev");
      }
    }
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
%}

%shared_ptr(BlackCapFloorEngine)
class BlackCapFloorEngine : public PricingEngine {
  public:
    BlackCapFloorEngine(const Handle<YieldTermStructure>& termStructure,
                        const Handle<Quote>& vol,
                        const DayCounter& dc = Actual365Fixed(),
                        Real displacement = 0.0);
    BlackCapFloorEngine(const Handle<YieldTermStructure>& termStructure,
                        const Handle<OptionletVolatilityStructure>& vol,
                        Real displacement = Null<Real>());
};

%{
using QuantLib::BachelierCapFloorEngine;
%}

%shared_ptr(BachelierCapFloorEngine)
class BachelierCapFloorEngine : public PricingEngine {
  public:
    BachelierCapFloorEngine(const Handle<YieldTermStructure>& termStructure,
                            const Handle<Quote>& vol);
    BachelierCapFloorEngine(const Handle<YieldTermStructure>& termStructure,
                            const Handle<OptionletVolatilityStructure>& vol);
};

#endif
