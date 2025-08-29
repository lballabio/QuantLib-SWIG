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
 <https://www.quantlib.org/license.shtml>.

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
%include swap.i
%include volatilities.i

%{
using QuantLib::CapFloor;
using QuantLib::Cap;
using QuantLib::Floor;
using QuantLib::Collar;
using QuantLib::MakeCapFloor;
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
    Cap(const std::vector<ext::shared_ptr<CashFlow> >& leg,
           const std::vector<Rate>& capRates);
};

%shared_ptr(Floor)
class Floor : public CapFloor {
  public:
    Floor(const std::vector<ext::shared_ptr<CashFlow> >& leg,
             const std::vector<Rate>& floorRates);
};

%shared_ptr(Collar)
class Collar : public CapFloor {
  public:
    Collar(const std::vector<ext::shared_ptr<CashFlow> >& leg,
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

#if defined(SWIGPYTHON)
%rename (_MakeCapFloor) MakeCapFloor;
#endif
class MakeCapFloor {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MakeCapFloor;
    #endif
    public:
        MakeCapFloor& withNominal(Real n);

        MakeCapFloor& withEffectiveDate(const Date&, bool firstCapletExcluded);
        MakeCapFloor& withTenor(const Period&);
        MakeCapFloor& withCalendar(const Calendar&);
        MakeCapFloor& withConvention(BusinessDayConvention bdc);
        MakeCapFloor& withTerminationDateConvention(BusinessDayConvention bdc);
        MakeCapFloor& withRule(DateGeneration::Rule r);
        MakeCapFloor& withEndOfMonth(bool flag = true);
        MakeCapFloor& withFirstDate(const Date&);
        MakeCapFloor& withNextToLastDate(const Date&);
        MakeCapFloor& withDayCount(const DayCounter&);

        MakeCapFloor& asOptionlet(bool b = true);

        MakeCapFloor& withPricingEngine(
                              const ext::shared_ptr<PricingEngine>& engine);

        MakeCapFloor(CapFloor::Type capFloorType,
                     const Period& capFloorTenor,
                     const ext::shared_ptr<IborIndex>& iborIndex,
                     doubleOrNull strike = Null<Rate>(),
                     const Period& forwardStart = 0*Days);
        
        %extend {
            ext::shared_ptr<CapFloor> makeCapFloor() {
                return (ext::shared_ptr<CapFloor>)(* $self);
            }
        }
};

#if defined(SWIGPYTHON)
%pythoncode {
_MAKECAPFLOOR_METHODS = {
    "nominal": "withNominal",
    "effectiveDate": "withEffectiveDate",
    "tenor": "withTenor",
    "calendar": "withCalendar",
    "convention": "withConvention",
    "terminationDateConvention": "withTerminationDateConvention",
    "rule": "withRule",
    "endOfMonth": "withEndOfMonth",
    "firstDate": "withFirstDate",
    "nextToLastDate": "withNextToLastDate",
    "dayCount": "withDayCount",
    "asOptionlet": "asOptionlet",
    "pricingEngine": "withPricingEngine",
}

def MakeCapFloor(capFloorType, capFloorTenor, iborIndex, strike=None, forwardStart=Period(0, Days), **kwargs):
    mv = _MakeCapFloor(capFloorType, capFloorTenor, iborIndex, strike, forwardStart)
    _apply_kwargs("MakeCapFloor", _MAKECAPFLOOR_METHODS, mv, kwargs)
    return mv.makeCapFloor()
}
#endif

#endif
