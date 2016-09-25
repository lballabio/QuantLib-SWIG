
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2016 Peter Caspers

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
using QuantLib::Settlement;
typedef boost::shared_ptr<Instrument> SwaptionPtr;
%}

struct Settlement {
   enum Type { Physical, Cash };
};

%rename(Swaption) SwaptionPtr;
class SwaptionPtr : public boost::shared_ptr<Instrument> {
  public:
    %extend {
        SwaptionPtr(const VanillaSwapPtr& simpleSwap,
                    const boost::shared_ptr<Exercise>& exercise,
                    Settlement::Type type = Settlement::Physical) {
            boost::shared_ptr<VanillaSwap> swap =
                 boost::dynamic_pointer_cast<VanillaSwap>(simpleSwap);
            QL_REQUIRE(swap, "simple swap required");
            return new SwaptionPtr(new Swaption(swap,exercise,type));
        }
    }
};


// pricing engines

%{
using QuantLib::BlackSwaptionEngine;
using QuantLib::BachelierSwaptionEngine;
typedef boost::shared_ptr<PricingEngine> BlackSwaptionEnginePtr;
typedef boost::shared_ptr<PricingEngine> BachelierSwaptionEnginePtr;
%}

%rename(BlackSwaptionEngine) BlackSwaptionEnginePtr;
class BlackSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        BlackSwaptionEnginePtr(
                           const Handle<YieldTermStructure> & discountCurve,
                           const Handle<Quote>& vol,
                           const DayCounter& dc = Actual365Fixed(),
                           Real displacement = 0.0) {
            return new BlackSwaptionEnginePtr(
                          new BlackSwaptionEngine(discountCurve, vol, dc, displacement));
        }
        BlackSwaptionEnginePtr(
                           const Handle<YieldTermStructure> & discountCurve,
                           const Handle<SwaptionVolatilityStructure>& v) {
            return new BlackSwaptionEnginePtr(
                                   new BlackSwaptionEngine(discountCurve, v));
        }
    }
};

%rename(BachelierSwaptionEngine) BachelierSwaptionEnginePtr;
class BachelierSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        BachelierSwaptionEnginePtr(
                           const Handle<YieldTermStructure> & discountCurve,
                           const Handle<Quote>& vol,
                           const DayCounter& dc = Actual365Fixed()) {
            return new BlackSwaptionEnginePtr(
                          new BachelierSwaptionEngine(discountCurve, vol, dc));
        }
        BachelierSwaptionEnginePtr(
                           const Handle<YieldTermStructure> & discountCurve,
                           const Handle<SwaptionVolatilityStructure>& v) {
            return new BlackSwaptionEnginePtr(
                                   new BachelierSwaptionEngine(discountCurve, v));
        }
    }
};

#endif

