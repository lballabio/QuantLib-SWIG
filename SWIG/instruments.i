
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

#ifndef quantlib_instruments_i
#define quantlib_instruments_i

%include common.i
%include types.i
%include marketelements.i
%include observer.i
%include stl.i

// pricing engine

%{
using QuantLib::PricingEngine;
%}

%shared_ptr(PricingEngine)
class PricingEngine : public Observable {
  private:
    PricingEngine();
};

// instrument

%{
using QuantLib::Instrument;
%}
    
%shared_ptr(Instrument)
class Instrument : public Observable {
  public:
    Real NPV() const;
    Real errorEstimate() const;
    bool isExpired() const;
    void setPricingEngine(const boost::shared_ptr<PricingEngine>&);
    void recalculate();
    void freeze();
    void unfreeze();
  private:
    Instrument();
};

#if defined(SWIGR)
%Rruntime %{
setMethod("summary", "_p_boost__shared_ptrTInstrument_t",
function(object) c(value=object$NPV()))

setMethod("print", "_p_boost__shared_ptrTInstrument_t",
function(x) print(summary(x)))
%}
#endif

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( boost::shared_ptr<Instrument> )
#endif
namespace std {
    %template(InstrumentVector) vector<boost::shared_ptr<Instrument> >;
}

// actual instruments

%{
using QuantLib::Stock;
%}

%shared_ptr(Stock)
class Stock : public Instrument {
  public:
    Stock(const Handle<Quote>& quote);
};


%{
using QuantLib::CompositeInstrument;
%}

%shared_ptr(CompositeInstrument)
class CompositeInstrument : public Instrument {
  public:
    CompositeInstrument();
    void add(const boost::shared_ptr<Instrument>& instrument,
             Real multiplier = 1.0);
    void subtract(const boost::shared_ptr<Instrument>& instrument,
                  Real multiplier = 1.0);
};


#endif
