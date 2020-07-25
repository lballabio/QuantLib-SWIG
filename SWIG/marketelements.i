
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier

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

#ifndef quantlib_market_elements_i
#define quantlib_market_elements_i

%include common.i
%include observer.i
%include functions.i

%{
using QuantLib::Quote;
%}

%shared_ptr(Quote)

class Quote : public Observable {
  private:
    Quote();
  public:
    Real value() const;
    bool isValid() const;
};

%template(QuoteHandle) Handle<Quote>;
%template(RelinkableQuoteHandle) RelinkableHandle<Quote>;

// actual quotes
%{
using QuantLib::SimpleQuote;
%}

%shared_ptr(SimpleQuote)

class SimpleQuote : public Quote {
  public:
    SimpleQuote(Real value);
    void setValue(Real value);
};


#if defined(SWIGPYTHON)
%{
typedef QuantLib::DerivedQuote<UnaryFunction> DerivedQuote;
typedef QuantLib::CompositeQuote<BinaryFunction> CompositeQuote;
%}

%shared_ptr(DerivedQuote)

class DerivedQuote : public Quote {
  public:
    %extend {
        DerivedQuote(const Handle<Quote>& h,
                     PyObject* function) {
            return new DerivedQuote(h,UnaryFunction(function));
        }
    }
};

%shared_ptr(CompositeQuote)

class CompositeQuote : public Quote {
  public:
    %extend {
        CompositeQuote(const Handle<Quote>& h1,
                       const Handle<Quote>& h2,
                       PyObject* function) {
            return new CompositeQuote(h1,h2,BinaryFunction(function));
        }
    }
};

#endif

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( boost::shared_ptr<Quote> )
SWIG_STD_VECTOR_ENHANCED( Handle<Quote> )
SWIG_STD_VECTOR_ENHANCED( RelinkableHandle<Quote> )
#endif
namespace std {
    %template(QuoteVector) vector<boost::shared_ptr<Quote> >;
    %template(QuoteVectorVector) vector<vector<boost::shared_ptr<Quote> > >;
    %template(QuoteHandleVector) vector<Handle<Quote> >;
    %template(QuoteHandleVectorVector) vector<vector<Handle<Quote> > >;
    %template(RelinkableQuoteHandleVector) vector<RelinkableHandle<Quote> >;
    %template(RelinkableQuoteHandleVectorVector)
                                  vector<vector<RelinkableHandle<Quote> > >;
}


#endif
