
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
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_market_elements_i
#define quantlib_market_elements_i

%include common.i
%include observer.i
%include functions.i
%include indexes.i
%include null.i

%{
using QuantLib::Quote;
using QuantLib::makeQuoteHandle;
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

RelinkableHandle<Quote> makeQuoteHandle(Real value);

#if defined(SWIGPYTHON)

%typecheck(SWIG_TYPECHECK_DOUBLE) std::variant<Real, Handle<Quote>>, const std::variant<Real, Handle<Quote>>& %{
    {
        int res = SWIG_AsVal_double($input, NULL);
        $1 = SWIG_CheckState(res);
    }
    if (!$1) {
        int res = SWIG_ConvertPtr($input, 0, $descriptor(Handle<Quote>*), SWIG_POINTER_NO_NULL | 0);
        $1 = SWIG_CheckState(res);;
    }
%}

%typemap(in) std::variant<Real, Handle<Quote>> (int res, Real val, void * argp) %{
    res = SWIG_AsVal_double($input, &val);
    if (SWIG_IsOK(res)) {
        $1 = val;
    } else {
        res = SWIG_ConvertPtr($input, &argp, $descriptor(Handle<Quote>*), SWIG_POINTER_NO_NULL);
        if (!SWIG_IsOK(res)) {
            SWIG_exception_fail(SWIG_ArgError(res), "in method '$symname', argument $argnum of type '$type'");
        } else if (!argp) {
            SWIG_exception_fail(SWIG_ValueError, "invalid null reference in method '$symname', argument $argnum of type '$type'");
        } else {
            Handle< Quote > * h = reinterpret_cast< Handle< Quote > * >(argp);
            $1 = *h;
            if (SWIG_IsNewObj(res)) delete h;
        }
    }
%}

%typemap(in) const std::variant<Real, Handle<Quote>>& (std::variant<Real, Handle<Quote>> temp, int res, Real val, void * argp) %{
    res = SWIG_AsVal_double($input, &val);
    if (SWIG_IsOK(res)) {
        temp = val;
        $1 = &temp;
    } else {
        res = SWIG_ConvertPtr($input, &argp, $descriptor(Handle<Quote>*), SWIG_POINTER_NO_NULL);
        if (!SWIG_IsOK(res)) {
            SWIG_exception_fail(SWIG_ArgError(res), "in method '$symname', argument $argnum of type '$type'");
        } else if (!argp) {
            SWIG_exception_fail(SWIG_ValueError, "invalid null reference in method '$symname', argument $argnum of type '$type'");
        } else {
            Handle< Quote > * h = reinterpret_cast< Handle< Quote > * >(argp);
            temp = *h;
            $1 = &temp;
            if (SWIG_IsNewObj(res)) delete h;
        }
    }
%}

#endif




// actual quotes
%{
using QuantLib::SimpleQuote;
using QuantLib::LastFixingQuote;
using QuantLib::FuturesConvAdjustmentQuote;
%}

%shared_ptr(SimpleQuote)

class SimpleQuote : public Quote {
  public:
    SimpleQuote(doubleOrNull value = Null<Real>());
    void setValue(doubleOrNull value = Null<Real>());
    void reset();
};

%shared_ptr(LastFixingQuote)

class LastFixingQuote : public Quote {
  public:
    LastFixingQuote(const ext::shared_ptr<Index>& index);
    ext::shared_ptr<Index> index() const;
    Date referenceDate() const;
};

%shared_ptr(FuturesConvAdjustmentQuote)

class FuturesConvAdjustmentQuote : public Quote {
  public:
    FuturesConvAdjustmentQuote(const ext::shared_ptr<IborIndex>& index,
                               const Date& futuresDate,
                               const Handle<Quote>& futuresQuote,
                               const Handle<Quote>& volatility,
                               const Handle<Quote>& meanReversion);
    FuturesConvAdjustmentQuote(const ext::shared_ptr<IborIndex>& index,
                               const std::string& immCode,
                               const Handle<Quote>& futuresQuote,
                               const Handle<Quote>& volatility,
                               const Handle<Quote>& meanReversion);
    Real futuresValue() const;
    Real volatility() const;
    Real meanReversion() const;
    Date immDate() const;
};

#if defined(SWIGPYTHON)
%{
using QuantLib::DerivedQuote;
using QuantLib::CompositeQuote;
%}

%shared_ptr(DerivedQuote<UnaryFunction>)

template <class F>
class DerivedQuote : public Quote {
  public:
    %extend {
        DerivedQuote(const Handle<Quote>& h,
                     PyObject* function) {
            return new DerivedQuote<F>(h,F(function));
        }
    }
};

%template(DerivedQuote) DerivedQuote<UnaryFunction>;

%shared_ptr(CompositeQuote<BinaryFunction>)

template <class F>
class CompositeQuote : public Quote {
  public:
    %extend {
        CompositeQuote(const Handle<Quote>& h1,
                       const Handle<Quote>& h2,
                       PyObject* function) {
            return new CompositeQuote<F>(h1,h2,F(function));
        }
    }
};

%template(CompositeQuote) CompositeQuote<BinaryFunction>;

#endif

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<Quote> )
SWIG_STD_VECTOR_ENHANCED( Handle<Quote> )
SWIG_STD_VECTOR_ENHANCED( RelinkableHandle<Quote> )
#endif
namespace std {
    %template(QuoteVector) vector<ext::shared_ptr<Quote> >;
    %template(QuoteVectorVector) vector<vector<ext::shared_ptr<Quote> > >;
    %template(QuoteHandleVector) vector<Handle<Quote> >;
    %template(QuoteHandleVectorVector) vector<vector<Handle<Quote> > >;
    %template(RelinkableQuoteHandleVector) vector<RelinkableHandle<Quote> >;
    %template(RelinkableQuoteHandleVectorVector)
                                  vector<vector<RelinkableHandle<Quote> > >;
}


#endif
