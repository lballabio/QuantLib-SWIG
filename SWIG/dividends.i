
/*
 Copyright (C) 2006 StatPro Italia srl

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

#ifndef quantlib_dividends_i
#define quantlib_dividends_i

%include cashflows.i

%{
using QuantLib::Dividend;
%}

%shared_ptr(Dividend)
class Dividend : public CashFlow {
  private:
    Dividend();
};

%{
using QuantLib::FixedDividend;
using QuantLib::FractionalDividend;
%}

%shared_ptr(FixedDividend)
class FixedDividend : public Dividend {
  public:
    FixedDividend(Real amount, const Date& date);
};

%shared_ptr(FractionalDividend)
class FractionalDividend : public Dividend {
  public:
    FractionalDividend(Rate rate, const Date& date);
};


%{
using QuantLib::DividendSchedule;
%}

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<Dividend> )
#endif
%template(DividendSchedule) std::vector<ext::shared_ptr<Dividend>>;
typedef std::vector<ext::shared_ptr<Dividend>> DividendSchedule;

#endif
