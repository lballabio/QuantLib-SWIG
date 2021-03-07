/*
 Copyright (C) 2008 StatPro Italia srl

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

#ifndef quantlib_futures_i
#define quantlib_futures_i

%include forward.i
%include indexes.i
%include options.i

%{
using QuantLib::Futures;
using QuantLib::OvernightIndexFuture;
%}

struct Futures {
    enum Type { IMM, ASX };
};

%shared_ptr(OvernightIndexFuture)
class OvernightIndexFuture : public Instrument {
  public:
    OvernightIndexFuture(
        ext::shared_ptr<OvernightIndex> overnightIndex,
        const Date& valueDate,
        const Date& maturityDate,
        Handle<Quote> convexityAdjustment = Handle<Quote>(),
        OvernightAveraging::Type averagingMethod = OvernightAveraging::Compound);

    Real convexityAdjustment() const;
};


#endif
