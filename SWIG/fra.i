/*
 Copyright (C) 2012 Tawanda Gwena
 Copyright (C) 2012 Francis Duffy
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

#ifndef quantlib_forward_rate_agreement_i
#define quantlib_forward_rate_agreement_i
 
%include instruments.i
%include termstructures.i
%include interestrate.i
%include forward.i

%{
using QuantLib::Position;
using QuantLib::ForwardRateAgreement;
%}
 
struct Position {
    enum Type { Long, Short };
};

%shared_ptr(ForwardRateAgreement)
class ForwardRateAgreement : public Instrument {
  public:
    ForwardRateAgreement(const Date& valueDate,
                         const Date& maturityDate,
                         Position::Type type,
                         Rate strikeForwardRate,
                         Real notionalAmount,
                         const ext::shared_ptr<IborIndex>& index,
                         const Handle<YieldTermStructure>& discountCurve = Handle<YieldTermStructure>(),
                         bool useIndexedCoupon = true);

    ForwardRateAgreement(const Date& valueDate,
                         Position::Type type,
                         Rate strikeForwardRate,
                         Real notionalAmount,
                         const ext::shared_ptr<IborIndex>& index,
                         const Handle<YieldTermStructure>& discountCurve = Handle<YieldTermStructure>());

    Real amount() const;
    Date fixingDate() const;
    InterestRate forwardRate() const;
};
 

#endif
