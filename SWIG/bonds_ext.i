/*
 Copyright (C) 2004, 2005, 2006, 2007 StatPro Italia srl
 Copyright (C) 2009 Joseph Malicki
 Copyright (C) 2011 Lluis Pujol Bajador
 Copyright (C) 2014 Simon Mazzucca
 Copyright (C) 2016 Gouthaman Balaraman

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

#ifndef quantlib_bonds_ext_i
#define quantlib_bonds_ext_i

%include bonds.i

%{
#include <qlext/instruments/bonds/ctbzerobond.hpp>
#include <qlext/instruments/bonds/ctbfixedbond.hpp>

using QuantLib::CTBZeroBond;
using QuantLib::CTBFixedBond;

typedef boost::shared_ptr<Instrument> CTBZeroBondPtr;
typedef boost::shared_ptr<Instrument> CTBFixedBondPtr;
%}

%rename(CTBZeroBond) CTBZeroBondPtr;
class CTBZeroBondPtr : public BondPtr {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CTBZeroBondPtr;
    #endif
  public:
    %extend {
        CTBZeroBondPtr(
                Natural settlementDays,
                const Calendar &calendar,
                Real faceAmount,
                Real issuePrice,
                const Date& issueDate,
                const DayCounter& accrualDayCounter,
                const Date& maturityDate,
                BusinessDayConvention paymentConvention = QuantLib::Following,
                Real redemption = 100.0) {
            return new CTBZeroBondPtr(
                new CTBZeroBond(settlementDays, calendar, faceAmount, issuePrice,
                                issueDate, accrualDayCounter, maturityDate,
                                paymentConvention, redemption));
        }
    }
};

%rename(CTBFixedBond) CTBFixedBondPtr;
class CTBFixedBondPtr : public BondPtr {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") CTBFixedBondPtr;
    #endif
  public:
    %extend {
        CTBFixedBondPtr(
                const Date& issueDate,
                Natural settlementDays,
                Real faceAmount,
                const Date& startDate,
                const Date& maturity,
                const Period& tenor,
                const std::vector<Real>& coupons,
                const DayCounter& accrualDayCounter = QuantLib::ActualActual(QuantLib::ActualActual::ISMA),
                bool endOfMonth = false,
                const Calendar& calendar = NullCalendar(),
                const Calendar& paymentCalendar = QuantLib::China(QuantLib::China::IB),
                BusinessDayConvention convention = Unadjusted,
                DateGeneration::Rule rule = DateGeneration::Backward,
                BusinessDayConvention paymentConvention = Unadjusted,
                Real redemption = 100.0,
                const Date& firstDate = Date(),
                const Date& nextToLastDate = Date()) {
            return new CTBFixedBondPtr(
                new CTBFixedBond(issueDate, settlementDays, faceAmount, startDate,
                                 maturity, tenor, coupons, accrualDayCounter, endOfMonth,
                                 calendar, paymentCalendar, convention, rule, paymentConvention,
                                 redemption, firstDate, nextToLastDate));
        }
    }
};


#endif
