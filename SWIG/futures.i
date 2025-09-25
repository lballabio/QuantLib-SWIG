/*
 Copyright (C) 2008 StatPro Italia srl
 Copyright (C) 2025 Hiroto Ogawa

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

#ifndef quantlib_futures_i
#define quantlib_futures_i

%include linearalgebra.i
%include forward.i
%include indexes.i
%include options.i

%{
using QuantLib::Futures;
using QuantLib::OvernightIndexFuture;
using QuantLib::PerpetualFutures;
using QuantLib::DiscountingPerpetualFuturesEngine;
%}

struct Futures {
    enum Type { IMM, ASX, Custom };
};

%shared_ptr(OvernightIndexFuture)
class OvernightIndexFuture : public Instrument {
  public:
    OvernightIndexFuture(
        ext::shared_ptr<OvernightIndex> overnightIndex,
        const Date& valueDate,
        const Date& maturityDate,
        Handle<Quote> convexityAdjustment = Handle<Quote>(),
        RateAveraging::Type averagingMethod = RateAveraging::Compound);

    Real convexityAdjustment() const;
};


%shared_ptr(PerpetualFutures)
class PerpetualFutures : public Instrument {
  public:
    enum PayoffType { Linear, Inverse };
    enum FundingType { FundingWithPreviousSpot, FundingWithCurrentSpot };
    PerpetualFutures(
        PerpetualFutures::PayoffType payoffType,
        PerpetualFutures::FundingType fundingType = FundingWithCurrentSpot,
        Period fundingFrequency = Period(8, Hours),
        Calendar cal = NullCalendar(),
        DayCounter dc = ActualActual(ActualActual::ISDA));
};

%shared_ptr(DiscountingPerpetualFuturesEngine)
class DiscountingPerpetualFuturesEngine : public PricingEngine {
  public:
    enum InterpolationType { PiecewiseConstant, Linear, CubicSpline };
    DiscountingPerpetualFuturesEngine(
        const Handle<YieldTermStructure>& domesticDiscountCurve,
        const Handle<YieldTermStructure>& foreignDiscountCurve,
        const Handle<Quote>& assetSpot,
        const std::vector<Time>& fundingTimes,
        const std::vector<Rate>& fundingRates,
        const std::vector<Spread>& interestRateDiffs,
        const InterpolationType fundingInterpType = PiecewiseConstant,
        const Real maxT = 60.);
};



#endif
