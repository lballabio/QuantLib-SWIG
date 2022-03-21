
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006 StatPro Italia srl
 Copyright (C) 2015 Matthias Groncki

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

#ifndef quantlib_discount_curve_i
#define quantlib_discount_curve_i

%include termstructures.i
%include interpolation.i

%{
using QuantLib::InterpolatedDiscountCurve;
%}

%shared_ptr(InterpolatedDiscountCurve<LogLinear>);
%shared_ptr(InterpolatedDiscountCurve<MonotonicLogCubic>);
%shared_ptr(InterpolatedDiscountCurve<SplineCubic>);
%shared_ptr(InterpolatedDiscountCurve<KrugerLog>);
%shared_ptr(InterpolatedDiscountCurve<SplineLogCubic>);

template <class Interpolator>
class InterpolatedDiscountCurve : public YieldTermStructure {
  public:
    InterpolatedDiscountCurve(const std::vector<Date>& dates,
                              const std::vector<DiscountFactor>& discounts,
                              const DayCounter& dayCounter,
                              const Calendar& calendar = Calendar(),
                              const Interpolator& i = Interpolator());
    const std::vector<Time>& times() const;
    const std::vector<Real>& data() const;
    const std::vector<Date>& dates() const;
    const std::vector<DiscountFactor>& discounts() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,DiscountFactor> > nodes() const;
    #endif
};

%template(DiscountCurve) InterpolatedDiscountCurve<LogLinear>;
%template(MonotonicLogCubicDiscountCurve) InterpolatedDiscountCurve<MonotonicLogCubic>;
%template(NaturalCubicDiscountCurve) InterpolatedDiscountCurve<SplineCubic>;
%template(KrugerLogDiscountCurve) InterpolatedDiscountCurve<KrugerLog>;
%template(NaturalLogCubicDiscountCurve) InterpolatedDiscountCurve<SplineLogCubic>;

#endif
