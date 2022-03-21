
/*
 Copyright (C) 2005, 2006 StatPro Italia srl
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

#ifndef quantlib_zero_curve_i
#define quantlib_zero_curve_i

%include termstructures.i
%include interpolation.i

%{
using QuantLib::InterpolatedZeroCurve;
%}

%shared_ptr(InterpolatedZeroCurve<Linear>);
%shared_ptr(InterpolatedZeroCurve<LogLinear>);
%shared_ptr(InterpolatedZeroCurve<Cubic>);
%shared_ptr(InterpolatedZeroCurve<SplineCubic>);
%shared_ptr(InterpolatedZeroCurve<DefaultLogCubic>);
%shared_ptr(InterpolatedZeroCurve<MonotonicCubic>);
%shared_ptr(InterpolatedZeroCurve<Kruger>);

template <class Interpolator>
class InterpolatedZeroCurve : public YieldTermStructure {
  public:
    InterpolatedZeroCurve(const std::vector<Date>& dates,
                          const std::vector<Rate>& yields,
                          const DayCounter& dayCounter,
                          const Calendar& calendar = Calendar(),
                          const Interpolator& i = Interpolator(),
                          Compounding compounding = Continuous,
                          Frequency frequency = Annual);
    const std::vector<Time>& times() const;
    const std::vector<Real>& data() const;
    const std::vector<Date>& dates() const;
    const std::vector<Rate>& zeroRates() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Rate> > nodes() const;
    #endif
};

%template(ZeroCurve) InterpolatedZeroCurve<Linear>;
%template(LogLinearZeroCurve) InterpolatedZeroCurve<LogLinear>;
%template(CubicZeroCurve) InterpolatedZeroCurve<Cubic>;
%template(NaturalCubicZeroCurve) InterpolatedZeroCurve<SplineCubic>;
%template(LogCubicZeroCurve) InterpolatedZeroCurve<DefaultLogCubic>;
%template(MonotonicCubicZeroCurve) InterpolatedZeroCurve<MonotonicCubic>;
%template(KrugerZeroCurve) InterpolatedZeroCurve<Kruger>;


#endif
