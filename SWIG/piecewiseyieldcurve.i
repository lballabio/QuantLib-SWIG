
/*
 Copyright (C) 2005, 2006, 2007, 2008 StatPro Italia srl
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

#ifndef quantlib_piecewise_yield_curve_i
#define quantlib_piecewise_yield_curve_i

%include termstructures.i
%include ratehelpers.i
%include interpolation.i
%include null.i

// bootstrap traits

%{
using QuantLib::Discount;
using QuantLib::ZeroYield;
using QuantLib::ForwardRate;
%}

struct Discount {};
struct ZeroYield {};
struct ForwardRate {};

// curve

%{
using QuantLib::PiecewiseYieldCurve;
%}

/* We have to resort to a macro, because the R implementation of shared_ptr
   can't take class templates with two or more template arguments. */

%{
struct IterativeBootstrap {
    double accuracy, minValue, maxValue;
    IterativeBootstrap(double accuracy = Null<double>(),
                       double minValue = Null<double>(),
                       double maxValue = Null<double>())
    : accuracy(accuracy), minValue(minValue), maxValue(maxValue) {}
};
%}

struct IterativeBootstrap {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") IterativeBootstrap;
    #endif
    IterativeBootstrap(doubleOrNull accuracy = Null<double>(),
                       doubleOrNull minValue = Null<double>(),
                       doubleOrNull maxValue = Null<double>());
};

%define export_piecewise_curve(Name,Traits,Interpolator)

%{
typedef PiecewiseYieldCurve<Traits, Interpolator> Name;
%}

%shared_ptr(Name);
class Name : public YieldTermStructure {
  public:
    %extend {
        Name(const Date& referenceDate,
             const std::vector<boost::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const std::vector<Handle<Quote> >& jumps = std::vector<Handle<Quote> >(),
             const std::vector<Date>& jumpDates = std::vector<Date>(),
             Real accuracy = 1.0e-12,
             const Interpolator& i = Interpolator(),
             const IterativeBootstrap& b = IterativeBootstrap()) {
            return new Name(referenceDate, instruments, dayCounter, jumps, jumpDates,
                            accuracy, i, Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
        Name(Integer settlementDays, const Calendar& calendar,
             const std::vector<boost::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const std::vector<Handle<Quote> >& jumps = std::vector<Handle<Quote> >(),
             const std::vector<Date>& jumpDates = std::vector<Date>(),
             Real accuracy = 1.0e-12,
             const Interpolator& i = Interpolator(),
             const IterativeBootstrap& b = IterativeBootstrap()) {
            return new Name(settlementDays, calendar, instruments, dayCounter,
                            jumps, jumpDates, accuracy, Interpolator(),
                            Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
        Name(const Date& referenceDate,
             const std::vector<boost::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const IterativeBootstrap& b) {
            return new Name(referenceDate, instruments, dayCounter, Interpolator(),
                            Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
        Name(Integer settlementDays, const Calendar& calendar,
             const std::vector<boost::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const IterativeBootstrap& b) {
            return new Name(settlementDays, calendar, instruments, dayCounter,
                            Interpolator(), Name::bootstrap_type(b.accuracy, b.minValue, b.maxValue));
        }
    }
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};

%enddef


export_piecewise_curve(PiecewiseFlatForward,ForwardRate,BackwardFlat);
export_piecewise_curve(PiecewiseLogLinearDiscount,Discount,LogLinear);
export_piecewise_curve(PiecewiseLinearForward,ForwardRate,Linear);
export_piecewise_curve(PiecewiseLinearZero,ZeroYield,Linear);
export_piecewise_curve(PiecewiseCubicZero,ZeroYield,Cubic);
export_piecewise_curve(PiecewiseLogCubicDiscount,Discount,MonotonicLogCubic);
export_piecewise_curve(PiecewiseSplineCubicDiscount,Discount,SplineCubic);
export_piecewise_curve(PiecewiseKrugerZero,ZeroYield,Kruger);
export_piecewise_curve(PiecewiseKrugerLogDiscount,Discount,KrugerLog);

#endif
