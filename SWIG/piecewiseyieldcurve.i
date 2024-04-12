
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
using QuantLib::SimpleZeroYield;
%}

struct Discount {};
struct ZeroYield {};
struct ForwardRate {};
struct SimpleZeroYield {};

// curve

%{
using QuantLib::PiecewiseYieldCurve;
%}

%{
using QuantLib::IterativeBootstrap;
using QuantLib::GlobalBootstrap;
%}

%{
struct _IterativeBootstrap {
    double accuracy, minValue, maxValue;
    Size maxAttempts;
    Real maxFactor, minFactor;
    bool dontThrow;
    Size dontThrowSteps, maxEvaluations;
    _IterativeBootstrap(double accuracy = Null<double>(),
                        double minValue = Null<double>(),
                        double maxValue = Null<double>(),
                        Size maxAttempts = 1,
                        Real maxFactor = 2.0,
                        Real minFactor = 2.0,
                        bool dontThrow = false,
                        Size dontThrowSteps = 10,
                        Size maxEvaluations = 100)
    : accuracy(accuracy), minValue(minValue), maxValue(maxValue),
      maxAttempts(maxAttempts), maxFactor(maxFactor), minFactor(minFactor),
      dontThrow(dontThrow), dontThrowSteps(dontThrowSteps),
      maxEvaluations(maxEvaluations) {}
};

template <class PiecewiseYieldCurve>
inline typename PiecewiseYieldCurve::bootstrap_type make_bootstrap(const _IterativeBootstrap& b) {
    return {
        b.accuracy, b.minValue, b.maxValue,
        b.maxAttempts, b.maxFactor, b.minFactor,
        b.dontThrow, b.dontThrowSteps,
        b.maxEvaluations
    };
}
%}

%rename(IterativeBootstrap) _IterativeBootstrap;
struct _IterativeBootstrap {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") _IterativeBootstrap;
    #endif
    _IterativeBootstrap(doubleOrNull accuracy = Null<double>(),
                        doubleOrNull minValue = Null<double>(),
                        doubleOrNull maxValue = Null<double>(),
                        Size maxAttempts = 1,
                        Real maxFactor = 2.0,
                        Real minFactor = 2.0,
                        bool dontThrow = false,
                        Size dontThrowSteps = 10,
                        Size maxEvaluations = 100);
};

// global boostrapper

%{
class AdditionalErrors {
    std::vector<ext::shared_ptr<RateHelper> > additionalHelpers_;
  public:
    AdditionalErrors(const std::vector<ext::shared_ptr<RateHelper> >& additionalHelpers)
    : additionalHelpers_(additionalHelpers) {}
    Array operator()() const {
        Array errors(additionalHelpers_.size() - 2);
        Real a = additionalHelpers_.front()->impliedQuote();
        Real b = additionalHelpers_.back()->impliedQuote();
        for (Size k = 0; k < errors.size(); ++k) {
            errors[k] = (static_cast<Real>(errors.size()-k) * a + static_cast<Real>(1+k) * b) / static_cast<Real>(errors.size()+1)
                - additionalHelpers_.at(1+k)->impliedQuote();
        }
        return errors;
    }
};

class AdditionalDates {
    std::vector<Date> additionalDates_;
  public:
    AdditionalDates(const std::vector<Date>& additionalDates)
    : additionalDates_(additionalDates) {}
    std::vector<Date> operator()() const {
        return additionalDates_;
    }
};

struct _GlobalBootstrap {
    std::vector<ext::shared_ptr<RateHelper> > additionalHelpers;
    std::vector<Date> additionalDates;
    double accuracy;
    _GlobalBootstrap(double accuracy = Null<double>())
    : accuracy(accuracy) {}
   _GlobalBootstrap(const std::vector<ext::shared_ptr<RateHelper> >& additionalHelpers,
                     const std::vector<Date>& additionalDates,
                     double accuracy = Null<double>())
    : additionalHelpers(additionalHelpers), additionalDates(additionalDates), accuracy(accuracy) {}
};

template <class PiecewiseYieldCurve>
inline typename PiecewiseYieldCurve::bootstrap_type make_bootstrap(const _GlobalBootstrap& b) {
    if (b.additionalHelpers.empty()) { 
        return { b.accuracy };
    } else {
        return { b.additionalHelpers,
                 AdditionalDates(b.additionalDates),
                 AdditionalErrors(b.additionalHelpers),
                 b.accuracy 
        };
    }
}

%}

%rename(GlobalBootstrap) _GlobalBootstrap;
struct _GlobalBootstrap {
    _GlobalBootstrap(doubleOrNull accuracy = Null<double>());
    _GlobalBootstrap(const std::vector<ext::shared_ptr<RateHelper> >& additionalHelpers,
                     const std::vector<Date>& additionalDates,
                     doubleOrNull accuracy = Null<double>());
};


/* We have to resort to a macro, because the R implementation of shared_ptr
   can't take class templates with two or more template arguments. */

%define export_piecewise_curve(Name,Traits,Interpolator,BootstrapTemplate,BootstrapProxy)

%{
typedef PiecewiseYieldCurve<Traits, Interpolator, BootstrapTemplate> Name;
%}

%shared_ptr(Name);
class Name : public YieldTermStructure {
  public:
    %extend {
        Name(const Date& referenceDate,
             const std::vector<ext::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const std::vector<Handle<Quote> >& jumps = std::vector<Handle<Quote> >(),
             const std::vector<Date>& jumpDates = std::vector<Date>(),
             const Interpolator& i = Interpolator(),
             const BootstrapProxy& b = BootstrapProxy()) {
            return new Name(referenceDate, instruments, dayCounter, jumps, jumpDates,
                            i, make_bootstrap<Name>(b));
        }
        Name(Integer settlementDays,
             const Calendar& calendar,
             const std::vector<ext::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const std::vector<Handle<Quote> >& jumps = std::vector<Handle<Quote> >(),
             const std::vector<Date>& jumpDates = std::vector<Date>(),
             const Interpolator& i = Interpolator(),
             const BootstrapProxy& b = BootstrapProxy()) {
            return new Name(settlementDays, calendar, instruments, dayCounter,
                            jumps, jumpDates, i, make_bootstrap<Name>(b));
        }
        Name(const Date& referenceDate,
             const std::vector<ext::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const BootstrapProxy& b,
             const Interpolator& i = Interpolator()) {
            return new Name(referenceDate, instruments, dayCounter, i,
                            make_bootstrap<Name>(b));
        }
        Name(Integer settlementDays, const Calendar& calendar,
             const std::vector<ext::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const BootstrapProxy& b,
             const Interpolator& i = Interpolator()) {
            return new Name(settlementDays, calendar, instruments, dayCounter,
                            i, make_bootstrap<Name>(b));
        }
    }
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif

    void recalculate();
    void freeze();
    void unfreeze();
};

%enddef

export_piecewise_curve(PiecewiseFlatForward,ForwardRate,BackwardFlat,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseLogLinearDiscount,Discount,LogLinear,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseLinearForward,ForwardRate,Linear,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseLinearZero,ZeroYield,Linear,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseCubicZero,ZeroYield,Cubic,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseLogCubicDiscount,Discount,MonotonicLogCubic,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseSplineCubicDiscount,Discount,SplineCubic,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseKrugerZero,ZeroYield,Kruger,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseKrugerLogDiscount,Discount,KrugerLog,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseConvexMonotoneForward,ForwardRate,ConvexMonotone,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseConvexMonotoneZero,ZeroYield,ConvexMonotone,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseNaturalCubicZero,ZeroYield,SplineCubic,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseNaturalLogCubicDiscount,Discount,SplineLogCubic,IterativeBootstrap,_IterativeBootstrap);
export_piecewise_curve(PiecewiseLogMixedLinearCubicDiscount,Discount,LogMixedLinearCubic,IterativeBootstrap,_IterativeBootstrap);

export_piecewise_curve(GlobalLinearSimpleZeroCurve,SimpleZeroYield,Linear,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseFlatForward,ForwardRate,BackwardFlat,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseLogLinearDiscount,Discount,LogLinear,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseLinearForward,ForwardRate,Linear,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseLinearZero,ZeroYield,Linear,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseCubicZero,ZeroYield,Cubic,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseLogCubicDiscount,Discount,MonotonicLogCubic,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseSplineCubicDiscount,Discount,SplineCubic,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseKrugerZero,ZeroYield,Kruger,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseKrugerLogDiscount,Discount,KrugerLog,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseConvexMonotoneForward,ForwardRate,ConvexMonotone,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseConvexMonotoneZero,ZeroYield,ConvexMonotone,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseNaturalCubicZero,ZeroYield,SplineCubic,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseNaturalLogCubicDiscount,Discount,SplineLogCubic,GlobalBootstrap,_GlobalBootstrap);
export_piecewise_curve(GlobalPiecewiseLogMixedLinearCubicDiscount,Discount,LogMixedLinearCubic,GlobalBootstrap,_GlobalBootstrap);

#endif
