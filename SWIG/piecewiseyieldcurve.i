
/*
 Copyright (C) 2005, 2006, 2007, 2008 StatPro Italia srl
 Copyright (C) 2018 Matthias Lungwitz

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

#ifndef quantlib_piecewise_yield_curve_i
#define quantlib_piecewise_yield_curve_i

%include termstructures.i
%include ratehelpers.i
%include interpolation.i
%include optimizers.i
%include null.i

%{
using QuantLib::Discount;
using QuantLib::ZeroYield;
using QuantLib::ForwardRate;
using QuantLib::PiecewiseYieldCurve;
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

/* We have to resort to a macro, because the R implementation of shared_ptr
   can't take class templates with two or more template arguments. */

%define export_piecewise_curve(Name,Traits,Interpolator)

%{
typedef PiecewiseYieldCurve<Traits, Interpolator> Name;
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
             const _IterativeBootstrap& b = _IterativeBootstrap()) {
            return new Name(referenceDate, instruments, dayCounter, jumps, jumpDates,
                            i, make_bootstrap<Name>(b));
        }
        Name(Integer settlementDays, const Calendar& calendar,
             const std::vector<ext::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const std::vector<Handle<Quote> >& jumps = std::vector<Handle<Quote> >(),
             const std::vector<Date>& jumpDates = std::vector<Date>(),
             const Interpolator& i = Interpolator(),
             const _IterativeBootstrap& b = _IterativeBootstrap()) {
            return new Name(settlementDays, calendar, instruments, dayCounter,
                            jumps, jumpDates, i, make_bootstrap<Name>(b));
        }
        Name(const Date& referenceDate,
             const std::vector<ext::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const _IterativeBootstrap& b,
             const Interpolator& i = Interpolator()) {
            return new Name(referenceDate, instruments, dayCounter, i,
                            make_bootstrap<Name>(b));
        }
        Name(Integer settlementDays, const Calendar& calendar,
             const std::vector<ext::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const _IterativeBootstrap& b,
             const Interpolator& i = Interpolator()) {
            return new Name(settlementDays, calendar, instruments, dayCounter,
                            i, make_bootstrap<Name>(b));
        }
    }
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    const std::vector<Real>& data() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif

    void recalculate();
    void freeze();
    void unfreeze();
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
export_piecewise_curve(PiecewiseConvexMonotoneForward,ForwardRate,ConvexMonotone);
export_piecewise_curve(PiecewiseConvexMonotoneZero,ZeroYield,ConvexMonotone);
export_piecewise_curve(PiecewiseNaturalCubicZero,ZeroYield,SplineCubic);
export_piecewise_curve(PiecewiseNaturalLogCubicDiscount,Discount,SplineLogCubic);
export_piecewise_curve(PiecewiseLogMixedLinearCubicDiscount,Discount,LogMixedLinearCubic);
export_piecewise_curve(PiecewiseParabolicCubicZero,ZeroYield,ParabolicCubic);
export_piecewise_curve(PiecewiseMonotonicParabolicCubicZero,ZeroYield,MonotonicParabolicCubic);
export_piecewise_curve(PiecewiseLogParabolicCubicDiscount,Discount,LogParabolicCubic);
export_piecewise_curve(PiecewiseMonotonicLogParabolicCubicDiscount,Discount,MonotonicLogParabolicCubic);


// global boostrapper
// hard-coded to linearly-interpolated, simply-compounded zero rates for now

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
    ext::shared_ptr<OptimizationMethod> optimizer;
    ext::shared_ptr<EndCriteria> endCriteria;
    _GlobalBootstrap(double accuracy = Null<double>(),
                     ext::shared_ptr<OptimizationMethod> optimizer = nullptr,
                     ext::shared_ptr<EndCriteria> endCriteria = nullptr)
    : accuracy(accuracy), optimizer(optimizer), endCriteria(endCriteria) {}
   _GlobalBootstrap(const std::vector<ext::shared_ptr<RateHelper> >& additionalHelpers,
                    const std::vector<Date>& additionalDates,
                    double accuracy = Null<double>(),
                    ext::shared_ptr<OptimizationMethod> optimizer = nullptr,
                    ext::shared_ptr<EndCriteria> endCriteria = nullptr)
   : additionalHelpers(additionalHelpers), additionalDates(additionalDates), accuracy(accuracy),
     optimizer(optimizer), endCriteria(endCriteria) {}
};
%}

%rename(GlobalBootstrap) _GlobalBootstrap;
struct _GlobalBootstrap {
    _GlobalBootstrap(doubleOrNull accuracy = Null<double>(),
                     ext::shared_ptr<OptimizationMethod> optimizer = nullptr,
                     ext::shared_ptr<EndCriteria> endCriteria = nullptr);
    _GlobalBootstrap(const std::vector<ext::shared_ptr<RateHelper> >& additionalHelpers,
                     const std::vector<Date>& additionalDates,
                     doubleOrNull accuracy = Null<double>(),
                     ext::shared_ptr<OptimizationMethod> optimizer = nullptr,
                     ext::shared_ptr<EndCriteria> endCriteria = nullptr);
};


%{
using QuantLib::SimpleZeroYield;
typedef PiecewiseYieldCurve<SimpleZeroYield, Linear, QuantLib::GlobalBootstrap>
    GlobalLinearSimpleZeroCurve;
%}

%shared_ptr(GlobalLinearSimpleZeroCurve);
class GlobalLinearSimpleZeroCurve : public YieldTermStructure {
  public:
    %extend {
        GlobalLinearSimpleZeroCurve(
             const Date& referenceDate,
             const std::vector<ext::shared_ptr<RateHelper> >& instruments,
             const DayCounter& dayCounter,
             const _GlobalBootstrap& b) {
            if (b.additionalHelpers.empty()) {
                return new GlobalLinearSimpleZeroCurve(
                    referenceDate, instruments, dayCounter, Linear(),
                    GlobalLinearSimpleZeroCurve::bootstrap_type(b.accuracy, b.optimizer, b.endCriteria));
            } else {
                return new GlobalLinearSimpleZeroCurve(
                    referenceDate, instruments, dayCounter, Linear(),
                    GlobalLinearSimpleZeroCurve::bootstrap_type(b.additionalHelpers,
                                                                AdditionalDates(b.additionalDates),
                                                                AdditionalErrors(b.additionalHelpers),
                                                                b.accuracy, b.optimizer, b.endCriteria));
            }
        }
    }
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif
};


%{
using QuantLib::PiecewiseSpreadYieldCurve;
%}

%define export_piecewise_spread_curve(Name,Traits,Interpolator)

%{
typedef PiecewiseSpreadYieldCurve<Traits, Interpolator> Name;
%}

%shared_ptr(Name);
class Name : public YieldTermStructure {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") Name;
    #endif
  public:
    %extend {
        Name(const Handle<YieldTermStructure>& baseCurve,
             const std::vector<ext::shared_ptr<RateHelper> >& instruments,
             const Interpolator& interpolator = Interpolator(),
             const _IterativeBootstrap& bootstrap = _IterativeBootstrap()) {
            return new Name(baseCurve, instruments, interpolator, make_bootstrap<Name>(bootstrap));
        }
    }
    const Handle<YieldTermStructure>& baseCurve() const;
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    const std::vector<Real>& data() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,Real> > nodes() const;
    #endif

    void recalculate();
    void freeze();
    void unfreeze();
};

%enddef

export_piecewise_spread_curve(PiecewiseLogLinearDiscountSpread,Discount,LogLinear);
export_piecewise_spread_curve(PiecewiseLogCubicDiscountSpread,Discount,MonotonicLogCubic);
export_piecewise_spread_curve(PiecewiseNaturalLogCubicDiscountSpread,Discount,SplineLogCubic);
export_piecewise_spread_curve(PiecewiseLogMixedLinearCubicDiscountSpread,Discount,LogMixedLinearCubic);

#endif
