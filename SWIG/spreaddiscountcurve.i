#ifndef quantlib_spread_discount_curve_i
#define quantlib_spread_discount_curve_i

%include termstructures.i
%include interpolation.i

%{
using QuantLib::InterpolatedSpreadDiscountCurve;
%}

%shared_ptr(InterpolatedSpreadDiscountCurve<LogLinear>);
%shared_ptr(InterpolatedSpreadDiscountCurve<MonotonicLogCubic>);
%shared_ptr(InterpolatedSpreadDiscountCurve<SplineLogCubic>);
%shared_ptr(InterpolatedSpreadDiscountCurve<LogMixedLinearCubic>);

template <class Interpolator>
class InterpolatedSpreadDiscountCurve : public YieldTermStructure {
  public:
    InterpolatedSpreadDiscountCurve(const Handle<YieldTermStructure>& baseCurve,
                                    const std::vector<Date>& dates,
                                    const std::vector<DiscountFactor>& dfs,
                                    const Interpolator& i = Interpolator());
    const Handle<YieldTermStructure>& baseCurve() const;
    const std::vector<Time>& times() const;
    const std::vector<Date>& dates() const;
    const std::vector<Real>& data() const;
    #if !defined(SWIGR)
    std::vector<std::pair<Date,DiscountFactor> > nodes() const;
    #endif
};

%template(SpreadDiscountCurve) InterpolatedSpreadDiscountCurve<LogLinear>;
%template(MonotonicLogCubicSpreadDiscountCurve) InterpolatedSpreadDiscountCurve<MonotonicLogCubic>;
%template(NaturalLogCubicSpreadDiscountCurve) InterpolatedSpreadDiscountCurve<SplineLogCubic>;
%template(LogMixedLinearCubicSpreadDiscountCurve) InterpolatedSpreadDiscountCurve<LogMixedLinearCubic>;

#endif
