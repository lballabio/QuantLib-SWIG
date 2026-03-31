
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2002, 2003 Ferdinando Ametrano
 Copyright (C) 2003, 2004, 2008 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2018, 2020 Matthias Lungwitz
 
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

#ifndef quantlib_interpolation_i
#define quantlib_interpolation_i

%include linearalgebra.i
%include optimizers.i

%{
// safe versions which copy their arguments
template <class I>
class SafeInterpolation {
  public:
    SafeInterpolation(const Array& x, const Array& y)
    : x_(x), y_(y), f_(x_.begin(), x_.end(), y_.begin()) {}
    Real operator()(Real x, bool allowExtrapolation = false) const {
        return f_(x, allowExtrapolation);
    }
    Real primitive(Real x, bool allowExtrapolation = false) const {
        return f_.primitive(x, allowExtrapolation);
    }
    Real derivative(Real x, bool allowExtrapolation = false) const {
        return f_.derivative(x, allowExtrapolation);
    }
    Real secondDerivative(Real x, bool allowExtrapolation = false) const {
        return f_.secondDerivative(x, allowExtrapolation);
    }
    Array x_, y_;
    I f_;
};
%}

%define make_safe_interpolation(T)
%{
typedef SafeInterpolation<QuantLib::T> Safe##T;
%}
%rename(T) Safe##T;
class Safe##T {
    #if defined(SWIGCSHARP)
    %rename(call) operator();
    #endif
  public:
    Safe##T(const Array& x, const Array& y);
    Real operator()(Real x, bool allowExtrapolation = false) const;
    Real derivative(Real x, bool extrapolate = false) const;
    Real secondDerivative(Real x, bool extrapolate = false) const;
    Real primitive(Real x, bool extrapolate = false) const;
};
%enddef

make_safe_interpolation(LinearInterpolation);
make_safe_interpolation(LogLinearInterpolation);

make_safe_interpolation(BackwardFlatInterpolation);
make_safe_interpolation(ForwardFlatInterpolation);

make_safe_interpolation(CubicNaturalSpline);
make_safe_interpolation(LogCubicNaturalSpline);
make_safe_interpolation(MonotonicCubicNaturalSpline);
make_safe_interpolation(MonotonicLogCubicNaturalSpline);

make_safe_interpolation(KrugerCubic);
make_safe_interpolation(KrugerLogCubic);

make_safe_interpolation(FritschButlandCubic);
make_safe_interpolation(FritschButlandLogCubic);

make_safe_interpolation(Parabolic);
make_safe_interpolation(LogParabolic);
make_safe_interpolation(MonotonicParabolic);
make_safe_interpolation(MonotonicLogParabolic);

make_safe_interpolation(LagrangeInterpolation);
%{
// safe versions which copy their arguments
template <class I>
class SafeInterpolation2D {
  public:
    SafeInterpolation2D(const Array& x, const Array& y, const Matrix& m)
    : x_(x), y_(y), m_(m), f_(x_.begin(), x_.end(), y_.begin(), y_.end(), m_) {}
    Real operator()(Real x, Real y, bool allowExtrapolation = false) const {
        return f_(x,y, allowExtrapolation);
    }
  protected:
    Array x_, y_;
    Matrix m_;
    I f_;
};
%}

%define make_safe_interpolation2d(T)
%{
typedef SafeInterpolation2D<QuantLib::T> Safe##T;
%}
%rename(T) Safe##T;
class Safe##T {
    #if defined(SWIGCSHARP)
    %rename(call) operator();
    #endif
  public:
    Safe##T(const Array& x, const Array& y, const Matrix& m);
    Real operator()(Real x, Real y, bool allowExtrapolation = false) const;
};
%enddef

make_safe_interpolation2d(BilinearInterpolation);
make_safe_interpolation2d(BicubicSpline);


// interpolation traits

%{
using QuantLib::CubicInterpolation;
using QuantLib::MixedInterpolation;
using QuantLib::BackwardFlat;
using QuantLib::ForwardFlat;
using QuantLib::Linear;
using QuantLib::LogLinear;
using QuantLib::Cubic;
using QuantLib::Bicubic;
using QuantLib::ConvexMonotone;
using QuantLib::DefaultLogCubic;

class MonotonicCubic : public Cubic {
  public:
    MonotonicCubic()
    : Cubic(CubicInterpolation::Spline, true,
            CubicInterpolation::SecondDerivative, 0.0,
            CubicInterpolation::SecondDerivative, 0.0) {}
};

class SplineCubic : public Cubic {
  public:
    SplineCubic()
    : Cubic(CubicInterpolation::Spline, false,
            CubicInterpolation::SecondDerivative, 0.0,
            CubicInterpolation::SecondDerivative, 0.0) {}
};

class Kruger : public Cubic {
  public:
    Kruger()
    : Cubic(CubicInterpolation::Kruger) {}
};

class LogCubic : public QuantLib::LogCubic {
  public:
    // We add defaults for all constructor arguments because wrappers for
    // InterpolatedDiscountCurve and PiecewiseYieldCurve assume that all
    // interpolators have default constructors.
    LogCubic(CubicInterpolation::DerivativeApprox da = CubicInterpolation::Spline,
             bool monotonic = true,
             CubicInterpolation::BoundaryCondition leftCondition
                 = CubicInterpolation::SecondDerivative,
             Real leftConditionValue = 0.0,
             CubicInterpolation::BoundaryCondition rightCondition
                 = CubicInterpolation::SecondDerivative,
             Real rightConditionValue = 0.0)
    : QuantLib::LogCubic(da, monotonic, leftCondition, leftConditionValue,
                         rightCondition, rightConditionValue) {}
};

class MonotonicLogCubic : public LogCubic {
  public:
    MonotonicLogCubic()
    : LogCubic(CubicInterpolation::Spline, true,
               CubicInterpolation::SecondDerivative, 0.0,
               CubicInterpolation::SecondDerivative, 0.0) {}
};

class KrugerLog : public LogCubic {
  public:
    KrugerLog()
    : LogCubic(CubicInterpolation::Kruger, false,
               CubicInterpolation::SecondDerivative, 0.0,
               CubicInterpolation::SecondDerivative, 0.0) {}
};

class SplineLogCubic : public LogCubic {
  public:
    SplineLogCubic()
    : LogCubic(CubicInterpolation::Spline, false,
               CubicInterpolation::SecondDerivative, 0.0,
               CubicInterpolation::SecondDerivative, 0.0) {}
};

class ParabolicCubic : public QuantLib::Cubic {
  public:
    ParabolicCubic()
    : QuantLib::Cubic(CubicInterpolation::Parabolic, false,
                      CubicInterpolation::SecondDerivative, 0.0,
                      CubicInterpolation::SecondDerivative, 0.0) {}
};

class MonotonicParabolicCubic : public QuantLib::Cubic {
  public:
    MonotonicParabolicCubic()
    : QuantLib::Cubic(CubicInterpolation::Parabolic, true,
                      CubicInterpolation::SecondDerivative, 0.0,
                      CubicInterpolation::SecondDerivative, 0.0) {}
};

class LogParabolicCubic : public LogCubic {
  public:
    LogParabolicCubic()
    : LogCubic(CubicInterpolation::Parabolic, false,
               CubicInterpolation::SecondDerivative, 0.0,
               CubicInterpolation::SecondDerivative, 0.0) {}
};

class MonotonicLogParabolicCubic : public LogCubic {
  public:
    MonotonicLogParabolicCubic()
    : LogCubic(CubicInterpolation::Parabolic, true,
               CubicInterpolation::SecondDerivative, 0.0,
               CubicInterpolation::SecondDerivative, 0.0) {}
};
%}

%nodefaultctor CubicInterpolation;
struct CubicInterpolation {
    enum DerivativeApprox {
        Spline,
        SplineOM1,
        SplineOM2,
        FourthOrder,
        Parabolic,
        FritschButland,
        Akima,
        Kruger,
        Harmonic,
    };
    enum BoundaryCondition {
        NotAKnot,
        FirstDerivative,
        SecondDerivative,
        Periodic,
        Lagrange,
    };
};

%nodefaultctor MixedInterpolation;
struct MixedInterpolation {
    enum Behavior { ShareRanges, SplitRanges };
};

struct BackwardFlat {};
struct ForwardFlat {};
struct Linear {};
struct LogLinear {};
struct Cubic {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") Cubic;
    #endif
    Cubic(CubicInterpolation::DerivativeApprox da = CubicInterpolation::Kruger,
          bool monotonic = false,
          CubicInterpolation::BoundaryCondition leftCondition
              = CubicInterpolation::SecondDerivative,
          doubleOrNull leftConditionValue = 0.0,
          CubicInterpolation::BoundaryCondition rightCondition
              = CubicInterpolation::SecondDerivative,
          doubleOrNull rightConditionValue = 0.0);
};
struct LogCubic {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") LogCubic;
    #endif
    LogCubic(CubicInterpolation::DerivativeApprox da = CubicInterpolation::Spline,
             bool monotonic = true,
             CubicInterpolation::BoundaryCondition leftCondition
                 = CubicInterpolation::SecondDerivative,
             doubleOrNull leftConditionValue = 0.0,
             CubicInterpolation::BoundaryCondition rightCondition
                 = CubicInterpolation::SecondDerivative,
             doubleOrNull rightConditionValue = 0.0);
};
struct Bicubic {};
struct MonotonicCubic : Cubic {};
struct DefaultLogCubic {};
struct MonotonicLogCubic : LogCubic {};
struct SplineCubic : Cubic {};
struct SplineLogCubic : LogCubic {};
struct Kruger : Cubic {};
struct KrugerLog : LogCubic {};
struct ConvexMonotone {
    ConvexMonotone(Real quadraticity = 0.3,
                   Real monotonicity = 0.7,
                   bool forcePositive = true);
};
struct ParabolicCubic : Cubic {};
struct MonotonicParabolicCubic : Cubic {};
struct LogParabolicCubic : LogCubic {};
struct MonotonicLogParabolicCubic : LogCubic {};

%define make_mixed_linear_cubic(T)
%{
class T : public QuantLib::T {
  public:
    // We add defaults for all constructor arguments because wrappers for
    // InterpolatedDiscountCurve and PiecewiseYieldCurve assume that all
    // interpolators have default constructors.
    T(Size n = 0,
      MixedInterpolation::Behavior behavior = MixedInterpolation::ShareRanges,
      CubicInterpolation::DerivativeApprox da = CubicInterpolation::Spline,
      bool monotonic = true,
      CubicInterpolation::BoundaryCondition leftCondition
          = CubicInterpolation::SecondDerivative,
      Real leftConditionValue = 0.0,
      CubicInterpolation::BoundaryCondition rightCondition
          = CubicInterpolation::SecondDerivative,
      Real rightConditionValue = 0.0)
    : QuantLib::T(n, behavior, da, monotonic, leftCondition, leftConditionValue,
                  rightCondition, rightConditionValue) {}
};
%}

struct T {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") T;
    #endif
    T(Size n = 0,
      MixedInterpolation::Behavior behavior = MixedInterpolation::ShareRanges,
      CubicInterpolation::DerivativeApprox da = CubicInterpolation::Spline,
      bool monotonic = true,
      CubicInterpolation::BoundaryCondition leftCondition
          = CubicInterpolation::SecondDerivative,
      doubleOrNull leftConditionValue = 0.0,
      CubicInterpolation::BoundaryCondition rightCondition
          = CubicInterpolation::SecondDerivative,
      doubleOrNull rightConditionValue = 0.0);
};
%enddef

make_mixed_linear_cubic(MixedLinearCubic);
make_mixed_linear_cubic(LogMixedLinearCubic);

%{
using QuantLib::RichardsonExtrapolation;
%}

class RichardsonExtrapolation {
  public:
    Real operator()(Real t=2.0) const;
    Real operator()(Real t, Real s) const;
    
#if defined(SWIGPYTHON)
    %extend {
        RichardsonExtrapolation(
            PyObject* fct, Real delta_h, Real n = Null<Real>()) {
        
            UnaryFunction f(fct);
            return new RichardsonExtrapolation(f, delta_h, n); 
        }
    }
#elif defined(SWIGJAVA) || defined(SWIGCSHARP)
    %extend {
        RichardsonExtrapolation(
            UnaryFunctionDelegate* fct, Real delta_h, Real n = Null<Real>()) {
        
            UnaryFunction f(fct);
            return new RichardsonExtrapolation(f, delta_h, n); 
        }
    }
#else
  private:
    RichardsonExtrapolation();
#endif
};


%{
class SafeConvexMonotoneInterpolation {
  public:
    SafeConvexMonotoneInterpolation(const Array& x, const Array& y,
                                    Real quadraticity = 0.3,
                                    Real monotonicity = 0.7,
                                    bool forcePositive = true)
    : x_(x), y_(y), f_(x_.begin(), x_.end(), y_.begin(),
                       quadraticity, monotonicity, forcePositive) {}
    Real operator()(Real x, bool allowExtrapolation=false) {
        return f_(x, allowExtrapolation);
    }
    Array x_, y_;
    QuantLib::ConvexMonotoneInterpolation<Array::const_iterator, Array::const_iterator> f_;
};
%}


%{
using QuantLib::ChebyshevInterpolation;
%}

class ChebyshevInterpolation {
    #if defined(SWIGCSHARP)
    %rename(call) operator();
    #endif

  public:
    enum PointsType {FirstKind, SecondKind};
    ChebyshevInterpolation(const Array& f, PointsType pointsType = SecondKind);
#if defined(SWIGPYTHON)
    %extend {
        ChebyshevInterpolation(
            Size n, PyObject* fct, PointsType pointsType = SecondKind) {
        
            UnaryFunction f(fct);
            return new ChebyshevInterpolation(n, f, pointsType); 
        }
    }
#elif defined(SWIGJAVA) || defined(SWIGCSHARP)
    %extend {
        ChebyshevInterpolation(
            Size n, UnaryFunctionDelegate* fct, PointsType pointsType = SecondKind) {
        
            UnaryFunction f(fct);
            return new ChebyshevInterpolation(n, f, pointsType); 
        }
    }
#endif
    
    Real operator()(Real z, bool allowExtrapolation=false) const;
    static Array nodes(Size n, PointsType pointsType);
};


%rename(ConvexMonotoneInterpolation) SafeConvexMonotoneInterpolation;
class SafeConvexMonotoneInterpolation {
    #if defined(SWIGCSHARP)
    %rename(call) operator();
    #endif
  public:
    SafeConvexMonotoneInterpolation(const Array& x, const Array& y,
                                    Real quadraticity = 0.3,
                                    Real monotonicity = 0.7,
                                    bool forcePositive = true);
    Real operator()(Real x, bool allowExtrapolation=false);
};


#endif
