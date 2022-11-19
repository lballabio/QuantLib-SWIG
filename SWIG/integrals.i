
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2013, 2022 Klaus Spanderen

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

#ifndef quantlib_segment_integral_i
#define quantlib_segment_integral_i

%include common.i
%include types.i
%include functions.i

%{
using QuantLib::SegmentIntegral;
using QuantLib::TrapezoidIntegral;
using QuantLib::Default;
using QuantLib::MidPoint;
using QuantLib::SimpsonIntegral;
using QuantLib::GaussKronrodAdaptive;
using QuantLib::GaussKronrodNonAdaptive;
using QuantLib::GaussLobattoIntegral;
using QuantLib::GaussLaguerreIntegration;
using QuantLib::GaussHermiteIntegration;
using QuantLib::GaussJacobiIntegration;
using QuantLib::GaussHyperbolicIntegration;
using QuantLib::GaussLegendreIntegration;
using QuantLib::GaussChebyshevIntegration;
using QuantLib::GaussChebyshev2ndIntegration;
using QuantLib::GaussGegenbauerIntegration;
using QuantLib::TanhSinhIntegral;
%}

%define INTEGRATION_METHODS
    %extend {
        #if defined(SWIGPYTHON)
        Real __call__(PyObject* pyFunction, Real a, Real b) {
            UnaryFunction f(pyFunction);
            return (*self)(f, a, b);
        }
        #elif defined(SWIGJAVA) || defined(SWIGCSHARP)
        Real calculate(UnaryFunctionDelegate* f, Real a, Real b) {
            return (*self)(UnaryFunction(f), a, b);		
        }
        #endif
    }
%enddef

class SegmentIntegral {
  public:
    SegmentIntegral(Size intervals);
    INTEGRATION_METHODS;
};


template <class IntegrationPolicy>
class TrapezoidIntegral {
  public:
    TrapezoidIntegral(Real accuracy, Size maxIterations);
    INTEGRATION_METHODS;
};

%template(TrapezoidIntegralDefault) TrapezoidIntegral<Default>;
%template(TrapezoidIntegralMidPoint) TrapezoidIntegral<MidPoint>;

class SimpsonIntegral {
  public:
    SimpsonIntegral(Real accuracy, Size maxIterations);
    INTEGRATION_METHODS;
};


class GaussKronrodAdaptive {
  public:
    GaussKronrodAdaptive(Real tolerance,
                         Size maxFunctionEvaluations = Null<Size>());
    INTEGRATION_METHODS;
};

class GaussKronrodNonAdaptive {
  public:
    GaussKronrodNonAdaptive(Real absoluteAccuracy,
                            Size maxEvaluations,
                            Real relativeAccuracy);
    INTEGRATION_METHODS;
};

class GaussLobattoIntegral {
  public:
    GaussLobattoIntegral(Size maxIterations,
                         Real absAccuracy,
                         Real relAccuracy = Null<Real>(),
                         bool useConvergenceEstimate = true);
    INTEGRATION_METHODS;
};

%{
using QuantLib::GaussianQuadrature;
%}

class GaussianQuadrature {
  private:
    GaussianQuadrature();
  public:
    Size order() const;
    %extend {
      Array weights() { 
        return self->weights(); 
      }      
      Array x() { 
        return self->x(); 
      }
      #if defined(SWIGPYTHON)
      Real __call__(PyObject* pyFunction) {
          UnaryFunction f(pyFunction);
          return (*self)(f);
      }
      #elif defined(SWIGJAVA) || defined(SWIGCSHARP)
      Real calculate(UnaryFunctionDelegate* f) {
          return (*self)(UnaryFunction(f));
      }
      #endif
    }
};

class GaussLaguerreIntegration: public GaussianQuadrature {
  public:
    GaussLaguerreIntegration(Size n, Real s = 0.0);
};

class GaussHermiteIntegration: public GaussianQuadrature {
  public:
    GaussHermiteIntegration(Size n, Real mu = 0.0);
};

class GaussJacobiIntegration: public GaussianQuadrature {
  public:
    GaussJacobiIntegration(Size n, Real alpha, Real beta);
};

class GaussHyperbolicIntegration: public GaussianQuadrature {
  public:
    GaussHyperbolicIntegration(Size n);
};

class GaussLegendreIntegration: public GaussianQuadrature {
  public:
    GaussLegendreIntegration(Size n);
};

class GaussChebyshevIntegration: public GaussianQuadrature {
  public:
    GaussChebyshevIntegration(Size n);
};

class GaussChebyshev2ndIntegration: public GaussianQuadrature {
  public:
    GaussChebyshev2ndIntegration(Size n);
};

class GaussGegenbauerIntegration: public GaussianQuadrature {
  public:
    GaussGegenbauerIntegration(Size n, Real lambda);
};

class TanhSinhIntegral {
  public:
    TanhSinhIntegral(
        Real relTolerance = std::sqrt(std::numeric_limits<Real>::epsilon()),
        Size maxRefinements = 15,
        Real minComplement = std::numeric_limits<Real>::min() * 4
    );
    INTEGRATION_METHODS;
};

#endif
