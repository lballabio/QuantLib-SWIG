/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2015 Klaus Spanderen

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

#ifndef quantlib_optimizers_i
#define quantlib_optimizers_i

%include functions.i
%include linearalgebra.i

// 1D Solvers

%{
using QuantLib::Bisection;
using QuantLib::Brent;
using QuantLib::FalsePosition;
using QuantLib::Newton;
using QuantLib::NewtonSafe;
using QuantLib::Ridder;
using QuantLib::Secant;
%}

#if defined(SWIGMZSCHEME)
%typecheck(SWIG_TYPECHECK_POINTER) Scheme_Object* {
    $1 = 1;
}
#elif defined(SWIGGUILE)
%typecheck(SWIG_TYPECHECK_POINTER) SCM {
    $1 = 1;
}
#endif

%define DeclareSolver(SolverName)
class SolverName {
    #if defined(SWIGRUBY)
    %rename("maxEvaluations=")      setMaxEvaluations;
    %rename("lowerBound=")          setLowerBound;
    %rename("upperBound=")          setUpperBound;
    #elif defined(SWIGMZSCHEME) || defined(SWIGGUILE)
    %rename("max-evaluations-set!") setMaxEvaluations;
    %rename("lower-bound-set!")     setLowerBound;
    %rename("upper-bound-set!")     setUpperBound;
    #endif
  public:
    void setMaxEvaluations(Size evaluations);
    void setLowerBound(Real lowerBound);
    void setUpperBound(Real upperBound);
    %extend {
        #if defined(SWIGPYTHON)
        Real solve(PyObject* function, Real xAccuracy,
                   Real guess, Real step) {
            UnaryFunction f(function);
            return self->solve(f, xAccuracy, guess, step);
        }
        Real solve(PyObject* function, Real xAccuracy,
                   Real guess, Real xMin, Real xMax) {
            UnaryFunction f(function);
            return self->solve(f, xAccuracy, guess, xMin, xMax);
        }
        #elif defined(SWIGRUBY)
        Real solve(Real xAccuracy, Real guess, Real step) {
            UnaryFunction f;
            return self->solve(f, xAccuracy, guess, step);
        }
        Real solve(Real xAccuracy, Real guess,
                   Real xMin, Real xMax) {
            UnaryFunction f;
            return self->solve(f, xAccuracy, guess, xMin, xMax);
        }
        #elif defined(SWIGMZSCHEME)
        Real solve(Scheme_Object* function, Real xAccuracy,
                   Real guess, Real step) {
            UnaryFunction f(function);
            return self->solve(f, xAccuracy, guess, step);
        }
        Real solve(Scheme_Object* function, Real xAccuracy,
                   Real guess, Real xMin, Real xMax) {
            UnaryFunction f(function);
            return self->solve(f, xAccuracy, guess, xMin, xMax);
        }
        #elif defined(SWIGGUILE)
        Real solve(SCM function, Real xAccuracy,
                   Real guess, Real step) {
            UnaryFunction f(function);
            return self->solve(f, xAccuracy, guess, step);
        }
        Real solve(SCM function, Real xAccuracy,
                   Real guess, Real xMin, Real xMax) {
            UnaryFunction f(function);
            return self->solve(f, xAccuracy, guess, xMin, xMax);
        }
        #elif defined(SWIGJAVA) || defined(SWIGCSHARP)
        Real solve(UnaryFunctionDelegate* function, Real xAccuracy,
                   Real guess, Real step) {
            UnaryFunction f(function);
            return self->solve(f, xAccuracy, guess, step);
        }
        Real solve(UnaryFunctionDelegate* function, Real xAccuracy,
                   Real guess, Real xMin, Real xMax) {
            UnaryFunction f(function);
            return self->solve(f, xAccuracy, guess, xMin, xMax);
        }
        #endif
    }
};
%enddef

// Keep this list in sync with bondfunctions.i yield solvers.
// Actual solvers
DeclareSolver(Brent);
DeclareSolver(Bisection);
DeclareSolver(FalsePosition);
DeclareSolver(Ridder);
DeclareSolver(Secant);

#if defined(SWIGPYTHON)
// these two need f.derivative()
DeclareSolver(Newton);
DeclareSolver(NewtonSafe);
#endif


// Optimizers

%{
using QuantLib::Constraint;
using QuantLib::BoundaryConstraint;
using QuantLib::NoConstraint;
using QuantLib::PositiveConstraint;
using QuantLib::CompositeConstraint;
using QuantLib::NonhomogeneousBoundaryConstraint;
%}

class Constraint {
    // prevent direct instantiation
  private:
    Constraint();
};

class BoundaryConstraint : public Constraint {
  public:
    BoundaryConstraint(Real lower, Real upper);
};

class NoConstraint : public Constraint {
  public:
    NoConstraint();
};

class PositiveConstraint : public Constraint {
  public:
    PositiveConstraint();
};

class CompositeConstraint : public Constraint {
  public:
    CompositeConstraint(const Constraint& c1, const Constraint& c2);
};

class NonhomogeneousBoundaryConstraint : public Constraint {
  public:
    NonhomogeneousBoundaryConstraint(const Array& l, const Array& u);
};

%{
using QuantLib::EndCriteria;
%}

class EndCriteria {
    #if defined(SWIGRUBY)
    %rename("setPositiveOptimization!") setPositiveOptimization;
    #elif defined(SWIGMZSCHEME) || defined(SWIGGUILE)
    %rename(call) operator();
    %rename("positive-optimization-set!") setPositiveOptimization;
    #elif defined(SWIGCSHARP) || defined(SWIGPERL)
    %rename(call) operator();
    #elif defined(SWIGPYTHON)
    %rename(NoCriteria) None;
    #endif
  public:
    enum Type {
        None,
        MaxIterations,
        StationaryPoint,
        StationaryFunctionValue,
        StationaryFunctionAccuracy,
        ZeroGradientNorm,
        Unknown
    };
    EndCriteria(Size maxIteration,
                Size maxStationaryStateIterations,
                Real rootEpsilon,
                Real functionEpsilon,
                Real gradientNormEpsilon);
    bool operator()(Size iteration,
                    Size &statState,
                    const bool positiveOptimization,
                    const Real fold,
                    const Real normgold,
                    const Real fnew,
                    const Real normgnewx,
                    EndCriteria::Type & ecType) const;
};


%{
using QuantLib::OptimizationMethod;
using QuantLib::ConjugateGradient;
using QuantLib::Simplex;
using QuantLib::SteepestDescent;
using QuantLib::BFGS;
using QuantLib::LevenbergMarquardt;
using QuantLib::DifferentialEvolution;
using QuantLib::SamplerGaussian;
using QuantLib::SamplerLogNormal;
using QuantLib::SamplerMirrorGaussian;
using QuantLib::ProbabilityBoltzmannDownhill;
using QuantLib::TemperatureExponential;
using QuantLib::ReannealingTrivial;
using QuantLib::GaussianSimulatedAnnealing;
using QuantLib::MirrorGaussianSimulatedAnnealing;
using QuantLib::LogNormalSimulatedAnnealing;

%}

class OptimizationMethod {
  private:
    // prevent direct instantiation
    OptimizationMethod();
};

class ConjugateGradient : public OptimizationMethod {
  public:
    ConjugateGradient();
};

class Simplex : public OptimizationMethod {
  public:
    Simplex(Real lambda);
};

class SteepestDescent : public OptimizationMethod {
  public:
    SteepestDescent();
};

class BFGS : public OptimizationMethod {
  public:
    BFGS();
};

class LevenbergMarquardt : public OptimizationMethod {
  public:
    LevenbergMarquardt(Real epsfcn = 1.0e-8,
                       Real xtol = 1.0e-8,
                       Real gtol = 1.0e-8);
};

class DifferentialEvolution : public OptimizationMethod {
  public:
    DifferentialEvolution();
};

class SamplerGaussian{
  public:
    SamplerGaussian(unsigned long seed = 0);
};

class SamplerLogNormal{
  public:
    SamplerLogNormal(unsigned long seed = 0);
};

class SamplerMirrorGaussian{
  public:
    SamplerMirrorGaussian(const Array& lower, const Array& upper, unsigned long seed = 0);
};

class ProbabilityBoltzmannDownhill{
  public:
    ProbabilityBoltzmannDownhill(unsigned long seed = 0);
};

class TemperatureExponential {
  public:
    TemperatureExponential(Real initialTemp, Size dimension, Real power = 0.95);
};

class ReannealingTrivial {
  public:
    ReannealingTrivial();
};

class GaussianSimulatedAnnealing : public OptimizationMethod {
  public:
    enum ResetScheme{
        NoResetScheme,
        ResetToBestPoint,
        ResetToOrigin
    };
    GaussianSimulatedAnnealing(const SamplerGaussian &sampler,
            const ProbabilityBoltzmannDownhill &probability,
            const TemperatureExponential &temperature,
            const ReannealingTrivial &reannealing = ReannealingTrivial(),
            Real startTemperature = 200.0,
            Real endTemperature = 0.01,
            Size reAnnealSteps = 50,
            ResetScheme resetScheme = ResetToBestPoint,
            Size resetSteps = 150);
};

class MirrorGaussianSimulatedAnnealing : public OptimizationMethod {
  public:
    enum ResetScheme{
        NoResetScheme,
        ResetToBestPoint,
        ResetToOrigin
    };
    MirrorGaussianSimulatedAnnealing(const SamplerMirrorGaussian &sampler,
            const ProbabilityBoltzmannDownhill &probability,
            const TemperatureExponential &temperature,
            const ReannealingTrivial &reannealing = ReannealingTrivial(),
            Real startTemperature = 200.0,
            Real endTemperature = 0.01,
            Size reAnnealSteps = 50,
            ResetScheme resetScheme = ResetToBestPoint,
            Size resetSteps = 150);
};

class LogNormalSimulatedAnnealing : public OptimizationMethod {
  public:
   enum ResetScheme{
        NoResetScheme,
        ResetToBestPoint,
        ResetToOrigin
    };
    LogNormalSimulatedAnnealing(const SamplerLogNormal &sampler,
            const ProbabilityBoltzmannDownhill &probability,
            const TemperatureExponential &temperature,
            const ReannealingTrivial &reannealing = ReannealingTrivial(),
            Real startTemperature = 10.0,
            Real endTemperature = 0.01,
            Size reAnnealSteps = 50,
            ResetScheme resetScheme = ResetToBestPoint,
            Size resetSteps = 150);
};

%{
using QuantLib::Problem;
%}

%inline %{
    class Optimizer {};
%}

#if defined(SWIGPYTHON)
%extend Optimizer {
    Array solve(PyObject* function, Constraint& c,
                OptimizationMethod& m, EndCriteria &e,
                Array &iv) {
        PyCostFunction f(function);
        Problem p(f,c,iv);
        m.minimize(p, e);
        return p.currentValue();
    }
}
#elif defined(SWIGRUBY)
%extend Optimizer {
    Array solve(Constraint& c, OptimizationMethod& m,
                EndCriteria &e,
                Array &iv) {
        RubyCostFunction f;
        Problem p(f,c,iv);
        m.minimize(p, e);
        return p.currentValue();
    }
}
#elif defined(SWIGMZSCHEME)
%extend Optimizer {
    Array solve(Scheme_Object* function, Constraint& c,
                OptimizationMethod& m,
                EndCriteria &e,
                Array &iv) {
        MzCostFunction f(function);
        Problem p(f,c,iv);
        m.minimize(p, e);
        return p.currentValue();
    }
}
#elif defined(SWIGGUILE)
%extend Optimizer {
    Array solve(SCM function, Constraint& c, OptimizationMethod& m,
                EndCriteria &e, Array &iv) {
        GuileCostFunction f(function);
        Problem p(f,c,iv);
        m.minimize(p, e);
        return p.currentValue();
    }
}
#elif defined(SWIGJAVA)
%extend Optimizer {
    Array solve(CostFunctionDelegate* function,
                Constraint& c, OptimizationMethod& m,
                EndCriteria &e, Array &iv) {
        JavaCostFunction f(function);
        Problem p(f,c,iv);
        m.minimize(p, e);
        return p.currentValue();
    }
}
#elif defined(SWIGCSHARP)
%extend Optimizer {
    Array solve(CostFunctionDelegate* function,
                Constraint& c, OptimizationMethod& m,
                EndCriteria &e, Array &iv) {
        DotNetCostFunction f(function);
        Problem p(f,c,iv);
        m.minimize(p, e);
        return p.currentValue();
    }
}
#endif
#endif
