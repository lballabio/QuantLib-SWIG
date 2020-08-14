/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2015 Klaus Spanderen
 Copyright (C) 2018, 2019 Matthias Lungwitz

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
%include stl.i

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

%define DeclareSolver(SolverName)
class SolverName {
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
#elif defined(SWIGJAVA) || defined(SWIGCSHARP)
%{
class NFunctAndDer {
  public:
    NFunctAndDer(const UnaryFunction& function,
                 const UnaryFunction& derivative)
    : f_(function), d_(derivative) {}
           
    Real operator()(Real x) const { return f_(x); }
    Real derivative(Real x) const { return d_(x); }
  private:          
    UnaryFunction f_, d_;
};
%}
%ignore NFunctAndDer;

%define DeclareJavaNewtonSolver(SolverName)
class SolverName {
  public:
    void setMaxEvaluations(Size evaluations);
    void setLowerBound(Real lowerBound);
    void setUpperBound(Real upperBound);
    %extend {
        Real solve(UnaryFunctionDelegate* function,
                   UnaryFunctionDelegate* derivative,
                   Real xAccuracy, Real guess, Real step) {
            UnaryFunction f(function), d(derivative);            
            return self->solve(NFunctAndDer(f, d), xAccuracy, guess, step);
        }
        Real solve(UnaryFunctionDelegate* function,
                   UnaryFunctionDelegate* derivative,         
                   Real xAccuracy, Real guess, Real xMin, Real xMax) {
            UnaryFunction f(function), d(derivative);            
            return self->solve(NFunctAndDer(f, d), xAccuracy, guess, xMin, xMax);
        }
    }
};
%enddef

DeclareJavaNewtonSolver(Newton);
DeclareJavaNewtonSolver(NewtonSafe);
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

%shared_ptr(Constraint)
class Constraint {
    // prevent direct instantiation
  private:
    Constraint();
};

%shared_ptr(BoundaryConstraint)
class BoundaryConstraint : public Constraint {
  public:
    BoundaryConstraint(Real lower, Real upper);
};

%shared_ptr(NoConstraint)
class NoConstraint : public Constraint {
  public:
    NoConstraint();
};

%shared_ptr(PositiveConstraint)
class PositiveConstraint : public Constraint {
  public:
    PositiveConstraint();
};

%shared_ptr(CompositeConstraint)
class CompositeConstraint : public Constraint {
  public:
    CompositeConstraint(const Constraint& c1, const Constraint& c2);
};

%shared_ptr(NonhomogeneousBoundaryConstraint)
class NonhomogeneousBoundaryConstraint : public Constraint {
  public:
    NonhomogeneousBoundaryConstraint(const Array& l, const Array& u);
};

%{
using QuantLib::EndCriteria;
%}

%shared_ptr(EndCriteria)
class EndCriteria {
    #if defined(SWIGCSHARP)
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

%shared_ptr(OptimizationMethod)
class OptimizationMethod {
  private:
    // prevent direct instantiation
    OptimizationMethod();
};

%shared_ptr(ConjugateGradient)
class ConjugateGradient : public OptimizationMethod {
  public:
    ConjugateGradient();
};

%shared_ptr(Simplex)
class Simplex : public OptimizationMethod {
  public:
    Simplex(Real lambda);
    #if defined(SWIGPYTHON)
    %rename(getLambda) lambda;
    #endif
    Real lambda();
};

%shared_ptr(SteepestDescent)
class SteepestDescent : public OptimizationMethod {
  public:
    SteepestDescent();
};

%shared_ptr(BFGS)
class BFGS : public OptimizationMethod {
  public:
    BFGS();
};

%shared_ptr(LevenbergMarquardt)
class LevenbergMarquardt : public OptimizationMethod {
  public:
    LevenbergMarquardt(Real epsfcn = 1.0e-8,
                       Real xtol = 1.0e-8,
                       Real gtol = 1.0e-8,
                       bool useCostFunctionsJacobian = false);
};

%shared_ptr(DifferentialEvolution)
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

%shared_ptr(GaussianSimulatedAnnealing)
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

%shared_ptr(MirrorGaussianSimulatedAnnealing)
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

%shared_ptr(LogNormalSimulatedAnnealing)
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
