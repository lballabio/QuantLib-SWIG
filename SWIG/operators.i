
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl

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

#ifndef quantlib_operators_i
#define quantlib_operators_i

%include common.i
%include linearalgebra.i

%{
typedef QuantLib::BoundaryCondition<QuantLib::TridiagonalOperator>
        DefaultBoundaryCondition;
%}

%shared_ptr(DefaultBoundaryCondition)
class DefaultBoundaryCondition {
    %rename(NoSide) None;
  private:
    DefaultBoundaryCondition();
  public:
    enum Side { None, Upper, Lower };
};

%{
using QuantLib::NeumannBC;
using QuantLib::DirichletBC;
%}

%shared_ptr(NeumannBC)
class NeumannBC : public DefaultBoundaryCondition {
  public:
    NeumannBC(Real value, DefaultBoundaryCondition::Side side);
};

%shared_ptr(DirichletBC)
class DirichletBC : public DefaultBoundaryCondition {
  public:
    DirichletBC(Real value, DefaultBoundaryCondition::Side side);
};



%{
using QuantLib::TridiagonalOperator;
%}

class TridiagonalOperator {
  public:
    // constructors
    TridiagonalOperator(const Array& low, const Array& mid, const Array& high);
    // operator interface
    Array solveFor(const Array& rhs) const;
    Array applyTo(const Array& v) const;
    // inspectors
    Size size() const;
    // modifiers
    void setFirstRow(Real, Real);
    void setMidRow(Size, Real, Real, Real);
    void setMidRows(Real, Real, Real);
    void setLastRow(Real, Real);
    // identity
    static TridiagonalOperator identity(Size size);
    #if defined(SWIGPYTHON)
    %extend {
        TridiagonalOperator __add__(const TridiagonalOperator& O) {
            return *self+O;
        }
        TridiagonalOperator __sub__(const TridiagonalOperator& O) {
            return *self-O;
        }
        TridiagonalOperator __mul__(Real a) {
            return *self*a;
        }
        TridiagonalOperator __div__(Real a) {
            return *self/a;
        }
        #if defined(SWIGPYTHON)
        TridiagonalOperator __iadd__(const TridiagonalOperator& O) {
            return *self+O;
        }
        TridiagonalOperator __isub__(const TridiagonalOperator& O) {
            return *self-O;
        }
        TridiagonalOperator __imul__(Real a) {
            return *self*a;
        }
        TridiagonalOperator __rmul__(Real a) {
            return *self*a;
        }
        TridiagonalOperator __idiv__(Real a) {
            return *self/a;
        }
        #endif
    }
    #endif
};


%{
using QuantLib::DPlus;
using QuantLib::DMinus;
using QuantLib::DZero;
using QuantLib::DPlusDMinus;
%}

class DPlus : public TridiagonalOperator {
  public:
    DPlus(Size gridPoints, Real h);
};
class DMinus : public TridiagonalOperator {
  public:
    DMinus(Size gridPoints, Real h);
};
class DZero : public TridiagonalOperator {
  public:
    DZero(Size gridPoints, Real h);
};
class DPlusDMinus : public TridiagonalOperator {
  public:
    DPlusDMinus(Size gridPoints, Real h);
};


#endif
