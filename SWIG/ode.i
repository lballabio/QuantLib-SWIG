/*
 Copyright (C) 2019 Klaus Spanderen

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

#ifndef quantlib_ode_i
#define quantlib_ode_i

%include common.i
%include functions.i

#if defined(SWIGJAVA) || defined(SWIGCSHARP)

%{
class OdeFctDelegate {
  public:
    virtual ~OdeFctDelegate() {}
    virtual std::vector<Real> value(
        Real x, const std::vector<Real>& y) const {
                
        QL_FAIL("implementation of OdeFctDelegate.value is missing");
    }
};

class OdeFct {
  public:
    OdeFct(OdeFctDelegate* delegate)
    : delegate_(delegate) { }

    virtual ~OdeFct() { }

    const std::vector<Real> operator()(Real x, const std::vector<Real>& y) const {
        std::vector<Real> retVal = delegate_->value(x, y);
        return retVal;
    }

  private:
    OdeFctDelegate* delegate_;
};
%}

%feature("director") OdeFctDelegate;

class OdeFctDelegate {
  public:
    virtual ~OdeFctDelegate();
    virtual std::vector<Real> value(Real x, const std::vector<Real>& y) const;
};

#elif defined(SWIGPYTHON)

%{
class OdeFct {
  public:
    OdeFct(PyObject* function)
    : function_(PyPtr::fromBorrowed(function)) {}

    const std::vector<Real> operator()(Real x, const std::vector<Real>& y) const {
        auto pyY = PyPtr::fromResult(PyList_New(y.size()), "failed to convert arguments");
        for (Size i=0; i < y.size(); ++i)
            PyList_SetItem(pyY.get(), i, PyFloat_FromDouble(y[i]));

        auto pyResult = PyPtr::fromResult(
            PyObject_CallFunction(function_.get(), "dO", x, pyY.get()),
            "failed to call Python function");

        QL_ENSURE(
            PyList_Check(pyResult.get()) && PyList_Size(pyResult.get()) == y.size(),
            "Python function did not return a list of correct size");

        std::vector<Real> retVal(y.size());
        for (Size i=0; i < y.size(); ++i)
            retVal[i] = PyFloat_AsDouble(PyList_GetItem(pyResult.get(), i));

        return retVal;
    }
  private:
    PyPtr function_;
};
%}

#endif


%{
using QuantLib::AdaptiveRungeKutta;
%}

template <class T = Real>
class AdaptiveRungeKutta {
  public:
    AdaptiveRungeKutta(const Real eps=1.0e-6,
                       const Real h1=1.0e-4,
                       const Real hmin=0.0);
                       
    %extend {
    
      #if defined(SWIGPYTHON)
      
        T operator()(PyObject* fct, T y1, Real x1, Real x2) {
            BinaryFunction f(fct);
            return self->operator()(f, y1, x1, x2);
        }
        
        std::vector<T> operator()(
            PyObject* fct, const std::vector<T>& y1, Real x1, Real x2) {
            OdeFct f(fct);
            return self->operator()(f, y1, x1, x2);
        }

      #elif defined(SWIGJAVA) || defined(SWIGCSHARP)
      
        T operator()(BinaryFunctionDelegate* fct, T y1, Real x1, Real x2) {
            BinaryFunction f(fct);
            return self->operator()(f, y1, x1, x2);         
        }

        std::vector<T> operator()(
            OdeFctDelegate* fct, const std::vector<T>& y1, Real x1, Real x2) {            
            OdeFct f(fct);
            return self->operator()(f, y1, x1, x2);
        }
        
      #endif     
    }                     
};

%template(RungeKutta) AdaptiveRungeKutta<Real>;

#endif
