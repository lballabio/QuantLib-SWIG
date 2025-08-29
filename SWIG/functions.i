
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2015 Klaus Spanderen
 
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

#ifndef quantlib_functions_i
#define quantlib_functions_i

%include linearalgebra.i
%include types.i

%{
using QuantLib::CostFunction;
%}

#if defined(SWIGPYTHON)

%{
class UnaryFunction {
  public:
    UnaryFunction(PyObject* function)
    : function_(PyPtr::fromBorrowed(function)) {}

    Real operator()(Real x) const {
        auto pyResult = PyPtr::fromResult(
            PyObject_CallFunction(function_.get(), "d", x),
            "failed to call Python function");
        return PyFloat_AsDouble(pyResult.get());
    }
    Real derivative(Real x) const {
        auto pyResult = PyPtr::fromResult(
            PyObject_CallMethod(function_.get(), "derivative", "d", x),
            "failed to call derivative() on Python object");
        return PyFloat_AsDouble(pyResult.get());
    }
  private:
    PyPtr function_;
};

class BinaryFunction {
  public:
    BinaryFunction(PyObject* function)
    : function_(PyPtr::fromBorrowed(function)) {}

    Real operator()(Real x, Real y) const {
        auto pyResult = PyPtr::fromResult(
            PyObject_CallFunction(function_.get(), "dd", x, y),
            "failed to call Python function");
        return PyFloat_AsDouble(pyResult.get());
    }
  private:
    PyPtr function_;
};

class PyCostFunction : public CostFunction {
  public:
    PyCostFunction(PyObject* function)
    : function_(PyPtr::fromBorrowed(function)) {}

    Real value(const Array& x) const {
        auto tuple = PyPtr::fromResult(PyTuple_New(x.size()), "failed to convert arguments");
        for (Size i=0; i<x.size(); i++)
            PyTuple_SetItem(tuple.get(), i, PyFloat_FromDouble(x[i]));

        auto pyResult = PyPtr::fromResult(
            PyObject_CallFunction(function_.get(), "O", tuple.get()),
            "failed to call Python function");
        return PyFloat_AsDouble(pyResult.get());
    }
    Array values(const Array& x) const {
        QL_FAIL("Not implemented");
        // Should be straight forward to copy from a python list
        // to an array
    }
  private:
    PyPtr function_;
};
%}

#elif defined(SWIGJAVA)

%{
class UnaryFunctionDelegate {
  public:
    virtual ~UnaryFunctionDelegate() {}
    virtual Real value(Real x) const {
        QL_FAIL("implementation of UnaryFunctionDelegate.value is missing");
    }
};

class UnaryFunction {
  public:
    UnaryFunction(UnaryFunctionDelegate* delegate)
    : delegate_(delegate) { }

    virtual ~UnaryFunction() { }

    Real operator()(Real x) const {
        return delegate_->value(x);
    }

  private:
    UnaryFunctionDelegate* delegate_;
};
%}

class UnaryFunction {
  public:
    UnaryFunction(UnaryFunctionDelegate*);
    Real operator()(Real x) const;
};

%feature("director") UnaryFunctionDelegate;

class UnaryFunctionDelegate {
  public:
    virtual ~UnaryFunctionDelegate();
    virtual Real value(Real x) const;
};

%{
class BinaryFunctionDelegate {
  public:
    virtual ~BinaryFunctionDelegate() {}
    virtual Real value(Real x, Real y) const {
    	QL_FAIL("implementation of BinaryFunctionDelegate.value is missing");
    }	
};

class BinaryFunction {
  public:
    BinaryFunction(BinaryFunctionDelegate* delegate)
    : delegate_(delegate) {}
    
    virtual ~BinaryFunction() {}
    
    Real operator()(Real x, Real y) const {
    	return delegate_->value(x, y);
    }
    
  private:
    BinaryFunctionDelegate* delegate_; 
};
%}

class BinaryFunction {
  public:
    BinaryFunction(BinaryFunctionDelegate*);
    Real operator()(Real, Real) const;
};

%feature("director") BinaryFunctionDelegate;

class BinaryFunctionDelegate {
  public:
    virtual ~BinaryFunctionDelegate();
    virtual Real value(Real, Real) const;
};

%{
class CostFunctionDelegate {
  public:
    virtual ~CostFunctionDelegate() {}
    virtual Real value(const Array& x) const {
      QL_FAIL("implementation of CostFunctionDelegate.value is missing");
    }

    virtual Array values(const Array& x) const {
      QL_FAIL("implementation of CostFunctionDelegate.values is missing");
    }
};

class JavaCostFunction : public CostFunction {
  public:
    JavaCostFunction(CostFunctionDelegate* delegate)
    : delegate_(delegate) { }

    virtual ~JavaCostFunction(){ }

    virtual Real value(const Array& x ) const{
      return delegate_->value(x);
    }

    virtual Array values(const Array& x) const {
      Array retVal = delegate_->values(x);
      return retVal;
    }

  private:
    CostFunctionDelegate* delegate_;
};
%}

class JavaCostFunction {
  public:
    JavaCostFunction(CostFunctionDelegate* delegate);

    virtual ~JavaCostFunction();
    virtual Real value(const Array& x ) const;
    virtual Array values(const Array& x) const;

  private:
    CostFunctionDelegate* delegate_;
};

%feature("director") CostFunctionDelegate;

class CostFunctionDelegate {
  public:
    virtual ~CostFunctionDelegate();

    virtual Real value(const Array& x) const;
    virtual Array values(const Array& x) const;
};

#elif defined(SWIGCSHARP)

%rename(call) operator();
%{
class UnaryFunctionDelegate {
  public:
    virtual ~UnaryFunctionDelegate() {}
    virtual Real value(Real x) const {
        QL_FAIL("implementation of UnaryFunctionDelegate.value is missing");
    };
};

class UnaryFunction {
  public:
    UnaryFunction(UnaryFunctionDelegate* delegate)
    : delegate_(delegate) { }

    virtual ~UnaryFunction() { }

    Real operator()(Real x) const {
        return delegate_->value(x);
    }

  private:
    UnaryFunctionDelegate* delegate_;
};
%}

class UnaryFunction {
  public:
    UnaryFunction(UnaryFunctionDelegate*);
    Real operator()(Real x) const;
};

%feature("director") UnaryFunctionDelegate;

class UnaryFunctionDelegate {
  public:
    virtual ~UnaryFunctionDelegate();
    virtual Real value(Real x) const;
};

%{
class BinaryFunctionDelegate {
  public:
    virtual ~BinaryFunctionDelegate() {}
    virtual Real value(Real x, Real y) const {
    	QL_FAIL("implementation of BinaryFunctionDelegate.value is missing");
    }	
};

class BinaryFunction {
  public:
    BinaryFunction(BinaryFunctionDelegate* delegate)
    : delegate_(delegate) {}
    
    virtual ~BinaryFunction() {}
    
    Real operator()(Real x, Real y) const {
    	return delegate_->value(x, y);
    }
    
  private:
    BinaryFunctionDelegate* delegate_; 
};
%}

class BinaryFunction {
  public:
    BinaryFunction(BinaryFunctionDelegate*);
    Real operator()(Real, Real) const;
};

%feature("director") BinaryFunctionDelegate;

class BinaryFunctionDelegate {
  public:
    virtual ~BinaryFunctionDelegate();
    virtual Real value(Real, Real) const;
};

%{
class CostFunctionDelegate {
  public:
    virtual ~CostFunctionDelegate() {}
    virtual Real value(const Array& x) const {
      QL_FAIL("implementation of CostFunctionDelegate.value is missing");
    }

    virtual Array values(const Array& x) const {
      QL_FAIL("implementation of CostFunctionDelegate.values is missing");
    }
};

class DotNetCostFunction : public CostFunction {
  public:
    DotNetCostFunction(CostFunctionDelegate* delegate)
    : delegate_(delegate) { }

    virtual ~DotNetCostFunction(){ }

    virtual Real value(const Array& x ) const{
      return delegate_->value(x);
    }

    virtual Array values(const Array& x) const {
      Array retVal = delegate_->values(x);
      return retVal;
    }

  private:
    CostFunctionDelegate* delegate_;
};
%}

class DotNetCostFunction {
  public:
    DotNetCostFunction(CostFunctionDelegate* delegate);

    virtual ~DotNetCostFunction();
    virtual Real value(const Array& x ) const;
    virtual Array values(const Array& x) const;

  private:
    CostFunctionDelegate* delegate_;
};

%feature("director") CostFunctionDelegate;

class CostFunctionDelegate {
  public:
    virtual ~CostFunctionDelegate();

    virtual Real value(const Array& x) const;
    virtual Array values(const Array& x) const;
};

#endif
#endif
