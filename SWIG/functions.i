
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
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

#ifndef quantlib_functions_i
#define quantlib_functions_i

%include linearalgebra.i
%include types.i

%{
using QuantLib::CostFunction;
using QuantLib::Disposable;
%}

#if defined(SWIGPYTHON)

%{
class UnaryFunction {
  public:
    UnaryFunction(PyObject* function) : function_(function) {
        Py_XINCREF(function_);
    }
    UnaryFunction(const UnaryFunction& f) : function_(f.function_) {
        Py_XINCREF(function_);
    }
    UnaryFunction& operator=(const UnaryFunction& f) {
        if ((this != &f) && (function_ != f.function_)) {
            Py_XDECREF(function_);
            function_ = f.function_;
            Py_XINCREF(function_);
        }
        return *this;
    }
    ~UnaryFunction() {
        Py_XDECREF(function_);
    }
    Real operator()(Real x) const {
        PyObject* pyResult = PyObject_CallFunction(function_,"d",x);
        QL_ENSURE(pyResult != NULL, "failed to call Python function");
        Real result = PyFloat_AsDouble(pyResult);
        Py_XDECREF(pyResult);
        return result;
    }
    Real derivative(Real x) const {
        PyObject* pyResult =
            PyObject_CallMethod(function_,"derivative","d",x);
        QL_ENSURE(pyResult != NULL,
                  "failed to call derivative() on Python object");
        Real result = PyFloat_AsDouble(pyResult);
        Py_XDECREF(pyResult);
        return result;
    }
  private:
    PyObject* function_;
};

class BinaryFunction {
  public:
    BinaryFunction(PyObject* function) : function_(function) {
        Py_XINCREF(function_);
    }
    BinaryFunction(const BinaryFunction& f)
    : function_(f.function_) {
        Py_XINCREF(function_);
    }
    BinaryFunction& operator=(const BinaryFunction& f) {
        if ((this != &f) && (function_ != f.function_)) {
            Py_XDECREF(function_);
            function_ = f.function_;
            Py_XINCREF(function_);
        }
        return *this;
    }
    ~BinaryFunction() {
        Py_XDECREF(function_);
    }
    Real operator()(Real x, Real y) const {
        PyObject* pyResult = PyObject_CallFunction(function_,"dd",x,y);
        QL_ENSURE(pyResult != NULL, "failed to call Python function");
        Real result = PyFloat_AsDouble(pyResult);
        Py_XDECREF(pyResult);
        return result;
    }
  private:
    PyObject* function_;
};

class PyCostFunction : public CostFunction {
  public:
    PyCostFunction(PyObject* function) : function_(function) {
        Py_XINCREF(function_);
    }
    PyCostFunction(const PyCostFunction& f)
    : function_(f.function_) {
        Py_XINCREF(function_);
    }
    PyCostFunction& operator=(const PyCostFunction& f) {
        if ((this != &f) && (function_ != f.function_)) {
            Py_XDECREF(function_);
            function_ = f.function_;
            Py_XINCREF(function_);
        }
        return *this;
    }
    ~PyCostFunction() {
        Py_XDECREF(function_);
    }
    Real value(const Array& x) const {
        PyObject* tuple = PyTuple_New(x.size());
        for (Size i=0; i<x.size(); i++)
            PyTuple_SetItem(tuple,i,PyFloat_FromDouble(x[i]));
        PyObject* pyResult = PyObject_CallObject(function_,tuple);
        Py_XDECREF(tuple);
        QL_ENSURE(pyResult != NULL, "failed to call Python function");
        Real result = PyFloat_AsDouble(pyResult);
        Py_XDECREF(pyResult);
        return result;
    }
    Disposable<Array> values(const Array& x) const {
        QL_FAIL("Not implemented");
        // Should be straight forward to copy from a python list
        // to an array
    }
  private:
    PyObject* function_;
};
%}

#elif defined(SWIGRUBY)

%{
class UnaryFunction {
  public:
    Real operator()(Real x) const {
        return NUM2DBL(rb_yield(rb_float_new(x)));
    }
};

class RubyCostFunction : public CostFunction {
  public:
    Real value(const Array& x) const {
        VALUE a = rb_ary_new2(x.size());
        for (Size i=0; i<x.size(); i++)
            rb_ary_store(a,i,rb_float_new(x[i]));
        return NUM2DBL(rb_yield(a));
    }
    Disposable<Array> values(const Array& x) const {
        QL_FAIL("Not implemented");
    }
};
%}

#elif defined(SWIGJAVA)

%{
class UnaryFunctionDelegate {
  public:
    virtual ~UnaryFunctionDelegate() {}
    virtual Real value(Real x) const {
        QL_FAIL("implementation of UnaryFunctionDelegate.value is missing");
    };
};

class UnaryFunction : public std::unary_function<Real, Real> {
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

    virtual Disposable<Array> values(const Array& x) const {
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
    virtual Disposable<Array> values(const Array& x) const;

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

class UnaryFunction : public std::unary_function<Real, Real> {
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

    virtual Disposable<Array> values(const Array& x) const {
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
    virtual Disposable<Array> values(const Array& x) const;

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
