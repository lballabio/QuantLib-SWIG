
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

#ifndef quantlib_observer_i
#define quantlib_observer_i

%include common.i
%include boost_shared_ptr.i

%{
using QuantLib::Observer;
using QuantLib::Observable;
%}

%shared_ptr(Observable);
class Observable {};


%extend Handle {
    ext::shared_ptr<Observable> asObservable() {
        return ext::shared_ptr<Observable>(*self);
    }
}


#if defined(SWIGPYTHON)

%typemap(in) const ObservableOrHandle& (PyObject* owned) %{
    owned = NULL;
    if (PyObject_HasAttrString($input, "asObservable")) {
        $input = owned = PyObject_CallMethod($input, "asObservable", NULL);
        if (!owned) SWIG_fail;
    }
    $typemap(in, const ext::shared_ptr<Observable>&)
%}
%typemap(freearg) const ObservableOrHandle& %{
    Py_XDECREF(owned$argnum);
%}

%{
typedef ext::shared_ptr<Observable> ObservableOrHandle;

// C++ wrapper for Python observer
class PyObserver : public Observer {
  public:
    PyObserver(PyObject* callback)
    : callback_(callback) {
        /* make sure the Python object stays alive
           as long as we need it */
        Py_XINCREF(callback_);
    }
    PyObserver(const PyObserver& o)
    : callback_(o.callback_) {
        /* make sure the Python object stays alive
           as long as we need it */
        Py_XINCREF(callback_);
    }
    PyObserver& operator=(const PyObserver& o) {
        if ((this != &o) && (callback_ != o.callback_)) {
            Py_XDECREF(callback_);
            callback_ = o.callback_;
            Py_XINCREF(callback_);
        }
        return *this;
    }
    ~PyObserver() {
        // now it can go as far as we are concerned
        Py_XDECREF(callback_);
    }
    void update() {
        PyObject* pyResult = PyObject_CallFunction(callback_,NULL);
        QL_ENSURE(pyResult != NULL, "failed to notify Python observer");
        Py_XDECREF(pyResult);
    }
  private:
    PyObject* callback_;
};
%}

// Python wrapper
%rename(Observer) PyObserver;
class PyObserver {
  public:
    PyObserver(PyObject* callback);
    void registerWith(const ObservableOrHandle&);
    void unregisterWith(const ObservableOrHandle&);
};

#endif


#endif
