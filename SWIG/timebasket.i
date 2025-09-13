
/*
 Copyright (C) 2000, 2001, 2002 RiskMap srl

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

#ifndef quantlib_timebasket_i
#define quantlib_timebasket_i

%include common.i
%include types.i
%include date.i

%{
using QuantLib::TimeBasket;
%}

class TimeBasket {
    #if defined (SWIGPYTHON)
    %rename(__len__) size;
    #endif
  public:
    TimeBasket();
    TimeBasket(const std::vector<Date>&, const std::vector<Real>&);
    Size size();
    TimeBasket rebin(const std::vector<Date>&) const;
    %extend {
        #if defined(SWIGPYTHON)
        Real __getitem__(const Date& d) {
            return (*self)[d];
        }
        void __setitem__(const Date& d, Real value) {
            (*self)[d] = value;
        }
        PyObject* items() {
            auto itemList = PyPtr::fromResult(PyList_New(self->size()),
                                              "failed to convert arguments");
            TimeBasket::iterator i;
            Size j;
            for (i=self->begin(), j=0; i!=self->end(); ++i, ++j) {
                Date* d = new Date(i->first);
                PyObject* item = PyTuple_New(2);
                PyTuple_SetItem(item,0,
                                SWIG_NewPointerObj((void *) d,
                                                   $descriptor(Date *),1));
                PyTuple_SetItem(item,1,PyFloat_FromDouble(i->second));
                PyList_SetItem(itemList.get(),j,item);
            }
            return itemList.release();
        }
        // Python 2.2 methods
        bool __contains__(const Date& d) {
            return self->hasDate(d);
        }
        PyObject* __iter__() {
            auto keyList = PyPtr::fromResult(PyList_New(self->size()),
                                             "failed to convert arguments");
            TimeBasket::iterator i;
            Size j;
            for (i=self->begin(), j=0; i!=self->end(); ++i, ++j) {
                Date* d = new Date(i->first);
                PyList_SetItem(keyList.get(), j,
                               SWIG_NewPointerObj((void *) d,
                                                  $descriptor(Date *),1));
            }
            return PyObject_GetIter(keyList.get());
        }
        #endif
    }
};


#endif
