
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005 StatPro Italia srl

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

#ifndef quantlib_null_values_i
#define quantlib_null_values_i

%include exception.i

%{
using QuantLib::Null;
typedef int intOrNull;
typedef double doubleOrNull;
%}

%inline%{
int nullInt() { return Null<int>(); }
double nullDouble() { return Null<double>(); }
%}

#if defined(SWIGPYTHON)

%typemap(in) intOrNull {
    if ($input == Py_None)
        $1 = Null<int>();
    else if (PyInt_Check($input))
        $1 = int(PyInt_AsLong($input));
    else
        SWIG_exception(SWIG_TypeError,"int expected");
}
%typecheck(SWIG_TYPECHECK_INTEGER) intOrNull {
    $1 = ($input == Py_None || PyInt_Check($input)) ? 1 : 0;
}
%typemap(out) intOrNull {
    if ($1 == Null<int>()) {
        Py_INCREF(Py_None);
        $result = Py_None;
    } else {
        $result = PyInt_FromLong(long($1));
    }
}

%typemap(in) doubleOrNull {
    if ($input == Py_None)
        $1 = Null<double>();
    else if (PyFloat_Check($input))
        $1 = PyFloat_AsDouble($input);
    else
        SWIG_exception(SWIG_TypeError,"double expected");
}
%typecheck(SWIG_TYPECHECK_DOUBLE) doubleOrNull {
    $1 = ($input == Py_None || PyFloat_Check($input)) ? 1 : 0;
}
%typemap(out) doubleOrNull {
    if ($1 == Null<double>()) {
        Py_INCREF(Py_None);
        $result = Py_None;
    } else {
        $result = PyFloat_FromDouble($1);
    }
}

#elif defined(SWIGJAVA)

typedef int intOrNull;
typedef double doubleOrNull;

#elif defined(SWIGCSHARP)

typedef int intOrNull;
typedef double doubleOrNull;

#elif defined(SWIGR)

%typemap(rtype) intOrNull "numeric";
%typemap(scoercein) intOrNull
   %{ $input = as($input, "integer"); %}
%typemap(scoerceout) intOrNull %{ %}

%typemap(in) intOrNull {
  $1 = ($1_ltype) INTEGER($input)[0];
  if ($1 == R_NaInt) $1 = Null<int>();
}

%typemap(out) intOrNull {
    if ($1 == Null<int>())
    $result = Rf_ScalarLogical(NA_LOGICAL)
    else
    $result = ScalarInteger($1);
}

%typemap(rtype) doubleOrNull "numeric";
%typemap(scoercein) doubleOrNull
    %{ $input = as($input, "numeric"); %}
%typemap(scoerceout) doubleOrNull %{ %}

%typemap(in) doubleOrNull {
  $1 = ($1_ltype) REAL($input)[0];
  if (R_IsNA($1)) $1 = Null<double>();
}

%typemap(out) doubleOrNull {
    if ($1 == Null<double>())
        $result = Rf_ScalarLogical(NA_LOGICAL);
    else
        $result = Rf_ScalarReal($1);
}

#endif


#endif
