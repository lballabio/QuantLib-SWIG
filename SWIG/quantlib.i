
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007 StatPro Italia srl
 Copyright (c) 2005 Dominic Thuillier

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

#if defined(SWIGCSHARP)
%module(directors="1") NQuantLibc
#elif defined(SWIGJAVA)
%module(directors="1") QuantLib
#else
%module QuantLib
#endif

%include exception.i

%exception {
    try {
        $action
    } catch (std::out_of_range& e) {
        SWIG_exception(SWIG_IndexError,const_cast<char*>(e.what()));
    } catch (std::exception& e) {
        SWIG_exception(SWIG_RuntimeError,const_cast<char*>(e.what()));
    } catch (...) {
        SWIG_exception(SWIG_UnknownError,"unknown error");
    }
}

#if defined(SWIGPYTHON)
%{
#include <ql/version.hpp>
const int    __hexversion__ = QL_HEX_VERSION;
const char* __version__    = QL_VERSION;
%}

const int __hexversion__;
%immutable;
const char* __version__;
%mutable;
#endif

#if defined(JAVA_AUTOLOAD)
// Automatically load the shared library for JAVA binding
%pragma(java) jniclasscode=%{
  /// Load the JNI library
  static {
    System.loadLibrary("QuantLibJNI");
  }
%}
#endif


#if defined(JAVA_AUTOCLOSEABLE)
%typemap(javaimports) SWIGTYPE %{
import java.lang.AutoCloseable;
%}
%typemap(javainterfaces) SWIGTYPE "AutoCloseable";
%typemap(javacode) SWIGTYPE %{
  @Override
  public void close() {
   this.delete();
  }
%}
#endif


#if !defined(JAVA_FINALIZER)
%typemap(javafinalize) SWIGTYPE %{%}
#endif

//#if defined(SWIGPYTHON)
//%feature("autodoc");
//#endif

%include ql.i
