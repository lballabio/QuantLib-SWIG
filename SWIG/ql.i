
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008 StatPro Italia srl

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

// Undefine symbols that are also used in quantlib

#if defined(SWIGPYTHON)
%{
#ifdef barrier
#undef barrier
#endif
%}
#endif

%{
#include <ql/quantlib.hpp>

#if QL_HEX_VERSION < 0x01410000
    #error at least QuantLib 1.41 required, please update
#endif

#if defined (SWIGJAVA) || defined (SWIGCSHARP) 
  #ifndef QL_ENABLE_THREAD_SAFE_OBSERVER_PATTERN
    #ifdef BOOST_MSVC
      #pragma message(\
          "Quantlib has not been compiled with the thread-safe "           \
          "observer pattern being enabled. This can lead to spurious "     \
          "crashes or pure virtual function calls within the JVM or .NET "  \
          "ecosystem due to the async garbage collector. Please consider " \
          "enabling QL_ENABLE_THREAD_SAFE_OBSERVER_PATTERN "               \
          "in ql/userconfig.hpp.")
    #else
      #warning \
Quantlib has not been compiled with the thread-safe \
observer pattern being enabled. This can lead to spurious \
crashes or pure virtual function calls within the JVM or .NET \
ecosystem due to the async garbage collector. Please consider \
passing --enable-thread-safe-observer-pattern when using the \
GNU autoconf configure script.
    #endif
  #endif
#endif


// add here SWIG version check

#if defined(_MSC_VER)         // Microsoft Visual C++ 6.0
// disable Swig-dependent warnings

// 'identifier1' has C-linkage specified,
// but returns UDT 'identifier2' which is incompatible with C
#pragma warning(disable: 4190)

// 'int' : forcing value to bool 'true' or 'false' (performance warning)
#pragma warning(disable: 4800)

// debug info too long etc etc
#pragma warning(disable: 4786)
#endif
%}

#ifdef SWIGJAVA
%include "enumtypesafe.swg"
#endif

// common name mappings
#if defined(SWIGJAVA)
%rename(add)           operator+;
%rename(add)           __add__;
%rename(subtract)      operator-;
%rename(subtract)      __sub__;
%rename(multiply)      operator*;
%rename(multiply)      __mul__;
%rename(divide)        operator/;
%rename(divide)        __div__;
%rename(getValue)      operator();
%rename(equals)        operator==;
%rename(unEquals)      operator!=;
%rename(compareTo)     __cmp__;
%javamethodmodifiers   __cmp__ "@Override public"
%rename(hashCode)      __hash__;
%rename(toString)      __str__;
%rename(repr)          __repr__;
#elif defined(SWIGCSHARP)
%rename(Add)           operator+;
%rename(Add)           __add__;
%rename(Subtract)      operator-;
%rename(Subtract)      __sub__;
%rename(Multiply)      operator*;
%rename(Multiply)      __mul__;
%rename(Divide)        operator/;
%rename(Divide)        __div__;
#endif

%{
// we do not want to see the deprecated warnings here
QL_DEPRECATED_DISABLE_WARNING
%}

%include common.i
%include vectors.i
%include tuple.i
%include asianoptions.i
%include barrieroptions.i
%include basketoptions.i
%include blackformula.i
%include bonds.i
%include bondfunctions.i
%include calendars.i
%include calibratedmodel.i
%include calibrationhelpers.i
%include capfloor.i
%include cashflows.i
%include cliquetoptions.i
%include convertiblebonds.i
%include credit.i
%include creditdefaultswap.i
%include currencies.i
%include date.i
%include daycounters.i
%include defaultprobability.i
%include discountcurve.i
%include distributions.i
%include dividends.i
%include exchangerates.i
%include exercise.i
%include fdm.i
%include fittedbondcurve.i
%include forward.i
%include forwardcurve.i
%include fra.i
%include functions.i
%include futures.i
%include gaussian1dmodel.i
%include grid.i
%include indexes.i
%include inflation.i
%include instruments.i
%include integrals.i
%include interestrate.i
%include interpolation.i
%include lazyobject.i
%include linearalgebra.i
%include localvolatilities.i
%include lmm.i
%include lookbackoptions.i
%include marketelements.i
%include money.i
%include montecarlo.i
%include null.i
%include observer.i
%include ode.i
%include operators.i
%include optimizers.i
%include parameter.i
%include options.i
%include payoffs.i
%include piecewiseyieldcurve.i
%include randomnumbers.i
%include ratehelpers.i
%include rounding.i
%include scheduler.i
%include settings.i
%include shortratemodels.i
%include slv.i
%include spreaddiscountcurve.i
%include spreadoption.i
%include statistics.i
%include stochasticprocess.i
%include swap.i
%include swaption.i
%include swingoption.i
%include termstructures.i
%include timebasket.i
%include timeseries.i
%include tracing.i
%include types.i
%include volatilities.i
%include volatilitymodels.i
%include zerocurve.i
%include old_volatility.i
