
/*
 Copyright (C) 2009 Cheng Li

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

#ifndef quantlib_rate_helpers_ext_i
#define quantlib_rate_helpers_ext_i

%include ratehelpers.i
%include swap_ext.i
%include indexes_ext.i

%{

#include <qlext/termstructures/yield/ratehelpers.hpp>

using QuantLib::ShiborSwapRateHelper;

typedef boost::shared_ptr<RateHelper> ShiborSwapRateHelperPtr;

%}

%rename(ShiborSwapRateHelper) ShiborSwapRateHelperPtr;
class ShiborSwapRateHelperPtr : public boost::shared_ptr<RateHelper> {
  public:
    %extend {
        ShiborSwapRateHelperPtr(
                const Handle<Quote>& rate,
                const Period& swapTenor,
                Frequency fixedFreq,
                const ShiborPtr& shiborIndex,
                const Period& fwdStart = 0*Days,
                const Handle<YieldTermStructure>& discountingCurve
                                            = Handle<YieldTermStructure>()) {
            boost::shared_ptr<Shibor> shibor =
                boost::dynamic_pointer_cast<Shibor>(shiborIndex);
            return new ShiborSwapRateHelperPtr(
                new ShiborSwapRateHelper(rate, swapTenor, fixedFreq,
                                   shibor, fwdStart,
                                   discountingCurve));
        }

        ShiborSwapPtr swap() {
            return boost::dynamic_pointer_cast<ShiborSwapRateHelper>(*self)->swap();
        }
    }
};

#endif