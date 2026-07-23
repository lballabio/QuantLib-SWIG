/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 QuantLib contributors

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

#ifndef quantlib_multicurve_i
#define quantlib_multicurve_i

%include piecewiseyieldcurve.i

%{
#include <ql/termstructures/multicurve.hpp>
using QuantLib::MultiCurve;
%}

%shared_ptr(MultiCurve);
class MultiCurve {
  public:
    explicit MultiCurve(Real accuracy);
    explicit MultiCurve(const ext::shared_ptr<OptimizationMethod>& optimizer = nullptr,
                        const ext::shared_ptr<EndCriteria>& endCriteria = nullptr);

    %extend {
        Handle<YieldTermStructure>
        addBootstrappedCurve(RelinkableHandle<YieldTermStructure>& internalHandle,
                             const ext::shared_ptr<YieldTermStructure>& curve) {
            ext::shared_ptr<YieldTermStructure> curveCopy = curve;
            return self->addBootstrappedCurve(internalHandle, std::move(curveCopy));
        }

        Handle<YieldTermStructure>
        addNonBootstrappedCurve(RelinkableHandle<YieldTermStructure>& internalHandle,
                                const ext::shared_ptr<YieldTermStructure>& curve) {
            ext::shared_ptr<YieldTermStructure> curveCopy = curve;
            return self->addNonBootstrappedCurve(internalHandle, std::move(curveCopy));
        }
    }
};

#endif
