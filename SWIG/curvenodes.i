/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2026 Kyrylo Protsenko

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

#ifndef quantlib_curve_nodes_i
#define quantlib_curve_nodes_i

%include common.i
%include types.i
%include date.i
%include vectors.i
%include termstructures.i

%{
#include <ql/termstructures/interpolatedcurve.hpp>
using QuantLib::InterpolatedNodes;
%}

// Node dates and values on the interpolated curve classes
%shared_ptr(InterpolatedNodes)
class InterpolatedNodes {
  private:
    InterpolatedNodes();
  public:
    const std::vector<Date>& dates() const;
    const std::vector<Time>& times() const;
    const std::vector<Real>& data() const;
};

// Declared after %shared_ptr(InterpolatedNodes) so that the C# binding returns an InterpolatedNodes proxy
%inline %{
    ext::shared_ptr<InterpolatedNodes> as_interpolated_nodes(const ext::shared_ptr<TermStructure>& curve) {
        return ext::dynamic_pointer_cast<InterpolatedNodes>(curve);
    }
%}

#endif
