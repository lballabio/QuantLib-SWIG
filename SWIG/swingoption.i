/*
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
 Copyright (C) 2018 Klaus Spanderen
 Copyright (C) 2018, 2019 Matthias Lungwitz

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

#ifndef quantlib_swing_options_i
#define quantlib_swing_options_i

%include options.i

%{
using QuantLib::VanillaSwingOption;
%}

%shared_ptr(VanillaSwingOption)
class VanillaSwingOption : public OneAssetOption {
  public:
    VanillaSwingOption(
        const ext::shared_ptr<Payoff>& payoff,
        const ext::shared_ptr<SwingExercise>& ex,
        Size minExerciseRights, Size maxExerciseRights);
};

%{
using QuantLib::FdSimpleBSSwingEngine;
using QuantLib::FdSimpleExtOUJumpSwingEngine;
%}

%shared_ptr(FdSimpleBSSwingEngine)
class FdSimpleBSSwingEngine : public PricingEngine {
  public:
    FdSimpleBSSwingEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
            Size tGrid = 50, Size xGrid = 100,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas());
};

%shared_ptr(FdSimpleExtOUJumpSwingEngine)
class FdSimpleExtOUJumpSwingEngine : public PricingEngine {
  public:
    %extend {
        FdSimpleExtOUJumpSwingEngine(
            const ext::shared_ptr<ExtOUWithJumpsProcess>& process,
            const ext::shared_ptr<YieldTermStructure>& rTS,
            Size tGrid = 50, Size xGrid = 200, Size yGrid=50,
            const std::vector<std::pair<Time,Real> >& shape =
                                         std::vector<std::pair<Time,Real> >(),
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer()) {

            ext::shared_ptr<FdSimpleExtOUJumpSwingEngine::Shape> curve(
                              new FdSimpleExtOUJumpSwingEngine::Shape(shape));

            return new FdSimpleExtOUJumpSwingEngine(
                    process, rTS, tGrid, xGrid, yGrid,
                    curve, schemeDesc);
        }
    }
};


#endif
