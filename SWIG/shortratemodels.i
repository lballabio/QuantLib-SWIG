
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2007, 2009 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2007 Luis Cota
 Copyright (C) 2015 Gouthaman Balaraman
 Copyright (C) 2018 Matthias Lungwitz 

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

#ifndef quantlib_short_rate_models_i
#define quantlib_short_rate_models_i

%include calibrationhelpers.i
%include grid.i
%include options.i

// the base class for models
%{
using QuantLib::ShortRateModel;
%}

%shared_ptr(ShortRateModel)
class ShortRateModel : public CalibratedModel {
  private:
    ShortRateModel();
};

%template(ShortRateModelHandle) Handle<ShortRateModel>;
%template(RelinkableShortRateModelHandle)
RelinkableHandle<ShortRateModel>;

// actual models

%{
using QuantLib::OneFactorAffineModel;
using QuantLib::Vasicek;
using QuantLib::HullWhite;
using QuantLib::BlackKarasinski;
using QuantLib::G2;
%}

%shared_ptr(OneFactorAffineModel)
class OneFactorAffineModel : public ShortRateModel {
  private:
    OneFactorAffineModel();
  public:
    virtual Real discountBond(Time now,
                          Time maturity,
                          Array factors) const;

    Real discountBond(Time now, Time maturity, Rate rate) const;
    DiscountFactor discount(Time t) const;  
};

%shared_ptr(Vasicek)
class Vasicek : public OneFactorAffineModel {
  public:
    Vasicek(Rate r0 = 0.05,
               Real a = 0.1,
               Real b = 0.05,
               Real sigma = 0.01,
               Real lambda = 0.0);
};


%shared_ptr(HullWhite)
class HullWhite : public Vasicek {
  public:
    HullWhite(const Handle<YieldTermStructure>& termStructure,
                 Real a = 0.1, Real sigma = 0.01);
                 
    static Rate convexityBias(Real futurePrice, Time t, Time T, 
                              Real sigma, Real a);

    // TermStructureConsistentModel
    const Handle<YieldTermStructure>& termStructure() const;
};

%shared_ptr(BlackKarasinski)
class BlackKarasinski : public ShortRateModel {
  public:
    BlackKarasinski(const Handle<YieldTermStructure>& termStructure,
                       Real a = 0.1, Real sigma = 0.1);

    // TermStructureConsistentModel
    const Handle<YieldTermStructure>& termStructure() const;
};

%shared_ptr(G2)
class G2 : public ShortRateModel {
  public:
    G2(const Handle<YieldTermStructure>& termStructure,
          Real a = 0.1, Real sigma = 0.01, Real b = 0.1,
          Real eta = 0.01, Real rho = -0.75);

    // TermStructureConsistentModel
    const Handle<YieldTermStructure>& termStructure() const;
};


// pricing engines for calibration helpers
%{
using QuantLib::JamshidianSwaptionEngine;
using QuantLib::TreeSwaptionEngine;
using QuantLib::AnalyticCapFloorEngine;
using QuantLib::TreeCapFloorEngine;
using QuantLib::G2SwaptionEngine;
using QuantLib::FdG2SwaptionEngine;
using QuantLib::FdHullWhiteSwaptionEngine;
%}

%shared_ptr(JamshidianSwaptionEngine)
class JamshidianSwaptionEngine : public PricingEngine {
  public:
    JamshidianSwaptionEngine(
                         const boost::shared_ptr<OneFactorAffineModel>& model,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>());
};

%shared_ptr(TreeSwaptionEngine)
class TreeSwaptionEngine : public PricingEngine {
  public:
    TreeSwaptionEngine(const boost::shared_ptr<ShortRateModel>& model,
                       Size timeSteps,
                       const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>());
    TreeSwaptionEngine(const boost::shared_ptr<ShortRateModel>& model,
                       const TimeGrid& grid,
                       const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>());
    TreeSwaptionEngine(const Handle<ShortRateModel>& model,
                       Size timeSteps,
                       const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>());
};

%shared_ptr(AnalyticCapFloorEngine)
class AnalyticCapFloorEngine : public PricingEngine {
  public:
    AnalyticCapFloorEngine(const boost::shared_ptr<OneFactorAffineModel>& model,
                           const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>());
};

%shared_ptr(TreeCapFloorEngine)
class TreeCapFloorEngine : public PricingEngine {
  public:
    TreeCapFloorEngine(const boost::shared_ptr<ShortRateModel>& model,
                       Size timeSteps,
                       const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>());
    TreeCapFloorEngine(const boost::shared_ptr<ShortRateModel>& model,
                       const TimeGrid& grid,
                       const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>());
};

%shared_ptr(G2SwaptionEngine)
class G2SwaptionEngine : public PricingEngine {
  public:
    G2SwaptionEngine(const boost::shared_ptr<G2>& model,
                     Real range, Size intervals);
};


%shared_ptr(FdG2SwaptionEngine)
class FdG2SwaptionEngine : public PricingEngine {
  public:
    FdG2SwaptionEngine(const boost::shared_ptr<G2>& model,
                       Size tGrid = 100, Size xGrid = 50, Size yGrid = 50,
                       Size dampingSteps = 0, Real invEps = 1e-5,
                       const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer());
};

%shared_ptr(FdHullWhiteSwaptionEngine)
class FdHullWhiteSwaptionEngine : public PricingEngine {
  public:
    FdHullWhiteSwaptionEngine(const boost::shared_ptr<HullWhite>& model,
                              Size tGrid = 100, Size xGrid = 100,
                              Size dampingSteps = 0, Real invEps = 1e-5,
                              const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas());
};

#endif
