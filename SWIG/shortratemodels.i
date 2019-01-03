
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
class HullWhite : public Vasicek, public TermStructureConsistentModel {
  public:
    HullWhite(const Handle<YieldTermStructure>& termStructure,
                 Real a = 0.1, Real sigma = 0.01);
                 
    static Rate convexityBias(Real futurePrice, Time t, Time T, 
                              Real sigma, Real a);
};

%shared_ptr(BlackKarasinski)
class BlackKarasinski : public ShortRateModel, public TermStructureConsistentModel {
  public:
    BlackKarasinski(const Handle<YieldTermStructure>& termStructure,
                       Real a = 0.1, Real sigma = 0.1);
};

%shared_ptr(G2)
class G2 : public ShortRateModel, public TermStructureConsistentModel {
  public:
    G2(const Handle<YieldTermStructure>& termStructure,
          Real a = 0.1, Real sigma = 0.01, Real b = 0.1,
          Real eta = 0.01, Real rho = -0.75);
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
typedef boost::shared_ptr<PricingEngine> JamshidianSwaptionEnginePtr;
typedef boost::shared_ptr<PricingEngine> TreeSwaptionEnginePtr;
typedef boost::shared_ptr<PricingEngine> AnalyticCapFloorEnginePtr;
typedef boost::shared_ptr<PricingEngine> TreeCapFloorEnginePtr;
typedef boost::shared_ptr<PricingEngine> G2SwaptionEnginePtr;
typedef boost::shared_ptr<PricingEngine> FdG2SwaptionEnginePtr;
typedef boost::shared_ptr<PricingEngine> FdHullWhiteSwaptionEnginePtr;
%}

%rename(JamshidianSwaptionEngine) JamshidianSwaptionEnginePtr;
class JamshidianSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        JamshidianSwaptionEnginePtr(
                         const boost::shared_ptr<OneFactorAffineModel>& model,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>()) {

            return new JamshidianSwaptionEnginePtr(
                               new JamshidianSwaptionEngine(model,termStructure));
        }
    }
};

%rename(TreeSwaptionEngine) TreeSwaptionEnginePtr;
class TreeSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        TreeSwaptionEnginePtr(
                         const boost::shared_ptr<ShortRateModel>& model,
                         Size timeSteps,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>()) {
            return new TreeSwaptionEnginePtr(
                       new TreeSwaptionEngine(model,timeSteps,termStructure));
        }
        TreeSwaptionEnginePtr(
                         const boost::shared_ptr<ShortRateModel>& model,
                         const TimeGrid& grid,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>()) {
            return new TreeSwaptionEnginePtr(
                            new TreeSwaptionEngine(model,grid,termStructure));
        }
        TreeSwaptionEnginePtr(
                         const Handle<ShortRateModel>& model,
                         Size timeSteps,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>()) {
            return new TreeSwaptionEnginePtr(
                       new TreeSwaptionEngine(model,timeSteps,termStructure));
        }
    }
};

%rename(AnalyticCapFloorEngine) AnalyticCapFloorEnginePtr;
class AnalyticCapFloorEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticCapFloorEnginePtr(
                         const boost::shared_ptr<OneFactorAffineModel>& model,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>()) {
            return new AnalyticCapFloorEnginePtr(
                                 new AnalyticCapFloorEngine(model,termStructure));
        }
    }
};

%rename(TreeCapFloorEngine) TreeCapFloorEnginePtr;
class TreeCapFloorEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        TreeCapFloorEnginePtr(
                         const boost::shared_ptr<ShortRateModel>& model,
                         Size timeSteps,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>()) {
            return new TreeCapFloorEnginePtr(
                       new TreeCapFloorEngine(model,timeSteps,termStructure));
        }
        TreeCapFloorEnginePtr(
                         const boost::shared_ptr<ShortRateModel>& model,
                         const TimeGrid& grid,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>()) {
            return new TreeCapFloorEnginePtr(
                            new TreeCapFloorEngine(model,grid,termStructure));
        }
    }
};

%rename(G2SwaptionEngine) G2SwaptionEnginePtr;
class G2SwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        G2SwaptionEnginePtr(const boost::shared_ptr<G2>& model,
                            Real range, Size intervals) {
            return new G2SwaptionEnginePtr(
                                    new G2SwaptionEngine(model,range,intervals));
        }
    }
};


%rename(FdG2SwaptionEngine) FdG2SwaptionEnginePtr;
class FdG2SwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdG2SwaptionEnginePtr(const boost::shared_ptr<G2>& model,
                            Size tGrid = 100, Size xGrid = 50, Size yGrid = 50,
							Size dampingSteps = 0, Real invEps = 1e-5,
							const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer()) {
            return new FdG2SwaptionEnginePtr(
                                    new FdG2SwaptionEngine(model,tGrid,xGrid,yGrid,dampingSteps,invEps,schemeDesc));
        }
    }
};

%rename(FdHullWhiteSwaptionEngine) FdHullWhiteSwaptionEnginePtr;
class FdHullWhiteSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdHullWhiteSwaptionEnginePtr(const boost::shared_ptr<HullWhite>& model,
                            Size tGrid = 100, Size xGrid = 100,
							Size dampingSteps = 0, Real invEps = 1e-5,
							const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas()) {
            return new FdHullWhiteSwaptionEnginePtr(
                                    new FdHullWhiteSwaptionEngine(model,tGrid,xGrid,dampingSteps,invEps,schemeDesc));
        }
    }
};

#endif
