
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2007, 2009 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2007 Luis Cota
 Copyright (C) 2015 Gouthaman Balaraman

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

%ignore ShortRateModel;
class ShortRateModel {
    #if defined(SWIGRUBY)
    %rename("calibrate!") calibrate;
    #elif defined(SWIGCSHARP)
    %rename("parameters") params;
    #endif
  public:
    Array params() const;
    void calibrate(
        const std::vector<boost::shared_ptr<CalibrationHelper> >&,
        OptimizationMethod&, const EndCriteria &,
        const Constraint& constraint = Constraint(),
        const std::vector<Real>& weights = std::vector<Real>(),
        const std::vector<bool> & fixParameters = std::vector<bool>());
    Real value(const Array& params,
               const std::vector<boost::shared_ptr<CalibrationHelper> >&);
    EndCriteria::Type endCriteria() const;
    void setParams(const Array& params);
};

%template(ShortRateModel) boost::shared_ptr<ShortRateModel>;
IsObservable(boost::shared_ptr<ShortRateModel>);

%template(ShortRateModelHandle) Handle<ShortRateModel>;
IsObservable(Handle<ShortRateModel>);
%template(RelinkableShortRateModelHandle)
RelinkableHandle<ShortRateModel>;

// actual models

%{
    using QuantLib::Vasicek;
    using QuantLib::HullWhite;
    using QuantLib::BlackKarasinski;
    using QuantLib::G2;
    typedef boost::shared_ptr<ShortRateModel> VasicekPtr;
    typedef boost::shared_ptr<ShortRateModel> HullWhitePtr;
    typedef boost::shared_ptr<ShortRateModel> BlackKarasinskiPtr;
    typedef boost::shared_ptr<ShortRateModel> G2Ptr;
%}

%rename(Vasicek) VasicekPtr;
class VasicekPtr : public boost::shared_ptr<ShortRateModel> {
  public:
    %extend {
        VasicekPtr(Rate r0 = 0.05,
                   Real a = 0.1,
                   Real b = 0.05,
                   Real sigma = 0.01,
                   Real lambda = 0.0) {
            return new VasicekPtr(new Vasicek(r0, a, b, sigma, lambda));
        }
        DiscountFactor discount(Time t) const {
            return boost::dynamic_pointer_cast<Vasicek>(*self)->discount(t);
        }
        
        Real discountBond(Time now, Time maturity, Rate rate) const {
            return boost::dynamic_pointer_cast<Vasicek>(*self)->discountBond(now, maturity, rate);
        }
    }
};


%rename(HullWhite) HullWhitePtr;
class HullWhitePtr : public boost::shared_ptr<ShortRateModel> {
  public:
    %extend {
        HullWhitePtr(const Handle<YieldTermStructure>& termStructure,
                     Real a = 0.1, Real sigma = 0.01) {
            return new HullWhitePtr(
                new HullWhite(termStructure, a, sigma));
        }
        DiscountFactor discount(Time t) const {
            return boost::dynamic_pointer_cast<HullWhite>(*self)->discount(t);
        }
        Real discountBond(Time now, Time maturity, Rate rate) const{
            return boost::dynamic_pointer_cast<HullWhite>(*self)->
                                            discountBond(now, maturity, rate);  
        }
        static Rate convexityBias(Real futurePrice, Time t, Time T, 
                                  Real sigma, Real a) {
            return QuantLib::HullWhite::convexityBias(futurePrice, t, T, 
                                                      sigma, a);
        }
    }
};

%rename(BlackKarasinski) BlackKarasinskiPtr;
class BlackKarasinskiPtr : public boost::shared_ptr<ShortRateModel> {
  public:
    %extend {
        BlackKarasinskiPtr(const Handle<YieldTermStructure>& termStructure,
                           Real a = 0.1, Real sigma = 0.1) {
            return new BlackKarasinskiPtr(
                new BlackKarasinski(termStructure, a, sigma));
        }
    }
};

%rename(G2) G2Ptr;
class G2Ptr : public boost::shared_ptr<ShortRateModel> {
  public:
    %extend {
        G2Ptr(const Handle<YieldTermStructure>& termStructure,
              Real a = 0.1, Real sigma = 0.01, Real b = 0.1,
              Real eta = 0.01, Real rho = -0.75) {
            return new G2Ptr(new G2(termStructure, a, sigma, b, eta, rho));
        }
        DiscountFactor discount(Time t) const {
            return boost::dynamic_pointer_cast<G2>(*self)->discount(t);
        }
        
        Real discountBond(Time t, Time T, Rate x, Rate y) const{
            return boost::dynamic_pointer_cast<G2>(*self)->discountBond(t, T, x, y);
        }
    }
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
                         const boost::shared_ptr<ShortRateModel>& model,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>()) {
            using QuantLib::OneFactorAffineModel;
            boost::shared_ptr<OneFactorAffineModel> m =
                 boost::dynamic_pointer_cast<OneFactorAffineModel>(model);
            QL_REQUIRE(model, "affine model required");
            return new JamshidianSwaptionEnginePtr(
                               new JamshidianSwaptionEngine(m,termStructure));
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
                         const boost::shared_ptr<ShortRateModel>& model,
                         const Handle<YieldTermStructure>& termStructure =
                                                Handle<YieldTermStructure>()) {
            using QuantLib::OneFactorAffineModel;
            boost::shared_ptr<OneFactorAffineModel> m =
                 boost::dynamic_pointer_cast<OneFactorAffineModel>(model);
            QL_REQUIRE(model, "affine model required");
            return new AnalyticCapFloorEnginePtr(
                                 new AnalyticCapFloorEngine(m,termStructure));
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
        G2SwaptionEnginePtr(const boost::shared_ptr<ShortRateModel>& model,
                            Real range, Size intervals) {
            boost::shared_ptr<G2> g2 =
                boost::dynamic_pointer_cast<G2>(model);
            return new G2SwaptionEnginePtr(
                                    new G2SwaptionEngine(g2,range,intervals));
        }
    }
};


%rename(FdG2SwaptionEngine) FdG2SwaptionEnginePtr;
class FdG2SwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdG2SwaptionEnginePtr(const boost::shared_ptr<ShortRateModel>& model,
                            Size tGrid = 100, Size xGrid = 50, Size yGrid = 50,
							Size dampingSteps = 0, Real invEps = 1e-5,
							const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer()) {
            boost::shared_ptr<G2> g2 =
                boost::dynamic_pointer_cast<G2>(model);
            return new FdG2SwaptionEnginePtr(
                                    new FdG2SwaptionEngine(g2,tGrid,xGrid,yGrid,dampingSteps,invEps,schemeDesc));
        }
    }
};

%rename(FdHullWhiteSwaptionEngine) FdHullWhiteSwaptionEnginePtr;
class FdHullWhiteSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdHullWhiteSwaptionEnginePtr(const boost::shared_ptr<ShortRateModel>& model,
                            Size tGrid = 100, Size xGrid = 100,
							Size dampingSteps = 0, Real invEps = 1e-5,
							const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas()) {
            boost::shared_ptr<HullWhite> hw =
                boost::dynamic_pointer_cast<HullWhite>(model);
            return new FdHullWhiteSwaptionEnginePtr(
                                    new FdHullWhiteSwaptionEngine(hw,tGrid,xGrid,dampingSteps,invEps,schemeDesc));
        }
    }
};

#endif
