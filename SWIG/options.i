
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2008 Tito Ingargiola
 Copyright (C) 2010, 2012, 2018 Klaus Spanderen
 Copyright (C) 2015 Thema Consulting SA
 Copyright (C) 2016 Gouthaman Balaraman
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

#ifndef quantlib_options_i
#define quantlib_options_i

%include common.i
%include exercise.i
%include stochasticprocess.i
%include instruments.i
%include stl.i
%include linearalgebra.i
%include calibrationhelpers.i
%include grid.i
%include parameter.i
%include vectors.i

// payoff

%{
using QuantLib::Payoff;
using QuantLib::TypePayoff;
using QuantLib::StrikedTypePayoff;
%}

%shared_ptr(Payoff);
class Payoff {
    #if defined(SWIGCSHARP)
    %rename(call) operator();
    #endif
  public:
    Real operator()(Real price) const;
  private:
    Payoff();
};

// option and barrier types
%{
using QuantLib::Option;
using QuantLib::Barrier;
%}

%shared_ptr(Option)
class Option : public Instrument {
  public:
    enum Type { Put = -1,
                Call = 1
    };
    boost::shared_ptr<Payoff> payoff();
    boost::shared_ptr<Exercise> exercise();
  private:
    Option();
};

struct Barrier {
    enum Type { DownIn, UpIn, DownOut, UpOut };
};

struct DoubleBarrier {
    enum Type { KnockIn, KnockOut, KIKO, KOKI };
};


#if defined(SWIGR)
%Rruntime %{
setMethod("summary", "_p_VanillaOption",
function(object) {object$freeze()
ans <- c(value=object$NPV(), delta=object$delta(),
gamma=object$gamma(), vega=object$vega(),
theta=object$theta(), rho=object$rho(),
divRho=object$dividendRho())
object$unfreeze()
ans
})

setMethod("summary", "_p_DividendVanillaOption",
function(object) {object$freeze()
ans <- c(value=object$NPV(), delta=object$delta(),
gamma=object$gamma(), vega=object$vega(),
theta=object$theta(), rho=object$rho(),
divRho=object$dividendRho())
object$unfreeze()
ans
})

%}
#endif


%shared_ptr(TypePayoff)
class TypePayoff : public Payoff {
  public:
    Option::Type optionType();
  private:
    TypePayoff();
};

%shared_ptr(StrikedTypePayoff)
class StrikedTypePayoff : public TypePayoff
{
  public:
    Real strike();
  private:
    StrikedTypePayoff();
};


// plain option and engines

%{
using QuantLib::OneAssetOption;
using QuantLib::VanillaOption;
using QuantLib::ForwardVanillaOption;
%}

%shared_ptr(OneAssetOption)
class OneAssetOption : public Option {
  public:
    OneAssetOption(const boost::shared_ptr<Payoff>&,
                   const boost::shared_ptr<Exercise>&);
    Real delta() const;
    Real deltaForward() const;
    Real elasticity() const;
    Real gamma() const;
    Real theta() const;
    Real thetaPerDay() const;
    Real vega() const;
    Real rho() const;
    Real dividendRho() const;
    Real strikeSensitivity() const;
    Real itmCashProbability() const;
};

%shared_ptr(VanillaOption)
class VanillaOption : public OneAssetOption {
  public:
    VanillaOption(
            const boost::shared_ptr<StrikedTypePayoff>& payoff,
            const boost::shared_ptr<Exercise>& exercise);

    Volatility impliedVolatility(
                         Real targetValue,
                         const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                         Real accuracy = 1.0e-4,
                         Size maxEvaluations = 100,
                         Volatility minVol = 1.0e-4,
                         Volatility maxVol = 4.0);
    %extend{                     
        SampledCurve priceCurve() {
            return self->result<SampledCurve>("priceCurve");
        }
    }
};


%{
using QuantLib::EuropeanOption;
%}


%shared_ptr(EuropeanOption)
class EuropeanOption : public VanillaOption {
  public:
    EuropeanOption(
            const boost::shared_ptr<StrikedTypePayoff>& payoff,
            const boost::shared_ptr<Exercise>& exercise);
};

// ForwardVanillaOption

%{
using QuantLib::ForwardVanillaOption;
%}

%shared_ptr(ForwardVanillaOption)
class ForwardVanillaOption : public OneAssetOption {
  public:
        ForwardVanillaOption(
                Real moneyness,
                Date resetDate,
                const boost::shared_ptr<StrikedTypePayoff>& payoff,
                const boost::shared_ptr<Exercise>& exercise);
};

// QuantoVanillaOption

%{
using QuantLib::QuantoVanillaOption;
%}

%shared_ptr(QuantoVanillaOption)
class QuantoVanillaOption : public OneAssetOption {
  public:
    QuantoVanillaOption(
            const boost::shared_ptr<StrikedTypePayoff>& payoff,
            const boost::shared_ptr<Exercise>& exercise);
    Real qvega();
    Real qrho();
    Real qlambda();
};

%{
using QuantLib::QuantoForwardVanillaOption;
%}

%shared_ptr(QuantoForwardVanillaOption)
class QuantoForwardVanillaOption : public ForwardVanillaOption {
  public:
    QuantoForwardVanillaOption(
            Real moneyness,
            Date resetDate,
            const boost::shared_ptr<StrikedTypePayoff>& payoff,
            const boost::shared_ptr<Exercise>& exercise);
};

%{
using QuantLib::MultiAssetOption;
%}
%shared_ptr(MultiAssetOption)
class MultiAssetOption : public Option {
  public:
    Real delta();
    Real gamma();
    Real theta();
    Real vega();
    Real rho();
    Real dividendRho();
};

// European engines

%{
using QuantLib::AnalyticEuropeanEngine;
typedef boost::shared_ptr<PricingEngine> AnalyticEuropeanEnginePtr;
%}

%rename(AnalyticEuropeanEngine) AnalyticEuropeanEnginePtr;
class AnalyticEuropeanEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticEuropeanEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new AnalyticEuropeanEnginePtr(
                                       new AnalyticEuropeanEngine(process));
        }
    }
};


%{
using QuantLib::HestonModel;
%}

%shared_ptr(HestonModel)
class HestonModel : public CalibratedModel {
  public:
    HestonModel(const boost::shared_ptr<HestonProcess>&  process);
    Real theta() const;
    Real kappa() const;
    Real sigma() const;
    Real rho() const;
    Real v0() const;
};



%{
using QuantLib::PiecewiseTimeDependentHestonModel;
%}

%shared_ptr(PiecewiseTimeDependentHestonModel)
class PiecewiseTimeDependentHestonModel : public CalibratedModel {
    public:
         PiecewiseTimeDependentHestonModel(
              const Handle<YieldTermStructure>& riskFreeRate,
              const Handle<YieldTermStructure>& dividendYield,
              const Handle<Quote>& s0,
              Real v0,
              const Parameter& theta,
              const Parameter& kappa,
              const Parameter& sigma,
              const Parameter& rho,
              const TimeGrid& timeGrid);
    
        Real theta(Time t) const;
        Real kappa(Time t) const;
        Real sigma(Time t) const;
        Real rho(Time t)   const;
        Real v0()          const;
        Real s0()          const;
        const TimeGrid& timeGrid() const;
        const Handle<YieldTermStructure>& dividendYield() const;
        const Handle<YieldTermStructure>& riskFreeRate() const;

};



%{
using QuantLib::AnalyticHestonEngine;
typedef boost::shared_ptr<PricingEngine> AnalyticHestonEnginePtr;
%}

%rename(AnalyticHestonEngine) AnalyticHestonEnginePtr;
class AnalyticHestonEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticHestonEnginePtr(const boost::shared_ptr<HestonModel>& model, 
                                Size integrationOrder = 144) {

            return new AnalyticHestonEnginePtr(
                new AnalyticHestonEngine(model, integrationOrder));
        }

        AnalyticHestonEnginePtr(const boost::shared_ptr<HestonModel>& model, 
                                Real relTolerance,
                                Size maxEvaluations) {

            return new AnalyticHestonEnginePtr(
                new AnalyticHestonEngine(model, relTolerance,maxEvaluations));
        }
    }
};

%{
using QuantLib::COSHestonEngine;
typedef boost::shared_ptr<PricingEngine> COSHestonEnginePtr;
%}

%rename(COSHestonEngine) COSHestonEnginePtr;
class COSHestonEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        COSHestonEnginePtr(const boost::shared_ptr<HestonModel>& model, 
                           Real L = 16, Size N=200) {

            return new COSHestonEnginePtr(
                new COSHestonEngine(model, L, N));
        }
    }
};


%{
using QuantLib::AnalyticPTDHestonEngine;
typedef boost::shared_ptr<PricingEngine> AnalyticPTDHestonEnginePtr;
%}

%rename(AnalyticPTDHestonEngine) AnalyticPTDHestonEnginePtr;
class AnalyticPTDHestonEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticPTDHestonEnginePtr(
            const boost::shared_ptr<PiecewiseTimeDependentHestonModel>& model,
            Real relTolerance, Size maxEvaluations) {
                return new AnalyticPTDHestonEnginePtr(
                    new AnalyticPTDHestonEngine(model, relTolerance, maxEvaluations)
            );
        }

        // Constructor using Laguerre integration
        // and Gatheral's version of complex log.
        AnalyticPTDHestonEnginePtr(
            const boost::shared_ptr<PiecewiseTimeDependentHestonModel>& model,
            Size integrationOrder = 144) {
                return new AnalyticPTDHestonEnginePtr(
                    new AnalyticPTDHestonEngine(model, integrationOrder)
            );
        }
    
    }
};


#if defined(SWIGPYTHON)
%rename(lambda_parameter) lambda;
#endif

%{
using QuantLib::BatesModel;
%}

%shared_ptr(BatesModel)
class BatesModel : public HestonModel {
  public:
    BatesModel(const boost::shared_ptr<BatesProcess>&  process);
    Real nu() const;
    Real delta() const;
    Real lambda() const;
};


%{
using QuantLib::BatesEngine;
typedef boost::shared_ptr<PricingEngine> BatesEnginePtr;
%}

%rename(BatesEngine) BatesEnginePtr;
class BatesEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        BatesEnginePtr(const boost::shared_ptr<BatesModel>& model,
                       Size integrationOrder = 144) {

            return new BatesEnginePtr(
                new BatesEngine(model, integrationOrder));
        }

        BatesEnginePtr(const boost::shared_ptr<BatesModel>& model, 
                       Real relTolerance,
                       Size maxEvaluations) {

            return new BatesEnginePtr(
                new BatesEngine(model, relTolerance,maxEvaluations));
        }
    }
};


%{
using QuantLib::IntegralEngine;
typedef boost::shared_ptr<PricingEngine> IntegralEnginePtr;
%}

%rename(IntegralEngine) IntegralEnginePtr;
class IntegralEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        IntegralEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new IntegralEnginePtr(new IntegralEngine(process));
        }
    }
};


%{
using QuantLib::FDBermudanEngine;
typedef boost::shared_ptr<PricingEngine> FDBermudanEnginePtr;
%}

%rename(FDBermudanEngine) FDBermudanEnginePtr;
class FDBermudanEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FDBermudanEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                            Size timeSteps = 100, Size gridPoints = 100,
                            bool timeDependent = false) {
            return new FDBermudanEnginePtr(
                            new FDBermudanEngine<>(process,timeSteps,
                                                   gridPoints,timeDependent));
        }
    }
};

%{
using QuantLib::FDEuropeanEngine;
typedef boost::shared_ptr<PricingEngine> FDEuropeanEnginePtr;
%}

%rename(FDEuropeanEngine) FDEuropeanEnginePtr;
class FDEuropeanEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FDEuropeanEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess> process,
                            Size timeSteps = 100, Size gridPoints = 100,
                            bool timeDependent = false) {
            return new FDEuropeanEnginePtr(
                            new FDEuropeanEngine<>(process,timeSteps,
                                                   gridPoints,timeDependent));
        }
    }
};

%{
using QuantLib::BinomialVanillaEngine;
using QuantLib::CoxRossRubinstein;
using QuantLib::JarrowRudd;
using QuantLib::AdditiveEQPBinomialTree;
using QuantLib::Trigeorgis;
using QuantLib::Tian;
using QuantLib::LeisenReimer;
using QuantLib::Joshi4;
typedef boost::shared_ptr<PricingEngine> BinomialVanillaEnginePtr;
%}

%rename(BinomialVanillaEngine) BinomialVanillaEnginePtr;
class BinomialVanillaEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        BinomialVanillaEnginePtr(
                             const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                             const std::string& type,
                             Size steps) {
            std::string s = boost::algorithm::to_lower_copy(type);
            if (s == "crr" || s == "coxrossrubinstein")
                return new BinomialVanillaEnginePtr(
                    new BinomialVanillaEngine<CoxRossRubinstein>(
                                                            process,steps));
            else if (s == "jr" || s == "jarrowrudd")
                return new BinomialVanillaEnginePtr(
                    new BinomialVanillaEngine<JarrowRudd>(process,steps));
            else if (s == "eqp" || s == "additiveeqpbinomialtree")
                return new BinomialVanillaEnginePtr(
                    new BinomialVanillaEngine<AdditiveEQPBinomialTree>(
                                                            process,steps));
            else if (s == "trigeorgis")
                return new BinomialVanillaEnginePtr(
                    new BinomialVanillaEngine<Trigeorgis>(process,steps));
            else if (s == "tian")
                return new BinomialVanillaEnginePtr(
                    new BinomialVanillaEngine<Tian>(process,steps));
            else if (s == "lr" || s == "leisenreimer")
                return new BinomialVanillaEnginePtr(
                    new BinomialVanillaEngine<LeisenReimer>(process,steps));
            else if (s == "j4" || s == "joshi4")
                return new BinomialVanillaEnginePtr(
                    new BinomialVanillaEngine<Joshi4>(process,steps));
            else
                QL_FAIL("unknown binomial engine type: "+s);
        }
    }
};


%{
using QuantLib::MCEuropeanEngine;
using QuantLib::MCAmericanEngine;
using QuantLib::PseudoRandom;
using QuantLib::LowDiscrepancy;
using QuantLib::LsmBasisSystem;
typedef boost::shared_ptr<PricingEngine> MCEuropeanEnginePtr;
typedef boost::shared_ptr<PricingEngine> MCAmericanEnginePtr;
%}

struct LsmBasisSystem {
    enum PolynomType  {Monomial, Laguerre, Hermite, Hyperbolic,
                           Legendre, Chebyshev, Chebyshev2nd };
};

%rename(MCEuropeanEngine) MCEuropeanEnginePtr;
class MCEuropeanEnginePtr : public boost::shared_ptr<PricingEngine> {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCEuropeanEnginePtr;
    #endif
  public:
    %extend {
        MCEuropeanEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                            const std::string& traits,
                            intOrNull timeSteps = Null<Size>(),
                            intOrNull timeStepsPerYear = Null<Size>(),
                            bool brownianBridge = false,
                            bool antitheticVariate = false,
                            intOrNull requiredSamples = Null<Size>(),
                            doubleOrNull requiredTolerance = Null<Real>(),
                            intOrNull maxSamples = Null<Size>(),
                            BigInteger seed = 0) {
            std::string s = boost::algorithm::to_lower_copy(traits);
            QL_REQUIRE(Size(timeSteps) != Null<Size>() ||
                       Size(timeStepsPerYear) != Null<Size>(),
                       "number of steps not specified");
            if (s == "pseudorandom" || s == "pr")
                return new MCEuropeanEnginePtr(
                         new MCEuropeanEngine<PseudoRandom>(process,
                                                            timeSteps,
                                                            timeStepsPerYear,
                                                            brownianBridge,
                                                            antitheticVariate,
                                                            requiredSamples,
                                                            requiredTolerance,
                                                            maxSamples,
                                                            seed));
            else if (s == "lowdiscrepancy" || s == "ld")
                return new MCEuropeanEnginePtr(
                       new MCEuropeanEngine<LowDiscrepancy>(process,
                                                            timeSteps,
                                                            timeStepsPerYear,
                                                            brownianBridge,
                                                            antitheticVariate,
                                                            requiredSamples,
                                                            requiredTolerance,
                                                            maxSamples,
                                                            seed));
            else
                QL_FAIL("unknown Monte Carlo engine type: "+s);
        }
    }
};

%rename(MCAmericanEngine) MCAmericanEnginePtr;
class MCAmericanEnginePtr : public boost::shared_ptr<PricingEngine> {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCAmericanEnginePtr;
    #endif
  public:
    %extend {
		static const VanillaSwap::Type Receiver = VanillaSwap::Receiver;
        static const VanillaSwap::Type Payer = VanillaSwap::Payer;
		
        MCAmericanEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                            const std::string& traits,
                            intOrNull timeSteps = Null<Size>(),
                            intOrNull timeStepsPerYear = Null<Size>(),
                            bool antitheticVariate = false,
                            bool controlVariate = false,							
                            intOrNull requiredSamples = Null<Size>(),
                            doubleOrNull requiredTolerance = Null<Real>(),
                            intOrNull maxSamples = Null<Size>(),
                            BigInteger seed = 0,
							intOrNull polynomOrder = 2, 
							LsmBasisSystem::PolynomType polynomType = LsmBasisSystem::Monomial,
							int nCalibrationSamples = 2048,
							boost::optional<bool> antitheticVariateCalibration = boost::none,
							BigNatural seedCalibration = Null<Size>()) {
            std::string s = boost::algorithm::to_lower_copy(traits);
            QL_REQUIRE(Size(timeSteps) != Null<Size>() ||
                       Size(timeStepsPerYear) != Null<Size>(),
                       "number of steps not specified");
            if (s == "pseudorandom" || s == "pr")
                return new MCAmericanEnginePtr(
                         new MCAmericanEngine<PseudoRandom>(process,
                                                            timeSteps,
                                                            timeStepsPerYear,
                                                            antitheticVariate,
															controlVariate,
                                                            requiredSamples,
                                                            requiredTolerance,
                                                            maxSamples,
                                                            seed,
															polynomOrder,
															polynomType,
															nCalibrationSamples,
															antitheticVariateCalibration,
															seedCalibration));
            else if (s == "lowdiscrepancy" || s == "ld")
                return new MCAmericanEnginePtr(
                       new MCAmericanEngine<LowDiscrepancy>(process,
                                                            timeSteps,
                                                            timeStepsPerYear,
                                                            antitheticVariate,
															controlVariate,
                                                            requiredSamples,
                                                            requiredTolerance,
                                                            maxSamples,
                                                            seed,
															polynomOrder,
															polynomType,
															nCalibrationSamples,
															antitheticVariateCalibration,
															seedCalibration));
            else
                QL_FAIL("unknown Monte Carlo engine type: "+s);
        }
    }
};

// American engines

%{
using QuantLib::FDAmericanEngine;
using QuantLib::FDShoutEngine;
typedef boost::shared_ptr<PricingEngine> FDAmericanEnginePtr;
typedef boost::shared_ptr<PricingEngine> FDShoutEnginePtr;
%}

%rename(FDAmericanEngine) FDAmericanEnginePtr;
class FDAmericanEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FDAmericanEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                            Size timeSteps = 100, Size gridPoints = 100,
                            bool timeDependent = false) {
            return new FDAmericanEnginePtr(
                            new FDAmericanEngine<>(process,timeSteps,
                                                   gridPoints,timeDependent));
        }
    }
};

%{
using QuantLib::ContinuousArithmeticAsianLevyEngine;
typedef boost::shared_ptr<PricingEngine> ContinuousArithmeticAsianLevyEnginePtr;
%}

%rename(ContinuousArithmeticAsianLevyEngine) ContinuousArithmeticAsianLevyEnginePtr;
class ContinuousArithmeticAsianLevyEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        ContinuousArithmeticAsianLevyEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                               const Handle<Quote>& runningAverage,
                                               const Date& startDate) {
            return new ContinuousArithmeticAsianLevyEnginePtr(
                            new ContinuousArithmeticAsianLevyEngine(process,runningAverage, startDate));
        }
    }
};

%{
using QuantLib::FdBlackScholesAsianEngine;
typedef boost::shared_ptr<PricingEngine> FdBlackScholesAsianEnginePtr;
%}

%rename(FdBlackScholesAsianEngine) FdBlackScholesAsianEnginePtr;
class FdBlackScholesAsianEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdBlackScholesAsianEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                     Size tGrid, Size xGrid, Size aGrid) {
            return new FdBlackScholesAsianEnginePtr(
                            new FdBlackScholesAsianEngine(process,tGrid, xGrid, aGrid));
        }
    }
};



%rename(FDShoutEngine) FDShoutEnginePtr;
class FDShoutEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FDShoutEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                         Size timeSteps = 100, Size gridPoints = 100,
                         bool timeDependent = false) {

            return new FDShoutEnginePtr(
                               new FDShoutEngine<>(process,timeSteps,
                                                   gridPoints,timeDependent));
        }
    }
};


%{
using QuantLib::BaroneAdesiWhaleyApproximationEngine;
typedef boost::shared_ptr<PricingEngine>
    BaroneAdesiWhaleyApproximationEnginePtr;
%}

%rename(BaroneAdesiWhaleyEngine) BaroneAdesiWhaleyApproximationEnginePtr;
class BaroneAdesiWhaleyApproximationEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        BaroneAdesiWhaleyApproximationEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new BaroneAdesiWhaleyApproximationEnginePtr(
                         new BaroneAdesiWhaleyApproximationEngine(process));
        }
    }
};


%{
using QuantLib::BjerksundStenslandApproximationEngine;
typedef boost::shared_ptr<PricingEngine>
    BjerksundStenslandApproximationEnginePtr;
%}

%rename(BjerksundStenslandEngine) BjerksundStenslandApproximationEnginePtr;
class BjerksundStenslandApproximationEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        BjerksundStenslandApproximationEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new BjerksundStenslandApproximationEnginePtr(
                        new BjerksundStenslandApproximationEngine(process));
        }
    }
};

%{
using QuantLib::AnalyticDigitalAmericanEngine;
typedef boost::shared_ptr<PricingEngine> AnalyticDigitalAmericanEnginePtr;
%}

%rename(AnalyticDigitalAmericanEngine) AnalyticDigitalAmericanEnginePtr;
class AnalyticDigitalAmericanEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticDigitalAmericanEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new AnalyticDigitalAmericanEnginePtr(
                                new AnalyticDigitalAmericanEngine(process));
        }
    }
};

%{
using QuantLib::AnalyticDigitalAmericanKOEngine;
typedef boost::shared_ptr<PricingEngine> AnalyticDigitalAmericanKOEnginePtr;
%}

%rename(AnalyticDigitalAmericanKOEngine) AnalyticDigitalAmericanKOEnginePtr;
class AnalyticDigitalAmericanKOEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticDigitalAmericanKOEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new AnalyticDigitalAmericanKOEnginePtr(
                                new AnalyticDigitalAmericanKOEngine(process));
        }
    }
};

// Dividend option

%{
using QuantLib::DividendVanillaOption;
%}


%shared_ptr(DividendVanillaOption)
class DividendVanillaOption : public OneAssetOption {
  public:
    DividendVanillaOption(
            const boost::shared_ptr<StrikedTypePayoff>& payoff,
            const boost::shared_ptr<Exercise>& exercise,
            const std::vector<Date>& dividendDates,
            const std::vector<Real>& dividends);
    Volatility impliedVolatility(
                         Real targetValue,
                         const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                         Real accuracy = 1.0e-4,
                         Size maxEvaluations = 100,
                         Volatility minVol = 1.0e-4,
                         Volatility maxVol = 4.0) {
        return boost::dynamic_pointer_cast<DividendVanillaOption>(*self)
            ->impliedVolatility(targetValue, process, accuracy,
                                maxEvaluations, minVol, maxVol);
    }
};


%{
using QuantLib::AnalyticDividendEuropeanEngine;
typedef boost::shared_ptr<PricingEngine> AnalyticDividendEuropeanEnginePtr;
%}

%rename(AnalyticDividendEuropeanEngine) AnalyticDividendEuropeanEnginePtr;
class AnalyticDividendEuropeanEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticDividendEuropeanEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new AnalyticDividendEuropeanEnginePtr(
                               new AnalyticDividendEuropeanEngine(process));
        }
    }
};

%{
using QuantLib::FDDividendEuropeanEngine;
using QuantLib::FDDividendAmericanEngine;
typedef boost::shared_ptr<PricingEngine> FDDividendEuropeanEnginePtr;
typedef boost::shared_ptr<PricingEngine> FDDividendAmericanEnginePtr;
%}

%rename(FDDividendEuropeanEngine) FDDividendEuropeanEnginePtr;
class FDDividendEuropeanEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FDDividendEuropeanEnginePtr(
                             const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                             Size timeSteps = 100,
                             Size gridPoints = 100,
                             bool timeDependent = false) {
            return new FDDividendEuropeanEnginePtr(
                   new FDDividendEuropeanEngine<>(process,timeSteps,
                                                  gridPoints, timeDependent));
        }
    }
};

%rename(FDDividendAmericanEngine) FDDividendAmericanEnginePtr;
class FDDividendAmericanEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FDDividendAmericanEnginePtr(
                             const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                             Size timeSteps = 100,
                             Size gridPoints = 100,
                             bool timeDependent = false) {
            return new FDDividendAmericanEnginePtr(
                   new FDDividendAmericanEngine<>(process,timeSteps,
                                                  gridPoints, timeDependent));
        }
    }
};


// Barrier option

%{
using QuantLib::BarrierOption;
%}

%shared_ptr(BarrierOption)
class BarrierOption : public OneAssetOption {
  public:
    BarrierOption(
               Barrier::Type barrierType,
               Real barrier,
               Real rebate,
               const boost::shared_ptr<StrikedTypePayoff>& payoff,
               const boost::shared_ptr<Exercise>& exercise);
    Volatility impliedVolatility(
                         Real targetValue,
                         const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                         Real accuracy = 1.0e-4,
                         Size maxEvaluations = 100,
                         Volatility minVol = 1.0e-4,
                         Volatility maxVol = 4.0) {
        return boost::dynamic_pointer_cast<BarrierOption>(*self)
             ->impliedVolatility(targetValue, process, accuracy,
                                 maxEvaluations, minVol, maxVol);
    }
};

// Barrier engines

%{
using QuantLib::AnalyticBarrierEngine;
typedef boost::shared_ptr<PricingEngine> AnalyticBarrierEnginePtr;
%}

%rename(AnalyticBarrierEngine) AnalyticBarrierEnginePtr;
class AnalyticBarrierEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticBarrierEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new AnalyticBarrierEnginePtr(
                                        new AnalyticBarrierEngine(process));
        }
    }
};

%{
using QuantLib::MCBarrierEngine;
typedef boost::shared_ptr<PricingEngine> MCBarrierEnginePtr;
%}

%rename(MCBarrierEngine) MCBarrierEnginePtr;
class MCBarrierEnginePtr : public boost::shared_ptr<PricingEngine> {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCBarrierEnginePtr;
    #endif
  public:
    %extend {
        MCBarrierEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                           const std::string& traits,
                           Size timeSteps = Null<Size>(),
                           Size timeStepsPerYear = Null<Size>(),
                           bool brownianBridge = false,
                           bool antitheticVariate = false,
                           intOrNull requiredSamples = Null<Size>(),
                           doubleOrNull requiredTolerance = Null<Real>(),
                           intOrNull maxSamples = Null<Size>(),
                           bool isBiased = false,
                           BigInteger seed = 0) {
            std::string s = boost::algorithm::to_lower_copy(traits);
            if (s == "pseudorandom" || s == "pr")
                return new MCBarrierEnginePtr(
                         new MCBarrierEngine<PseudoRandom>(process,
                                                           timeSteps,
                                                           timeStepsPerYear,
                                                           brownianBridge,
                                                           antitheticVariate,
                                                           requiredSamples,
                                                           requiredTolerance,
                                                           maxSamples,
                                                           isBiased,
                                                           seed));
            else if (s == "lowdiscrepancy" || s == "ld")
                return new MCBarrierEnginePtr(
                       new MCBarrierEngine<LowDiscrepancy>(process,
                                                           timeSteps,
                                                           timeStepsPerYear,
                                                           brownianBridge,
                                                           antitheticVariate,
                                                           requiredSamples,
                                                           requiredTolerance,
                                                           maxSamples,
                                                           isBiased,
                                                           seed));
            else
                QL_FAIL("unknown Monte Carlo engine type: "+s);
        }
    }
};

%{
using QuantLib::FdmSchemeDesc;
%}

struct FdmSchemeDesc {
  enum FdmSchemeType { HundsdorferType, DouglasType, 
                       CraigSneydType, ModifiedCraigSneydType, 
                       ImplicitEulerType, ExplicitEulerType };
  
  FdmSchemeDesc(FdmSchemeType type, Real theta, Real mu);
  
  const FdmSchemeType type;
  const Real theta, mu;

  // some default scheme descriptions
  static FdmSchemeDesc Douglas();
  static FdmSchemeDesc ImplicitEuler();
  static FdmSchemeDesc ExplicitEuler();
  static FdmSchemeDesc CraigSneyd();
  static FdmSchemeDesc ModifiedCraigSneyd(); 
  static FdmSchemeDesc Hundsdorfer();
  static FdmSchemeDesc ModifiedHundsdorfer();
};

%{
using QuantLib::FdBlackScholesVanillaEngine;
typedef boost::shared_ptr<PricingEngine> FdBlackScholesVanillaEnginePtr;
%}

%rename(FdBlackScholesVanillaEngine) FdBlackScholesVanillaEnginePtr;
class FdBlackScholesVanillaEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdBlackScholesVanillaEnginePtr(
        	const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>()) {
            return new FdBlackScholesVanillaEnginePtr(
                new FdBlackScholesVanillaEngine( 
                    process,tGrid, xGrid, dampingSteps, 
                    schemeDesc, localVol, illegalLocalVolOverwrite));
        }
    }
};

%{
using QuantLib::FdBatesVanillaEngine;
typedef boost::shared_ptr<PricingEngine> FdBatesVanillaEnginePtr;
%}

%rename(FdBatesVanillaEngine) FdBatesVanillaEnginePtr;
class FdBatesVanillaEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdBatesVanillaEnginePtr(
            const boost::shared_ptr<BatesModel>& model,
            Size tGrid = 100, Size xGrid = 100,
            Size vGrid=50, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer()) {
            
            return new FdBatesVanillaEnginePtr(
                new FdBatesVanillaEngine(
                    model, tGrid, xGrid, vGrid, dampingSteps, schemeDesc));
        }
    }
};

%{
using QuantLib::FdHestonVanillaEngine;
typedef boost::shared_ptr<PricingEngine> FdHestonVanillaEnginePtr;
%}

%rename(FdHestonVanillaEngine) FdHestonVanillaEnginePtr;
class FdHestonVanillaEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdHestonVanillaEnginePtr(
            const boost::shared_ptr<HestonModel>& model,
            Size tGrid = 100, Size xGrid = 100,
            Size vGrid = 50, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer()) {
            
            return new FdHestonVanillaEnginePtr(
                new FdHestonVanillaEngine(model, tGrid, xGrid,
                                          vGrid, dampingSteps, schemeDesc));
        }
    }
};

%{
using QuantLib::FdBlackScholesBarrierEngine;
typedef boost::shared_ptr<PricingEngine> FdBlackScholesBarrierEnginePtr;
%}

%rename(FdBlackScholesBarrierEngine) FdBlackScholesBarrierEnginePtr;
class FdBlackScholesBarrierEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdBlackScholesBarrierEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                       Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
                                       const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
                                       bool localVol = false, 
                                       Real illegalLocalVolOverwrite = -Null<Real>()) {
            return new FdBlackScholesBarrierEnginePtr(
                            new FdBlackScholesBarrierEngine(process,
                                                              tGrid,  xGrid, dampingSteps, 
                                                              schemeDesc, localVol,
                                                              illegalLocalVolOverwrite));
        }
  }
};


%{
using QuantLib::AnalyticBinaryBarrierEngine;
typedef boost::shared_ptr<PricingEngine> AnalyticBinaryBarrierEnginePtr;
%}

%rename(AnalyticBinaryBarrierEngine) AnalyticBinaryBarrierEnginePtr;
class AnalyticBinaryBarrierEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticBinaryBarrierEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new AnalyticBinaryBarrierEnginePtr(
                                        new AnalyticBinaryBarrierEngine(process));
        }
    }
};


%{
using QuantLib::BinomialBarrierEngine;
using QuantLib::DiscretizedDermanKaniBarrierOption;
typedef boost::shared_ptr<PricingEngine> BinomialBarrierEnginePtr;
%}

#if defined(SWIGPYTHON)
%feature("docstring") BinomialBarrierEnginePtr "Binomial Engine for barrier options.
Features different binomial models, selected by the type parameters.
Uses Boyle-Lau adjustment for optimize steps and Derman-Kani optimization to speed
up convergence.
Type values:
    crr or coxrossrubinstein:        Cox-Ross-Rubinstein model
    jr  or jarrowrudd:               Jarrow-Rudd model
    eqp or additiveeqpbinomialtree:  Additive EQP model
    trigeorgis:                      Trigeorgis model
    tian:                            Tian model
    lr  or leisenreimer              Leisen-Reimer model
    j4  or joshi4:                   Joshi 4th (smoothed) model

Boyle-Lau adjustment is controlled by parameter max_steps.
If max_steps is equal to steps Boyle-Lau is disabled.
Il max_steps is 0 (default value), max_steps is calculated by capping it to 
5*steps when Boyle-Lau would need more than 1000 steps.
If max_steps is specified, it would limit binomial steps to this value.
"
#endif
%rename(BinomialBarrierEngine) BinomialBarrierEnginePtr;
class BinomialBarrierEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        BinomialBarrierEnginePtr(
                             const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                             const std::string& type,
                             Size steps,
                             Size max_steps = 0) {
            std::string s = boost::algorithm::to_lower_copy(type);
            if (s == "crr" || s == "coxrossrubinstein")
                return new BinomialBarrierEnginePtr(
                    new BinomialBarrierEngine<CoxRossRubinstein,
                                          DiscretizedDermanKaniBarrierOption>(
                                                  process,steps,max_steps));
            else if (s == "jr" || s == "jarrowrudd")
                return new BinomialBarrierEnginePtr(
                    new BinomialBarrierEngine<JarrowRudd,
                                          DiscretizedDermanKaniBarrierOption>(
                                                  process,steps,max_steps));
            else if (s == "eqp" || s == "additiveeqpbinomialtree")
                return new BinomialBarrierEnginePtr(
                    new BinomialBarrierEngine<AdditiveEQPBinomialTree,
                                          DiscretizedDermanKaniBarrierOption>(
                                                  process,steps,max_steps));
            else if (s == "trigeorgis")
                return new BinomialBarrierEnginePtr(
                    new BinomialBarrierEngine<Trigeorgis,
                                          DiscretizedDermanKaniBarrierOption>(
                                                  process,steps,max_steps));
            else if (s == "tian")
                return new BinomialBarrierEnginePtr(
                    new BinomialBarrierEngine<Tian,
                                          DiscretizedDermanKaniBarrierOption>(
                                                  process,steps,max_steps));
            else if (s == "lr" || s == "leisenreimer")
                return new BinomialBarrierEnginePtr(
                    new BinomialBarrierEngine<LeisenReimer,
                                          DiscretizedDermanKaniBarrierOption>(
                                                  process,steps,max_steps));
            else if (s == "j4" || s == "joshi4")
                return new BinomialBarrierEnginePtr(
                    new BinomialBarrierEngine<Joshi4,
                                          DiscretizedDermanKaniBarrierOption>(
                                                  process,steps,max_steps));
            else
                QL_FAIL("unknown binomial barrier engine type: "+s);
        }
    }
};

%{
using QuantLib::QuantoEngine;
using QuantLib::ForwardVanillaEngine;
typedef boost::shared_ptr<PricingEngine> ForwardEuropeanEnginePtr;
typedef boost::shared_ptr<PricingEngine> QuantoEuropeanEnginePtr;
typedef boost::shared_ptr<PricingEngine> QuantoForwardEuropeanEnginePtr;
%}


%rename(ForwardEuropeanEngine) ForwardEuropeanEnginePtr;
class ForwardEuropeanEnginePtr: public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        ForwardEuropeanEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new ForwardEuropeanEnginePtr(
                 new ForwardVanillaEngine<AnalyticEuropeanEngine>(process));
        }
    }
};


%rename(QuantoEuropeanEngine) QuantoEuropeanEnginePtr;
class QuantoEuropeanEnginePtr: public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        QuantoEuropeanEnginePtr(
                  const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                  const Handle<YieldTermStructure>& foreignRiskFreeRate,
                  const Handle<BlackVolTermStructure>& exchangeRateVolatility,
                  const Handle<Quote>& correlation) {
            return new QuantoEuropeanEnginePtr(
                new QuantoEngine<VanillaOption,AnalyticEuropeanEngine>(
                                                       process,
                                                       foreignRiskFreeRate,
                                                       exchangeRateVolatility,
                                                       correlation));
        }
    }
};

%rename(QuantoForwardEuropeanEngine) QuantoForwardEuropeanEnginePtr;
class QuantoForwardEuropeanEnginePtr: public boost::shared_ptr<PricingEngine> {
public:
    %extend {
        QuantoForwardEuropeanEnginePtr(
                  const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                  const Handle<YieldTermStructure>& foreignRiskFreeRate,
                  const Handle<BlackVolTermStructure>& exchangeRateVolatility,
                  const Handle<Quote>& correlation) {
            return new QuantoForwardEuropeanEnginePtr(
                new QuantoEngine<ForwardVanillaOption,AnalyticEuropeanEngine>(
                                                       process,
                                                       foreignRiskFreeRate,
                                                       exchangeRateVolatility,
                                                       correlation));
        }
    }
};

%{
using QuantLib::BlackCalculator;
%}

class BlackCalculator {
  public:
    %extend {
        BlackCalculator (
                   const boost::shared_ptr<Payoff>& payoff,
                   Real forward,
                   Real stdDev,
                   Real discount = 1.0) {

            boost::shared_ptr<StrikedTypePayoff> stPayoff =
                boost::dynamic_pointer_cast<StrikedTypePayoff>(payoff);

            QL_REQUIRE(stPayoff, "wrong payoff given");

            return new BlackCalculator(stPayoff,forward,stdDev,discount);
        }
    }
    Real value() const;
    Real deltaForward() const;
    Real delta(Real spot) const;
    Real elasticityForward() const;
    Real elasticity(Real spot) const;
    Real gammaForward() const;
    Real gamma(Real spot) const;
    Real theta(Real spot, Time maturity) const;
    Real thetaPerDay(Real spot, Time maturity) const;
    Real vega(Time maturity) const;
    Real rho(Time maturity) const;
    Real dividendRho(Time maturity) const;
    Real itmCashProbability() const;
    Real itmAssetProbability() const;
    Real strikeSensitivity() const;
    Real alpha() const;
    Real beta() const;
};



// Asian options

%{
using QuantLib::Average;
using QuantLib::ContinuousAveragingAsianOption;
using QuantLib::DiscreteAveragingAsianOption;
%}

struct Average {
    enum Type { Arithmetic, Geometric };
};

%shared_ptr(ContinuousAveragingAsianOption)
class ContinuousAveragingAsianOption : public OneAssetOption {
  public:
    ContinuousAveragingAsianOption(
            Average::Type averageType,
            const boost::shared_ptr<StrikedTypePayoff>& payoff,
            const boost::shared_ptr<Exercise>& exercise);
};

%shared_ptr(DiscreteAveragingAsianOption)
class DiscreteAveragingAsianOption : public OneAssetOption {
  public:
    DiscreteAveragingAsianOption(
            Average::Type averageType,
            Real runningAccumulator,
            Size pastFixings,
            const std::vector<Date>& fixingDates,
            const boost::shared_ptr<StrikedTypePayoff>& payoff,
            const boost::shared_ptr<Exercise>& exercise);
};

// Asian engines


%{
using QuantLib::AnalyticContinuousGeometricAveragePriceAsianEngine;
typedef boost::shared_ptr<PricingEngine>
    AnalyticContinuousGeometricAveragePriceAsianEnginePtr;
%}

%rename(AnalyticContinuousGeometricAveragePriceAsianEngine)
        AnalyticContinuousGeometricAveragePriceAsianEnginePtr;
class AnalyticContinuousGeometricAveragePriceAsianEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticContinuousGeometricAveragePriceAsianEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new AnalyticContinuousGeometricAveragePriceAsianEnginePtr(
                new AnalyticContinuousGeometricAveragePriceAsianEngine(
                                                                  process));
        }
    }
};


%{
using QuantLib::AnalyticDiscreteGeometricAveragePriceAsianEngine;
typedef boost::shared_ptr<PricingEngine>
    AnalyticDiscreteGeometricAveragePriceAsianEnginePtr;
%}

%rename(AnalyticDiscreteGeometricAveragePriceAsianEngine)
        AnalyticDiscreteGeometricAveragePriceAsianEnginePtr;
class AnalyticDiscreteGeometricAveragePriceAsianEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticDiscreteGeometricAveragePriceAsianEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new AnalyticDiscreteGeometricAveragePriceAsianEnginePtr(
                new AnalyticDiscreteGeometricAveragePriceAsianEngine(
                                                                  process));
        }
    }
};


%{
using QuantLib::AnalyticDiscreteGeometricAverageStrikeAsianEngine;
typedef boost::shared_ptr<PricingEngine>
    AnalyticDiscreteGeometricAverageStrikeAsianEnginePtr;
%}

%rename(AnalyticDiscreteGeometricAverageStrikeAsianEngine)
        AnalyticDiscreteGeometricAverageStrikeAsianEnginePtr;
class AnalyticDiscreteGeometricAverageStrikeAsianEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticDiscreteGeometricAverageStrikeAsianEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {
            return new AnalyticDiscreteGeometricAverageStrikeAsianEnginePtr(
                new AnalyticDiscreteGeometricAverageStrikeAsianEngine(
                                                                  process));
        }
    }
};



%{
using QuantLib::MCDiscreteArithmeticAPEngine;
typedef boost::shared_ptr<PricingEngine> MCDiscreteArithmeticAPEnginePtr;
%}

%rename(MCDiscreteArithmeticAPEngine) MCDiscreteArithmeticAPEnginePtr;
class MCDiscreteArithmeticAPEnginePtr
    : public boost::shared_ptr<PricingEngine> {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCDiscreteArithmeticAPEnginePtr;
    #endif
  public:
    %extend {
        MCDiscreteArithmeticAPEnginePtr(
                            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                            const std::string& traits,
                            bool brownianBridge = false,
                            bool antitheticVariate = false,
                            bool controlVariate = false,
                            intOrNull requiredSamples = Null<Size>(),
                            doubleOrNull requiredTolerance = Null<Real>(),
                            intOrNull maxSamples = Null<Size>(),
                            BigInteger seed = 0) {

                 


            std::string s = boost::algorithm::to_lower_copy(traits);
            if (s == "pseudorandom" || s == "pr")
                return new MCDiscreteArithmeticAPEnginePtr(
                         new MCDiscreteArithmeticAPEngine<PseudoRandom>(
                                                            process,
                                                            brownianBridge,
                                                            antitheticVariate,
                                                            controlVariate,
                                                            requiredSamples,
                                                            requiredTolerance,
                                                            maxSamples,
                                                            seed));
            else if (s == "lowdiscrepancy" || s == "ld")
                return new MCDiscreteArithmeticAPEnginePtr(
                       new MCDiscreteArithmeticAPEngine<LowDiscrepancy>(
                                                            process,
                                                            brownianBridge,
                                                            antitheticVariate,
                                                            controlVariate,
                                                            requiredSamples,
                                                            requiredTolerance,
                                                            maxSamples,
                                                            seed));
            else
                QL_FAIL("unknown Monte Carlo engine type: "+s);
        }
    }
};


%{
using QuantLib::MCDiscreteArithmeticASEngine;
typedef boost::shared_ptr<PricingEngine> MCDiscreteArithmeticASEnginePtr;
%}

%rename(MCDiscreteArithmeticASEngine) MCDiscreteArithmeticASEnginePtr;
class MCDiscreteArithmeticASEnginePtr
    : public boost::shared_ptr<PricingEngine> {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCDiscreteArithmeticASEnginePtr;
    #endif
  public:
    %extend {
        MCDiscreteArithmeticASEnginePtr(
                            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                            const std::string& traits,
                            bool brownianBridge = false,
                            bool antitheticVariate = false,
                            intOrNull requiredSamples = Null<Size>(),
                            doubleOrNull requiredTolerance = Null<Real>(),
                            intOrNull maxSamples = Null<Size>(),
                            BigInteger seed = 0) {

                 


            std::string s = boost::algorithm::to_lower_copy(traits);
            if (s == "pseudorandom" || s == "pr")
                return new MCDiscreteArithmeticASEnginePtr(
                         new MCDiscreteArithmeticASEngine<PseudoRandom>(
                                                            process,
                                                            brownianBridge,
                                                            antitheticVariate,
                                                            requiredSamples,
                                                            requiredTolerance,
                                                            maxSamples,
                                                            seed));
            else if (s == "lowdiscrepancy" || s == "ld")
                return new MCDiscreteArithmeticASEnginePtr(
                       new MCDiscreteArithmeticASEngine<LowDiscrepancy>(
                                                            process,
                                                            brownianBridge,
                                                            antitheticVariate,
                                                            requiredSamples,
                                                            requiredTolerance,
                                                            maxSamples,
                                                            seed));
            else
                QL_FAIL("unknown Monte Carlo engine type: "+s);
        }
    }
};


%{
using QuantLib::MCDiscreteGeometricAPEngine;
typedef boost::shared_ptr<PricingEngine> MCDiscreteGeometricAPEnginePtr;
%}

%rename(MCDiscreteGeometricAPEngine) MCDiscreteGeometricAPEnginePtr;
class MCDiscreteGeometricAPEnginePtr
    : public boost::shared_ptr<PricingEngine> {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCDiscreteGeometricAPEnginePtr;
    #endif
  public:
    %extend {
        MCDiscreteGeometricAPEnginePtr(
                            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                            const std::string& traits,
                            bool brownianBridge = false,
                            bool antitheticVariate = false,
                            intOrNull requiredSamples = Null<Size>(),
                            doubleOrNull requiredTolerance = Null<Real>(),
                            intOrNull maxSamples = Null<Size>(),
                            BigInteger seed = 0) {

                 


            std::string s = boost::algorithm::to_lower_copy(traits);
            if (s == "pseudorandom" || s == "pr")
                return new MCDiscreteGeometricAPEnginePtr(
                         new MCDiscreteGeometricAPEngine<PseudoRandom>(
                                                            process,
                                                            brownianBridge,
                                                            antitheticVariate,
                                                            requiredSamples,
                                                            requiredTolerance,
                                                            maxSamples,
                                                            seed));
            else if (s == "lowdiscrepancy" || s == "ld")
                return new MCDiscreteGeometricAPEnginePtr(
                       new MCDiscreteGeometricAPEngine<LowDiscrepancy>(
                                                            process,
                                                            brownianBridge,
                                                            antitheticVariate,
                                                            requiredSamples,
                                                            requiredTolerance,
                                                            maxSamples,
                                                            seed));
            else
                QL_FAIL("unknown Monte Carlo engine type: "+s);
        }
    }
};

%{
using QuantLib::VarianceGammaEngine;
typedef boost::shared_ptr<PricingEngine>
    VarianceGammaEnginePtr;
%}

%rename(VarianceGammaEngine) VarianceGammaEnginePtr;
class VarianceGammaEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        VarianceGammaEnginePtr(const boost::shared_ptr<VarianceGammaProcess>& process) {
            return new VarianceGammaEnginePtr(new VarianceGammaEngine(process));
        }
    }
};

%{
using QuantLib::FFTVarianceGammaEngine;
typedef boost::shared_ptr<PricingEngine>
    FFTVarianceGammaEnginePtr;
%}

%rename(FFTVarianceGammaEngine) FFTVarianceGammaEnginePtr;
class FFTVarianceGammaEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FFTVarianceGammaEnginePtr(const boost::shared_ptr<VarianceGammaProcess>& process, Real logStrikeSpacing = 0.001) {
            return new FFTVarianceGammaEnginePtr(new FFTVarianceGammaEngine(process, logStrikeSpacing));
        }
        void precalculate(const std::vector<boost::shared_ptr<Instrument> >& optionList)
        {
            boost::dynamic_pointer_cast<FFTVarianceGammaEngine>(*self)->precalculate(optionList);
        }
    }
};

// Double barrier options
%{
using QuantLib::DoubleBarrierOption;
using QuantLib::DoubleBarrier;
%}

%{
using QuantLib::DoubleBarrierOption;
%}

%shared_ptr(DoubleBarrierOption)
class DoubleBarrierOption : public OneAssetOption {
  public:
    DoubleBarrierOption(
               DoubleBarrier::Type barrierType,
               Real barrier_lo,
               Real barrier_hi,
               Real rebate,
               const boost::shared_ptr<StrikedTypePayoff>& payoff,
               const boost::shared_ptr<Exercise>& exercise);
};

// QuantoVanillaOption

%{
using QuantLib::QuantoDoubleBarrierOption;
%}

%shared_ptr(QuantoDoubleBarrierOption)
class QuantoDoubleBarrierOption : public DoubleBarrierOption {
  public:
    QuantoDoubleBarrierOption(
            DoubleBarrier::Type barrierType,
            Real barrier_lo,
            Real barrier_hi,
            Real rebate,
            const boost::shared_ptr<StrikedTypePayoff>& payoff,
            const boost::shared_ptr<Exercise>& exercise);
    Real qvega();
    Real qrho();
    Real qlambda();
};

// Double Barrier engines

%{
using QuantLib::AnalyticDoubleBarrierEngine;
typedef boost::shared_ptr<PricingEngine> AnalyticDoubleBarrierEnginePtr;
%}

#if defined(SWIGPYTHON)
%feature("docstring") AnalyticDoubleBarrierEnginePtr "Double barrier engine implementing Ikeda-Kunitomo series."
#endif
%rename(AnalyticDoubleBarrierEngine) AnalyticDoubleBarrierEnginePtr;
class AnalyticDoubleBarrierEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticDoubleBarrierEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                           int series = 5) {

                 


            return new AnalyticDoubleBarrierEnginePtr(
                            new AnalyticDoubleBarrierEngine(process, series));
        }
    }
};

%{
using QuantLib::WulinYongDoubleBarrierEngine;
typedef boost::shared_ptr<PricingEngine> WulinYongDoubleBarrierEnginePtr;
%}

%rename(WulinYongDoubleBarrierEngine) WulinYongDoubleBarrierEnginePtr;
class WulinYongDoubleBarrierEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        WulinYongDoubleBarrierEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                           int series = 5) {

                 


            return new WulinYongDoubleBarrierEnginePtr(
                            new WulinYongDoubleBarrierEngine(process, series));
        }
    }
};

%{
using QuantLib::VannaVolgaDoubleBarrierEngine;
using QuantLib::VannaVolgaBarrierEngine;
using QuantLib::DeltaVolQuote;
typedef boost::shared_ptr<PricingEngine> VannaVolgaDoubleBarrierEnginePtr;
typedef boost::shared_ptr<PricingEngine> VannaVolgaBarrierEnginePtr;
%}

%shared_ptr(DeltaVolQuote)

class DeltaVolQuote : public Quote {
  public:
    enum DeltaType { Spot, Fwd, PaSpot, PaFwd };
    enum AtmType { AtmNull, AtmSpot, AtmFwd, AtmDeltaNeutral,
                   AtmVegaMax, AtmGammaMax, AtmPutCall50 };
    DeltaVolQuote(Real delta,
                  const Handle<Quote>& vol,
                  Time maturity,
                  DeltaVolQuote::DeltaType deltaType);
    DeltaVolQuote(const Handle<Quote>& vol,
                  DeltaVolQuote::DeltaType deltaType,
                  Time maturity,
                  DeltaVolQuote::AtmType atmType);
};

%template(DeltaVolQuoteHandle) Handle<DeltaVolQuote>;
%template(RelinkableDeltaVolQuoteHandle)
RelinkableHandle<DeltaVolQuote>;

#if defined(SWIGPYTHON)
%feature("docstring") VannaVolgaDoubleBarrierEnginePtr "
Vanna-Volga engine for double barrier options.
Supports different double barrier engines, selected by the type parameters.
Type values:
    ik or analytic:  Ikeda-Kunitomo standard engine (default)
    wo:              Wulin-Yong engine
"
#endif
%rename(VannaVolgaDoubleBarrierEngine) VannaVolgaDoubleBarrierEnginePtr;
class VannaVolgaDoubleBarrierEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        VannaVolgaDoubleBarrierEnginePtr(
                           const Handle<DeltaVolQuote> atmVol,
                           const Handle<DeltaVolQuote> vol25Put,
                           const Handle<DeltaVolQuote> vol25Call,
                           const Handle<Quote> spotFX,
                           const Handle<YieldTermStructure> domesticTS,
                           const Handle<YieldTermStructure> foreignTS,
                           const std::string& type = "ik",
                           const bool adaptVanDelta = false,
                           const Real bsPriceWithSmile = 0.0,
                           int series = 5) {
            std::string s = boost::algorithm::to_lower_copy(type);
            if (s == "ik" || s == "analytic")
                return new VannaVolgaDoubleBarrierEnginePtr(
                   new VannaVolgaDoubleBarrierEngine<AnalyticDoubleBarrierEngine>(
                                        atmVol, vol25Put, vol25Call, spotFX, 
                                        domesticTS, foreignTS, adaptVanDelta, 
                                        bsPriceWithSmile, series));
            else if (s == "wo")
                return new VannaVolgaDoubleBarrierEnginePtr(
                   new VannaVolgaDoubleBarrierEngine<WulinYongDoubleBarrierEngine>(
                                        atmVol, vol25Put, vol25Call, spotFX, 
                                        domesticTS, foreignTS, adaptVanDelta, 
                                        bsPriceWithSmile, series));
            else
                QL_FAIL("unknown binomial engine type: "+s);
        }
    }
};

%rename(VannaVolgaBarrierEngine) VannaVolgaBarrierEnginePtr;
class VannaVolgaBarrierEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        VannaVolgaBarrierEnginePtr(
                const Handle<DeltaVolQuote>& atmVol,
                const Handle<DeltaVolQuote>& vol25Put,
                const Handle<DeltaVolQuote>& vol25Call,
                const Handle<Quote>& spotFX,
                const Handle<YieldTermStructure>& domesticTS,
                const Handle<YieldTermStructure>& foreignTS,
                const bool adaptVanDelta = false,
                const Real bsPriceWithSmile = 0.0) {
                return new VannaVolgaBarrierEnginePtr(
                   new VannaVolgaBarrierEngine(
                                        atmVol, vol25Put, vol25Call, spotFX, 
                                        domesticTS, foreignTS, adaptVanDelta, 
                                        bsPriceWithSmile));
        }
    }
};

%{
using QuantLib::AnalyticDoubleBarrierBinaryEngine;
typedef boost::shared_ptr<PricingEngine> AnalyticDoubleBarrierBinaryEnginePtr;
%}

%rename(AnalyticDoubleBarrierBinaryEngine) AnalyticDoubleBarrierBinaryEnginePtr;
class AnalyticDoubleBarrierBinaryEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        AnalyticDoubleBarrierBinaryEnginePtr(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process) {

                 


            return new AnalyticDoubleBarrierBinaryEnginePtr(
                            new AnalyticDoubleBarrierBinaryEngine(process));
        }
    }
};

%{
using QuantLib::BinomialDoubleBarrierEngine;
using QuantLib::DiscretizedDermanKaniDoubleBarrierOption;
typedef boost::shared_ptr<PricingEngine> BinomialDoubleBarrierEnginePtr;
%}

#if defined(SWIGPYTHON)
%feature("docstring") BinomialDoubleBarrierEnginePtr "Binomial Engine for double barrier options.
Features different binomial models, selected by the type parameters.
Uses Derman-Kani optimization to speed up convergence.
Type values:
    crr or coxrossrubinstein:        Cox-Ross-Rubinstein model
    jr  or jarrowrudd:               Jarrow-Rudd model
    eqp or additiveeqpbinomialtree:  Additive EQP model
    trigeorgis:                      Trigeorgis model
    tian:                            Tian model
    lr  or leisenreimer              Leisen-Reimer model
    j4  or joshi4:                   Joshi 4th (smoothed) model
"
#endif
%rename(BinomialDoubleBarrierEngine) BinomialDoubleBarrierEnginePtr;
class BinomialDoubleBarrierEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        BinomialDoubleBarrierEnginePtr(
                             const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                             const std::string& type,
                             Size steps) {

                 


            std::string s = boost::algorithm::to_lower_copy(type);
            if (s == "crr" || s == "coxrossrubinstein")
                return new BinomialDoubleBarrierEnginePtr(
                    new BinomialDoubleBarrierEngine<CoxRossRubinstein,
                                     DiscretizedDermanKaniDoubleBarrierOption>(
                                                  process,steps));
            else if (s == "jr" || s == "jarrowrudd")
                return new BinomialDoubleBarrierEnginePtr(
                    new BinomialDoubleBarrierEngine<JarrowRudd,
                                     DiscretizedDermanKaniDoubleBarrierOption>(
                                                  process,steps));
            else if (s == "eqp" || s == "additiveeqpbinomialtree")
                return new BinomialDoubleBarrierEnginePtr(
                    new BinomialDoubleBarrierEngine<AdditiveEQPBinomialTree,
                                     DiscretizedDermanKaniDoubleBarrierOption>(
                                                  process,steps));
            else if (s == "trigeorgis")
                return new BinomialDoubleBarrierEnginePtr(
                    new BinomialDoubleBarrierEngine<Trigeorgis,
                                     DiscretizedDermanKaniDoubleBarrierOption>(
                                                  process,steps));
            else if (s == "tian")
                return new BinomialDoubleBarrierEnginePtr(
                    new BinomialDoubleBarrierEngine<Tian,
                                     DiscretizedDermanKaniDoubleBarrierOption>(
                                                  process,steps));
            else if (s == "lr" || s == "leisenreimer")
                return new BinomialDoubleBarrierEnginePtr(
                    new BinomialDoubleBarrierEngine<LeisenReimer,
                                     DiscretizedDermanKaniDoubleBarrierOption>(
                                                  process,steps));
            else if (s == "j4" || s == "joshi4")
                return new BinomialDoubleBarrierEnginePtr(
                    new BinomialDoubleBarrierEngine<Joshi4,
                                     DiscretizedDermanKaniDoubleBarrierOption>(
                                                  process,steps));
            else
                QL_FAIL("unknown binomial double barrier engine type: "+s);
        }
    }
};


// Swing option

%{
using QuantLib::VanillaSwingOption;
%}

%shared_ptr(VanillaSwingOption)
class VanillaSwingOption : public OneAssetOption {
  public:
    VanillaSwingOption(
        const boost::shared_ptr<Payoff>& payoff,
        const boost::shared_ptr<SwingExercise>& ex,
        Size minExerciseRights, Size maxExerciseRights);
};

// Swing engines

%{
using QuantLib::FdSimpleBSSwingEngine;
using QuantLib::FdSimpleExtOUJumpSwingEngine;
typedef boost::shared_ptr<PricingEngine> FdSimpleExtOUJumpSwingEnginePtr;
typedef boost::shared_ptr<PricingEngine> FdSimpleBSSwingEnginePtr;
%}

%rename(FdSimpleBSSwingEngine) FdSimpleBSSwingEnginePtr;
class FdSimpleBSSwingEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdSimpleBSSwingEnginePtr(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
            Size tGrid = 50, Size xGrid = 100,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas()) {

                 


            
            return new FdSimpleBSSwingEnginePtr(
                new FdSimpleBSSwingEngine(process,tGrid, xGrid, schemeDesc));
        }
    }
};

%rename(FdSimpleExtOUJumpSwingEngine) FdSimpleExtOUJumpSwingEnginePtr;
class FdSimpleExtOUJumpSwingEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdSimpleExtOUJumpSwingEnginePtr(            
        	const boost::shared_ptr<ExtOUWithJumpsProcess>& process,
            const boost::shared_ptr<YieldTermStructure>& rTS,
            Size tGrid = 50, Size xGrid = 200, Size yGrid=50,
            const std::vector<std::pair<Time,Real> >& shape =
                                         std::vector<std::pair<Time,Real> >(),
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer()) {

            boost::shared_ptr<FdSimpleExtOUJumpSwingEngine::Shape> curve(
                              new FdSimpleExtOUJumpSwingEngine::Shape(shape));
            
            return new FdSimpleExtOUJumpSwingEnginePtr(
                new FdSimpleExtOUJumpSwingEngine(
                    process, rTS, tGrid, xGrid, yGrid, 
                    curve, schemeDesc));
        }
    }
};


#endif
