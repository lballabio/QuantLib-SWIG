/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2008 Tito Ingargiola
 Copyright (C) 2010, 2012, 2018, 2019 Klaus Spanderen
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
  private:
    OneAssetOption();
  public:
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
%}

%shared_ptr(AnalyticEuropeanEngine)
class AnalyticEuropeanEngine : public PricingEngine {
  public:
    AnalyticEuropeanEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>&);
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
%}

%shared_ptr(AnalyticHestonEngine)
class AnalyticHestonEngine : public PricingEngine {
  public:
    AnalyticHestonEngine(const boost::shared_ptr<HestonModel>& model,
                         Size integrationOrder = 144);
    AnalyticHestonEngine(const boost::shared_ptr<HestonModel>& model,
                         Real relTolerance,
                         Size maxEvaluations);
};

%{
using QuantLib::COSHestonEngine;
%}

%shared_ptr(COSHestonEngine)
class COSHestonEngine : public PricingEngine {
  public:
    COSHestonEngine(const boost::shared_ptr<HestonModel>& model,
                    Real L = 16, Size N = 200);
};


%{
using QuantLib::AnalyticPTDHestonEngine;
%}

%shared_ptr(AnalyticPTDHestonEngine)
class AnalyticPTDHestonEngine : public PricingEngine {
  public:
    AnalyticPTDHestonEngine(
            const boost::shared_ptr<PiecewiseTimeDependentHestonModel>& model,
            Real relTolerance, Size maxEvaluations);
    // Constructor using Laguerre integration
    // and Gatheral's version of complex log.
    AnalyticPTDHestonEngine(
            const boost::shared_ptr<PiecewiseTimeDependentHestonModel>& model,
            Size integrationOrder = 144);
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
%}

%shared_ptr(BatesEngine)
class BatesEngine : public PricingEngine {
  public:
    BatesEngine(const boost::shared_ptr<BatesModel>& model,
                Size integrationOrder = 144);
    BatesEngine(const boost::shared_ptr<BatesModel>& model,
                Real relTolerance,
                Size maxEvaluations);
};


%{
using QuantLib::IntegralEngine;
%}

%shared_ptr(IntegralEngine)
class IntegralEngine : public PricingEngine {
  public:
    IntegralEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>&);
};


%{
using QuantLib::CrankNicolson;
typedef QuantLib::FDBermudanEngine<CrankNicolson> FDBermudanEngine;
%}

%shared_ptr(FDBermudanEngine)
class FDBermudanEngine : public PricingEngine {
  public:
    FDBermudanEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                     Size timeSteps = 100, Size gridPoints = 100,
                     bool timeDependent = false);
};

%{
typedef QuantLib::FDEuropeanEngine<CrankNicolson> FDEuropeanEngine;
%}

%shared_ptr(FDEuropeanEngine)
class FDEuropeanEngine : public PricingEngine {
  public:
    FDEuropeanEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess> process,
                     Size timeSteps = 100, Size gridPoints = 100,
                     bool timeDependent = false);
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
%}

%shared_ptr(BinomialVanillaEngine<CoxRossRubinstein>)
%shared_ptr(BinomialVanillaEngine<JarrowRudd>)
%shared_ptr(BinomialVanillaEngine<AdditiveEQPBinomialTree>)
%shared_ptr(BinomialVanillaEngine<Trigeorgis>)
%shared_ptr(BinomialVanillaEngine<Tian>)
%shared_ptr(BinomialVanillaEngine<LeisenReimer>)
%shared_ptr(BinomialVanillaEngine<Joshi4>)

template <class T>
class BinomialVanillaEngine : public PricingEngine {
  public:
    BinomialVanillaEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>&,
                          Size steps);
};

%template(BinomialCRRVanillaEngine) BinomialVanillaEngine<CoxRossRubinstein>;
%template(BinomialJRVanillaEngine) BinomialVanillaEngine<JarrowRudd>;
%template(BinomialEQPVanillaEngine) BinomialVanillaEngine<AdditiveEQPBinomialTree>;
%template(BinomialTrigeorgisVanillaEngine) BinomialVanillaEngine<Trigeorgis>;
%template(BinomialTianVanillaEngine) BinomialVanillaEngine<Tian>;
%template(BinomialLRVanillaEngine) BinomialVanillaEngine<LeisenReimer>;
%template(BinomialJ4VanillaEngine) BinomialVanillaEngine<Joshi4>;


#if defined(SWIGPYTHON)
%pythoncode %{
    def BinomialVanillaEngine(process, type, steps):
        type = type.lower()
        if type == "crr" or type == "coxrossrubinstein":
            cls = BinomialCRRVanillaEngine
        elif type == "jr" or type == "jarrowrudd":
            cls = BinomialJRVanillaEngine
        elif type == "eqp":
            cls = BinomialEQPVanillaEngine
        elif type == "trigeorgis":
            cls = BinomialTrigeorgisVanillaEngine
        elif type == "tian":
            cls = BinomialTianVanillaEngine
        elif type == "lr" or type == "leisenreimer":
            cls = BinomialLRVanillaEngine
        elif type == "j4" or type == "joshi4":
            cls = BinomialJ4VanillaEngine
        else:
            raise RuntimeError("unknown binomial engine type: %s" % type);
        return cls(process, steps)
%}
#endif

%{
using QuantLib::MCEuropeanEngine;
using QuantLib::MCEuropeanHestonEngine;
using QuantLib::MCAmericanEngine;
using QuantLib::PseudoRandom;
using QuantLib::LowDiscrepancy;
using QuantLib::LsmBasisSystem;
%}

struct LsmBasisSystem {
    enum PolynomType  {Monomial, Laguerre, Hermite, Hyperbolic,
                           Legendre, Chebyshev, Chebyshev2nd };
};

%shared_ptr(MCEuropeanEngine<PseudoRandom>);
%shared_ptr(MCEuropeanEngine<LowDiscrepancy>);

template <class RNG>
class MCEuropeanEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCEuropeanEngine;
    #endif
  public:
    %extend {
        MCEuropeanEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                         intOrNull timeSteps = Null<Size>(),
                         intOrNull timeStepsPerYear = Null<Size>(),
                         bool brownianBridge = false,
                         bool antitheticVariate = false,
                         intOrNull requiredSamples = Null<Size>(),
                         doubleOrNull requiredTolerance = Null<Real>(),
                         intOrNull maxSamples = Null<Size>(),
                         BigInteger seed = 0) {
            QL_REQUIRE(Size(timeSteps) != Null<Size>() ||
                       Size(timeStepsPerYear) != Null<Size>(),
                       "number of steps not specified");
            return new MCEuropeanEngine<RNG>(process,
                                             timeSteps,
                                             timeStepsPerYear,
                                             brownianBridge,
                                             antitheticVariate,
                                             requiredSamples,
                                             requiredTolerance,
                                             maxSamples,
                                             seed);
        }
    }
};

%template(MCPREuropeanEngine) MCEuropeanEngine<PseudoRandom>;
%template(MCLDEuropeanEngine) MCEuropeanEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCEuropeanEngine(process,
                         traits,
                         timeSteps=None,
                         timeStepsPerYear=None,
                         brownianBridge=False,
                         antitheticVariate=False,
                         requiredSamples=None,
                         requiredTolerance=None,
                         maxSamples=None,
                         seed=0):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPREuropeanEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDEuropeanEngine
        else:
            raise RuntimeError("unknown MC traits: %s" % traits);
        return cls(process,
                   timeSteps,
                   timeStepsPerYear,
                   brownianBridge,
                   antitheticVariate,
                   requiredSamples,
                   requiredTolerance,
                   maxSamples,
                   seed)
%}
#endif


%shared_ptr(MCAmericanEngine<PseudoRandom>);
%shared_ptr(MCAmericanEngine<LowDiscrepancy>);

template <class RNG>
class MCAmericanEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCAmericanEngine;
    #endif
  public:
    %extend {
        MCAmericanEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
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
            return new MCAmericanEngine<RNG>(process,
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
                                             seedCalibration);
        }
    }
};

%template(MCPRAmericanEngine) MCAmericanEngine<PseudoRandom>;
%template(MCLDAmericanEngine) MCAmericanEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCAmericanEngine(process,
                         traits,
                         timeSteps=None,
                         timeStepsPerYear=None,
                         antitheticVariate=False,
                         controlVariate=False,
                         requiredSamples=None,
                         requiredTolerance=None,
                         maxSamples=None,
                         seed=0,
                         polynomOrder=2,
                         polynomType=LsmBasisSystem.Monomial,
                         nCalibrationSamples=2048,
                         antitheticVariateCalibration=None,
                         seedCalibration=None):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPRAmericanEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDAmericanEngine
        else:
            raise RuntimeError("unknown MC traits: %s" % traits);
        return cls(process,
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
                   seedCalibration)
%}
#endif


%shared_ptr(MCEuropeanHestonEngine<PseudoRandom>);
%shared_ptr(MCEuropeanHestonEngine<LowDiscrepancy>);

template <class RNG>
class MCEuropeanHestonEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCEuropeanHestonEngine;
    #endif
  public:
    %extend {
        MCEuropeanHestonEngine(const boost::shared_ptr<HestonProcess>& process,
                               intOrNull timeSteps = Null<Size>(),
                               intOrNull timeStepsPerYear = Null<Size>(),
                               bool antitheticVariate = false,
                               intOrNull requiredSamples = Null<Size>(),
                               doubleOrNull requiredTolerance = Null<Real>(),
                               intOrNull maxSamples = Null<Size>(),
                               BigInteger seed = 0) {
            QL_REQUIRE(Size(timeSteps) != Null<Size>() ||
                       Size(timeStepsPerYear) != Null<Size>(),
                       "number of steps not specified");
            return new MCEuropeanHestonEngine<RNG>(process,
                                                   timeSteps,
                                                   timeStepsPerYear,
                                                   antitheticVariate,
                                                   requiredSamples,
                                                   requiredTolerance,
                                                   maxSamples,
                                                   seed);
        }
    }
};

%template(MCPREuropeanHestonEngine) MCEuropeanHestonEngine<PseudoRandom>;
%template(MCLDEuropeanHestonEngine) MCEuropeanHestonEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCEuropeanHestonEngine(process,
                               traits,
                               timeSteps=None,
                               timeStepsPerYear=None,
                               antitheticVariate=False,
                               requiredSamples=None,
                               requiredTolerance=None,
                               maxSamples=None,
                               seed=0):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPREuropeanHestonEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDEuropeanHestonEngine
        else:
            raise RuntimeError("unknown MC traits: %s" % traits);
        return cls(process,
                   timeSteps,
                   timeStepsPerYear,
                   antitheticVariate,
                   requiredSamples,
                   requiredTolerance,
                   maxSamples,
                   seed)
%}
#endif


// American engines

%{
typedef QuantLib::FDAmericanEngine<CrankNicolson> FDAmericanEngine;
typedef QuantLib::FDShoutEngine<CrankNicolson> FDShoutEngine;
%}

%shared_ptr(FDAmericanEngine)
class FDAmericanEngine : public PricingEngine {
  public:
    FDAmericanEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                     Size timeSteps = 100, Size gridPoints = 100,
                     bool timeDependent = false);
};

%shared_ptr(FDShoutEngine)
class FDShoutEngine : public PricingEngine {
  public:
    FDShoutEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                  Size timeSteps = 100, Size gridPoints = 100,
                  bool timeDependent = false);
};


%{
using QuantLib::ContinuousArithmeticAsianLevyEngine;
%}

%shared_ptr(ContinuousArithmeticAsianLevyEngine)
class ContinuousArithmeticAsianLevyEngine : public PricingEngine {
  public:
    ContinuousArithmeticAsianLevyEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                        const Handle<Quote>& runningAverage,
                                        const Date& startDate);
};

%{
using QuantLib::FdBlackScholesAsianEngine;
%}

%shared_ptr(FdBlackScholesAsianEngine)
class FdBlackScholesAsianEngine : public PricingEngine {
  public:
    FdBlackScholesAsianEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                              Size tGrid, Size xGrid, Size aGrid);
};

%{
using QuantLib::BaroneAdesiWhaleyApproximationEngine;
%}

%shared_ptr(BaroneAdesiWhaleyApproximationEngine);
%rename(BaroneAdesiWhaleyEngine) BaroneAdesiWhaleyApproximationEngine;
class BaroneAdesiWhaleyApproximationEngine : public PricingEngine {
  public:
    BaroneAdesiWhaleyApproximationEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


%{
using QuantLib::BjerksundStenslandApproximationEngine;
%}

%shared_ptr(BjerksundStenslandApproximationEngine);
%rename(BjerksundStenslandEngine) BjerksundStenslandApproximationEngine;
class BjerksundStenslandApproximationEngine : public PricingEngine {
  public:
    BjerksundStenslandApproximationEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


%{
using QuantLib::JuQuadraticApproximationEngine;
%}

%shared_ptr(JuQuadraticApproximationEngine);
class JuQuadraticApproximationEngine : public PricingEngine {
  public:
    JuQuadraticApproximationEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%{
using QuantLib::AnalyticDigitalAmericanEngine;
%}

%shared_ptr(AnalyticDigitalAmericanEngine)
class AnalyticDigitalAmericanEngine : public PricingEngine {
  public:
    AnalyticDigitalAmericanEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%{
using QuantLib::AnalyticDigitalAmericanKOEngine;
%}

%shared_ptr(AnalyticDigitalAmericanKOEngine)
class AnalyticDigitalAmericanKOEngine : public PricingEngine {
  public:
    AnalyticDigitalAmericanKOEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process);
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
%}

%shared_ptr(AnalyticDividendEuropeanEngine)
class AnalyticDividendEuropeanEngine : public PricingEngine {
  public:
    AnalyticDividendEuropeanEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%{
using QuantLib::FDDividendEuropeanEngine;
using QuantLib::FDDividendAmericanEngine;
%}

%shared_ptr(FDDividendEuropeanEngine<CrankNicolson>)

%rename(FDDividendEuropeanEngineT) FDDividendEuropeanEngine;
template <class S>
class FDDividendEuropeanEngine : public PricingEngine {
  public:
    FDDividendEuropeanEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                             Size timeSteps = 100,
                             Size gridPoints = 100,
                             bool timeDependent = false);
};

%template(FDDividendEuropeanEngine) FDDividendEuropeanEngine<CrankNicolson>;


%shared_ptr(FDDividendAmericanEngine<CrankNicolson>)

%rename(FDDividendAmericanEngineT) FDDividendAmericanEngine;
template <class S>
class FDDividendAmericanEngine : public PricingEngine {
  public:
    FDDividendAmericanEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                             Size timeSteps = 100,
                             Size gridPoints = 100,
                             bool timeDependent = false);
};

%template(FDDividendAmericanEngine) FDDividendAmericanEngine<CrankNicolson>;


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
using QuantLib::MCBarrierEngine;
%}

%shared_ptr(AnalyticBarrierEngine)
class AnalyticBarrierEngine : public PricingEngine {
  public:
    AnalyticBarrierEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>&);
};


%shared_ptr(MCBarrierEngine<PseudoRandom>);
%shared_ptr(MCBarrierEngine<LowDiscrepancy>);

template <class RNG>
class MCBarrierEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCBarrierEngine;
    #endif
  public:
    %extend {
        MCBarrierEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                        intOrNull timeSteps = Null<Size>(),
                        intOrNull timeStepsPerYear = Null<Size>(),
                        bool brownianBridge = false,
                        bool antitheticVariate = false,
                        intOrNull requiredSamples = Null<Size>(),
                        doubleOrNull requiredTolerance = Null<Real>(),
                        intOrNull maxSamples = Null<Size>(),
                        bool isBiased = false,
                        BigInteger seed = 0) {
            return new MCBarrierEngine<RNG>(process,
                                            timeSteps,
                                            timeStepsPerYear,
                                            brownianBridge,
                                            antitheticVariate,
                                            requiredSamples,
                                            requiredTolerance,
                                            maxSamples,
                                            isBiased,
                                            seed);
        }
    }
};

%template(MCPRBarrierEngine) MCBarrierEngine<PseudoRandom>;
%template(MCLDBarrierEngine) MCBarrierEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCBarrierEngine(process,
                        traits,
                        timeSteps=None,
                        timeStepsPerYear=None,
                        brownianBridge=False,
                        antitheticVariate=False,
                        requiredSamples=None,
                        requiredTolerance=None,
                        maxSamples=None,
                        isBiased=False,
                        seed=0):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPRBarrierEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDBarrierEngine
        else:
            raise RuntimeError("unknown MC traits: %s" % traits);
        return cls(process,
                   timeSteps,
                   timeStepsPerYear,
                   brownianBridge,
                   antitheticVariate,
                   requiredSamples,
                   requiredTolerance,
                   maxSamples,
                   isBiased,
                   seed)
%}
#endif


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
using QuantLib::FdBatesVanillaEngine;
using QuantLib::FdHestonVanillaEngine;
%}

%shared_ptr(FdBlackScholesVanillaEngine)
class FdBlackScholesVanillaEngine : public PricingEngine {
  public:
    FdBlackScholesVanillaEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
            Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
            bool localVol = false,
            Real illegalLocalVolOverwrite = -Null<Real>());
};

%shared_ptr(FdBatesVanillaEngine)
class FdBatesVanillaEngine : public PricingEngine {
  public:
    FdBatesVanillaEngine(
            const boost::shared_ptr<BatesModel>& model,
            Size tGrid = 100, Size xGrid = 100,
            Size vGrid=50, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer());
};

%shared_ptr(FdHestonVanillaEngine)
class FdHestonVanillaEngine : public PricingEngine {
  public:
    FdHestonVanillaEngine(
            const boost::shared_ptr<HestonModel>& model,
            Size tGrid = 100, Size xGrid = 100,
            Size vGrid = 50, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
            const boost::shared_ptr<LocalVolTermStructure>& leverageFct
                = boost::shared_ptr<LocalVolTermStructure>());
};


%{
using QuantLib::AnalyticCEVEngine;
using QuantLib::FdCEVVanillaEngine;
%}

%shared_ptr(AnalyticCEVEngine);
class AnalyticCEVEngine : public PricingEngine {
  public:
    AnalyticCEVEngine(
            Real f0, Real alpha, Real beta,
            const Handle<YieldTermStructure>& rTS);
};

%shared_ptr(FdCEVVanillaEngine);
class FdCEVVanillaEngine : public PricingEngine {
  public:
    FdCEVVanillaEngine(
            Real f0, Real alpha, Real beta,
            const Handle<YieldTermStructure>& rTS,
            Size tGrid = 50, Size xGrid = 400,
            Size dampingSteps = 0,
            Real scalingFactor = 1.0, Real eps = 1e-4,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas());
};


%{
using QuantLib::FdSabrVanillaEngine;
%}

%shared_ptr(FdSabrVanillaEngine);
class FdSabrVanillaEngine : public PricingEngine {
  public:
    FdSabrVanillaEngine(
            Real f0, Real alpha, Real beta, Real nu, Real rho,
            const Handle<YieldTermStructure>& rTS,
            Size tGrid = 50, Size fGrid = 400, Size xGrid = 50,
            Size dampingSteps = 0,
            Real scalingFactor = 1.0, Real eps = 1e-4,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer());
};


%{
using QuantLib::FdBlackScholesBarrierEngine;
%}

%shared_ptr(FdBlackScholesBarrierEngine)
class FdBlackScholesBarrierEngine : public PricingEngine {
  public:
    FdBlackScholesBarrierEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
                                const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
                                bool localVol = false,
                                Real illegalLocalVolOverwrite = -Null<Real>());
};


%{
using QuantLib::AnalyticBinaryBarrierEngine;
%}

%shared_ptr(AnalyticBinaryBarrierEngine)
class AnalyticBinaryBarrierEngine : public PricingEngine {
  public:
    AnalyticBinaryBarrierEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


%{
using QuantLib::BinomialBarrierEngine;
using QuantLib::DiscretizedDermanKaniBarrierOption;
%}

#if defined(SWIGPYTHON)
%feature("docstring") BinomialBarrierEngine "Binomial Engine for barrier options.
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

#if !defined(SWIGR)

%shared_ptr(BinomialBarrierEngine<CoxRossRubinstein, DiscretizedDermanKaniBarrierOption>);
%shared_ptr(BinomialBarrierEngine<JarrowRudd, DiscretizedDermanKaniBarrierOption>);
%shared_ptr(BinomialBarrierEngine<AdditiveEQPBinomialTree, DiscretizedDermanKaniBarrierOption>);
%shared_ptr(BinomialBarrierEngine<Trigeorgis, DiscretizedDermanKaniBarrierOption>);
%shared_ptr(BinomialBarrierEngine<Tian, DiscretizedDermanKaniBarrierOption>);
%shared_ptr(BinomialBarrierEngine<LeisenReimer, DiscretizedDermanKaniBarrierOption>);
%shared_ptr(BinomialBarrierEngine<Joshi4, DiscretizedDermanKaniBarrierOption>);

template <class T, class U>
class BinomialBarrierEngine : public PricingEngine {
  public:
    BinomialBarrierEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>&,
                          Size steps,
                          Size max_steps = 0);
};

%template(BinomialCRRBarrierEngine) BinomialBarrierEngine<CoxRossRubinstein, DiscretizedDermanKaniBarrierOption>;
%template(BinomialJRBarrierEngine) BinomialBarrierEngine<JarrowRudd, DiscretizedDermanKaniBarrierOption>;
%template(BinomialEQPBarrierEngine) BinomialBarrierEngine<AdditiveEQPBinomialTree, DiscretizedDermanKaniBarrierOption>;
%template(BinomialTrigeorgisBarrierEngine) BinomialBarrierEngine<Trigeorgis, DiscretizedDermanKaniBarrierOption>;
%template(BinomialTianBarrierEngine) BinomialBarrierEngine<Tian, DiscretizedDermanKaniBarrierOption>;
%template(BinomialLRBarrierEngine) BinomialBarrierEngine<LeisenReimer, DiscretizedDermanKaniBarrierOption>;
%template(BinomialJ4BarrierEngine) BinomialBarrierEngine<Joshi4, DiscretizedDermanKaniBarrierOption>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def BinomialBarrierEngine(process, type, steps):
        type = type.lower()
        if type == "crr" or type == "coxrossrubinstein":
            cls = BinomialCRRBarrierEngine
        elif type == "jr" or type == "jarrowrudd":
            cls = BinomialJRBarrierEngine
        elif type == "eqp":
            cls = BinomialEQPBarrierEngine
        elif type == "trigeorgis":
            cls = BinomialTrigeorgisBarrierEngine
        elif type == "tian":
            cls = BinomialTianBarrierEngine
        elif type == "lr" or type == "leisenreimer":
            cls = BinomialLRBarrierEngine
        elif type == "j4" or type == "joshi4":
            cls = BinomialJ4BarrierEngine
        else:
            raise RuntimeError("unknown binomial engine type: %s" % type);
        return cls(process, steps)
%}
#endif

#endif


%{
using QuantLib::ForwardVanillaEngine;
using QuantLib::QuantoEngine;
typedef ForwardVanillaEngine<AnalyticEuropeanEngine> ForwardEuropeanEngine;
typedef QuantoEngine<VanillaOption,AnalyticEuropeanEngine> QuantoEuropeanEngine;
typedef QuantoEngine<ForwardVanillaOption,AnalyticEuropeanEngine> QuantoForwardEuropeanEngine;
%}


%shared_ptr(ForwardEuropeanEngine)
class ForwardEuropeanEngine: public PricingEngine {
  public:
    ForwardEuropeanEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>&);
};

%shared_ptr(QuantoEuropeanEngine)
class QuantoEuropeanEngine : public PricingEngine {
  public:
    QuantoEuropeanEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                         const Handle<YieldTermStructure>& foreignRiskFreeRate,
                         const Handle<BlackVolTermStructure>& exchangeRateVolatility,
                         const Handle<Quote>& correlation);
};

%shared_ptr(QuantoForwardEuropeanEngine)
class QuantoForwardEuropeanEngine : public PricingEngine {
  public:
    QuantoForwardEuropeanEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                const Handle<YieldTermStructure>& foreignRiskFreeRate,
                                const Handle<BlackVolTermStructure>& exchangeRateVolatility,
                                const Handle<Quote>& correlation);
};


%{
using QuantLib::BlackCalculator;
%}

class BlackCalculator {
  public:
    BlackCalculator(const boost::shared_ptr<StrikedTypePayoff>& payoff,
                    Real forward,
                    Real stdDev,
                    Real discount = 1.0);
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
%}

%shared_ptr(AnalyticContinuousGeometricAveragePriceAsianEngine)
class AnalyticContinuousGeometricAveragePriceAsianEngine : public PricingEngine {
  public:
    AnalyticContinuousGeometricAveragePriceAsianEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


%{
using QuantLib::AnalyticDiscreteGeometricAveragePriceAsianEngine;
%}

%shared_ptr(AnalyticDiscreteGeometricAveragePriceAsianEngine)
class AnalyticDiscreteGeometricAveragePriceAsianEngine : public PricingEngine {
  public:
    AnalyticDiscreteGeometricAveragePriceAsianEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


%{
using QuantLib::AnalyticDiscreteGeometricAverageStrikeAsianEngine;
%}

%shared_ptr(AnalyticDiscreteGeometricAverageStrikeAsianEngine)
class AnalyticDiscreteGeometricAverageStrikeAsianEngine : public PricingEngine {
  public:
    AnalyticDiscreteGeometricAverageStrikeAsianEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


%{
using QuantLib::MCDiscreteArithmeticAPEngine;
using QuantLib::MCDiscreteArithmeticASEngine;
using QuantLib::MCDiscreteGeometricAPEngine;
%}

%shared_ptr(MCDiscreteArithmeticAPEngine<PseudoRandom>);
%shared_ptr(MCDiscreteArithmeticAPEngine<LowDiscrepancy>);

template <class RNG>
class MCDiscreteArithmeticAPEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCDiscreteArithmeticAPEngine;
    #endif
  public:
    %extend {
        MCDiscreteArithmeticAPEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                     bool brownianBridge = false,
                                     bool antitheticVariate = false,
                                     bool controlVariate = false,
                                     intOrNull requiredSamples = Null<Size>(),
                                     doubleOrNull requiredTolerance = Null<Real>(),
                                     intOrNull maxSamples = Null<Size>(),
                                     BigInteger seed = 0) {
            return new MCDiscreteArithmeticAPEngine<RNG>(process,
                                                         brownianBridge,
                                                         antitheticVariate,
                                                         controlVariate,
                                                         requiredSamples,
                                                         requiredTolerance,
                                                         maxSamples,
                                                         seed);
        }
    }
};

%template(MCPRDiscreteArithmeticAPEngine) MCDiscreteArithmeticAPEngine<PseudoRandom>;
%template(MCLDDiscreteArithmeticAPEngine) MCDiscreteArithmeticAPEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCDiscreteArithmeticAPEngine(process,
                                     traits,
                                     brownianBridge=False,
                                     antitheticVariate=False,
                                     controlVariate=False,
                                     requiredSamples=None,
                                     requiredTolerance=None,
                                     maxSamples=None,
                                     seed=0):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPRDiscreteArithmeticAPEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDDiscreteArithmeticAPEngine
        else:
            raise RuntimeError("unknown MC traits: %s" % traits);
        return cls(process,
                   brownianBridge,
                   antitheticVariate,
                   controlVariate,
                   requiredSamples,
                   requiredTolerance,
                   maxSamples,
                   seed)
%}
#endif


%shared_ptr(MCDiscreteArithmeticASEngine<PseudoRandom>);
%shared_ptr(MCDiscreteArithmeticASEngine<LowDiscrepancy>);

template <class RNG>
class MCDiscreteArithmeticASEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCDiscreteArithmeticASEngine;
    #endif
  public:
    %extend {
        MCDiscreteArithmeticASEngine(
                            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                            bool brownianBridge = false,
                            bool antitheticVariate = false,
                            intOrNull requiredSamples = Null<Size>(),
                            doubleOrNull requiredTolerance = Null<Real>(),
                            intOrNull maxSamples = Null<Size>(),
                            BigInteger seed = 0) {
            return new MCDiscreteArithmeticASEngine<RNG>(process,
                                                         brownianBridge,
                                                         antitheticVariate,
                                                         requiredSamples,
                                                         requiredTolerance,
                                                         maxSamples,
                                                         seed);
        }
    }
};

%template(MCPRDiscreteArithmeticASEngine) MCDiscreteArithmeticASEngine<PseudoRandom>;
%template(MCLDDiscreteArithmeticASEngine) MCDiscreteArithmeticASEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCDiscreteArithmeticASEngine(process,
                                     traits,
                                     brownianBridge=False,
                                     antitheticVariate=False,
                                     requiredSamples=None,
                                     requiredTolerance=None,
                                     maxSamples=None,
                                     seed=0):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPRDiscreteArithmeticASEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDDiscreteArithmeticASEngine
        else:
            raise RuntimeError("unknown MC traits: %s" % traits);
        return cls(process,
                   brownianBridge,
                   antitheticVariate,
                   requiredSamples,
                   requiredTolerance,
                   maxSamples,
                   seed)
%}
#endif


%shared_ptr(MCDiscreteGeometricAPEngine<PseudoRandom>);
%shared_ptr(MCDiscreteGeometricAPEngine<LowDiscrepancy>);

template <class RNG>
class MCDiscreteGeometricAPEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCDiscreteGeometricAPEngine;
    #endif
  public:
    %extend {
        MCDiscreteGeometricAPEngine(
                            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                            bool brownianBridge = false,
                            bool antitheticVariate = false,
                            intOrNull requiredSamples = Null<Size>(),
                            doubleOrNull requiredTolerance = Null<Real>(),
                            intOrNull maxSamples = Null<Size>(),
                            BigInteger seed = 0) {
            return new MCDiscreteGeometricAPEngine<RNG>(process,
                                                        brownianBridge,
                                                        antitheticVariate,
                                                        requiredSamples,
                                                        requiredTolerance,
                                                        maxSamples,
                                                        seed);
        }
    }
};


%template(MCPRDiscreteGeometricAPEngine) MCDiscreteGeometricAPEngine<PseudoRandom>;
%template(MCLDDiscreteGeometricAPEngine) MCDiscreteGeometricAPEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCDiscreteGeometricAPEngine(process,
                                     traits,
                                     brownianBridge=False,
                                     antitheticVariate=False,
                                     requiredSamples=None,
                                     requiredTolerance=None,
                                     maxSamples=None,
                                     seed=0):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPRDiscreteGeometricAPEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDDiscreteGeometricAPEngine
        else:
            raise RuntimeError("unknown MC traits: %s" % traits);
        return cls(process,
                   brownianBridge,
                   antitheticVariate,
                   requiredSamples,
                   requiredTolerance,
                   maxSamples,
                   seed)
%}
#endif


%{
using QuantLib::VarianceGammaEngine;
%}

%shared_ptr(VarianceGammaEngine)
class VarianceGammaEngine : public PricingEngine {
  public:
    VarianceGammaEngine(const boost::shared_ptr<VarianceGammaProcess>& process);
};

%{
using QuantLib::FFTVarianceGammaEngine;
%}

%shared_ptr(FFTVarianceGammaEngine)
class FFTVarianceGammaEngine : public PricingEngine {
  public:
    FFTVarianceGammaEngine(const boost::shared_ptr<VarianceGammaProcess>& process,
                           Real logStrikeSpacing = 0.001);
    void precalculate(const std::vector<boost::shared_ptr<Instrument> >& optionList);
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
%}

#if defined(SWIGPYTHON)
%feature("docstring") AnalyticDoubleBarrierEngine "Double barrier engine implementing Ikeda-Kunitomo series."
#endif
%shared_ptr(AnalyticDoubleBarrierEngine)
class AnalyticDoubleBarrierEngine : public PricingEngine {
  public:
    AnalyticDoubleBarrierEngine(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                           int series = 5);
};

%{
using QuantLib::WulinYongDoubleBarrierEngine;
%}

%shared_ptr(WulinYongDoubleBarrierEngine)
class WulinYongDoubleBarrierEngine : public PricingEngine {
  public:
    WulinYongDoubleBarrierEngine(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
                           int series = 5);
};

%{
using QuantLib::VannaVolgaBarrierEngine;
using QuantLib::DeltaVolQuote;
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
%template(RelinkableDeltaVolQuoteHandle) RelinkableHandle<DeltaVolQuote>;

#if defined(SWIGPYTHON)
%feature("docstring") VannaVolgaDoubleBarrierEngine "
Vanna-Volga engine for double barrier options.
Supports different double barrier engines, selected by the type parameters.
Type values:
    ik or analytic:  Ikeda-Kunitomo standard engine (default)
    wo:              Wulin-Yong engine
"
#endif

%{
using QuantLib::VannaVolgaDoubleBarrierEngine;
%}

%shared_ptr(VannaVolgaDoubleBarrierEngine<AnalyticDoubleBarrierEngine>);
%shared_ptr(VannaVolgaDoubleBarrierEngine<WulinYongDoubleBarrierEngine>);

template <class E>
class VannaVolgaDoubleBarrierEngine : public PricingEngine {
  public:
    VannaVolgaDoubleBarrierEngine(
                           const Handle<DeltaVolQuote> atmVol,
                           const Handle<DeltaVolQuote> vol25Put,
                           const Handle<DeltaVolQuote> vol25Call,
                           const Handle<Quote> spotFX,
                           const Handle<YieldTermStructure> domesticTS,
                           const Handle<YieldTermStructure> foreignTS,
                           const bool adaptVanDelta = false,
                           const Real bsPriceWithSmile = 0.0,
                           int series = 5);
};

%template(VannaVolgaIKDoubleBarrierEngine) VannaVolgaDoubleBarrierEngine<AnalyticDoubleBarrierEngine>;
%template(VannaVolgaWODoubleBarrierEngine) VannaVolgaDoubleBarrierEngine<WulinYongDoubleBarrierEngine>;


%shared_ptr(VannaVolgaBarrierEngine)
class VannaVolgaBarrierEngine : public PricingEngine {
  public:
    VannaVolgaBarrierEngine(
                const Handle<DeltaVolQuote>& atmVol,
                const Handle<DeltaVolQuote>& vol25Put,
                const Handle<DeltaVolQuote>& vol25Call,
                const Handle<Quote>& spotFX,
                const Handle<YieldTermStructure>& domesticTS,
                const Handle<YieldTermStructure>& foreignTS,
                const bool adaptVanDelta = false,
                const Real bsPriceWithSmile = 0.0);
};

%{
using QuantLib::AnalyticDoubleBarrierBinaryEngine;
%}

%shared_ptr(AnalyticDoubleBarrierBinaryEngine)
class AnalyticDoubleBarrierBinaryEngine : public PricingEngine {
  public:
    AnalyticDoubleBarrierBinaryEngine(
                           const boost::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%{
using QuantLib::BinomialDoubleBarrierEngine;
using QuantLib::DiscretizedDermanKaniDoubleBarrierOption;
%}

#if defined(SWIGPYTHON)
%feature("docstring") BinomialDoubleBarrierEngine "Binomial Engine for double barrier options.
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

#if !defined(SWIGR)

%shared_ptr(BinomialDoubleBarrierEngine<CoxRossRubinstein, DiscretizedDermanKaniDoubleBarrierOption>);
%shared_ptr(BinomialDoubleBarrierEngine<JarrowRudd, DiscretizedDermanKaniDoubleBarrierOption>);
%shared_ptr(BinomialDoubleBarrierEngine<AdditiveEQPBinomialTree, DiscretizedDermanKaniDoubleBarrierOption>);
%shared_ptr(BinomialDoubleBarrierEngine<Trigeorgis, DiscretizedDermanKaniDoubleBarrierOption>);
%shared_ptr(BinomialDoubleBarrierEngine<Tian, DiscretizedDermanKaniDoubleBarrierOption>);
%shared_ptr(BinomialDoubleBarrierEngine<LeisenReimer, DiscretizedDermanKaniDoubleBarrierOption>);
%shared_ptr(BinomialDoubleBarrierEngine<Joshi4, DiscretizedDermanKaniDoubleBarrierOption>);

template <class T, class U>
class BinomialDoubleBarrierEngine : public PricingEngine {
  public:
    BinomialDoubleBarrierEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>&,
                                Size steps);
};

%template(BinomialCRRDoubleBarrierEngine) BinomialDoubleBarrierEngine<CoxRossRubinstein, DiscretizedDermanKaniDoubleBarrierOption>;
%template(BinomialJRDoubleBarrierEngine) BinomialDoubleBarrierEngine<JarrowRudd, DiscretizedDermanKaniDoubleBarrierOption>;
%template(BinomialEQPDoubleBarrierEngine) BinomialDoubleBarrierEngine<AdditiveEQPBinomialTree, DiscretizedDermanKaniDoubleBarrierOption>;
%template(BinomialTrigeorgisDoubleBarrierEngine) BinomialDoubleBarrierEngine<Trigeorgis, DiscretizedDermanKaniDoubleBarrierOption>;
%template(BinomialTianDoubleBarrierEngine) BinomialDoubleBarrierEngine<Tian, DiscretizedDermanKaniDoubleBarrierOption>;
%template(BinomialLRDoubleBarrierEngine) BinomialDoubleBarrierEngine<LeisenReimer, DiscretizedDermanKaniDoubleBarrierOption>;
%template(BinomialJ4DoubleBarrierEngine) BinomialDoubleBarrierEngine<Joshi4, DiscretizedDermanKaniDoubleBarrierOption>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def BinomialDoubleBarrierEngine(process, type, steps):
        type = type.lower()
        if type == "crr" or type == "coxrossrubinstein":
            cls = BinomialCRRDoubleBarrierEngine
        elif type == "jr" or type == "jarrowrudd":
            cls = BinomialJRDoubleBarrierEngine
        elif type == "eqp":
            cls = BinomialEQPDoubleBarrierEngine
        elif type == "trigeorgis":
            cls = BinomialTrigeorgisDoubleBarrierEngine
        elif type == "tian":
            cls = BinomialTianDoubleBarrierEngine
        elif type == "lr" or type == "leisenreimer":
            cls = BinomialLRDoubleBarrierEngine
        elif type == "j4" or type == "joshi4":
            cls = BinomialJ4DoubleBarrierEngine
        else:
            raise RuntimeError("unknown binomial engine type: %s" % type);
        return cls(process, steps)
%}
#endif

#endif


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
%}

%shared_ptr(FdSimpleBSSwingEngine)
class FdSimpleBSSwingEngine : public PricingEngine {
  public:
    FdSimpleBSSwingEngine(
            const boost::shared_ptr<GeneralizedBlackScholesProcess>& process,
            Size tGrid = 50, Size xGrid = 100,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas());
};

%shared_ptr(FdSimpleExtOUJumpSwingEngine)
class FdSimpleExtOUJumpSwingEngine : public PricingEngine {
  public:
    %extend {
        FdSimpleExtOUJumpSwingEngine(
            const boost::shared_ptr<ExtOUWithJumpsProcess>& process,
            const boost::shared_ptr<YieldTermStructure>& rTS,
            Size tGrid = 50, Size xGrid = 200, Size yGrid=50,
            const std::vector<std::pair<Time,Real> >& shape =
                                         std::vector<std::pair<Time,Real> >(),
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer()) {

            boost::shared_ptr<FdSimpleExtOUJumpSwingEngine::Shape> curve(
                              new FdSimpleExtOUJumpSwingEngine::Shape(shape));

            return new FdSimpleExtOUJumpSwingEngine(
                    process, rTS, tGrid, xGrid, yGrid,
                    curve, schemeDesc);
        }
    }
};


#endif
