/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2008 Tito Ingargiola
 Copyright (C) 2010, 2012, 2018, 2019 Klaus Spanderen
 Copyright (C) 2015 Thema Consulting SA
 Copyright (C) 2016 Gouthaman Balaraman
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
    ext::shared_ptr<Payoff> payoff();
    ext::shared_ptr<Exercise> exercise();
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
            const ext::shared_ptr<StrikedTypePayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise);

    Volatility impliedVolatility(
                         Real targetValue,
                         const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
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

%template(CalibrationPair) std::pair< ext::shared_ptr<VanillaOption>, ext::shared_ptr<Quote> >;
%template(CalibrationSet) std::vector<std::pair< ext::shared_ptr<VanillaOption>, ext::shared_ptr<Quote> > >;

%{
using QuantLib::EuropeanOption;
%}


%shared_ptr(EuropeanOption)
class EuropeanOption : public VanillaOption {
  public:
    EuropeanOption(
            const ext::shared_ptr<StrikedTypePayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise);
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
                const ext::shared_ptr<StrikedTypePayoff>& payoff,
                const ext::shared_ptr<Exercise>& exercise);
};

// QuantoVanillaOption

%{
using QuantLib::QuantoVanillaOption;
%}

%shared_ptr(QuantoVanillaOption)
class QuantoVanillaOption : public OneAssetOption {
  public:
    QuantoVanillaOption(
            const ext::shared_ptr<StrikedTypePayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise);
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
            const ext::shared_ptr<StrikedTypePayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise);
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
    AnalyticEuropeanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>&);
};


%{
using QuantLib::HestonModel;
%}

%shared_ptr(HestonModel)
class HestonModel : public CalibratedModel {
  public:
    HestonModel(const ext::shared_ptr<HestonProcess>&  process);
    Real theta() const;
    Real kappa() const;
    Real sigma() const;
    Real rho() const;
    Real v0() const;
};

%template(HestonModelHandle) Handle<HestonModel>;

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
%rename (AnalyticHestonEngine_Integration) AnalyticHestonEngine::Integration;
%shared_ptr(AnalyticHestonEngine)
#if !defined(SWIGCSHARP)
%feature ("flatnested") AnalyticHestonEngine::Integration;
#endif

class AnalyticHestonEngine : public PricingEngine {
  public:
    class Integration
    {
    public:
        // non adaptive integration algorithms based on Gaussian quadrature
        static Integration gaussLaguerre    (Size integrationOrder = 128);
        static Integration gaussLegendre    (Size integrationOrder = 128);
        static Integration gaussChebyshev   (Size integrationOrder = 128);
        static Integration gaussChebyshev2nd(Size integrationOrder = 128);

        // for an adaptive integration algorithm Gatheral's version has to
        // be used.Be aware: using a too large number for maxEvaluations might
        // result in a stack overflow as the these integrations are based on
        // recursive algorithms.
        static Integration gaussLobatto(Real relTolerance, Real absTolerance,
                                        Size maxEvaluations = 1000);

        // usually these routines have a poor convergence behavior.
        static Integration gaussKronrod(Real absTolerance,
                                        Size maxEvaluations = 1000);
        static Integration simpson(Real absTolerance,
                                   Size maxEvaluations = 1000);
        static Integration trapezoid(Real absTolerance,
                                     Size maxEvaluations = 1000);
        static Integration discreteSimpson(Size evaluation = 1000);
        static Integration discreteTrapezoid(Size evaluation = 1000);

        static Real andersenPiterbargIntegrationLimit(
            Real c_inf, Real epsilon, Real v0, Real t);

        Real calculate(Real c_inf,
                       const ext::function<Real(Real)>& f,
                       doubleOrNull maxBound = Null<Real>()) const;

        Size numberOfEvaluations() const;
        bool isAdaptiveIntegration() const;

    private:
      enum Algorithm
        { GaussLobatto, GaussKronrod, Simpson, Trapezoid,
          DiscreteTrapezoid, DiscreteSimpson,
          GaussLaguerre, GaussLegendre,
          GaussChebyshev, GaussChebyshev2nd };

      Integration(Algorithm intAlgo,
                const ext::shared_ptr<GaussianQuadrature>& quadrature);

      Integration(Algorithm intAlgo,
                const ext::shared_ptr<Integrator>& integrator);
    };
    enum ComplexLogFormula { 
        Gatheral, BranchCorrection, AndersenPiterbarg, 
        AndersenPiterbargOptCV, AsymptoticChF, OptimalCV
    };
    AnalyticHestonEngine(const ext::shared_ptr<HestonModel>& model,
                         Size integrationOrder = 144);
    AnalyticHestonEngine(const ext::shared_ptr<HestonModel>& model,
                         Real relTolerance,
                         Size maxEvaluations);
    AnalyticHestonEngine(const ext::shared_ptr<HestonModel>& model,
                     ComplexLogFormula cpxLog, const AnalyticHestonEngine::Integration& itg,
                     Real andersenPiterbargEpsilon = 1e-8);

    %extend {                     
        std::pair<Real, Real> chF(Real real, Real imag, Time t) const {
            const std::complex<Real> tmp 
                = self->chF(std::complex<Real>(real, imag), t);
            return std::pair<Real, Real>(tmp.real(), tmp.imag());
        }
    }
};

%{
using QuantLib::COSHestonEngine;
%}

%shared_ptr(COSHestonEngine)
class COSHestonEngine : public PricingEngine {
  public:
    COSHestonEngine(const ext::shared_ptr<HestonModel>& model,
                    Real L = 16, Size N = 200);
};

%{
using QuantLib::ExponentialFittingHestonEngine;
%}

%shared_ptr(ExponentialFittingHestonEngine)
class ExponentialFittingHestonEngine : public PricingEngine {
  public:
    enum ControlVariate { AndersenPiterbarg, AndersenPiterbargOptCV,
                          AsymptoticChF, OptimalCV };
    
    ExponentialFittingHestonEngine(
        const ext::shared_ptr<HestonModel>& model,
        ControlVariate cv = AndersenPiterbargOptCV,
        doubleOrNull scaling = Null<Real>());
};


%{
using QuantLib::AnalyticPTDHestonEngine;
%}

%shared_ptr(AnalyticPTDHestonEngine)
class AnalyticPTDHestonEngine : public PricingEngine {
  public:
    enum ComplexLogFormula { Gatheral, AndersenPiterbarg };
    typedef AnalyticHestonEngine::Integration Integration;
    
    AnalyticPTDHestonEngine(
            const ext::shared_ptr<PiecewiseTimeDependentHestonModel>& model,
            Real relTolerance, Size maxEvaluations);
    // Constructor using Laguerre integration
    // and Gatheral's version of complex log.
    AnalyticPTDHestonEngine(
            const ext::shared_ptr<PiecewiseTimeDependentHestonModel>& model,
            Size integrationOrder = 144);
            
    // Constructor giving full control over Fourier integration algorithm
    AnalyticPTDHestonEngine(
        const ext::shared_ptr<PiecewiseTimeDependentHestonModel>& model,
        ComplexLogFormula cpxLog,
        const Integration& itg,
        Real andersenPiterbargEpsilon = 1e-8);            
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
    BatesModel(const ext::shared_ptr<BatesProcess>&  process);
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
    BatesEngine(const ext::shared_ptr<BatesModel>& model,
                Size integrationOrder = 144);
    BatesEngine(const ext::shared_ptr<BatesModel>& model,
                Real relTolerance,
                Size maxEvaluations);
};


%{
using QuantLib::IntegralEngine;
%}

%shared_ptr(IntegralEngine)
class IntegralEngine : public PricingEngine {
  public:
    IntegralEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>&);
};


%{
using QuantLib::CrankNicolson;
using QuantLib::FDBermudanEngine;
using QuantLib::FDEuropeanEngine;
%}

%shared_ptr(FDBermudanEngine<CrankNicolson>);

template <class S>
class FDBermudanEngine : public PricingEngine {
    #if defined(SWIGPYTHON)
    %pythonprepend FDBermudanEngine %{
        from warnings import warn
        warn("FDBermudanEngine is deprecated; use FdBlackScholesVanillaEngine")
    %}
    #endif
  public:
    FDBermudanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                     Size timeSteps = 100, Size gridPoints = 100,
                     bool timeDependent = false);
};

%template(FDBermudanEngine) FDBermudanEngine<CrankNicolson>;

%shared_ptr(FDEuropeanEngine<CrankNicolson>);

template <class S>
class FDEuropeanEngine : public PricingEngine {
    #if defined(SWIGPYTHON)
    %pythonprepend FDEuropeanEngine %{
        from warnings import warn
        warn("FDEuropeanEngine is deprecated; use FdBlackScholesVanillaEngine")
    %}
    #endif
  public:
    FDEuropeanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess> process,
                     Size timeSteps = 100, Size gridPoints = 100,
                     bool timeDependent = false);
};

%template(FDEuropeanEngine) FDEuropeanEngine<CrankNicolson>;


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
    BinomialVanillaEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>&,
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
        MCEuropeanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
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
        MCAmericanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
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
                   seedCalibration if seedCalibration is not None else nullInt())
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
        MCEuropeanHestonEngine(const ext::shared_ptr<HestonProcess>& process,
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
using QuantLib::FDAmericanEngine;
using QuantLib::FDShoutEngine;
%}

%shared_ptr(FDAmericanEngine<CrankNicolson>);

template <class S>
class FDAmericanEngine : public PricingEngine {
    #if defined(SWIGPYTHON)
    %pythonprepend FDAmericanEngine %{
        from warnings import warn
        warn("FDAmericanEngine is deprecated; use FdBlackScholesVanillaEngine")
    %}
    #endif
  public:
    FDAmericanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                     Size timeSteps = 100, Size gridPoints = 100,
                     bool timeDependent = false);
};

%template(FDAmericanEngine) FDAmericanEngine<CrankNicolson>;


%shared_ptr(FDShoutEngine<CrankNicolson>);

template <class S>
class FDShoutEngine : public PricingEngine {
  public:
    FDShoutEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                  Size timeSteps = 100, Size gridPoints = 100,
                  bool timeDependent = false);
};

%template(FDShoutEngine) FDShoutEngine<CrankNicolson>;


%{
using QuantLib::ContinuousArithmeticAsianLevyEngine;
%}

%shared_ptr(ContinuousArithmeticAsianLevyEngine)
class ContinuousArithmeticAsianLevyEngine : public PricingEngine {
  public:
    ContinuousArithmeticAsianLevyEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                        const Handle<Quote>& runningAverage,
                                        const Date& startDate);
};

%{
using QuantLib::FdBlackScholesAsianEngine;
%}

%shared_ptr(FdBlackScholesAsianEngine)
class FdBlackScholesAsianEngine : public PricingEngine {
  public:
    FdBlackScholesAsianEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                              Size tGrid, Size xGrid, Size aGrid);
};

%{
using QuantLib::BaroneAdesiWhaleyApproximationEngine;
%}

%shared_ptr(BaroneAdesiWhaleyApproximationEngine);
class BaroneAdesiWhaleyApproximationEngine : public PricingEngine {
  public:
    BaroneAdesiWhaleyApproximationEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};
deprecate_feature(BaroneAdesiWhaleyEngine, BaroneAdesiWhaleyApproximationEngine);


%{
using QuantLib::BjerksundStenslandApproximationEngine;
%}

%shared_ptr(BjerksundStenslandApproximationEngine);
class BjerksundStenslandApproximationEngine : public PricingEngine {
  public:
    BjerksundStenslandApproximationEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};
deprecate_feature(BjerksundStenslandEngine, BjerksundStenslandApproximationEngine);


%{
using QuantLib::JuQuadraticApproximationEngine;
%}

%shared_ptr(JuQuadraticApproximationEngine);
class JuQuadraticApproximationEngine : public PricingEngine {
  public:
    JuQuadraticApproximationEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%{
using QuantLib::AnalyticDigitalAmericanEngine;
%}

%shared_ptr(AnalyticDigitalAmericanEngine)
class AnalyticDigitalAmericanEngine : public PricingEngine {
  public:
    AnalyticDigitalAmericanEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%{
using QuantLib::AnalyticDigitalAmericanKOEngine;
%}

%shared_ptr(AnalyticDigitalAmericanKOEngine)
class AnalyticDigitalAmericanKOEngine : public PricingEngine {
  public:
    AnalyticDigitalAmericanKOEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

// Dividend option

%{
using QuantLib::DividendVanillaOption;
%}


%shared_ptr(DividendVanillaOption)
class DividendVanillaOption : public OneAssetOption {
  public:
    DividendVanillaOption(
            const ext::shared_ptr<StrikedTypePayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise,
            const std::vector<Date>& dividendDates,
            const std::vector<Real>& dividends);
    Volatility impliedVolatility(
                         Real targetValue,
                         const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                         Real accuracy = 1.0e-4,
                         Size maxEvaluations = 100,
                         Volatility minVol = 1.0e-4,
                         Volatility maxVol = 4.0);
};


%{
using QuantLib::AnalyticDividendEuropeanEngine;
%}

%shared_ptr(AnalyticDividendEuropeanEngine)
class AnalyticDividendEuropeanEngine : public PricingEngine {
  public:
    AnalyticDividendEuropeanEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%{
using QuantLib::FDDividendEuropeanEngine;
using QuantLib::FDDividendAmericanEngine;
%}

%shared_ptr(FDDividendEuropeanEngine<CrankNicolson>)

%rename(FDDividendEuropeanEngineT) FDDividendEuropeanEngine;
template <class S>
class FDDividendEuropeanEngine : public PricingEngine {
    #if defined(SWIGPYTHON)
    %pythonprepend FDDividendEuropeanEngine %{
        from warnings import warn
        warn("FDDividendEuropeanEngine is deprecated; use FdBlackScholesVanillaEngine")
    %}
    #endif
  public:
    FDDividendEuropeanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                             Size timeSteps = 100,
                             Size gridPoints = 100,
                             bool timeDependent = false);
};

%template(FDDividendEuropeanEngine) FDDividendEuropeanEngine<CrankNicolson>;


%shared_ptr(FDDividendAmericanEngine<CrankNicolson>)

%rename(FDDividendAmericanEngineT) FDDividendAmericanEngine;
template <class S>
class FDDividendAmericanEngine : public PricingEngine {
    #if defined(SWIGPYTHON)
    %pythonprepend FDDividendAmericanEngine %{
        from warnings import warn
        warn("FDDividendAmericanEngine is deprecated; use FdBlackScholesVanillaEngine")
    %}
    #endif
  public:
    FDDividendAmericanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                             Size timeSteps = 100,
                             Size gridPoints = 100,
                             bool timeDependent = false);
};

%template(FDDividendAmericanEngine) FDDividendAmericanEngine<CrankNicolson>;


// Barrier option

%{
using QuantLib::BarrierOption;
using QuantLib::DividendBarrierOption;
%}

%shared_ptr(BarrierOption)
class BarrierOption : public OneAssetOption {
  public:
    BarrierOption(
               Barrier::Type barrierType,
               Real barrier,
               Real rebate,
               const ext::shared_ptr<StrikedTypePayoff>& payoff,
               const ext::shared_ptr<Exercise>& exercise);
    Volatility impliedVolatility(
                         Real targetValue,
                         const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                         Real accuracy = 1.0e-4,
                         Size maxEvaluations = 100,
                         Volatility minVol = 1.0e-4,
                         Volatility maxVol = 4.0);
};

%shared_ptr(DividendBarrierOption)
class DividendBarrierOption : public BarrierOption {
  public:
    DividendBarrierOption(Barrier::Type barrierType,
                          Real barrier,
                          Real rebate,
                          const ext::shared_ptr<StrikedTypePayoff>& payoff,
                          const ext::shared_ptr<Exercise>& exercise,
                          const std::vector<Date>& dividendDates,
                          const std::vector<Real>& dividends);
};


// Barrier engines

%{
using QuantLib::AnalyticBarrierEngine;
using QuantLib::MCBarrierEngine;
%}

%shared_ptr(AnalyticBarrierEngine)
class AnalyticBarrierEngine : public PricingEngine {
  public:
    AnalyticBarrierEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>&);
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
        MCBarrierEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
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
                       ImplicitEulerType, ExplicitEulerType,
                       MethodOfLinesType, TrBDF2Type, 
                       CrankNicolsonType };

  FdmSchemeDesc(FdmSchemeType type, Real theta, Real mu);

  const FdmSchemeType type;
  const Real theta, mu;

  // some default scheme descriptions
  static FdmSchemeDesc Douglas();
  static FdmSchemeDesc CrankNicolson();
  static FdmSchemeDesc ImplicitEuler();
  static FdmSchemeDesc ExplicitEuler();
  static FdmSchemeDesc CraigSneyd();
  static FdmSchemeDesc ModifiedCraigSneyd();
  static FdmSchemeDesc Hundsdorfer();
  static FdmSchemeDesc ModifiedHundsdorfer();
  static FdmSchemeDesc MethodOfLines(
      Real eps=0.001, Real relInitStepSize=0.01);
  static FdmSchemeDesc TrBDF2();
};

%{
using QuantLib::FdmQuantoHelper;
%}

%shared_ptr(FdmQuantoHelper)
class FdmQuantoHelper {
  public:
    FdmQuantoHelper(
        const ext::shared_ptr<YieldTermStructure>& rTS,
        const ext::shared_ptr<YieldTermStructure>& fTS,
        const ext::shared_ptr<BlackVolTermStructure>& fxVolTS,
        Real equityFxCorrelation,
        Real exchRateATMlevel);
};

%{
using QuantLib::LocalVolTermStructure;
using QuantLib::FdBlackScholesVanillaEngine;
using QuantLib::FdOrnsteinUhlenbeckVanillaEngine;
using QuantLib::FdBatesVanillaEngine;
using QuantLib::FdHestonVanillaEngine;
%}

%shared_ptr(FdBlackScholesVanillaEngine)
class FdBlackScholesVanillaEngine : public PricingEngine {
  public:
    enum CashDividendModel { Spot, Escrowed };

    FdBlackScholesVanillaEngine(
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
        bool localVol = false,
        Real illegalLocalVolOverwrite = -Null<Real>(),
        CashDividendModel cashDividendModel = Spot);

    FdBlackScholesVanillaEngine(
        const ext::shared_ptr<GeneralizedBlackScholesProcess>&,
        const ext::shared_ptr<FdmQuantoHelper>& quantoHelper,
        Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
        bool localVol = false,
        Real illegalLocalVolOverwrite = -Null<Real>(),
        CashDividendModel cashDividendModel = Spot);

    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") make;
    %extend {
        static ext::shared_ptr<FdBlackScholesVanillaEngine> make(
                    const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                    const ext::shared_ptr<FdmQuantoHelper>& quantoHelper
                        = ext::shared_ptr<FdmQuantoHelper>(),
                    Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
                    const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
                    bool localVol = false,
                    Real illegalLocalVolOverwrite = -Null<Real>(),
                    CashDividendModel cashDividendModel = Spot) {
            return ext::shared_ptr<FdBlackScholesVanillaEngine>(
                new FdBlackScholesVanillaEngine(process, quantoHelper, tGrid, xGrid,
                                                dampingSteps, schemeDesc,
                                                localVol, illegalLocalVolOverwrite,
                                                cashDividendModel));
        }
    }
    #endif
};

%shared_ptr(FdOrnsteinUhlenbeckVanillaEngine)
class FdOrnsteinUhlenbeckVanillaEngine : public PricingEngine {
  public:
    #if !defined(SWIGPYTHON)
    %feature("kwargs") FdOrnsteinUhlenbeckVanillaEngine;
    #endif
    FdOrnsteinUhlenbeckVanillaEngine(
        const ext::shared_ptr<OrnsteinUhlenbeckProcess>&,
        const ext::shared_ptr<YieldTermStructure>& rTS,
        Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
        Real epsilon = 0.0001,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas());
};

%shared_ptr(FdBatesVanillaEngine)
class FdBatesVanillaEngine : public PricingEngine {
  public:
    FdBatesVanillaEngine(
            const ext::shared_ptr<BatesModel>& model,
            Size tGrid = 100, Size xGrid = 100,
            Size vGrid=50, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer());
};

%shared_ptr(FdHestonVanillaEngine)
class FdHestonVanillaEngine : public PricingEngine {
  public:
    FdHestonVanillaEngine(
        const ext::shared_ptr<HestonModel>& model,
        Size tGrid = 100, Size xGrid = 100,
        Size vGrid = 50, Size dampingSteps = 0,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
        const ext::shared_ptr<LocalVolTermStructure>& leverageFct
            = ext::shared_ptr<LocalVolTermStructure>(),
        const Real mixingFactor = 1.0);

    FdHestonVanillaEngine(
        const ext::shared_ptr<HestonModel>& model,
        const ext::shared_ptr<FdmQuantoHelper>& quantoHelper,
        Size tGrid = 100, 
        Size xGrid = 100,
        Size vGrid = 50, 
        Size dampingSteps = 0,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
        const ext::shared_ptr<LocalVolTermStructure>& leverageFct
            = ext::shared_ptr<LocalVolTermStructure>(),
        const Real mixingFactor = 1.0);

    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") make;
    %extend {
        static ext::shared_ptr<FdHestonVanillaEngine> make(
                    const ext::shared_ptr<HestonModel>& model,
                    const ext::shared_ptr<FdmQuantoHelper>& quantoHelper
                        = ext::shared_ptr<FdmQuantoHelper>(),
                    Size tGrid = 100, Size xGrid = 100, Size vGrid = 50,
                    Size dampingSteps = 0,
                    const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
                    const ext::shared_ptr<LocalVolTermStructure>& leverageFct
                        = ext::shared_ptr<LocalVolTermStructure>(),
                    const Real mixingFactor = 1.0) {
            return ext::shared_ptr<FdHestonVanillaEngine>(
                new FdHestonVanillaEngine(model, quantoHelper, tGrid, xGrid, vGrid,
                                          dampingSteps, schemeDesc, leverageFct, mixingFactor));
        }
    }
    #endif
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
using QuantLib::FdBlackScholesRebateEngine;
using QuantLib::FdHestonBarrierEngine;
using QuantLib::FdHestonRebateEngine;
%}

%shared_ptr(FdBlackScholesBarrierEngine)
class FdBlackScholesBarrierEngine : public PricingEngine {
  public:
    FdBlackScholesBarrierEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
                                const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
                                bool localVol = false,
                                Real illegalLocalVolOverwrite = -Null<Real>());
};

%shared_ptr(FdBlackScholesRebateEngine)
class FdBlackScholesRebateEngine : public PricingEngine {
  public:
    FdBlackScholesRebateEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                               Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
                               const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
                               bool localVol = false,
                               Real illegalLocalVolOverwrite = -Null<Real>());
};

%shared_ptr(FdHestonBarrierEngine)
class FdHestonBarrierEngine : public PricingEngine {
  public:
    FdHestonBarrierEngine(const ext::shared_ptr<HestonModel>& model,
                          Size tGrid = 100, Size xGrid = 100, Size vGrid = 50, Size dampingSteps = 0,
                          const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
                          const ext::shared_ptr<LocalVolTermStructure>& leverageFct
                              = ext::shared_ptr<LocalVolTermStructure>(),
                          const Real mixingFactor = 1.0);
};

%shared_ptr(FdHestonRebateEngine)
class FdHestonRebateEngine : public PricingEngine {
  public:
    FdHestonRebateEngine(const ext::shared_ptr<HestonModel>& model,
                         Size tGrid = 100, Size xGrid = 100, Size vGrid = 50, Size dampingSteps = 0,
                         const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
                         const ext::shared_ptr<LocalVolTermStructure>& leverageFct
                             = ext::shared_ptr<LocalVolTermStructure>(),
                         const Real mixingFactor = 1.0);
};


%{
using QuantLib::AnalyticBinaryBarrierEngine;
%}

%shared_ptr(AnalyticBinaryBarrierEngine)
class AnalyticBinaryBarrierEngine : public PricingEngine {
  public:
    AnalyticBinaryBarrierEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
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
    BinomialBarrierEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>&,
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
class ForwardEuropeanEngine : public PricingEngine {
  public:
    ForwardEuropeanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>&);
};

%shared_ptr(QuantoEuropeanEngine)
class QuantoEuropeanEngine : public PricingEngine {
  public:
    QuantoEuropeanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                         const Handle<YieldTermStructure>& foreignRiskFreeRate,
                         const Handle<BlackVolTermStructure>& exchangeRateVolatility,
                         const Handle<Quote>& correlation);
};

%shared_ptr(QuantoForwardEuropeanEngine)
class QuantoForwardEuropeanEngine : public PricingEngine {
  public:
    QuantoForwardEuropeanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                const Handle<YieldTermStructure>& foreignRiskFreeRate,
                                const Handle<BlackVolTermStructure>& exchangeRateVolatility,
                                const Handle<Quote>& correlation);
};


%{
using QuantLib::AnalyticHestonForwardEuropeanEngine;
using QuantLib::MCForwardEuropeanBSEngine;
using QuantLib::MCForwardEuropeanHestonEngine;
%}

%shared_ptr(AnalyticHestonForwardEuropeanEngine)
class AnalyticHestonForwardEuropeanEngine : public PricingEngine {
  public:
    AnalyticHestonForwardEuropeanEngine(const ext::shared_ptr<HestonProcess>& process,
                                        Size integrationOrder = 144);
};

%shared_ptr(MCForwardEuropeanBSEngine<PseudoRandom>);
%shared_ptr(MCForwardEuropeanBSEngine<LowDiscrepancy>);

template <class RNG>
class MCForwardEuropeanBSEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCForwardEuropeanBSEngine;
    #endif
  public:
    %extend {
        MCForwardEuropeanBSEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                  intOrNull timeSteps = Null<Size>(),
                                  intOrNull timeStepsPerYear = Null<Size>(),
                                  bool brownianBridge = false,
                                  bool antitheticVariate = false,
                                  intOrNull requiredSamples = Null<Size>(),
                                  doubleOrNull requiredTolerance = Null<Real>(),
                                  intOrNull maxSamples = Null<Size>(),
                                  BigInteger seed = 0) {
            return new MCForwardEuropeanBSEngine<RNG>(process,
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

%template(MCPRForwardEuropeanBSEngine) MCForwardEuropeanBSEngine<PseudoRandom>;
%template(MCLDForwardEuropeanBSEngine) MCForwardEuropeanBSEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCForwardEuropeanBSEngine(process,
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
            cls = MCPRForwardEuropeanBSEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDForwardEuropeanBSEngine
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


%shared_ptr(MCForwardEuropeanHestonEngine<PseudoRandom>);
%shared_ptr(MCForwardEuropeanHestonEngine<LowDiscrepancy>);

template <class RNG>
class MCForwardEuropeanHestonEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCForwardEuropeanHestonEngine;
    #endif
  public:
    %extend {
        MCForwardEuropeanHestonEngine(const ext::shared_ptr<HestonProcess>& process,
                                      intOrNull timeSteps = Null<Size>(),
                                      intOrNull timeStepsPerYear = Null<Size>(),
                                      bool antitheticVariate = false,
                                      intOrNull requiredSamples = Null<Size>(),
                                      doubleOrNull requiredTolerance = Null<Real>(),
                                      intOrNull maxSamples = Null<Size>(),
                                      BigInteger seed = 0,
                                      bool controlVariate = false) {
            return new MCForwardEuropeanHestonEngine<RNG>(process,
                                                          timeSteps,
                                                          timeStepsPerYear,
                                                          antitheticVariate,
                                                          requiredSamples,
                                                          requiredTolerance,
                                                          maxSamples,
                                                          seed,
                                                          controlVariate);
        }
    }
};

%template(MCPRForwardEuropeanHestonEngine) MCForwardEuropeanHestonEngine<PseudoRandom>;
%template(MCLDForwardEuropeanHestonEngine) MCForwardEuropeanHestonEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCForwardEuropeanHestonEngine(process,
                                      traits,
                                      timeSteps=None,
                                      timeStepsPerYear=None,
                                      antitheticVariate=False,
                                      requiredSamples=None,
                                      requiredTolerance=None,
                                      maxSamples=None,
                                      seed=0,
                                      controlVariate=False):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPRForwardEuropeanHestonEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDForwardEuropeanHestonEngine
        else:
            raise RuntimeError("unknown MC traits: %s" % traits);
        return cls(process,
                   timeSteps,
                   timeStepsPerYear,
                   antitheticVariate,
                   requiredSamples,
                   requiredTolerance,
                   maxSamples,
                   seed,
                   controlVariate)
%}
#endif




%{
using QuantLib::BlackCalculator;
%}

class BlackCalculator {
  public:
    BlackCalculator(const ext::shared_ptr<StrikedTypePayoff>& payoff,
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
            const ext::shared_ptr<StrikedTypePayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise);
};

%shared_ptr(DiscreteAveragingAsianOption)
class DiscreteAveragingAsianOption : public OneAssetOption {
  public:
    DiscreteAveragingAsianOption(
            Average::Type averageType,
            Real runningAccumulator,
            Size pastFixings,
            const std::vector<Date>& fixingDates,
            const ext::shared_ptr<StrikedTypePayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise);
    %extend {
        TimeGrid timeGrid() {
            return self->result<TimeGrid>("TimeGrid");
        }
    }
};

// Asian engines


%{
using QuantLib::AnalyticContinuousGeometricAveragePriceAsianEngine;
using QuantLib::AnalyticContinuousGeometricAveragePriceAsianHestonEngine;
using QuantLib::AnalyticDiscreteGeometricAveragePriceAsianEngine;
using QuantLib::AnalyticDiscreteGeometricAveragePriceAsianHestonEngine;
using QuantLib::AnalyticDiscreteGeometricAverageStrikeAsianEngine;
%}

%shared_ptr(AnalyticContinuousGeometricAveragePriceAsianEngine)
class AnalyticContinuousGeometricAveragePriceAsianEngine : public PricingEngine {
  public:
    AnalyticContinuousGeometricAveragePriceAsianEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%shared_ptr(AnalyticContinuousGeometricAveragePriceAsianHestonEngine)
class AnalyticContinuousGeometricAveragePriceAsianHestonEngine : public PricingEngine {
  public:
    AnalyticContinuousGeometricAveragePriceAsianHestonEngine(
            const ext::shared_ptr<HestonProcess>& process,
            Size summationCutoff = 50,
            Real xiRightLimit = 100.0);
};

%shared_ptr(AnalyticDiscreteGeometricAveragePriceAsianEngine)
class AnalyticDiscreteGeometricAveragePriceAsianEngine : public PricingEngine {
  public:
    AnalyticDiscreteGeometricAveragePriceAsianEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%shared_ptr(AnalyticDiscreteGeometricAveragePriceAsianHestonEngine)
class AnalyticDiscreteGeometricAveragePriceAsianHestonEngine : public PricingEngine {
  public:
    AnalyticDiscreteGeometricAveragePriceAsianHestonEngine(
            const ext::shared_ptr<HestonProcess>& process,
            Real xiRightLimit = 100.0);
};

%shared_ptr(AnalyticDiscreteGeometricAverageStrikeAsianEngine)
class AnalyticDiscreteGeometricAverageStrikeAsianEngine : public PricingEngine {
  public:
    AnalyticDiscreteGeometricAverageStrikeAsianEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


%{
using QuantLib::MCDiscreteArithmeticAPEngine;
using QuantLib::MCDiscreteArithmeticAPHestonEngine;
using QuantLib::MCDiscreteArithmeticASEngine;
using QuantLib::MCDiscreteGeometricAPEngine;
using QuantLib::MCDiscreteGeometricAPHestonEngine;
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
        MCDiscreteArithmeticAPEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
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

%shared_ptr(MCDiscreteArithmeticAPHestonEngine<PseudoRandom>);
%shared_ptr(MCDiscreteArithmeticAPHestonEngine<LowDiscrepancy>);

template <class RNG>
class MCDiscreteArithmeticAPHestonEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCDiscreteArithmeticAPHestonEngine;
    #endif
  public:
    %extend {
        MCDiscreteArithmeticAPHestonEngine(const ext::shared_ptr<HestonProcess>& process,
                                           bool antitheticVariate = false,
                                           intOrNull requiredSamples = Null<Size>(),
                                           doubleOrNull requiredTolerance = Null<Real>(),
                                           intOrNull maxSamples = Null<Size>(),
                                           BigInteger seed = 0,
                                           intOrNull timeSteps = Null<Size>(),
                                           intOrNull timeStepsPerYear = Null<Size>(),
                                           bool controlVariate = false) {
            return new MCDiscreteArithmeticAPHestonEngine<RNG>(process,
                                                               antitheticVariate,
                                                               requiredSamples,
                                                               requiredTolerance,
                                                               maxSamples,
                                                               seed,
                                                               timeSteps,
                                                               timeStepsPerYear,
                                                               controlVariate);
        }
    }
};

%template(MCPRDiscreteArithmeticAPHestonEngine) MCDiscreteArithmeticAPHestonEngine<PseudoRandom>;
%template(MCLDDiscreteArithmeticAPHestonEngine) MCDiscreteArithmeticAPHestonEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCDiscreteArithmeticAPHestonEngine(process,
                                           traits,
                                           antitheticVariate=False,
                                           requiredSamples=None,
                                           requiredTolerance=None,
                                           maxSamples=None,
                                           seed=0,
                                           timeSteps=None,
                                           timeStepsPerYear=None,
                                           controlVariate=False):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPRDiscreteArithmeticAPHestonEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDDiscreteArithmeticAPHestonEngine
        else:
            raise RuntimeError("unknown MC traits: %s" % traits);
        return cls(process,
                   antitheticVariate,
                   requiredSamples,
                   requiredTolerance,
                   maxSamples,
                   seed,
                   timeSteps,
                   timeStepsPerYear,
                   controlVariate)
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
                            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
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
                            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
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

%shared_ptr(MCDiscreteGeometricAPHestonEngine<PseudoRandom>);
%shared_ptr(MCDiscreteGeometricAPHestonEngine<LowDiscrepancy>);

template <class RNG>
class MCDiscreteGeometricAPHestonEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCDiscreteGeometricAPHestonEngine;
    #endif
  public:
    %extend {
        MCDiscreteGeometricAPHestonEngine(const ext::shared_ptr<HestonProcess>& process,
                                          bool antitheticVariate = false,
                                          intOrNull requiredSamples = Null<Size>(),
                                          doubleOrNull requiredTolerance = Null<Real>(),
                                          intOrNull maxSamples = Null<Size>(),
                                          BigInteger seed = 0,
                                          intOrNull timeSteps = Null<Size>(),
                                          intOrNull timeStepsPerYear = Null<Size>()) {
            return new MCDiscreteGeometricAPHestonEngine<RNG>(process,
                                                              antitheticVariate,
                                                              requiredSamples,
                                                              requiredTolerance,
                                                              maxSamples,
                                                              seed,
                                                              timeSteps,
                                                              timeStepsPerYear);
        }
    }
};

%template(MCPRDiscreteGeometricAPHestonEngine) MCDiscreteGeometricAPHestonEngine<PseudoRandom>;
%template(MCLDDiscreteGeometricAPHestonEngine) MCDiscreteGeometricAPHestonEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCDiscreteGeometricAPHestonEngine(process,
                                          traits,
                                          antitheticVariate=False,
                                          requiredSamples=None,
                                          requiredTolerance=None,
                                          maxSamples=None,
                                          seed=0,
                                          timeSteps=None,
                                          timeStepsPerYear=None):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPRDiscreteGeometricAPHestonEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDDiscreteGeometricAPHestonEngine
        else:
            raise RuntimeError("unknown MC traits: %s" % traits);
        return cls(process,
                   antitheticVariate,
                   requiredSamples,
                   requiredTolerance,
                   maxSamples,
                   seed,
                   timeSteps,
                   timeStepsPerYear)
%}
#endif


%{
using QuantLib::VarianceGammaEngine;
%}

%shared_ptr(VarianceGammaEngine)
class VarianceGammaEngine : public PricingEngine {
  public:
    VarianceGammaEngine(const ext::shared_ptr<VarianceGammaProcess>& process);
};

%{
using QuantLib::FFTVarianceGammaEngine;
%}

%shared_ptr(FFTVarianceGammaEngine)
class FFTVarianceGammaEngine : public PricingEngine {
  public:
    FFTVarianceGammaEngine(const ext::shared_ptr<VarianceGammaProcess>& process,
                           Real logStrikeSpacing = 0.001);
    void precalculate(const std::vector<ext::shared_ptr<Instrument> >& optionList);
};

// Double barrier options
%{
using QuantLib::DoubleBarrierOption;
using QuantLib::DoubleBarrier;
%}

%shared_ptr(DoubleBarrierOption)
class DoubleBarrierOption : public OneAssetOption {
  public:
    DoubleBarrierOption(
               DoubleBarrier::Type barrierType,
               Real barrier_lo,
               Real barrier_hi,
               Real rebate,
               const ext::shared_ptr<StrikedTypePayoff>& payoff,
               const ext::shared_ptr<Exercise>& exercise);
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
            const ext::shared_ptr<StrikedTypePayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise);
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
                           const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                           int series = 5);
};

%{
using QuantLib::FdHestonDoubleBarrierEngine;
%}

%shared_ptr(FdHestonDoubleBarrierEngine);
class FdHestonDoubleBarrierEngine : public PricingEngine {
  public:
    FdHestonDoubleBarrierEngine(
            const ext::shared_ptr<HestonModel>& model,
            Size tGrid = 100, Size xGrid = 100,
            Size vGrid = 50, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
            const ext::shared_ptr<LocalVolTermStructure>& leverageFct
                = ext::shared_ptr<LocalVolTermStructure>(),
            const Real mixingFactor = 1.0);
};

%{
using QuantLib::WulinYongDoubleBarrierEngine;
%}

%shared_ptr(WulinYongDoubleBarrierEngine)
class WulinYongDoubleBarrierEngine : public PricingEngine {
  public:
    WulinYongDoubleBarrierEngine(
                           const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
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

    Real delta() const;
    Time maturity() const;
    AtmType atmType() const;
    DeltaType deltaType() const;
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
                           const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
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
    BinomialDoubleBarrierEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>&,
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
        const ext::shared_ptr<Payoff>& payoff,
        const ext::shared_ptr<SwingExercise>& ex,
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


%{
using QuantLib::GJRGARCHModel;
%}

%shared_ptr(GJRGARCHModel)
class GJRGARCHModel : public CalibratedModel {
      public:
        GJRGARCHModel(const ext::shared_ptr<GJRGARCHProcess>& process);
        Real omega() const;
        Real alpha() const;
        Real beta() const;
        Real gamma() const;
        Real lambda() const;
        Real v0() const;
};


%{
using QuantLib::AnalyticGJRGARCHEngine;
%}

%shared_ptr(AnalyticGJRGARCHEngine)
class AnalyticGJRGARCHEngine : public PricingEngine {
  public:
    AnalyticGJRGARCHEngine(const ext::shared_ptr<GJRGARCHModel>& process);
};

%{
using QuantLib::MCEuropeanGJRGARCHEngine;
%}

%shared_ptr(MCEuropeanGJRGARCHEngine<PseudoRandom>);
%shared_ptr(MCEuropeanGJRGARCHEngine<LowDiscrepancy>);

template <class RNG>
class MCEuropeanGJRGARCHEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCEuropeanGJRGARCHEngine;
    #endif
  public:
    %extend {
        MCEuropeanGJRGARCHEngine(const ext::shared_ptr<GJRGARCHProcess>& process,
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
            return new MCEuropeanGJRGARCHEngine<RNG>(process,
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

%template(MCPREuropeanGJRGARCHEngine) MCEuropeanGJRGARCHEngine<PseudoRandom>;
%template(MCLDEuropeanGJRGARCHEngine) MCEuropeanGJRGARCHEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCEuropeanGJRGARCHEngine(process,
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
            cls = MCPREuropeanGJRGARCHEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDEuropeanGJRGARCHEngine
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


%{
using QuantLib::SpreadOption;
using QuantLib::KirkSpreadOptionEngine;
%}

%shared_ptr(SpreadOption);
class SpreadOption : public MultiAssetOption {
public:
  SpreadOption(const ext::shared_ptr<PlainVanillaPayoff>& payoff,
               const ext::shared_ptr<Exercise>& exercise);
};

%shared_ptr(KirkSpreadOptionEngine);
class KirkSpreadOptionEngine : public PricingEngine {
public:
  KirkSpreadOptionEngine(const ext::shared_ptr<BlackProcess>& process1,
                         const ext::shared_ptr<BlackProcess>& process2,
                         const Handle<Quote>& correlation);
};

#endif
