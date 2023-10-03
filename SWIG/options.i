/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2008 Tito Ingargiola
 Copyright (C) 2010, 2012, 2018, 2019 Klaus Spanderen
 Copyright (C) 2015 Thema Consulting SA
 Copyright (C) 2016 Gouthaman Balaraman
 Copyright (C) 2018, 2019 Matthias Lungwitz
 Copyright (C) 2019 Pedro Coelho

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
%include dividends.i
%include exercise.i
%include stochasticprocess.i
%include instruments.i
%include stl.i
%include linearalgebra.i
%include calibratedmodel.i
%include grid.i
%include parameter.i
%include vectors.i

// payoff

%{
using QuantLib::Payoff;
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

// option types
%{
using QuantLib::Option;
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


%{
using QuantLib::TypePayoff;
using QuantLib::FloatingTypePayoff;
using QuantLib::StrikedTypePayoff;
%}

%shared_ptr(TypePayoff)
class TypePayoff : public Payoff {
  public:
    Option::Type optionType();
  private:
    TypePayoff();
};

%shared_ptr(FloatingTypePayoff)
class FloatingTypePayoff : public TypePayoff {
  public:
    FloatingTypePayoff(Option::Type type);
    Real operator()(Real price, Real strike) const;
    Real operator()(Real price) const;
};

%shared_ptr(StrikedTypePayoff)
class StrikedTypePayoff : public TypePayoff
{
  public:
    Real strike();
  private:
    StrikedTypePayoff();
};


%{
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
%}
#endif

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
    Volatility impliedVolatility(
                         Real targetValue,
                         const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                         const DividendSchedule& dividends,
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
    AnalyticEuropeanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
    AnalyticEuropeanEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                           const Handle<YieldTermStructure>& discountCurve);
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
using QuantLib::MCDigitalEngine;
using QuantLib::PseudoRandom;
using QuantLib::LowDiscrepancy;
using QuantLib::LsmBasisSystem;
%}

struct LsmBasisSystem {
    enum PolynomialType { Monomial, Laguerre, Hermite, Hyperbolic,
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
                         LsmBasisSystem::PolynomialType polynomType = LsmBasisSystem::Monomial,
                         int nCalibrationSamples = 2048,
                         ext::optional<bool> antitheticVariateCalibration = ext::nullopt,
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



%shared_ptr(MCDigitalEngine<PseudoRandom>);
%shared_ptr(MCDigitalEngine<LowDiscrepancy>);

template <class RNG>
class MCDigitalEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCDigitalEngine;
    #endif
  public:
    %extend {
        MCDigitalEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
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
            return new MCDigitalEngine<RNG>(process,
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

%template(MCPRDigitalEngine) MCDigitalEngine<PseudoRandom>;
%template(MCLDDigitalEngine) MCDigitalEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCDigitalEngine(process,
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
            cls = MCPRDigitalEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDDigitalEngine
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


// American engines

%{
using QuantLib::BaroneAdesiWhaleyApproximationEngine;
%}

%shared_ptr(BaroneAdesiWhaleyApproximationEngine);
class BaroneAdesiWhaleyApproximationEngine : public PricingEngine {
  public:
    BaroneAdesiWhaleyApproximationEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


%{
using QuantLib::BjerksundStenslandApproximationEngine;
%}

%shared_ptr(BjerksundStenslandApproximationEngine);
class BjerksundStenslandApproximationEngine : public PricingEngine {
  public:
    BjerksundStenslandApproximationEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


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


#if defined(SWIGR)
%Rruntime %{
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
    AnalyticDividendEuropeanEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
            DividendSchedule dividends);
};

%{
using QuantLib::QdPlusAmericanEngine;
%}

%shared_ptr(QdPlusAmericanEngine)
class QdPlusAmericanEngine: public PricingEngine {
  public:
    enum SolverType {Brent, Newton, Ridder, Halley, SuperHalley};

    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") QdPlusAmericanEngine;
    #endif
    explicit QdPlusAmericanEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        Size interpolationPoints = 8,
        SolverType solverType = Halley,
        Real eps = 1e-6,
        Size maxIter = Null<Size>());        
};

%{
using QuantLib::QdFpLegendreScheme;
using QuantLib::QdFpIterationScheme;
using QuantLib::QdFpLegendreTanhSinhScheme;
using QuantLib::QdFpTanhSinhIterationScheme;
using QuantLib::QdFpAmericanEngine;
%}

%shared_ptr(QdFpIterationScheme)
class QdFpIterationScheme {
  private:
    QdFpIterationSchem();
};

%shared_ptr(QdFpLegendreScheme)
class QdFpLegendreScheme: public QdFpIterationScheme {
  public:
    QdFpLegendreScheme(Size l, Size m, Size n, Size p);
};

%shared_ptr(QdFpLegendreTanhSinhScheme)
class QdFpLegendreTanhSinhScheme: public QdFpLegendreScheme {
  public:
    QdFpLegendreTanhSinhScheme(Size l, Size m, Size n, Real eps);
};

%shared_ptr(QdFpTanhSinhIterationScheme)
class QdFpTanhSinhIterationScheme : public QdFpIterationScheme {
  public:
    QdFpTanhSinhIterationScheme(Size m, Size n, Real eps);
};


%shared_ptr(QdFpAmericanEngine)
class QdFpAmericanEngine : public PricingEngine {
  public:
    enum FixedPointEquation { FP_A, FP_B, Auto };

    explicit QdFpAmericanEngine(
      ext::shared_ptr<GeneralizedBlackScholesProcess> bsProcess,
      ext::shared_ptr<QdFpIterationScheme> iterationScheme =
          accurateScheme(),
      FixedPointEquation fpEquation = Auto);
      
    static ext::shared_ptr<QdFpIterationScheme> fastScheme();
    static ext::shared_ptr<QdFpIterationScheme> accurateScheme();
    static ext::shared_ptr<QdFpIterationScheme> highPrecisionScheme();      
};


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
using QuantLib::FdBlackScholesShoutEngine;
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

    FdBlackScholesVanillaEngine(
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        DividendSchedule dividends,
        Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
        bool localVol = false,
        Real illegalLocalVolOverwrite = -Null<Real>(),
        CashDividendModel cashDividendModel = Spot);

    FdBlackScholesVanillaEngine(
        const ext::shared_ptr<GeneralizedBlackScholesProcess>&,
        DividendSchedule dividends,
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
                    const DividendSchedule& dividends = {},
                    const ext::shared_ptr<FdmQuantoHelper>& quantoHelper = {},
                    Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
                    const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas(),
                    bool localVol = false,
                    Real illegalLocalVolOverwrite = -Null<Real>(),
                    CashDividendModel cashDividendModel = Spot) {
            if (dividends.empty()) {
                return ext::make_shared<FdBlackScholesVanillaEngine>(
                                                process, quantoHelper, tGrid, xGrid,
                                                dampingSteps, schemeDesc,
                                                localVol, illegalLocalVolOverwrite,
                                                cashDividendModel);
            } else {
                return ext::make_shared<FdBlackScholesVanillaEngine>(
                                                process, dividends, quantoHelper, tGrid, xGrid,
                                                dampingSteps, schemeDesc,
                                                localVol, illegalLocalVolOverwrite,
                                                cashDividendModel);
            }
        }
    }
    #endif
};

%shared_ptr(FdBlackScholesShoutEngine)
class FdBlackScholesShoutEngine : public PricingEngine {
  public:
    FdBlackScholesShoutEngine(
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas());
    FdBlackScholesShoutEngine(
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        DividendSchedule dividends,
        Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas());
};

%shared_ptr(FdOrnsteinUhlenbeckVanillaEngine)
class FdOrnsteinUhlenbeckVanillaEngine : public PricingEngine {
  public:
    FdOrnsteinUhlenbeckVanillaEngine(
        const ext::shared_ptr<OrnsteinUhlenbeckProcess>&,
        const ext::shared_ptr<YieldTermStructure>& rTS,
        Size tGrid = 100, Size xGrid = 100, Size dampingSteps = 0,
        Real epsilon = 0.0001,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Douglas());
    FdOrnsteinUhlenbeckVanillaEngine(
        const ext::shared_ptr<OrnsteinUhlenbeckProcess>&,
        const ext::shared_ptr<YieldTermStructure>& rTS,
        DividendSchedule dividends,
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
    FdBatesVanillaEngine(
            const ext::shared_ptr<BatesModel>& model,
            DividendSchedule dividends,
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
        const ext::shared_ptr<LocalVolTermStructure>& leverageFct = {},
        const Real mixingFactor = 1.0);

    FdHestonVanillaEngine(
        const ext::shared_ptr<HestonModel>& model,
        const ext::shared_ptr<FdmQuantoHelper>& quantoHelper,
        Size tGrid = 100, 
        Size xGrid = 100,
        Size vGrid = 50, 
        Size dampingSteps = 0,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
        const ext::shared_ptr<LocalVolTermStructure>& leverageFct = {},
        const Real mixingFactor = 1.0);

    FdHestonVanillaEngine(
        const ext::shared_ptr<HestonModel>& model,
        DividendSchedule dividends,
        Size tGrid = 100, Size xGrid = 100,
        Size vGrid = 50, Size dampingSteps = 0,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
        const ext::shared_ptr<LocalVolTermStructure>& leverageFct = {},
        const Real mixingFactor = 1.0);

    FdHestonVanillaEngine(
        const ext::shared_ptr<HestonModel>& model,
        DividendSchedule dividends,
        const ext::shared_ptr<FdmQuantoHelper>& quantoHelper,
        Size tGrid = 100, 
        Size xGrid = 100,
        Size vGrid = 50, 
        Size dampingSteps = 0,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
        const ext::shared_ptr<LocalVolTermStructure>& leverageFct = {},
        const Real mixingFactor = 1.0);

    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") make;
    %extend {
        static ext::shared_ptr<FdHestonVanillaEngine> make(
                    const ext::shared_ptr<HestonModel>& model,
                    const DividendSchedule& dividends = {},
                    const ext::shared_ptr<FdmQuantoHelper>& quantoHelper = {},
                    Size tGrid = 100, Size xGrid = 100, Size vGrid = 50,
                    Size dampingSteps = 0,
                    const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
                    const ext::shared_ptr<LocalVolTermStructure>& leverageFct = {},
                    const Real mixingFactor = 1.0) {
            if (dividends.empty()) {
                return ext::make_shared<FdHestonVanillaEngine>(
                    model, quantoHelper, tGrid, xGrid, vGrid,
                    dampingSteps, schemeDesc, leverageFct, mixingFactor);
            } else {
                return ext::make_shared<FdHestonVanillaEngine>(
                    model, dividends, quantoHelper, tGrid, xGrid, vGrid,
                    dampingSteps, schemeDesc, leverageFct, mixingFactor);
            }
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
using QuantLib::FdHestonHullWhiteVanillaEngine;
%}

%shared_ptr(FdHestonHullWhiteVanillaEngine);
class FdHestonHullWhiteVanillaEngine : public PricingEngine {
  public:
    FdHestonHullWhiteVanillaEngine(
        const ext::shared_ptr<HestonModel>& model,
        ext::shared_ptr<HullWhiteProcess> hwProcess,
        Real corrEquityShortRate,
        Size tGrid = 50,
        Size xGrid = 100,
        Size vGrid = 40,
        Size rGrid = 20,
        Size dampingSteps = 0,
        bool controlVariate = true,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer());    
    FdHestonHullWhiteVanillaEngine(
        const ext::shared_ptr<HestonModel>& model,
        ext::shared_ptr<HullWhiteProcess> hwProcess,
        DividendSchedule dividends,
        Real corrEquityShortRate,
        Size tGrid = 50,
        Size xGrid = 100,
        Size vGrid = 40,
        Size rGrid = 20,
        Size dampingSteps = 0,
        bool controlVariate = true,
        const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer());    
};


%{
using QuantLib::AnalyticHestonHullWhiteEngine;
%}

%shared_ptr(AnalyticHestonHullWhiteEngine);
class AnalyticHestonHullWhiteEngine : public PricingEngine {
  public:
    AnalyticHestonHullWhiteEngine(const ext::shared_ptr<HestonModel>& hestonModel,
                                  ext::shared_ptr<HullWhite> hullWhiteModel,
                                  Size integrationOrder = 144);

    AnalyticHestonHullWhiteEngine(const ext::shared_ptr<HestonModel>& model,
                                  ext::shared_ptr<HullWhite> hullWhiteModel,
                                  Real relTolerance,
                                  Size maxEvaluations);
};


%{
using QuantLib::AnalyticH1HWEngine;
%}

%shared_ptr(AnalyticH1HWEngine);
class AnalyticH1HWEngine : public PricingEngine {
  public:
    AnalyticH1HWEngine(const ext::shared_ptr<HestonModel>& hestonModel,
                       const ext::shared_ptr<HullWhite>& hullWhiteModel,
                       Real rhoSr, Size integrationOrder = 144);

    AnalyticH1HWEngine(const ext::shared_ptr<HestonModel>& model,
                       const ext::shared_ptr<HullWhite>& hullWhiteModel,
                       Real rhoSr,
                       Real relTolerance,
                       Size maxEvaluations);
};


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
    Real strikeGamma() const;
    Real alpha() const;
    Real beta() const;
};




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
using QuantLib::MargrabeOption;
using QuantLib::AnalyticEuropeanMargrabeEngine;
using QuantLib::AnalyticAmericanMargrabeEngine;
%}

%shared_ptr(MargrabeOption)
class MargrabeOption : public MultiAssetOption {
  public:
    MargrabeOption(Integer Q1,
                   Integer Q2,
                   const ext::shared_ptr<Exercise>&);
    Real delta1() const;
    Real delta2() const;
    Real gamma1() const;
    Real gamma2() const;
};

%shared_ptr(AnalyticEuropeanMargrabeEngine)
class AnalyticEuropeanMargrabeEngine : public PricingEngine {
  public:
    AnalyticEuropeanMargrabeEngine(ext::shared_ptr<GeneralizedBlackScholesProcess> process1,
                                   ext::shared_ptr<GeneralizedBlackScholesProcess> process2,
                                   Real correlation);
};

%shared_ptr(AnalyticAmericanMargrabeEngine)
class AnalyticAmericanMargrabeEngine : public PricingEngine {
  public:
    AnalyticAmericanMargrabeEngine(ext::shared_ptr<GeneralizedBlackScholesProcess> process1,
                                   ext::shared_ptr<GeneralizedBlackScholesProcess> process2,
                                   Real correlation);
};


%{
using QuantLib::CompoundOption;
using QuantLib::AnalyticCompoundOptionEngine;
%}

%shared_ptr(CompoundOption)
class CompoundOption : public OneAssetOption {
  public:
    CompoundOption(const ext::shared_ptr<StrikedTypePayoff>& motherPayoff,
                   const ext::shared_ptr<Exercise>& motherExercise,
                   ext::shared_ptr<StrikedTypePayoff> daughterPayoff,
                   ext::shared_ptr<Exercise> daughterExercise);
};

%shared_ptr(AnalyticCompoundOptionEngine)
class AnalyticCompoundOptionEngine : public PricingEngine {
  public:
    AnalyticCompoundOptionEngine(ext::shared_ptr<GeneralizedBlackScholesProcess> process);
};


%{
using QuantLib::SimpleChooserOption;
using QuantLib::AnalyticSimpleChooserEngine;
%}

%shared_ptr(SimpleChooserOption)
class SimpleChooserOption : public OneAssetOption {
  public:
    SimpleChooserOption(Date choosingDate,
                        Real strike,
                        const ext::shared_ptr<Exercise>& exercise);
};

%shared_ptr(AnalyticSimpleChooserEngine)
class AnalyticSimpleChooserEngine : public PricingEngine {
  public:
    AnalyticSimpleChooserEngine(ext::shared_ptr<GeneralizedBlackScholesProcess> process);
};


%{
using QuantLib::ComplexChooserOption;
using QuantLib::AnalyticComplexChooserEngine;
%}

%shared_ptr(ComplexChooserOption)
class ComplexChooserOption : public OneAssetOption {
  public:
    ComplexChooserOption(Date choosingDate,
                         Real strikeCall,
                         Real strikePut,
                         const ext::shared_ptr<Exercise>& exerciseCall,
                         const ext::shared_ptr<Exercise>& exercisePut);
};

%shared_ptr(AnalyticComplexChooserEngine)
class AnalyticComplexChooserEngine : public PricingEngine {
  public:
    AnalyticComplexChooserEngine(ext::shared_ptr<GeneralizedBlackScholesProcess> process);
};


#endif
