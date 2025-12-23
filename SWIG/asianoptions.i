/*
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
 Copyright (C) 2010, 2012, 2018, 2019 Klaus Spanderen
 Copyright (C) 2018, 2019 Matthias Lungwitz
 Copyright (C) 2020, 2022 Jack Gillett

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

#ifndef quantlib_asian_options_i
#define quantlib_asian_options_i

%include options.i
%include grid.i

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
    ContinuousAveragingAsianOption(
            Average::Type averageType,
            Date startDate,
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
            std::vector<Date> fixingDates,
            const ext::shared_ptr<StrikedTypePayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise);
    DiscreteAveragingAsianOption(
            Average::Type averageType,
            std::vector<Date> fixingDates,
            const ext::shared_ptr<StrikedTypePayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise,
            std::vector<Real> allPastFixings = std::vector<Real>());
    %extend {
        TimeGrid timeGrid() {
            return self->result<TimeGrid>("TimeGrid");
        }
    }
};


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
using QuantLib::ContinuousArithmeticAsianLevyEngine;
%}

%shared_ptr(ContinuousArithmeticAsianLevyEngine)
class ContinuousArithmeticAsianLevyEngine : public PricingEngine {
  public:
    ContinuousArithmeticAsianLevyEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                                        const Handle<Quote>& runningAverage);
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
using QuantLib::ChoiAsianEngine;
%}

%shared_ptr(ChoiAsianEngine)
class ChoiAsianEngine : public PricingEngine {
  public:
    ChoiAsianEngine(
        ext::shared_ptr<GeneralizedBlackScholesProcess> process,
        Real lambda = 15,
        Size maxNrIntegrationSteps = 2 << 21);
};

%{
using QuantLib::TurnbullWakemanAsianEngine;
%}

%shared_ptr(TurnbullWakemanAsianEngine)
class TurnbullWakemanAsianEngine : public PricingEngine {
  public:
    TurnbullWakemanAsianEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


#endif
