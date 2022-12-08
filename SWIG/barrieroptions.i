/*
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
 Copyright (C) 2015 Thema Consulting SA
 Copyright (C) 2018, 2019 Matthias Lungwitz
 Copyright (C) 2022 Ignacio Anguita

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

#ifndef quantlib_barrier_options_i
#define quantlib_barrier_options_i

%include options.i

%{
using QuantLib::Barrier;
%}

struct Barrier {
    enum Type { DownIn, UpIn, DownOut, UpOut };
};

%{
using QuantLib::BarrierOption;
using QuantLib::DividendBarrierOption;
using QuantLib::QuantoBarrierOption;
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

%shared_ptr(QuantoBarrierOption)
class QuantoBarrierOption : public BarrierOption {
  public:
    QuantoBarrierOption(
               Barrier::Type barrierType,
               Real barrier,
               Real rebate,
               const ext::shared_ptr<StrikedTypePayoff>& payoff,
               const ext::shared_ptr<Exercise>& exercise);
};

%{
using QuantLib::PartialBarrier;
%}

struct PartialBarrier : public Barrier {
    enum Range { Start, End, EndB1, EndB2 };
};

%{
using QuantLib::PartialTimeBarrierOption;
%}

%shared_ptr(PartialTimeBarrierOption)
class PartialTimeBarrierOption : public OneAssetOption {
      public:
        PartialTimeBarrierOption(PartialBarrier::Type barrierType,
            PartialBarrier::Range barrierRange,
            Real barrier,
            Real rebate,
            Date coverEventDate,
            const ext::shared_ptr<StrikedTypePayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise);
}; 

%{
using QuantLib::AnalyticPartialTimeBarrierOptionEngine;
%}

#if defined(SWIGPYTHON)
%feature("docstring") AnalyticPartialTimeBarrierOptionEngine "Partial Time Barrier Option Engine"
#endif
%shared_ptr(AnalyticPartialTimeBarrierOptionEngine )
class AnalyticPartialTimeBarrierOptionEngine : public PricingEngine {
  public:
    AnalyticPartialTimeBarrierOptionEngine (
                           const ext::shared_ptr<GeneralizedBlackScholesProcess>& process
                 );
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
using QuantLib::QuantoEngine;
typedef QuantoEngine<BarrierOption,AnalyticBarrierEngine> QuantoBarrierEngine;
%}

%shared_ptr(QuantoBarrierEngine);
class QuantoBarrierEngine : public PricingEngine {
      public:
        QuantoBarrierEngine(ext::shared_ptr<GeneralizedBlackScholesProcess>,
                            Handle<YieldTermStructure> foreignRiskFreeRate,
                            Handle<BlackVolTermStructure> exchangeRateVolatility,
                            Handle<Quote> correlation);
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
using QuantLib::VannaVolgaBarrierEngine;
%}

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
using QuantLib::DoubleBarrierOption;
using QuantLib::DoubleBarrier;
%}

struct DoubleBarrier {
    enum Type { KnockIn, KnockOut, KIKO, KOKI };
};

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
using QuantLib::SuoWangDoubleBarrierEngine;
%}

%shared_ptr(SuoWangDoubleBarrierEngine)
class SuoWangDoubleBarrierEngine : public PricingEngine {
  public:
    SuoWangDoubleBarrierEngine(
                           const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                           int series = 5);
};

deprecate_feature(WulinYongDoubleBarrierEngine, SuoWangDoubleBarrierEngine);


%{
using QuantLib::VannaVolgaDoubleBarrierEngine;
%}

%shared_ptr(VannaVolgaDoubleBarrierEngine<AnalyticDoubleBarrierEngine>);
%shared_ptr(VannaVolgaDoubleBarrierEngine<SuoWangDoubleBarrierEngine>);

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
%template(VannaVolgaWODoubleBarrierEngine) VannaVolgaDoubleBarrierEngine<SuoWangDoubleBarrierEngine>;

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

#endif
