
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2007 Joseph Wang
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

#ifndef quantlib_basket_options_i
#define quantlib_basket_options_i

%include date.i
%include options.i
%include payoffs.i

%{
using QuantLib::BasketOption;
using QuantLib::BasketPayoff;
using QuantLib::MinBasketPayoff;
using QuantLib::MaxBasketPayoff;
using QuantLib::AverageBasketPayoff;
using QuantLib::SpreadBasketPayoff;
%}

%shared_ptr(BasketPayoff)
class BasketPayoff : public Payoff {
  private:
    BasketPayoff();
};

%shared_ptr(MinBasketPayoff)
class MinBasketPayoff : public BasketPayoff  {
  public:
    MinBasketPayoff(const ext::shared_ptr<Payoff> p);
};

%shared_ptr(MaxBasketPayoff)
class MaxBasketPayoff : public BasketPayoff  {
  public:
    MaxBasketPayoff(const ext::shared_ptr<Payoff> p);
};

%shared_ptr(AverageBasketPayoff)
class AverageBasketPayoff :
      public BasketPayoff  {
  public:
    AverageBasketPayoff(const ext::shared_ptr<Payoff> p,
                           const Array &a);
    AverageBasketPayoff(const ext::shared_ptr<Payoff> p,
                           Size n);
};

%shared_ptr(SpreadBasketPayoff)
class SpreadBasketPayoff : public BasketPayoff  {
  public:
    SpreadBasketPayoff(const ext::shared_ptr<Payoff> p);
};

%shared_ptr(BasketOption)
class BasketOption : public MultiAssetOption {
  public:
    BasketOption(
            const ext::shared_ptr<BasketPayoff>& payoff,
            const ext::shared_ptr<Exercise>& exercise);
};


%{
using QuantLib::MCEuropeanBasketEngine;
using QuantLib::MCAmericanBasketEngine;
%}

%shared_ptr(MCEuropeanBasketEngine<PseudoRandom>);
%shared_ptr(MCEuropeanBasketEngine<LowDiscrepancy>);

template <class RNG>
class MCEuropeanBasketEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCEuropeanBasketEngine;
    #endif
  public:
    %extend {
        MCEuropeanBasketEngine(const ext::shared_ptr<StochasticProcessArray>& process,
                               intOrNull timeSteps = Null<Size>(),
                               intOrNull timeStepsPerYear = Null<Size>(),
                               bool brownianBridge = false,
                               bool antitheticVariate = false,
                               intOrNull requiredSamples = Null<Size>(),
                               doubleOrNull requiredTolerance = Null<Real>(),
                               intOrNull maxSamples = Null<Size>(),
                               BigInteger seed = 0) {
            return new MCEuropeanBasketEngine<RNG>(process,
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

%template(MCPREuropeanBasketEngine) MCEuropeanBasketEngine<PseudoRandom>;
%template(MCLDEuropeanBasketEngine) MCEuropeanBasketEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCEuropeanBasketEngine(process,
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
            cls = MCPREuropeanBasketEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDEuropeanBasketEngine
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


%shared_ptr(MCAmericanBasketEngine<PseudoRandom>);
%shared_ptr(MCAmericanBasketEngine<LowDiscrepancy>);

template <class RNG>
class MCAmericanBasketEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCAmericanBasketEngine;
    #endif
  public:
    %extend {
        MCAmericanBasketEngine(const ext::shared_ptr<StochasticProcessArray>& process,
                               intOrNull timeSteps = Null<Size>(),
                               intOrNull timeStepsPerYear = Null<Size>(),
                               bool brownianBridge = false,
                               bool antitheticVariate = false,
                               intOrNull requiredSamples = Null<Size>(),
                               doubleOrNull requiredTolerance = Null<Real>(),
                               intOrNull maxSamples = Null<Size>(),
                               BigInteger seed = 0,
                               Size nCalibrationSamples = Null<Size>(),
                               Size polynomOrder = 2,
                               LsmBasisSystem::PolynomialType polynomType 
                                  	= LsmBasisSystem::Monomial) {
            return new MCAmericanBasketEngine<RNG>(process,
                                                   timeSteps,
                                                   timeStepsPerYear,
                                                   brownianBridge,
                                                   antitheticVariate,
                                                   requiredSamples,
                                                   requiredTolerance,
                                                   maxSamples,
                                                   seed,
                                                   nCalibrationSamples,
                                                   polynomOrder,
                                                   polynomType);
        }
    }
};

%template(MCPRAmericanBasketEngine) MCAmericanBasketEngine<PseudoRandom>;
%template(MCLDAmericanBasketEngine) MCAmericanBasketEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCAmericanBasketEngine(process,
                               traits,
                               timeSteps=None,
                               timeStepsPerYear=None,
                               brownianBridge=False,
                               antitheticVariate=False,
                               requiredSamples=None,
                               requiredTolerance=None,
                               maxSamples=None,
                               seed=0,
                               nCalibrationSamples=2048,
                               polynomOrder=2,
                               polynomType=LsmBasisSystem.Monomial):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPRAmericanBasketEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDAmericanBasketEngine
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
                   seed,
                   nCalibrationSamples,
                   polynomOrder,
                   polynomType)
%}
#endif

%{
using QuantLib::StulzEngine;
using QuantLib::KirkEngine;
using QuantLib::Fd2dBlackScholesVanillaEngine;
%}

%shared_ptr(StulzEngine)
class StulzEngine : public PricingEngine {
  public:
    StulzEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process1,
                const ext::shared_ptr<GeneralizedBlackScholesProcess>& process2,
                Real correlation);
};

%shared_ptr(KirkEngine)
class KirkEngine : public PricingEngine {
  public:
    KirkEngine(const ext::shared_ptr<BlackProcess>& process1,
               const ext::shared_ptr<BlackProcess>& process2,
               Real correlation);
};

%shared_ptr(Fd2dBlackScholesVanillaEngine)
class Fd2dBlackScholesVanillaEngine : public PricingEngine {
  public:
    Fd2dBlackScholesVanillaEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& p1,
                const ext::shared_ptr<GeneralizedBlackScholesProcess>& p2,
                Real correlation,
                Size xGrid = 100, Size yGrid = 100, 
                Size tGrid = 50, Size dampingSteps = 0,
                const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
                bool localVol = false,
                Real illegalLocalVolOverwrite = -Null<Real>());
};

%{
using QuantLib::EverestOption;
using QuantLib::MCEverestEngine;
%}

%shared_ptr(EverestOption)
class EverestOption : public MultiAssetOption {
  public:
    EverestOption(Real notional,
                     Rate guarantee,
                     const ext::shared_ptr<Exercise>& exercise);
};

%shared_ptr(MCEverestEngine<PseudoRandom>);
%shared_ptr(MCEverestEngine<LowDiscrepancy>);

template <class RNG>
class MCEverestEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCEverestEngine;
    #endif
  public:
    %extend {
        MCEverestEngine(const ext::shared_ptr<StochasticProcessArray>& process,
                        Size timeSteps = Null<Size>(),
                        Size timeStepsPerYear = Null<Size>(),
                        bool brownianBridge = false,
                        bool antitheticVariate = false,
                        intOrNull requiredSamples = Null<Size>(),
                        doubleOrNull requiredTolerance = Null<Real>(),
                        intOrNull maxSamples = Null<Size>(),
                        BigInteger seed = 0) {
            return new MCEverestEngine<RNG>(process,
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

%template(MCPREverestEngine) MCEverestEngine<PseudoRandom>;
%template(MCLDEverestEngine) MCEverestEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCEverestEngine(process,
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
            cls = MCPREverestEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDEverestEngine
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


%{
using QuantLib::HimalayaOption;
using QuantLib::MCHimalayaEngine;
%}

%shared_ptr(HimalayaOption)
class HimalayaOption : public MultiAssetOption {
  public:
    HimalayaOption(const std::vector<Date>& fixingDates,
                      Real strike);
};

%shared_ptr(MCHimalayaEngine<PseudoRandom>);
%shared_ptr(MCHimalayaEngine<LowDiscrepancy>);

template <class RNG>
class MCHimalayaEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCHimalayaEngine;
    #endif
  public:
    %extend {
        MCHimalayaEngine(const ext::shared_ptr<StochasticProcessArray>& process,
                         bool brownianBridge = false,
                         bool antitheticVariate = false,
                         intOrNull requiredSamples = Null<Size>(),
                         doubleOrNull requiredTolerance = Null<Real>(),
                         intOrNull maxSamples = Null<Size>(),
                         BigInteger seed = 0) {
            return new MCHimalayaEngine<RNG>(process,
                                             brownianBridge,
                                             antitheticVariate,
                                             requiredSamples,
                                             requiredTolerance,
                                             maxSamples,
                                             seed);
        }
    }
};

%template(MCPRHimalayaEngine) MCHimalayaEngine<PseudoRandom>;
%template(MCLDHimalayaEngine) MCHimalayaEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCHimalayaEngine(process,
                         traits,
                         brownianBridge=False,
                         antitheticVariate=False,
                         requiredSamples=None,
                         requiredTolerance=None,
                         maxSamples=None,
                         seed=0):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPRHimalayaEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDHimalayaEngine
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


#endif
