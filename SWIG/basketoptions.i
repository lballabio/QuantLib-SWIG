
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2007 Joseph Wang
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
%}

%shared_ptr(BasketPayoff)
class BasketPayoff : public Payoff {
  private:
    BasketPayoff();
};

%shared_ptr(MinBasketPayoff)
class MinBasketPayoff : public BasketPayoff  {
  public:
    MinBasketPayoff(const boost::shared_ptr<Payoff> p);
};

%shared_ptr(MaxBasketPayoff)
class MaxBasketPayoff : public BasketPayoff  {
  public:
    MaxBasketPayoff(const boost::shared_ptr<Payoff> p);
};

%shared_ptr(AverageBasketPayoff)
class AverageBasketPayoff :
      public BasketPayoff  {
  public:
    AverageBasketPayoff(const boost::shared_ptr<Payoff> p,
                           const Array &a);
    AverageBasketPayoff(const boost::shared_ptr<Payoff> p,
                           Size n);
};


%shared_ptr(BasketOption)
class BasketOption : public MultiAssetOption {
  public:
    BasketOption(
            const boost::shared_ptr<BasketPayoff>& payoff,
            const boost::shared_ptr<Exercise>& exercise);
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
        MCEuropeanBasketEngine(const boost::shared_ptr<StochasticProcessArray>& process,
                               Size timeSteps = Null<Size>(),
                               Size timeStepsPerYear = Null<Size>(),
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


%shared_ptr(MCAmericanBasketEngine<PseudoRandom>);
%shared_ptr(MCAmericanBasketEngine<LowDiscrepancy>);

template <class RNG>
class MCAmericanBasketEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCAmericanBasketEngine;
    #endif
  public:
    %extend {
        MCAmericanBasketEngine(const boost::shared_ptr<StochasticProcessArray>& process,
                               Size timeSteps = Null<Size>(),
                               Size timeStepsPerYear = Null<Size>(),
                               bool brownianBridge = false,
                               bool antitheticVariate = false,
                               intOrNull requiredSamples = Null<Size>(),
                               doubleOrNull requiredTolerance = Null<Real>(),
                               intOrNull maxSamples = Null<Size>(),
                               BigInteger seed = 0,
                               Size nCalibrationSamples = Null<Size>(),
                               Size polynomOrder = 2,
                               LsmBasisSystem::PolynomType polynomType 
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


%{
using QuantLib::StulzEngine;
%}

%shared_ptr(StulzEngine)
class StulzEngine : public PricingEngine {
  public:
    StulzEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process1,
                const boost::shared_ptr<GeneralizedBlackScholesProcess>& process2,
                Real correlation);
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
                     const boost::shared_ptr<Exercise>& exercise);
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
        MCEverestEngine(const boost::shared_ptr<StochasticProcessArray>& process,
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
        MCHimalayaEngine(const boost::shared_ptr<StochasticProcessArray>& process,
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


#endif
