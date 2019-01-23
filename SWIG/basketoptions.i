
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
typedef boost::shared_ptr<PricingEngine> MCEuropeanBasketEnginePtr;
%}

%rename(MCEuropeanBasketEngine) MCEuropeanBasketEnginePtr;
class MCEuropeanBasketEnginePtr : public boost::shared_ptr<PricingEngine> {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCEuropeanBasketEnginePtr;
    #endif
  public:
    %extend {
        MCEuropeanBasketEnginePtr(const boost::shared_ptr<StochasticProcessArray>& process,
                                  const std::string& traits,
                                  Size timeSteps = Null<Size>(),
                                  Size timeStepsPerYear = Null<Size>(),
                                  bool brownianBridge = false,
                                  bool antitheticVariate = false,
                                  intOrNull requiredSamples = Null<Size>(),
                                  doubleOrNull requiredTolerance = Null<Real>(),
                                  intOrNull maxSamples = Null<Size>(),
                                  BigInteger seed = 0) {

            std::string s = boost::algorithm::to_lower_copy(traits);
            if (s == "pseudorandom" || s == "pr")
                return new MCEuropeanBasketEnginePtr(
                   new MCEuropeanBasketEngine<PseudoRandom>(process,
                                                            timeSteps,
                                                            timeStepsPerYear,
                                                            brownianBridge,
                                                            antitheticVariate,
                                                            requiredSamples,
                                                            requiredTolerance,
                                                            maxSamples,
                                                            seed));
            else if (s == "lowdiscrepancy" || s == "ld")
                return new MCEuropeanBasketEnginePtr(
                   new MCEuropeanBasketEngine<LowDiscrepancy>(process,
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

%{
using QuantLib::MCAmericanBasketEngine;
typedef boost::shared_ptr<PricingEngine> MCAmericanBasketEnginePtr;
%}

%rename(MCAmericanBasketEngine) MCAmericanBasketEnginePtr;
class MCAmericanBasketEnginePtr : public boost::shared_ptr<PricingEngine> {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCAmericanBasketEnginePtr;
    #endif
  public:
    %extend {
        MCAmericanBasketEnginePtr(const boost::shared_ptr<StochasticProcessArray>& process,
                                  const std::string& traits,
                                  Size timeSteps = Null<Size>(),
                                  Size timeStepsPerYear = Null<Size>(),
                                  bool brownianBridge = false,
                                  bool antitheticVariate = false,
                                  intOrNull requiredSamples = Null<Size>(),
                                  doubleOrNull requiredTolerance = Null<Real>(),
                                  intOrNull maxSamples = Null<Size>(),
                                  BigInteger seed = 0,
                                  Size polynomOrder = 2,
                                  LsmBasisSystem::PolynomType polynomType 
                                  	= LsmBasisSystem::Monomial,
                                  Size nCalibrationSamples = Null<Size>()) {
            boost::shared_ptr<StochasticProcessArray> processes =
                 boost::dynamic_pointer_cast<StochasticProcessArray>(process);
            QL_REQUIRE(processes, "stochastic-process array required");
            std::string s = boost::algorithm::to_lower_copy(traits);
            if (s == "pseudorandom" || s == "pr")
                  return new MCAmericanBasketEnginePtr(
                  new MCAmericanBasketEngine<PseudoRandom>(process,
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
                                                           polynomType));
            else if (s == "lowdiscrepancy" || s == "ld")
                return new MCAmericanBasketEnginePtr(
                new MCAmericanBasketEngine<LowDiscrepancy>(process,
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
                                                           polynomType));
            else
                QL_FAIL("unknown Monte Carlo engine type: "+s);
        }
    }
};


%{
using QuantLib::StulzEngine;
typedef boost::shared_ptr<PricingEngine> StulzEnginePtr;
%}

%rename(StulzEngine) StulzEnginePtr;
class StulzEnginePtr
    : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        StulzEnginePtr(const boost::shared_ptr<GeneralizedBlackScholesProcess>& process1,
                       const boost::shared_ptr<GeneralizedBlackScholesProcess>& process2,
                       Real correlation) {
            return new StulzEnginePtr(
                          new StulzEngine(process1,process2,correlation));
        }
    }
};


%{
using QuantLib::EverestOption;
using QuantLib::MCEverestEngine;
typedef boost::shared_ptr<PricingEngine> MCEverestEnginePtr;
%}

%shared_ptr(EverestOption)
class EverestOption : public MultiAssetOption {
  public:
    EverestOption(Real notional,
                     Rate guarantee,
                     const boost::shared_ptr<Exercise>& exercise);
};

%rename(MCEverestEngine) MCEverestEnginePtr;
class MCEverestEnginePtr : public boost::shared_ptr<PricingEngine> {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCEverestEnginePtr;
    #endif
  public:
    %extend {
        MCEverestEnginePtr(const boost::shared_ptr<StochasticProcessArray>& process,
                           const std::string& traits,
                           Size timeSteps = Null<Size>(),
                           Size timeStepsPerYear = Null<Size>(),
                           bool brownianBridge = false,
                           bool antitheticVariate = false,
                           intOrNull requiredSamples = Null<Size>(),
                           doubleOrNull requiredTolerance = Null<Real>(),
                           intOrNull maxSamples = Null<Size>(),
                           BigInteger seed = 0) {
            std::string s = boost::algorithm::to_lower_copy(traits);
            if (s == "pseudorandom" || s == "pr")
                return new MCEverestEnginePtr(
                        new MCEverestEngine<PseudoRandom>(process,
                                                          timeSteps,
                                                          timeStepsPerYear,
                                                          brownianBridge,
                                                          antitheticVariate,
                                                          requiredSamples,
                                                          requiredTolerance,
                                                          maxSamples,
                                                          seed));
            else if (s == "lowdiscrepancy" || s == "ld")
                return new MCEverestEnginePtr(
                      new MCEverestEngine<LowDiscrepancy>(process,
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


%{
using QuantLib::HimalayaOption;
using QuantLib::MCHimalayaEngine;
typedef boost::shared_ptr<PricingEngine> MCHimalayaEnginePtr;
%}

%shared_ptr(HimalayaOption)
class HimalayaOption : public MultiAssetOption {
  public:
    HimalayaOption(const std::vector<Date>& fixingDates,
                      Real strike);
};

%rename(MCHimalayaEngine) MCHimalayaEnginePtr;
class MCHimalayaEnginePtr : public boost::shared_ptr<PricingEngine> {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCHimalayaEnginePtr;
    #endif
  public:
    %extend {
        MCHimalayaEnginePtr(const boost::shared_ptr<StochasticProcessArray>& process,
                            const std::string& traits,
                            bool brownianBridge = false,
                            bool antitheticVariate = false,
                            intOrNull requiredSamples = Null<Size>(),
                            doubleOrNull requiredTolerance = Null<Real>(),
                            intOrNull maxSamples = Null<Size>(),
                            BigInteger seed = 0) {
            std::string s = boost::algorithm::to_lower_copy(traits);
            if (s == "pseudorandom" || s == "pr")
                return new MCHimalayaEnginePtr(
                       new MCHimalayaEngine<PseudoRandom>(process,
                                                          brownianBridge,
                                                          antitheticVariate,
                                                          requiredSamples,
                                                          requiredTolerance,
                                                          maxSamples,
                                                          seed));
            else if (s == "lowdiscrepancy" || s == "ld")
                return new MCHimalayaEnginePtr(
                     new MCHimalayaEngine<LowDiscrepancy>(process,
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

#endif
