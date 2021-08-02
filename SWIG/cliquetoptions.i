/*
 Copyright (C) 2021 Jack Gillett

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

#ifndef quantlib_cliquet_options_i
#define quantlib_cliquet_options_i

%include options.i

%{
using QuantLib::CliquetOption;
%}

%shared_ptr(CliquetOption)
class CliquetOption : public OneAssetOption {
  public:
    CliquetOption(const ext::shared_ptr<PercentageStrikePayoff>& payoff,
                  const ext::shared_ptr<EuropeanExercise>& maturity,
                  std::vector<Date> resetDates);
};


%{
using QuantLib::AnalyticCliquetEngine;
using QuantLib::AnalyticPerformanceEngine;
using QuantLib::MCPerformanceEngine;
%}


%shared_ptr(AnalyticCliquetEngine)
class AnalyticCliquetEngine : public PricingEngine {
  public:
    AnalyticCliquetEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


%shared_ptr(AnalyticPerformanceEngine)
class AnalyticPerformanceEngine : public PricingEngine {
  public:
    AnalyticPerformanceEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess> process);
};


%shared_ptr(MCPerformanceEngine<PseudoRandom>);
%shared_ptr(MCPerformanceEngine<LowDiscrepancy>);

template <class RNG>
class MCPerformanceEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCPerformanceEngine;
    #endif
  public:
    %extend {
        MCPerformanceEngine(ext::shared_ptr<GeneralizedBlackScholesProcess> process,
                            bool brownianBridge = false,
                            bool antitheticVariate = false,
                            intOrNull requiredSamples = Null<Size>(),
                            doubleOrNull requiredTolerance = Null<Real>(),
                            intOrNull maxSamples = Null<Size>(),
                            BigInteger seed = 0) {
            return new MCPerformanceEngine<RNG>(process,
                                                brownianBridge,
                                                antitheticVariate,
                                                requiredSamples,
                                                requiredTolerance,
                                                maxSamples,
                                                seed);
        }
    }
};

%template(MCPRPerformanceEngine) MCPerformanceEngine<PseudoRandom>;
%template(MCLDPerformanceEngine) MCPerformanceEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCPerformanceEngine(process,
                            traits,
                            brownianBridge=False,
                            antitheticVariate=False,
                            requiredSamples=None,
                            requiredTolerance=None,
                            maxSamples=None,
                            seed=0):
        traits = traits.lower()
        if traits == "pr" or traits == "pseudorandom":
            cls = MCPRPerformanceEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDPerformanceEngine
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
