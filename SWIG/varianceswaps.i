/*
 Copyright (C) 2026 Mahimn Patel

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

#ifndef quantlib_variance_swaps_i
#define quantlib_variance_swaps_i

%include date.i
%include fra.i
%include options.i

%{
using QuantLib::VarianceSwap;
using QuantLib::ReplicatingVarianceSwapEngine;
using QuantLib::MCVarianceSwapEngine;
%}

%shared_ptr(VarianceSwap)
class VarianceSwap : public Instrument {
  public:
    VarianceSwap(Position::Type position,
                 Real strike,
                 Real notional,
                 const Date& startDate,
                 const Date& maturityDate);
    Real strike() const;
    Position::Type position() const;
    Date startDate() const;
    Date maturityDate() const;
    Real notional() const;
    Real variance() const;
};

%shared_ptr(ReplicatingVarianceSwapEngine)
class ReplicatingVarianceSwapEngine : public PricingEngine {
  public:
    ReplicatingVarianceSwapEngine(
        const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
        Real dk = 5.0,
        const std::vector<Real>& callStrikes = std::vector<Real>(),
        const std::vector<Real>& putStrikes = std::vector<Real>());
};

%shared_ptr(MCVarianceSwapEngine<PseudoRandom>);
%shared_ptr(MCVarianceSwapEngine<LowDiscrepancy>);

template <class RNG>
class MCVarianceSwapEngine : public PricingEngine {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") MCVarianceSwapEngine;
    #endif
  public:
    %extend {
        MCVarianceSwapEngine(const ext::shared_ptr<GeneralizedBlackScholesProcess>& process,
                             intOrNull timeSteps = Null<Size>(),
                             intOrNull timeStepsPerYear = Null<Size>(),
                             bool brownianBridge = false,
                             bool antitheticVariate = false,
                             intOrNull requiredSamples = Null<Size>(),
                             doubleOrNull requiredTolerance = Null<Real>(),
                             intOrNull maxSamples = Null<Size>(),
                             BigInteger seed = 0) {
            return new MCVarianceSwapEngine<RNG>(process,
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

%template(MCPRVarianceSwapEngine) MCVarianceSwapEngine<PseudoRandom>;
%template(MCLDVarianceSwapEngine) MCVarianceSwapEngine<LowDiscrepancy>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def MCVarianceSwapEngine(process,
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
            cls = MCPRVarianceSwapEngine
        elif traits == "ld" or traits == "lowdiscrepancy":
            cls = MCLDVarianceSwapEngine
        else:
            raise RuntimeError("unknown MC traits: %s" % traits)
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

#endif
