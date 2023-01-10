/*
 Copyright (C) 2022 StatPro Italia srl
 
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

#ifndef quantlib_lmm_i
#define quantlib_lmm_i

%include common.i
%include types.i
%include vectors.i
%include linearalgebra.i
%include randomnumbers.i

%{
using QuantLib::EvolutionDescription;
%}

class EvolutionDescription {
  public:
    %extend {
        EvolutionDescription(
                const std::vector<Time>& rateTimes,
                const std::vector<Time>& evolutionTimes = {},
                const std::vector<std::pair<unsigned int, unsigned int> >& relevanceRates = {}) {
            return new EvolutionDescription(rateTimes, evolutionTimes,
                                            to_vector<std::pair<Size, Size>>(relevanceRates));
        }
    }
    const std::vector<Time>& rateTimes() const;
    const std::vector<Time>& rateTaus() const;
    const std::vector<Time>& evolutionTimes() const;
    %extend {
        std::vector<unsigned int> firstAliveRate() const {
            return to_vector<unsigned int>($self->firstAliveRate());
        }
        std::vector<std::pair<unsigned int, unsigned int> > relevanceRates() const {
            return to_vector<std::pair<unsigned int, unsigned int>>($self->relevanceRates());
        }
    }
    Size numberOfRates() const;
    Size numberOfSteps() const;
};

%rename(checkCompatibility) _checkCompatibility;
%rename(isInTerminalMeasure) _isInTerminalMeasure;
%rename(isInMoneyMarketPlusMeasure) _isInMoneyMarketPlusMeasure;
%rename(isInMoneyMarketMeasure) _isInMoneyMarketMeasure;
%rename(terminalMeasure) _terminalMeasure;
%rename(moneyMarketPlusMeasure) _moneyMarketPlusMeasure;
%rename(moneyMarketMeasure) _moneyMarketMeasure;

%inline %{

    void _checkCompatibility(const EvolutionDescription& evolution,
                             const std::vector<unsigned int>& numeraires) {
        QuantLib::checkCompatibility(evolution, to_vector<Size>(numeraires));
    }

    bool _isInTerminalMeasure(const EvolutionDescription& evolution,
                              const std::vector<unsigned int>& numeraires) {
        return QuantLib::isInTerminalMeasure(evolution, to_vector<Size>(numeraires));
    }

    bool _isInMoneyMarketPlusMeasure(const EvolutionDescription& evolution,
                                     const std::vector<unsigned int>& numeraires,
                                     Size offset = 1) {
        return QuantLib::isInMoneyMarketPlusMeasure(evolution, to_vector<Size>(numeraires), offset);
    }

    bool _isInMoneyMarketMeasure(const EvolutionDescription& evolution,
                                 const std::vector<unsigned int>& numeraires) {
        return QuantLib::isInMoneyMarketMeasure(evolution, to_vector<Size>(numeraires));
    }

    std::vector<unsigned int> _terminalMeasure(const EvolutionDescription& evolution) {
        return to_vector<unsigned int>(QuantLib::terminalMeasure(evolution));
    }

    std::vector<unsigned int> _moneyMarketPlusMeasure(const EvolutionDescription& evolution,
                                                      Size offset = 1) {
        return to_vector<unsigned int>(QuantLib::moneyMarketPlusMeasure(evolution, offset));
    }

    std::vector<unsigned int> _moneyMarketMeasure(const EvolutionDescription& evolution) {
        return to_vector<unsigned int>(QuantLib::moneyMarketMeasure(evolution));
    }

%}


%{
using QuantLib::MarketModel;
using QuantLib::MarketModelFactory;
%}

%shared_ptr(MarketModel)
class MarketModel {
  public:
    const std::vector<Rate>& initialRates() const;
    const std::vector<Spread>& displacements() const;
    const EvolutionDescription& evolution() const;
    Size numberOfRates() const;
    Size numberOfFactors() const;
    Size numberOfSteps() const;
    const Matrix& pseudoRoot(Size i) const;
    const Matrix& covariance(Size i) const;
    const Matrix& totalCovariance(Size endIndex) const;
    std::vector<Volatility> timeDependentVolatility(Size i) const;
  private:
    MarketModel();
};

class MarketModelFactory {
  public:
    ext::shared_ptr<MarketModel> create(const EvolutionDescription&,
                                        Size numberOfFactors) const;
  private:
    MarketModelFactory();
};


%{
using QuantLib::PiecewiseConstantCorrelation;
using QuantLib::ExponentialForwardCorrelation;
%}

%template() std::vector<Matrix>;

%shared_ptr(PiecewiseConstantCorrelation)
class PiecewiseConstantCorrelation {
  public:
    const std::vector<Time>& times() const;
    const std::vector<Time>& rateTimes() const;
    const std::vector<Matrix>& correlations() const;
    const Matrix& correlation(Size i) const;
    Size numberOfRates() const;
  private:
    PiecewiseConstantCorrelation();
};

%shared_ptr(ExponentialForwardCorrelation)
class ExponentialForwardCorrelation : public PiecewiseConstantCorrelation {
  public:
    ExponentialForwardCorrelation(const std::vector<Time>& rateTimes,
                                  Real longTermCorr = 0.5,
                                  Real beta = 0.2,
                                  Real gamma = 1.0,
                                  std::vector<Time> times = {});
};


%{
using QuantLib::CurveState;
using QuantLib::LMMCurveState;
%}

class CurveState {
  public:
    Size numberOfRates() const;

    const std::vector<Time>& rateTimes() const;
    const std::vector<Time>& rateTaus() const;

    Real discountRatio(Size i, Size j) const;
    Rate forwardRate(Size i) const;
    Rate coterminalSwapAnnuity(Size numeraire, Size i) const;
    Rate coterminalSwapRate(Size i) const;
    Rate cmSwapAnnuity(Size numeraire, Size i, Size spanningForwards) const;
    Rate cmSwapRate(Size i, Size spanningForwards) const;

    const std::vector<Rate>& forwardRates() const;
    const std::vector<Rate>& coterminalSwapRates() const;
    const std::vector<Rate>& cmSwapRates(Size spanningForwards) const;

    Rate swapRate(Size begin, Size end) const;
  private:
    CurveState();
};

class LMMCurveState : public CurveState {
  public:
    LMMCurveState(const std::vector<Time>& rateTimes);

    void setOnForwardRates(const std::vector<Rate>& fwdRates,
                           Size firstValidIndex = 0);
    void setOnDiscountRatios(const std::vector<DiscountFactor>& discRatios,
                             Size firstValidIndex = 0);
};


%{
using QuantLib::LMMDriftCalculator;
%}

class LMMDriftCalculator {
  public:
    LMMDriftCalculator(const Matrix& pseudo,
                       const std::vector<Spread>& displacements,
                       const std::vector<Time>& taus,
                       Size numeraire,
                       Size alive);

    void compute(const LMMCurveState& cs, std::vector<Real>& drifts) const;
    void compute(const std::vector<Rate>& fwds, std::vector<Real>& drifts) const;

    void computePlain(const LMMCurveState& cs, std::vector<Real>& drifts) const;
    void computePlain(const std::vector<Rate>& fwds, std::vector<Real>& drifts) const;

    void computeReduced(const LMMCurveState& cs, std::vector<Real>& drifts) const;
    void computeReduced(const std::vector<Rate>& fwds, std::vector<Real>& drifts) const;
};


%{
using QuantLib::MarketModelEvolver;
using QuantLib::LogNormalFwdRateIpc;
%}

%shared_ptr(MarketModelEvolver)
class MarketModelEvolver {
  public:
    %extend {
        std::vector<unsigned int> numeraires() const {
            return to_vector<unsigned int>($self->numeraires());
        }
    }
    Real startNewPath();
    Real advanceStep();
    Size currentStep() const;
    const CurveState& currentState() const;
    void setInitialState(const CurveState&);
  private:
    MarketModelEvolver();
};

%shared_ptr(LogNormalFwdRateIpc)
class LogNormalFwdRateIpc : public MarketModelEvolver {
  public:
    %extend {
        LogNormalFwdRateIpc(const ext::shared_ptr<MarketModel>& model,
                            const BrownianGeneratorFactory& factory,
                            const std::vector<unsigned int>& numeraires,
                            Size initialStep = 0) {
            return new LogNormalFwdRateIpc(model, factory,
                                           to_vector<Size>(numeraires),
                                           initialStep);
        }
    }
};


%{
using QuantLib::AbcdVol;
%}

%shared_ptr(AbcdVol)
class AbcdVol : public MarketModel {
  public:
    AbcdVol(Real a,
            Real b,
            Real c,
            Real d,
            const std::vector<Real>& ks,
            const ext::shared_ptr<PiecewiseConstantCorrelation>& corr,
            const EvolutionDescription& evolution,
            Size numberOfFactors,
            const std::vector<Rate>& initialRates,
            const std::vector<Spread>& displacements);
};


%{
using QuantLib::AbcdMathFunction;
using QuantLib::AbcdFunction;
%}

class AbcdMathFunction {
  public:
    AbcdMathFunction(Real a = 0.002,
                     Real b = 0.001, 
                     Real c = 0.16,
                     Real d = 0.0005);
    AbcdMathFunction(std::vector<Real> abcd);

    Real operator()(Time t) const;

    Time maximumLocation() const;
    Real maximumValue() const;
    Real longTermValue() const;

    Real derivative(Time t) const;
    Real primitive(Time t) const;
    Real definiteIntegral(Time t1, Time t2) const;

    Real a() const;
    Real b() const;
    Real c() const;
    Real d() const;
    const std::vector<Real>& coefficients();
    const std::vector<Real>& derivativeCoefficients();
    std::vector<Real> definiteIntegralCoefficients(Time t, Time t2) const;
    std::vector<Real> definiteDerivativeCoefficients(Time t, Time t2) const;

    static void validate(Real a, Real b, Real c, Real d);
};

class AbcdFunction : public AbcdMathFunction {
  public:
    AbcdFunction(Real a = -0.06,
                 Real b =  0.17,
                 Real c =  0.54,
                 Real d =  0.17);

    Real maximumVolatility() const;
    Real shortTermVolatility() const;
    Real longTermVolatility() const;

    Real covariance(Time t, Time T, Time S) const;
    Real covariance(Time t1, Time t2, Time T, Time S) const;
    Real volatility(Time tMin, Time tMax, Time T) const;
    Real variance(Time tMin, Time tMax, Time T) const;
        
    Real instantaneousVolatility(Time t, Time T) const;
    Real instantaneousVariance(Time t, Time T) const;
    Real instantaneousCovariance(Time u, Time T, Time S) const;

    Real primitive(Time t, Time T, Time S) const;
};

#endif
