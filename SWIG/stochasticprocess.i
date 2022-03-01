
/*
 Copyright (C) 2004, 2005, 2007, 2008 StatPro Italia srl
 Copyright (C) 2010 Klaus Spanderen
 Copyright (C) 2015 Matthias Groncki
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

#ifndef quantlib_stochastic_process_i
#define quantlib_stochastic_process_i

%include marketelements.i
%include termstructures.i
%include volatilities.i
%include observer.i

%{
using QuantLib::StochasticProcess;
%}

%shared_ptr(StochasticProcess)
class StochasticProcess : public Observable {
  private:
    StochasticProcess();
  public:
    Size size() const;
    Size factors() const;
    Array initialValues() const;
    Array drift(Time t, const Array& x) const;
    Matrix diffusion(Time t, const Array& x) const;
    Array expectation(Time t0, const Array& x0, Time dt) const;
    Matrix stdDeviation(Time t0, const Array& x0, Time dt) const;
    Matrix covariance(Time t0, const Array& x0, Time dt) const;
    Array evolve(Time t0, const Array& x0,
                 Time dt, const Array& dw) const;
};

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<StochasticProcess> )
#endif
%template(StochasticProcessVector)
std::vector<ext::shared_ptr<StochasticProcess> >;


%{
using QuantLib::StochasticProcess1D;
%}

%shared_ptr(StochasticProcess1D)
class StochasticProcess1D
    : public StochasticProcess {
  public:
      Real x0() const;
      Real drift(Time t, Real x) const;
      Real diffusion(Time t, Real x) const;
      Real expectation(Time t0, Real x0, Time dt) const;
      Real stdDeviation(Time t0, Real x0, Time dt) const;
      Real variance(Time t0, Real x0, Time dt) const;
      Real evolve(Time t0, Real x0, Time dt, Real dw) const;
      Real apply(Real x0, Real dx) const;
};

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<StochasticProcess1D> )
#endif
%template(StochasticProcess1DVector)
std::vector<ext::shared_ptr<StochasticProcess1D> >;


%{
using QuantLib::GeneralizedBlackScholesProcess;
%}

%shared_ptr(GeneralizedBlackScholesProcess)
class GeneralizedBlackScholesProcess : public StochasticProcess1D {
  public:
      GeneralizedBlackScholesProcess(
                             const Handle<Quote>& s0,
                             const Handle<YieldTermStructure>& dividendTS,
                             const Handle<YieldTermStructure>& riskFreeTS,
                             const Handle<BlackVolTermStructure>& volTS);

      GeneralizedBlackScholesProcess(
            const Handle<Quote>& x0,
            const Handle<YieldTermStructure>& dividendTS,
            const Handle<YieldTermStructure>& riskFreeTS,
            const Handle<BlackVolTermStructure>& blackVolTS,
            const Handle<LocalVolTermStructure>& localVolTS);

      Handle<Quote> stateVariable();
      Handle<YieldTermStructure> dividendYield();
      Handle<YieldTermStructure> riskFreeRate();
      Handle<BlackVolTermStructure> blackVolatility();
      Handle<LocalVolTermStructure> localVolatility();
};

%{
using QuantLib::BlackScholesProcess;
%}

%shared_ptr(BlackScholesProcess)
class BlackScholesProcess : public GeneralizedBlackScholesProcess {
  public:
  BlackScholesProcess(const Handle<Quote>& s0,
                           const Handle<YieldTermStructure>& riskFreeTS,
                           const Handle<BlackVolTermStructure>& volTS);
};

%{
using QuantLib::BlackScholesMertonProcess;
%}

%shared_ptr(BlackScholesMertonProcess)
class BlackScholesMertonProcess : public GeneralizedBlackScholesProcess {
  public:
      BlackScholesMertonProcess(
                             const Handle<Quote>& s0,
                             const Handle<YieldTermStructure>& dividendTS,
                             const Handle<YieldTermStructure>& riskFreeTS,
                             const Handle<BlackVolTermStructure>& volTS);
};

%{
using QuantLib::BlackProcess;
%}

%shared_ptr(BlackProcess)
class BlackProcess : public GeneralizedBlackScholesProcess {
  public:
      BlackProcess(const Handle<Quote>& s0,
                      const Handle<YieldTermStructure>& riskFreeTS,
                      const Handle<BlackVolTermStructure>& volTS);
};

%{
using QuantLib::GarmanKohlagenProcess;
%}

%shared_ptr(GarmanKohlagenProcess)
class GarmanKohlagenProcess : public GeneralizedBlackScholesProcess {
  public:
      GarmanKohlagenProcess(
                         const Handle<Quote>& s0,
                         const Handle<YieldTermStructure>& foreignRiskFreeTS,
                         const Handle<YieldTermStructure>& domesticRiskFreeTS,
                         const Handle<BlackVolTermStructure>& volTS);
};



%{
using QuantLib::Merton76Process;
%}

%shared_ptr(Merton76Process)
class Merton76Process : public StochasticProcess1D {
  public:
      Merton76Process(const Handle<Quote>& stateVariable,
                         const Handle<YieldTermStructure>& dividendTS,
                         const Handle<YieldTermStructure>& riskFreeTS,
                         const Handle<BlackVolTermStructure>& volTS,
                         const Handle<Quote>& jumpIntensity,
                         const Handle<Quote>& meanLogJump,
                         const Handle<Quote>& jumpVolatility);
};

%{
using QuantLib::StochasticProcessArray;
%}

%shared_ptr(StochasticProcessArray)
class StochasticProcessArray : public StochasticProcess {
  public:
      StochasticProcessArray(
               const std::vector<ext::shared_ptr<StochasticProcess1D> >&array,
               const Matrix &correlation);
};


%{
using QuantLib::GeometricBrownianMotionProcess;
%}

%shared_ptr(GeometricBrownianMotionProcess)
class GeometricBrownianMotionProcess : public StochasticProcess1D {
  public:
      GeometricBrownianMotionProcess(Real initialValue,
                                        Real mu,
                                        Real sigma);
};

%{
using QuantLib::VarianceGammaProcess;
%}

%shared_ptr(VarianceGammaProcess)
class VarianceGammaProcess : public StochasticProcess1D {
  public:
      VarianceGammaProcess(const Handle<Quote>& s0,
            const Handle<YieldTermStructure>& dividendYield,
            const Handle<YieldTermStructure>& riskFreeRate,
            Real sigma, Real nu, Real theta);
};


%{
using QuantLib::HestonProcess;
%}

%shared_ptr(HestonProcess)
class HestonProcess : public StochasticProcess {
  public:
        enum Discretization { PartialTruncation,
                    FullTruncation,
                    Reflection,
                    NonCentralChiSquareVariance,
                    QuadraticExponential,
                    QuadraticExponentialMartingale,
                    BroadieKayaExactSchemeLobatto,
                    BroadieKayaExactSchemeLaguerre,
                    BroadieKayaExactSchemeTrapezoidal };

        HestonProcess(const Handle<YieldTermStructure>& riskFreeTS,
					   const Handle<YieldTermStructure>& dividendTS,
					   const Handle<Quote>& s0,
					   Real v0, Real kappa,
                       Real theta, Real sigma, Real rho,
                       Discretization d = QuadraticExponentialMartingale);

        Handle<Quote> s0();
      Handle<YieldTermStructure> dividendYield();
      Handle<YieldTermStructure> riskFreeRate();
};

%{
using QuantLib::BatesProcess;
%}

%shared_ptr(BatesProcess)
class BatesProcess : public HestonProcess {
  public:
      BatesProcess(const Handle<YieldTermStructure>& riskFreeRate,
                      const Handle<YieldTermStructure>& dividendYield,
                      const Handle<Quote>& s0,
                      Real v0, Real kappa,
                      Real theta, Real sigma, Real rho,
                      Real lambda, Real nu, Real delta);
};

%{
using QuantLib::HullWhiteProcess;
%}

%shared_ptr(HullWhiteProcess)
class HullWhiteProcess : public StochasticProcess1D {
  public:
      HullWhiteProcess(const Handle<YieldTermStructure>& riskFreeTS,
                          Real a, Real sigma);
};

%{
using QuantLib::HullWhiteForwardProcess;
%}

%shared_ptr(HullWhiteForwardProcess)
class HullWhiteForwardProcess : public StochasticProcess1D {
  public:
    HullWhiteForwardProcess(const Handle<YieldTermStructure>& riskFreeTS,
                             Real a,
                             Real sigma);
    Real alpha(Time t) const;
    Real M_T(Real s, Real t, Real T) const;
    Real B(Time t, Time T) const;
    void setForwardMeasureTime(Time t);
};

%{
using QuantLib::G2Process;
%}

%shared_ptr(G2Process)
class G2Process : public StochasticProcess {
  public:
    G2Process(Real a, Real sigma, Real b, Real eta, Real rho);
};

%{
using QuantLib::G2ForwardProcess;
%}

%shared_ptr(G2ForwardProcess)
class G2ForwardProcess : public StochasticProcess {
  public:
    G2ForwardProcess(Real a, Real sigma, Real b, Real eta, Real rho);
    void setForwardMeasureTime(Time t);
};

%{
using QuantLib::GsrProcess;
%}

%shared_ptr(GsrProcess)
class GsrProcess : public StochasticProcess1D {
    public:
    GsrProcess(const Array &times, const Array &vols,
               const Array &reversions, const Real T = 60.0);
    Real sigma(Time t);
    Real reversion(Time t);
    Real y(Time t);
    Real G(Time t, Time T, Real x);
    void setForwardMeasureTime(Time t);
};

%inline %{
    const ext::shared_ptr<GsrProcess> as_gsr_process(
                           const ext::shared_ptr<StochasticProcess>& proc) {
        return ext::dynamic_pointer_cast<GsrProcess>(proc);
    }
%}


%{
using QuantLib::OrnsteinUhlenbeckProcess;
%}

%shared_ptr(OrnsteinUhlenbeckProcess)
class OrnsteinUhlenbeckProcess : public StochasticProcess1D {
  public:
    OrnsteinUhlenbeckProcess(
    	Real speed, Volatility vol, Real x0 = 0.0, Real level = 0.0);
    	    	
    Real speed() const;
    Real volatility() const;
    Real level() const;
};


%{
using QuantLib::KlugeExtOUProcess;
using QuantLib::ExtendedOrnsteinUhlenbeckProcess;
using QuantLib::ExtOUWithJumpsProcess;
%}

%shared_ptr(ExtendedOrnsteinUhlenbeckProcess)
class ExtendedOrnsteinUhlenbeckProcess : public StochasticProcess1D {
    public:
        enum Discretization { MidPoint, Trapezodial, GaussLobatto };

        ExtendedOrnsteinUhlenbeckProcess(
                                Real speed, Volatility sigma, Real x0,
                                const ext::function<Real (Real)>& b,
                                Discretization discretization = MidPoint,
                                Real intEps = 1e-4);
    %extend{                            
        #if defined(SWIGPYTHON)    
        ExtendedOrnsteinUhlenbeckProcess(
            Real speed, Volatility sigma, Real x0, 
            PyObject* function, 
            Real intEps = 1e-4) {
            
            const UnaryFunction f(function);
            return new ExtendedOrnsteinUhlenbeckProcess(
            	    speed, sigma, x0, f, 
            	    ExtendedOrnsteinUhlenbeckProcess::MidPoint, intEps);
        }
        #elif defined(SWIGJAVA) || defined(SWIGCSHARP)
        ExtendedOrnsteinUhlenbeckProcess(
            Real speed, Volatility sigma, Real x0, 
            UnaryFunctionDelegate* function,
            Real intEps = 1e-4) {
            
            const UnaryFunction f(function);
            return new ExtendedOrnsteinUhlenbeckProcess(
            	    speed, sigma, x0, f, 
            	    ExtendedOrnsteinUhlenbeckProcess::MidPoint, intEps);
        }
		#endif
    }
};

%shared_ptr(ExtOUWithJumpsProcess)
class ExtOUWithJumpsProcess : public StochasticProcess {
  public:
    ExtOUWithJumpsProcess(
            const ext::shared_ptr<ExtendedOrnsteinUhlenbeckProcess>& process,
            Real Y0, Real beta, Real jumpIntensity, Real eta);
};

%shared_ptr(KlugeExtOUProcess)
class KlugeExtOUProcess : public StochasticProcess {
  public:
    KlugeExtOUProcess(
            Real rho,
            const ext::shared_ptr<ExtOUWithJumpsProcess>& kluge,
            const ext::shared_ptr<ExtendedOrnsteinUhlenbeckProcess>& extOU);
};

%{
using QuantLib::GJRGARCHProcess;
%}

%shared_ptr(GJRGARCHProcess)
class GJRGARCHProcess : public StochasticProcess {
  public:
      enum Discretization { PartialTruncation, FullTruncation, Reflection};

      GJRGARCHProcess(const Handle<YieldTermStructure>& riskFreeRate,
                      const Handle<YieldTermStructure>& dividendYield,
                      const Handle<Quote>& s0,
                      Real v0, Real omega, Real alpha, Real beta,
                      Real gamma, Real lambda, Real daysPerYear = 252.0,
                      Discretization d = FullTruncation);

  Handle<Quote> s0();
  Handle<YieldTermStructure> dividendYield();
  Handle<YieldTermStructure> riskFreeRate();
};



#endif
