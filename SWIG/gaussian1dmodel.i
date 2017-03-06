
/*
 Copyright (C) 2014 Matthias Groncki
 Copyright (C) 2017 Matthias Lungwitz

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

#ifndef quantlib_gaussian1dmodel_i
#define quantlib_gaussian1dmodel_i

%include stochasticprocess.i
%include date.i
%include options.i
%include indexes.i
%include optimizers.i

%{
using QuantLib::Gaussian1dModel;

%}

%ignore Gaussian1dModel;
class Gaussian1dModel {
	public:
		const StochasticProcess1DPtr stateProcess() const;

		const Real numeraire(const Time t, const Real y = 0.0,
                             const Handle<YieldTermStructure> &yts =
                                 Handle<YieldTermStructure>()) const;

		const Real zerobond(const Time T, const Time t = 0.0,
                            const Real y = 0.0,
                            const Handle<YieldTermStructure> &yts =
                                Handle<YieldTermStructure>());

		 const Real numeraire(const Date &referenceDate, const Real y = 0.0,
                             const Handle<YieldTermStructure> &yts =
                                 Handle<YieldTermStructure>()) const;

        const Real zerobond(const Date &maturity,
                            const Date &referenceDate = Null<Date>(),
                            const Real y = 0.0,
                            const Handle<YieldTermStructure> &yts =
                                Handle<YieldTermStructure>()) const;

		const Real zerobondOption(
            const Option::Type &type, const Date &expiry, const Date &valueDate,
            const Date &maturity, const Rate strike,
            const Date &referenceDate = Null<Date>(), const Real y = 0.0,
            const Handle<YieldTermStructure> &yts =
                Handle<YieldTermStructure>(),
            const Real yStdDevs = 7.0, const Size yGridPoints = 64,
            const bool extrapolatePayoff = true,
            const bool flatPayoffExtrapolation = false) const;

        const Real forwardRate(const Date &fixing,
                               const Date &referenceDate = Null<Date>(),
                               const Real y = 0.0,
                               boost::shared_ptr<IborIndex> iborIdx =
                                   boost::shared_ptr<IborIndex>()) const;

        const Real swapRate(const Date &fixing, const Period &tenor,
                            const Date &referenceDate = Null<Date>(),
                            const Real y = 0.0,
                            boost::shared_ptr<SwapIndex> swapIdx =
                                boost::shared_ptr<SwapIndex>()) const;

        const Real swapAnnuity(const Date &fixing, const Period &tenor,
                               const Date &referenceDate = Null<Date>(),
                               const Real y = 0.0,
                               boost::shared_ptr<SwapIndex> swapIdx =
                                   boost::shared_ptr<SwapIndex>()) const;
};

%template(Gaussian1dModel) boost::shared_ptr<Gaussian1dModel>;

%{
using QuantLib::Gsr;
using QuantLib::MarkovFunctional;

typedef boost::shared_ptr<Gaussian1dModel> GsrPtr;
typedef boost::shared_ptr<MarkovFunctional> MarkovFunctionalPtr;
%}


%rename(Gsr) GsrPtr;
class GsrPtr : public boost::shared_ptr<Gaussian1dModel> {
  public:
    %extend {

	GsrPtr(const Handle<YieldTermStructure> &termStructure,
            const std::vector<Date> &volstepdates,
            const std::vector<Handle<Quote> > &volatilities,
            const std::vector<Handle<Quote> > &reversions, const Real T = 60.0) {
			return new GsrPtr(new Gsr(termStructure, volstepdates, volatilities, reversions, T));
		}
	
	void calibrateVolatilitiesIterative(
            const std::vector<boost::shared_ptr<CalibrationHelper> > &helpers,
            OptimizationMethod &method, const EndCriteria &endCriteria,
            const Constraint &constraint = Constraint(),
            const std::vector<Real> &weights = std::vector<Real>()) {
				boost::dynamic_pointer_cast<Gsr>(*self)->calibrateVolatilitiesIterative(helpers, method, 
						endCriteria, constraint, weights);
            }
	
	const Array &reversion() const {
		return boost::dynamic_pointer_cast<Gsr>(*self)->reversion();
	}

	const Array &volatility() const {
		return boost::dynamic_pointer_cast<Gsr>(*self)->volatility();
	}
	
    }

};


%rename(MarkovFunctional) MarkovFunctionalPtr;
class MarkovFunctionalPtr : public boost::shared_ptr<Gaussian1dModel> {
  public:
    %extend {

	MarkovFunctionalPtr(const Handle<YieldTermStructure> &termStructure,
			const Real reversion,
			const std::vector<Date> &volstepdates,
			const std::vector<Real> &volatilities,
			const Handle<SwaptionVolatilityStructure> &swaptionVol,
			const std::vector<Date> &swaptionExpiries,
			const std::vector<Period> &swaptionTenors,
			const SwapIndexPtr& swapIndexBase) {
			const boost::shared_ptr<SwapIndex> swi =
                boost::dynamic_pointer_cast<SwapIndex>(swapIndexBase);
			return new MarkovFunctionalPtr(new MarkovFunctional(termStructure, reversion, volstepdates, volatilities, swaptionVol,
			swaptionExpiries, swaptionTenors, swi));
		}

	
	void calibrate(
		const std::vector<boost::shared_ptr<CalibrationHelper> > &helper,
            OptimizationMethod &method, const EndCriteria &endCriteria,
            const Constraint &constraint = Constraint(),
            const std::vector<Real> &weights = std::vector<Real>(),
            const std::vector<bool> &fixParameters = std::vector<bool>()) {
			boost::dynamic_pointer_cast<MarkovFunctional>(*self)->calibrate(helper, method, 
					endCriteria, constraint, weights, fixParameters);
		}
		
	const Array &volatility() const {
		return boost::dynamic_pointer_cast<MarkovFunctional>(*self)->volatility();
	}
	
    }
};

// Pricing Engines

%{
using QuantLib::Gaussian1dSwaptionEngine;
using QuantLib::Gaussian1dNonstandardSwaptionEngine;
using QuantLib::Gaussian1dFloatFloatSwaptionEngine;
 
typedef boost::shared_ptr<PricingEngine> Gaussian1dSwaptionEnginePtr;
typedef boost::shared_ptr<PricingEngine> Gaussian1dNonstandardSwaptionEnginePtr;
typedef boost::shared_ptr<PricingEngine> Gaussian1dFloatFloatSwaptionEnginePtr;
%}

%rename(Gaussian1dSwaptionEngine) Gaussian1dSwaptionEnginePtr;
class Gaussian1dSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {

	Gaussian1dSwaptionEnginePtr(const boost::shared_ptr<Gaussian1dModel> &model,
            const int integrationPoints = 64, const Real stddevs = 7.0,
            const bool extrapolatePayoff = true,
            const bool flatPayoffExtrapolation = false,
            const Handle<YieldTermStructure> &discountCurve =
                Handle<YieldTermStructure>()) {
			return new Gaussian1dSwaptionEnginePtr(new Gaussian1dSwaptionEngine(model, integrationPoints, 
					stddevs, extrapolatePayoff, flatPayoffExtrapolation, discountCurve));
		}
	
    }

};

%rename(Gaussian1dNonstandardSwaptionEngine) Gaussian1dNonstandardSwaptionEnginePtr;
class Gaussian1dNonstandardSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {

	Gaussian1dNonstandardSwaptionEnginePtr(
            const boost::shared_ptr<Gaussian1dModel> &model,
            const int integrationPoints = 64, const Real stddevs = 7.0,
            const bool extrapolatePayoff = true,
            const bool flatPayoffExtrapolation = false,
            const Handle<Quote> &oas = Handle<Quote>(), // continuously
                                                        // compounded w.r.t. yts
                                                        // daycounter
            const Handle<YieldTermStructure> &discountCurve =
                Handle<YieldTermStructure>()) {
			return new Gaussian1dNonstandardSwaptionEnginePtr(new Gaussian1dNonstandardSwaptionEngine(model, integrationPoints, 
					stddevs, extrapolatePayoff, flatPayoffExtrapolation, oas, discountCurve));
		}
	
    }

};

#if defined(SWIGJAVA) || defined(SWIGCSHARP)
%rename(_Gaussian1dFloatFloatSwaptionEngine) Gaussian1dFloatFloatSwaptionEngine;
#else
%ignore Gaussian1dFloatFloatSwaptionEngine;
#endif
class Gaussian1dFloatFloatSwaptionEngine {
  public:
    enum Probabilities {  None,
            Naive,
            Digital };
#if defined(SWIGJAVA) || defined(SWIGCSHARP)
  private:
    Gaussian1dFloatFloatSwaptionEngine();
#endif
};

%rename(Gaussian1dFloatFloatSwaptionEngine) Gaussian1dFloatFloatSwaptionEnginePtr;
class Gaussian1dFloatFloatSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {

	static const Gaussian1dFloatFloatSwaptionEngine::Probabilities None = Gaussian1dFloatFloatSwaptionEngine::None;
    static const Gaussian1dFloatFloatSwaptionEngine::Probabilities Naive = Gaussian1dFloatFloatSwaptionEngine::Naive;
	static const Gaussian1dFloatFloatSwaptionEngine::Probabilities Digital = Gaussian1dFloatFloatSwaptionEngine::Digital;

	Gaussian1dFloatFloatSwaptionEnginePtr(const boost::shared_ptr<Gaussian1dModel> &model,
            const int integrationPoints = 64, const Real stddevs = 7.0,
            const bool extrapolatePayoff = true,
            const bool flatPayoffExtrapolation = false,
			const Handle<Quote> &oas = Handle<Quote>(),
            const Handle<YieldTermStructure> &discountCurve =
                Handle<YieldTermStructure>(),
			const bool includeTodaysExercise = false,
			const Gaussian1dFloatFloatSwaptionEngine::Probabilities probabilities = Gaussian1dFloatFloatSwaptionEngine::None) {
			return new Gaussian1dFloatFloatSwaptionEnginePtr(new Gaussian1dFloatFloatSwaptionEngine(model, integrationPoints,
					stddevs, extrapolatePayoff, flatPayoffExtrapolation, oas, discountCurve, includeTodaysExercise, probabilities));
		}

	}

};

#endif
