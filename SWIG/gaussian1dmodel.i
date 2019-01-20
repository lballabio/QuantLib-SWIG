
/*
 Copyright (C) 2014 Matthias Groncki
 Copyright (C) 2017, 2018 Matthias Lungwitz

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
%include calibrationhelpers.i
%include observer.i

%{
using QuantLib::Gaussian1dModel;
%}

%shared_ptr(Gaussian1dModel)
class Gaussian1dModel : public TermStructureConsistentModel {
    public:
        const boost::shared_ptr<StochasticProcess1D> stateProcess() const;

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

%{
using QuantLib::Gsr;
using QuantLib::MarkovFunctional;
%}


%shared_ptr(Gsr)
class Gsr : public Gaussian1dModel {
  public:
    Gsr(const Handle<YieldTermStructure> &termStructure,
           const std::vector<Date> &volstepdates,
           const std::vector<Handle<Quote> > &volatilities,
           const std::vector<Handle<Quote> > &reversions, const Real T = 60.0);
    
    void calibrateVolatilitiesIterative(
            const std::vector<boost::shared_ptr<BlackCalibrationHelper> > &helpers,
            OptimizationMethod &method, const EndCriteria &endCriteria,
            const Constraint &constraint = Constraint(),
            const std::vector<Real> &weights = std::vector<Real>());

    const Array &reversion() const;

    const Array &volatility() const;

    // Calibrated Model functions
    %extend{
        Array params() const { return self->params();}
        void calibrate(
            const std::vector<boost::shared_ptr<CalibrationHelperBase> >& instruments,
            OptimizationMethod& method, const EndCriteria& endCriteria,
            const Constraint& constraint = Constraint(),
            const std::vector<Real>& weights = std::vector<Real>(),
            const std::vector<bool>& fixParameters = std::vector<bool>()) { self->calibrate(instruments, method, endCriteria, constraint, weights, fixParameters);}

        void setParams(const Array& params) {self->setParams(params);}
        Real value(const Array& params,
                   const std::vector<boost::shared_ptr<CalibrationHelperBase> >& instruments) {return self->value(params, instruments);}
        const boost::shared_ptr<Constraint>& constraint() const {return self->constraint();}
        EndCriteria::Type endCriteria() const {return self->endCriteria();}
        const Array& problemValues() const {return self->problemValues();}
        Integer functionEvaluation() const {return self->functionEvaluation();}
    }
};


%rename (MarkovFunctionalSettings) MarkovFunctional::ModelSettings;
%feature ("flatnested") ModelSettings;

%shared_ptr(MarkovFunctional)
class MarkovFunctional : public Gaussian1dModel {
  public:
  
    struct ModelSettings {
       enum Adjustments {
            AdjustNone = 0,
            AdjustDigitals = 1 << 0,
            AdjustYts = 1 << 1,
            ExtrapolatePayoffFlat = 1 << 2,
            NoPayoffExtrapolation = 1 << 3,
            KahaleSmile = 1 << 4,
            SmileExponentialExtrapolation = 1 << 5,
            KahaleInterpolation = 1 << 6,
            SmileDeleteArbitragePoints = 1 << 7,
            SabrSmile = 1 << 8
        };

        ModelSettings();
        
        ModelSettings(Size yGridPoints, Real yStdDevs, Size gaussHermitePoints,
                      Real digitalGap, Real marketRateAccuracy,
                      Real lowerRateBound, Real upperRateBound,
                      int adjustments,
                      const std::vector<Real>& smileMoneyCheckpoints = std::vector<Real>());
        void validate();
    };

    // Constructor for a swaption smile calibrated model
    MarkovFunctional(const Handle<YieldTermStructure> &termStructure,
                     const Real reversion,
                     const std::vector<Date> &volstepdates,
                     const std::vector<Real> &volatilities,
                     const Handle<SwaptionVolatilityStructure> &swaptionVol,
                     const std::vector<Date> &swaptionExpiries,
                     const std::vector<Period> &swaptionTenors,
                     const boost::shared_ptr<SwapIndex> &swapIndexBase,
                     const MarkovFunctional::ModelSettings &modelSettings =
                         ModelSettings());
    
    // Constructor for a caplet smile calibrated model
    MarkovFunctional(const Handle<YieldTermStructure> &termStructure,
                     const Real reversion,
                     const std::vector<Date> &volstepdates,
                     const std::vector<Real> &volatilities,
                     const Handle<OptionletVolatilityStructure> &capletVol,
                     const std::vector<Date> &capletExpiries,
                     const boost::shared_ptr<IborIndex> &iborIndex,
                     const MarkovFunctional::ModelSettings &modelSettings =
                         ModelSettings());
                         
    const Array &volatility();

    void calibrate(
        const std::vector<boost::shared_ptr<BlackCalibrationHelper> > &helper,
        OptimizationMethod &method, const EndCriteria &endCriteria,
        const Constraint &constraint = Constraint(),
        const std::vector<Real> &weights = std::vector<Real>(),
        const std::vector<bool> &fixParameters = std::vector<bool>());  

    //  Calibrated Model functions
    %extend{
        Array params() const { return self->params();}
        void setParams(const Array& params) {self->setParams(params);}
        Real value(const Array& params,
                   const std::vector<boost::shared_ptr<CalibrationHelperBase> >& instruments) {return self->value(params, instruments);}
        const boost::shared_ptr<Constraint>& constraint() const {return self->constraint();}
        EndCriteria::Type endCriteria() const {return self->endCriteria();}
        const Array& problemValues() const {return self->problemValues();}
        Integer functionEvaluation() const {return self->functionEvaluation();}
    }
    
};

// Pricing Engines

%{
using QuantLib::Gaussian1dSwaptionEngine;
using QuantLib::Gaussian1dJamshidianSwaptionEngine;
using QuantLib::Gaussian1dNonstandardSwaptionEngine;
using QuantLib::Gaussian1dFloatFloatSwaptionEngine;
 
typedef boost::shared_ptr<PricingEngine> Gaussian1dSwaptionEnginePtr;
typedef boost::shared_ptr<PricingEngine> Gaussian1dJamshidianSwaptionEnginePtr;
typedef boost::shared_ptr<PricingEngine> Gaussian1dNonstandardSwaptionEnginePtr;
typedef boost::shared_ptr<PricingEngine> Gaussian1dFloatFloatSwaptionEnginePtr;
%}

#if defined(SWIGJAVA) || defined(SWIGCSHARP)
%rename(_Gaussian1dSwaptionEngine) Gaussian1dSwaptionEngine;
#else
%ignore Gaussian1dSwaptionEngine;
#endif
class Gaussian1dSwaptionEngine {
  public:
    enum Probabilities { None, Naive, Digital };
#if defined(SWIGJAVA) || defined(SWIGCSHARP)
  private:
    Gaussian1dSwaptionEngine();
#endif
};

%rename(Gaussian1dSwaptionEngine) Gaussian1dSwaptionEnginePtr;
class Gaussian1dSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
	static const Gaussian1dSwaptionEngine::Probabilities None = Gaussian1dSwaptionEngine::None;
	static const Gaussian1dSwaptionEngine::Probabilities Naive = Gaussian1dSwaptionEngine::Naive;
	static const Gaussian1dSwaptionEngine::Probabilities Digital = Gaussian1dSwaptionEngine::Digital;

    Gaussian1dSwaptionEnginePtr(const boost::shared_ptr<Gaussian1dModel> &model,
            const int integrationPoints = 64, const Real stddevs = 7.0,
            const bool extrapolatePayoff = true,
            const bool flatPayoffExtrapolation = false,
            const Handle<YieldTermStructure> &discountCurve =
                Handle<YieldTermStructure>(),
			const Gaussian1dSwaptionEngine::Probabilities probabilities = Gaussian1dSwaptionEngine::None) {
            return new Gaussian1dSwaptionEnginePtr(new Gaussian1dSwaptionEngine(model, integrationPoints, 
                    stddevs, extrapolatePayoff, flatPayoffExtrapolation, discountCurve, probabilities));
        }
    
    }
};

%rename(Gaussian1dJamshidianSwaptionEngine) Gaussian1dJamshidianSwaptionEnginePtr;
class Gaussian1dJamshidianSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {

    Gaussian1dJamshidianSwaptionEnginePtr(const boost::shared_ptr<Gaussian1dModel> &model) {
            return new Gaussian1dJamshidianSwaptionEnginePtr(new Gaussian1dJamshidianSwaptionEngine(model));
        }
    
    }
};

#if defined(SWIGJAVA) || defined(SWIGCSHARP)
%rename(_Gaussian1dNonstandardSwaptionEngine) Gaussian1dNonstandardSwaptionEngine;
#else
%ignore Gaussian1dNonstandardSwaptionEngine;
#endif
class Gaussian1dNonstandardSwaptionEngine {
  public:
    enum Probabilities { None, Naive, Digital };
#if defined(SWIGJAVA) || defined(SWIGCSHARP)
  private:
    Gaussian1dNonstandardSwaptionEngine();
#endif
};

%rename(Gaussian1dNonstandardSwaptionEngine) Gaussian1dNonstandardSwaptionEnginePtr;
class Gaussian1dNonstandardSwaptionEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {

	static const Gaussian1dNonstandardSwaptionEngine::Probabilities None = Gaussian1dNonstandardSwaptionEngine::None;
	static const Gaussian1dNonstandardSwaptionEngine::Probabilities Naive = Gaussian1dNonstandardSwaptionEngine::Naive;
	static const Gaussian1dNonstandardSwaptionEngine::Probabilities Digital = Gaussian1dNonstandardSwaptionEngine::Digital;

    Gaussian1dNonstandardSwaptionEnginePtr(
            const boost::shared_ptr<Gaussian1dModel> &model,
            const int integrationPoints = 64, const Real stddevs = 7.0,
            const bool extrapolatePayoff = true,
            const bool flatPayoffExtrapolation = false,
            const Handle<Quote> &oas = Handle<Quote>(), // continuously
                                                        // compounded w.r.t. yts
                                                        // daycounter
            const Handle<YieldTermStructure> &discountCurve =
                Handle<YieldTermStructure>(),
				const Gaussian1dNonstandardSwaptionEngine::Probabilities probabilities = Gaussian1dNonstandardSwaptionEngine::None) {
            return new Gaussian1dNonstandardSwaptionEnginePtr(new Gaussian1dNonstandardSwaptionEngine(model, integrationPoints, 
                    stddevs, extrapolatePayoff, flatPayoffExtrapolation, oas, discountCurve, probabilities));
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

        static const Gaussian1dFloatFloatSwaptionEngine::Probabilities None
            = Gaussian1dFloatFloatSwaptionEngine::None;
        static const Gaussian1dFloatFloatSwaptionEngine::Probabilities Naive
            = Gaussian1dFloatFloatSwaptionEngine::Naive;
        static const Gaussian1dFloatFloatSwaptionEngine::Probabilities Digital
            = Gaussian1dFloatFloatSwaptionEngine::Digital;

        Gaussian1dFloatFloatSwaptionEnginePtr(
                const boost::shared_ptr<Gaussian1dModel> &model,
                const int integrationPoints = 64, const Real stddevs = 7.0,
                const bool extrapolatePayoff = true,
                const bool flatPayoffExtrapolation = false,
                const Handle<Quote> &oas = Handle<Quote>(),
                const Handle<YieldTermStructure> &discountCurve =
                                                 Handle<YieldTermStructure>(),
                const bool includeTodaysExercise = false,
                const Gaussian1dFloatFloatSwaptionEngine::Probabilities probabilities = Gaussian1dFloatFloatSwaptionEngine::None) {
            return new Gaussian1dFloatFloatSwaptionEnginePtr(
                new Gaussian1dFloatFloatSwaptionEngine(model,
                                                       integrationPoints,
                                                       stddevs,
                                                       extrapolatePayoff,
                                                       flatPayoffExtrapolation,
                                                       oas,
                                                       discountCurve,
                                                       includeTodaysExercise,
                                                       probabilities));
        }

    }

};

#endif
