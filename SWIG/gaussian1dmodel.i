
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
        const ext::shared_ptr<StochasticProcess1D> stateProcess() const;

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
                               ext::shared_ptr<IborIndex> iborIdx =
                                   ext::shared_ptr<IborIndex>()) const;

        const Real swapRate(const Date &fixing, const Period &tenor,
                            const Date &referenceDate = Null<Date>(),
                            const Real y = 0.0,
                            ext::shared_ptr<SwapIndex> swapIdx =
                                ext::shared_ptr<SwapIndex>()) const;

        const Real swapAnnuity(const Date &fixing, const Period &tenor,
                               const Date &referenceDate = Null<Date>(),
                               const Real y = 0.0,
                               ext::shared_ptr<SwapIndex> swapIdx =
                                   ext::shared_ptr<SwapIndex>()) const;
};

%{
using QuantLib::Gsr;
using QuantLib::MarkovFunctional;
%}


%shared_ptr(Gsr)
class Gsr : public Gaussian1dModel {
    #if defined(SWIGCSHARP)
    %rename("parameters") params;
    #endif
  public:
    Gsr(const Handle<YieldTermStructure> &termStructure,
           const std::vector<Date> &volstepdates,
           const std::vector<Handle<Quote> > &volatilities,
           const std::vector<Handle<Quote> > &reversions, const Real T = 60.0);
    
    void calibrateVolatilitiesIterative(
            const std::vector<ext::shared_ptr<BlackCalibrationHelper> > &helpers,
            OptimizationMethod &method, const EndCriteria &endCriteria,
            const Constraint &constraint = Constraint(),
            const std::vector<Real> &weights = std::vector<Real>());

    const Array &reversion() const;

    const Array &volatility() const;

    // Calibrated Model functions
    Array params() const;
    void calibrate(
            const std::vector<ext::shared_ptr<CalibrationHelper> >& instruments,
            OptimizationMethod& method, const EndCriteria& endCriteria,
            const Constraint& constraint = Constraint(),
            const std::vector<Real>& weights = std::vector<Real>(),
            const std::vector<bool>& fixParameters = std::vector<bool>());
    void setParams(const Array& params);
    Real value(const Array& params,
               const std::vector<ext::shared_ptr<CalibrationHelper> >& instruments);
    const ext::shared_ptr<Constraint>& constraint() const;
    EndCriteria::Type endCriteria() const;
    const Array& problemValues() const;
    Integer functionEvaluation() const;
};


%rename (MarkovFunctionalSettings) MarkovFunctional::ModelSettings;
%feature ("flatnested") ModelSettings;

%shared_ptr(MarkovFunctional)
class MarkovFunctional : public Gaussian1dModel {
    #if defined(SWIGCSHARP)
    %rename("parameters") params;
    #endif
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
                     const ext::shared_ptr<SwapIndex> &swapIndexBase,
                     const MarkovFunctional::ModelSettings &modelSettings =
                         ModelSettings());
    
    // Constructor for a caplet smile calibrated model
    MarkovFunctional(const Handle<YieldTermStructure> &termStructure,
                     const Real reversion,
                     const std::vector<Date> &volstepdates,
                     const std::vector<Real> &volatilities,
                     const Handle<OptionletVolatilityStructure> &capletVol,
                     const std::vector<Date> &capletExpiries,
                     const ext::shared_ptr<IborIndex> &iborIndex,
                     const MarkovFunctional::ModelSettings &modelSettings =
                         ModelSettings());
                         
    const Array &volatility();

    void calibrate(
        const std::vector<ext::shared_ptr<CalibrationHelper> > &helper,
        OptimizationMethod &method, const EndCriteria &endCriteria,
        const Constraint &constraint = Constraint(),
        const std::vector<Real> &weights = std::vector<Real>(),
        const std::vector<bool> &fixParameters = std::vector<bool>());  

    //  Calibrated Model functions
    Array params() const;
    void setParams(const Array& params);
    Real value(const Array& params,
               const std::vector<ext::shared_ptr<CalibrationHelper> >& instruments);
    const ext::shared_ptr<Constraint>& constraint() const;
    EndCriteria::Type endCriteria() const;
    const Array& problemValues() const;
    Integer functionEvaluation() const;
};

// Pricing Engines

%{
using QuantLib::Gaussian1dCapFloorEngine;
using QuantLib::Gaussian1dSwaptionEngine;
using QuantLib::Gaussian1dJamshidianSwaptionEngine;
using QuantLib::Gaussian1dNonstandardSwaptionEngine;
using QuantLib::Gaussian1dFloatFloatSwaptionEngine;
%}

%shared_ptr(Gaussian1dCapFloorEngine)
class Gaussian1dCapFloorEngine : public PricingEngine {
  public:
    Gaussian1dCapFloorEngine(const ext::shared_ptr<Gaussian1dModel> &model,
                             const int integrationPoints = 64, const Real stddevs = 7.0,
                             const bool extrapolatePayoff = true,
                             const bool flatPayoffExtrapolation = false,
                             const Handle<YieldTermStructure> &discountCurve =
                                                           Handle<YieldTermStructure>());
};

%shared_ptr(Gaussian1dSwaptionEngine)
class Gaussian1dSwaptionEngine : public PricingEngine {
    #if defined(SWIGPYTHON)
    %rename(NoProb) None;
    #endif
  public:
    enum Probabilities { None, Naive, Digital };
    Gaussian1dSwaptionEngine(const ext::shared_ptr<Gaussian1dModel> &model,
                             const int integrationPoints = 64, const Real stddevs = 7.0,
                             const bool extrapolatePayoff = true,
                             const bool flatPayoffExtrapolation = false,
                             const Handle<YieldTermStructure> &discountCurve =
                                                           Handle<YieldTermStructure>(),
                             const Gaussian1dSwaptionEngine::Probabilities probabilities =
                                                        Gaussian1dSwaptionEngine::None);
};

%shared_ptr(Gaussian1dJamshidianSwaptionEngine)
class Gaussian1dJamshidianSwaptionEngine : public PricingEngine {
  public:
    Gaussian1dJamshidianSwaptionEngine(const ext::shared_ptr<Gaussian1dModel>& model);
};

%shared_ptr(Gaussian1dNonstandardSwaptionEngine)
class Gaussian1dNonstandardSwaptionEngine : public PricingEngine {
    #if defined(SWIGPYTHON)
    %rename(NoProb) None;
    #endif
  public:
    enum Probabilities { None, Naive, Digital };
    Gaussian1dNonstandardSwaptionEngine(
            const ext::shared_ptr<Gaussian1dModel> &model,
            const int integrationPoints = 64, const Real stddevs = 7.0,
            const bool extrapolatePayoff = true,
            const bool flatPayoffExtrapolation = false,
            const Handle<Quote> &oas = Handle<Quote>(), // continuously
                                                        // compounded w.r.t. yts
                                                        // daycounter
            const Handle<YieldTermStructure> &discountCurve =
                                                 Handle<YieldTermStructure>(),
            const Gaussian1dNonstandardSwaptionEngine::Probabilities probabilities =
                                   Gaussian1dNonstandardSwaptionEngine::None);
};

%shared_ptr(Gaussian1dFloatFloatSwaptionEngine)
class Gaussian1dFloatFloatSwaptionEngine : public PricingEngine {
    #if defined(SWIGPYTHON)
    %rename(NoProb) None;
    #endif
  public:
    enum Probabilities {  None,
                          Naive,
                          Digital };
    Gaussian1dFloatFloatSwaptionEngine(
                const ext::shared_ptr<Gaussian1dModel> &model,
                const int integrationPoints = 64, const Real stddevs = 7.0,
                const bool extrapolatePayoff = true,
                const bool flatPayoffExtrapolation = false,
                const Handle<Quote> &oas = Handle<Quote>(),
                const Handle<YieldTermStructure> &discountCurve =
                                                 Handle<YieldTermStructure>(),
                const bool includeTodaysExercise = false,
                const Gaussian1dFloatFloatSwaptionEngine::Probabilities probabilities =
                                    Gaussian1dFloatFloatSwaptionEngine::None);
};

#endif
