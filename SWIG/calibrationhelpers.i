/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2007, 2009 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2007 Luis Cota
 Copyright (C) 2016 Gouthaman Balaraman
 Copyright (C) 2016 Peter Caspers
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

#ifndef quantlib_calibration_helpers_i
#define quantlib_calibration_helpers_i

%include date.i
%include calendars.i
%include daycounters.i
%include cashflows.i
%include marketelements.i
%include termstructures.i
%include optimizers.i
%include options.i
%include linearalgebra.i
%include types.i
%include vectors.i

%{
using QuantLib::CalibrationHelperBase;
using QuantLib::BlackCalibrationHelper;
using QuantLib::SwaptionHelper;
using QuantLib::CapHelper;
using QuantLib::HestonModelHelper;
%}

// calibration helpers
%shared_ptr(CalibrationHelperBase)
class CalibrationHelperBase {
  public:
	Real calibrationError();
  private:
    CalibrationHelperBase();
};

%shared_ptr(BlackCalibrationHelper)
class BlackCalibrationHelper : public CalibrationHelperBase {
  public:
    enum CalibrationErrorType { RelativePriceError, PriceError, ImpliedVolError };
                       
    void setPricingEngine(const boost::shared_ptr<PricingEngine>& engine);
    Real marketValue() const;
    virtual Real modelValue() const;
    Volatility impliedVolatility(Real targetValue,
                                 Real accuracy, Size maxEvaluations,
                                 Volatility minVol, Volatility maxVol);
    Real blackPrice(Volatility volatility);
    Handle<Quote> volatility() const;
    VolatilityType volatilityType() const;
    Real calibrationError();
  private:
    BlackCalibrationHelper();

};

%inline %{
    boost::shared_ptr<BlackCalibrationHelper> as_black_helper(const boost::shared_ptr<CalibrationHelperBase>& h) {
        return boost::dynamic_pointer_cast<BlackCalibrationHelper>(h);
    }
    boost::shared_ptr<SwaptionHelper> as_swaption_helper(const boost::shared_ptr<BlackCalibrationHelper>& h) {
        return boost::dynamic_pointer_cast<SwaptionHelper>(h);
    }
%}

%shared_ptr(SwaptionHelper)
class SwaptionHelper : public BlackCalibrationHelper {
  public:
    SwaptionHelper(const Period& maturity, const Period& length,
                      const Handle<Quote>& volatility,
                      const boost::shared_ptr<IborIndex>& index,
                      const Period& fixedLegTenor,
                      const DayCounter& fixedLegDayCounter,
                      const DayCounter& floatingLegDayCounter,
                      const Handle<YieldTermStructure>& termStructure,
                      BlackCalibrationHelper::CalibrationErrorType errorType
                                = BlackCalibrationHelper::RelativePriceError,
                      const Real strike = Null<Real>(),
                      const Real nominal = 1.0,
                      const VolatilityType type = ShiftedLognormal,
                      const Real shift = 0.0);
    
    SwaptionHelper(const Date& exerciseDate, const Period& length,
                      const Handle<Quote>& volatility,
                      const boost::shared_ptr<IborIndex>& index,
                      const Period& fixedLegTenor,
                      const DayCounter& fixedLegDayCounter,
                      const DayCounter& floatingLegDayCounter,
                      const Handle<YieldTermStructure>& termStructure,
                      BlackCalibrationHelper::CalibrationErrorType errorType
                                = BlackCalibrationHelper::RelativePriceError,
                      const Real strike = Null<Real>(),
                      const Real nominal = 1.0,
                      const VolatilityType type = ShiftedLognormal,
                      const Real shift = 0.0);
    
    SwaptionHelper(const Date& exerciseDate, const Date& endDate,
                      const Handle<Quote>& volatility,
                      const boost::shared_ptr<IborIndex>& index,
                      const Period& fixedLegTenor,
                      const DayCounter& fixedLegDayCounter,
                      const DayCounter& floatingLegDayCounter,
                      const Handle<YieldTermStructure>& termStructure,
                      BlackCalibrationHelper::CalibrationErrorType errorType
                                = BlackCalibrationHelper::RelativePriceError,
                      const Real strike = Null<Real>(),
                      const Real nominal = 1.0,
                      const VolatilityType type = ShiftedLognormal,
                      const Real shift = 0.0);

    boost::shared_ptr<VanillaSwap> underlyingSwap() const;
    boost::shared_ptr<Swaption> swaption() const;
                      
    %extend {
        std::vector<Time> times() {
            std::list<Time> l;
            self->addTimesTo(l);
            std::vector<Time> v;
            std::copy(l.begin(),l.end(),std::back_inserter(v));
            return v;
        }
        Date swaptionExpiryDate() {
            return self->swaption()->exercise()->date(0);
        }
        Real swaptionStrike() {
            return self->swaption()->underlyingSwap()->fixedRate();
        }
        Real swaptionNominal() {
            return self->swaption()->underlyingSwap()->nominal();
        }
        Date swaptionMaturityDate() {
            return self->swaption()->underlyingSwap()->fixedSchedule().dates().back();
        }
    }
};

%shared_ptr(CapHelper)
class CapHelper : public BlackCalibrationHelper {
  public:
    CapHelper(const Period& length,
                 const Handle<Quote>& volatility,
                 const boost::shared_ptr<IborIndex>& index,
                 Frequency fixedLegFrequency,
                 const DayCounter& fixedLegDayCounter,
                 bool includeFirstSwaplet,
                 const Handle<YieldTermStructure>& termStructure,
                 BlackCalibrationHelper::CalibrationErrorType errorType
                                = BlackCalibrationHelper::RelativePriceError);
    %extend {        
        std::vector<Time> times() {
            std::list<Time> l;
            self->addTimesTo(l);
            std::vector<Time> v;
            std::copy(l.begin(),l.end(),std::back_inserter(v));
            return v;
        }
    }
};

%shared_ptr(HestonModelHelper)
class HestonModelHelper : public BlackCalibrationHelper {
  public:
    HestonModelHelper(const Period& maturity,
                      const Calendar& calendar,
                      const Real s0,
                      const Real strikePrice,
                      const Handle<Quote>& volatility,
                      const Handle<YieldTermStructure>& riskFreeRate,
                      const Handle<YieldTermStructure>& dividendYield,
                      BlackCalibrationHelper::CalibrationErrorType errorType
                          = BlackCalibrationHelper::RelativePriceError);
};

// allow use of vectors of helpers
#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( boost::shared_ptr<CalibrationHelperBase> )
SWIG_STD_VECTOR_ENHANCED( boost::shared_ptr<BlackCalibrationHelper> )
#endif
namespace std {
    %template(CalibrationHelperVector)
        vector<boost::shared_ptr<CalibrationHelperBase> >;
    %template(BlackCalibrationHelperVector)
        vector<boost::shared_ptr<BlackCalibrationHelper> >;        
}

// the base class for calibrated models
%{
using QuantLib::CalibratedModel;
using QuantLib::TermStructureConsistentModel;
%}

%shared_ptr(CalibratedModel)
class CalibratedModel : public virtual Observable {
    #if defined(SWIGCSHARP)
    %rename("parameters") params;
    #endif
  public:
    Array params() const;
    virtual void calibrate(
        const std::vector<boost::shared_ptr<CalibrationHelperBase> >&,
        OptimizationMethod&, const EndCriteria &,
        const Constraint& constraint = Constraint(),
        const std::vector<Real>& weights = std::vector<Real>(),
        const std::vector<bool>& fixParameters = std::vector<bool>());
     
    void setParams(const Array& params);
    Real value(const Array& params,
               const std::vector<boost::shared_ptr<CalibrationHelperBase> >&);
    const boost::shared_ptr<Constraint>& constraint() const;
    EndCriteria::Type endCriteria() const;
    const Array& problemValues() const;
    Integer functionEvaluation() const;
  private:
    CalibratedModel();
};

%shared_ptr(TermStructureConsistentModel)
class TermStructureConsistentModel : public virtual Observable{
  public:
    const Handle<YieldTermStructure>& termStructure() const;
  private:
    TermStructureConsistentModel();
};

%template(CalibratedModelHandle) Handle<CalibratedModel>;
%template(RelinkableCalibratedModelHandle)
RelinkableHandle<CalibratedModel>;

#endif
