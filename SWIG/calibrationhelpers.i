/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2007, 2009 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2007 Luis Cota
 Copyright (C) 2016 Gouthaman Balaraman
 Copyright (C) 2016 Peter Caspers

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
typedef boost::shared_ptr<CalibrationHelperBase> BlackCalibrationHelperPtr;
typedef boost::shared_ptr<CalibrationHelperBase> SwaptionHelperPtr;
typedef boost::shared_ptr<CalibrationHelperBase> CapHelperPtr;
typedef boost::shared_ptr<CalibrationHelperBase> HestonModelHelperPtr;
%}

// calibration helpers

%ignore CalibrationHelperBase;
class CalibrationHelperBase {
  public:
	Real calibrationError();
};

%template(CalibrationHelperBase) boost::shared_ptr<CalibrationHelperBase>;

#if defined(SWIGJAVA) || defined(SWIGCSHARP)
%rename(_BlackCalibrationHelper) BlackCalibrationHelper;
#else
%ignore BlackCalibrationHelper;
#endif
class BlackCalibrationHelper : public CalibrationHelperBase {
  private:
    BlackCalibrationHelper();
  public:
    enum CalibrationErrorType { RelativePriceError, PriceError, ImpliedVolError };
};

%rename(BlackCalibrationHelper) BlackCalibrationHelperPtr;
class BlackCalibrationHelperPtr : public boost::shared_ptr<CalibrationHelperBase> {
    #if defined(SWIGRUBY)
    %rename("pricingEngine=")      setPricingEngine;
    #endif
  private:
    BlackCalibrationHelperPtr();
  public:
    %extend {
        static const BlackCalibrationHelper::CalibrationErrorType RelativePriceError =
            BlackCalibrationHelper::RelativePriceError;
        static const BlackCalibrationHelper::CalibrationErrorType PriceError =
            BlackCalibrationHelper::PriceError;
        static const BlackCalibrationHelper::CalibrationErrorType ImpliedVolError =
            BlackCalibrationHelper::ImpliedVolError;
        void setPricingEngine(const boost::shared_ptr<PricingEngine>& engine) {
            boost::dynamic_pointer_cast<BlackCalibrationHelper>(*self)
                ->setPricingEngine(engine);
        }
        Real marketValue() const {
            return boost::dynamic_pointer_cast<BlackCalibrationHelper>(*self)
                ->marketValue();
        }
        Real modelValue() const {
            return boost::dynamic_pointer_cast<BlackCalibrationHelper>(*self)
                ->modelValue();
        }
        Volatility impliedVolatility(Real targetValue,
                                     Real accuracy, Size maxEvaluations,
                                     Volatility minVol, Volatility maxVol) {
            return boost::dynamic_pointer_cast<BlackCalibrationHelper>(*self)
                ->impliedVolatility(targetValue, accuracy, maxEvaluations,
                                    minVol, maxVol);
        }
        Real blackPrice(Volatility volatility) {
            return boost::dynamic_pointer_cast<BlackCalibrationHelper>(*self)
                ->blackPrice(volatility);
        }
        Handle<Quote> volatility() const {
            return boost::dynamic_pointer_cast<BlackCalibrationHelper>(*self)
                ->volatility();
        }
        VolatilityType volatilityType() const {
            return boost::dynamic_pointer_cast<BlackCalibrationHelper>(*self)
                ->volatilityType();
        }

        Date swaptionExpiryDate() {
            boost::shared_ptr<SwaptionHelper> s =
                boost::dynamic_pointer_cast<SwaptionHelper>(*self);
            return s ? s->swaption()->exercise()->date(0) : Null<Date>();
        }
        Real swaptionStrike() {
            boost::shared_ptr<SwaptionHelper> s =
                boost::dynamic_pointer_cast<SwaptionHelper>(*self);
            return s ? s->swaption()->underlyingSwap()->fixedRate() : Null<Real>();
        }
        Real swaptionNominal() {
            boost::shared_ptr<SwaptionHelper> s =
                boost::dynamic_pointer_cast<SwaptionHelper>(*self);
            return s ? s->swaption()->underlyingSwap()->nominal() : Null<Real>();
        }
        Date swaptionMaturityDate() {
            boost::shared_ptr<SwaptionHelper> s =
                boost::dynamic_pointer_cast<SwaptionHelper>(*self);
            return s ? s->swaption()->underlyingSwap()->fixedSchedule().dates().back() : Null<Date>();
        }
    }
};

deprecate_feature(CalibrationHelper,
                  BlackCalibrationHelper);

%inline %{
    BlackCalibrationHelperPtr as_black_helper(const boost::shared_ptr<CalibrationHelperBase>& h) {
        return boost::dynamic_pointer_cast<BlackCalibrationHelper>(h);
    }
%}

%rename(SwaptionHelper) SwaptionHelperPtr;
class SwaptionHelperPtr : public BlackCalibrationHelperPtr {
  public:
    %extend {
        SwaptionHelperPtr(const Period& maturity, const Period& length,
                          const Handle<Quote>& volatility,
                          const IborIndexPtr& index,
                          const Period& fixedLegTenor,
                          const DayCounter& fixedLegDayCounter,
                          const DayCounter& floatingLegDayCounter,
                          const Handle<YieldTermStructure>& termStructure,
                          BlackCalibrationHelper::CalibrationErrorType errorType
                                    = BlackCalibrationHelper::RelativePriceError,
                          const Real strike = Null<Real>(),
                          const Real nominal = 1.0,
                          const VolatilityType type = ShiftedLognormal,
                          const Real shift = 0.0) {
            boost::shared_ptr<IborIndex> libor =
                boost::dynamic_pointer_cast<IborIndex>(index);
            return new SwaptionHelperPtr(
                new SwaptionHelper(maturity,length,volatility,
                                   libor,fixedLegTenor,
                                   fixedLegDayCounter,
                                   floatingLegDayCounter,
                                   termStructure,
                                   errorType,
                                   strike, nominal, 
                                   type, shift));
        }
        
        SwaptionHelperPtr(const Date& exerciseDate, const Period& length,
                          const Handle<Quote>& volatility,
                          const IborIndexPtr& index,
                          const Period& fixedLegTenor,
                          const DayCounter& fixedLegDayCounter,
                          const DayCounter& floatingLegDayCounter,
                          const Handle<YieldTermStructure>& termStructure,
                          BlackCalibrationHelper::CalibrationErrorType errorType
                                    = BlackCalibrationHelper::RelativePriceError,
                          const Real strike = Null<Real>(),
                          const Real nominal = 1.0,
                          const VolatilityType type = ShiftedLognormal,
                          const Real shift = 0.0) {
            boost::shared_ptr<IborIndex> libor =
                boost::dynamic_pointer_cast<IborIndex>(index);
            return new SwaptionHelperPtr(
                new SwaptionHelper(exerciseDate,length,volatility,
                                   libor,fixedLegTenor,
                                   fixedLegDayCounter,
                                   floatingLegDayCounter,
                                   termStructure,
                                   errorType,
                                   strike, nominal, 
                                   type, shift));
        }
        
        SwaptionHelperPtr(const Date& exerciseDate, const Date& endDate,
                          const Handle<Quote>& volatility,
                          const IborIndexPtr& index,
                          const Period& fixedLegTenor,
                          const DayCounter& fixedLegDayCounter,
                          const DayCounter& floatingLegDayCounter,
                          const Handle<YieldTermStructure>& termStructure,
                          BlackCalibrationHelper::CalibrationErrorType errorType
                                    = BlackCalibrationHelper::RelativePriceError,
                          const Real strike = Null<Real>(),
                          const Real nominal = 1.0,
                          const VolatilityType type = ShiftedLognormal,
                          const Real shift = 0.0) {
            boost::shared_ptr<IborIndex> libor =
                boost::dynamic_pointer_cast<IborIndex>(index);
            return new SwaptionHelperPtr(
                new SwaptionHelper(exerciseDate,endDate,volatility,
                                   libor,fixedLegTenor,
                                   fixedLegDayCounter,
                                   floatingLegDayCounter,
                                   termStructure,
                                   errorType,
                                   strike, nominal, 
                                   type, shift));
        }
        
        std::vector<Time> times() {
            std::list<Time> l;
            boost::dynamic_pointer_cast<BlackCalibrationHelper>(*self)->addTimesTo(l);
            std::vector<Time> v;
            std::copy(l.begin(),l.end(),std::back_inserter(v));
            return v;
        }
    }
};

%rename(CapHelper) CapHelperPtr;
class CapHelperPtr : public BlackCalibrationHelperPtr {
  public:
    %extend {
        CapHelperPtr(const Period& length,
                     const Handle<Quote>& volatility,
                     const IborIndexPtr& index,
                     Frequency fixedLegFrequency,
                     const DayCounter& fixedLegDayCounter,
                     bool includeFirstSwaplet,
                     const Handle<YieldTermStructure>& termStructure,
                     BlackCalibrationHelper::CalibrationErrorType errorType
                                    = BlackCalibrationHelper::RelativePriceError) {
            boost::shared_ptr<IborIndex> libor =
                boost::dynamic_pointer_cast<IborIndex>(index);
            return new CapHelperPtr(
                new CapHelper(length,volatility,libor,fixedLegFrequency,
                              fixedLegDayCounter,includeFirstSwaplet,
                              termStructure));
        }
        std::vector<Time> times() {
            std::list<Time> l;
            boost::dynamic_pointer_cast<BlackCalibrationHelper>(*self)->addTimesTo(l);
            std::vector<Time> v;
            std::copy(l.begin(),l.end(),std::back_inserter(v));
            return v;
        }
    }
};

%rename(HestonModelHelper) HestonModelHelperPtr;
class HestonModelHelperPtr : public BlackCalibrationHelperPtr {
  public:
	%extend {
		HestonModelHelperPtr(const Period& maturity,
                             const Calendar& calendar,
                             const Real s0,
                             const Real strikePrice,
                             const Handle<Quote>& volatility,
                             const Handle<YieldTermStructure>& riskFreeRate,
                             const Handle<YieldTermStructure>& dividendYield,
                             BlackCalibrationHelper::CalibrationErrorType errorType
								 = BlackCalibrationHelper::RelativePriceError) {
			return new HestonModelHelperPtr(
				new HestonModelHelper(maturity, calendar, s0, strikePrice,
									  volatility, riskFreeRate, dividendYield,
									  errorType)); 
		}
	}
};

// allow use of vectors of helpers
#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( boost::shared_ptr<CalibrationHelperBase> )
#endif
namespace std {
    %template(CalibrationHelperVector)
        vector<boost::shared_ptr<CalibrationHelperBase> >;
}

// the base class for calibrated models
%{
using QuantLib::CalibratedModel;
%}

%ignore CalibratedModel;
class CalibratedModel {
    #if defined(SWIGRUBY)
    %rename("calibrate!") calibrate;
    #elif defined(SWIGCSHARP)
    %rename("parameters") params;
    #endif
  public:
    Array params() const;
    void calibrate(
        const std::vector<boost::shared_ptr<CalibrationHelperBase> >&,
        OptimizationMethod&, const EndCriteria &,
        const Constraint& constraint = Constraint(),
        const std::vector<Real>& weights = std::vector<Real>());
    void setParams(const Array& params);
    Real value(const Array& params,
               const std::vector<boost::shared_ptr<CalibrationHelperBase> >&);
};


%template(CalibratedModel) boost::shared_ptr<CalibratedModel>;
IsObservable(boost::shared_ptr<CalibratedModel>);

%template(CalibratedModelHandle) Handle<CalibratedModel>;
IsObservable(Handle<CalibratedModel>);
%template(RelinkableCalibratedModelHandle)
RelinkableHandle<CalibratedModel>;

#endif
