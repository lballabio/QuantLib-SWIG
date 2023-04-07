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

%include calibratedmodel.i
%include date.i
%include calendars.i
%include daycounters.i
%include cashflows.i
%include marketelements.i
%include termstructures.i
%include linearalgebra.i
%include types.i
%include vectors.i

%{
using QuantLib::VanillaSwap;
using QuantLib::Swaption;
using QuantLib::BlackCalibrationHelper;
using QuantLib::SwaptionHelper;
using QuantLib::CapHelper;
using QuantLib::HestonModelHelper;
%}

%shared_ptr(BlackCalibrationHelper)
class BlackCalibrationHelper : public CalibrationHelper {
  public:
    enum CalibrationErrorType { RelativePriceError, PriceError, ImpliedVolError };

    void setPricingEngine(const ext::shared_ptr<PricingEngine>& engine);
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
    ext::shared_ptr<BlackCalibrationHelper> as_black_helper(const ext::shared_ptr<CalibrationHelper>& h) {
        return ext::dynamic_pointer_cast<BlackCalibrationHelper>(h);
    }
    ext::shared_ptr<SwaptionHelper> as_swaption_helper(const ext::shared_ptr<BlackCalibrationHelper>& h) {
        return ext::dynamic_pointer_cast<SwaptionHelper>(h);
    }
%}

%shared_ptr(SwaptionHelper)
class SwaptionHelper : public BlackCalibrationHelper {
  public:
    SwaptionHelper(const Period& maturity, const Period& length,
                      const Handle<Quote>& volatility,
                      const ext::shared_ptr<IborIndex>& index,
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
                      const ext::shared_ptr<IborIndex>& index,
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
                      const ext::shared_ptr<IborIndex>& index,
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

    ext::shared_ptr<VanillaSwap> underlyingSwap() const;
    ext::shared_ptr<Swaption> swaption() const;

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
              const ext::shared_ptr<IborIndex>& index,
              Frequency fixedLegFrequency,
              const DayCounter& fixedLegDayCounter,
              bool includeFirstSwaplet,
              const Handle<YieldTermStructure>& termStructure,
              BlackCalibrationHelper::CalibrationErrorType errorType
                                = BlackCalibrationHelper::RelativePriceError,
              const VolatilityType type = ShiftedLognormal,
              const Real shift = 0.0);
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
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<BlackCalibrationHelper> )
#endif
namespace std {
    %template(BlackCalibrationHelperVector)
        vector<ext::shared_ptr<BlackCalibrationHelper> >;
}

#endif
