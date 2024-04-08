/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2008, 2023 StatPro Italia srl
 Copyright (C) 2011 Lluis Pujol Bajador
 Copyright (C) 2015 Matthias Groncki
 Copyright (C) 2016 Peter Caspers
 Copyright (C) 2018, 2019, 2020 Matthias Lungwitz
 Copyright (C) 2022 Skandinaviska Enskilda Banken AB (publ)

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

#ifndef quantlib_local_volatilities_i
#define quantlib_local_volatilities_i

%include volatilities.i


// constant local vol term structure
%{
using QuantLib::LocalConstantVol;
%}

%shared_ptr(LocalConstantVol);
class LocalConstantVol : public LocalVolTermStructure {
  public:
    LocalConstantVol(const Date& referenceDate, Volatility volatility,
                     const DayCounter& dayCounter);
    LocalConstantVol(const Date& referenceDate,
                     const Handle<Quote>& volatility,
                     const DayCounter& dayCounter);
    LocalConstantVol(Integer settlementDays, const Calendar& calendar,
                     Volatility volatility,
                     const DayCounter& dayCounter);
    LocalConstantVol(Integer settlementDays, const Calendar& calendar,
                     const Handle<Quote>& volatility,
                     const DayCounter& dayCounter);
};



// local vol surface
%{
using QuantLib::LocalVolSurface;
%}

%shared_ptr(LocalVolSurface);
class LocalVolSurface : public LocalVolTermStructure {
  public:
    LocalVolSurface(const Handle<BlackVolTermStructure>& blackTS,
                    const Handle<YieldTermStructure>& riskFreeTS,
                    const Handle<YieldTermStructure>& dividendTS,
                    const Handle<Quote>& underlying);
    LocalVolSurface(const Handle<BlackVolTermStructure>& blackTS,
                    const Handle<YieldTermStructure>& riskFreeTS,
                    const Handle<YieldTermStructure>& dividendTS,
                    Real underlying);
};



// no except local vol surface (override bad points - use with care)
%{
using QuantLib::NoExceptLocalVolSurface;
%}

%shared_ptr(NoExceptLocalVolSurface);
class NoExceptLocalVolSurface : public LocalVolSurface {
  public:
    NoExceptLocalVolSurface(const Handle<BlackVolTermStructure>& blackTS,
                            const Handle<YieldTermStructure>& riskFreeTS,
                            const Handle<YieldTermStructure>& dividendTS,
                            const Handle<Quote>& underlying,
                            Real illegalLocalVolOverwrite);
    NoExceptLocalVolSurface(const Handle<BlackVolTermStructure>& blackTS,
                            const Handle<YieldTermStructure>& riskFreeTS,
                            const Handle<YieldTermStructure>& dividendTS,
                            Real underlying,
                            Real illegalLocalVolOverwrite);
};


%{
using QuantLib::FixedLocalVolSurface;
using QuantLib::GridModelLocalVolSurface;
%}

%shared_ptr(FixedLocalVolSurface);
class FixedLocalVolSurface : public LocalVolTermStructure {
    %warnfilter(509) FixedLocalVolSurface;
  public:
    enum Extrapolation { ConstantExtrapolation,
                         InterpolatorDefaultExtrapolation };
    %extend {
        FixedLocalVolSurface(const Date& referenceDate,
                             const std::vector<Date>& dates,
                             const std::vector<Real>& strikes,
                             Matrix localVolMatrix,
                             const DayCounter& dayCounter,
                             Extrapolation lowerExtrapolation = ConstantExtrapolation,
                             Extrapolation upperExtrapolation = ConstantExtrapolation) {
            ext::shared_ptr<Matrix> pMatrix = ext::make_shared<Matrix>(localVolMatrix);
            return new FixedLocalVolSurface(referenceDate, dates, strikes, pMatrix,
                                            dayCounter, lowerExtrapolation, upperExtrapolation);
        }
        FixedLocalVolSurface(const Date& referenceDate,
                             const std::vector<Time>& times,
                             const std::vector<Real>& strikes,
                             Matrix localVolMatrix,
                             const DayCounter& dayCounter,
                             Extrapolation lowerExtrapolation = ConstantExtrapolation,
                             Extrapolation upperExtrapolation = ConstantExtrapolation) {
            ext::shared_ptr<Matrix> pMatrix = ext::make_shared<Matrix>(localVolMatrix);
            return new FixedLocalVolSurface(referenceDate, times, strikes, pMatrix,
                                            dayCounter, lowerExtrapolation, upperExtrapolation);
        }
        FixedLocalVolSurface(const Date& referenceDate,
                             const std::vector<Time>& times,
                             const std::vector<std::vector<Real>>& strikes,
                             Matrix localVolMatrix,
                             const DayCounter& dayCounter,
                             Extrapolation lowerExtrapolation = ConstantExtrapolation,
                             Extrapolation upperExtrapolation = ConstantExtrapolation) {
            std::vector<ext::shared_ptr<std::vector<Real> > > pStrikes(strikes.size());
            for (Size i=0; i<strikes.size(); ++i)
                pStrikes[i] = ext::make_shared<std::vector<Real>>(strikes[i]);
            ext::shared_ptr<Matrix> pMatrix = ext::make_shared<Matrix>(localVolMatrix);
            return new FixedLocalVolSurface(referenceDate, times, pStrikes, pMatrix,
                                            dayCounter, lowerExtrapolation, upperExtrapolation);
        }
                
        void setInterpolation(const std::string& interpolator = "") {
            const std::string s = boost::to_lower_copy(interpolator);
            if (s == "" || s == "linear") {
                self->setInterpolation<QuantLib::Linear>();
            } else if (s == "cubic") {
                self->setInterpolation<QuantLib::Cubic>();
            } else {
                QL_FAIL("Unknown interpolator: " << interpolator);
            }
        }        
    }
};

%shared_ptr(GridModelLocalVolSurface);
class GridModelLocalVolSurface
#if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    : public LocalVolTermStructure,
      public CalibratedModel
#else
    // multiple inheritance not allowed
    : public LocalVolTermStructure
#endif
{
  public:
    typedef FixedLocalVolSurface::Extrapolation Extrapolation;

    %extend {
        GridModelLocalVolSurface(
                    const Date& referenceDate,
                    const std::vector<Date>& dates,
                    const std::vector<std::vector<Real>>& strikes,
                    const DayCounter& dayCounter,
                    Extrapolation lowerExtrapolation = FixedLocalVolSurface::ConstantExtrapolation,
                    Extrapolation upperExtrapolation = FixedLocalVolSurface::ConstantExtrapolation) {
            std::vector<ext::shared_ptr<std::vector<Real> > > pStrikes(strikes.size());
            for (Size i=0; i<strikes.size(); ++i)
                pStrikes[i] = ext::make_shared<std::vector<Real>>(strikes[i]);
            return new GridModelLocalVolSurface(referenceDate, dates, pStrikes,
                                                dayCounter, lowerExtrapolation, upperExtrapolation);
        }
    }

    #if defined(SWIGCSHARP) || defined(SWIGJAVA)
    // not inheriting from CalibratedModel, add minimal set of methods
    #if defined(SWIGCSHARP)
    %rename("parameters") params;
    #endif
    Array params() const;
    virtual void calibrate(
        const std::vector<ext::shared_ptr<CalibrationHelper> >&,
        OptimizationMethod&, const EndCriteria &,
        const Constraint& constraint = Constraint(),
        const std::vector<Real>& weights = std::vector<Real>(),
        const std::vector<bool>& fixParameters = std::vector<bool>());
    EndCriteria::Type endCriteria() const;
    #endif

};


#endif

