/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2008 StatPro Italia srl
 Copyright (C) 2011 Lluis Pujol Bajador
 Copyright (C) 2015 Matthias Groncki
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

#ifndef quantlib_volatilities_i
#define quantlib_volatilities_i

%include common.i
%include date.i
%include daycounters.i
%include types.i
%include currencies.i
%include observer.i
%include marketelements.i
%include interpolation.i
%include indexes.i
%include optimizers.i

%{
using QuantLib::VolatilityType;
using QuantLib::ShiftedLognormal;
using QuantLib::Normal;
%}

enum VolatilityType { ShiftedLognormal, Normal};

%{
using QuantLib::BlackVolTermStructure;
using QuantLib::LocalVolTermStructure;
using QuantLib::OptionletVolatilityStructure;
using QuantLib::SwaptionVolatilityStructure;
%}

%ignore BlackVolTermStructure;
class BlackVolTermStructure : public Extrapolator {
  public:
    Date referenceDate() const;
    DayCounter dayCounter() const;
    Calendar calendar() const;
    Date maxDate() const;
    Time maxTime() const;
    Real minStrike() const;
    Real maxStrike() const;
    Volatility blackVol(const Date&, Real strike,
                        bool extrapolate = false) const;
    Volatility blackVol(Time, Real strike,
                        bool extrapolate = false) const;
    Real blackVariance(const Date&, Real strike,
                       bool extrapolate = false) const;
    Real blackVariance(Time, Real strike,
                       bool extrapolate = false) const;
    Volatility blackForwardVol(const Date&, const Date&,
                               Real strike, bool extrapolate = false) const;
    Volatility blackForwardVol(Time, Time, Real strike,
                               bool extrapolate = false) const;
    Real blackForwardVariance(const Date&, const Date&,
                              Real strike, bool extrapolate = false) const;
    Real blackForwardVariance(Time, Time, Real strike,
                              bool extrapolate = false) const;
};

%template(BlackVolTermStructure) boost::shared_ptr<BlackVolTermStructure>;
IsObservable(boost::shared_ptr<BlackVolTermStructure>);

%template(BlackVolTermStructureHandle) Handle<BlackVolTermStructure>;
IsObservable(Handle<BlackVolTermStructure>);
%template(RelinkableBlackVolTermStructureHandle)
RelinkableHandle<BlackVolTermStructure>;


%ignore LocalVolTermStructure;
class LocalVolTermStructure : public Extrapolator {
  public:
    Date referenceDate() const;
    DayCounter dayCounter() const;
    Calendar calendar() const;
    Date maxDate() const;
    Time maxTime() const;
    Real minStrike() const;
    Real maxStrike() const;
    Volatility localVol(const Date&, Real u,
                        bool extrapolate = false) const;
    Volatility localVol(Time, Real u,
                        bool extrapolate = false) const;
};

%template(LocalVolTermStructure) boost::shared_ptr<LocalVolTermStructure>;
IsObservable(boost::shared_ptr<LocalVolTermStructure>);

%template(LocalVolTermStructureHandle) Handle<LocalVolTermStructure>;
IsObservable(Handle<LocalVolTermStructure>);
%template(RelinkableLocalVolTermStructureHandle)
RelinkableHandle<LocalVolTermStructure>;


%ignore OptionletVolatilityStructure;
class OptionletVolatilityStructure : public Extrapolator {
  public:
    Date referenceDate() const;
    DayCounter dayCounter() const;
    Calendar calendar() const;
    Date maxDate() const;
    Time maxTime() const;
    Real minStrike() const;
    Real maxStrike() const;
    Volatility volatility(const Date&, Real strike,
                          bool extrapolate = false) const;
    Volatility volatility(Time, Real strike,
                          bool extrapolate = false) const;
    Real blackVariance(const Date&, Rate strike,
                       bool extrapolate = false) const ;
    Real blackVariance(Time, Rate strike,
                       bool extrapolate = false) const;
};

%template(OptionletVolatilityStructure)
boost::shared_ptr<OptionletVolatilityStructure>;
IsObservable(boost::shared_ptr<OptionletVolatilityStructure>);

%template(OptionletVolatilityStructureHandle)
Handle<OptionletVolatilityStructure>;
IsObservable(Handle<OptionletVolatilityStructure>);

%template(RelinkableOptionletVolatilityStructureHandle)
RelinkableHandle<OptionletVolatilityStructure>;


%{
using QuantLib::SwaptionVolatilityStructure;
%}

%ignore SwaptionVolatilityStructure;
class SwaptionVolatilityStructure : public Extrapolator {
  public:
    Date referenceDate() const;
    DayCounter dayCounter() const;
    Calendar calendar() const;
    Period maxSwapTenor() const;
    Time maxSwapLength() const;
    Real minStrike() const;
    Real maxStrike() const;
    Volatility volatility(const Date& start, const Period& length,
                          Rate strike, bool extrapolate = false) const;
    Volatility volatility(Time start, Time length,
                          Rate strike, bool extrapolate = false) const;
    Real blackVariance(const Date& start, const Period& length,
                       Rate strike, bool extrapolate = false) const;
    Real blackVariance(Time start, Time length,
                       Rate strike, bool extrapolate = false) const;
    Date optionDateFromTenor(const Period& p) const;
};

%template(SwaptionVolatilityStructure)
    boost::shared_ptr<SwaptionVolatilityStructure>;
IsObservable(boost::shared_ptr<SwaptionVolatilityStructure>);

%template(SwaptionVolatilityStructureHandle)
Handle<SwaptionVolatilityStructure>;
IsObservable(Handle<SwaptionVolatilityStructure>);
%template(RelinkableSwaptionVolatilityStructureHandle)
RelinkableHandle<SwaptionVolatilityStructure>;



// actual term structures below

// constant Black vol term structure
%{
using QuantLib::BlackConstantVol;
typedef boost::shared_ptr<BlackVolTermStructure> BlackConstantVolPtr;
%}

%rename(BlackConstantVol) BlackConstantVolPtr;
class BlackConstantVolPtr : public boost::shared_ptr<BlackVolTermStructure> {
  public:
    %extend {
        BlackConstantVolPtr(const Date& referenceDate,
                            const Calendar & c,
                            Volatility volatility,
                            const DayCounter& dayCounter) {
            return new BlackConstantVolPtr(
                new BlackConstantVol(referenceDate, c,
                                     volatility, dayCounter));
        }
        BlackConstantVolPtr(const Date& referenceDate,
                            const Calendar &c,
                            const Handle<Quote>& volatility,
                            const DayCounter& dayCounter) {
            return new BlackConstantVolPtr(
                new BlackConstantVol(referenceDate, c,
                                     volatility, dayCounter));
        }
        BlackConstantVolPtr(Natural settlementDays, const Calendar& calendar,
                            Volatility volatility,
                            const DayCounter& dayCounter) {
            return new BlackConstantVolPtr(
                new BlackConstantVol(settlementDays, calendar,
                                     volatility, dayCounter));
        }
        BlackConstantVolPtr(Natural settlementDays, const Calendar& calendar,
                            const Handle<Quote>& volatility,
                            const DayCounter& dayCounter) {
            return new BlackConstantVolPtr(
                new BlackConstantVol(settlementDays, calendar,
                                     volatility, dayCounter));
        }
    }
};

// Black ATM curve

%{
using QuantLib::BlackVarianceCurve;
typedef boost::shared_ptr<BlackVolTermStructure> BlackVarianceCurvePtr;
%}

%rename(BlackVarianceCurve) BlackVarianceCurvePtr;
class BlackVarianceCurvePtr : public boost::shared_ptr<BlackVolTermStructure> {
  public:
    %extend {
        BlackVarianceCurvePtr(const Date& referenceDate,
                              const std::vector<Date>& dates,
                              const std::vector<Real>& volatilities,
                              const DayCounter& dayCounter,
                              bool forceMonotoneVariance = true) {
            return new BlackVarianceCurvePtr(
                new BlackVarianceCurve(referenceDate,
                                       dates, volatilities,
                                       dayCounter, forceMonotoneVariance));
        }
    }
};



// Black smiled surface
%{
using QuantLib::BlackVarianceSurface;
typedef boost::shared_ptr<BlackVolTermStructure> BlackVarianceSurfacePtr;
%}

#if defined(SWIGJAVA) || defined(SWIGCSHARP)
%rename(_BlackVarianceSurface) BlackVarianceSurface;
#else
%ignore BlackVarianceSurface;
#endif
class BlackVarianceSurface {
  public:
    enum Extrapolation { ConstantExtrapolation,
                         InterpolatorDefaultExtrapolation };
#if defined(SWIGJAVA) || defined(SWIGCSHARP)
  private:
    BlackVarianceSurface();
#endif
};

%rename(BlackVarianceSurface) BlackVarianceSurfacePtr;
class BlackVarianceSurfacePtr
    : public boost::shared_ptr<BlackVolTermStructure> {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") BlackVarianceSurfacePtr;
    #endif
  public:
    %extend {
        BlackVarianceSurfacePtr(
                const Date& referenceDate,
                const Calendar & cal,
                const std::vector<Date>& dates,
                const std::vector<Real>& strikes,
                const Matrix& blackVols,
                const DayCounter& dayCounter,
                BlackVarianceSurface::Extrapolation lower =
                    BlackVarianceSurface::InterpolatorDefaultExtrapolation,
                BlackVarianceSurface::Extrapolation upper =
                    BlackVarianceSurface::InterpolatorDefaultExtrapolation,
                const std::string& interpolator = "") {
            BlackVarianceSurface* surf =
                new BlackVarianceSurface(referenceDate,cal,
                                         dates,strikes,
                                         blackVols,dayCounter,lower,upper);
            std::string s = boost::algorithm::to_lower_copy(interpolator);
            if (s == "" || s == "bilinear") {
                surf->setInterpolation<QuantLib::Bilinear>();
            } else if (s == "bicubic") {
                surf->setInterpolation<QuantLib::Bicubic>();
            } else {
                QL_FAIL("Unknown interpolator: " << interpolator);
            }
            return new BlackVarianceSurfacePtr(surf);
        }
        void setInterpolation(const std::string& interpolator = "")
        {
            std::string s = boost::algorithm::to_lower_copy(interpolator);
            boost::shared_ptr<BlackVarianceSurface> surf =
                boost::dynamic_pointer_cast<BlackVarianceSurface>(*self);
            if (s == "" || s == "bilinear") {
                surf->setInterpolation<QuantLib::Bilinear>();
            } else if (s == "bicubic") {
                surf->setInterpolation<QuantLib::Bicubic>();
            } else {
                QL_FAIL("Unknown interpolator: " << interpolator);
            }
            
        }
        static const BlackVarianceSurface::Extrapolation
            ConstantExtrapolation =
            BlackVarianceSurface::ConstantExtrapolation;
        static const BlackVarianceSurface::Extrapolation
            InterpolatorDefaultExtrapolation =
            BlackVarianceSurface::InterpolatorDefaultExtrapolation;
    }
};



// constant local vol term structure
%{
using QuantLib::LocalConstantVol;
typedef boost::shared_ptr<LocalVolTermStructure> LocalConstantVolPtr;
%}

%rename(LocalConstantVol) LocalConstantVolPtr;
class LocalConstantVolPtr : public boost::shared_ptr<LocalVolTermStructure> {
  public:
    %extend {
        LocalConstantVolPtr(
                const Date& referenceDate, Volatility volatility,
                const DayCounter& dayCounter) {
            return new LocalConstantVolPtr(
                new LocalConstantVol(referenceDate, volatility, dayCounter));
        }
        LocalConstantVolPtr(
                const Date& referenceDate,
                const Handle<Quote>& volatility,
                const DayCounter& dayCounter) {
            return new LocalConstantVolPtr(
                new LocalConstantVol(referenceDate, volatility, dayCounter));
        }
        LocalConstantVolPtr(
                Integer settlementDays, const Calendar& calendar,
                Volatility volatility,
                const DayCounter& dayCounter) {
            return new LocalConstantVolPtr(
                new LocalConstantVol(settlementDays, calendar,
                                     volatility, dayCounter));
        }
        LocalConstantVolPtr(
                Integer settlementDays, const Calendar& calendar,
                const Handle<Quote>& volatility,
                const DayCounter& dayCounter) {
            return new LocalConstantVolPtr(
                new LocalConstantVol(settlementDays, calendar,
                                     volatility, dayCounter));
        }
    }
};



// constant local vol term structure
%{
using QuantLib::LocalVolSurface;
typedef boost::shared_ptr<LocalVolTermStructure> LocalVolSurfacePtr;
%}
%rename(LocalVolSurface) LocalVolSurfacePtr;
class LocalVolSurfacePtr : public boost::shared_ptr<LocalVolTermStructure> {
    public:
     %extend {
        LocalVolSurfacePtr(const Handle<BlackVolTermStructure>& blackTS,
                        const Handle<YieldTermStructure>& riskFreeTS,
                        const Handle<YieldTermStructure>& dividendTS,
                        const Handle<Quote>& underlying) {
            return new LocalVolSurfacePtr(
                new LocalVolSurface(blackTS, riskFreeTS, 
                    dividendTS, underlying));
        
        }
        LocalVolSurfacePtr(const Handle<BlackVolTermStructure>& blackTS,
                        const Handle<YieldTermStructure>& riskFreeTS,
                        const Handle<YieldTermStructure>& dividendTS,
                        Real underlying) {
            return new LocalVolSurfacePtr(
                new LocalVolSurface(blackTS, riskFreeTS, 
                    dividendTS, underlying));
        }
    }
};


// constant caplet constant term structure
%{
using QuantLib::ConstantOptionletVolatility;
typedef boost::shared_ptr<OptionletVolatilityStructure>
    ConstantOptionletVolatilityPtr;
%}

%rename(ConstantOptionletVolatility) ConstantOptionletVolatilityPtr;
class ConstantOptionletVolatilityPtr
    : public boost::shared_ptr<OptionletVolatilityStructure> {
  public:
    %extend {
        ConstantOptionletVolatilityPtr(const Date& referenceDate,
                                       const Calendar &cal,
                                       BusinessDayConvention bdc,
                                       Volatility volatility,
                                       const DayCounter& dayCounter,
                                       const VolatilityType type = ShiftedLognormal,
                                       const Real shift = 0.0) {
            return new ConstantOptionletVolatilityPtr(
                new ConstantOptionletVolatility(referenceDate,
                                                cal, bdc, volatility,
                                                dayCounter, type, shift));
        }
        ConstantOptionletVolatilityPtr(const Date& referenceDate,
                                       const Calendar &cal,
                                       BusinessDayConvention bdc,
                                       const Handle<Quote>& volatility,
                                       const DayCounter& dayCounter,
                                       const VolatilityType type = ShiftedLognormal,
                                       const Real shift = 0.0) {
            return new ConstantOptionletVolatilityPtr(
                new ConstantOptionletVolatility(referenceDate,
                                                cal, bdc, volatility,
                                                dayCounter, type, shift));
        }
        ConstantOptionletVolatilityPtr(Natural settlementDays,
                                       const Calendar &cal,
                                       BusinessDayConvention bdc,
                                       Volatility volatility,
                                       const DayCounter& dayCounter,
                                       const VolatilityType type = ShiftedLognormal,
                                       const Real shift = 0.0) {
            return new ConstantOptionletVolatilityPtr(
                new ConstantOptionletVolatility(settlementDays,
                                                cal, bdc, volatility,
                                                dayCounter, type, shift));
        }
        ConstantOptionletVolatilityPtr(Natural settlementDays,
                                       const Calendar &cal,
                                       BusinessDayConvention bdc,
                                       const Handle<Quote>& volatility,
                                       const DayCounter& dayCounter,
                                       const VolatilityType type = ShiftedLognormal,
                                       const Real shift = 0.0) {
            return new ConstantOptionletVolatilityPtr(
                new ConstantOptionletVolatility(settlementDays,
                                                cal, bdc, volatility,
                                                dayCounter, type, shift));
        }
    }
};



%{
using QuantLib::ConstantSwaptionVolatility;
typedef boost::shared_ptr<SwaptionVolatilityStructure>
    ConstantSwaptionVolatilityPtr;
%}

%rename(ConstantSwaptionVolatility) ConstantSwaptionVolatilityPtr;
class ConstantSwaptionVolatilityPtr
    : public boost::shared_ptr<SwaptionVolatilityStructure> {
  public:
    %extend {
        ConstantSwaptionVolatilityPtr(Natural settlementDays,
                                      const Calendar& cal,
                                      BusinessDayConvention bdc,
                                      const Handle<Quote>& volatility,
                                      const DayCounter& dc,
                                      const VolatilityType type = ShiftedLognormal,
                                      const Real shift = 0.0) {
            return new ConstantSwaptionVolatilityPtr(
                new ConstantSwaptionVolatility(settlementDays, cal, bdc,
                                               volatility, dc, type, shift));
        }
        ConstantSwaptionVolatilityPtr(const Date& referenceDate,
                                      const Calendar& cal,
                                      BusinessDayConvention bdc,
                                      const Handle<Quote>& volatility,
                                      const DayCounter& dc,
                                      const VolatilityType type = ShiftedLognormal,
                                      const Real shift = 0.0) {
            return new ConstantSwaptionVolatilityPtr(
                new ConstantSwaptionVolatility(referenceDate, cal, bdc,
                                               volatility, dc, type, shift));
        }
        ConstantSwaptionVolatilityPtr(Natural settlementDays,
                                      const Calendar& cal,
                                      BusinessDayConvention bdc,
                                      Volatility volatility,
                                      const DayCounter& dc,
                                      const VolatilityType type = ShiftedLognormal,
                                      const Real shift = 0.0) {
            return new ConstantSwaptionVolatilityPtr(
                new ConstantSwaptionVolatility(settlementDays, cal, bdc,
                                               volatility, dc, type, shift));
        }
        ConstantSwaptionVolatilityPtr(const Date& referenceDate,
                                      const Calendar& cal,
                                      BusinessDayConvention bdc,
                                      Volatility volatility,
                                      const DayCounter& dc,
                                      const VolatilityType type = ShiftedLognormal,
                                      const Real shift = 0.0) {
            return new ConstantSwaptionVolatilityPtr(
                new ConstantSwaptionVolatility(referenceDate, cal, bdc,
                                               volatility, dc, type, shift));
        }
    }
};

%{
using QuantLib::SwaptionVolatilityMatrix;
typedef boost::shared_ptr<SwaptionVolatilityStructure>
    SwaptionVolatilityMatrixPtr;
%}

%rename(SwaptionVolatilityMatrix) SwaptionVolatilityMatrixPtr;
class SwaptionVolatilityMatrixPtr
    : public boost::shared_ptr<SwaptionVolatilityStructure> {
  public:
    %extend {
        SwaptionVolatilityMatrixPtr(const Date& referenceDate,
                                    const std::vector<Date>& dates,
                                    const std::vector<Period>& lengths,
                                    const Matrix& vols,
                                    const DayCounter& dayCounter,
                                    const bool flatExtrapolation = false,
                                    const VolatilityType type = ShiftedLognormal,
                                    const Matrix& shifts = Matrix()) {
            return new SwaptionVolatilityMatrixPtr(
                new SwaptionVolatilityMatrix(referenceDate,dates,lengths,
                                             vols,dayCounter,
                                             flatExtrapolation, type, shifts));
        }
        SwaptionVolatilityMatrixPtr(
                        const Calendar& calendar,
                        BusinessDayConvention bdc,
                        const std::vector<Period>& optionTenors,
                        const std::vector<Period>& swapTenors,
                        const std::vector<std::vector<Handle<Quote> > >& vols,
                        const DayCounter& dayCounter,
                        const bool flatExtrapolation = false,
                        const VolatilityType type = ShiftedLognormal,
                        const std::vector<std::vector<Real> >& shifts =
                                          std::vector<std::vector<Real> >()) {
            return new SwaptionVolatilityMatrixPtr(
                new SwaptionVolatilityMatrix(calendar,bdc,optionTenors,
                                             swapTenors,vols,dayCounter,
                                             flatExtrapolation, type, shifts));
        }
        SwaptionVolatilityMatrixPtr(const Calendar& calendar,
                                    BusinessDayConvention bdc,
                                    const std::vector<Period>& optionTenors,
                                    const std::vector<Period>& swapTenors,
                                    const Matrix& vols,
                                    const DayCounter& dayCounter,
                                    const bool flatExtrapolation = false,
                                    const VolatilityType type = ShiftedLognormal,
                                    const Matrix& shifts = Matrix()) {
            return new SwaptionVolatilityMatrixPtr(
                new SwaptionVolatilityMatrix(calendar,bdc,optionTenors,
                                             swapTenors,vols,dayCounter,
                                             flatExtrapolation, type, shifts));
        }
    }
};

%{
using QuantLib::SwaptionVolCube1;
using QuantLib::SwaptionVolCube2;
typedef boost::shared_ptr<SwaptionVolatilityStructure> SwaptionVolCube1Ptr;
typedef boost::shared_ptr<SwaptionVolatilityStructure> SwaptionVolCube2Ptr;
%}

%rename(SwaptionVolCube1) SwaptionVolCube1Ptr;
class SwaptionVolCube1Ptr
    : public boost::shared_ptr<SwaptionVolatilityStructure> {
  public:
    %extend {
        SwaptionVolCube1Ptr(
             const Handle<SwaptionVolatilityStructure>& atmVolStructure,
             const std::vector<Period>& optionTenors,
             const std::vector<Period>& swapTenors,
             const std::vector<Spread>& strikeSpreads,
             const std::vector<std::vector<Handle<Quote> > >& volSpreads,
             const SwapIndexPtr& swapIndexBase,
             const SwapIndexPtr& shortSwapIndexBase,
             bool vegaWeightedSmileFit,
             const std::vector<std::vector<Handle<Quote> > >& parametersGuess,
             const std::vector<bool>& isParameterFixed,
             bool isAtmCalibrated,
             const boost::shared_ptr<EndCriteria>& endCriteria
                                           = boost::shared_ptr<EndCriteria>(),
             Real maxErrorTolerance = Null<Real>(),
             const boost::shared_ptr<OptimizationMethod>& optMethod
                                  = boost::shared_ptr<OptimizationMethod>()) {
            const boost::shared_ptr<SwapIndex> swi =
                boost::dynamic_pointer_cast<SwapIndex>(swapIndexBase);
            const boost::shared_ptr<SwapIndex> shortSwi =
                boost::dynamic_pointer_cast<SwapIndex>(shortSwapIndexBase);
            return new SwaptionVolCube1Ptr(
                new SwaptionVolCube1(
                    atmVolStructure,optionTenors,swapTenors, strikeSpreads,
                    volSpreads, swi, shortSwi, vegaWeightedSmileFit,
                    parametersGuess,isParameterFixed,isAtmCalibrated,
                    endCriteria,maxErrorTolerance,optMethod));
        }

        Matrix sparseSabrParameters() const {
            return boost::dynamic_pointer_cast<SwaptionVolCube1>(*self)
                ->sparseSabrParameters();
        }

        Matrix denseSabrParameters() const {
            return boost::dynamic_pointer_cast<SwaptionVolCube1>(*self)
                ->denseSabrParameters();
        }

        Matrix marketVolCube() const {
            return boost::dynamic_pointer_cast<SwaptionVolCube1>(*self)
                ->marketVolCube();
        }

        Matrix volCubeAtmCalibrated() const {
            return boost::dynamic_pointer_cast<SwaptionVolCube1>(*self)
                ->volCubeAtmCalibrated();
        }
    }
};

%rename(SwaptionVolCube2) SwaptionVolCube2Ptr;
class SwaptionVolCube2Ptr
    : public boost::shared_ptr<SwaptionVolatilityStructure> {
  public:
    %extend {
        SwaptionVolCube2Ptr(
                   const Handle<SwaptionVolatilityStructure>& atmVolStructure,
                   const std::vector<Period>& optionTenors,
                   const std::vector<Period>& swapTenors,
                   const std::vector<Spread>& strikeSpreads,
                   const std::vector<std::vector<Handle<Quote> > >& volSpreads,
                   const SwapIndexPtr& swapIndexBase,
                   const SwapIndexPtr& shortSwapIndexBase,
                   bool vegaWeightedSmileFit) {
            const boost::shared_ptr<SwapIndex> swi =
                boost::dynamic_pointer_cast<SwapIndex>(swapIndexBase);
            const boost::shared_ptr<SwapIndex> shortSwi =
                boost::dynamic_pointer_cast<SwapIndex>(shortSwapIndexBase);
            return new SwaptionVolCube2Ptr(
                new SwaptionVolCube2(
                    atmVolStructure,optionTenors,swapTenors,strikeSpreads,
                    volSpreads, swi, shortSwi, vegaWeightedSmileFit));
        }
    }
};

#endif
