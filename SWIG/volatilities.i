/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2008 StatPro Italia srl
 Copyright (C) 2011 Lluis Pujol Bajador
 Copyright (C) 2015 Matthias Groncki
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
%include options.i

%define QL_TYPECHECK_VOLATILITYTYPE       8210    %enddef

%{
using QuantLib::VolatilityType;
using QuantLib::ShiftedLognormal;
using QuantLib::Normal;
%}

enum VolatilityType { ShiftedLognormal, Normal};

#if defined(SWIGPYTHON)
%typemap(in) boost::optional<VolatilityType> %{
    if($input == Py_None)
        $1 = boost::none;
    else if (PyInt_Check($input))
        $1 = (VolatilityType) PyInt_AsLong($input);
    else
        $1 = (VolatilityType) PyLong_AsLong($input);
%}
%typecheck (QL_TYPECHECK_VOLATILITYTYPE) boost::optional<VolatilityType> {
if (PyInt_Check($input) || PyLong_Check($input) || Py_None == $input)
    $1 = 1;
else
    $1 = 0;
}
#endif

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

%{
using QuantLib::SmileSection;
%}

%ignore SmileSection;
class SmileSection{
  public:
    SmileSection(const Date& d,
                 const DayCounter& dc = DayCounter(),
                 const Date& referenceDate = Date(),
                 const VolatilityType type = ShiftedLognormal,
                 const Rate shift = 0.0);
    SmileSection(Time exerciseTime,
                 const DayCounter& dc = DayCounter(),
                 const VolatilityType type = ShiftedLognormal,
                 const Rate shift = 0.0);
    SmileSection() {}

    Real variance(Rate strike) const;
    Volatility volatility(Rate strike) const;
    virtual const Date& exerciseDate() const;
    virtual VolatilityType volatilityType() const;
    virtual Rate shift() const;
    virtual const Date& referenceDate() const;
    virtual Time exerciseTime() const;
    virtual const DayCounter& dayCounter();
    virtual Real optionPrice(Rate strike,
                             Option::Type type = Option::Call,
                             Real discount=1.0) const;
    virtual Real digitalOptionPrice(Rate strike,
                                    Option::Type type = Option::Call,
                                    Real discount=1.0,
                                    Real gap=1.0e-5) const;
    virtual Real vega(Rate strike,
                      Real discount=1.0) const;
    virtual Real density(Rate strike,
                         Real discount=1.0,
                         Real gap=1.0E-4) const;
    Volatility volatility(Rate strike, VolatilityType type, Real shift=0.0) const;
};

%template(SmileSection) boost::shared_ptr<SmileSection>;
IsObservable(boost::shared_ptr<SmileSection>);

%{
using QuantLib::FlatSmileSection;
typedef boost::shared_ptr<SmileSection> FlatSmileSectionPtr;
%}

%rename(FlatSmileSection) FlatSmileSectionPtr;
class FlatSmileSectionPtr
    : public boost::shared_ptr<SmileSection> {
  public:
    %extend {
        FlatSmileSectionPtr(const Date& d,
                         Volatility vol,
                         const DayCounter& dc,
                         const Date& referenceDate = Date(),
                         Real atmLevel = Null<Rate>(),
                         VolatilityType type = ShiftedLognormal,
                         Real shift = 0.0) {
            return new FlatSmileSectionPtr(
                new FlatSmileSection(
                    d,
                    vol,
                    dc,
                    referenceDate,
                    atmLevel,
                    type,
                    shift
                )
            );
        }
        FlatSmileSectionPtr(Time exerciseTime,
                         Volatility vol,
                         const DayCounter& dc,
                         Real atmLevel = Null<Rate>(),
                         VolatilityType type = ShiftedLognormal,
                         Real shift = 0.0) {
            return new FlatSmileSectionPtr(
                new FlatSmileSection(
                    exerciseTime,
                    vol,
                    dc,
                    atmLevel,
                    type,
                    shift
                )
            );
        }
    }
};

%{
using QuantLib::InterpolatedSmileSection;
using QuantLib::Actual365Fixed;
%}

%define export_smileinterpolation_curve(Name,Interpolator)

%{
typedef boost::shared_ptr<SmileSection> Name##Ptr;
%}

%rename(Name) Name##Ptr;
class Name##Ptr : public boost::shared_ptr<SmileSection> {
  public:
    %extend {
        Name##Ptr(
               Time expiryTime,
               const std::vector<Rate>& strikes,
               const std::vector<Handle<Quote> >& stdDevHandles,
               const Handle<Quote>& atmLevel,
               const Interpolator& interpolator = Interpolator(),
               const DayCounter& dc = Actual365Fixed(),
               const VolatilityType type = ShiftedLognormal,
               const Real shift = 0.0) {
            return new Name##Ptr(
                new InterpolatedSmileSection<Interpolator>(
                          expiryTime,strikes,stdDevHandles,atmLevel,interpolator,dc,type,shift));
        }
        Name##Ptr(
               Time expiryTime,
               const std::vector<Rate>& strikes,
               const std::vector<Real>& stdDevs,
               Real atmLevel,
               const Interpolator& interpolator = Interpolator(),
               const DayCounter& dc = Actual365Fixed(),
               const VolatilityType type = ShiftedLognormal,
               const Real shift = 0.0) {
            return new Name##Ptr(
                new InterpolatedSmileSection<Interpolator>(
                          expiryTime,strikes,stdDevs,atmLevel,interpolator,dc,type,shift));
        }
        Name##Ptr(
               const Date& d,
               const std::vector<Rate>& strikes,
               const std::vector<Handle<Quote> >& stdDevHandles,
               const Handle<Quote>& atmLevel,
               const DayCounter& dc = Actual365Fixed(),               
               const Interpolator& interpolator = Interpolator(),
               const Date& referenceDate = Date(),
               const VolatilityType type = ShiftedLognormal,
               const Real shift = 0.0) {
            return new Name##Ptr(
                new InterpolatedSmileSection<Interpolator>(
                          d,strikes,stdDevHandles,atmLevel,dc,interpolator,referenceDate,type,shift));
        }
        Name##Ptr(
               const Date& d,
               const std::vector<Rate>& strikes,
               const std::vector<Real>& stdDevs,
               Real atmLevel,
               const DayCounter& dc = Actual365Fixed(),
               const Interpolator& interpolator = Interpolator(),
               const Date& referenceDate = Date(),
               const VolatilityType type = ShiftedLognormal,
               const Real shift = 0.0) {
            return new Name##Ptr(
                new InterpolatedSmileSection<Interpolator>(
                          d,strikes,stdDevs,atmLevel,dc,interpolator,referenceDate,type,shift));
        }
    }
};

%enddef

export_smileinterpolation_curve(LinearInterpolatedSmileSection, Linear);
export_smileinterpolation_curve(CubicInterpolatedSmileSection, Cubic);
export_smileinterpolation_curve(MonotonicCubicInterpolatedSmileSection, MonotonicCubic);
export_smileinterpolation_curve(SplineCubicInterpolatedSmileSection, SplineCubic);

%{
using QuantLib::SabrSmileSection;
typedef boost::shared_ptr<SmileSection> SabrSmileSectionPtr;
%}

%rename(SabrSmileSection) SabrSmileSectionPtr;
class SabrSmileSectionPtr
    : public boost::shared_ptr<SmileSection> {
  public:
    %extend {
        SabrSmileSectionPtr(const Date& d,
                         Rate forward,
                         const std::vector<Real>& sabrParameters,
                         const DayCounter& dc = Actual365Fixed(),
                         Real shift = 0.0) {
            return new SabrSmileSectionPtr(
                new SabrSmileSection(
                    d,
                    forward,
                    sabrParameters,
                    dc,
                    shift
                )
            );
        }
        SabrSmileSectionPtr(Time timeToExpiry,
                         Rate forward,
                         const std::vector<Real>& sabrParameters,
                         Real shift = 0.0) {
            return new SabrSmileSectionPtr(
                new SabrSmileSection(
                    timeToExpiry,
                    forward,
                    sabrParameters,
                    shift
                )
            );
        }
    }
};

%{
using QuantLib::KahaleSmileSection;
typedef boost::shared_ptr<SmileSection> KahaleSmileSectionPtr;
%}

%rename(KahaleSmileSection) KahaleSmileSectionPtr;
class KahaleSmileSectionPtr
    : public boost::shared_ptr<SmileSection> {
  public:
    %extend {
        KahaleSmileSectionPtr(const boost::shared_ptr<SmileSection> source,
                           const Real atm = Null<Real>(),
                           const bool interpolate = false,
                           const bool exponentialExtrapolation = false,
                           const bool deleteArbitragePoints = false,
                           const std::vector<Real> &moneynessGrid =
                               std::vector<Real>(),
                           const Real gap = 1.0E-5,
                           const int forcedLeftIndex = -1,
                           const int forcedRightIndex = QL_MAX_INTEGER) {
            return new KahaleSmileSectionPtr(
                new KahaleSmileSection(
                    source,
                    atm,
                    interpolate,
                    exponentialExtrapolation,
                    deleteArbitragePoints,
                    moneynessGrid,
                    gap,
                    forcedLeftIndex,
                    QL_MAX_INTEGER
                )
            );
        }
    }
};

%{
using QuantLib::ZabrShortMaturityLognormal;
using QuantLib::ZabrShortMaturityNormal;
using QuantLib::ZabrLocalVolatility;
using QuantLib::ZabrFullFd;
using QuantLib::ZabrSmileSection;
using QuantLib::ZabrInterpolatedSmileSection;
using QuantLib::NoArbSabrSmileSection;
using QuantLib::NoArbSabrInterpolatedSmileSection;
using QuantLib::Option;
%}

struct ZabrShortMaturityLognormal {};
struct ZabrShortMaturityNormal {};
struct ZabrLocalVolatility {};
struct ZabrFullFd {};

%define export_zabrsmilesection_curve(Name,Evaluation)

%{
typedef boost::shared_ptr<SmileSection> Name##Ptr;
%}

%rename(Name) Name##Ptr;
class Name##Ptr : public boost::shared_ptr<SmileSection> {
  public:
    %extend {
        Name##Ptr(
               Time timeToExpiry, 
               Rate forward,
               const std::vector<Real> &zabrParameters,
               const std::vector<Real> &moneyness = std::vector<Real>(),
               const Size fdRefinement = 5) {
            return new Name##Ptr(
                new ZabrSmileSection<Evaluation>(
                          timeToExpiry,forward,zabrParameters,moneyness,fdRefinement));
        }
        Name##Ptr(
               const Date &d, 
               Rate forward,
               const std::vector<Real> &zabrParameters,
               const DayCounter &dc = Actual365Fixed(),
               const std::vector<Real> &moneyness = std::vector<Real>(),
               const Size fdRefinement = 5) {
            return new Name##Ptr(
                new ZabrSmileSection<Evaluation>(
                          d,forward,zabrParameters,dc,moneyness,fdRefinement));
        }
    }
};

%enddef

export_zabrsmilesection_curve(ZabrShortMaturityLognormalSmileSection, ZabrShortMaturityLognormal);
export_zabrsmilesection_curve(ZabrShortMaturityNormalSmileSection, ZabrShortMaturityNormal);
export_zabrsmilesection_curve(ZabrLocalVolatilitySmileSection, ZabrLocalVolatility);
export_zabrsmilesection_curve(ZabrFullFdSmileSection, ZabrFullFd);

%define export_zabrinterpolatedsmilesection_curve(Name,Evaluation)

%{
typedef boost::shared_ptr<SmileSection> Name##Ptr;
%}

%rename(Name) Name##Ptr;
class Name##Ptr : public boost::shared_ptr<SmileSection> {
  public:
    %extend {
        Name##Ptr(
               const Date &optionDate, const Handle<Quote> &forward,
               const std::vector<Rate> &strikes, bool hasFloatingStrikes,
               const Handle<Quote> &atmVolatility,
               const std::vector<Handle<Quote> > &volHandles, Real alpha, Real beta,
               Real nu, Real rho, Real gamma, bool isAlphaFixed = false,
               bool isBetaFixed = false, bool isNuFixed = false,
               bool isRhoFixed = false, bool isGammaFixed = false,
               bool vegaWeighted = true,
               const boost::shared_ptr<EndCriteria> &endCriteria =
               boost::shared_ptr<EndCriteria>(),
               const boost::shared_ptr<OptimizationMethod> &method =
               boost::shared_ptr<OptimizationMethod>(),
               const DayCounter &dc = Actual365Fixed()) {
            return new Name##Ptr(
                new ZabrInterpolatedSmileSection<Evaluation>(
                          optionDate, forward, strikes, hasFloatingStrikes, atmVolatility,
                          volHandles, alpha, beta, nu, rho, gamma, isAlphaFixed,
                          isBetaFixed, isNuFixed, isRhoFixed, isGammaFixed, vegaWeighted,
                          endCriteria, method, dc));
        }
        Name##Ptr(
               const Date &optionDate, const Rate &forward,
               const std::vector<Rate> &strikes, bool hasFloatingStrikes,
               const Volatility &atmVolatility, const std::vector<Volatility> &vols,
               Real alpha, Real beta, Real nu, Real rho, Real gamma,
               bool isAlphaFixed = false, bool isBetaFixed = false,
               bool isNuFixed = false, bool isRhoFixed = false,
               bool isGammaFixed = false, bool vegaWeighted = true,
               const boost::shared_ptr<EndCriteria> &endCriteria =
               boost::shared_ptr<EndCriteria>(),
               const boost::shared_ptr<OptimizationMethod> &method =
               boost::shared_ptr<OptimizationMethod>(),
               const DayCounter &dc = Actual365Fixed()) {
            return new Name##Ptr(
                new ZabrInterpolatedSmileSection<Evaluation>(
                          optionDate, forward, strikes, hasFloatingStrikes, atmVolatility,
                          vols, alpha, beta, nu, rho, gamma, isAlphaFixed,
                          isBetaFixed, isNuFixed, isRhoFixed, isGammaFixed, vegaWeighted,
                          endCriteria, method, dc));
        }
        Real alpha() const {
            return boost::dynamic_pointer_cast<ZabrInterpolatedSmileSection<Evaluation> >(*self)
                ->alpha();
        }
        Real beta() const {
            return boost::dynamic_pointer_cast<ZabrInterpolatedSmileSection<Evaluation> >(*self)
                ->beta();
        }
        Real nu() const {
            return boost::dynamic_pointer_cast<ZabrInterpolatedSmileSection<Evaluation> >(*self)
                ->nu();
        }
        Real rho() const {
            return boost::dynamic_pointer_cast<ZabrInterpolatedSmileSection<Evaluation> >(*self)
                ->rho();
        }
        Real rmsError() const {
            return boost::dynamic_pointer_cast<ZabrInterpolatedSmileSection<Evaluation> >(*self)
                ->rmsError();
        }
        Real maxError() const {
            return boost::dynamic_pointer_cast<ZabrInterpolatedSmileSection<Evaluation> >(*self)
                ->maxError();
        }
        EndCriteria::Type endCriteria() const {
            return boost::dynamic_pointer_cast<ZabrInterpolatedSmileSection<Evaluation> >(*self)
                ->endCriteria();
        }
    }
};

%enddef

export_zabrinterpolatedsmilesection_curve(ZabrShortMaturityLognormalInterpolatedSmileSection, ZabrShortMaturityLognormal);
export_zabrinterpolatedsmilesection_curve(ZabrShortMaturityNormalInterpolatedSmileSection, ZabrShortMaturityNormal);
export_zabrinterpolatedsmilesection_curve(ZabrLocalVolatilityInterpolatedSmileSection, ZabrLocalVolatility);
export_zabrinterpolatedsmilesection_curve(ZabrFullFdInterpolatedSmileSection, ZabrFullFd);

%{
typedef boost::shared_ptr<SmileSection> NoArbSabrSmileSectionPtr;
typedef boost::shared_ptr<SmileSection> NoArbSabrInterpolatedSmileSectionPtr;
%}

%rename(NoArbSabrSmileSection) NoArbSabrSmileSectionPtr;
class NoArbSabrSmileSectionPtr : public boost::shared_ptr<SmileSection> {
  public:
    %extend {
        NoArbSabrSmileSectionPtr(
               Time timeToExpiry, 
               Rate forward,
               const std::vector<Real> &sabrParameters,
               const Real shift = 0.0) {
            return new NoArbSabrSmileSectionPtr(
                new NoArbSabrSmileSection(
                          timeToExpiry,forward,sabrParameters,shift));
        }
        NoArbSabrSmileSectionPtr(
               const Date &d, 
               Rate forward,
               const std::vector<Real> &sabrParameters,
               const DayCounter &dc = Actual365Fixed(),
               const Real shift = 0.0) {
            return new NoArbSabrSmileSectionPtr(
                new NoArbSabrSmileSection(
                          d,forward,sabrParameters,dc,shift));
        }       
    }
};

%rename(NoArbSabrInterpolatedSmileSection) NoArbSabrInterpolatedSmileSectionPtr;
class NoArbSabrInterpolatedSmileSectionPtr : public boost::shared_ptr<SmileSection> {
  public:
    %extend {
        NoArbSabrInterpolatedSmileSectionPtr(
               const Date &optionDate, const Handle<Quote> &forward,
               const std::vector<Rate> &strikes, bool hasFloatingStrikes,
               const Handle<Quote> &atmVolatility,
               const std::vector<Handle<Quote> > &volHandles, Real alpha, Real beta,
               Real nu, Real rho, bool isAlphaFixed = false,
               bool isBetaFixed = false, bool isNuFixed = false,
               bool isRhoFixed = false,
               bool vegaWeighted = true,
               const boost::shared_ptr<EndCriteria> &endCriteria =
               boost::shared_ptr<EndCriteria>(),
               const boost::shared_ptr<OptimizationMethod> &method =
               boost::shared_ptr<OptimizationMethod>(),
               const DayCounter &dc = Actual365Fixed()) {
            return new NoArbSabrInterpolatedSmileSectionPtr(
                new NoArbSabrInterpolatedSmileSection(
                          optionDate, forward, strikes, hasFloatingStrikes, atmVolatility,
                          volHandles, alpha, beta, nu, rho, isAlphaFixed,
                          isBetaFixed, isNuFixed, isRhoFixed, vegaWeighted,
                          endCriteria, method, dc));
        }
        NoArbSabrInterpolatedSmileSectionPtr(
               const Date &optionDate, const Rate &forward,
               const std::vector<Rate> &strikes, bool hasFloatingStrikes,
               const Volatility &atmVolatility, const std::vector<Volatility> &vols,
               Real alpha, Real beta, Real nu, Real rho,
               bool isAlphaFixed = false, bool isBetaFixed = false,
               bool isNuFixed = false, bool isRhoFixed = false,
               bool vegaWeighted = true,
               const boost::shared_ptr<EndCriteria> &endCriteria =
               boost::shared_ptr<EndCriteria>(),
               const boost::shared_ptr<OptimizationMethod> &method =
               boost::shared_ptr<OptimizationMethod>(),
               const DayCounter &dc = Actual365Fixed()) {
            return new NoArbSabrInterpolatedSmileSectionPtr(
                new NoArbSabrInterpolatedSmileSection(
                          optionDate, forward, strikes, hasFloatingStrikes, atmVolatility,
                          vols, alpha, beta, nu, rho, isAlphaFixed,
                          isBetaFixed, isNuFixed, isRhoFixed, vegaWeighted,
                          endCriteria, method, dc));
        }
        Real alpha() const {
            return boost::dynamic_pointer_cast<NoArbSabrInterpolatedSmileSection>(*self)
                ->alpha();
        }
        Real beta() const {
            return boost::dynamic_pointer_cast<NoArbSabrInterpolatedSmileSection>(*self)
                ->beta();
        }
        Real nu() const {
            return boost::dynamic_pointer_cast<NoArbSabrInterpolatedSmileSection>(*self)
                ->nu();
        }
        Real rho() const {
            return boost::dynamic_pointer_cast<NoArbSabrInterpolatedSmileSection>(*self)
                ->rho();
        }
        Real rmsError() const {
            return boost::dynamic_pointer_cast<NoArbSabrInterpolatedSmileSection>(*self)
                ->rmsError();
        }
        Real maxError() const {
            return boost::dynamic_pointer_cast<NoArbSabrInterpolatedSmileSection>(*self)
                ->maxError();
        }
        EndCriteria::Type endCriteria() const {
            return boost::dynamic_pointer_cast<NoArbSabrInterpolatedSmileSection>(*self)
                ->endCriteria();
        }
    }
};

#endif
