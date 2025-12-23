/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2008 StatPro Italia srl
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
 <https://www.quantlib.org/license.shtml>.

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
%include termstructures.i
%include vectors.i
%include tuple.i

%define QL_TYPECHECK_VOLATILITYTYPE       8210    %enddef

%{
using QuantLib::VolatilityType;
using QuantLib::ShiftedLognormal;
using QuantLib::Normal;
%}

enum VolatilityType { ShiftedLognormal, Normal };

#if defined(SWIGPYTHON)
%typemap(in) ext::optional<VolatilityType> %{
    if ($input == Py_None)
        $1 = ext::nullopt;
    else if (PyLong_Check($input))
        $1 = (VolatilityType)PyLong_AsLong($input);
    else
        SWIG_exception(SWIG_TypeError, "int expected");
%}
%typecheck (QL_TYPECHECK_VOLATILITYTYPE) ext::optional<VolatilityType> %{
    $1 = (PyLong_Check($input) || $input == Py_None) ? 1 : 0;
%}
#endif

%{
using QuantLib::VolatilityTermStructure;
using QuantLib::BlackVolTermStructure;
using QuantLib::LocalVolTermStructure;
using QuantLib::OptionletVolatilityStructure;
using QuantLib::SwaptionVolatilityStructure;
using QuantLib::YoYOptionletVolatilitySurface;
%}

%shared_ptr(VolatilityTermStructure);
class VolatilityTermStructure : public TermStructure {
  private:
    VolatilityTermStructure();
  public:
    Real minStrike() const;
    Real maxStrike() const;
};


%shared_ptr(BlackVolTermStructure);
class BlackVolTermStructure : public VolatilityTermStructure {
  private:
    BlackVolTermStructure();
  public:
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

%template(BlackVolTermStructureHandle) Handle<BlackVolTermStructure>;
%template(RelinkableBlackVolTermStructureHandle) RelinkableHandle<BlackVolTermStructure>;


%shared_ptr(LocalVolTermStructure);
class LocalVolTermStructure : public VolatilityTermStructure {
  private:
    LocalVolTermStructure();
  public:
    Volatility localVol(const Date&, Real u,
                        bool extrapolate = false) const;
    Volatility localVol(Time, Real u,
                        bool extrapolate = false) const;
};

%template(LocalVolTermStructureHandle) Handle<LocalVolTermStructure>;
%template(RelinkableLocalVolTermStructureHandle) RelinkableHandle<LocalVolTermStructure>;


%shared_ptr(OptionletVolatilityStructure);
class OptionletVolatilityStructure : public VolatilityTermStructure {
  private:
    OptionletVolatilityStructure();
  public:
    Volatility volatility(const Date&, Real strike,
                          bool extrapolate = false) const;
    Volatility volatility(Time, Real strike,
                          bool extrapolate = false) const;
    Real blackVariance(const Date&, Rate strike,
                       bool extrapolate = false) const ;
    Real blackVariance(Time, Rate strike,
                       bool extrapolate = false) const;
};

%template(OptionletVolatilityStructureHandle) Handle<OptionletVolatilityStructure>;
%template(RelinkableOptionletVolatilityStructureHandle) RelinkableHandle<OptionletVolatilityStructure>;


%shared_ptr(YoYOptionletVolatilitySurface)
class YoYOptionletVolatilitySurface : public VolatilityTermStructure {
  private:
    YoYOptionletVolatilitySurface();
  public:
    Period observationLag() const;
    Real frequency() const;
    bool indexIsInterpolated() const;
    Date baseDate() const;
    Time timeFromBase(const Date& date,
                      const Period& obsLag = Period(-1,Days)) const;
    Real minStrike() const;
    Real maxStrike() const;
    Volatility baseLevel() const;
    Volatility volatility(const Date& maturityDate, Real strike,
                          const Period& obsLag = Period(-1,Days),
                          bool extrapolate = false) const;
    Volatility volatility(const Period& optionTenor, Real strike,
                          const Period& obsLag = Period(-1,Days),
                          bool extrapolate = false) const;
    Real totalVariance(const Date& exerciseDate, Rate strike,
                       const Period& obsLag = Period(-1,Days),
                       bool extrapolate = false) const ;
    Real totalVariance(const Period& optionTenor, Rate strike,
                       const Period& obsLag = Period(-1,Days),
                       bool extrapolate = false) const;
};

%template(YoYOptionletVolatilitySurfaceHandle) Handle<YoYOptionletVolatilitySurface>;
%template(RelinkableYoYOptionletVolatilitySurfaceHandle) RelinkableHandle<YoYOptionletVolatilitySurface>;


%{
using QuantLib::SmileSection;
%}

%shared_ptr(SmileSection);
class SmileSection : public Observable {
  private:
    SmileSection();
  public:
    Real minStrike() const;
    Real maxStrike() const;
    Real atmLevel() const;
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

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<SmileSection> )
#endif
%template(SmileSectionVector) std::vector<ext::shared_ptr<SmileSection> >;


%{
using QuantLib::SwaptionVolatilityStructure;
%}

%shared_ptr(SwaptionVolatilityStructure);
class SwaptionVolatilityStructure : public VolatilityTermStructure {
  private:
    SwaptionVolatilityStructure();
  public:
    Volatility volatility(const Date& start, const Period& length,
                          Rate strike, bool extrapolate = false) const;
    Volatility volatility(Time start, Time length,
                          Rate strike, bool extrapolate = false) const;
    Real blackVariance(const Date& start, const Period& length,
                       Rate strike, bool extrapolate = false) const;
    Real blackVariance(Time start, Time length,
                       Rate strike, bool extrapolate = false) const;
    Date optionDateFromTenor(const Period& p) const;
    Real shift(const Period& optionTenor,
               const Period& swapTenor,
               bool extrapolate = false) const;
    Real shift(const Date& optionDate,
               const Period& swapTenor,
               bool extrapolate = false) const;
    Real shift(Time optionTime,
               const Period& swapTenor,
               bool extrapolate = false) const;
    Real shift(const Period& optionTenor,
               Time swapLength,
               bool extrapolate = false) const;
    Real shift(const Date& optionDate,
               Time swapLength,
               bool extrapolate = false) const;
    Real shift(Time optionTime,
               Time swapLength,
               bool extrapolate = false) const;
    ext::shared_ptr<SmileSection> smileSection(const Period& optionTenor,
                                                 const Period& swapTenor,
                                                 bool extr = false) const;
    ext::shared_ptr<SmileSection> smileSection(const Date& optionDate,
                                                 const Period& swapTenor,
                                                 bool extr = false) const;
    ext::shared_ptr<SmileSection> smileSection(Time optionTime,
                                                 const Period& swapTenor,
                                                 bool extr = false) const;
    ext::shared_ptr<SmileSection> smileSection(const Period& optionTenor,
                                                 Time swapLength,
                                                 bool extr = false) const;
    ext::shared_ptr<SmileSection> smileSection(const Date& optionDate,
                                                 Time swapLength,
                                                 bool extr = false) const;
    ext::shared_ptr<SmileSection> smileSection(Time optionTime,
                                                 Time swapLength,
                                                 bool extr = false) const;
};

%template(SwaptionVolatilityStructureHandle) Handle<SwaptionVolatilityStructure>;
%template(RelinkableSwaptionVolatilityStructureHandle) RelinkableHandle<SwaptionVolatilityStructure>;



%{
using QuantLib::sabrVolatility;
using QuantLib::shiftedSabrVolatility;
using QuantLib::sabrFlochKennedyVolatility;
%}

Real sabrVolatility(Rate strike,
                    Rate forward,
                    Time expiryTime,
                    Real alpha,
                    Real beta,
                    Real nu,
                    Real rho,
                    VolatilityType volatilityType = VolatilityType::ShiftedLognormal);

Real shiftedSabrVolatility(Rate strike,
                           Rate forward,
                           Time expiryTime,
                           Real alpha,
                           Real beta,
                           Real nu,
                           Real rho,
                           Real shift,
                           VolatilityType volatilityType = VolatilityType::ShiftedLognormal);

Real sabrFlochKennedyVolatility(Rate strike,
                                Rate forward,
                                Time expiryTime,
                                Real alpha,
                                Real beta,
                                Real nu,
                                Real rho);

%rename(sabrGuess) _sabrGuess;
%inline %{
    std::vector<Real> _sabrGuess(Real k_m, Volatility vol_m,
                                 Real k_0, Volatility vol_0,
                                 Real k_p, Volatility vol_p,
                                 Rate forward,
                                 Time expiryTime,
                                 Real beta,
                                 Real shift,
                                 VolatilityType volatilityType) {
        auto [alpha, beta_1, nu, rho] =
            QuantLib::sabrGuess(k_m, vol_m, k_0, vol_0, k_p, vol_p,
                                forward, expiryTime, beta, shift, volatilityType);
        return { alpha, beta_1, nu, rho };
    }
%}


%{
// safe version which copies its arguments
class SafeSABRInterpolation {
  public:
    SafeSABRInterpolation(const Array& strikes,
                          const Array& volatilities,
                          Time expiryTime,
                          Real forward,
                          Real alpha,
                          Real beta,
                          Real nu,
                          Real rho,
                          bool alphaIsFixed = false,
                          bool betaIsFixed = false,
                          bool nuIsFixed = false,
                          bool rhoIsFixed = false,
                          bool vegaWeighted = true,
                          const ext::shared_ptr<EndCriteria>& endCriteria = {},
                          const ext::shared_ptr<OptimizationMethod>& optMethod = {},
                          Real errorAccept=0.0020,
                          bool useMaxError=false,
                          Size maxGuesses=50,
                          Real shift = 0.0,
                          VolatilityType volatilityType = VolatilityType::ShiftedLognormal)
    : x_(strikes), y_(volatilities), forward_(forward),
      f_(x_.begin(),x_.end(),y_.begin(),
         expiryTime, forward_, alpha, beta, nu, rho,
         alphaIsFixed, betaIsFixed,
         nuIsFixed, rhoIsFixed,
         vegaWeighted, endCriteria, optMethod,
         errorAccept, useMaxError, maxGuesses,
         shift, volatilityType) {f_.update();}
    Real operator()(Real x, bool allowExtrapolation=false) const {
        return f_(x, allowExtrapolation);
    }
    Real alpha() const {return f_.alpha();}
    Real beta() const {return f_.beta();}
    Real rho() const {return f_.rho();}
    Real nu() const {return f_.nu();}

  private:
    Array x_, y_;  // passed via iterators, need to stay alive
    Real forward_; // passed by reference, same
    QuantLib::SABRInterpolation f_;
};
%}

%rename(SABRInterpolation) SafeSABRInterpolation;
class SafeSABRInterpolation {
    #if defined(SWIGCSHARP)
    %rename(call) operator();
    #endif
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") SafeSABRInterpolation;
    #endif
  public:
    SafeSABRInterpolation(const Array& strikes,
                          const Array& volatilities,
                          Time expiryTime,
                          Real forward,
                          Real alpha,
                          Real beta,
                          Real nu,
                          Real rho,
                          bool alphaIsFixed = false,
                          bool betaIsFixed = false,
                          bool nuIsFixed = false,
                          bool rhoIsFixed = false,
                          bool vegaWeighted = true,
                          const ext::shared_ptr<EndCriteria>& endCriteria = {},
                          const ext::shared_ptr<OptimizationMethod>& optMethod = {},
                          const Real errorAccept=0.0020,
                          const bool useMaxError=false,
                          const Size maxGuesses=50,
                          const Real shift = 0.0,
                          VolatilityType volatilityType = VolatilityType::ShiftedLognormal);
    Real operator()(Real x, bool allowExtrapolation=false) const;
    Real alpha() const;
    Real beta() const;
    Real rho() const;
    Real nu() const;
};


// actual term structures below

// constant Black vol term structure
%{
using QuantLib::BlackConstantVol;
%}

%shared_ptr(BlackConstantVol);
class BlackConstantVol : public BlackVolTermStructure {
  public:
    BlackConstantVol(const Date& referenceDate,
                     const Calendar & c,
                     Volatility volatility,
                     const DayCounter& dayCounter);
    BlackConstantVol(const Date& referenceDate,
                     const Calendar &c,
                     const Handle<Quote>& volatility,
                     const DayCounter& dayCounter);
    BlackConstantVol(Natural settlementDays, const Calendar& calendar,
                     Volatility volatility,
                     const DayCounter& dayCounter);
    BlackConstantVol(Natural settlementDays, const Calendar& calendar,
                     const Handle<Quote>& volatility,
                     const DayCounter& dayCounter);
};

// Black ATM curve

%{
using QuantLib::BlackVarianceCurve;
%}

%shared_ptr(BlackVarianceCurve);
class BlackVarianceCurve : public BlackVolTermStructure {
  public:
    BlackVarianceCurve(const Date& referenceDate,
                       const std::vector<Date>& dates,
                       const std::vector<Real>& volatilities,
                       const DayCounter& dayCounter,
                       bool forceMonotoneVariance = true);
    %extend {
      void setInterpolation(const std::string& interpolator = "") {
          std::string s = boost::to_lower_copy(interpolator);
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



// Black smiled surface
%{
using QuantLib::BlackVarianceSurface;
%}

%shared_ptr(BlackVarianceSurface);
class BlackVarianceSurface : public BlackVolTermStructure {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") BlackVarianceSurface;
    #endif
  public:
    enum Extrapolation { ConstantExtrapolation,
                         InterpolatorDefaultExtrapolation };
    %extend {
        BlackVarianceSurface(
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
            BlackVarianceSurface* surface =
                new BlackVarianceSurface(referenceDate,cal,
                                         dates,strikes,
                                         blackVols,dayCounter,lower,upper);
            std::string s = boost::to_lower_copy(interpolator);
            if (s == "" || s == "bilinear") {
                surface->setInterpolation<QuantLib::Bilinear>();
            } else if (s == "bicubic") {
                surface->setInterpolation<QuantLib::Bicubic>();
            } else {
                QL_FAIL("Unknown interpolator: " << interpolator);
            }
            return surface;
        }
        void setInterpolation(const std::string& interpolator = "") {
            std::string s = boost::to_lower_copy(interpolator);
            if (s == "" || s == "bilinear") {
                self->setInterpolation<QuantLib::Bilinear>();
            } else if (s == "bicubic") {
                self->setInterpolation<QuantLib::Bicubic>();
            } else {
                QL_FAIL("Unknown interpolator: " << interpolator);
            }
        }
    }
};



// constant caplet constant term structure
%{
using QuantLib::ConstantOptionletVolatility;
%}

%shared_ptr(ConstantOptionletVolatility);
class ConstantOptionletVolatility : public OptionletVolatilityStructure {
  public:
    ConstantOptionletVolatility(const Date& referenceDate,
                                const Calendar &cal,
                                BusinessDayConvention bdc,
                                Volatility volatility,
                                const DayCounter& dayCounter,
                                const VolatilityType type = ShiftedLognormal,
                                const Real shift = 0.0);
    ConstantOptionletVolatility(const Date& referenceDate,
                                const Calendar &cal,
                                BusinessDayConvention bdc,
                                const Handle<Quote>& volatility,
                                const DayCounter& dayCounter,
                                const VolatilityType type = ShiftedLognormal,
                                const Real shift = 0.0);
    ConstantOptionletVolatility(Natural settlementDays,
                                const Calendar &cal,
                                BusinessDayConvention bdc,
                                Volatility volatility,
                                const DayCounter& dayCounter,
                                const VolatilityType type = ShiftedLognormal,
                                const Real shift = 0.0);
    ConstantOptionletVolatility(Natural settlementDays,
                                const Calendar &cal,
                                BusinessDayConvention bdc,
                                const Handle<Quote>& volatility,
                                const DayCounter& dayCounter,
                                const VolatilityType type = ShiftedLognormal,
                                const Real shift = 0.0);
};



%{
using QuantLib::ConstantSwaptionVolatility;
%}

%shared_ptr(ConstantSwaptionVolatility);
class ConstantSwaptionVolatility : public SwaptionVolatilityStructure {
  public:
    ConstantSwaptionVolatility(Natural settlementDays,
                               const Calendar& cal,
                               BusinessDayConvention bdc,
                               const Handle<Quote>& volatility,
                               const DayCounter& dc,
                               const VolatilityType type = ShiftedLognormal,
                               const Real shift = 0.0);
    ConstantSwaptionVolatility(const Date& referenceDate,
                               const Calendar& cal,
                               BusinessDayConvention bdc,
                               const Handle<Quote>& volatility,
                               const DayCounter& dc,
                               const VolatilityType type = ShiftedLognormal,
                               const Real shift = 0.0);
    ConstantSwaptionVolatility(Natural settlementDays,
                               const Calendar& cal,
                               BusinessDayConvention bdc,
                               Volatility volatility,
                               const DayCounter& dc,
                               const VolatilityType type = ShiftedLognormal,
                               const Real shift = 0.0);
    ConstantSwaptionVolatility(const Date& referenceDate,
                               const Calendar& cal,
                               BusinessDayConvention bdc,
                               Volatility volatility,
                               const DayCounter& dc,
                               const VolatilityType type = ShiftedLognormal,
                               const Real shift = 0.0);
};

%{
using QuantLib::SwaptionVolatilityMatrix;
using QuantLib::SwaptionVolatilityDiscrete;
%}

%shared_ptr(SwaptionVolatilityDiscrete);
class SwaptionVolatilityDiscrete : public SwaptionVolatilityStructure {
    private:
        SwaptionVolatilityDiscrete();
    public:
        const std::vector<Period>& optionTenors() const;
        const std::vector<Date>& optionDates() const;
        const std::vector<Time>& optionTimes() const;
        const std::vector<Period>& swapTenors() const;
        const std::vector<Time>& swapLengths() const;
        const Date optionDateFromTime(Time optionTime) const;
};

%shared_ptr(SwaptionVolatilityMatrix);
class SwaptionVolatilityMatrix : public SwaptionVolatilityDiscrete {
  public:
    SwaptionVolatilityMatrix(const Date& referenceDate,
                             const Calendar& calendar,
                             BusinessDayConvention bdc,
                             const std::vector<Date>& dates,
                             const std::vector<Period>& lengths,
                             const Matrix& vols,
                             const DayCounter& dayCounter,
                             const bool flatExtrapolation = false,
                             const VolatilityType type = ShiftedLognormal,
                             const Matrix& shifts = Matrix());
    SwaptionVolatilityMatrix(const Calendar& calendar,
                             BusinessDayConvention bdc,
                             const std::vector<Period>& optionTenors,
                             const std::vector<Period>& swapTenors,
                             const std::vector<std::vector<Handle<Quote> > >& vols,
                             const DayCounter& dayCounter,
                             const bool flatExtrapolation = false,
                             const VolatilityType type = ShiftedLognormal,
                             const std::vector<std::vector<Real> >& shifts =
                                          std::vector<std::vector<Real> >());
    SwaptionVolatilityMatrix(const Calendar& calendar,
                             BusinessDayConvention bdc,
                             const std::vector<Period>& optionTenors,
                             const std::vector<Period>& swapTenors,
                             const Matrix& vols,
                             const DayCounter& dayCounter,
                             const bool flatExtrapolation = false,
                             const VolatilityType type = ShiftedLognormal,
                             const Matrix& shifts = Matrix());
    %extend {
        SwaptionVolatilityMatrix(const Date& referenceDate,
                                 const std::vector<Date>& dates,
                                 const std::vector<Period>& lengths,
                                 const Matrix& vols,
                                 const DayCounter& dayCounter,
                                 const bool flatExtrapolation = false,
                                 const VolatilityType type = ShiftedLognormal,
                                 const Matrix& shifts = Matrix()) {
            return new SwaptionVolatilityMatrix(referenceDate, QuantLib::NullCalendar(), Following,
                                                dates, lengths, vols, dayCounter,
                                                flatExtrapolation, type, shifts);
        }
    }
    
    %extend {
        std::pair<unsigned int, unsigned int> locate(const Date& optionDate,
                                                     const Period& swapTenor) const {
            auto sizes = self->locate(optionDate, swapTenor);
            return { (unsigned int)sizes.first, (unsigned int)sizes.second };
        }
        std::pair<unsigned int, unsigned int> locate(Time optionTime,
                                                     Time swapLength) const {
            auto sizes = self->locate(optionTime, swapLength);
            return { (unsigned int)sizes.first, (unsigned int)sizes.second };
        }
    }

    VolatilityType volatilityType() const;
};


%{
using QuantLib::SabrSmileSection;
%}

%shared_ptr(SabrSmileSection)
class SabrSmileSection : public SmileSection {
  public:
    SabrSmileSection(const Date& d,
                     Rate forward,
                     const std::vector<Real>& sabrParameters,
                     const Date& referenceDate = Date(),
                     const DayCounter& dc = Actual365Fixed(),
                     Real shift = 0.0,
                     VolatilityType volatilityType = VolatilityType::ShiftedLognormal);
    SabrSmileSection(Time timeToExpiry,
                     Rate forward,
                     const std::vector<Real>& sabrParameters,
                     Real shift = 0.0,
                     VolatilityType volatilityType = VolatilityType::ShiftedLognormal);
    Real alpha() const;
    Real beta() const;
    Real nu() const;
    Real rho() const;
};


%{
using QuantLib::SviSmileSection;
%}

%shared_ptr(SviSmileSection)
class SviSmileSection : public SmileSection {
  public:
    SviSmileSection(Time timeToExpiry, Rate forward, std::vector<Real> sviParameters);
    SviSmileSection(const Date& d,
                    Rate forward,
                    std::vector<Real> sviParameters,
                    const DayCounter& dc = Actual365Fixed());
};

%{
using QuantLib::SviInterpolatedSmileSection;
%}

%shared_ptr(SviInterpolatedSmileSection)
class SviInterpolatedSmileSection : public SmileSection {
  public:
    SviInterpolatedSmileSection(
        const Date& optionDate,
        Handle<Quote> forward,
        const std::vector<Rate>& strikes,
        bool hasFloatingStrikes,
        Handle<Quote> atmVolatility,
        const std::vector<Handle<Quote> >& volHandles,
        Real a,
        Real b,
        Real sigma,
        Real rho,
        Real m,
        bool aIsFixed,
        bool bIsFixed,
        bool sigmaIsFixed,
        bool rhoIsFixed,
        bool mIsFixed,
        bool vegaWeighted = true,
        ext::shared_ptr<EndCriteria> endCriteria = ext::shared_ptr<EndCriteria>(),
        ext::shared_ptr<OptimizationMethod> method = ext::shared_ptr<OptimizationMethod>(),
        const DayCounter& dc = Actual365Fixed());
    SviInterpolatedSmileSection(
        const Date& optionDate,
        const Rate& forward,
        const std::vector<Rate>& strikes,
        bool hasFloatingStrikes,
        const Volatility& atmVolatility,
        const std::vector<Volatility>& vols,
        Real a,
        Real b,
        Real sigma,
        Real rho,
        Real m,
        bool isAFixed,
        bool isBFixed,
        bool isSigmaFixed,
        bool isRhoFixed,
        bool isMFixed,
        bool vegaWeighted = true,
        ext::shared_ptr<EndCriteria> endCriteria = ext::shared_ptr<EndCriteria>(),
        ext::shared_ptr<OptimizationMethod> method = ext::shared_ptr<OptimizationMethod>(),
        const DayCounter& dc = Actual365Fixed());
    Real a() const;
    Real b() const;
    Real sigma() const;
    Real rho() const;
    Real m() const;
    Real rmsError() const;
    Real maxError() const;
    EndCriteria::Type endCriteria() const;
};


%{
using QuantLib::SwaptionVolatilityCube;
using QuantLib::SabrSwaptionVolatilityCube;
using QuantLib::InterpolatedSwaptionVolatilityCube;
%}

%shared_ptr(SwaptionVolatilityCube);
class SwaptionVolatilityCube : public SwaptionVolatilityDiscrete {
    private:
        SwaptionVolatilityCube();
    public:
        Rate atmStrike(const Date& optionDate,
                       const Period& swapTenor) const;
};

%shared_ptr(SabrSwaptionVolatilityCube);
class SabrSwaptionVolatilityCube : public SwaptionVolatilityCube {
  public:
    SabrSwaptionVolatilityCube(
             const Handle<SwaptionVolatilityStructure>& atmVolStructure,
             const std::vector<Period>& optionTenors,
             const std::vector<Period>& swapTenors,
             const std::vector<Spread>& strikeSpreads,
             const std::vector<std::vector<Handle<Quote> > >& volSpreads,
             const ext::shared_ptr<SwapIndex>& swapIndex,
             const ext::shared_ptr<SwapIndex>& shortSwapIndex,
             bool vegaWeightedSmileFit,
             const std::vector<std::vector<Handle<Quote> > >& parametersGuess,
             const std::vector<bool>& isParameterFixed,
             bool isAtmCalibrated,
             const ext::shared_ptr<EndCriteria>& endCriteria
                                           = ext::shared_ptr<EndCriteria>(),
             Real maxErrorTolerance = Null<Real>(),
             const ext::shared_ptr<OptimizationMethod>& optMethod
                                  = ext::shared_ptr<OptimizationMethod>(),
             const Real errorAccept = Null<Real>(),
             const bool useMaxError = false,
             const Size maxGuesses = 50,
             const bool backwardFlat = false,
             const Real cutoffStrike = 0.0001);
    Matrix sparseSabrParameters() const;
    Matrix denseSabrParameters() const;
    Matrix marketVolCube() const;
    Matrix volCubeAtmCalibrated() const;
    %extend {
        ext::shared_ptr<SabrSmileSection> smileSection(Time optionTime, Time swapLength, bool extr = false) const {
            auto base = dynamic_cast<const SwaptionVolatilityStructure*>($self);
            return ext::dynamic_pointer_cast<SabrSmileSection>(base->smileSection(optionTime, swapLength, extr));
        }
        ext::shared_ptr<SabrSmileSection> smileSection(const Period& optionTenor, const Period& swapTenor, bool extr = false) const {
            auto base = dynamic_cast<const SwaptionVolatilityStructure*>($self);
            return ext::dynamic_pointer_cast<SabrSmileSection>(base->smileSection(optionTenor, swapTenor, extr));
        }
    }
};


%shared_ptr(InterpolatedSwaptionVolatilityCube);
class InterpolatedSwaptionVolatilityCube : public SwaptionVolatilityCube {
  public:
    InterpolatedSwaptionVolatilityCube(const Handle<SwaptionVolatilityStructure>& atmVolStructure,
                                       const std::vector<Period>& optionTenors,
                                       const std::vector<Period>& swapTenors,
                                       const std::vector<Spread>& strikeSpreads,
                                       const std::vector<std::vector<Handle<Quote> > >& volSpreads,
                                       const ext::shared_ptr<SwapIndex>& swapIndex,
                                       const ext::shared_ptr<SwapIndex>& shortSwapIndex,
                                       bool vegaWeightedSmileFit);
};


%{
using QuantLib::SpreadedSwaptionVolatility;
%}

%shared_ptr(SpreadedSwaptionVolatility)
class SpreadedSwaptionVolatility : public SwaptionVolatilityStructure {
  public:
    SpreadedSwaptionVolatility(const Handle<SwaptionVolatilityStructure>&,
                               Handle<Quote> spread);
};



%{
using QuantLib::ConstantYoYOptionletVolatility;
%}

%shared_ptr(ConstantYoYOptionletVolatility)
class ConstantYoYOptionletVolatility : public YoYOptionletVolatilitySurface {
  public:
    ConstantYoYOptionletVolatility(Volatility volatility,
                                   Natural settlementDays,
                                   const Calendar &cal,
                                   BusinessDayConvention bdc,
                                   const DayCounter& dc,
                                   const Period& observationLag,
                                   Frequency frequency,
                                   bool indexIsInterpolated,
                                   Real minStrike = -1.0,
                                   Real maxStrike = 100.0);
};


%{
using QuantLib::FlatSmileSection;
%}

%shared_ptr(FlatSmileSection)

class FlatSmileSection : public SmileSection {
  public:
    FlatSmileSection(const Date& d,
                     Volatility vol,
                     const DayCounter& dc,
                     const Date& referenceDate = Date(),
                     Real atmLevel = Null<Rate>(),
                     VolatilityType type = ShiftedLognormal,
                     Real shift = 0.0);
    FlatSmileSection(Time exerciseTime,
                     Volatility vol,
                     const DayCounter& dc,
                     Real atmLevel = Null<Rate>(),
                     VolatilityType type = ShiftedLognormal,
                     Real shift = 0.0);
};

%{
using QuantLib::InterpolatedSmileSection;
using QuantLib::Actual365Fixed;
%}

template<class Interpolator>
class InterpolatedSmileSection : public SmileSection {
  public:
    InterpolatedSmileSection(
               Time expiryTime,
               const std::vector<Rate>& strikes,
               const std::vector<Handle<Quote> >& stdDevHandles,
               const Handle<Quote>& atmLevel,
               const Interpolator& interpolator = Interpolator(),
               const DayCounter& dc = Actual365Fixed(),
               const VolatilityType type = ShiftedLognormal,
               const Real shift = 0.0);
    InterpolatedSmileSection(
               Time expiryTime,
               const std::vector<Rate>& strikes,
               const std::vector<Real>& stdDevs,
               Real atmLevel,
               const Interpolator& interpolator = Interpolator(),
               const DayCounter& dc = Actual365Fixed(),
               const VolatilityType type = ShiftedLognormal,
               const Real shift = 0.0);
    InterpolatedSmileSection(
               const Date& d,
               const std::vector<Rate>& strikes,
               const std::vector<Handle<Quote> >& stdDevHandles,
               const Handle<Quote>& atmLevel,
               const DayCounter& dc = Actual365Fixed(),               
               const Interpolator& interpolator = Interpolator(),
               const Date& referenceDate = Date(),
               const VolatilityType type = ShiftedLognormal,
               const Real shift = 0.0);
    InterpolatedSmileSection(
               const Date& d,
               const std::vector<Rate>& strikes,
               const std::vector<Real>& stdDevs,
               Real atmLevel,
               const DayCounter& dc = Actual365Fixed(),
               const Interpolator& interpolator = Interpolator(),
               const Date& referenceDate = Date(),
               const VolatilityType type = ShiftedLognormal,
               const Real shift = 0.0);
};

%define export_smileinterpolation_curve(Name,Interpolator)
%shared_ptr(InterpolatedSmileSection<Interpolator>)
%template(Name) InterpolatedSmileSection<Interpolator>;
%enddef

export_smileinterpolation_curve(LinearInterpolatedSmileSection, Linear);
export_smileinterpolation_curve(CubicInterpolatedSmileSection, Cubic);
export_smileinterpolation_curve(MonotonicCubicInterpolatedSmileSection, MonotonicCubic);
export_smileinterpolation_curve(SplineCubicInterpolatedSmileSection, SplineCubic);

%{
using QuantLib::KahaleSmileSection;
%}

%shared_ptr(KahaleSmileSection)

class KahaleSmileSection : public SmileSection {
  public:
    KahaleSmileSection(const ext::shared_ptr<SmileSection> source,
                       const Real atm = Null<Real>(),
                       const bool interpolate = false,
                       const bool exponentialExtrapolation = false,
                       const bool deleteArbitragePoints = false,
                       const std::vector<Real> &moneynessGrid = std::vector<Real>(),
                       const Real gap = 1.0E-5,
                       const int forcedLeftIndex = -1,
                       const int forcedRightIndex = QL_MAX_INTEGER);
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

template <class Evaluation>
class ZabrSmileSection : public SmileSection {
  public:
    ZabrSmileSection(
               Time timeToExpiry, 
               Rate forward,
               const std::vector<Real> &zabrParameters,
               const std::vector<Real> &moneyness = std::vector<Real>(),
               const Size fdRefinement = 5);
    ZabrSmileSection(
               const Date &d, 
               Rate forward,
               const std::vector<Real> &zabrParameters,
               const DayCounter &dc = Actual365Fixed(),
               const std::vector<Real> &moneyness = std::vector<Real>(),
               const Size fdRefinement = 5);
};

%define export_zabrsmilesection_curve(Name,Evaluation)
%shared_ptr(ZabrSmileSection<Evaluation>)
%template(Name) ZabrSmileSection<Evaluation>;
%enddef

export_zabrsmilesection_curve(ZabrShortMaturityLognormalSmileSection, ZabrShortMaturityLognormal);
export_zabrsmilesection_curve(ZabrShortMaturityNormalSmileSection, ZabrShortMaturityNormal);
export_zabrsmilesection_curve(ZabrLocalVolatilitySmileSection, ZabrLocalVolatility);
export_zabrsmilesection_curve(ZabrFullFdSmileSection, ZabrFullFd);


template <class Evaluation>
class ZabrInterpolatedSmileSection : public SmileSection {
  public:
    ZabrInterpolatedSmileSection(
               const Date &optionDate, const Handle<Quote> &forward,
               const std::vector<Rate> &strikes, bool hasFloatingStrikes,
               const Handle<Quote> &atmVolatility,
               const std::vector<Handle<Quote> > &volHandles, Real alpha, Real beta,
               Real nu, Real rho, Real gamma, bool isAlphaFixed = false,
               bool isBetaFixed = false, bool isNuFixed = false,
               bool isRhoFixed = false, bool isGammaFixed = false,
               bool vegaWeighted = true,
               const ext::shared_ptr<EndCriteria> &endCriteria =
               ext::shared_ptr<EndCriteria>(),
               const ext::shared_ptr<OptimizationMethod> &method =
               ext::shared_ptr<OptimizationMethod>(),
               const DayCounter &dc = Actual365Fixed());
    ZabrInterpolatedSmileSection(
               const Date &optionDate, const Rate &forward,
               const std::vector<Rate> &strikes, bool hasFloatingStrikes,
               const Volatility &atmVolatility, const std::vector<Volatility> &vols,
               Real alpha, Real beta, Real nu, Real rho, Real gamma,
               bool isAlphaFixed = false, bool isBetaFixed = false,
               bool isNuFixed = false, bool isRhoFixed = false,
               bool isGammaFixed = false, bool vegaWeighted = true,
               const ext::shared_ptr<EndCriteria> &endCriteria =
               ext::shared_ptr<EndCriteria>(),
               const ext::shared_ptr<OptimizationMethod> &method =
               ext::shared_ptr<OptimizationMethod>(),
               const DayCounter &dc = Actual365Fixed());
    Real alpha() const;
    Real beta() const;
    Real nu() const;
    Real rho() const;
    Real rmsError() const;
    Real maxError() const;
    EndCriteria::Type endCriteria() const;
};

%define export_zabrinterpolatedsmilesection_curve(Name,Evaluation)
%shared_ptr(ZabrInterpolatedSmileSection<Evaluation>)
%template(Name) ZabrInterpolatedSmileSection<Evaluation>;
%enddef

export_zabrinterpolatedsmilesection_curve(ZabrShortMaturityLognormalInterpolatedSmileSection, ZabrShortMaturityLognormal);
export_zabrinterpolatedsmilesection_curve(ZabrShortMaturityNormalInterpolatedSmileSection, ZabrShortMaturityNormal);
export_zabrinterpolatedsmilesection_curve(ZabrLocalVolatilityInterpolatedSmileSection, ZabrLocalVolatility);
export_zabrinterpolatedsmilesection_curve(ZabrFullFdInterpolatedSmileSection, ZabrFullFd);


%shared_ptr(NoArbSabrSmileSection)

class NoArbSabrSmileSection : public SmileSection {
  public:
    NoArbSabrSmileSection(
               Time timeToExpiry, 
               Rate forward,
               const std::vector<Real> &sabrParameters,
               const Real shift = 0.0,
               VolatilityType volatilityType = VolatilityType::ShiftedLognormal);
    NoArbSabrSmileSection(
               const Date &d, 
               Rate forward,
               const std::vector<Real> &sabrParameters,
               const DayCounter &dc = Actual365Fixed(),
               const Real shift = 0.0,
               VolatilityType volatilityType = VolatilityType::ShiftedLognormal);
};

%shared_ptr(NoArbSabrInterpolatedSmileSection)

class NoArbSabrInterpolatedSmileSection : public SmileSection {
  public:
    NoArbSabrInterpolatedSmileSection(
               const Date &optionDate, const Handle<Quote> &forward,
               const std::vector<Rate> &strikes, bool hasFloatingStrikes,
               const Handle<Quote> &atmVolatility,
               const std::vector<Handle<Quote> > &volHandles, Real alpha, Real beta,
               Real nu, Real rho, bool isAlphaFixed = false,
               bool isBetaFixed = false, bool isNuFixed = false,
               bool isRhoFixed = false,
               bool vegaWeighted = true,
               const ext::shared_ptr<EndCriteria> &endCriteria =
               ext::shared_ptr<EndCriteria>(),
               const ext::shared_ptr<OptimizationMethod> &method =
               ext::shared_ptr<OptimizationMethod>(),
               const DayCounter &dc = Actual365Fixed());
    NoArbSabrInterpolatedSmileSection(
               const Date &optionDate, const Rate &forward,
               const std::vector<Rate> &strikes, bool hasFloatingStrikes,
               const Volatility &atmVolatility, const std::vector<Volatility> &vols,
               Real alpha, Real beta, Real nu, Real rho,
               bool isAlphaFixed = false, bool isBetaFixed = false,
               bool isNuFixed = false, bool isRhoFixed = false,
               bool vegaWeighted = true,
               const ext::shared_ptr<EndCriteria> &endCriteria =
               ext::shared_ptr<EndCriteria>(),
               const ext::shared_ptr<OptimizationMethod> &method =
               ext::shared_ptr<OptimizationMethod>(),
               const DayCounter &dc = Actual365Fixed());
    Real alpha() const;
    Real beta() const;
    Real nu() const;
    Real rho() const;
    Real rmsError() const;
    Real maxError() const;
    EndCriteria::Type endCriteria() const;
};

%{
using QuantLib::AndreasenHugeVolatilityInterpl;
using QuantLib::AndreasenHugeVolatilityAdapter;
using QuantLib::AndreasenHugeLocalVolAdapter;
using QuantLib::HestonBlackVolSurface;
%}

%template(CalibrationErrorTuple) std::tuple<Real, Real, Real>;

%shared_ptr(AndreasenHugeVolatilityInterpl)
class AndreasenHugeVolatilityInterpl : public Observable {
  public:
        enum InterpolationType {PiecewiseConstant, Linear, CubicSpline};
        enum CalibrationType {
            // we specify values directly to work around a problem in
            // the SWIG C# module
            Call = 1, // Option::Call,
            Put = -1, // Option::Put,
            CallPut};

        typedef std::vector<std::pair<
            ext::shared_ptr<VanillaOption>, ext::shared_ptr<Quote> > >
          CalibrationSet;

        AndreasenHugeVolatilityInterpl(
            const CalibrationSet& calibrationSet,
            const Handle<Quote>& spot,
            const Handle<YieldTermStructure>& rTS,
            const Handle<YieldTermStructure>& qTS,
            InterpolationType interpolationType = CubicSpline,
            CalibrationType calibrationType = Call,
            Size nGridPoints = 500,
            Real minStrike = Null<Real>(),
            Real maxStrike = Null<Real>(),
            const ext::shared_ptr<OptimizationMethod>& optimizationMethod =
                ext::shared_ptr<OptimizationMethod>(new LevenbergMarquardt),
            const EndCriteria& endCriteria =
                EndCriteria(500, 100, 1e-12, 1e-10, 1e-10));

        Date maxDate() const;
        Real minStrike() const;
        Real maxStrike() const;

        Real fwd(Time t) const;
        const Handle<YieldTermStructure>& riskFreeRate() const;

        // returns min, max and average error in volatility units
        std::tuple<Real, Real, Real> calibrationError() const;

        // returns the option price of the calibration type. In case
        // of CallPut it return the call option price
        Real optionPrice(Time t, Real strike, Option::Type optionType) const;

        Volatility localVol(Time t, Real strike) const;
};

%shared_ptr(AndreasenHugeVolatilityAdapter)
class AndreasenHugeVolatilityAdapter : public BlackVolTermStructure {
  public:
    AndreasenHugeVolatilityAdapter(
        const ext::shared_ptr<AndreasenHugeVolatilityInterpl>& volInterpl,
        Real eps = 1e-6);
};

%shared_ptr(AndreasenHugeLocalVolAdapter)
class AndreasenHugeLocalVolAdapter : public LocalVolTermStructure {
  public:
    explicit AndreasenHugeLocalVolAdapter(
        const ext::shared_ptr<AndreasenHugeVolatilityInterpl>& localVol);
};

%shared_ptr(HestonBlackVolSurface)
class HestonBlackVolSurface : public BlackVolTermStructure {
  public:
    explicit HestonBlackVolSurface(
        const Handle<HestonModel>& hestonModel,
        const AnalyticHestonEngine::ComplexLogFormula cpxLogFormula
            = AnalyticHestonEngine::Gatheral,
        const AnalyticHestonEngine::Integration& integration =
            AnalyticHestonEngine::Integration::gaussLaguerre(164));
};

%{
using QuantLib::CmsMarket;
%}

%shared_ptr(CmsMarket)
class CmsMarket{
  public:       
    CmsMarket(
        const std::vector<Period>& swapLengths,
        const std::vector<ext::shared_ptr<SwapIndex> >& swapIndexes,
        const ext::shared_ptr<IborIndex>& iborIndex,
        const std::vector<std::vector<Handle<Quote> > >& bidAskSpreads,
        const std::vector<ext::shared_ptr<CmsCouponPricer> >& pricers,
        const Handle<YieldTermStructure>& discountingTS);

        void reprice(const Handle<SwaptionVolatilityStructure>& volStructure,
                     Real meanReversion);

        const std::vector<Period>& swapTenors() const;
        const std::vector<Period>& swapLengths() const;
        const Matrix& impliedCmsSpreads();
        const Matrix& spreadErrors();
        Matrix browse() const;

        Real weightedSpreadError(const Matrix& weights);
        Real weightedSpotNpvError(const Matrix& weights);
        Real weightedFwdNpvError(const Matrix& weights);
        Array weightedSpreadErrors(const Matrix& weights);
        Array weightedSpotNpvErrors(const Matrix& weights);
        Array weightedFwdNpvErrors(const Matrix& weights);
};

%{
using QuantLib::CmsMarketCalibration;
%}

class CmsMarketCalibration {
  public:
    enum CalibrationType {OnSpread, OnPrice, OnForwardCmsPrice };

    CmsMarketCalibration(
        Handle<SwaptionVolatilityStructure>& volCube,
        ext::shared_ptr<CmsMarket>& cmsMarket,
        const Matrix& weights,
        CalibrationType calibrationType);

    Array compute(const ext::shared_ptr<EndCriteria>& endCriteria,
              const ext::shared_ptr<OptimizationMethod>& method,
              const Array& guess,
              bool isMeanReversionFixed);

    Matrix compute(const ext::shared_ptr<EndCriteria>& endCriteria,
                  const ext::shared_ptr<OptimizationMethod>& method,
                  const Matrix& guess,
                  bool isMeanReversionFixed,
                  const Real meanReversionGuess = Null<Real>());


    Matrix computeParametric(const ext::shared_ptr<EndCriteria> &endCriteria,
                      const ext::shared_ptr<OptimizationMethod> &method,
                      const Matrix &guess, bool isMeanReversionFixed,
                      const Real meanReversionGuess = Null<Real>());

    Real error();
    EndCriteria::Type endCriteria();
};

#endif
