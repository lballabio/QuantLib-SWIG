
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2007 StatPro Italia srl
 Copyright (C) 2015 Matthias Groncki

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

#ifndef quantlib_old_volatility_i
#define quantlib_old_volatility_i

%include common.i
%include date.i
%include calendars.i
%include daycounters.i
%include vectors.i
%include linearalgebra.i
%include interpolation.i
%include observer.i
%include volatilities.i


// cap/floor volatilities

%{
using QuantLib::CapFloorTermVolatilityStructure;
%}

%shared_ptr(CapFloorTermVolatilityStructure);
class CapFloorTermVolatilityStructure : public VolatilityTermStructure {
  private:
    CapFloorTermVolatilityStructure();
  public:
    Volatility volatility(const Period& length, Rate strike,
                          bool extrapolate = false);
    Volatility volatility(const Date& end, Rate strike,
                          bool extrapolate = false);
    Volatility volatility(Time end, Rate strike,
                          bool extrapolate = false);
};

%template(CapFloorTermVolatilityStructureHandle) Handle<CapFloorTermVolatilityStructure>;
%template(RelinkableCapFloorTermVolatilityStructureHandle)
RelinkableHandle<CapFloorTermVolatilityStructure>;

%{
using QuantLib::CapFloorTermVolCurve;
%}

%shared_ptr(CapFloorTermVolCurve);
class CapFloorTermVolCurve : public CapFloorTermVolatilityStructure {
  public:
    CapFloorTermVolCurve(const Date& referenceDate,
                         const Calendar& calendar,
                         BusinessDayConvention bdc,
                         const std::vector<Period>& lengths,
                         const std::vector<Volatility>& vols,
                         const DayCounter& dc =
                                           QuantLib::Actual365Fixed());
    CapFloorTermVolCurve(Natural settlementDays,
                         const Calendar& calendar,
                         BusinessDayConvention bdc,
                         const std::vector<Period>& lengths,
                         const std::vector<Volatility>& vols,
                         const DayCounter& dc =
                                            QuantLib::Actual365Fixed());
};

%{
using QuantLib::CapFloorTermVolSurface;
%}

%shared_ptr(CapFloorTermVolSurface);
class CapFloorTermVolSurface : public CapFloorTermVolatilityStructure {
  public:
    CapFloorTermVolSurface(Natural settlementDays,
                           const Calendar& calendar,
                           BusinessDayConvention bdc,
                           const std::vector<Period>& optionTenors,
                           const std::vector<Rate>& strikes,
                           const std::vector<std::vector<Handle<Quote> > >& quotes,
                           const DayCounter& dc = QuantLib::Actual365Fixed());
    CapFloorTermVolSurface(const Date& settlementDate,
                           const Calendar& calendar,
                           BusinessDayConvention bdc,
                           const std::vector<Period>& optionTenors,
                           const std::vector<Rate>& strikes,
                           const std::vector<std::vector<Handle<Quote> > >& quotes,
                           const DayCounter& dc = QuantLib::Actual365Fixed());
    CapFloorTermVolSurface(const Date& settlementDate,
                           const Calendar& calendar,
                           BusinessDayConvention bdc,
                           const std::vector<Period>& optionTenors,
                           const std::vector<Rate>& strikes,
                           const Matrix& volatilities,
                           const DayCounter& dc = QuantLib::Actual365Fixed());
    CapFloorTermVolSurface(Natural settlementDays,
                           const Calendar& calendar,
                           BusinessDayConvention bdc,
                           const std::vector<Period>& optionTenors,
                           const std::vector<Rate>& strikes,
                           const Matrix& volatilities,
                           const DayCounter& dc = QuantLib::Actual365Fixed());
    Date maxDate() const;
    Real minStrike() const;
    Real maxStrike() const;
    const std::vector<Period>& optionTenors() const;
    const std::vector<Date>& optionDates() const;
    const std::vector<Time>& optionTimes() const;
    const std::vector<Rate>& strikes() const;
};

%{
using QuantLib::StrippedOptionletBase;
using QuantLib::VolatilityType;
%}

%shared_ptr(StrippedOptionletBase);
class StrippedOptionletBase {
  private:
    StrippedOptionletBase();
  public:
    const std::vector<Rate>& optionletStrikes(Size i);
    const std::vector<Volatility>& optionletVolatilities(Size i);
    const std::vector<Date>& optionletFixingDates();
    const std::vector<Time>& optionletFixingTimes();
    Size optionletMaturities();
    const std::vector<Rate>& atmOptionletRates();
    DayCounter dayCounter();
    Calendar calendar();
    Natural settlementDays();
    BusinessDayConvention businessDayConvention();
    VolatilityType volatilityType();
    Real displacement();
};

%{
using QuantLib::StrippedOptionlet;
%}

%shared_ptr(StrippedOptionlet)
class StrippedOptionlet : public StrippedOptionletBase {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") StrippedOptionlet;
    #endif
  public:
    StrippedOptionlet(Natural settlementDays,
                      const Calendar& calendar,
                      BusinessDayConvention bdc,
                      ext::shared_ptr<IborIndex> iborIndex,
                      const std::vector<Date>& optionletDates,
                      const std::vector<Rate>& strikes,
                      std::vector<std::vector<Handle<Quote> > > volatilities,
                      DayCounter dc,
                      VolatilityType type = ShiftedLognormal,
                      Real displacement = 0.0);
};

%{
using QuantLib::OptionletStripper1;
%}

%shared_ptr(OptionletStripper1);
class OptionletStripper1 : public StrippedOptionletBase {
    #if !defined(SWIGJAVA) && !defined(SWIGCSHARP)
    %feature("kwargs") OptionletStripper1;
    #endif
  public:
    %extend {
        OptionletStripper1(const ext::shared_ptr<CapFloorTermVolSurface>& parVolSurface,
                           const ext::shared_ptr<IborIndex> &index,
                           Rate switchStrikes = Null<Rate>(),
                           Real accuracy = 1.0e-6, Natural maxIter = 100,
                           const Handle<YieldTermStructure> &discount = {},
                           VolatilityType type = ShiftedLognormal,
                           Real displacement = 0.0,
                           bool dontThrow = false,
                           Period optionletFrequency = Period()) {
            ext::optional<Period> frequency = ext::nullopt;
            if (optionletFrequency != Period())
                frequency = optionletFrequency;
            return new OptionletStripper1(parVolSurface, index, switchStrikes, accuracy, maxIter,
                                          discount, type, displacement, dontThrow, frequency);
        }
    }
    const Matrix& capFloorPrices() const;
    const Matrix& capFloorVolatilities() const;
    const Matrix& optionletPrices() const;
    Rate switchStrike() const;
};


%{
using QuantLib::StrippedOptionletAdapter;
%}

%shared_ptr(StrippedOptionletAdapter);
class StrippedOptionletAdapter : public OptionletVolatilityStructure {
  public:
    StrippedOptionletAdapter(const ext::shared_ptr<StrippedOptionletBase>&);
};
 
#endif
