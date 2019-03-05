
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
 <http://quantlib.org/license.shtml>.

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


// eventually the classes exported here will be redesigned or deprecated

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
    OptionletStripper1(const boost::shared_ptr<CapFloorTermVolSurface>& parVolSurface,
                       const boost::shared_ptr<IborIndex> &index,
                       Rate switchStrikes = Null<Rate>(),
                       Real accuracy = 1.0e-6, Natural maxIter = 100,
                       const Handle<YieldTermStructure> &discount =
                              Handle<YieldTermStructure>(),
                       VolatilityType type = ShiftedLognormal,
                       Real displacement = 0.0,
                       bool dontThrow = false);
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
    StrippedOptionletAdapter(const boost::shared_ptr<StrippedOptionletBase>&);
};
 
#endif
