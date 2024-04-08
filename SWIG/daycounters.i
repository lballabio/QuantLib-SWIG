
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005 StatPro Italia srl
 Copyright (C) 2005 Johan Witters
 Copyright (C) 2022 Ignacio Anguita

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

#ifndef quantlib_day_counters_i
#define quantlib_day_counters_i

%include common.i
%include calendars.i
%include date.i
%include types.i
%include stl.i
%include null.i

%{
using QuantLib::DayCounter;
%}

class DayCounter {
  protected:
    DayCounter();
  public:
    BigInteger dayCount(const Date& d1, const Date& d2) const;
    Time yearFraction(const Date& d1, const Date& d2,
                      const Date& startRef = Date(),
                      const Date& endRef = Date()) const;
    std::string name() const;
    bool empty();
    %extend {
        std::string __str__() {
            return self->name()+" day counter";
        }
        #if defined(SWIGPYTHON) || defined(SWIGJAVA)
        bool operator==(const DayCounter& other) {
            return (*self) == other;
        }
        bool operator!=(const DayCounter& other) {
            return (*self) != other;
        }
        hash_t __hash__() {
            return self->empty() ? 0 : std::hash<std::string>()(self->name());
        }
        #endif
    }
};

namespace QuantLib {

    class Actual360 : public DayCounter {
      public:
        Actual360(const bool includeLastDay = false);
    };
    class Actual366 : public DayCounter {
      public:
        Actual366(const bool includeLastDay = false);
    };
    class Actual36525 : public DayCounter {
      public:
        Actual36525(const bool includeLastDay = false);
    };
    class Actual364 : public DayCounter {};
    class Actual365Fixed : public DayCounter {
      public:
        enum Convention { Standard, Canadian, NoLeap };
        Actual365Fixed(Convention c = Standard);
    };
    class Thirty360 : public DayCounter {
      public:
        enum Convention { USA, BondBasis, European, EurobondBasis, Italian, German, ISMA, ISDA, NASD };
        Thirty360(Convention c, const Date& terminationDate = Date());
    };
    class Thirty365 : public DayCounter {};
    class ActualActual : public DayCounter {
      public:
        enum Convention { ISMA, Bond, ISDA, Historical, Actual365, AFB, Euro };
        ActualActual(Convention c, const Schedule& schedule = Schedule());
    };
    class OneDayCounter : public DayCounter {};
    class SimpleDayCounter : public DayCounter {};
    class Business252 : public DayCounter {
      public:
        Business252(Calendar c = Brazil());
    };

}

%{
using QuantLib::yearFractionToDate;
%}

Date yearFractionToDate(
    const DayCounter& dayCounter, const Date& referenceDate, Time t);

#endif
