
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007 StatPro Italia srl
 Copyright (C) 2005 Johan Witters
 Copyright (C) 2018 Matthias Groncki
 Copyright (C) 2023 Skandinaviska Enskilda Banken AB (publ)

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

#ifndef quantlib_calendar_i
#define quantlib_calendar_i

%include common.i
%include date.i
%include stl.i

%define QL_TYPECHECK_BUSINESSDAYCONVENTION       6210    %enddef

%{
using QuantLib::Calendar;
%}

%{
using QuantLib::BusinessDayConvention;
using QuantLib::Following;
using QuantLib::ModifiedFollowing;
using QuantLib::Preceding;
using QuantLib::ModifiedPreceding;
using QuantLib::Unadjusted;
using QuantLib::HalfMonthModifiedFollowing;
using QuantLib::Nearest;
%}

enum BusinessDayConvention {
    Following,
    ModifiedFollowing,
    Preceding,
    ModifiedPreceding,
    Unadjusted,
    HalfMonthModifiedFollowing,
    Nearest
};

%{
using QuantLib::JointCalendarRule;
using QuantLib::JoinHolidays;
using QuantLib::JoinBusinessDays;
%}

enum JointCalendarRule { JoinHolidays, JoinBusinessDays };

#if defined(SWIGPYTHON)
%typemap(in) ext::optional<BusinessDayConvention> %{
    if ($input == Py_None)
        $1 = ext::nullopt;
    else if (PyLong_Check($input))
        $1 = (BusinessDayConvention)PyLong_AsLong($input);
    else
        SWIG_exception(SWIG_TypeError, "int expected");
%}
%typecheck (QL_TYPECHECK_BUSINESSDAYCONVENTION) ext::optional<BusinessDayConvention> %{
    $1 = (PyLong_Check($input) || $input == Py_None) ? 1 : 0;
%}
#endif

class Calendar {
  protected:
    Calendar();
  public:
    bool isWeekend(Weekday w);
    Date startOfMonth(const Date&);
    Date endOfMonth(const Date&);
    bool isBusinessDay(const Date&);
    bool isHoliday(const Date&);
    bool isEndOfMonth(const Date&);
    bool isStartOfMonth(const Date&);
    void addHoliday(const Date&);
    void removeHoliday(const Date&);
    void resetAddedAndRemovedHolidays();

    Date adjust(const Date& d,
                BusinessDayConvention convention = QuantLib::Following);
    Date advance(const Date& d, Integer n, TimeUnit unit,
                 BusinessDayConvention convention = QuantLib::Following,
                 bool endOfMonth = false);
    Date advance(const Date& d, const Period& period,
                 BusinessDayConvention convention = QuantLib::Following,
                 bool endOfMonth = false);
    BigInteger businessDaysBetween(const Date& from,
                                   const Date& to,
                                   bool includeFirst = true,
                                   bool includeLast = false);
    std::vector<Date> holidayList(const Date& from,
                                  const Date& to,
                                  bool includeWeekEnds = false);
    std::vector<Date> businessDayList(const Date& from,
                                      const Date& to);
    std::string name();
    bool empty();
    %extend {
        std::string __str__() {
            return self->name()+" calendar";
        }
        #if defined(SWIGPYTHON) || defined(SWIGJAVA)
        bool operator==(const Calendar& other) {
            return (*self) == other;
        }
        bool operator!=(const Calendar& other) {
            return (*self) != other;
        }
        hash_t __hash__() {
            return self->empty() ? 0 : std::hash<std::string>()(self->name());
        }
        #endif
    }
};

namespace std {
    %template(CalendarVector) vector<Calendar>;
}

namespace QuantLib {

    class Argentina : public Calendar {
      public:
        enum Market { Merval };
        Argentina(Market m = Merval);
    };

    class Australia : public Calendar {
      public:
        enum Market { Settlement, ASX };
        Australia(Market market = Settlement);
    };

    class Austria : public Calendar {
      public:
        enum Market { Settlement, Exchange };
        Austria(Market m = Settlement);
    };

    class Botswana : public Calendar {};

    class Brazil : public Calendar {
      public:
        enum Market { Settlement, Exchange };
        Brazil(Market m = Settlement);
    };

    class Canada : public Calendar {
      public:
        enum Market { Settlement, TSX };
        Canada(Market m = Settlement);
    };

    class Chile : public Calendar {
      public:
        enum Market { SSE };
        Chile(Market m = SSE);
    };

    class China : public Calendar {
      public:
        enum Market { SSE, IB };
        China(Market m = SSE);
    };

    class CzechRepublic : public Calendar {
      public:
        enum Market { PSE };
        CzechRepublic(Market m = PSE);
    };

    class Denmark : public Calendar {};
    class Finland : public Calendar {};

    class France : public Calendar {
      public:
        enum Market { Settlement, Exchange };
        France(Market m = Settlement);
    };

    class Germany : public Calendar {
      public:
        enum Market { Settlement, FrankfurtStockExchange, Xetra, Eurex };
        Germany(Market m = FrankfurtStockExchange);
    };

    class HongKong : public Calendar {
      public:
        enum Market { HKEx };
        HongKong(Market m = HKEx);
    };

    class Hungary : public Calendar {};

    class Iceland : public Calendar {
      public:
        enum Market { ICEX };
        Iceland(Market m = ICEX);
    };

    class India : public Calendar {
      public:
        enum Market { NSE };
        India(Market m = NSE);
    };

    class Indonesia : public Calendar {
      public:
        enum Market { BEJ, JSX };
        Indonesia(Market m = BEJ);
    };

    class Israel : public Calendar {
      public:
        enum Market { Settlement, TASE, SHIR };
        Israel(Market m = Settlement);
    };

    class Italy : public Calendar {
      public:
        enum Market { Settlement, Exchange };
        Italy(Market m = Settlement);
    };

    class Japan : public Calendar {};

    class Mexico : public Calendar {
      public:
        enum Market { BMV };
        Mexico(Market m = BMV);
    };

    class NewZealand : public Calendar {
      public:
        enum Market { Wellington, Auckland };
        NewZealand(Market m = Wellington);
    };

    class Norway : public Calendar {};

    class Poland : public Calendar {
      public:
        enum Market { Settlement, WSE };
        Poland(Market m = Settlement);
    };

    class Romania : public Calendar {
      public:
        enum Market { Public, BVB };
        Romania(Market m = BVB);
    };

    class Russia : public Calendar {
      public:
        enum Market { Settlement, MOEX };
        Russia(Market m = Settlement);
    };

    class SaudiArabia : public Calendar {
      public:
        enum Market { Tadawul };
        SaudiArabia(Market m = Tadawul);
    };

    class Singapore : public Calendar {
      public:
        enum Market { SGX };
        Singapore(Market m = SGX);
    };

    class Slovakia : public Calendar {
      public:
        enum Market { BSSE };
        Slovakia(Market m = BSSE);
    };

    class SouthAfrica : public Calendar {};

    class SouthKorea : public Calendar {
      public:
        enum Market { Settlement, KRX };
        SouthKorea(Market m = KRX);
    };

    class Sweden : public Calendar {};
    class Switzerland : public Calendar {};

    class Taiwan : public Calendar {
      public:
        enum Market { TSEC };
        Taiwan(Market m = TSEC);
    };

    class TARGET : public Calendar {};
    class Thailand : public Calendar {};
    class Turkey : public Calendar {};

    class Ukraine : public Calendar {
      public:
        enum Market { USE };
        Ukraine(Market m = USE);
    };

    class UnitedKingdom : public Calendar {
      public:
        enum Market { Settlement, Exchange, Metals };
        UnitedKingdom(Market m = Settlement);
    };

    class UnitedStates : public Calendar {
      public:
        enum Market { Settlement, NYSE, GovernmentBond,
                      NERC, LiborImpact, FederalReserve, SOFR };
        UnitedStates(Market m);
    };

    // others

    class NullCalendar : public Calendar {};

    class WeekendsOnly : public Calendar {};

    class JointCalendar : public Calendar {
      public:
        JointCalendar(const Calendar&, const Calendar&,
                      JointCalendarRule rule = QuantLib::JoinHolidays);
        JointCalendar(const Calendar&, const Calendar&, const Calendar&,
                      JointCalendarRule rule = QuantLib::JoinHolidays);
        JointCalendar(const Calendar&, const Calendar&,
                      const Calendar&, const Calendar&,
                      JointCalendarRule rule = QuantLib::JoinHolidays);
        explicit JointCalendar(const std::vector<Calendar>&,
                               JointCalendarRule = QuantLib::JoinHolidays);
    };

    class BespokeCalendar : public Calendar {
      public:
        BespokeCalendar(const std::string& name);
        void addWeekend(Weekday);
    };

}


#endif

