
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2007, 2008 StatPro Italia srl

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

#ifndef quantlib_scheduler_i
#define quantlib_scheduler_i

%include date.i
%include calendars.i
%include types.i

%define QL_TYPECHECK_DATEGENERATION       7210    %enddef

%{
using QuantLib::Schedule;
using QuantLib::DateGeneration;
using QuantLib::MakeSchedule;
%}

struct DateGeneration {
    enum Rule { Backward, Forward,
                Zero, ThirdWednesday,
                Twentieth, TwentiethIMM,
                OldCDS, CDS, CDS2015 };
};
#if defined(SWIGPYTHON)
%typemap(in) boost::optional<DateGeneration::Rule> %{
    if($input == Py_None)
        $1 = boost::none;
    else if (PyInt_Check($input))
        $1 = (DateGeneration::Rule) PyInt_AsLong($input);
    else
        $1 = (DateGeneration::Rule) PyLong_AsLong($input);
%}
%typecheck (QL_TYPECHECK_DATEGENERATION) boost::optional<DateGeneration::Rule> {
if (PyInt_Check($input) || PyLong_Check($input) || Py_None == $input)
    $1 = 1;
else
    $1 = 0;
}
#endif

#if defined(SWIGRUBY)
%mixin Schedule "Enumerable";
#endif
class Schedule {
    #if defined(SWIGPYTHON) || defined(SWIGRUBY)
    %rename(__len__)       size;
    %ignore                date;
    #endif
    #if defined(SWIGRUBY)
    %rename("isRegular?")  isRegular;
    #endif
  public:
    #if defined(SWIGPYTHON)
    Schedule(const std::vector<Date>&,
             const Calendar& calendar = NullCalendar(),
             const BusinessDayConvention convention = Unadjusted,
             boost::optional<BusinessDayConvention>
             terminationDateConvention = boost::none,
             const boost::optional<Period> tenor = boost::none,
             boost::optional<DateGeneration::Rule> rule = boost::none,
             boost::optional<bool> endOfMonth = boost::none,
             const std::vector<bool>& isRegular = std::vector<bool>(0));
    #else
    Schedule(const std::vector<Date>&,
         const Calendar& calendar = NullCalendar(),
         const BusinessDayConvention convention = Unadjusted);
    #endif
    Schedule(const Date& effectiveDate,
             const Date& terminationDate,
             const Period& tenor,
             const Calendar& calendar,
             BusinessDayConvention convention,
             BusinessDayConvention terminationDateConvention,
             DateGeneration::Rule rule,
             bool endOfMonth,
             const Date& firstDate = Date(),
             const Date& nextToLastDate = Date());
    Schedule();
    Size size() const;
    Date date(Size i) const;
    bool isRegular(Size i) const;
    Schedule until(Date truncationDate) const;
    %extend {
        #if defined(SWIGPYTHON) || defined(SWIGRUBY)
        Date __getitem__(Integer i) {
            Integer size_ = static_cast<Integer>(self->size());
            if (i>=0 && i<size_) {
                return self->date(i);
            } else if (i<0 && -i<=size_) {
                return self->date(size_+i);
            } else {
                throw std::out_of_range("schedule index out of range");
            }
        }
        #endif
        #if defined(SWIGRUBY)
        void each() {
            for (Size i=0; i<self->size(); i++) {
                Date* d = new Date(self->date(i));
                rb_yield(SWIG_NewPointerObj((void *) d,
                                            $descriptor(Date *), 1));
            }
        }
        #endif
    }
};

#if defined(SWIGPYTHON)
%rename (_MakeSchedule) MakeSchedule;
#endif

/*! This class provides a more comfortable interface to the
    argument list of Schedule's constructor.
*/
class MakeSchedule {
  public:
    MakeSchedule();
#if defined(SWIGPYTHON)
    // 'from' can't be overridden in python
    %rename("fromDate") from;
#endif
    MakeSchedule& from(const Date& effectiveDate);
    MakeSchedule& to(const Date& terminationDate);
    MakeSchedule& withTenor(const Period&);
    MakeSchedule& withFrequency(Frequency);
    MakeSchedule& withCalendar(const Calendar&);
    MakeSchedule& withConvention(BusinessDayConvention);
    MakeSchedule& withTerminationDateConvention(BusinessDayConvention);
    MakeSchedule& withRule(DateGeneration::Rule);
    MakeSchedule& forwards();
    MakeSchedule& backwards();
    MakeSchedule& endOfMonth(bool flag=true);
    MakeSchedule& withFirstDate(const Date& d);
    MakeSchedule& withNextToLastDate(const Date& d);

    %extend{
      Schedule schedule(){
        return (Schedule)(* $self);
      }
    }
};

#if defined(SWIGPYTHON)
%pythoncode{
def MakeSchedule(effectiveDate=None,terminationDate=None,tenor=None,
    frequency=None,calendar=None,convention=None,terminalDateConvention=None,
    rule=None,forwards=False,backwards=False,
    endOfMonth=None,firstDate=None,nextToLastDate=None):
    ms = _MakeSchedule()
    if effectiveDate is not None:
        ms.fromDate(effectiveDate)
    if terminationDate is not None:
        ms.to(terminationDate)
    if tenor is not None:
        ms.withTenor(tenor)
    if frequency is not None:
        ms.withFrequency(frequency)
    if calendar is not None:
        ms.withCalendar(calendar)
    if convention is not None:
        ms.withConvention(convention)
    if terminalDateConvention is not None:
        ms.withTerminationDateConvention(terminalDateConvention)
    if rule is not None:
        ms.withRule(rule)
    if forwards:
        ms.forwards()
    if backwards:
        ms.backwards()
    if endOfMonth is not None:
        ms.endOfMonth(endOfMonth)
    if firstDate is not None:
        ms.withFirstDate(firstDate)
    if nextToLastDate is not None:
        ms.withNextToLastDate(nextToLastDate)
    return ms.schedule()
}
#endif

#endif
