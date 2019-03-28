
/*
 Copyright (C) 2003 StatPro Italia srl
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

#ifndef quantlib_exercise_i
#define quantlib_exercise_i

%include common.i
%include boost_shared_ptr.i

// exercise conditions

%{
using QuantLib::Exercise;
%}

%shared_ptr(Exercise)
class Exercise {
  public:
    enum Type {
        American, Bermudan, European
    };
    explicit Exercise(Type type);
    Type type() const;
    Date date(Size index);
    Date dateAt(Size index);
    const std::vector<Date>& dates();
    Date lastDate() const;
    #if defined(SWIGJAVA)
    // Scala can't use "type" as a method name
    %extend {
        Exercise::Type exerciseType() {
            return self->type();
        }
    }
    #endif
};

%{
using QuantLib::EuropeanExercise;
using QuantLib::AmericanExercise;
using QuantLib::BermudanExercise;
using QuantLib::RebatedExercise;
using QuantLib::SwingExercise;
%}

%shared_ptr(EuropeanExercise)
class EuropeanExercise : public Exercise {
  public:
    EuropeanExercise(const Date& date);
};

%shared_ptr(AmericanExercise)
class AmericanExercise : public Exercise {
  public:
    AmericanExercise(const Date& earliestDate,
                        const Date& latestDate,
                        bool payoffAtExpiry = false);
};

%shared_ptr(BermudanExercise)
class BermudanExercise : public Exercise {
  public:
    BermudanExercise(const std::vector<Date>& dates,
                        bool payoffAtExpiry = false);
};

%{
using QuantLib::BusinessDayConvention;
using QuantLib::Calendar;
using QuantLib::Following;
using QuantLib::NullCalendar;
%}

%shared_ptr(RebatedExercise)
class RebatedExercise : public Exercise {
  public:
    RebatedExercise(const Exercise &exercise,
                       const std::vector<Real> rebates,
                       Natural rebateSettlementDays = 0,
                       const Calendar &rebatePaymentCalendar = NullCalendar(),
                       const BusinessDayConvention rebatePaymentConvention = Following);
};


%shared_ptr(SwingExercise)
class SwingExercise : public Exercise {
  public:
    SwingExercise(const std::vector<Date>& dates);
};

#endif
