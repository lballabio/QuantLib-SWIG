
/*
 Copyright (C) 2006, 2007 StatPro Italia srl

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

#ifndef quantlib_callability_i
#define quantlib_callability_i

%include date.i
%include vectors.i
%include common.i

%{
using QuantLib::Callability;
using QuantLib::SoftCallability;
typedef Callability::Price CallabilityPrice;
using QuantLib::CallabilitySchedule;
%}

class CallabilityPrice {
  public:
    enum Type { Dirty, Clean };
    CallabilityPrice(Real amount, Type type);
    Real amount() const;
    Type type() const;
};

%shared_ptr(Callability)

class Callability {
  public:
    enum Type { Call, Put };
    Callability(const CallabilityPrice& price,
                Type type,
                const Date& date);
    const CallabilityPrice& price() const;
    Type type() const;
    Date date() const;
};

%shared_ptr(SoftCallability)

class SoftCallability : public Callability {
  public:
    SoftCallability(const CallabilityPrice& price,
                    const Date& date,
                    Real trigger);
};


#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<Callability> )
#endif
namespace std {
    %template(CallabilitySchedule) vector<ext::shared_ptr<Callability> >;
}

#endif
