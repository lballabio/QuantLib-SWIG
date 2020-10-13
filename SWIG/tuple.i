/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007 StatPro Italia srl
 Copyright (C) 2019 Matthias Lungwitz

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


#ifndef quantlib_tuple_i
#define quantlib_tuple_i

%include common.i

namespace ext {
  template <typename T1=void, typename T2=void, typename T3=void>
  struct tuple;

  template <>
  struct tuple<void,void,void> {
  };

  template <typename T1>
  struct tuple<T1, void, void> {
    tuple(T1);
    %extend {
      T1 first() const {
        return ext::get<0>(*$self);
      }
    }
  };

  template <typename T1, typename T2>
  struct tuple <T1, T2, void> {
    tuple(T1,T2);
    %extend {
      T1 first() const {
        return ext::get<0>(*$self);
      }
      T2 second() const {
        return ext::get<1>(*$self);
      }
    }
  };

  template <typename T1, typename T2, typename T3>
  struct tuple <T1,T2,T3> {
    tuple(T1,T2,T3);
    %extend {
      T1 first() const {
        return ext::get<0>(*$self);
      }
      T2 second() const {
        return ext::get<1>(*$self);
      }
      T3 third() const {
        return ext::get<2>(*$self);
      }
    }
  };
}

#endif
