
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005 StatPro Italia srl

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

#ifndef quantlib_money_i
#define quantlib_money_i

%include currencies.i

%{
using QuantLib::Money;
%}

class Money {
    #if defined(SWIGRUBY)
    %rename("conversionType=") setConversionType;
    %rename("baseCurrency=")   setBaseCurrency;
    #elif defined(SWIGJAVA)
    %rename("compare") __cmp__;
    #endif
  public:
    Money(const Currency& currency, Decimal value);
    Money(Decimal value, const Currency& currency);
    const Currency& currency() const;
    Decimal value() const;
    Money rounded() const;

    #if defined(SWIGPYTHON) || defined(SWIGRUBY) || defined(SWIGJAVA)
    Money operator+() const;
    Money operator-() const;
    %extend {
        Money operator+(const Money& m) { return *self+m; }
        Money operator-(const Money& m) { return *self-m; }
        Money operator*(Decimal x) { return *self*x; }
        Money operator/(Decimal x) { return *self/x; }
        Decimal operator/(const Money& m) { return *self/m; }
        #if defined(SWIGPYTHON)
        Money __rmul__(Decimal x) { return *self*x; }
        bool __lt__(const Money& other) {
            return *self < other;
        }
        bool __gt__(const Money& other) {
            return other < *self;
        }
        bool __le__(const Money& other) {
            return !(other < *self);
        }
        bool __ge__(const Money& other) {
            return !(*self < other);
        }
        #endif
        int __cmp__(const Money& other) {
            if (*self < other)
                return -1;
            else if (*self == other)
                return 0;
            else
                return 1;
        }
        std::string __str__() {
            std::ostringstream out;
            out << *self;
            return out.str();
        }
    }
    #endif

    enum ConversionType { NoConversion,
                          BaseCurrencyConversion,
                          AutomatedConversion };
    %extend {
        static void setConversionType(ConversionType type) {
            Money::conversionType = type;
        }
        static void setBaseCurrency(const Currency& c) {
            Money::baseCurrency = c;
        }
    }
};


#endif
