
/*
 Copyright (C) 2006, 2007 StatPro Italia srl
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

#ifndef quantlib_convertible_bonds_i
#define quantlib_convertible_bonds_i

%include bonds.i
%include callability.i
%include dividends.i
%include options.i

%{
using QuantLib::ConvertibleZeroCouponBond;
using QuantLib::ConvertibleFixedCouponBond;
using QuantLib::ConvertibleFloatingRateBond;
using QuantLib::BinomialConvertibleEngine;
%}

%shared_ptr(ConvertibleZeroCouponBond)
class ConvertibleZeroCouponBond : public Bond {
  public:
    ConvertibleZeroCouponBond(
          const boost::shared_ptr<Exercise>& exercise,
          Real conversionRatio,
          const std::vector<boost::shared_ptr<Dividend> >& dividends,
          const std::vector<boost::shared_ptr<Callability> >& callability,
          const Handle<Quote>& creditSpread,
          const Date& issueDate,
          Integer settlementDays,
          const DayCounter& dayCounter,
          const Schedule& schedule,
          Real redemption = 100.0);
};


%shared_ptr(ConvertibleFixedCouponBond)
class ConvertibleFixedCouponBond : public Bond {
  public:
    ConvertibleFixedCouponBond(
          const boost::shared_ptr<Exercise>& exercise,
          Real conversionRatio,
          const std::vector<boost::shared_ptr<Dividend> >& dividends,
          const std::vector<boost::shared_ptr<Callability> >& callability,
          const Handle<Quote>& creditSpread,
          const Date& issueDate,
          Integer settlementDays,
          const std::vector<Rate>& coupons,
          const DayCounter& dayCounter,
          const Schedule& schedule,
          Real redemption = 100.0,
          const Period& exCouponPeriod = Period(),
          const Calendar& exCouponCalendar = Calendar(),
          const BusinessDayConvention exCouponConvention = Unadjusted,
          bool exCouponEndOfMonth = false);
};


%shared_ptr(ConvertibleFloatingRateBond)
class ConvertibleFloatingRateBond : public Bond {
  public:
    ConvertibleFloatingRateBond(
          const boost::shared_ptr<Exercise>& exercise,
          Real conversionRatio,
          const std::vector<boost::shared_ptr<Dividend> >& dividends,
          const std::vector<boost::shared_ptr<Callability> >& callability,
          const Handle<Quote>& creditSpread,
          const Date& issueDate,
          Integer settlementDays,
          const boost::shared_ptr<IborIndex>& index,
          Integer fixingDays,
          const std::vector<Spread>& spreads,
          const DayCounter& dayCounter,
          const Schedule& schedule,
          Real redemption = 100.0,
          const Period& exCouponPeriod = Period(),
          const Calendar& exCouponCalendar = Calendar(),
          const BusinessDayConvention exCouponConvention = Unadjusted,
          bool exCouponEndOfMonth = false);
};


%shared_ptr(BinomialConvertibleEngine<CoxRossRubinstein>)
%shared_ptr(BinomialConvertibleEngine<JarrowRudd>)
%shared_ptr(BinomialConvertibleEngine<AdditiveEQPBinomialTree>)
%shared_ptr(BinomialConvertibleEngine<Trigeorgis>)
%shared_ptr(BinomialConvertibleEngine<Tian>)
%shared_ptr(BinomialConvertibleEngine<LeisenReimer>)
%shared_ptr(BinomialConvertibleEngine<Joshi4>)

template <class T>
class BinomialConvertibleEngine : public PricingEngine {
  public:
    BinomialConvertibleEngine(const boost::shared_ptr<GeneralizedBlackScholesProcess>&,
                              Size steps);
};

%template(BinomialCRRConvertibleEngine) BinomialConvertibleEngine<CoxRossRubinstein>;
%template(BinomialJRConvertibleEngine) BinomialConvertibleEngine<JarrowRudd>;
%template(BinomialEQPConvertibleEngine) BinomialConvertibleEngine<AdditiveEQPBinomialTree>;
%template(BinomialTrigeorgisConvertibleEngine) BinomialConvertibleEngine<Trigeorgis>;
%template(BinomialTianConvertibleEngine) BinomialConvertibleEngine<Tian>;
%template(BinomialLRConvertibleEngine) BinomialConvertibleEngine<LeisenReimer>;
%template(BinomialJ4ConvertibleEngine) BinomialConvertibleEngine<Joshi4>;

#if defined(SWIGPYTHON)
%pythoncode %{
    def BinomialConvertibleEngine(process, type, steps):
        type = type.lower()
        if type == "crr" or type == "coxrossrubinstein":
            cls = BinomialCRRConvertibleEngine
        elif type == "jr" or type == "jarrowrudd":
            cls = BinomialJRConvertibleEngine
        elif type == "eqp":
            cls = BinomialEQPConvertibleEngine
        elif type == "trigeorgis":
            cls = BinomialTrigeorgisConvertibleEngine
        elif type == "tian":
            cls = BinomialTianConvertibleEngine
        elif type == "lr" or type == "leisenreimer":
            cls = BinomialLRConvertibleEngine
        elif type == "j4" or type == "joshi4":
            cls = BinomialJ4ConvertibleEngine
        else:
            raise RuntimeError("unknown binomial engine type: %s" % type);
        return cls(process, steps)
%}
#endif

#endif
