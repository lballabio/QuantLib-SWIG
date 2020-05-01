/*
  Copyright (C) 2020 Gorazd Brumen

  Spread Option SWIG part.
 */


#ifndef quantlib_spread_option_i
#define quantlib_spread_option_i

%include common.i
/* %include types.i */
/* %include interestrate.i */
/* %include date.i */
/* %include calendars.i */
/* %include daycounters.i */
/* %include currencies.i */
/* %include observer.i */
/* %include marketelements.i */
/* %include interpolation.i */


%{
  using QuantLib::SpreadOption;
  using QuantLib::MultiAssetOption;
  using QuantLib::PlainVanillaPayoff;
  using QuantLib::PricingEngine;
%}

//! Spread option on two assets
%shared_ptr(SpreadOption);
class SpreadOption : public MultiAssetOption {
public:
  SpreadOption(const boost::shared_ptr<PlainVanillaPayoff>& payoff,
               const boost::shared_ptr<Exercise>& exercise);
  void setPricingEngine(const boost::shared_ptr<PricingEngine>&);
};

#endif
