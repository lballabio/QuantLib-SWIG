
/*
 Copyright (C) 2003 StatPro Italia srl
 Copyright (C) 2018, 2019 Matthias Lungwitz

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

#ifndef quantlib_payoffs_i
#define quantlib_payoffs_i

%include options.i

// payoffs

%{
using QuantLib::PlainVanillaPayoff;
using QuantLib::PercentageStrikePayoff;
using QuantLib::CashOrNothingPayoff;
using QuantLib::AssetOrNothingPayoff;
using QuantLib::SuperSharePayoff;
using QuantLib::GapPayoff;
using QuantLib::VanillaForwardPayoff;
%}


%shared_ptr(PlainVanillaPayoff)
class PlainVanillaPayoff : public StrikedTypePayoff {
  public:
    PlainVanillaPayoff(Option::Type type,
                          Real strike);
};

%inline %{
    const ext::shared_ptr<PlainVanillaPayoff> as_plain_vanilla_payoff(
                           const ext::shared_ptr<Payoff>& payoff) {
        return ext::dynamic_pointer_cast<PlainVanillaPayoff>(payoff);
    }
%}

%shared_ptr(PercentageStrikePayoff)
class PercentageStrikePayoff : public StrikedTypePayoff {
  public:
    PercentageStrikePayoff(Option::Type type,
                              Real moneyness);
};

%shared_ptr(CashOrNothingPayoff)
class CashOrNothingPayoff : public StrikedTypePayoff {
  public:
    CashOrNothingPayoff(Option::Type type,
                           Real strike,
                           Real payoff);
};

%shared_ptr(AssetOrNothingPayoff)
class AssetOrNothingPayoff : public StrikedTypePayoff {
  public:
    AssetOrNothingPayoff(Option::Type type,
                            Real strike);
};

%shared_ptr(SuperSharePayoff)
class SuperSharePayoff : public StrikedTypePayoff {
  public:
    SuperSharePayoff(Option::Type type,
                        Real strike,
                        Real increment);
};

%shared_ptr(GapPayoff)
class GapPayoff : public StrikedTypePayoff {
  public:
    GapPayoff(Option::Type type,
                        Real strike,
                        Real strikePayoff);
};

%shared_ptr(VanillaForwardPayoff)
class VanillaForwardPayoff : public StrikedTypePayoff {
  public:
    VanillaForwardPayoff(Option::Type type, Real strike);
};


#endif
