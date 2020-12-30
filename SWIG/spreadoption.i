/*
 Copyright (C) 2020 Gorazd Brumen

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

#ifndef quantlib_spread_options_i
#define quantlib_spread_options_i

%include options.i

%{
using QuantLib::SpreadOption;
using QuantLib::KirkSpreadOptionEngine;
%}

%shared_ptr(SpreadOption);
class SpreadOption : public MultiAssetOption {
public:
  SpreadOption(const ext::shared_ptr<PlainVanillaPayoff>& payoff,
               const ext::shared_ptr<Exercise>& exercise);
};

%shared_ptr(KirkSpreadOptionEngine);
class KirkSpreadOptionEngine : public PricingEngine {
public:
  KirkSpreadOptionEngine(const ext::shared_ptr<BlackProcess>& process1,
                         const ext::shared_ptr<BlackProcess>& process2,
                         const Handle<Quote>& correlation);
};

#endif
