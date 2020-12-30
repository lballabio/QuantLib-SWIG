/*
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 StatPro Italia srl

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

#ifndef quantlib_lookback_options_i
#define quantlib_lookback_options_i

%include options.i

%{
using QuantLib::ContinuousFloatingLookbackOption;
using QuantLib::ContinuousFixedLookbackOption;
using QuantLib::ContinuousPartialFloatingLookbackOption;
using QuantLib::ContinuousPartialFixedLookbackOption;
%}

%shared_ptr(ContinuousFloatingLookbackOption)
class ContinuousFloatingLookbackOption : public OneAssetOption {
  public:
    ContinuousFloatingLookbackOption(
                          Real currentMinmax,
                          const ext::shared_ptr<TypePayoff>& payoff,
                          const ext::shared_ptr<Exercise>& exercise);
};

%shared_ptr(ContinuousFixedLookbackOption)
class ContinuousFixedLookbackOption : public OneAssetOption {
  public:
    ContinuousFixedLookbackOption(
                          Real currentMinmax,
                          const ext::shared_ptr<StrikedTypePayoff>& payoff,
                          const ext::shared_ptr<Exercise>& exercise);
};

%shared_ptr(ContinuousPartialFloatingLookbackOption)
class ContinuousPartialFloatingLookbackOption : public ContinuousFloatingLookbackOption {
  public:
    ContinuousPartialFloatingLookbackOption(
                          Real currentMinmax,
                          Real lambda,
                          Date lookbackPeriodEnd,
                          const ext::shared_ptr<TypePayoff>& payoff,
                          const ext::shared_ptr<Exercise>& exercise);
};

%shared_ptr(ContinuousPartialFixedLookbackOption)
class ContinuousPartialFixedLookbackOption : public ContinuousFixedLookbackOption {
  public:
    ContinuousPartialFixedLookbackOption(
                          Date lookbackPeriodStart,
                          const ext::shared_ptr<StrikedTypePayoff>& payoff,
                          const ext::shared_ptr<Exercise>& exercise);
};


%{
using QuantLib::AnalyticContinuousFloatingLookbackEngine;
using QuantLib::AnalyticContinuousFixedLookbackEngine;
using QuantLib::AnalyticContinuousPartialFloatingLookbackEngine;
using QuantLib::AnalyticContinuousPartialFixedLookbackEngine;
%}

%shared_ptr(AnalyticContinuousFloatingLookbackEngine)
class AnalyticContinuousFloatingLookbackEngine : public PricingEngine {
  public:
        AnalyticContinuousFloatingLookbackEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%shared_ptr(AnalyticContinuousFixedLookbackEngine)
class AnalyticContinuousFixedLookbackEngine : public PricingEngine {
  public:
        AnalyticContinuousFixedLookbackEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%shared_ptr(AnalyticContinuousPartialFloatingLookbackEngine)
class AnalyticContinuousPartialFloatingLookbackEngine : public PricingEngine {
  public:
        AnalyticContinuousPartialFloatingLookbackEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};

%shared_ptr(AnalyticContinuousPartialFixedLookbackEngine)
class AnalyticContinuousPartialFixedLookbackEngine : public PricingEngine {
  public:
        AnalyticContinuousPartialFixedLookbackEngine(
            const ext::shared_ptr<GeneralizedBlackScholesProcess>& process);
};


#endif
