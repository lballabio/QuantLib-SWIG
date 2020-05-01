/*
  Copyright (C) 2020 Gorazd Brumen

  Kirk Spread Option engine SWIG part.
*/

#ifndef quantlib_kirk_spread_option_engine_hpp
#define quantlib_kirk_spread_option_engine_hpp

%include common.i

%{
  using QuantLib::KirkSpreadOptionEngine;
  using QuantLib::BlackProcess;
  using QuantLib::Handle;
  using QuantLib::Quote;
%}


// Kirk spread option engine

%shared_ptr(KirkSpreadOptionEngine);
//class KirkSpreadOptionEngine : public SpreadOption::engine {
class KirkSpreadOptionEngine : public PricingEngine {
public:
  KirkSpreadOptionEngine(
                         const boost::shared_ptr<BlackProcess>& process1,
                         const boost::shared_ptr<BlackProcess>& process2,
                         const Handle<Quote>& correlation);
  void calculate() const;
};

#endif
