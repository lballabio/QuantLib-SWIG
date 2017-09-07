
/*
 Copyright (C) 2017 Wojciech Ślusarski

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

#ifndef quantlib_black_formula_i
#define quantlib_black_formula_i


%{
using QuantLib::blackFormula;
using QuantLib::blackFormulaImpliedStdDev;
using QuantLib::blackFormulaCashItmProbability;
using QuantLib::bachelierBlackFormula;
using QuantLib::bachelierBlackFormulaImpliedVol;
%}


Real blackFormula(Option::Type optionType,
                    Real strike,
                    Real forward,
                    Real stdDev,
                    Real discount = 1.0,
                    Real displacement = 0.0);

Real blackFormulaImpliedStdDev(Option::Type optionType,
                                Real strike,
                                Real forward,
                                Real blackPrice,
                                Real discount = 1.0,
                                Real displacement = 0.0,
                                Real guess = Null<Real>(),
                                Real accuracy = 1.0e-6,
                                Natural maxIterations = 100);


Real blackFormulaCashItmProbability(Option::Type optionType,
                                    Real strike,
                                    Real forward,
                                    Real stdDev,
                                    Real displacement = 0.0);

Real bachelierBlackFormula(Option::Type optionType,
                            Real strike,
                            Real forward,
                            Real stdDev,
                            Real discount = 1.0);

Real bachelierBlackFormulaImpliedVol(Option::Type optionType,
                                Real strike,
                                Real forward,
                                Real tte,
                                Real bachelierPrice,
                                Real discount = 1.0);


%{
using QuantLib::BlackDeltaCalculator;
%}

class BlackDeltaCalculator{
  public:
    BlackDeltaCalculator(
        Option::Type ot,
        DeltaVolQuote::DeltaType dt,
        Real spot,
        DiscountFactor dDiscount,
        DiscountFactor fDiscount,
        Real stDev);

    Real deltaFromStrike(Real strike) const;
    Real strikeFromDelta(Real delta) const;
    Real atmStrike(DeltaVolQuote::AtmType atmT) const;
};


#endif

