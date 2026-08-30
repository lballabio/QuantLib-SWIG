# coding=utf-8-unix
"""
 Copyright (C) 2017 Wojciech Ślusarski
 Copyright (C) 2026 Chirag Desai
 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
"""
import math
import pytest
import QuantLib as ql


@pytest.fixture
def black_formula_params():
    """Fixture for Black formula test parameters"""
    option_type = ql.Option.Call
    spot = 100.0
    strike = 100.0
    risk_free_rate = 0.05
    expiry = 1.0
    forward = spot * math.exp(risk_free_rate * expiry)
    df = math.exp(-risk_free_rate * expiry)
    vol = 0.2 * math.sqrt(expiry)
    displacement = 0.0

    return {
        'option_type': option_type,
        'spot': spot,
        'strike': strike,
        'risk_free_rate': risk_free_rate,
        'expiry': expiry,
        'forward': forward,
        'df': df,
        'vol': vol,
        'displacement': displacement
    }


def test_blackFormula(black_formula_params):
    """Testing blackFormula in a simple Black-Scholes World..."""
    expected = 10.4506
    res = ql.blackFormula(black_formula_params['option_type'],
                         black_formula_params['strike'],
                         black_formula_params['forward'],
                         black_formula_params['vol'],
                         black_formula_params['df'],
                         black_formula_params['displacement'])
    assert res == pytest.approx(expected, abs=1e-4), \
        "Failed to calculate simple Black-Scholes-Merton price rounded to four decimal places."


def test_black_formula_implied_stdev(black_formula_params):
    """Testing implied volatility calculator"""
    expected = 0.2 * math.sqrt(black_formula_params['expiry'])
    black_price = 10.4506
    res = ql.blackFormulaImpliedStdDev(black_formula_params['option_type'],
                                      black_formula_params['strike'],
                                      black_formula_params['forward'],
                                      black_price,
                                      black_formula_params['df'])
    assert res == pytest.approx(expected, abs=1e-4), \
        "Failed to determine Implied Vol rounded to a single vol bps."


@pytest.fixture
def black_delta_params():
    """Fixture for BlackDeltaCalculator test parameters"""
    todaysDate = ql.Date(5, ql.September, 2017)
    ql.Settings.instance().evaluationDate = todaysDate
    spotDate = ql.Date(7, ql.September, 2017)
    domestic_rate = ql.FlatForward(spotDate, 0.017, ql.Actual365Fixed())
    foreign_rate = ql.FlatForward(spotDate, 0.013, ql.Actual365Fixed())

    yield {
        'todaysDate': todaysDate,
        'spotDate': spotDate,
        'domestic_rate': domestic_rate,
        'foreign_rate': foreign_rate
    }

    ql.Settings.instance().evaluationDate = ql.Date()


def test_single_spot_delta(black_delta_params):
    """Test for a single strike for call spot delta 75"""
    volatility = 0.2
    expiry = 2
    spot_price = 3.6
    domDf = black_delta_params['domestic_rate'].discount(expiry)
    forDf = black_delta_params['foreign_rate'].discount(expiry)
    forward = spot_price * forDf / domDf

    spot_delta_level = 0.75
    stDev = volatility * expiry ** 0.5

    inv_norm_dist = ql.InverseCumulativeNormal()
    expected_strike = inv_norm_dist(spot_delta_level / forDf)
    expected_strike *= stDev
    expected_strike -= 0.5 * stDev ** 2
    expected_strike = math.exp(expected_strike) / forward
    expected_strike = 1 / expected_strike

    option_type = ql.Option.Call
    delta_type = ql.DeltaVolQuote.Spot

    black_calculator = ql.BlackDeltaCalculator(option_type,
                                               delta_type,
                                               spot_price,
                                               domDf,
                                               forDf,
                                               stDev)

    strike = black_calculator.strikeFromDelta(spot_delta_level)

    assert strike == pytest.approx(expected_strike, abs=1e-4)


def test_spot_atm_delta_calculator(black_delta_params):
    """Test for 0-delta straddle strike"""
    volatility = 0.2
    expiry = 2
    spot_price = 3.6
    domDf = black_delta_params['domestic_rate'].discount(expiry)
    forDf = black_delta_params['foreign_rate'].discount(expiry)
    forward = spot_price * forDf / domDf
    expected_strike = forward * math.exp(-0.5 * volatility ** 2 * expiry)

    option_type = ql.Option.Call
    delta_type = ql.DeltaVolQuote.AtmDeltaNeutral
    stDev = volatility * expiry ** 0.5

    black_calculator = ql.BlackDeltaCalculator(option_type,
                                               delta_type,
                                               spot_price,
                                               domDf,
                                               forDf,
                                               stDev)

    strike = black_calculator.atmStrike(ql.DeltaVolQuote.AtmDeltaNeutral)

    assert strike == pytest.approx(expected_strike, abs=1e-4)


if __name__ == '__main__':
    print("testing QuantLib", ql.__version__)
    pytest.main([__file__, '-v'])
