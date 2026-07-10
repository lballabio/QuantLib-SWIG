"""
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

'''
TODO 
The Heston MC tests (testMCDiscreteGeometricAveragePriceHeston, testMCDiscreteArithmeticAveragePriceHeston, testDiscreteGeometricAveragePriceHestonPastFixings, testChoiAsianEngineVsMC, testChoiAsianEngineSpecialCases, testPastFixingsModelDependency) are 
omitted because the corresponding Python bindings 
(MakeMCDiscreteGeometricAPHestonEngine, ChoiAsianEngine, Heston MC arithmetic engine are not exposed

'''
import math
import pytest
import QuantLib as ql

# ---------------------------------------------------------------------------
# Helpers mirroring utilities.hpp flatRate / flatVol
# ---------------------------------------------------------------------------

def flat_rate(forward, dc, today=None):
    if today is not None:
        return ql.FlatForward(today, forward, dc)
    return ql.FlatForward(0, ql.NullCalendar(), ql.QuoteHandle(ql.SimpleQuote(forward)), dc)


def flat_vol(vol, dc, today=None):
    if today is not None:
        return ql.BlackConstantVol(today, ql.NullCalendar(), vol, dc)
    return ql.BlackConstantVol(0, ql.NullCalendar(), ql.QuoteHandle(ql.SimpleQuote(vol)), dc)


def make_bsm_process(spot_val, q_val, r_val, vol_val, dc, today):
    spot = ql.SimpleQuote(spot_val)
    qTS = ql.YieldTermStructureHandle(flat_rate(q_val, dc, today))
    rTS = ql.YieldTermStructureHandle(flat_rate(r_val, dc, today))
    volTS = ql.BlackVolTermStructureHandle(flat_vol(vol_val, dc, today))
    return ql.BlackScholesMertonProcess(
        ql.QuoteHandle(spot), qTS, rTS, volTS
    )

# ---------------------------------------------------------------------------
# Test 1 – Analytic continuous geometric average-price
# ---------------------------------------------------------------------------
def test_analytic_continuous_geometric_average_price():
    """Analytic continuous geometric average-price Asians 
     // data from "Option Pricing Formulas 2nd Edition", Haug, page 183"""
    dc = ql.Actual360()
    today = ql.Settings.instance().evaluationDate

    spot = ql.SimpleQuote(80.0)
    qRate = ql.SimpleQuote(-0.03)
    qTS = ql.YieldTermStructureHandle(ql.FlatForward(today, ql.QuoteHandle(qRate), dc))
    rRate = ql.SimpleQuote(0.05)
    rTS = ql.YieldTermStructureHandle(ql.FlatForward(today, ql.QuoteHandle(rRate), dc))
    vol = ql.SimpleQuote(0.20)
    volTS = ql.BlackVolTermStructureHandle(ql.BlackConstantVol(today, ql.NullCalendar(), ql.QuoteHandle(vol), dc))

    process = ql.BlackScholesMertonProcess(
        ql.QuoteHandle(spot), qTS, rTS, volTS
    )

    engine = ql.AnalyticContinuousGeometricAveragePriceAsianEngine(process)

    payoff = ql.PlainVanillaPayoff(ql.Option.Put, strike =  85.0)
    exercise = ql.EuropeanExercise(today + 90)
    average_type = ql.Average.Geometric
    
    option = ql.ContinuousAveragingAsianOption( average_type, payoff, exercise)
    option.setPricingEngine(engine)

    calculated = option.NPV()
    expected = 4.6922
    tolerance = 1.0e-4
    assert abs(calculated - expected) <= tolerance, (
        f"Continuous geometric Asian: expected {expected}, got {calculated}"
    )

    # Approximate continuous version with discrete version
    running_accumulator = 1.0
    past_fixings = 0
    fixing_dates = [today + i for i in range(exercise.lastDate() - today + 1)]

    engine2 = ql.AnalyticDiscreteGeometricAveragePriceAsianEngine(process)
    option2 = ql.DiscreteAveragingAsianOption(
        ql.Average.Geometric, running_accumulator, past_fixings, fixing_dates, payoff, exercise
    )
    option2.setPricingEngine(engine2)

    calculated = option2.NPV()
    tolerance = 3.0e-3
    assert abs(calculated - expected) <= tolerance, (
        f"Discrete approx of continuous geometric Asian: expected {expected}, got {calculated}"
    )

# ---------------------------------------------------------------------------
# Test 2 – Analytic continuous geometric average-price Greeks
# ---------------------------------------------------------------------------
def test_analytic_continuous_geometric_average_price_greeks():
    """Analytic continuous geometric average-price Asian greeks."""
    tol = {
        "delta": 1.0e-5,
        "gamma": 1.0e-5,
        "theta": 1.0e-5,
        "rho": 1.0e-5,
        "divRho": 1.0e-5,
        "vega": 1.0e-5,
    }

    types = [ql.Option.Call, ql.Option.Put]
    underlyings = [100.0]
    strikes = [90.0, 100.0, 110.0]
    q_rates = [0.04, 0.05, 0.06]
    r_rates = [0.01, 0.05, 0.15]
    lengths = [1, 2]
    vols = [0.11, 0.50, 1.20]

    dc = ql.Actual360()
    today = ql.Settings.instance().evaluationDate

    spot = ql.SimpleQuote(0.0)
    qRate = ql.SimpleQuote(0.0)
    qTS = ql.YieldTermStructureHandle(
        ql.FlatForward(0, ql.NullCalendar(), ql.QuoteHandle(qRate), dc)
    )
    rRate = ql.SimpleQuote(0.0)
    rTS = ql.YieldTermStructureHandle(
        ql.FlatForward(0, ql.NullCalendar(), ql.QuoteHandle(rRate), dc)
    )
    vol_q = ql.SimpleQuote(0.0)
    volTS = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(0, ql.NullCalendar(), ql.QuoteHandle(vol_q), dc)
    )

    process = ql.BlackScholesMertonProcess(ql.QuoteHandle(spot), qTS, rTS, volTS)

    for opt_type in types:
        for strike in strikes:
            for length in lengths:
                maturity = ql.EuropeanExercise(today + ql.Period(length, ql.Years))
                payoff = ql.PlainVanillaPayoff(opt_type, strike)
                engine = ql.AnalyticContinuousGeometricAveragePriceAsianEngine(process)
                option = ql.ContinuousAveragingAsianOption(ql.Average.Geometric, payoff, maturity)
                option.setPricingEngine(engine)

                for u in underlyings:
                    for q in q_rates:
                        for r in r_rates:
                            for v in vols:
                                spot.setValue(u)
                                qRate.setValue(q)
                                rRate.setValue(r)
                                vol_q.setValue(v)

                                value = option.NPV()
                                if value <= u * 1.0e-5:
                                    continue

                                calc = {
                                    "delta": option.delta(),
                                    "gamma": option.gamma(),
                                    "theta": option.theta(),
                                    "rho": option.rho(),
                                    "divRho": option.dividendRho(),
                                    "vega": option.vega(),
                                }
                                # perturb spot and get delta and gamma
                                du = u * 1.0e-4
                                spot.setValue(u + du)
                                vp = option.NPV(); dp = option.delta()
                                spot.setValue(u - du)
                                vm = option.NPV(); dm = option.delta()
                                spot.setValue(u)
                                exp = {
                                    "delta": (vp - vm) / (2 * du),
                                    "gamma": (dp - dm) / (2 * du),
                                }
                                # perturb rates and get rho and dividend rho
                                dr = r * 1.0e-4
                                rRate.setValue(r + dr); vp = option.NPV()
                                rRate.setValue(r - dr); vm = option.NPV()
                                rRate.setValue(r)
                                exp["rho"] = (vp - vm) / (2 * dr)

                                dq = q * 1.0e-4
                                qRate.setValue(q + dq); vp = option.NPV()
                                qRate.setValue(q - dq); vm = option.NPV()
                                qRate.setValue(q)
                                exp["divRho"] = (vp - vm) / (2 * dq)
                                
                                # perturb volatility and get vega
                                dv = v * 1.0e-4
                                vol_q.setValue(v + dv); vp = option.NPV()
                                vol_q.setValue(v - dv); vm = option.NPV()
                                vol_q.setValue(v)
                                exp["vega"] = (vp - vm) / (2 * dv)

                                # // perturb date and get theta
                                dT = dc.yearFraction(today - 1, today + 1)
                                ql.Settings.instance().evaluationDate = today - 1
                                vm = option.NPV()
                                ql.Settings.instance().evaluationDate = today + 1
                                vp = option.NPV()
                                ql.Settings.instance().evaluationDate = today
                                exp["theta"] = (vp - vm) / dT

                                for greek, calcl in calc.items():
                                    expct = exp[greek]
                                    denom = u if u != 0 else 1.0
                                    error = abs(expct - calcl) / denom
                                    assert error <= tol[greek], (
                                        f"Continuous geom Asian greek '{greek}': "
                                        f"expected {expct}, got {calcl}, error {error}"
                                    )
# ---------------------------------------------------------------------------
# Test 3 – Analytic discrete geometric average-price
# ---------------------------------------------------------------------------
def test_analytic_discrete_geometric_average_price():
    """Analytic discrete geometric average-price Asians (Clewlow & Strickland p.118-123)."""
    dc = ql.Actual360()
    today = ql.Settings.instance().evaluationDate

    spot = ql.SimpleQuote(100.0)
    qTS = ql.YieldTermStructureHandle(ql.FlatForward(today, 0.03, dc))
    rTS = ql.YieldTermStructureHandle(ql.FlatForward(today, 0.06, dc))
    volTS = ql.BlackVolTermStructureHandle(ql.BlackConstantVol(today, ql.NullCalendar(), 0.20, dc))

    process = ql.BlackScholesMertonProcess(ql.QuoteHandle(spot), qTS, rTS, volTS)
    engine = ql.AnalyticDiscreteGeometricAveragePriceAsianEngine(process)

    future_fixings = 10 #assume 10 equally spaced fixings over the next year
    dt = round(360.0 / future_fixings)
    fixing_dates = [today + int(i * dt) for i in range(1, future_fixings + 1)]
    exercise = ql.EuropeanExercise(today + 360)
    payoff = ql.PlainVanillaPayoff(ql.Option.Call, strike = 100.0)

    option = ql.DiscreteAveragingAsianOption(
        ql.Average.Geometric, 1.0, 0, fixing_dates, payoff, exercise
    )
    option.setPricingEngine(engine)

    calculated = option.NPV()
    expected = 5.3425606635
    tolerance = 1e-10
    assert abs(calculated - expected) <= tolerance, (
        f"Discrete geom avg price: expected {expected}, got {calculated}"
    )

# ---------------------------------------------------------------------------
# Test 4 – Analytic discrete geometric average-strike
 # data from "Implementing Derivatives Model",
 #Clewlow, Strickland, p.118-123
# ---------------------------------------------------------------------------
def test_analytic_discrete_geometric_average_strike():
    """Analytic discrete geometric average-strike Asians."""
    dc = ql.Actual360()
    today = ql.Settings.instance().evaluationDate

    spot = ql.SimpleQuote(100.0)
    qTS = ql.YieldTermStructureHandle(ql.FlatForward(today, 0.03, dc))
    rTS = ql.YieldTermStructureHandle(ql.FlatForward(today, 0.06, dc))
    volTS = ql.BlackVolTermStructureHandle(ql.BlackConstantVol(today, ql.NullCalendar(), 0.20, dc))

    process = ql.BlackScholesMertonProcess(ql.QuoteHandle(spot), qTS, rTS, volTS)
    engine = ql.AnalyticDiscreteGeometricAverageStrikeAsianEngine(process)

    future_fixings = 10
    dt = round(360.0 / future_fixings)
    fixing_dates = [today + int(i * dt) for i in range(1, future_fixings + 1)]
    exercise = ql.EuropeanExercise(today + 360)
    payoff = ql.PlainVanillaPayoff(ql.Option.Call, 100.0)

    option = ql.DiscreteAveragingAsianOption(
        ql.Average.Geometric, 1.0, 0, fixing_dates, payoff, exercise
    )
    option.setPricingEngine(engine)

    calculated = option.NPV()
    expected = 4.97109
    tolerance = 1e-5
    assert abs(calculated - expected) <= tolerance, (
        f"Discrete geom avg strike: expected {expected}, got {calculated}"
    )

# ---------------------------------------------------------------------------
# Test 6 – Analytic discrete geometric average-price Greeks
# ---------------------------------------------------------------------------
def test_analytic_discrete_geometric_average_price_greeks():
    """Discrete-averaging geometric Asian greeks."""
    tol = {
        "delta": 1.0e-5,
        "gamma": 1.0e-5,
        "theta": 1.0e-5,
        "rho": 1.0e-5,
        "divRho": 1.0e-5,
        "vega": 1.0e-5,
    }

    types = [ql.Option.Call, ql.Option.Put]
    underlyings = [100.0]
    strikes = [90.0, 100.0, 110.0]
    q_rates = [0.04, 0.05, 0.06]
    r_rates = [0.01, 0.05, 0.15]
    lengths = [1, 2]
    vols = [0.11, 0.50, 1.20]

    dc = ql.Actual360()
    today = ql.Settings.instance().evaluationDate

    spot = ql.SimpleQuote(0.0)
    qRate = ql.SimpleQuote(0.0)
    qTS = ql.YieldTermStructureHandle(
        ql.FlatForward(0, ql.NullCalendar(), ql.QuoteHandle(qRate), dc)
    )
    rRate = ql.SimpleQuote(0.0)
    rTS = ql.YieldTermStructureHandle(
        ql.FlatForward(0, ql.NullCalendar(), ql.QuoteHandle(rRate), dc)
    )
    vol_q = ql.SimpleQuote(0.0)
    volTS = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(0, ql.NullCalendar(), ql.QuoteHandle(vol_q), dc)
    )
    process = ql.BlackScholesMertonProcess(ql.QuoteHandle(spot), qTS, rTS, volTS)

    for opt_type in types:
        for strike in strikes:
            for length in lengths:
                maturity = ql.EuropeanExercise(today + ql.Period(length, ql.Years))
                payoff = ql.PlainVanillaPayoff(opt_type, strike)
                running_average = 120.0
                past_fixings = 1
                fixing_dates = []
                d = today + ql.Period(3, ql.Months)
                while d <= maturity.lastDate():
                    fixing_dates.append(d)
                    d += ql.Period(3, ql.Months)

                engine = ql.AnalyticDiscreteGeometricAveragePriceAsianEngine(process)
                option = ql.DiscreteAveragingAsianOption(
                    ql.Average.Geometric, running_average, past_fixings,
                    fixing_dates, payoff, maturity
                )
                option.setPricingEngine(engine)

                for u in underlyings:
                    for q in q_rates:
                        for r in r_rates:
                            for v in vols:
                                spot.setValue(u)
                                qRate.setValue(q)
                                rRate.setValue(r)
                                vol_q.setValue(v)

                                value = option.NPV()
                                if value <= u * 1.0e-5:
                                    continue

                                calc = {
                                    "delta": option.delta(),
                                    "gamma": option.gamma(),
                                    "theta": option.theta(),
                                    "rho": option.rho(),
                                    "divRho": option.dividendRho(),
                                    "vega": option.vega(),
                                }

                                du = u * 1.0e-4
                                spot.setValue(u + du)
                                vp = option.NPV(); dp = option.delta()
                                spot.setValue(u - du)
                                vm = option.NPV(); dm = option.delta()
                                spot.setValue(u)
                                exp = {
                                    "delta": (vp - vm) / (2 * du),
                                    "gamma": (dp - dm) / (2 * du),
                                }

                                dr = r * 1.0e-4
                                rRate.setValue(r + dr); vp = option.NPV()
                                rRate.setValue(r - dr); vm = option.NPV()
                                rRate.setValue(r)
                                exp["rho"] = (vp - vm) / (2 * dr)

                                dq = q * 1.0e-4
                                qRate.setValue(q + dq); vp = option.NPV()
                                qRate.setValue(q - dq); vm = option.NPV()
                                qRate.setValue(q)
                                exp["divRho"] = (vp - vm) / (2 * dq)

                                dv = v * 1.0e-4
                                vol_q.setValue(v + dv); vp = option.NPV()
                                vol_q.setValue(v - dv); vm = option.NPV()
                                vol_q.setValue(v)
                                exp["vega"] = (vp - vm) / (2 * dv)

                                dT = dc.yearFraction(today - 1, today + 1)
                                ql.Settings.instance().evaluationDate = today - 1
                                vm = option.NPV()
                                ql.Settings.instance().evaluationDate = today + 1
                                vp = option.NPV()
                                ql.Settings.instance().evaluationDate = today
                                exp["theta"] = (vp - vm) / dT

                                for greek, calcl in calc.items():
                                    expct = exp[greek]
                                    denom = u if u != 0 else 1.0
                                    error = abs(expct - calcl) / denom
                                    assert error <= tol[greek], (
                                        f"Discrete geom Asian greek '{greek}': "
                                        f"expected {expct}, got {calcl}, error {error}"
                                    )


# ---------------------------------------------------------------------------
# Test 7 – MC discrete arithmetic average-price
#TODO --- expose MakeMCDiscreteGeometricAPEngine
# ---------------------------------------------------------------------------

DISCRETE_ARITH_AP_CASES = [
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 0.0,       11.0/12.0, 2,    0.13, True,  1.3942835683),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 0.0,       11.0/12.0, 4,    0.13, True,  1.5852442983),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 0.0,       11.0/12.0, 8,    0.13, True,  1.66970673),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 0.0,       11.0/12.0, 12,   0.13, True,  1.6980019214),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 0.0,       11.0/12.0, 26,   0.13, True,  1.7255070456),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 0.0,       11.0/12.0, 52,   0.13, True,  1.7401553533),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 0.0,       11.0/12.0, 100,  0.13, True,  1.7478303712),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 0.0,       11.0/12.0, 250,  0.13, True,  1.7490291943),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 0.0,       11.0/12.0, 500,  0.13, True,  1.7515113291),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 0.0,       11.0/12.0, 1000, 0.13, True,  1.7537344885),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0,  11.0/12.0, 2,    0.13, True,  1.8496053697),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0,  11.0/12.0, 4,    0.13, True,  2.0111495205),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0,  11.0/12.0, 8,    0.13, True,  2.0852138818),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0,  11.0/12.0, 12,   0.13, True,  2.1105094397),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0,  11.0/12.0, 26,   0.13, True,  2.1346526695),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0,  11.0/12.0, 52,   0.13, True,  2.147489651),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0,  11.0/12.0, 100,  0.13, True,  2.154728109),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0,  11.0/12.0, 250,  0.13, True,  2.1564276565),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0,  11.0/12.0, 500,  0.13, True,  2.1594238588),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 1.0/12.0,  11.0/12.0, 1000, 0.13, True,  2.1595367326),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0,  11.0/12.0, 2,    0.13, True,  2.63315092584),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0,  11.0/12.0, 4,    0.13, True,  2.76723962361),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0,  11.0/12.0, 8,    0.13, True,  2.83124836881),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0,  11.0/12.0, 12,   0.13, True,  2.84290301412),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0,  11.0/12.0, 26,   0.13, True,  2.88179560417),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0,  11.0/12.0, 52,   0.13, True,  2.88447044543),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0,  11.0/12.0, 100,  0.13, True,  2.89985329603),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0,  11.0/12.0, 250,  0.13, True,  2.90047296063),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0,  11.0/12.0, 500,  0.13, True,  2.89813412160),
    (ql.Option.Put, 90.0, 87.0, 0.06, 0.025, 3.0/12.0,  11.0/12.0, 1000, 0.13, True,  2.89703362437),
]

def _time_to_days(t, basis=360):
    return int(round(t * basis))

# ---------------------------------------------------------------------------
# Helper for Levy engine data
# ---------------------------------------------------------------------------
LEVY_CASES = [
    (ql.Option.Call, 6.80,  6.80,  6.90,  0.09, 0.07, 0.14, 180, 0,   0.0944),
    (ql.Option.Put,  6.80,  6.80,  6.90,  0.09, 0.07, 0.14, 180, 0,   0.2237),
    (ql.Option.Call, 100.0, 100.0, 95.0,  0.05, 0.1,  0.15, 270, 0,   7.0544),
    (ql.Option.Call, 100.0, 100.0, 95.0,  0.05, 0.1,  0.15, 270, 90,  5.6731),
    (ql.Option.Call, 100.0, 100.0, 95.0,  0.05, 0.1,  0.15, 270, 180, 5.0806),
    (ql.Option.Call, 100.0, 100.0, 95.0,  0.05, 0.1,  0.35, 270, 0,   10.1213),
    (ql.Option.Call, 100.0, 100.0, 95.0,  0.05, 0.1,  0.35, 270, 90,  6.9705),
    (ql.Option.Call, 100.0, 100.0, 95.0,  0.05, 0.1,  0.35, 270, 180, 5.1411),
    (ql.Option.Call, 100.0, 100.0, 100.0, 0.05, 0.1,  0.15, 270, 0,   3.7845),
    (ql.Option.Call, 100.0, 100.0, 100.0, 0.05, 0.1,  0.15, 270, 90,  1.9964),
    (ql.Option.Call, 100.0, 100.0, 100.0, 0.05, 0.1,  0.15, 270, 180, 0.6722),
    (ql.Option.Call, 100.0, 100.0, 100.0, 0.05, 0.1,  0.35, 270, 0,   7.5038),
    (ql.Option.Call, 100.0, 100.0, 100.0, 0.05, 0.1,  0.35, 270, 90,  4.0687),
    (ql.Option.Call, 100.0, 100.0, 100.0, 0.05, 0.1,  0.35, 270, 180, 1.4222),
    (ql.Option.Call, 100.0, 100.0, 105.0, 0.05, 0.1,  0.15, 270, 0,   1.6729),
    (ql.Option.Call, 100.0, 100.0, 105.0, 0.05, 0.1,  0.15, 270, 90,  0.3565),
    (ql.Option.Call, 100.0, 100.0, 105.0, 0.05, 0.1,  0.15, 270, 180, 0.0004),
    (ql.Option.Call, 100.0, 100.0, 105.0, 0.05, 0.1,  0.35, 270, 0,   5.4071),
    (ql.Option.Call, 100.0, 100.0, 105.0, 0.05, 0.1,  0.35, 270, 90,  2.1359),
    (ql.Option.Call, 100.0, 100.0, 105.0, 0.05, 0.1,  0.35, 270, 180, 0.1552),
]


# ---------------------------------------------------------------------------
# Test 11 – Levy engine
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "opt_type,spot_val,current_avg,strike,div_yield,rfr,vol,length,elapsed,result",
    LEVY_CASES
)
def test_levy_engine(opt_type, spot_val, current_avg, strike, div_yield, rfr, vol,
                     length, elapsed, result):
    """Levy engine for Asian options (Haug p.99-100)."""
    dc = ql.Actual360()
    today = ql.Settings.instance().evaluationDate

    spot = ql.SimpleQuote(spot_val)
    qTS = ql.YieldTermStructureHandle(ql.FlatForward(today, div_yield, dc))
    rTS = ql.YieldTermStructureHandle(ql.FlatForward(today, rfr, dc))
    volTS = ql.BlackVolTermStructureHandle(ql.BlackConstantVol(today, ql.NullCalendar(), vol, dc))
    process = ql.BlackScholesMertonProcess(ql.QuoteHandle(spot), qTS, rTS, volTS)

    average = ql.SimpleQuote(current_avg)
    payoff = ql.PlainVanillaPayoff(opt_type, strike)
    start_date = today - elapsed
    maturity = start_date + length
    exercise = ql.EuropeanExercise(maturity)

    engine = ql.ContinuousArithmeticAsianLevyEngine(process, ql.QuoteHandle(average))
    option = ql.ContinuousAveragingAsianOption(ql.Average.Arithmetic, start_date, payoff, exercise)
    option.setPricingEngine(engine)

    calculated = option.NPV()
    tolerance = 1.0e-4
    assert abs(calculated - result) <= tolerance, (
        f"Levy engine: expected {result}, got {calculated}"
    )


# ---------------------------------------------------------------------------
# Test 13 – Analytic continuous geometric average-price under Heston
# ---------------------------------------------------------------------------

def test_analytic_continuous_geometric_average_price_heston():
    """Analytic continuous geometric Asians under Heston (Kim & Wee 2011)."""
    days   = [73, 73, 73, 73, 73, 548, 548, 548, 548, 548, 1095, 1095, 1095, 1095, 1095]
    strikes = [90.0, 95.0, 100.0, 105.0, 110.0] * 3
    prices = [
        10.6571, 6.5871, 3.4478, 1.4552, 0.4724,
        16.5030, 13.7625, 11.3374, 9.2245, 7.4122,
        20.5102, 18.3060, 16.2895, 14.4531, 12.7882,
    ]
    prices_2 = [
        10.6425, 6.4362, 3.1578, 1.1936, 0.3609,
        14.9955, 11.6707, 8.7767, 6.3818, 4.5118,
        18.1219, 15.2009, 12.5707, 10.2539, 8.2611,
    ]

    tolerance = 1.0e-2
    dc = ql.Actual365Fixed()
    today = ql.Settings.instance().evaluationDate

    spot = ql.QuoteHandle(ql.SimpleQuote(100.0))
    qTS = ql.YieldTermStructureHandle(ql.FlatForward(today, 0.0, dc))
    rTS = ql.YieldTermStructureHandle(ql.FlatForward(today, 0.05, dc))

    # Table 1 parameters
    process1 = ql.HestonProcess(rTS, qTS, spot, 0.09, 1.15, 0.348, 0.39, -0.64)
    engine1 = ql.AnalyticContinuousGeometricAveragePriceAsianHestonEngine(process1)

    for i, (day, strike, expected) in enumerate(zip(days, strikes, prices)):
        expiry = ql.EuropeanExercise(today + day)
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, strike)
        option = ql.ContinuousAveragingAsianOption(ql.Average.Geometric, payoff, expiry)
        option.setPricingEngine(engine1)
        calculated = option.NPV()
        assert abs(calculated - expected) <= tolerance, (
            f"Heston cont geom [table1] i={i}: expected {expected}, got {calculated}"
        )

    # Table 4 parameters
    process2 = ql.HestonProcess(rTS, qTS, spot, 0.09, 2.0, 0.09, 1.0, -0.3)
    engine2 = ql.AnalyticContinuousGeometricAveragePriceAsianHestonEngine(process2)

    for i, (day, strike, expected) in enumerate(zip(days, strikes, prices_2)):
        expiry = ql.EuropeanExercise(today + day)
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, strike)
        option = ql.ContinuousAveragingAsianOption(ql.Average.Geometric, payoff, expiry)
        option.setPricingEngine(engine2)
        calculated = option.NPV()
        assert abs(calculated - expected) <= tolerance, (
            f"Heston cont geom [table4] i={i}: expected {expected}, got {calculated}"
        )

    # Kim et al 2016 data
    days_3   = [30,91,182,365,730,1095, 30,91,182,365,730,1095, 30,91,182,365,730,1095]
    strikes_3 = [90]*6 + [100]*6 + [110]*6
    tol_3     = [2e-2,1e-2,1e-2,1e-2,1e-2,1e-2]*3
    prices_3  = [
        10.1513, 10.8175, 11.8664, 13.5931, 16.0988, 17.9475,
        2.0472,  3.5735,  5.0588,  7.1132,  9.9139,  11.9959,
        0.0350,  0.4869,  1.3376,  2.8569,  5.2804,  7.2682,
    ]
    process3 = ql.HestonProcess(rTS, qTS, spot, 0.09, 1.15, 0.0348, 0.39, -0.64)
    engine3 = ql.AnalyticContinuousGeometricAveragePriceAsianHestonEngine(process3)

    for day, strike, expected, tol_i in zip(days_3, strikes_3, prices_3, tol_3):
        expiry = ql.EuropeanExercise(today + day)
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, strike)
        option = ql.ContinuousAveragingAsianOption(ql.Average.Geometric, payoff, expiry)
        option.setPricingEngine(engine3)
        calculated = option.NPV()
        assert abs(calculated - expected) <= tol_i, (
            f"Heston cont geom [2016] day={day} K={strike}: expected {expected}, got {calculated}"
        )


# ---------------------------------------------------------------------------
# Test 14 – Analytic discrete geometric average-price under Heston
# ---------------------------------------------------------------------------

def _discrete_geom_heston_data():
    days = [30,91,182,365,730,1095]*3
    strikes = [90]*6 + [100]*6 + [110]*6
    prices = [
        10.2732, 10.9554, 11.9916, 13.6950, 16.1773, 18.0146,
        2.4389,  3.7881,  5.2132,  7.2243,  9.9948,  12.0639,
        0.1012,  0.5949,  1.4444,  2.9479,  5.3531,  7.3315,
    ]
    tols = [
        3e-2,2e-2,2e-2,2e-2,3e-2,4e-2,
        8e-2,1e-2,2e-2,3e-2,3e-2,4e-2,
        2e-2,1e-2,1e-2,2e-2,3e-2,4e-2,
    ]
    return days, strikes, prices, tols


def test_analytic_discrete_geometric_average_price_heston():
    """Analytic discrete geometric average-price Asians under Heston."""
    days, strikes, prices, tols = _discrete_geom_heston_data()

    dc = ql.Actual365Fixed()
    today = ql.Settings.instance().evaluationDate
    spot = ql.QuoteHandle(ql.SimpleQuote(100.0))
    qTS = ql.YieldTermStructureHandle(ql.FlatForward(today, 0.0, dc))
    rTS = ql.YieldTermStructureHandle(ql.FlatForward(today, 0.05, dc))

    process = ql.HestonProcess(rTS, qTS, spot, 0.09, 1.15, 0.0348, 0.39, -0.64)
    engine = ql.AnalyticDiscreteGeometricAveragePriceAsianHestonEngine(process)

    for day, strike, expected, tol_i in zip(days, strikes, prices, tols):
        future_fixings = int(math.floor(day / 7.0))
        expiry_date = today + day
        fixing_dates = [expiry_date - i * 7 for i in range(future_fixings - 1, -1, -1)]
        exercise = ql.EuropeanExercise(expiry_date)
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, strike)
        option = ql.DiscreteAveragingAsianOption(
            ql.Average.Geometric, 1.0, 0, fixing_dates, payoff, exercise
        )
        option.setPricingEngine(engine)
        calculated = option.NPV()
        assert abs(calculated - expected) <= tol_i, (
            f"Discrete geom Heston day={day} K={strike}: expected {expected}, got {calculated}"
        )


# ---------------------------------------------------------------------------
# Test 15 – Turnbull-Wakeman engine with term structure
# ---------------------------------------------------------------------------

TW_CASES = [
    (ql.Option.Call, 100, 80,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "flat", 19.5152),
    (ql.Option.Call, 100, 80,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "up",   19.5063),
    (ql.Option.Call, 100, 80,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "down", 19.5885),
    (ql.Option.Put,  100, 80,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "flat", 0.0090),
    (ql.Option.Put,  100, 80,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "up",   0.0001),
    (ql.Option.Put,  100, 80,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "down", 0.0823),
    (ql.Option.Call, 100, 90,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "flat", 10.1437),
    (ql.Option.Call, 100, 90,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "up",   9.8313),
    (ql.Option.Call, 100, 90,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "down", 10.7062),
    (ql.Option.Put,  100, 90,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "flat", 0.3906),
    (ql.Option.Put,  100, 90,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "up",   0.0782),
    (ql.Option.Put,  100, 90,  0, 0.05, 1.0/52, 0.5, 26, 0.2, "down", 0.9531),
    (ql.Option.Call, 100, 100, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "flat", 3.2700),
    (ql.Option.Call, 100, 100, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "up",   2.2819),
    (ql.Option.Call, 100, 100, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "down", 4.3370),
    (ql.Option.Put,  100, 100, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "flat", 3.2700),
    (ql.Option.Put,  100, 100, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "up",   2.2819),
    (ql.Option.Put,  100, 100, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "down", 4.3370),
    (ql.Option.Call, 100, 110, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "flat", 0.5515),
    (ql.Option.Call, 100, 110, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "up",   0.1314),
    (ql.Option.Call, 100, 110, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "down", 1.2429),
    (ql.Option.Put,  100, 110, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "flat", 10.3046),
    (ql.Option.Put,  100, 110, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "up",   9.8845),
    (ql.Option.Put,  100, 110, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "down", 10.9960),
    (ql.Option.Call, 100, 120, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "flat", 0.0479),
    (ql.Option.Call, 100, 120, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "up",   0.0016),
    (ql.Option.Call, 100, 120, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "down", 0.2547),
    (ql.Option.Put,  100, 120, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "flat", 19.5541),
    (ql.Option.Put,  100, 120, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "up",   19.5078),
    (ql.Option.Put,  100, 120, 0, 0.05, 1.0/52, 0.5, 26, 0.2, "down", 19.7609),
]


def _time_to_days_360(t):
    return int(round(t * 360))


@pytest.mark.parametrize(
    "opt_type,underlying,strike,b,rfr,first,expiry,fixings,base_vol,slope,result",
    TW_CASES
)
def test_turnbull_wakeman_engine(opt_type, underlying, strike, b, rfr,
                                  first, expiry, fixings, base_vol, slope, result):
    """Turnbull-Wakeman engine with term structure (Haug Table 4-28)."""
    dc = ql.Actual360()
    today = ql.Settings.instance().evaluationDate

    dt = (expiry - first) / (fixings - 1)
    fixing_dates = [today + _time_to_days_360(i * dt + first) for i in range(fixings)]

    spot = ql.SimpleQuote(underlying)
    qTS = ql.YieldTermStructureHandle(ql.FlatForward(today, rfr - b, dc))
    rTS = ql.YieldTermStructureHandle(ql.FlatForward(today, rfr, dc))

    vol_slope = 0.005
    if slope == "flat":
        volTS = ql.BlackVolTermStructureHandle(
            ql.BlackConstantVol(today, ql.NullCalendar(), base_vol, dc)
        )
    elif slope == "up":
        vols = [base_vol - (fixings - 1) * vol_slope + i * vol_slope for i in range(fixings)]
        volTS = ql.BlackVolTermStructureHandle(
            ql.BlackVarianceCurve(today, fixing_dates, vols, dc, True)
        )
    else:  # down
        vols = [base_vol + (fixings - 1) * vol_slope - i * vol_slope for i in range(fixings)]
        volTS = ql.BlackVolTermStructureHandle(
            ql.BlackVarianceCurve(today, fixing_dates, vols, dc, False)
        )

    process = ql.BlackScholesMertonProcess(ql.QuoteHandle(spot), qTS, rTS, volTS)
    maturity = today + _time_to_days_360(expiry)
    exercise = ql.EuropeanExercise(maturity)
    payoff = ql.PlainVanillaPayoff(opt_type, strike)

    option = ql.DiscreteAveragingAsianOption(
        ql.Average.Arithmetic, 0.0, 0, fixing_dates, payoff, exercise
    )
    option.setPricingEngine(ql.TurnbullWakemanAsianEngine(process))

    calculated = option.NPV()
    tolerance = 2.5e-3
    assert abs(calculated - result) <= tolerance, (
        f"TW engine [{slope}]: expected {result}, got {calculated}"
    )

    # Check analytical delta and gamma against bump
    dS = 0.001
    delta = option.delta()
    gamma = option.gamma()

    spot_up = ql.SimpleQuote(underlying + dS)
    spot_down = ql.SimpleQuote(underlying - dS)
    proc_up = ql.BlackScholesMertonProcess(ql.QuoteHandle(spot_up), qTS, rTS, volTS)
    proc_dn = ql.BlackScholesMertonProcess(ql.QuoteHandle(spot_down), qTS, rTS, volTS)

    option.setPricingEngine(ql.TurnbullWakemanAsianEngine(proc_up))
    npv_up = option.NPV()
    option.setPricingEngine(ql.TurnbullWakemanAsianEngine(proc_dn))
    npv_dn = option.NPV()

    delta_bump = (npv_up - npv_dn) / (2 * dS)
    gamma_bump = (npv_up + npv_dn - 2 * calculated) / (dS * dS)

    tol_greek = 1.0e-6
    assert abs(delta_bump - delta) <= tol_greek, (
        f"TW delta [{slope}]: bump={delta_bump}, analytic={delta}"
    )
    assert abs(gamma_bump - gamma) <= tol_greek, (
        f"TW gamma [{slope}]: bump={gamma_bump}, analytic={gamma}"
    )


# ---------------------------------------------------------------------------
# Test 16 – Seasoned continuous Asian options
# ---------------------------------------------------------------------------

def test_continuous_seasoned_asian_options():
    """Seasoned continuous averaging Asian options."""
    dc = ql.Actual365Fixed()
    today = ql.Date(15, ql.November, 2025)
    settlement_date = ql.Date(17, ql.November, 2025)
    ql.Settings.instance().evaluationDate = today

    spot = 100.0
    div_yield = 0.03
    rfr = 0.06
    vol = 0.20
    maturity = ql.Date(17, ql.November, 2026)
    start_date = ql.Date(17, ql.August, 2025)
    strike = 100.0

    underlying_h = ql.QuoteHandle(ql.SimpleQuote(spot))
    rTS = ql.YieldTermStructureHandle(ql.FlatForward(settlement_date, rfr, dc))
    qTS = ql.YieldTermStructureHandle(ql.FlatForward(settlement_date, div_yield, dc))
    volTS = ql.BlackVolTermStructureHandle(
        ql.BlackConstantVol(settlement_date, ql.TARGET(), vol, dc)
    )
    process = ql.BlackScholesMertonProcess(underlying_h, qTS, rTS, volTS)

    payoff = ql.PlainVanillaPayoff(ql.Option.Put, strike)
    exercise = ql.EuropeanExercise(maturity)

    # Fresh (unseasoned) option
    fresh_option = ql.ContinuousAveragingAsianOption(
        ql.Average.Arithmetic, settlement_date, payoff, exercise
    )
    fresh_option.setPricingEngine(
        ql.ContinuousArithmeticAsianLevyEngine(process, ql.QuoteHandle(ql.SimpleQuote(0.0)))
    )
    fresh_npv = fresh_option.NPV()

    # Seasoned option with low average
    current_average = 98.5
    seasoned_option = ql.ContinuousAveragingAsianOption(
        ql.Average.Arithmetic, start_date, payoff, exercise
    )
    seasoned_option.setPricingEngine(
        ql.ContinuousArithmeticAsianLevyEngine(
            process, ql.QuoteHandle(ql.SimpleQuote(current_average))
        )
    )
    seasoned_npv = seasoned_option.NPV()

    assert seasoned_npv < fresh_npv, (
        f"Seasoned put ({seasoned_npv}) should be < fresh put ({fresh_npv}) "
        f"when avg ({current_average}) < strike ({strike})"
    )

    # Seasoned option with high average
    high_average = 102.0
    seasoned_high_option = ql.ContinuousAveragingAsianOption(
        ql.Average.Arithmetic, start_date, payoff, exercise
    )
    seasoned_high_option.setPricingEngine(
        ql.ContinuousArithmeticAsianLevyEngine(
            process, ql.QuoteHandle(ql.SimpleQuote(high_average))
        )
    )
    seasoned_high_npv = seasoned_high_option.NPV()

    assert seasoned_high_npv < seasoned_npv, (
        f"Higher average ({high_average}, NPV={seasoned_high_npv}) should give "
        f"lower NPV than lower average ({current_average}, NPV={seasoned_npv})"
    )

    # Continuous geometric should raise for seasoned options
    seasoned_geom = ql.ContinuousAveragingAsianOption(
        ql.Average.Geometric, start_date, payoff, exercise
    )
    seasoned_geom.setPricingEngine(
        ql.AnalyticContinuousGeometricAveragePriceAsianEngine(process)
    )
    with pytest.raises(Exception):
        seasoned_geom.NPV()

    ql.Settings.instance().evaluationDate = ql.Date.todaysDate()
