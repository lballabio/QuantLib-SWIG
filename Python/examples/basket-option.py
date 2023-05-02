# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Basket options
#
# Copyright (&copy;) 2004, 2005, 2006 StatPro Italia srl
#
# This file is part of QuantLib, a free-software/open-source library
# for financial quantitative analysts and developers - https://www.quantlib.org/
#
# QuantLib is free software: you can redistribute it and/or modify it under the
# terms of the QuantLib license.  You should have received a copy of the
# license along with this program; if not, please email
# <quantlib-dev@lists.sf.net>. The license is also available online at
# <https://www.quantlib.org/license.shtml>.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the license for more details.

import QuantLib as ql

# ### Global data

todaysDate = ql.Date(15, ql.May, 1998)
ql.Settings.instance().evaluationDate = todaysDate
settlementDate = ql.Date(17, ql.May, 1998)
riskFreeRate = ql.FlatForward(settlementDate, 0.05, ql.Actual365Fixed())

# ### Option parameters

exercise = ql.EuropeanExercise(ql.Date(17, ql.May, 1999))
payoff = ql.PlainVanillaPayoff(ql.Option.Call, 8.0)

# ### Market data


underlying1 = ql.SimpleQuote(7.0)
volatility1 = ql.BlackConstantVol(todaysDate, ql.TARGET(), 0.10, ql.Actual365Fixed())
dividendYield1 = ql.FlatForward(settlementDate, 0.05, ql.Actual365Fixed())
underlying2 = ql.SimpleQuote(7.0)
volatility2 = ql.BlackConstantVol(todaysDate, ql.TARGET(), 0.10, ql.Actual365Fixed())
dividendYield2 = ql.FlatForward(settlementDate, 0.05, ql.Actual365Fixed())


process1 = ql.BlackScholesMertonProcess(
    ql.QuoteHandle(underlying1),
    ql.YieldTermStructureHandle(dividendYield1),
    ql.YieldTermStructureHandle(riskFreeRate),
    ql.BlackVolTermStructureHandle(volatility1),
)

process2 = ql.BlackScholesMertonProcess(
    ql.QuoteHandle(underlying2),
    ql.YieldTermStructureHandle(dividendYield2),
    ql.YieldTermStructureHandle(riskFreeRate),
    ql.BlackVolTermStructureHandle(volatility2),
)

matrix = ql.Matrix(2, 2)
matrix[0][0] = 1.0
matrix[1][1] = 1.0
matrix[0][1] = 0.5
matrix[1][0] = 0.5

process = ql.StochasticProcessArray([process1, process2], matrix)

# ### Pricing

basketoption = ql.BasketOption(ql.MaxBasketPayoff(payoff), exercise)
basketoption.setPricingEngine(
    ql.MCEuropeanBasketEngine(process, "pseudorandom", timeStepsPerYear=1, requiredTolerance=0.02, seed=42)
)
print("Maximum Basket Payoff: ", basketoption.NPV())

basketoption = ql.BasketOption(ql.MinBasketPayoff(payoff), exercise)
basketoption.setPricingEngine(
    ql.MCEuropeanBasketEngine(process, "pseudorandom", timeStepsPerYear=1, requiredTolerance=0.02, seed=42)
)
print("Minimum Basket Payoff: ", basketoption.NPV())

basketoption = ql.BasketOption(ql.AverageBasketPayoff(payoff, 2), exercise)
basketoption.setPricingEngine(
    ql.MCEuropeanBasketEngine(process, "pseudorandom", timeStepsPerYear=1, requiredTolerance=0.02, seed=42)
)
print("Average Basket Payoff: ", basketoption.NPV())

americanExercise = ql.AmericanExercise(settlementDate, ql.Date(17, ql.May, 1999))
americanbasketoption = ql.BasketOption(ql.MaxBasketPayoff(payoff), americanExercise)
americanbasketoption.setPricingEngine(
    ql.MCAmericanBasketEngine(
        process,
        "pseudorandom",
        timeSteps=10,
        requiredTolerance=0.02,
        seed=42,
        polynomOrder=5,
        polynomType=ql.LsmBasisSystem.Hermite,
    )
)
print("Basket American Exercise: ", americanbasketoption.NPV())
