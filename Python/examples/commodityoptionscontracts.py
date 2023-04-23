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

# # Treasury Futures Contracts
#
#
# This file is part of QuantLib, a free-software/open-source library
# for financial quantitative analysts and developers - https://www.quantlib.org/
#
# QuantLib is free software: you can redistribute it and/or modify it
# under the terms of the QuantLib license.  You should have received a
# # copy of the license along with this program; if not, please email
# <quantlib-dev@lists.sf.net>. The license is also available online at
# <https://www.quantlib.org/license.shtml>.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the license for more details.

#    This example shows how to set up a term structure and then price
#    some simple bonds. The last part is dedicated to peripherical
#    computations such as "Yield to Price" or "Price to Yield"

# ### Imports
import QuantLib as ql
import math

# ### Global Data
calendar = ql.UnitedStates()
businessConvention = ql.ModifiedFollowing
settlementDays = 0
daysCount = ql.ActualActual(ql.ActualActual.Bond)

# ### Option on Treasury Futures Contracts

# ------------------------------------------------------- #

# ------------------------------------------------------- #

# ### 10-Year Note Options (Put)

calendar = ql.UnitedStates()
business_convention = ql.ModifiedFollowing
settlement_days = 0
dayCount = ql.ActualActual(ql.ActualActual.Bond)

interestRate = 0.0011
calcDate = ql.Date(1, 1, 2022)
yieldCurve = ql.FlatForward(calcDate, interestRate, dayCount, ql.Compounded, ql.Continuous)

ql.Settings.instance().evaluationDate = calcDate

maturity = ql.Date(28, 6, 2022)

strike = 120
spot = 125
volatility = 0.20
type = ql.Option.Put

discount = yieldCurve.discount(maturity)

time = yieldCurve.dayCounter().yearFraction(calcDate, maturity)

stddev = volatility*math.sqrt(time)

strikepayoff = ql.PlainVanillaPayoff(type, strike)

black = ql.BlackCalculator(strikepayoff, spot, stddev, discount)

print("%-20s: %4.4f" %("Option values on Treasury Futures Put", black.value()))
print("%-20s: %4.4f" %("Delta", black.delta(spot)))
print("%-20s: %4.4f" %("Gamma", black.gamma(spot)))
print("%-20s: %4.4f" %("Theta", black.theta(spot, T)))
print("%-20s: %4.4f" %("Vega", black.vega(T)))
print("%-20s: %4.4f" %("Rho", black.rho( T)))

# ------------------------------------------------------- #

# ### Natural Silver Futures Options (Call)

interestRate = 0.05
calcDate = ql.Date(1, 1, 2022)
yieldCurve = ql.FlatForward(calcDate, interestRate, dayCount, ql.Compounded, ql.Continuous)

ql.Settings.instance().evaluationDate = calcDate
T = 31/365.

strike = 120
spot = 125
volatility = 0.20
type = ql.Option.Put

discount = yieldCurve.discount(T)
stddev = volatility*math.sqrt(T)
strikepayoff = ql.PlainVanillaPayoff(type, strike)
black = ql.BlackCalculator(strikepayoff, spot, stddev, discount)

print("%-20s: %4.4f" %("Option Commodity Price for Put", black.value()))
print("%-20s: %4.4f" %("Delta", black.delta(spot)))
print("%-20s: %4.4f" %("Gamma", black.gamma(spot)))
print("%-20s: %4.4f" %("Theta", black.theta(spot, T)))
print("%-20s: %4.4f" %("Vega", black.vega(T)))
print("%-20s: %4.4f" %("Rho", black.rho( T)))
