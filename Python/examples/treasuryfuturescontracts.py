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

# ### 2-Year Note 

interestRate = 0.003
calcDate = ql.Date(1,12,2023)
yieldCurve = ql.FlatForward(calcDate, interestRate, daysCount, ql.Compounded, ql.Continuous)

ql.Settings.instance().evaluationDate = calcDate
optionMaturityDate = ql.Date(24,12,2025)
strike = 100
spot = 125 # spot price is the futures price
volatility = 20/100.
optionType = ql.Option.Call

discount = yieldCurve.discount(optionMaturityDate)
strikepayoff = ql.PlainVanillaPayoff(optionType, strike)
T = yieldCurve.dayCounter().yearFraction(calcDate, optionMaturityDate)

stddev = volatility*math.sqrt(T)

black = ql.BlackCalculator(strikepayoff, spot, stddev, discount)

print("%-20s: %4.4f" %("Option Price for Treasury Futures Contract (2-Year Note)", black.value() )) 
print("%-20s: %4.4f" %("Delta", black.delta(spot)))
print("%-20s: %4.4f" %("Gamma", black.gamma(spot)))
print("%-20s: %4.4f" %("Theta", black.theta(spot, T))) 
print("%-20s: %4.4f" %("Vega", black.vega(T)))
print("%-20s: %4.4f" %("Rho", black.rho( T)))

# ------------------------------------------------------- #

# ### 3-Year Note

interestRate = 0.003
calcDate = ql.Date(1,12,2023)
yieldCurve = ql.FlatForward(calcDate, interestRate, daysCount, ql.Compounded, ql.Continuous)

ql.Settings.instance().evaluationDate = calcDate
optionMaturityDate = ql.Date(24,12,2026)
strike = 100
spot = 125 # spot price is the futures price
volatility = 20/100.
optionType = ql.Option.Call

discount = yieldCurve.discount(optionMaturityDate)
strikepayoff = ql.PlainVanillaPayoff(optionType, strike)
T = yieldCurve.dayCounter().yearFraction(calcDate, optionMaturityDate)

stddev = volatility*math.sqrt(T)
black = ql.BlackCalculator(strikepayoff, spot, stddev, discount)

print("%-20s: %4.4f" %("Option Price for Treasury Futures Contract (3-Year Note)", black.value() )) 
print("%-20s: %4.4f" %("Delta", black.delta(spot)))
print("%-20s: %4.4f" %("Gamma", black.gamma(spot)))
print("%-20s: %4.4f" %("Theta", black.theta(spot, T))) 
print("%-20s: %4.4f" %("Vega", black.vega(T)))
print("%-20s: %4.4f" %("Rho", black.rho( T)))

# ------------------------------------------------------- #

# ### 5-Year Note

interestRate = 0.003
calcDate = ql.Date(1,12,2023)
yieldCurve = ql.FlatForward(calcDate, interestRate, daysCount, ql.Compounded, ql.Continuous)

ql.Settings.instance().evaluationDate = calcDate
optionMaturityDate = ql.Date(24,12,2028)
strike = 100
spot = 125 # spot price is the futures price
volatility = 20/100.
optionType = ql.Option.Call

discount = yieldCurve.discount(optionMaturityDate)
strikepayoff = ql.PlainVanillaPayoff(optionType, strike)
T = yieldCurve.dayCounter().yearFraction(calcDate, optionMaturityDate)

stddev = volatility*math.sqrt(T)

black = ql.BlackCalculator(strikepayoff, spot, stddev, discount)

print("%-20s: %4.4f" %("Option Price for Treasury Futures Contract (5-Year Note)", black.value() )) 
print("%-20s: %4.4f" %("Delta", black.delta(spot)))
print("%-20s: %4.4f" %("Gamma", black.gamma(spot)))
print("%-20s: %4.4f" %("Theta", black.theta(spot, T))) 
print("%-20s: %4.4f" %("Vega", black.vega(T)))
print("%-20s: %4.4f" %("Rho", black.rho( T)))

# ------------------------------------------------------- #

# ### 10-Year Note

interestRate = 0.003
calcDate = ql.Date(1,12,2023)
yieldCurve = ql.FlatForward(calcDate, interestRate, daysCount, ql.Compounded, ql.Continuous)

ql.Settings.instance().evaluationDate = calcDate
optionMaturityDate = ql.Date(24,12,2033)
strike = 100
spot = 125 # futures price
volatility = 20/100.
optionType = ql.Option.Call

discount = yieldCurve.discount(optionMaturityDate)
strikepayoff = ql.PlainVanillaPayoff(optionType, strike)
T = yieldCurve.dayCounter().yearFraction(calcDate, optionMaturityDate)

stddev = volatility*math.sqrt(T)
black = ql.BlackCalculator(strikepayoff, spot, stddev, discount)

print("%-20s: %4.4f" %("Option Price for Treasury Futures Contract (10-Year Note)", black.value() )) 
print("%-20s: %4.4f" %("Delta", black.delta(spot)))
print("%-20s: %4.4f" %("Gamma", black.gamma(spot)))
print("%-20s: %4.4f" %("Theta", black.theta(spot, T))) 
print("%-20s: %4.4f" %("Vega", black.vega(T)))
print("%-20s: %4.4f" %("Rho", black.rho( T)))

# ------------------------------------------------------- #

