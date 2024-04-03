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

# # Bermudan swaptions
#
# Copyright (&copy;) 2004, 2005, 2006, 2007 StatPro Italia srl
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
import pandas as pd

# ### Setup

todaysDate = ql.Date(15, ql.February, 2002)
ql.Settings.instance().evaluationDate = todaysDate
calendar = ql.TARGET()
settlementDate = ql.Date(19, ql.February, 2002)


def calibrate(model, helpers, l, name):
    print("Model: %s" % name)

    method = ql.Simplex(l)
    model.calibrate(helpers, method, ql.EndCriteria(1000, 250, 1e-7, 1e-7, 1e-7))

    print("Parameters: %s" % model.params())

    totalError = 0.0
    data = []
    for swaption, helper in zip(swaptionVols, helpers):
        maturity, length, vol = swaption
        NPV = helper.modelValue()
        implied = helper.impliedVolatility(NPV, 1.0e-4, 1000, 0.05, 0.50)
        error = implied - vol
        totalError += abs(error)
        data.append((maturity, length, vol, implied, error))
    averageError = totalError / len(helpers)

    print(pd.DataFrame(data, columns=["maturity", "length", "volatility", "implied", "error"]))

    print("Average error: %.4f" % averageError)


# ### Market data

swaptionVols = [
    # maturity,          length,             volatility
    (ql.Period(1, ql.Years), ql.Period(5, ql.Years), 0.1148),
    (ql.Period(2, ql.Years), ql.Period(4, ql.Years), 0.1108),
    (ql.Period(3, ql.Years), ql.Period(3, ql.Years), 0.1070),
    (ql.Period(4, ql.Years), ql.Period(2, ql.Years), 0.1021),
    (ql.Period(5, ql.Years), ql.Period(1, ql.Years), 0.1000),
]

# This is a flat yield term structure implying a 1x5 swap at 5%.

rate = ql.makeQuoteHandle(0.04875825)
termStructure = ql.YieldTermStructureHandle(ql.FlatForward(settlementDate, rate, ql.Actual365Fixed()))

# Define the ATM/OTM/ITM swaps:

swapEngine = ql.DiscountingSwapEngine(termStructure)

fixedLegFrequency = ql.Annual
fixedLegTenor = ql.Period(1, ql.Years)
fixedLegConvention = ql.Unadjusted
floatingLegConvention = ql.ModifiedFollowing
fixedLegDayCounter = ql.Thirty360(ql.Thirty360.European)
floatingLegFrequency = ql.Semiannual
floatingLegTenor = ql.Period(6, ql.Months)

payFixed = ql.Swap.Payer
fixingDays = 2
index = ql.Euribor6M(termStructure)
floatingLegDayCounter = index.dayCounter()

swapStart = calendar.advance(settlementDate, 1, ql.Years, floatingLegConvention)
swapEnd = calendar.advance(swapStart, 5, ql.Years, floatingLegConvention)

fixedSchedule = ql.Schedule(
    swapStart,
    swapEnd,
    fixedLegTenor,
    calendar,
    fixedLegConvention,
    fixedLegConvention,
    ql.DateGeneration.Forward,
    False,
)
floatingSchedule = ql.Schedule(
    swapStart,
    swapEnd,
    floatingLegTenor,
    calendar,
    floatingLegConvention,
    floatingLegConvention,
    ql.DateGeneration.Forward,
    False,
)

dummy = ql.VanillaSwap(
    payFixed, 100.0, fixedSchedule, 0.0, fixedLegDayCounter, floatingSchedule, index, 0.0, floatingLegDayCounter
)
dummy.setPricingEngine(swapEngine)
atmRate = dummy.fairRate()

atmSwap = ql.VanillaSwap(
    payFixed, 1000.0, fixedSchedule, atmRate, fixedLegDayCounter,
    floatingSchedule, index, 0.0, floatingLegDayCounter
)

otmSwap = ql.VanillaSwap(
    payFixed, 1000.0, fixedSchedule, atmRate * 1.2, fixedLegDayCounter,
    floatingSchedule, index, 0.0, floatingLegDayCounter
)

itmSwap = ql.VanillaSwap(
    payFixed, 1000.0, fixedSchedule, atmRate * 0.8, fixedLegDayCounter,
    floatingSchedule, index, 0.0, floatingLegDayCounter
)

atmSwap.setPricingEngine(swapEngine)
otmSwap.setPricingEngine(swapEngine)
itmSwap.setPricingEngine(swapEngine)

helpers = [
    ql.SwaptionHelper(
        maturity,
        length,
        ql.makeQuoteHandle(vol),
        index,
        index.tenor(),
        index.dayCounter(),
        index.dayCounter(),
        termStructure,
    )
    for maturity, length, vol in swaptionVols
]

times = {}
for h in helpers:
    for t in h.times():
        times[t] = 1
times = sorted(times.keys())

grid = ql.TimeGrid(times, 30)

G2model = ql.G2(termStructure)
HWmodel = ql.HullWhite(termStructure)
HWmodel2 = ql.HullWhite(termStructure)
BKmodel = ql.BlackKarasinski(termStructure)

# ### Calibrations

for h in helpers:
    h.setPricingEngine(ql.G2SwaptionEngine(G2model, 6.0, 16))
calibrate(G2model, helpers, 0.05, "G2 (analytic formulae)")

for h in helpers:
    h.setPricingEngine(ql.JamshidianSwaptionEngine(HWmodel))
calibrate(HWmodel, helpers, 0.05, "Hull-White (analytic formulae)")

for h in helpers:
    h.setPricingEngine(ql.TreeSwaptionEngine(HWmodel2, grid))
calibrate(HWmodel2, helpers, 0.05, "Hull-White (numerical calibration)")

for h in helpers:
    h.setPricingEngine(ql.TreeSwaptionEngine(BKmodel, grid))
calibrate(BKmodel, helpers, 0.05, "Black-Karasinski (numerical calibration)")


# ### Price Bermudan swaptions on defined swaps

bermudanDates = [d for d in fixedSchedule][:-1]
exercise = ql.BermudanExercise(bermudanDates)

atmSwaption = ql.Swaption(atmSwap, exercise)
otmSwaption = ql.Swaption(otmSwap, exercise)
itmSwaption = ql.Swaption(itmSwap, exercise)

data = []

# +
atmSwaption.setPricingEngine(ql.TreeSwaptionEngine(G2model, 50))
otmSwaption.setPricingEngine(ql.TreeSwaptionEngine(G2model, 50))
itmSwaption.setPricingEngine(ql.TreeSwaptionEngine(G2model, 50))

data.append(("G2 analytic", itmSwaption.NPV(), atmSwaption.NPV(), otmSwaption.NPV()))

# +
atmSwaption.setPricingEngine(ql.TreeSwaptionEngine(HWmodel, 50))
otmSwaption.setPricingEngine(ql.TreeSwaptionEngine(HWmodel, 50))
itmSwaption.setPricingEngine(ql.TreeSwaptionEngine(HWmodel, 50))

data.append(("HW analytic", itmSwaption.NPV(), atmSwaption.NPV(), otmSwaption.NPV()))

# +
atmSwaption.setPricingEngine(ql.TreeSwaptionEngine(HWmodel2, 50))
otmSwaption.setPricingEngine(ql.TreeSwaptionEngine(HWmodel2, 50))
itmSwaption.setPricingEngine(ql.TreeSwaptionEngine(HWmodel2, 50))

data.append(("HW numerical", itmSwaption.NPV(), atmSwaption.NPV(), otmSwaption.NPV()))

# +
atmSwaption.setPricingEngine(ql.TreeSwaptionEngine(BKmodel, 50))
otmSwaption.setPricingEngine(ql.TreeSwaptionEngine(BKmodel, 50))
itmSwaption.setPricingEngine(ql.TreeSwaptionEngine(BKmodel, 50))

data.append(("BK numerical", itmSwaption.NPV(), atmSwaption.NPV(), otmSwaption.NPV()))
# -

print(pd.DataFrame(data, columns=["model", "in-the-money", "at-the-money", "out-of-the-money"]))
