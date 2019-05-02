# Copyright (C) 2004, 2005, 2006, 2007 StatPro Italia srl
#
# This file is part of QuantLib, a free-software/open-source library
# for financial quantitative analysts and developers - http://quantlib.org/
#
# QuantLib is free software: you can redistribute it and/or modify it under the
# terms of the QuantLib license.  You should have received a copy of the
# license along with this program; if not, please email
# <quantlib-dev@lists.sf.net>. The license is also available online at
# <http://quantlib.org/license.shtml>.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the license for more details.

import QuantLib as ql

# global data
calendar = ql.TARGET()
todaysDate = ql.Date(6, ql.November, 2001)
ql.Settings.instance().evaluationDate = todaysDate
settlementDate = ql.Date(8, ql.November, 2001)

# market quotes
deposits = {
    (1, ql.Weeks): 0.0382,
    (1, ql.Months): 0.0372,
    (3, ql.Months): 0.0363,
    (6, ql.Months): 0.0353,
    (9, ql.Months): 0.0348,
    (1, ql.Years): 0.0345,
}

FRAs = {(3, 6): 0.037125, (6, 9): 0.037125, (9, 12): 0.037125}

futures = {
    ql.Date(19, 12, 2001): 96.2875,
    ql.Date(20, 3, 2002): 96.7875,
    ql.Date(19, 6, 2002): 96.9875,
    ql.Date(18, 9, 2002): 96.6875,
    ql.Date(18, 12, 2002): 96.4875,
    ql.Date(19, 3, 2003): 96.3875,
    ql.Date(18, 6, 2003): 96.2875,
    ql.Date(17, 9, 2003): 96.0875,
}

swaps = {
    (2, ql.Years): 0.037125,
    (3, ql.Years): 0.0398,
    (5, ql.Years): 0.0443,
    (10, ql.Years): 0.05165,
    (15, ql.Years): 0.055175,
}

# convert them to Quote objects
for n, unit in deposits.keys():
    deposits[(n, unit)] = ql.SimpleQuote(deposits[(n, unit)])
for n, m in FRAs.keys():
    FRAs[(n, m)] = ql.SimpleQuote(FRAs[(n, m)])
for d in futures.keys():
    futures[d] = ql.SimpleQuote(futures[d])
for n, unit in swaps.keys():
    swaps[(n, unit)] = ql.SimpleQuote(swaps[(n, unit)])

# build rate helpers

dayCounter = ql.Actual360()
settlementDays = 2
depositHelpers = [
    ql.DepositRateHelper(
        ql.QuoteHandle(deposits[(n, unit)]),
        ql.Period(n, unit),
        settlementDays,
        calendar,
        ql.ModifiedFollowing,
        False,
        dayCounter,
    )
    for n, unit in [(1, ql.Weeks), (1, ql.Months), (3, ql.Months), (6, ql.Months), (9, ql.Months), (1, ql.Years)]
]

dayCounter = ql.Actual360()
settlementDays = 2
fraHelpers = [
    ql.FraRateHelper(
        ql.QuoteHandle(FRAs[(n, m)]), n, m, settlementDays, calendar, ql.ModifiedFollowing, False, dayCounter
    )
    for n, m in FRAs.keys()
]

dayCounter = ql.Actual360()
months = 3
futuresHelpers = [
    ql.FuturesRateHelper(
        ql.QuoteHandle(futures[d]),
        d,
        months,
        calendar,
        ql.ModifiedFollowing,
        True,
        dayCounter,
        ql.QuoteHandle(ql.SimpleQuote(0.0)),
    )
    for d in futures.keys()
]

settlementDays = 2
fixedLegFrequency = ql.Annual
fixedLegTenor = ql.Period(1, ql.Years)
fixedLegAdjustment = ql.Unadjusted
fixedLegDayCounter = ql.Thirty360()
floatingLegFrequency = ql.Semiannual
floatingLegTenor = ql.Period(6, ql.Months)
floatingLegAdjustment = ql.ModifiedFollowing
swapHelpers = [
    ql.SwapRateHelper(
        ql.QuoteHandle(swaps[(n, unit)]),
        ql.Period(n, unit),
        calendar,
        fixedLegFrequency,
        fixedLegAdjustment,
        fixedLegDayCounter,
        ql.Euribor6M(),
    )
    for n, unit in swaps.keys()
]

# term structure handles

discountTermStructure = ql.RelinkableYieldTermStructureHandle()
forecastTermStructure = ql.RelinkableYieldTermStructureHandle()

# term-structure construction

helpers = depositHelpers[:2] + futuresHelpers + swapHelpers[1:]
depoFuturesSwapCurve = ql.PiecewiseFlatForward(settlementDate, helpers, ql.Actual360())

helpers = depositHelpers[:3] + fraHelpers + swapHelpers
depoFraSwapCurve = ql.PiecewiseFlatForward(settlementDate, helpers, ql.Actual360())

# swaps to be priced

swapEngine = ql.DiscountingSwapEngine(discountTermStructure)

nominal = 1000000
length = 5
maturity = calendar.advance(settlementDate, length, ql.Years)
payFixed = True

fixedLegFrequency = ql.Annual
fixedLegAdjustment = ql.Unadjusted
fixedLegDayCounter = ql.Thirty360()
fixedRate = 0.04

floatingLegFrequency = ql.Semiannual
spread = 0.0
fixingDays = 2
index = ql.Euribor6M(forecastTermStructure)
floatingLegAdjustment = ql.ModifiedFollowing
floatingLegDayCounter = index.dayCounter()

fixedSchedule = ql.Schedule(
    settlementDate,
    maturity,
    fixedLegTenor,
    calendar,
    fixedLegAdjustment,
    fixedLegAdjustment,
    ql.DateGeneration.Forward,
    False,
)
floatingSchedule = ql.Schedule(
    settlementDate,
    maturity,
    floatingLegTenor,
    calendar,
    floatingLegAdjustment,
    floatingLegAdjustment,
    ql.DateGeneration.Forward,
    False,
)

spot = ql.VanillaSwap(
    ql.VanillaSwap.Payer,
    nominal,
    fixedSchedule,
    fixedRate,
    fixedLegDayCounter,
    floatingSchedule,
    index,
    spread,
    floatingLegDayCounter,
)
spot.setPricingEngine(swapEngine)

forwardStart = calendar.advance(settlementDate, 1, ql.Years)
forwardEnd = calendar.advance(forwardStart, length, ql.Years)
fixedSchedule = ql.Schedule(
    forwardStart,
    forwardEnd,
    fixedLegTenor,
    calendar,
    fixedLegAdjustment,
    fixedLegAdjustment,
    ql.DateGeneration.Forward,
    False,
)
floatingSchedule = ql.Schedule(
    forwardStart,
    forwardEnd,
    floatingLegTenor,
    calendar,
    floatingLegAdjustment,
    floatingLegAdjustment,
    ql.DateGeneration.Forward,
    False,
)

forward = ql.VanillaSwap(
    ql.VanillaSwap.Payer,
    nominal,
    fixedSchedule,
    fixedRate,
    fixedLegDayCounter,
    floatingSchedule,
    index,
    spread,
    floatingLegDayCounter,
)
forward.setPricingEngine(swapEngine)

# price on the bootstrapped curves


def formatPrice(p, digits=2):
    format = "%%.%df" % digits
    return format % p


def formatRate(r, digits=2):
    format = "%%.%df %%%%" % digits
    return format % (r * 100)


headers = ("term structure", "net present value", "fair spread", "fair fixed rate")
separator = " | "

format = ""
width = 0
for h in headers[:-1]:
    format += "%%%ds" % len(h)
    format += separator
    width += len(h) + len(separator)
format += "%%%ds" % len(headers[-1])
width += len(headers[-1])

rule = "-" * width
dblrule = "=" * width
tab = " " * 8


def report(swap, name):
    print(format % (name, formatPrice(swap.NPV(), 2), formatRate(swap.fairSpread(), 4), formatRate(swap.fairRate(), 4)))


print(dblrule)
print("5-year market swap-rate = %s" % formatRate(swaps[(5, ql.Years)].value()))
print(dblrule)

# price on two different term structures

print(tab + "5-years swap paying %s" % formatRate(fixedRate))
print(separator.join(headers))
print(rule)

discountTermStructure.linkTo(depoFuturesSwapCurve)
forecastTermStructure.linkTo(depoFuturesSwapCurve)
report(spot, "depo-fut-swap")

discountTermStructure.linkTo(depoFraSwapCurve)
forecastTermStructure.linkTo(depoFraSwapCurve)
report(spot, "depo-FRA-swap")

print(rule)

# price the 1-year forward swap

print(tab + "5-years, 1-year forward swap paying %s" % formatRate(fixedRate))
print(rule)

discountTermStructure.linkTo(depoFuturesSwapCurve)
forecastTermStructure.linkTo(depoFuturesSwapCurve)
report(forward, "depo-fut-swap")

discountTermStructure.linkTo(depoFraSwapCurve)
forecastTermStructure.linkTo(depoFraSwapCurve)
report(forward, "depo-FRA-swap")

# modify the 5-years swap rate and reprice

swaps[(5, ql.Years)].setValue(0.046)

print(dblrule)
print("5-year market swap-rate = %s" % formatRate(swaps[(5, ql.Years)].value()))
print(dblrule)

print(tab + "5-years swap paying %s" % formatRate(fixedRate))
print(separator.join(headers))
print(rule)

discountTermStructure.linkTo(depoFuturesSwapCurve)
forecastTermStructure.linkTo(depoFuturesSwapCurve)
report(spot, "depo-fut-swap")

discountTermStructure.linkTo(depoFraSwapCurve)
forecastTermStructure.linkTo(depoFraSwapCurve)
report(spot, "depo-FRA-swap")

print(rule)

print(tab + "5-years, 1-year forward swap paying %s" % formatRate(fixedRate))
print(rule)

discountTermStructure.linkTo(depoFuturesSwapCurve)
forecastTermStructure.linkTo(depoFuturesSwapCurve)
report(forward, "depo-fut-swap")

discountTermStructure.linkTo(depoFraSwapCurve)
forecastTermStructure.linkTo(depoFraSwapCurve)
report(forward, "depo-FRA-swap")
