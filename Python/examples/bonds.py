# Copyright (C) 2008 Florent Grenier
# Copyright (C) 2010 Lluis Pujol Bajador
#
# This file is part of QuantLib, a free-software/open-source library
# for financial quantitative analysts and developers - http://quantlib.org/
#
# QuantLib is free software: you can redistribute it and/or modify it
# under the terms of the QuantLib license.  You should have received a
# copy of the license along with this program; if not, please email
# <quantlib-dev@lists.sf.net>. The license is also available online at
# <http://quantlib.org/license.shtml>.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the license for more details.

#    This example shows how to set up a term structure and then price
#    some simple bonds. The last part is dedicated to peripherical
#    computations such as "Yield to Price" or "Price to Yield"

import QuantLib as ql

# global data
calendar = ql.TARGET()
settlementDate = ql.Date(18, ql.September, 2008)
settlementDate = calendar.adjust(settlementDate)

fixingDays = 3
settlementDays = 3

todaysDate = calendar.advance(settlementDate, -fixingDays, ql.Days)
ql.Settings.instance().evaluationDate = todaysDate

print("Today: " + str(todaysDate))
print("Settlement Date: " + str(settlementDate))

# market quotes

# constructing bond yield curve

zcQuotes = [(0.0096, ql.Period(3, ql.Months)), (0.0145, ql.Period(6, ql.Months)), (0.0194, ql.Period(1, ql.Years))]

zcBondsDayCounter = ql.Actual365Fixed()

zcHelpers = [
    ql.DepositRateHelper(
        ql.QuoteHandle(ql.SimpleQuote(r)), tenor, fixingDays, calendar, ql.ModifiedFollowing, True, zcBondsDayCounter
    )
    for (r, tenor) in zcQuotes
]

# setup bonds

redemption = 100.0
numberOfBonds = 5

bondQuotes = [
    (ql.Date(15, ql.March, 2005), ql.Date(31, ql.August, 2010), 0.02375, 100.390625),
    (ql.Date(15, ql.June, 2005), ql.Date(31, ql.August, 2011), 0.04625, 106.21875),
    (ql.Date(30, ql.June, 2006), ql.Date(31, ql.August, 2013), 0.03125, 100.59375),
    (ql.Date(15, ql.November, 2002), ql.Date(15, ql.August, 2018), 0.04000, 101.6875),
    (ql.Date(15, ql.May, 1987), ql.Date(15, ql.May, 2038), 0.04500, 102.140625),
]

# Definition of the rate helpers

bondsHelpers = []

for issueDate, maturity, couponRate, marketQuote in bondQuotes:
    schedule = ql.Schedule(
        issueDate,
        maturity,
        ql.Period(ql.Semiannual),
        ql.UnitedStates(ql.UnitedStates.GovernmentBond),
        ql.Unadjusted,
        ql.Unadjusted,
        ql.DateGeneration.Backward,
        False,
    )
    bondsHelpers.append(
        ql.FixedRateBondHelper(
            ql.QuoteHandle(ql.SimpleQuote(marketQuote)),
            settlementDays,
            100.0,
            schedule,
            [couponRate],
            ql.ActualActual(ql.ActualActual.Bond),
            ql.Unadjusted,
            redemption,
            issueDate,
        )
    )

###################################
####  **  CURVE BUILDING  **  #####
###################################

termStructureDayCounter = ql.ActualActual(ql.ActualActual.ISDA)

# not needed as defined in the interface file:  tolerance = 1.0e-15

bondInstruments = zcHelpers + bondsHelpers

bondDiscountingTermStructure = ql.PiecewiseFlatForward(settlementDate, bondInstruments, termStructureDayCounter)

# Building of the Libor forecasting curve
# deposits
dQuotes = [
    (0.043375, ql.Period(1, ql.Weeks)),
    (0.031875, ql.Period(1, ql.Months)),
    (0.0320375, ql.Period(3, ql.Months)),
    (0.03385, ql.Period(6, ql.Months)),
    (0.0338125, ql.Period(9, ql.Months)),
    (0.0335125, ql.Period(1, ql.Years)),
]
sQuotes = [
    (0.0295, ql.Period(2, ql.Years)),
    (0.0323, ql.Period(3, ql.Years)),
    (0.0359, ql.Period(5, ql.Years)),
    (0.0412, ql.Period(10, ql.Years)),
    (0.0433, ql.Period(15, ql.Years)),
]

# deposits

depositDayCounter = ql.Actual360()
depositHelpers = [
    ql.DepositRateHelper(
        ql.QuoteHandle(ql.SimpleQuote(rate)), tenor, fixingDays, calendar, ql.ModifiedFollowing, True, depositDayCounter
    )
    for rate, tenor in dQuotes
]

# swaps

swFixedLegFrequency = ql.Annual
swFixedLegConvention = ql.Unadjusted
swFixedLegDayCounter = ql.Thirty360(ql.Thirty360.European)
swFloatingLegIndex = ql.Euribor6M()
forwardStart = ql.Period(1, ql.Days)
swapHelpers = [
    ql.SwapRateHelper(
        ql.QuoteHandle(ql.SimpleQuote(rate)),
        tenor,
        calendar,
        swFixedLegFrequency,
        swFixedLegConvention,
        swFixedLegDayCounter,
        swFloatingLegIndex,
        ql.QuoteHandle(),
        forwardStart,
    )
    for rate, tenor in sQuotes
]

depoSwapInstruments = depositHelpers + swapHelpers

depoSwapTermStructure = ql.PiecewiseFlatForward(settlementDate, depoSwapInstruments, termStructureDayCounter)

# Term structures that will be used for pricing:
# the one used for discounting cash flows

discountingTermStructure = ql.RelinkableYieldTermStructureHandle()

# the one used for forward rate forecasting

forecastingTermStructure = ql.RelinkableYieldTermStructureHandle()

#######################################
#        BONDS TO BE PRICED           #
#######################################

# common data

faceAmount = 100

# pricing engine
bondEngine = ql.DiscountingBondEngine(discountingTermStructure)

# zero coupon bond

zeroCouponBond = ql.ZeroCouponBond(
    settlementDays,
    ql.UnitedStates(ql.UnitedStates.GovernmentBond),
    faceAmount,
    ql.Date(15, ql.August, 2013),
    ql.Following,
    116.92,
    ql.Date(15, ql.August, 2003),
)

zeroCouponBond.setPricingEngine(bondEngine)

# fixed 4.5% US Treasury note

fixedBondSchedule = ql.Schedule(
    ql.Date(15, ql.May, 2007),
    ql.Date(15, ql.May, 2017),
    ql.Period(ql.Semiannual),
    ql.UnitedStates(ql.UnitedStates.GovernmentBond),
    ql.Unadjusted,
    ql.Unadjusted,
    ql.DateGeneration.Backward,
    False,
)

fixedRateBond = ql.FixedRateBond(
    settlementDays,
    faceAmount,
    fixedBondSchedule,
    [0.045],
    ql.ActualActual(ql.ActualActual.Bond),
    ql.ModifiedFollowing,
    100.0,
    ql.Date(15, ql.May, 2007),
)

fixedRateBond.setPricingEngine(bondEngine)

# Floating rate bond (3M USD Libor + 0.1%)
# Should and will be priced on another curve later...

liborTermStructure = ql.RelinkableYieldTermStructureHandle()

libor3m = ql.USDLibor(ql.Period(3, ql.Months), liborTermStructure)
libor3m.addFixing(ql.Date(17, ql.July, 2008), 0.0278625)

floatingBondSchedule = ql.Schedule(
    ql.Date(21, ql.October, 2005),
    ql.Date(21, ql.October, 2010),
    ql.Period(ql.Quarterly),
    ql.UnitedStates(ql.UnitedStates.NYSE),
    ql.Unadjusted,
    ql.Unadjusted,
    ql.DateGeneration.Backward,
    True,
)

floatingRateBond = ql.FloatingRateBond(
    settlementDays,
    faceAmount,
    floatingBondSchedule,
    libor3m,
    ql.Actual360(),
    ql.ModifiedFollowing,
    spreads=[0.001],
    inArrears=True,
    issueDate=ql.Date(21, ql.October, 2005),
)

floatingRateBond.setPricingEngine(bondEngine)

# coupon pricers

pricer = ql.BlackIborCouponPricer()

# optionlet volatilities
volatility = 0.0
vol = ql.ConstantOptionletVolatility(settlementDays, calendar, ql.ModifiedFollowing, volatility, ql.Actual365Fixed())

pricer.setCapletVolatility(ql.OptionletVolatilityStructureHandle(vol))
ql.setCouponPricer(floatingRateBond.cashflows(), pricer)

# Yield curve bootstrapping
forecastingTermStructure.linkTo(depoSwapTermStructure)
discountingTermStructure.linkTo(bondDiscountingTermStructure)

# We are using the depo & swap curve to estimate the future Libor rates
liborTermStructure.linkTo(depoSwapTermStructure)

#############################
#       BOND PRICING        #
#############################

# write column headings
def formatPrice(p, digits=2):
    format = "%%.%df" % digits
    return format % p


def formatRate(r, digits=2):
    format = "%%.%df %%%%" % digits
    return format % (r * 100)


def report(Info, Zc, Fix, Frn, format):
    if format == "Price":
        Zc = formatPrice(Zc)
        Fix = formatPrice(Fix)
        Frn = formatPrice(Frn)
    else:
        if Info.find("coupon") == -1:
            Zc = formatRate(Zc)
        else:
            Zc = "N/A"
        Fix = formatRate(Fix)
        Frn = formatRate(Frn)

    print("%19s" % Info + " |" + " |".join(["%10s" % y for y in [Zc, Fix, Frn]]))


headers = ["ZC", "Fixed", "Floating"]
print("")
print("%19s" % "" + " |" + " |".join(["%10s" % y for y in headers]))

separator = " | "
widths = [18, 10, 10, 10]
width = widths[0] + widths[1] + widths[2] + widths[3] + widths[3]
rule = "-" * width
dblrule = "=" * width
tab = " " * 8

print(rule)
report("Net present value", zeroCouponBond.NPV(), fixedRateBond.NPV(), floatingRateBond.NPV(), "Price")
report("Clean price", zeroCouponBond.cleanPrice(), fixedRateBond.cleanPrice(), floatingRateBond.cleanPrice(), "Price")
report("Dirty price", zeroCouponBond.dirtyPrice(), fixedRateBond.dirtyPrice(), floatingRateBond.dirtyPrice(), "Price")
report(
    "Accrued coupon",
    zeroCouponBond.accruedAmount(),
    fixedRateBond.accruedAmount(),
    floatingRateBond.accruedAmount(),
    "Price",
)
report("Previous coupon", 0, fixedRateBond.previousCouponRate(), floatingRateBond.previousCouponRate(), "Rate")
report("Next coupon", 0, fixedRateBond.nextCouponRate(), floatingRateBond.nextCouponRate(), "Rate")
report(
    "Yield",
    zeroCouponBond.bondYield(ql.Actual360(), ql.Compounded, ql.Annual),
    fixedRateBond.bondYield(ql.Actual360(), ql.Compounded, ql.Annual),
    floatingRateBond.bondYield(ql.Actual360(), ql.Compounded, ql.Annual),
    "Rate",
)
print("")

# Other computations

print("Sample indirect computations (for the floating rate bond): ")
print(rule)
print(
    "Yield to Clean Price: "
    + formatPrice(
        floatingRateBond.cleanPrice(
            floatingRateBond.bondYield(ql.Actual360(), ql.Compounded, ql.Annual),
            ql.Actual360(),
            ql.Compounded,
            ql.Annual,
            settlementDate,
        )
    )
)
print(
    "Clean Price to Yield: "
    + formatRate(
        floatingRateBond.bondYield(
            floatingRateBond.cleanPrice(), ql.Actual360(), ql.Compounded, ql.Annual, settlementDate
        )
    )
)
