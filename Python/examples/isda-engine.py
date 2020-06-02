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

# # ISDA CDS engine
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

interactive = 'get_ipython' in globals()

tradeDate = ql.Date(21,5,2009)
ql.Settings.instance().setEvaluationDate(tradeDate)

ql.IborCoupon.createAtParCoupons()

dep_tenors = [1,2,3,6,9,12]
dep_quotes = [0.003081,0.005525,0.007163,0.012413,0.014,0.015488]
isdaRateHelpers = [ql.DepositRateHelper(dep_quotes[i],
                                        dep_tenors[i]*ql.Period(ql.Monthly),
                                        2,ql.WeekendsOnly(),
                                        ql.ModifiedFollowing,
                                        False,ql.Actual360())
                   for i in range(len(dep_tenors))]

swap_tenors = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30]
swap_quotes = [0.011907,
               0.01699,
               0.021198,
               0.02444,
               0.026937,
               0.028967,
               0.030504,
               0.031719,
               0.03279,
               0.034535,
               0.036217,
               0.036981,
               0.037246,
               0.037605]

isda_ibor = ql.IborIndex('IsdaIbor',3*ql.Period(ql.Monthly),2,
                         ql.USDCurrency(),ql.WeekendsOnly(),
                         ql.ModifiedFollowing,False,ql.Actual360())
isdaRateHelpers = isdaRateHelpers + [
    ql.SwapRateHelper(swap_quotes[i],swap_tenors[i]*ql.Period(ql.Annual),
                      ql.WeekendsOnly(),ql.Semiannual,ql.ModifiedFollowing,
                      ql.Thirty360(),isda_ibor)
    for i in range(len(swap_tenors))]

discountCurve = ql.RelinkableYieldTermStructureHandle()
discountCurve.linkTo(ql.PiecewiseLogLinearDiscount(0,ql.WeekendsOnly(),
                                                   isdaRateHelpers,
                                                   ql.Actual365Fixed()))

probabilityCurve = ql.RelinkableDefaultProbabilityTermStructureHandle()

termDates = [ql.Date(20, 6, 2010),
             ql.Date(20, 6, 2011),
             ql.Date(20, 6, 2012),
             ql.Date(20, 6, 2016),
             ql.Date(20, 6, 2019)]

spreads = [0.001, 0.1]
recoveries = [0.2, 0.4]

markitValues = [97798.29358, #0.001
                97776.11889, #0.001
                -914971.5977, #0.1
                -894985.6298, #0.1
                186921.3594, #0.001
                186839.8148, #0.001
                -1646623.672, #0.1
                -1579803.626, #0.1
                274298.9203,
                274122.4725,
                -2279730.93,
                -2147972.527,
                592420.2297,
                591571.2294,
                -3993550.206,
                -3545843.418,
                797501.1422,
                795915.9787,
                -4702034.688,
                -4042340.999]

tolerance = 1.0e-6


l = 0;
distance = 0

# +
data = []
for termDate in termDates:
    for spread in spreads:
        for recovery in recoveries:

            cdsSchedule = ql.Schedule(tradeDate+1, termDate,
                                      3*ql.Period(ql.Monthly),
                                      ql.WeekendsOnly(),
                                      ql.Following, ql.Unadjusted,
                                      ql.DateGeneration.CDS, False)

            quotedTrade = ql.CreditDefaultSwap(
                ql.Protection.Buyer,10000000,0,spread,cdsSchedule,
                ql.Following,ql.Actual360(),True,True,tradeDate+1,
                ql.WeekendsOnly().advance(tradeDate,3*ql.Period(ql.Daily)),
                ql.FaceValueClaim(), ql.Actual360(True))

            h = quotedTrade.impliedHazardRate(0,discountCurve,ql.Actual365Fixed(),
                                              recovery,1e-10,
                                              ql.CreditDefaultSwap.ISDA)

            probabilityCurve.linkTo(
                ql.FlatHazardRate(0,ql.WeekendsOnly(),
                                  ql.QuoteHandle(ql.SimpleQuote(h)),
                                  ql.Actual365Fixed()))

            engine = ql.IsdaCdsEngine(probabilityCurve,recovery,discountCurve)
            conventionalTrade = ql.CreditDefaultSwap(
                ql.Protection.Buyer,10000000,0,0.01,cdsSchedule,
                ql.Following,ql.Actual360(),True,True,tradeDate+1,
                ql.WeekendsOnly().advance(tradeDate,3*ql.Period(ql.Daily)),
                ql.FaceValueClaim(), ql.Actual360(True))
            conventionalTrade.setPricingEngine(engine)

            upfront = conventionalTrade.notional() * conventionalTrade.fairUpfront()

            data.append(
                (h,
                 upfront,
                 abs(upfront-markitValues[l]),
                 abs(upfront-markitValues[l])<tolerance)
            )
            distance = distance + abs(upfront-markitValues[l])

            l = l + 1

df = pd.DataFrame(data, columns=["Hazard", "Upfront", "Distance", "Within tolerance"])
if not interactive:
    print(df)
df
# -

print('total distance:',distance)
