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

# # Credit default swaps
#
# Copyright (&copy;) 2014 Thema Consulting SA
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

# ### Setup

import QuantLib as ql

calendar = ql.TARGET()

todaysDate = ql.Date(15, ql.May, 2007)
ql.Settings.instance().evaluationDate = todaysDate

risk_free_rate = ql.YieldTermStructureHandle(ql.FlatForward(todaysDate, 0.01, ql.Actual365Fixed()))

# ### CDS parameters

recovery_rate = 0.5
quoted_spreads = [0.0150, 0.0150, 0.0150, 0.0150]
tenors = [ql.Period(3, ql.Months), ql.Period(6, ql.Months), ql.Period(1, ql.Years), ql.Period(2, ql.Years)]
maturities = [calendar.adjust(todaysDate + x, ql.Following) for x in tenors]

instruments = [
    ql.SpreadCdsHelper(
        ql.makeQuoteHandle(s),
        tenor,
        0,
        calendar,
        ql.Quarterly,
        ql.Following,
        ql.DateGeneration.TwentiethIMM,
        ql.Actual365Fixed(),
        recovery_rate,
        risk_free_rate,
    )
    for s, tenor in zip(quoted_spreads, tenors)
]

hazard_curve = ql.PiecewiseFlatHazardRate(todaysDate, instruments, ql.Actual365Fixed())
print("Calibrated hazard rate values: ")
for x in hazard_curve.nodes():
    print("hazard rate on %s is %.7f" % x)

print("Some survival probability values: ")
print(
    "1Y survival probability: %.4g, \n\t\texpected %.4g"
    % (hazard_curve.survivalProbability(todaysDate + ql.Period("1Y")), 0.9704)
)
print(
    "2Y survival probability: %.4g, \n\t\texpected %.4g"
    % (hazard_curve.survivalProbability(todaysDate + ql.Period("2Y")), 0.9418)
)

# ### Reprice instruments

nominal = 1000000.0
probability = ql.DefaultProbabilityTermStructureHandle(hazard_curve)

# We'll create a cds for every maturity:

all_cds = []
for maturity, s in zip(maturities, quoted_spreads):
    schedule = ql.Schedule(
        todaysDate,
        maturity,
        ql.Period(ql.Quarterly),
        calendar,
        ql.Following,
        ql.Unadjusted,
        ql.DateGeneration.TwentiethIMM,
        False,
    )
    cds = ql.CreditDefaultSwap(ql.Protection.Seller, nominal, s, schedule, ql.Following, ql.Actual365Fixed())
    engine = ql.MidPointCdsEngine(probability, recovery_rate, risk_free_rate)
    cds.setPricingEngine(engine)
    all_cds.append(cds)

print("Repricing of quoted CDSs employed for calibration: ")
for cds, tenor in zip(all_cds, tenors):
    print("%s fair spread: %.7g" % (tenor, cds.fairSpread()))
    print("   NPV: %g" % cds.NPV())
    print("   default leg: %.7g" % cds.defaultLegNPV())
    print("   coupon leg: %.7g" % cds.couponLegNPV())
    print("")
