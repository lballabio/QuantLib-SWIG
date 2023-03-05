# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Callable bonds (Calls)
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
import numpy as np
calcDate = ql.Date(16, 8, 2006) 
ql.Settings.instance().evaluationDate = calcDate

dayCount = ql.ActualActual(ql.ActualActual.Bond)
rate = 0.0465
termStructure = ql.FlatForward(calcDate, rate, dayCount, ql.Compounded, ql.Semiannual)
term_Structure_Handle = ql.RelinkableYieldTermStructureHandle(termStructure)

callabilitySchedule = ql.CallabilitySchedule()
callPrice = 100.0
callDate = ql.Date(15, ql.September, 2006); 
nc = ql.NullCalendar()

# Number of calldates is 24
for i in range(0, 24):
    callabilityPrice  = ql.BondPrice(callPrice, ql.BondPrice.Clean)
    callabilitySchedule.append(ql.Callability(callabilityPrice, ql.Callability.Call, callDate))
    callDate = nc.advance(callDate, 3, ql.Months)

issueDate = ql.Date(16, ql.September, 2004)
maturityDate = ql.Date(15, ql.September, 2012)
calendar = ql.UnitedStates(ql.UnitedStates.GovernmentBond)
tenor = ql.Period(ql.Quarterly)
accrualConvention = ql.Unadjusted
schedule = ql.Schedule(issueDate, maturityDate, tenor, calendar, accrualConvention, accrualConvention, ql.DateGeneration.Backward, False)

settlement_days = 3
faceAmount = 100
accrual_daycount = ql.ActualActual(ql.ActualActual.Bond)
coupon = 0.025
bond = ql.CallableFixedRateBond(settlement_days, faceAmount, schedule, [coupon], accrual_daycount, ql.Following, faceAmount, issueDate, callabilitySchedule)

maxIterations = 1000
accuracy = 1e-8
gridIntervals = 40
reversionParameter = .03

def value_bond(a, s, grid_points, bond): 
    model = ql.HullWhite(term_Structure_Handle, a, s)
    engine = ql.TreeCallableFixedRateBondEngine(model, grid_points) 
    bond.setPricingEngine(engine)
    return bond

def value_bond2(reversionParameter, s, gridIntervals, bond): 
    model2 = ql.HullWhite(term_Structure_Handle, reversionParameter, s)
    engine2 = ql.TreeCallableFixedRateBondEngine(model2, gridIntervals) 
    bond.setPricingEngine(engine2)
    return bond

# 6% mean reversion and 20% volatility
value_bond2(0.06, 0.20, 40, bond) 
print("Bond price using clean price: ", bond.cleanPrice())
print("Bond price using net present value: ", bond.NPV())

value_bond(0.03, 0.15, 40, bond) 
print("Bond price using clean price: ", bond.cleanPrice())
print("Bond price using net present value: ", bond.NPV())
