# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # Caps and Floors
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

calcDate = ql.Date(14, 6, 2016)
ql.Settings.instance().evaluationDate = calcDate

dates = [ql.Date(14,6,2016), ql.Date(14,9,2016),
         ql.Date(14,12,2016), ql.Date(14,6,2017),
         ql.Date(14,6,2019), ql.Date(14,6,2021),
         ql.Date(15,6,2026), ql.Date(16,6,2031),
         ql.Date(16,6,2036), ql.Date(14,6,2046)]

yields = [0.000000, 0.006616, 0.007049, 0.007795,
          0.009599, 0.011203, 0.015068, 0.017583,
          0.018998, 0.020080]

dayCount = ql.ActualActual(ql.ActualActual.Bond)
calendar = ql.UnitedStates(ql.UnitedStates.GovernmentBond)
interpolation = ql.Linear()
compounding = ql.Compounded
compoundingFrequency = ql.Annual
term_structure = ql.ZeroCurve(dates, yields, dayCount, calendar, interpolation, compounding, compoundingFrequency)
ts_handle = ql.YieldTermStructureHandle(term_structure)

start_date = ql.Date(14, 6, 2016)
end_date = ql.Date(14, 6 , 2026)
period = ql.Period(3, ql.Months)
buss_convention = ql.ModifiedFollowing
rule = ql.DateGeneration.Forward
end_of_month = False
schedule = ql.Schedule(start_date, end_date, period, calendar, buss_convention, buss_convention, rule, end_of_month)

iborIndex = ql.USDLibor(ql.Period(3, ql.Months), ts_handle)
iborIndex.addFixing(ql.Date(10,6,2016), 0.0065560)
ibor_leg = ql.IborLeg([1000000], schedule, iborIndex)

strike = 0.02
cap = ql.Cap(ibor_leg, [strike])
vols = ql.makeQuoteHandle(0.547295)
engine = ql.BlackCapFloorEngine(ts_handle, vols)
cap.setPricingEngine(engine)
print("Value of Caps given constant volatility:", cap.NPV())
