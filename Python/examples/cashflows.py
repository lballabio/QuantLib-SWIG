# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Cash-flow analysis
#
# Copyright (&copy;) 2020 StatPro Italia srl
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

# %%
import QuantLib as ql
import pandas as pd

# %%
interactive = "get_ipython" in globals()

# %%
today = ql.Date(19, ql.October, 2020)
ql.Settings.instance().evaluationDate = today

# %% [markdown]
# ### Term structure construction

# %%
dates = [
    ql.Date(19,10,2020),
    ql.Date(19,11,2020),
    ql.Date(19, 1,2021),
    ql.Date(19, 4,2021),
    ql.Date(19,10,2021),
    ql.Date(19, 4,2022),
    ql.Date(19,10,2022),
    ql.Date(19,10,2023),
    ql.Date(19,10,2025),
    ql.Date(19,10,2030),
    ql.Date(19,10,2035),
    ql.Date(19,10,2040),
]

rates = [
    -0.004,
    -0.002,
    0.001,
    0.005,
    0.009,
    0.010,
    0.010,
    0.012,
    0.017,
    0.019,
    0.028,
    0.032,
]

forecast_curve = ql.ZeroCurve(dates, rates, ql.Actual365Fixed())

# %%
forecast_handle = ql.YieldTermStructureHandle(forecast_curve)

# %% [markdown]
# ### Swap construction
#
# We'll use an overnight swap as an example.  We're keeping the initialization simple, but the analysis work in the same way for more complex ones, as well as for other kinds of swaps and bonds (once we extract the cashflows from them using the proper methods).

# %%
swap = ql.MakeOIS(swapTenor=ql.Period(5, ql.Years),
                  overnightIndex=ql.Eonia(forecast_handle),
                  fixedRate=0.002)

# %% [markdown]
# ### Cash-flow analysis
#
# The fixed-rate coupons can be extracted from the swap using the `fixedLeg` method.  They are returned as instances of the base `Cashflow` class, so the only methods we have directly available are from that class interface:

# %%
fixed_leg = swap.fixedLeg()

# %%
df = pd.DataFrame([(c.date(), c.amount()) for c in fixed_leg if c.date() > today],
                  columns=['date', 'amount'])
df

# %% [markdown]
# The following displays the results when this is run as a Python script (in which case the cell above is not displayed).

# %%
if not interactive:
    print(df)

# %% [markdown]
# If we want to extract more information, we need to upcast the coupons to a more specific class.  This can be done by using the `as_fixed_rate_coupon` method.  In this case, the upcast works by construction; but in the general case we might have cashflows for which the upcast fails (e.g., the redemption for a bond) so we have to check for nulls.

# %%
coupons = []
for cf in fixed_leg:
    c = ql.as_fixed_rate_coupon(cf)
    if c:
        coupons.append(c)

# %% [markdown]
# We can now access methods from the coupon class.

# %%
df = pd.DataFrame([(c.date(), c.amount(), c.rate(), c.accrualStartDate(), c.accrualEndDate(), c.accrualPeriod())
                   for c in coupons if c.date() > today],
                  columns=['payment date', 'amount', 'rate', 'start date', 'end date', 'accrual period'])
df

# %%
if not interactive:
    print(df)

# %% [markdown]
# The same goes for the floating leg: in this case, we need to upcast to floating-rate coupons in order to access the specific methods we'll need.

# %%
floating_leg = swap.overnightLeg()

# %%
coupons = []
for cf in floating_leg:
    c = ql.as_floating_rate_coupon(cf)
    if c:
        coupons.append(c)

# %%
df = pd.DataFrame([(c.date(), c.amount(), c.rate(), c.accrualStartDate(), c.accrualEndDate(), c.accrualPeriod())
                   for c in coupons if c.date() > today],
                  columns=['payment date', 'amount', 'rate', 'start date', 'end date', 'accrual period'])
df

# %%
if not interactive:
    print(df)
