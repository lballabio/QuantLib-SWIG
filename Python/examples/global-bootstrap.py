# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Global curve bootstrap
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

# %% [markdown]
# ### Setup

# %%
today = ql.Date(26, 9, 2019)
spot = ql.TARGET().advance(today, 2, ql.Days)

# %%
ql.Settings.instance().evaluationDate = today

# %% [markdown]
# ### Data
#
# We'll use the following data as input:

# %%
refMktRates = [
    -0.373,
    -0.388,
    -0.402,
    -0.418,
    -0.431,
    -0.441,
    -0.45,
    -0.457,
    -0.463,
    -0.469,
    -0.461,
    -0.463,
    -0.479,
    -0.4511,
    -0.45418,
    -0.439,
    -0.4124,
    -0.37703,
    -0.3335,
    -0.28168,
    -0.22725,
    -0.1745,
    -0.12425,
    -0.07746,
    0.0385,
    0.1435,
    0.17525,
    0.17275,
    0.1515,
    0.1225,
    0.095,
    0.0644,
]

# %% [markdown]
# ### Market instruments

# %%
index = ql.Euribor6M()

# %% [markdown]
# The first market rate is for the 6-months deposit...

# %%
helpers = [
    ql.DepositRateHelper(
        refMktRates[0] / 100.0, ql.Period(6, ql.Months), 2, ql.TARGET(), ql.ModifiedFollowing, True, ql.Actual360()
    )
]

# %% [markdown]
# ...the next 12 are for FRAs...

# %%
helpers += [ql.FraRateHelper(r / 100.0, i + 1, index) for i, r in enumerate(refMktRates[1:13])]

# %% [markdown]
# ...and the others are swap rates.

# %%
swapTenors = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20, 25, 30, 35, 40, 45, 50]
helpers += [
    ql.SwapRateHelper(
        r / 100.0, ql.Period(T, ql.Years), ql.TARGET(), ql.Annual, ql.ModifiedFollowing, ql.Thirty360(ql.Thirty360.BondBasis), index
    )
    for r, T in zip(refMktRates[13:32], swapTenors)
]

# %% [markdown]
# We'll also add a few synthetic helpers:

# %%
additional_helpers = [ql.FraRateHelper(-0.004, 12 + i, index) for i in range(7)]
additional_dates = [ql.TARGET().advance(spot, 1 + i, ql.Months) for i in range(5)]


# %% [markdown]
# ### Global bootstrap
#
# This curve takes into account the market instruments, as well as the passed additional ones.

# %%
curve = ql.GlobalLinearSimpleZeroCurve(
    spot, helpers, ql.Actual365Fixed(), ql.GlobalBootstrap(additional_helpers, additional_dates, 1.0e-12)
)
curve.enableExtrapolation()


# %% [markdown]
# ### Report

# %%
data = []
for i, h in enumerate(helpers):
    pillar = h.pillarDate()

    if i < 13:
        day_counter = ql.Actual360()
        compounding = ql.Simple
    else:
        day_counter = ql.Thirty360(ql.Thirty360.BondBasis)
        compounding = ql.SimpleThenCompounded

    r = curve.zeroRate(pillar, day_counter, compounding, ql.Annual).rate()
    data.append((pillar.to_date(), r * 100))

# %%
df = pd.DataFrame(data, columns=["pillar", "zero rate"])
if not interactive:
    print(df)
df
