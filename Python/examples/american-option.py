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
# # American options
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

# %%
import QuantLib as ql
import pandas as pd

# %% [markdown]
# ### Global parameters

# %%
todaysDate = ql.Date(15, ql.May, 1998)
ql.Settings.instance().evaluationDate = todaysDate

# %%
interactive = "get_ipython" in globals()

# %% [markdown]
# ### Option construction

# %%
exercise = ql.AmericanExercise(todaysDate, ql.Date(17, ql.May, 1999))
payoff = ql.PlainVanillaPayoff(ql.Option.Put, 40.0)

# %%
option = ql.VanillaOption(payoff, exercise)

# %% [markdown]
# ### Market data

# %%
underlying = ql.SimpleQuote(36.0)
dividendYield = ql.FlatForward(todaysDate, 0.00, ql.Actual365Fixed())
volatility = ql.BlackConstantVol(todaysDate, ql.TARGET(), 0.20, ql.Actual365Fixed())
riskFreeRate = ql.FlatForward(todaysDate, 0.06, ql.Actual365Fixed())

# %%
process = ql.BlackScholesMertonProcess(
    ql.QuoteHandle(underlying),
    ql.YieldTermStructureHandle(dividendYield),
    ql.YieldTermStructureHandle(riskFreeRate),
    ql.BlackVolTermStructureHandle(volatility),
)

# %% [markdown]
# ### Pricing
#
# We'll collect tuples of method name, option value, and estimated error from the analytic formula.

# %%
results = []

# %% [markdown]
# #### Analytic approximations

# %%
option.setPricingEngine(ql.BaroneAdesiWhaleyApproximationEngine(process))
results.append(("Barone-Adesi-Whaley", option.NPV()))

# %%
option.setPricingEngine(ql.BjerksundStenslandApproximationEngine(process))
results.append(("Bjerksund-Stensland", option.NPV()))

# %% [markdown]
# #### Finite-difference method

# %%
timeSteps = 801
gridPoints = 800

# %%
option.setPricingEngine(ql.FdBlackScholesVanillaEngine(process, timeSteps, gridPoints))
results.append(("finite differences", option.NPV()))


# %% [markdown]
# #### Li, M. QD+ American engine

# %%
option.setPricingEngine(ql.QdPlusAmericanEngine(process))
results.append(("QD+", option.NPV()))


# %% [markdown]
# #### Leif Andersen, Mark Lake and Dimitri Offengenden high performance American engine

# %%
option.setPricingEngine(
    ql.QdFpAmericanEngine(process, ql.QdFpAmericanEngine.accurateScheme())
)
results.append(("QD+ fixed point", option.NPV()))


# %% [markdown]
# #### Binomial method

# %%
timeSteps = 801

# %%
for tree in ["JR", "CRR", "EQP", "Trigeorgis", "Tian", "LR", "Joshi4"]:
    option.setPricingEngine(ql.BinomialVanillaEngine(process, tree, timeSteps))
    results.append(("Binomial (%s)" % tree, option.NPV()))

# %% [markdown]
# ### Results

# %%
df = pd.DataFrame(results, columns=["Method", "Option value"])
df.style.hide(axis="index")

# %% [markdown]
# The following displays the results when this is run as a Python script (in which case the cell above is not displayed).

# %%
if not interactive:
    print(df)
