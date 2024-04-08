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
# # Leverage function for the Heston SLV model
#
# Copyright (&copy;) 2019 Klaus Spanderen
#
# This file is part of QuantLib, a free-software/open-source library for financial quantitative analysts and developers - https://www.quantlib.org/
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

# %% [markdown]
# This notebook only works with Python 3, at least on Travis.

# %%
import sys

if sys.version_info.major < 3:
    sys.exit()

# %%
import QuantLib as ql
from matplotlib import pyplot as plt
import numpy as np
import math

# %matplotlib inline

# %%
is_interactive = 'get_ipython' in globals()

# %%
todaysDate = ql.Date(15, ql.May, 2019)
ql.Settings.instance().evaluationDate = todaysDate

# %%
settlementDate = todaysDate + ql.Period(2, ql.Days)
exerciseDate = todaysDate + ql.Period(4, ql.Years)

# %%
dc = ql.Actual365Fixed()

spot = 100
underlying = ql.makeQuoteHandle(spot)

riskFreeRate = ql.YieldTermStructureHandle(ql.FlatForward(settlementDate, 0.05, dc))
dividendYield = ql.YieldTermStructureHandle(ql.FlatForward(settlementDate, 0.025, dc))

vol = 0.30
blackVol = ql.BlackVolTermStructureHandle(ql.BlackConstantVol(settlementDate, ql.TARGET(), vol, dc))

# %%
localVol = ql.LocalVolSurface(
    blackVol,
    riskFreeRate,
    dividendYield,
    underlying,
)

hestonProcess = ql.HestonProcess(riskFreeRate, dividendYield, underlying, 0.09, 1.0, 0.06, 0.4, -0.75)

hestonModel = ql.HestonModel(hestonProcess)

# %%
leverageFct = ql.HestonSLVMCModel(
    localVol, hestonModel, ql.MTBrownianGeneratorFactory(1234), exerciseDate, 91
).leverageFunction()

# %%
tSteps = 40
uSteps = 30

tv = np.linspace(0.1, dc.yearFraction(settlementDate, exerciseDate), tSteps)

t = np.empty(tSteps * uSteps)
s = np.empty(tSteps * uSteps)
z = np.empty(tSteps * uSteps)

for i in range(0, tSteps):
    scale = min(4, math.exp(3 * math.sqrt(tv[i]) * vol))
    sv = np.linspace(spot / scale, spot * scale, uSteps)

    for j in range(0, uSteps):
        idx = i * uSteps + j
        t[idx] = tv[i]
        s[idx] = math.log(sv[j])
        z[idx] = leverageFct.localVol(t[idx], sv[j])

# %%
fig = plt.figure(figsize=(12,8))
ax = plt.axes(projection="3d")

surf = ax.plot_trisurf(s, t, z, cmap=plt.cm.viridis, linewidth=0, antialiased=False, edgecolor="none")
ax.view_init(30, -120)

ax.set_xlabel("ln(S)")
ax.set_ylabel("Time")
ax.text2D(0.225, 0.985, "Leverage Function with $\eta=1.0$", transform=ax.transAxes)

fig.colorbar(surf, shrink=0.75, aspect=14)

plt.show(block=False)

# %% [markdown]
# When this is run as a Python script (i.e., from Travis), we need to close the figure in order to terminate.

# %%
if not is_interactive:
    plt.pause(3)
    plt.close()
