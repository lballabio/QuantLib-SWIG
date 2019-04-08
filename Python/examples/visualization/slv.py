# Copyright (C) 2019 Klaus Spanderen
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

import math

from QuantLib import *

from mpl_toolkits import mplot3d
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

todaysDate = Date(15,May,2019)
exerciseDate = todaysDate + Period(4, Years)
Settings.instance().evaluationDate = todaysDate
settlementDate = todaysDate + Period(2, Days)
dc = Actual365Fixed()
riskFreeRate = YieldTermStructureHandle(FlatForward(settlementDate, 0.05, dc))
dividendYield = YieldTermStructureHandle(FlatForward(settlementDate, 0.025, dc))

spot = 100
underlying = QuoteHandle(SimpleQuote(spot))

vol = 0.30

localVol = LocalVolSurface(
    BlackVolTermStructureHandle(
        BlackConstantVol(settlementDate, TARGET(), vol, dc)), 
    riskFreeRate, dividendYield, underlying)

hestonProcess = HestonProcess(riskFreeRate, dividendYield, underlying, 
                              0.09, 1.0, 0.06, 0.4, -0.75)

hestonModel = HestonModel(hestonProcess)

leverageFct = HestonSLVMCModel(
    localVol, hestonModel, MTBrownianGeneratorFactory(1234),
    exerciseDate, 91).leverageFunction()


#uncomment this if you want to use PDE calibration instead of Monte-Carlo
#
#fdmParams = HestonSLVFokkerPlanckFdmParams(
#    201, 101, 200, 30, 2.0, 0, 2, 0.01, 1e-4, 10000,
#    1e-5, 1e-5, 0.0000025, 1.0, 0.1, 0.9, 1e-4,
#    FdmHestonGreensFct.Gaussian, FdmSquareRootFwdOp.Log,
#    FdmSchemeDesc.Hundsdorfer()
#)
#leverageFct = HestonSLVFDMModel(
#    localVol, hestonModel, exerciseDate, fdmParams).leverageFunction()

tSteps = 40
uSteps = 30

tv = np.linspace(0.1, dc.yearFraction(settlementDate, exerciseDate), tSteps)

t = np.empty(tSteps*uSteps)
s = np.empty(tSteps*uSteps)
z = np.empty(tSteps*uSteps)

for i in range(0, tSteps): 
    scale = min(4, math.exp(3*math.sqrt(tv[i])*vol))
    sv = np.linspace(spot/scale, spot*scale, uSteps)    
    
    for j in range(0, uSteps):
        idx = i*uSteps + j
        t[idx] = tv[i]
        s[idx] = math.log(sv[j])
        z[idx] = leverageFct.localVol(t[idx], sv[j])


fig = plt.figure()
ax = plt.axes(projection='3d')

surf = ax.plot_trisurf(s, t, z, cmap=plt.cm.viridis, 
                       linewidth=0, antialiased=False, edgecolor='none')
ax.view_init(30, -120)

ax.set_xlabel("ln(S)")
ax.set_ylabel("Time")
ax.text2D(0.225, 0.985, "Leverage Function with $\eta=1.0$", transform=ax.transAxes)

fig.colorbar(surf, shrink=0.75, aspect=14)

plt.show()
