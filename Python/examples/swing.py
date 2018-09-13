# Copyright (C) 2018 Klaus Spanderen
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

from QuantLib import *
import math

todaysDate = Date(30,September,2018)
Settings.instance().evaluationDate = todaysDate
settlementDate = todaysDate
riskFreeRate = FlatForward(settlementDate, 0.0, Actual365Fixed())
dividendYield = FlatForward(settlementDate, 0.0, Actual365Fixed())
underlying = SimpleQuote(30.0)
volatility = BlackConstantVol(todaysDate, TARGET(), 0.20, Actual365Fixed())


exerciseDates = [ Date(1, January, 2019) + i
                  for i in range(31) ]

swingOption = VanillaSwingOption(VanillaForwardPayoff(Option.Call, underlying.value()),
                                 SwingExercise(exerciseDates), 0, len(exerciseDates))

bsProcess = BlackScholesMertonProcess(
    QuoteHandle(underlying),
    YieldTermStructureHandle(dividendYield),
    YieldTermStructureHandle(riskFreeRate),
    BlackVolTermStructureHandle(volatility))

swingOption.setPricingEngine(FdSimpleBSSwingEngine(bsProcess))

print("Black Scholes Price: %f" % swingOption.NPV())

x0 = 0.0
x1 = 0.0

beta = 4.0
eta = 4.0
jumpIntensity = 1.0
speed = 1.0
volatility = 0.1

curveShape = []
for d in exerciseDates:
    t = Actual365Fixed().yearFraction(todaysDate, d)
    gs = math.log(underlying.value())                                   \
        - volatility*volatility/(4*speed)*(1-math.exp(-2*speed*t))      \
        - jumpIntensity/beta*math.log((eta-math.exp(-beta*t))/(eta-1.0))
    curveShape.append((t, gs))

ouProcess = ExtendedOrnsteinUhlenbeckProcess(speed, volatility, x0, lambda x: x0)
jProcess = ExtOUWithJumpsProcess(ouProcess, x1, beta, jumpIntensity, eta)

swingOption.setPricingEngine(FdSimpleExtOUJumpSwingEngine(
    jProcess, riskFreeRate, 25, 25, 200, curveShape))

print("Kluge Model Price  : %f" % swingOption.NPV())
