# Copyright (C) 2004, 2005, 2006, 2007 StatPro Italia srl
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

import QuantLib as ql

# global data
todaysDate = ql.Date(15, ql.May, 1998)
ql.Settings.instance().evaluationDate = todaysDate
settlementDate = ql.Date(17, ql.May, 1998)
riskFreeRate = ql.FlatForward(settlementDate, 0.05, ql.Actual365Fixed())

# option parameters
exercise = ql.EuropeanExercise(ql.Date(17, ql.May, 1999))
payoff = ql.PlainVanillaPayoff(ql.Option.Call, 8.0)

# market data
underlying = ql.SimpleQuote(7.0)
volatility = ql.BlackConstantVol(settlementDate, ql.TARGET(), 0.10, ql.Actual365Fixed())
dividendYield = ql.FlatForward(settlementDate, 0.05, ql.Actual365Fixed())

# report
header = " |".join(["%17s" % tag for tag in ["method", "value", "estimated error", "actual error"]])
print("")
print(header)
print("-" * len(header))

refValue = None


def report(method, x, dx=None):
    e = "%.4f" % abs(x - refValue)
    x = "%.5f" % x
    if dx:
        dx = "%.4f" % dx
    else:
        dx = "n/a"
    print(" |".join(["%17s" % y for y in [method, x, dx, e]]))


# good to go

process = ql.BlackScholesMertonProcess(
    ql.QuoteHandle(underlying),
    ql.YieldTermStructureHandle(dividendYield),
    ql.YieldTermStructureHandle(riskFreeRate),
    ql.BlackVolTermStructureHandle(volatility),
)

hestonProcess = ql.HestonProcess(
    ql.YieldTermStructureHandle(riskFreeRate),
    ql.YieldTermStructureHandle(dividendYield),
    ql.QuoteHandle(underlying),
    0.1 * 0.1,
    1.0,
    0.1 * 0.1,
    0.0001,
    0.0,
)
hestonModel = ql.HestonModel(hestonProcess)

option = ql.VanillaOption(payoff, exercise)

# method: analytic
option.setPricingEngine(ql.AnalyticEuropeanEngine(process))
value = option.NPV()
refValue = value
report("analytic", value)

# method: Heston semi-analytic
option.setPricingEngine(ql.AnalyticHestonEngine(hestonModel))
report("Heston analytic", option.NPV())

# method: Heston COS method
option.setPricingEngine(ql.COSHestonEngine(hestonModel))
report("Heston COS Method", option.NPV())

# method: Heston Exponentially-Fitted Gauss-Laguerre quadrature rule
option.setPricingEngine(ql.ExponentialFittingHestonEngine(hestonModel))
report("Heston ExpFitting", option.NPV())

# method: integral
option.setPricingEngine(ql.IntegralEngine(process))
report("integral", option.NPV())

# method: finite differences
timeSteps = 801
gridPoints = 800

option.setPricingEngine(ql.FdBlackScholesVanillaEngine(process, timeSteps, gridPoints))
report("finite diff.", option.NPV())

# method: binomial
timeSteps = 801

option.setPricingEngine(ql.BinomialVanillaEngine(process, "JR", timeSteps))
report("binomial (JR)", option.NPV())

option.setPricingEngine(ql.BinomialVanillaEngine(process, "CRR", timeSteps))
report("binomial (CRR)", option.NPV())

option.setPricingEngine(ql.BinomialVanillaEngine(process, "EQP", timeSteps))
report("binomial (EQP)", option.NPV())

option.setPricingEngine(ql.BinomialVanillaEngine(process, "Trigeorgis", timeSteps))
report("bin. (Trigeorgis)", option.NPV())

option.setPricingEngine(ql.BinomialVanillaEngine(process, "Tian", timeSteps))
report("binomial (Tian)", option.NPV())

option.setPricingEngine(ql.BinomialVanillaEngine(process, "LR", timeSteps))
report("binomial (LR)", option.NPV())

option.setPricingEngine(ql.BinomialVanillaEngine(process, "Joshi4", timeSteps))
report("binomial (Joshi)", option.NPV())

# method: finite differences
# not yet implemented

# method: Monte Carlo
option.setPricingEngine(ql.MCEuropeanEngine(process, "pseudorandom", timeSteps=1, requiredTolerance=0.02, seed=42))
report("MC (crude)", option.NPV(), option.errorEstimate())

option.setPricingEngine(ql.MCEuropeanEngine(process, "lowdiscrepancy", timeSteps=1, requiredSamples=32768))
report("MC (Sobol)", option.NPV())
