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
riskFreeRate = ql.FlatForward(settlementDate, 0.06, ql.Actual365Fixed())

# option parameters
exercise = ql.AmericanExercise(settlementDate, ql.Date(17, ql.May, 1999))
payoff = ql.PlainVanillaPayoff(ql.Option.Put, 40.0)

# market data
underlying = ql.SimpleQuote(36.0)
volatility = ql.BlackConstantVol(todaysDate, ql.TARGET(), 0.20, ql.Actual365Fixed())
dividendYield = ql.FlatForward(settlementDate, 0.00, ql.Actual365Fixed())

# report
header = "%19s" % "method" + " |" + " |".join(["%17s" % tag for tag in ["value", "estimated error", "actual error"]])
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
    print("%19s" % method + " |" + " |".join(["%17s" % y for y in [x, dx, e]]))


# good to go

process = ql.BlackScholesMertonProcess(
    ql.QuoteHandle(underlying),
    ql.YieldTermStructureHandle(dividendYield),
    ql.YieldTermStructureHandle(riskFreeRate),
    ql.BlackVolTermStructureHandle(volatility),
)

option = ql.VanillaOption(payoff, exercise)

refValue = 4.48667344
report("reference value", refValue)

# method: analytic

option.setPricingEngine(ql.BaroneAdesiWhaleyEngine(process))
report("Barone-Adesi-Whaley", option.NPV())

option.setPricingEngine(ql.BjerksundStenslandEngine(process))
report("Bjerksund-Stensland", option.NPV())

# method: finite differences
timeSteps = 801
gridPoints = 800

option.setPricingEngine(ql.FdBlackScholesVanillaEngine(process, timeSteps, gridPoints))
report("finite differences", option.NPV())

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
