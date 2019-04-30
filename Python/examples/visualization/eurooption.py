# Example of option baskets
# Distributed under BSD license

import QuantLib as ql

todaysDate = ql.Date(15, ql.May, 1998)
ql.Settings.instance().evaluationDate = todaysDate
settlementDate = ql.Date(17, ql.May, 1998)
riskFreeQuote = ql.SimpleQuote(0.05)
riskFreeRate = ql.FlatForward(settlementDate, ql.QuoteHandle(riskFreeQuote), ql.Actual365Fixed())

# option parameters
exercise1 = ql.AmericanExercise(settlementDate, ql.Date(17, ql.May, 1999))
exercise2 = ql.EuropeanExercise(settlementDate)
payoff = ql.PlainVanillaPayoff(ql.Option.Call, 40.0)

# market data
underlying = ql.SimpleQuote(36.0)
volatilityQuote = ql.SimpleQuote(0.05)
volatility = ql.BlackConstantVol(todaysDate, ql.QuoteHandle(volatilityQuote), ql.Actual365Fixed())
dividendYield = ql.FlatForward(settlementDate, 0.00, ql.Actual365Fixed())

# good to go
process = ql.BlackScholesMertonProcess(
    ql.QuoteHandle(underlying),
    ql.YieldTermStructureHandle(dividendYield),
    ql.YieldTermStructureHandle(riskFreeRate),
    ql.BlackVolTermStructureHandle(volatility),
)
option1 = ql.VanillaOption(process, payoff, exercise1)
option1.setPricingEngine(ql.BaroneAdesiWhaleyEngine())


def f(x, y):
    underlying.setValue(x)
    volatilityQuote.setValue(y)
    return option1.NPV()


def setQuote(r):
    riskFreeQuote.setValue(r)
