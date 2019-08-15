# Example of option baskets
# Distributed under BSD License

import QuantLib as ql


class BasketOptionClass:
    def __init__(self, btype):
        # global data
        self.todaysDate = ql.Date(15, ql.May, 1998)
        ql.Settings.instance().evaluationDate = self.todaysDate
        self.settlementDate = ql.Date(17, ql.May, 1998)
        self.riskFreeQuote = ql.SimpleQuote(0.05)
        self.riskFreeRate = ql.FlatForward(self.settlementDate, ql.QuoteHandle(self.riskFreeQuote), ql.Actual365Fixed())

        # option parameters
        self.exercise = ql.EuropeanExercise(ql.Date(17, ql.May, 1999))
        self.payoff = ql.PlainVanillaPayoff(ql.Option.Call, 8.0)

        # market data
        self.underlying1 = ql.SimpleQuote(10.0)
        self.volatility1 = ql.BlackConstantVol(self.todaysDate, 0.20, ql.Actual365Fixed())
        self.dividendYield1 = ql.FlatForward(self.settlementDate, 0.05, ql.Actual365Fixed())
        self.underlying2 = ql.SimpleQuote(7.0)
        self.volatility2 = ql.BlackConstantVol(self.todaysDate, 0.10, ql.Actual365Fixed())
        self.dividendYield2 = ql.FlatForward(self.settlementDate, 0.05, ql.Actual365Fixed())

        self.process1 = ql.BlackScholesMertonProcess(
            ql.QuoteHandle(self.underlying1),
            ql.YieldTermStructureHandle(self.dividendYield1),
            ql.YieldTermStructureHandle(self.riskFreeRate),
            ql.BlackVolTermStructureHandle(self.volatility1),
        )

        self.process2 = ql.BlackScholesMertonProcess(
            ql.QuoteHandle(self.underlying2),
            ql.YieldTermStructureHandle(self.dividendYield2),
            ql.YieldTermStructureHandle(self.riskFreeRate),
            ql.BlackVolTermStructureHandle(self.volatility2),
        )

        self.procs = ql.StochasticProcessVector()
        self.procs.push_back(self.process1)
        self.procs.push_back(self.process2)

        self.matrix = ql.Matrix(2, 2)
        self.matrix[0][0] = 1.0
        self.matrix[1][1] = 1.0
        self.matrix[0][1] = 0.25
        self.matrix[1][0] = 0.25

        self.process = ql.StochasticProcessArray(self.procs, self.matrix)
        self.engine = ql.MCBasketEngine("lowdiscrepancy", timeStepsPerYear=1, requiredTolerance=0.02, seed=42)
        if btype == "min":
            bpayoff = ql.MinBasketPayoff(self.payoff)
        elif btype == "max":
            bpayoff = ql.MaxBasketPayoff(self.payoff)
        else:
            bpayoff = ql.AverageBasketPayoff(self.payoff, 2)

        self.option = ql.BasketOption(self.process, bpayoff, self.exercise, self.engine)

    def npv(self, x, y):
        self.underlying1.setValue(x)
        self.underlying2.setValue(y)
        return self.option.NPV()

    def setQuote(self, r):
        self.riskFreeQuote.setValue(r)
