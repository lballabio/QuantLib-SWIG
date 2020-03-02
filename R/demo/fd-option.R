
suppressMessages(library(QuantLib))

quote <- function(x) { sapply(x, function(x1) {
  todaysDate <- Date(15, "May", 1998)
  invisible(Settings_instance()$setEvaluationDate(d=todaysDate))
  settlementDate <- Date(17, "May", 1998)
  riskFreeRate <- FlatForward(settlementDate, 0.05, Actual365Fixed())
  exercise <- EuropeanExercise(Date(17, "May", 1999))
  payoff <- PlainVanillaPayoff("Call", 8.0)
  underlying <- SimpleQuote(x1)
  volatility <- BlackConstantVol(todaysDate, TARGET(), 0.10, Actual365Fixed())
  dividendYield <- FlatForward(settlementDate, 0.05, Actual365Fixed())
  process <- BlackScholesMertonProcess(QuoteHandle(underlying),
                                       YieldTermStructureHandle(dividendYield),
                                       YieldTermStructureHandle(riskFreeRate),
                                       BlackVolTermStructureHandle(volatility))
  option <- VanillaOption(payoff, exercise)
  invisible(option$setPricingEngine(s_arg2=FdBlackScholesVanillaEngine(process, 25)))
  option$NPV()
})}

curve(quote,xlim=c(1,20))
