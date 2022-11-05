
suppressMessages(library(QuantLib))

## global data
calendar <- TARGET()

settlementDate <- Date(15, "May", 1998)
settlementDate <- Calendar_adjust(calendar, settlementDate)

fixingDays <- 3
settlementDays <- 3
todaysDate <- Calendar_advance(calendar, settlementDate, -fixingDays, "Days")
invisible(Settings_instance()$setEvaluationDate(d=todaysDate))

cat('Today          : ', todaysDate$`__str__`(), "\n")
cat('Settlement Date: ', settlementDate$`__str__`(), "\n")


## option
exercise = AmericanExercise(todaysDate, Date(17, "May", 1999))
payoff = PlainVanillaPayoff("Put", 40.0)
option = VanillaOption(payoff, exercise)


# ### Market data

# %%
underlying = SimpleQuote(36.0)
dividendYield = FlatForward(todaysDate, 0.00, Actual365Fixed())
volatility = BlackConstantVol(todaysDate, calendar, 0.20, Actual365Fixed())
riskFreeRate = FlatForward(todaysDate, 0.06, Actual365Fixed())

# %%
process = BlackScholesMertonProcess(
    QuoteHandle(underlying),
    YieldTermStructureHandle(dividendYield),
    YieldTermStructureHandle(riskFreeRate),
    BlackVolTermStructureHandle(volatility)
)

# ### Pricing
# We'll collect tuples of method name, option value, and estimated error from the analytic formula.

# %%
results = NULL

# #### Analytic approximations

option$setPricingEngine(BaroneAdesiWhaleyApproximationEngine(process))
results = rbind(results, data.frame("method"="Barone-Adesi-Whaley", "NPV"=option$NPV()))

option$setPricingEngine(BjerksundStenslandApproximationEngine(process))
results = rbind(results, data.frame("method"="Bjerksund-Stensland", "NPV"=option$NPV()))


# #### Finite-difference method

timeSteps = 801
gridPoints = 800

option$setPricingEngine(FdBlackScholesVanillaEngine(process, timeSteps, gridPoints))
results = rbind(results, data.frame("method"="finite differences", "NPV"=option$NPV()))


# #### Binomial method

timeSteps = 801

# %%
listEngineFuncsNames = lsf.str("package:QuantLib")
listEngineFuncsNames = listEngineFuncsNames[grepl(pattern = "^Binomial", x = listEngineFuncsNames)]
listEngineFuncsNames = listEngineFuncsNames[grepl(pattern = "VanillaEngine$", x = listEngineFuncsNames)]

for (engineFuncName in listEngineFuncsNames) {
  eval(parse(text=paste0("option$setPricingEngine(", engineFuncName, "(process, timeSteps))")))
  results = rbind(results, data.frame("method"=engineFuncName, "NPV"=option$NPV()))
}


# ### Results
print(results)
