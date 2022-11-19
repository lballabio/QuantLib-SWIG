
# inspired by python example with the same name

suppressMessages(library(QuantLib))

## global data
calendar <- TARGET()

todaysDate <- Date(15, "May", 1998)
invisible(Settings_instance()$setEvaluationDate(d=todaysDate))

settlementDays <- 3
settlementDate <- Calendar_advance(calendar, todaysDate, settlementDays, "Days")

cat('Today          : ', todaysDate$`__str__`(), "\n")
cat('Settlement Date: ', settlementDate$`__str__`(), "\n")

riskFreeRate = FlatForward(settlementDate, 0.05, Actual365Fixed())


# ### Option parameters
exerciseDate = Date(17, "May", 1999)
exercise = EuropeanExercise(exerciseDate)
payoff = PlainVanillaPayoff("Call", 8.0)


# ### Market data
underlying1 = SimpleQuote(7.0)
volatility1 = BlackConstantVol(todaysDate, calendar, 0.10, Actual365Fixed())
dividendYield1 = FlatForward(settlementDate, 0.05, Actual365Fixed())
underlying2 = SimpleQuote(7.0)
volatility2 = BlackConstantVol(todaysDate, calendar, 0.10, Actual365Fixed())
dividendYield2 = FlatForward(settlementDate, 0.05, Actual365Fixed())


process1 = BlackScholesMertonProcess(
    QuoteHandle(underlying1),
    YieldTermStructureHandle(dividendYield1),
    YieldTermStructureHandle(riskFreeRate),
    BlackVolTermStructureHandle(volatility1)
)

process2 = BlackScholesMertonProcess(
    QuoteHandle(underlying2),
    YieldTermStructureHandle(dividendYield2),
    YieldTermStructureHandle(riskFreeRate),
    BlackVolTermStructureHandle(volatility2)
)

corrMatrix = Matrix(2, 2)
invisible(Matrix_setitem(corrMatrix, 0, 0, 1.0))
invisible(Matrix_setitem(corrMatrix, 1, 1, 1.0))
invisible(Matrix_setitem(corrMatrix, 1, 0, 0.5))
invisible(Matrix_setitem(corrMatrix, 0, 1, 0.5))
cat(corrMatrix$`__str__`())

procVector = StochasticProcess1DVector(2)
invisible(StochasticProcess1DVector___setitem__(self = procVector, i = 0, x = process1))
invisible(StochasticProcess1DVector___setitem__(self = procVector, i = 1, x = process2))

process = StochasticProcessArray(array = procVector, correlation = corrMatrix)

# ### Pricing - European

basketoption = BasketOption(MaxBasketPayoff(payoff), exercise)
basketoption$setPricingEngine(
  MCPREuropeanBasketEngine(process=process, timeSteps=NA, timeStepsPerYear=1, brownianBridge=F, antitheticVariate=F, requiredSamples=NA, requiredTolerance=0.02, maxSamples=10000, seed=42)
)
print(basketoption$NPV())

basketoption = BasketOption(MinBasketPayoff(payoff), exercise)
basketoption$setPricingEngine(
  MCPREuropeanBasketEngine(process=process, timeSteps=NA, timeStepsPerYear=1, brownianBridge=F, antitheticVariate=F, requiredSamples=NA, requiredTolerance=0.02, maxSamples=10000, seed=42)
)
print(basketoption$NPV())

basketoption = BasketOption(AverageBasketPayoff(payoff, 2), exercise)
basketoption$setPricingEngine(
  MCPREuropeanBasketEngine(process=process, timeSteps=NA, timeStepsPerYear=1, brownianBridge=F, antitheticVariate=F, requiredSamples=NA, requiredTolerance=0.02, maxSamples=10000, seed=42)
)
print(basketoption$NPV())


# ### Pricing - American

americanExercise = AmericanExercise(settlementDate, exerciseDate)
americanbasketoption = BasketOption(MaxBasketPayoff(payoff), americanExercise)
americanbasketoption$setPricingEngine(
  MCPRAmericanBasketEngine(
        process=process,
        timeSteps=10,
        timeStepsPerYear=NA,
        brownianBridge=F,
        antitheticVariate=F,
        requiredSamples=50000,
        requiredTolerance=0.02,
        maxSamples=100000,
        seed=42,
        nCalibrationSamples=5000,
        polynomOrder=5,
        polynomType=LsmBasisSystem_Hermite_get()
    )
)
print(americanbasketoption$NPV())
