
# inspired by python example with the same name

suppressMessages(library(QuantLib))

## global data
calendar <- TARGET()

todaysDate <- DateParser_parseISO("2014-04-30")
invisible(Settings_instance()$setEvaluationDate(d=todaysDate))

settlementDays <- 3
settlementDate <- Calendar_advance(calendar, todaysDate, settlementDays, "Days")

cat('Today          : ', Date_ISO(todaysDate), "\n")
cat('Settlement Date: ', Date_ISO(settlementDate), "\n")

refDate = todaysDate

# ### Calculations

# This exercise tries to replicate the Quantlib C++ `Gaussian1dModel` example on how to use the GSR and Markov Functional model.

# We assume a multicurve setup, for simplicity with flat yield term structures.
#
# The discounting curve is an Eonia curve at a level of 2% and the forwarding curve is an Euribor 6m curve at a level of 2.5%.
#
# For the volatility we assume a flat swaption volatility at 20%.


forward6mQuote = QuoteHandle(SimpleQuote(0.025))
oisQuote = QuoteHandle(SimpleQuote(0.02))
volQuote = QuoteHandle(SimpleQuote(0.2))

# %%
dc = Actual365Fixed()
yts6m = FlatForward(refDate, forward6mQuote, dc)
ytsOis = FlatForward(refDate, oisQuote, dc)
yts6m$enableExtrapolation()
ytsOis$enableExtrapolation()

hyts6m = RelinkableYieldTermStructureHandle(yts6m)
t0_curve = YieldTermStructureHandle(yts6m)
t0_Ois = YieldTermStructureHandle(ytsOis)
euribor6m = Euribor6M(hyts6m)
swaptionVol = ConstantSwaptionVolatility(0, calendar, "ModifiedFollowing", volQuote, Actual365Fixed())


# %%
effectiveDate = Calendar_advance(calendar, refDate, 2, "Days")
print(Date_ISO(effectiveDate))

maturityDate = Calendar_advance(calendar, effectiveDate, 10, "Years")
print(Date_ISO(maturityDate))

# %%
fixedSchedule = Schedule(effectiveDate,
                            maturityDate,
                            Period(1, "Years"),
                            calendar,
                            "ModifiedFollowing",
                            "ModifiedFollowing",
                            DateGeneration_Forward_get(), F)

# %%
floatSchedule = Schedule(effectiveDate,
                            maturityDate,
                            Period(6, "Months"),
                            calendar,
                            "ModifiedFollowing",
                            "ModifiedFollowing",
                            DateGeneration_Forward_get(), F)


# We consider a standard 10-years Bermudan payer swaption with yearly exercises at a strike of 4%.

# %%
fixedNominal    = rep(1, fixedSchedule$size()-1)
floatingNominal = rep(1, floatSchedule$size()-1)
strike          = rep(0.04, fixedSchedule$size()-1)
gearing         = rep(1, floatSchedule$size()-1)
spread          = rep(0, floatSchedule$size()-1)

# %%
underlying = NonstandardSwap(
    Swap_Payer_get(),
    fixedNominal, floatingNominal, fixedSchedule, strike,
    Thirty360(Thirty360_BondBasis_get()), floatSchedule,
    euribor6m, gearing, spread, Actual360(), F, F, "ModifiedFollowing")

# %%
exerciseDates = Schedule_dates(self = fixedSchedule, .copy = T)
for (i in seq(0, exerciseDates$size()-1)) { # C Style vector numbering
  DateVector___setitem__(exerciseDates, i, Calendar_advance(calendar, DateVector___getitem__(exerciseDates, i), -2, "Days"))
}

DateVector___delitem__(exerciseDates, 0)
DateVector___delitem__(exerciseDates, DateVector___len__(exerciseDates)-1)
for (i in seq(1, exerciseDates$size())) { # R style vector numbering
  print(exerciseDates[i][[1]])
}

exercise = BermudanExercise(exerciseDates)
swaption = NonstandardSwaption(underlying,exercise,Settlement_Physical_get())

# The model is a one factor Hull White model with piecewise volatility adapted to our exercise dates.
#
# The reversion is just kept constant at a level of 1%.

stepDates = DateVector(exerciseDates)
DateVector___delitem__(stepDates, DateVector___len__(stepDates)-1)

sigmas = QuoteHandleVector()
for (i in seq(1, 9)) {
  QuoteHandleVector_append(sigmas, QuoteHandle(SimpleQuote(0.01)))
}

reversion = QuoteHandleVector()
QuoteHandleVector_append(reversion, QuoteHandle(SimpleQuote(0.01)))

#
# The model's curve is set to the 6m forward curve. 
# Note that the model adapts automatically to other curves where appropriate 
# (e.g. if an index requires a different forwarding curve) or where explicitly specified (e.g. in a swaption pricing engine).


gsr = Gsr(t0_curve, stepDates, sigmas, reversion)
swaptionEngine = Gaussian1dSwaptionEngine(gsr, 64, 7.0, T, F, t0_Ois)
nonstandardSwaptionEngine = Gaussian1dNonstandardSwaptionEngine(
    gsr, 64, 7.0, T, F, QuoteHandle(SimpleQuote(0)), t0_Ois)

# %%
swaption$setPricingEngine(nonstandardSwaptionEngine)

# %%
swapBase = EuriborSwapIsdaFixA(Period(10, 'Years'), t0_curve, t0_Ois)
basket = swaption$calibrationBasket(swapBase, swaptionVol, 'Naive')

# %%
for (i in seq(1, basket$size())) {
  basket_i = basket[i][[1]]
  BlackCalibrationHelper_setPricingEngine(basket_i, swaptionEngine)
}
    

# %%
method = LevenbergMarquardt()
ec = EndCriteria(1000, 10, 1e-8, 1e-8, 1e-8)

# %%
gsr$calibrateVolatilitiesIterative(basket, method, ec)


# %% [markdown]
# The engine can generate a calibration basket in two modes.
#
# The first one is called Naive and generates ATM swaptions adapted to the exercise dates of the swaption and its maturity date. The resulting basket looks as follows:

# %%
basket_data = function(basket) {
  df = NULL
  for (i in seq(1, basket$size())) {
    basket_i = basket[i][[1]]
    h = as_swaption_helper(basket_i)
    hopt = SwaptionHelper_swaption(h)
    df = rbind(df, data.frame(expiry=Date_ISO(SwaptionHelper_swaptionExpiryDate(h)),
                              maturity=Date_ISO(SwaptionHelper_swaptionMaturityDate(h)),
                              nominal=SwaptionHelper_swaptionNominal(h),
                              strike=SwaptionHelper_swaptionStrike(h),
                              optType=Swaption_type(hopt)
                              
    ))
  }
  df
}  

print(basket_data(basket))




#
# Let's calibrate our model to this basket. 
# We use a specialized calibration method calibrating the sigma function one by one to the calibrating vanilla swaptions. The result of this is as follows:

# %%
calibration_data = function(basket, volatilities) {
  # volatilities = gsr$volatility()
  df = NULL
  for (i in seq(1, basket$size())) {
    basket_i = basket[i][[1]]
    vola_i = volatilities[i][[1]]
    h = as_swaption_helper(basket_i)
    hopt = SwaptionHelper_swaption(h)
    modelPrice = basket_i$modelValue()
    modelImpVol = BlackCalibrationHelper_impliedVolatility(self = basket_i, targetValue = modelPrice, accuracy = 1e-6, maxEvaluations = 1000, minVol = 0.0, maxVol = 2.0)
    marketPrice = basket_i$marketValue()
    
    df = rbind(df, data.frame(expiry=Date_ISO(SwaptionHelper_swaptionExpiryDate(h)),
                              modelSigma=vola_i,
                              modelPrice=modelPrice,
                              marketPrice=marketPrice,
                              modelImpVol=modelImpVol,
                              marketImpVol=QuoteHandle_value(basket_i$volatility())
    ))
  }  
                              
  df
}

print(calibration_data(basket, gsr$volatility()))

# %% [markdown]
# Bermudan swaption NPV (ATM calibrated GSR):

# %%
print(swaption$NPV())


#
# There is another mode to generate a calibration basket called `MaturityStrikeByDeltaGamma`.
# This means that the maturity, the strike and the nominal of the calibrating swaptions are obtained matching the NPV, 
# first derivative and second derivative of the swap you will exercise into at at each bermudan call date.
# The derivatives are taken with respect to the model's state variable.
#
# Let's try this in our case.

# %%
basket = swaption$calibrationBasket(swapBase, swaptionVol, 'MaturityStrikeByDeltaGamma')
print(basket_data(basket))

# %%
for (i in seq(1, basket$size())) {
  basket_i = basket[i][[1]]
  BlackCalibrationHelper_setPricingEngine(basket_i, swaptionEngine)
}


# %% [markdown]
# The calibrated nominal is close to the exotics nominal. The expiries and maturity dates of the vanillas are the same as in the case above. The difference is the strike which is now equal to the exotics strike.
#
# Let's see how this affects the exotics NPV. The recalibrated model is:

# %%
gsr$calibrateVolatilitiesIterative(basket, method, ec)
print(calibration_data(basket, gsr$volatility()))

# %% [markdown]
# Bermudan swaption NPV (deal strike calibrated GSR):

# %%
print(swaption$NPV())

# 
# We can do more complicated things.  
# Let's e.g. modify the nominal schedule to be linear amortizing and see what the effect on the generated calibration basket is:

# %%
for (i in seq(1,fixedSchedule$size()-1)) {
  tmp = 1.0 - i/ (fixedSchedule$size())
  fixedNominal[i]        = tmp
  floatingNominal[i*2-1]   = tmp
  floatingNominal[i*2] = tmp
}

# %%
underlying2 = NonstandardSwap(Swap_Payer_get(),
                            fixedNominal, floatingNominal, fixedSchedule, strike,
                            Thirty360(Thirty360_BondBasis_get()), floatSchedule,
                            euribor6m, gearing, spread, Actual360(), F, F, "ModifiedFollowing")

# %%
swaption2 = NonstandardSwaption(underlying2, exercise, Settlement_Physical_get())

# %%
swaption2$setPricingEngine(nonstandardSwaptionEngine)
basket = swaption2$calibrationBasket(swapBase, swaptionVol, 'MaturityStrikeByDeltaGamma')

# %%
print(basket_data(basket))

#
# The notional is weighted over the underlying exercised into and the maturity is adjusted downwards. The rate, on the other hand, is not affected.

#
# You can also price exotic bond's features. 
# If you have e.g. a Bermudan callable fixed bond you can set up the call right as a swaption 
# to enter into a one leg swap with notional reimbursement at maturity. 
# The exercise should then be written as a rebated exercise paying the notional in case of exercise. The calibration basket looks like this:

# %%
fixedNominal2    = rep(1, fixedSchedule$size()-1)
floatingNominal2 = rep(0, fixedSchedule$size()*2-2) #null the second leg

# %%
underlying3 = NonstandardSwap(Swap_Receiver_get(),
                            fixedNominal2, floatingNominal2, fixedSchedule, strike,
                            Thirty360(Thirty360_BondBasis_get()), floatSchedule,
                            euribor6m, gearing, spread, Actual360(), F, T, "ModifiedFollowing")

# %%
rebateAmount = rep(-1, exerciseDates$size())
exercise2 = RebatedExercise(exercise, rebateAmount, 2, calendar)
swaption3 = NonstandardSwaption(underlying3, exercise2, Settlement_Physical_get())

# %%
oas0 = SimpleQuote(0)
oas100 = SimpleQuote(0.01)
oas = RelinkableQuoteHandle(oas0)

# %%
nonstandardSwaptionEngine2 = Gaussian1dNonstandardSwaptionEngine(
    gsr, 64, 7.0, T, F, oas, t0_curve) # Change discounting to 6m

# %%
swaption3$setPricingEngine(nonstandardSwaptionEngine2)
basket = swaption3$calibrationBasket(swapBase, swaptionVol, 'MaturityStrikeByDeltaGamma')

# %%
print(basket_data(basket))

# %% [markdown]
# Note that nominals are not exactly 1.0 here. This is because we do our bond discounting on 6m level while the swaptions are still discounted on OIS level. (You can try this by changing the OIS level to the 6m level, which will produce nominals near 1.0).
#
# The NPV of the call right is (after recalibrating the model):

# %%
for (i in seq(1, basket$size())) {
  basket_i = basket[i][[1]]
  BlackCalibrationHelper_setPricingEngine(basket_i, swaptionEngine)
}

# %%
gsr$calibrateVolatilitiesIterative(basket, method, ec)

# %%
print(swaption3$NPV())

# %% [markdown]
# Up to now, no credit spread is included in the pricing. We can do so by specifying an oas in the pricing engine. Let's set the spread level to 100bp and regenerate the calibration basket.

# %%
oas$linkTo(oas100)
basket = swaption3$calibrationBasket(swapBase, swaptionVol, 'MaturityStrikeByDeltaGamma')
print(basket_data(basket))

# %% [markdown]
# The adjusted basket takes the credit spread into account. This is consistent to a hedge where you would have a margin on the float leg around 100bp,too.

# %%
for (i in seq(1, basket$size())) {
  basket_i = basket[i][[1]]
  BlackCalibrationHelper_setPricingEngine(basket_i, swaptionEngine)
}

# %%
gsr$calibrateVolatilitiesIterative(basket, method, ec)

# %%
print(swaption3$NPV())

