
# inspired by python example with the same name

suppressMessages(library(QuantLib))

## global data
calendar <- TARGET()

todaysDate <- DateParser_parseISO("2001-09-06")
invisible(Settings_instance()$setEvaluationDate(d=todaysDate))

settlementDays = 2
settlementDate = Calendar_advance(calendar, todaysDate, settlementDays, "Days")

cat('Today          : ', Date_ISO(todaysDate), "\n")
cat('Settlement Date: ', Date_ISO(settlementDate), "\n")


# ### Market quotes

# %%
deposits = list(
  "3M" = 0.0363
)

# %%
FRAs = list(
  "3x6" = 0.037125,
  "6x9" = 0.037125,
  "9x12" =0.037125
)


# %%
futures = list(
  "2001-12-19" = 96.2875,
  "2002-03-20" = 96.7875,
  "2002-06-19" = 96.9875,
  "2002-09-18" = 96.6875,
  "2002-12-18" = 96.4875,
  "2003-03-19" = 96.3875,
  "2003-06-18" = 96.2875,
  "2003-09-17" = 96.0875
)

# %%
swaps = list(
#  "2Y" = 0.037125,
  "3Y" = 0.0398,
  "5Y" = 0.0443,
  "10Y" = 0.05165,
  "15Y" = 0.055175
  
)

allHelpers_DepoFutSwap = RateHelperVector()
allHelpers_DepoFraSwap = RateHelperVector()

# Rate Helpers - Depos

for (it in names(deposits)) {
  itVal = deposits[[it]]  
  
  dayCounter = Actual360()
  
  dh = DepositRateHelper(
      QuoteHandle(SimpleQuote(itVal)),
      PeriodParser_parse(it),
      settlementDays,
      calendar,
      "ModifiedFollowing",
      F,
      dayCounter
    )  
  
  RateHelperVector_append(allHelpers_DepoFutSwap, dh)
  RateHelperVector_append(allHelpers_DepoFraSwap, dh)
}

# Rate Helpers - FRAs

for (it in names(FRAs)) {
  itVal = FRAs[[it]]  
  
  itSplit = strsplit(x = it, split = "x", fixed = T)[[1]]
  n = as.numeric(itSplit[1])
  m = as.numeric(itSplit[2])
  
  dayCounter = Actual360()
  months = 3

  dh = FraRateHelper(
    QuoteHandle(SimpleQuote(itVal)), n, m, settlementDays, calendar, "ModifiedFollowing", F, dayCounter
  )
  
  RateHelperVector_append(allHelpers_DepoFraSwap, dh)
}
  
  
# Rate Helpers - Futures

for (it in names(futures)) {
  itVal = futures[[it]]  

  dayCounter = Actual360()
  months = 3

  dh = FuturesRateHelper(
    QuoteHandle(SimpleQuote(itVal)),
    DateParser_parseISO(it),
    months,
    calendar,
    "ModifiedFollowing",
    T,
    dayCounter,
    QuoteHandle(SimpleQuote(0.0))    
  )
  
  RateHelperVector_append(allHelpers_DepoFutSwap, dh)
}


# 
# The discount curve for the swaps will come from elsewhere. 
# -> A real application would use some kind of risk-free curve; here we're using a flat one for convenience.
#

# %%
discountTermStructure = YieldTermStructureHandle(FlatForward(settlementDate, 0.04, Actual360()))

# %%

fixedLegFrequency = "Annual"
fixedLegTenor = PeriodParser_parse("1Y") 
fixedLegAdjustment = "Unadjusted"
fixedLegDayCounter = Thirty360(Thirty360_BondBasis_get())
floatingLegFrequency = "Quarterly"
floatingLegTenor = PeriodParser_parse("3M")
floatingLegAdjustment = "ModifiedFollowing"

for (it in names(swaps)) {
  itVal = swaps[[it]]  
  
  dh = SwapRateHelper(
    QuoteHandle(SimpleQuote(itVal)),
    PeriodParser_parse(it),
    calendar,
    fixedLegFrequency,
    fixedLegAdjustment,
    fixedLegDayCounter,
    Euribor3M(),
    QuoteHandle(),
    PeriodParser_parse("0D"),
    discountTermStructure
  )
  
  RateHelperVector_append(allHelpers_DepoFutSwap, dh)
  RateHelperVector_append(allHelpers_DepoFraSwap, dh)
}

# ### Term structure construction

printCurve <- function(curve, helpers) {
  df = NULL
  for (i in seq(1, helpers$size())) {
    h = helpers[i][[1]]
    pillar = RateHelper_pillarDate(h)

    day_counter = Actual365Fixed()
    compounding = "Continuous"
    
    r = YieldTermStructure_zeroRate(curve, pillar, day_counter, compounding, "Annual")
    disc = YieldTermStructure_discount(curve, pillar)
    
    df = rbind(df, data.frame(
      pillar=Date_ISO(pillar), zeroRate=InterestRate_rate(r)*100.0, df=disc, impliedRate=RateHelper_impliedQuote(h)
    ))
  }
  print(df)
}


# %%
forecastTermStructure = RelinkableYieldTermStructureHandle()

# %%
depoFuturesSwapCurve = PiecewiseFlatForward(settlementDate, allHelpers_DepoFutSwap, Actual360())
printCurve(depoFuturesSwapCurve, allHelpers_DepoFutSwap)


# %%
depoFraSwapCurve = PiecewiseFlatForward(settlementDate, allHelpers_DepoFraSwap, Actual360())
printCurve(depoFraSwapCurve, allHelpers_DepoFraSwap)



# ### Swap pricing

# %%
swapEngine = DiscountingSwapEngine(discountTermStructure)

# %%
nominal = 1000000
length = 5
maturity = Calendar_advance(calendar, settlementDate, length, "Years")
payFixed = T

# %%
fixedLegFrequency = "Annual"
fixedLegAdjustment = "Unadjusted"
fixedLegDayCounter = Thirty360(Thirty360_BondBasis_get())
fixedRate = swaps[["5Y"]] # 0.04

# %%
floatingLegFrequency = "Quarterly"
spread = 0.0
fixingDays = 2
index = Euribor3M(forecastTermStructure)
floatingLegAdjustment = "ModifiedFollowing"
floatingLegDayCounter = index$dayCounter()

# %%
fixedSchedule = Schedule(
    settlementDate,
    maturity,
    fixedLegTenor,
    calendar,
    fixedLegAdjustment,
    fixedLegAdjustment,
    DateGeneration_Forward_get(),
    F
)
floatingSchedule = Schedule(
    settlementDate,
    maturity,
    floatingLegTenor,
    calendar,
    floatingLegAdjustment,
    floatingLegAdjustment,
    DateGeneration_Forward_get(),
    F
)

# We'll build a 5-years swap starting spot...

# %%
spot = VanillaSwap(
    Swap_Payer_get(),
    nominal,
    fixedSchedule,
    fixedRate,
    fixedLegDayCounter,
    floatingSchedule,
    index,
    spread,
    floatingLegDayCounter
)
spot$setPricingEngine(swapEngine)



# ...and one starting 1 year forward.

# %%
forwardStart = Calendar_advance(calendar, settlementDate, 1, "Years")
forwardEnd = Calendar_advance(calendar, forwardStart, length, "Years")
fixedSchedule = Schedule(
    forwardStart,
    forwardEnd,
    fixedLegTenor,
    calendar,
    fixedLegAdjustment,
    fixedLegAdjustment,
    DateGeneration_Forward_get(),
    F
)
floatingSchedule = Schedule(
    forwardStart,
    forwardEnd,
    floatingLegTenor,
    calendar,
    floatingLegAdjustment,
    floatingLegAdjustment,
    DateGeneration_Forward_get(),
    F
)

# %%
forward = VanillaSwap(
    Swap_Payer_get(),
    nominal,
    fixedSchedule,
    fixedRate,
    fixedLegDayCounter,
    floatingSchedule,
    index,
    spread,
    floatingLegDayCounter
)
forward$setPricingEngine(swapEngine)


# We'll price them both on the bootstrapped curves.
#
# This is the quoted 5-years market rate; we expect the fair rate of the spot swap to match it.

# %%
print(swaps["5Y"])

showSwap <- function(swap) {
  print(paste("NPV         = ", swap$NPV()))
  print(paste("Fair spread = ", (swap$fairSpread()*100)))
  print(paste("Fair rate   = ", (swap$fairRate()*100)))
}


# %% [markdown]
# These are the results for the 5-years spot swap on the deposit/futures/swap curve...

# %%
forecastTermStructure$linkTo(depoFuturesSwapCurve)
showSwap(spot)

# ...and these are on the deposit/fra/swap curve.

# %%
forecastTermStructure$linkTo(depoFraSwapCurve)
showSwap(spot)

# %% [markdown]
# The same goes for the 1-year forward swap, except for the fair rate not matching the spot rate.

# %%
forecastTermStructure$linkTo(depoFuturesSwapCurve)
showSwap(forward)

# %%
forecastTermStructure$linkTo(depoFraSwapCurve)
showSwap(forward)

