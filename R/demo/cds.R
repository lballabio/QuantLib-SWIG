
# inspired by python example with the same name

suppressMessages(library(QuantLib))

## global data
calendar <- TARGET()

todaysDate <- DateParser_parseISO("2007-05-15")
invisible(Settings_instance()$setEvaluationDate(d=todaysDate))

settlementDays <- 3
settlementDate <- Calendar_advance(calendar, todaysDate, settlementDays, "Days")

cat('Today          : ', Date_ISO(todaysDate), "\n")
cat('Settlement Date: ', Date_ISO(settlementDate), "\n")

risk_free_rate = YieldTermStructureHandle(FlatForward(todaysDate, 0.01, Actual365Fixed()))

# ### CDS parameters

recovery_rate = 0.5
quoted_spreads = c(0.0150, 0.0150, 0.0150, 0.0150)

tenors = PeriodVector()
tenors$append(Period(3, "Months"))
tenors$append(Period(6, "Months"))
tenors$append(Period(1, "Years"))
tenors$append(Period(2, "Years"))

maturities = DateVector()
for (i in seq(1, tenors$size())) {
  maturities$append(
    Calendar_adjust(calendar,
                    Calendar_advance(calendar, todaysDate, tenors[i][[1]]),
                    "Following")
  )  
}

instruments = DefaultProbabilityHelperVector()
for (i in seq(1, length(quoted_spreads))) {
  s = quoted_spreads[i]
  tenor = tenors[i][[1]]
  
  it = SpreadCdsHelper(
    QuoteHandle(SimpleQuote(s)),
    tenor,
    0,
    calendar,
    "Quarterly",
    "Following",
    DateGeneration_TwentiethIMM_get(),
    Actual365Fixed(),
    recovery_rate,
    risk_free_rate
  )
  
  instruments$append(it) 
}

hazard_curve = PiecewiseFlatHazardRate(todaysDate, instruments, Actual365Fixed())

print("Calibrated hazard rate values: ")
for (i in seq(1, hazard_curve$dates()$size())) {
  print(paste("hazard rate on ", 
              Date_ISO(hazard_curve$dates()[i][[1]]), 
              "is ", 
              HazardRateCurve_hazardRates(hazard_curve)[i]))
}

print(paste("1Y survival probability:",
            DefaultProbabilityTermStructure_survivalProbability(hazard_curve, Calendar_advance(calendar, todaysDate, Period(1, "Years"))),
            "expected 0.9704 "
))

print(paste("2Y survival probability:",
            DefaultProbabilityTermStructure_survivalProbability(hazard_curve, Calendar_advance(calendar, todaysDate, Period(2, "Years"))),
            "expected 0.9418 "
))

# ### Reprice instruments

nominal = 1000000.0
probability = DefaultProbabilityTermStructureHandle(hazard_curve)

# We'll create a cds for every maturity:

print("Repricing of quoted CDSs employed for calibration: ")
for (i in seq(1, length(quoted_spreads))) {
  maturity = maturities[i][[1]]
  s = quoted_spreads[i]
  tenor = tenors[i][[1]]
  
  schedule = Schedule(
    todaysDate,
    maturity,
    Period(3, "Months"), # "Quarterly",
    calendar,
    "Following",
    "Unadjusted",
    DateGeneration_TwentiethIMM_get(),
    F
  )
  
  cds = CreditDefaultSwap(Protection_Seller_get(), nominal, s, schedule, "Following", Actual365Fixed())
  engine = MidPointCdsEngine(probability, recovery_rate, risk_free_rate)
  cds$setPricingEngine(engine)
  
  print(paste("fair spread: ", Period___str__(tenor), cds$fairSpread()))
  print(paste("   NPV: ", cds$NPV()))
  print(paste("   default leg: ", cds$defaultLegNPV()))
  print(paste("   coupon leg: ", cds$couponLegNPV()))
  print("")
  
}

