
# inspired by python example with the same name

suppressMessages(library(QuantLib))

## global data
calendar <- TARGET()

todaysDate <- DateParser_parseISO("2019-09-26")
invisible(Settings_instance()$setEvaluationDate(d=todaysDate))

settlementDays <- 2
settlementDate <- Calendar_advance(calendar, todaysDate, settlementDays, "Days")

cat('Today          : ', Date_ISO(todaysDate), "\n")
cat('Settlement Date: ', Date_ISO(settlementDate), "\n")

spot = settlementDate

# %% [markdown]
# ### Data
#
# We'll use the following data as input:

# %%
refMktRates = c(
    -0.373,
    -0.388,
    -0.402,
    -0.418,
    -0.431,
    -0.441,
    -0.45,
    -0.457,
    -0.463,
    -0.469,
    -0.461,
    -0.463,
    -0.479,
    -0.4511,
    -0.45418,
    -0.439,
    -0.4124,
    -0.37703,
    -0.3335,
    -0.28168,
    -0.22725,
    -0.1745,
    -0.12425,
    -0.07746,
    0.0385,
    0.1435,
    0.17525,
    0.17275,
    0.1515,
    0.1225,
    0.095,
    0.0644
)

# ### Market instruments

# %%
index = Euribor6M()

# The first market rate is for the 6-months deposit...

# %%
helpers = RateHelperVector()

helpers$append(
    DepositRateHelper(
        refMktRates[1] / 100.0, Period(6, "Months"), 2, calendar, "ModifiedFollowing", T, Actual360()
    )
)

# ...the next 12 are for FRAs...

# %%
for (i in seq(2, 13)) {
  helpers$append(
    FraRateHelper(refMktRates[i] / 100.0, i + 1, index)
  )
}

# ...and the others are swap rates.

# %%
swapTenors = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20, 25, 30, 35, 40, 45, 50)

for (i in seq(14, length(refMktRates))) {
  r = refMktRates[i]
  tenor = swapTenors[i-13]
  helpers$append(
    SwapRateHelper(
      r / 100.0, Period(tenor, "Years"), calendar, "Annual", "ModifiedFollowing", Thirty360(Thirty360_BondBasis_get()), index
    )
  )
}

# We'll also add a few synthetic helpers:

# %%
additional_helpers = RateHelperVector()
for (i in seq(0, 6)){
  additional_helpers$append(
    FraRateHelper(-0.004, 12 + i, index) 
  )
}

# %%
additional_dates = DateVector()
for (i in seq(0, 4)){
  additional_dates$append(
    Calendar_advance(calendar, spot, Period(i+1, "Months"))
  )
}



# ### Global bootstrap
#
# This curve takes into account the market instruments, as well as the passed additional ones.

# %%
curve = GlobalLinearSimpleZeroCurve(
    spot, helpers, Actual365Fixed(), GlobalBootstrap(additional_helpers, additional_dates, 1.0e-4)
)
curve$enableExtrapolation()


# ### Report
df = NULL
for (i in seq(1, helpers$size())) {
  h = helpers[i][[1]]
  pillar = RateHelper_pillarDate(h)
  
  day_counter = Actual360()
  compounding = "Simple"
  if (i > 13) {
    day_counter = Thirty360(Thirty360_BondBasis_get())
    compounding = "SimpleThenCompounded"
  }

  r = YieldTermStructure_zeroRate(curve, pillar, day_counter, compounding, "Annual")
  
  df = rbind(df, data.frame(
    pillar=Date_ISO(pillar), zeroRate=InterestRate_rate(r)*100.0
  ))
}
print(df)

