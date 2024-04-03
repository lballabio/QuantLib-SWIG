
# inspired by python example with the same name

suppressMessages(library(QuantLib))

## global data
calendar <- TARGET()

todaysDate <- Date(19, "October", 2020)
invisible(Settings_instance()$setEvaluationDate(d=todaysDate))

settlementDays <- 3
settlementDate <- Calendar_advance(calendar, todaysDate, settlementDays, "Days")

cat('Today          : ', todaysDate$`__str__`(), "\n")
cat('Settlement Date: ', settlementDate$`__str__`(), "\n")

# ### Term structure construction

# %%
dates = DateVector()
DateVector_append(dates, DateParser_parseISO("2020-10-19"))
DateVector_append(dates, DateParser_parseISO("2020-11-19"))

DateVector_append(dates, DateParser_parseISO("2021-01-19"))
DateVector_append(dates, DateParser_parseISO("2021-04-19"))
DateVector_append(dates, DateParser_parseISO("2021-10-19"))

DateVector_append(dates, DateParser_parseISO("2022-04-19"))
DateVector_append(dates, DateParser_parseISO("2022-10-19"))

DateVector_append(dates, DateParser_parseISO("2023-10-19"))
DateVector_append(dates, DateParser_parseISO("2025-10-19"))
DateVector_append(dates, DateParser_parseISO("2030-10-19"))
DateVector_append(dates, DateParser_parseISO("2035-10-19"))
DateVector_append(dates, DateParser_parseISO("2040-10-19"))

DateVector_size(dates)

rates = c(
    -0.004,
    -0.002,
    0.001,
    0.005,
    0.009,
    0.010,
    0.010,
    0.012,
    0.017,
    0.019,
    0.028,
    0.032
)

length(rates)

forecast_curve = ZeroCurve(dates, rates, Actual365Fixed())
forecast_handle = YieldTermStructureHandle(forecast_curve)

# ### Swap construction
#
# We'll use an overnight swap as an example.  
# We're keeping the initialization simple, 
# but the analysis work in the same way for more complex ones, 
# as well as for other kinds of swaps and bonds (once we extract the cashflows from them using the proper methods).

# %%
swapBuilder = MakeOIS(swapTenor=Period(5, "Years"),
                  overnightIndex=Eonia(forecast_handle),
                  fixedRate=0.002)

swap = MakeOIS_makeOIS(swapBuilder)

# ### Cash-flow analysis
#
# The fixed-rate coupons can be extracted from the swap using the `fixedLeg` method.  
#

fixed_leg = swap$fixedLeg()
CashFlows_maturityDate(fixed_leg)

df = NULL
for (i in seq(1, fixed_leg$size())) {
  cfl = fixed_leg[i][[1]]
  df = rbind(df, data.frame(date=Date_ISO(CashFlow_date(cfl)), amount=CashFlow_amount(cfl)))
}
print(df)


#
# If we want to extract more information, we need to upcast the coupons to a more specific class.  
# This can be done by using the `as_coupon` or the `as_fixed_rate_coupon` method. 
#

df = NULL
for (i in seq(1, fixed_leg$size())) {
  cfl = fixed_leg[i][[1]]
  cflCoupon = as_coupon(cfl)
  cflCouponFix = as_fixed_rate_coupon(cfl)
  cflCouponFixIr = FixedRateCoupon_interestRate(cflCouponFix)
  
  df = rbind(df, data.frame(date=Date_ISO(CashFlow_date(cfl)), amount=CashFlow_amount(cfl), 
                            accrualStartDate=Date_ISO(Coupon_accrualStartDate(cflCoupon)),
                            accrualEndDate=Date_ISO(Coupon_accrualEndDate(cflCoupon)),
                            accrualPeriod=Coupon_accrualPeriod(cflCoupon),
                            accrualDays=Coupon_accrualDays(cflCoupon),
                            accrualDayCounter=DayCounter_name(Coupon_dayCounter(cflCoupon)),
                            accruedAmount=Coupon_accruedAmount(self = cflCoupon, todaysDate),
                            rate=InterestRate_rate(cflCouponFixIr), rateDayCount=DayCounter_name(InterestRate_dayCounter(cflCouponFixIr))
                            ))
}
print(df)



# %%
floating_leg = swap$overnightLeg()

# %%
df = NULL
for (i in seq(1, floating_leg$size())) {
  cfl = floating_leg[i][[1]]
  cflCoupon = as_coupon(cfl)
  cflCouponFloat = as_floating_rate_coupon(cfl)
  cflCouponFloatIdx = FloatingRateCoupon_index(cflCouponFloat)
  
  
  df = rbind(df, data.frame(date=Date_ISO(CashFlow_date(cfl)), amount=CashFlow_amount(cfl), 
                            accrualStartDate=Date_ISO(Coupon_accrualStartDate(cflCoupon)),
                            accrualEndDate=Date_ISO(Coupon_accrualEndDate(cflCoupon)),
                            accrualPeriod=Coupon_accrualPeriod(cflCoupon),
                            accrualDays=Coupon_accrualDays(cflCoupon),
                            accrualDayCounter=DayCounter_name(Coupon_dayCounter(cflCoupon)),
                            accruedAmount=Coupon_accruedAmount(self = cflCoupon, todaysDate),
                            spread=FloatingRateCoupon_spread(cflCouponFloat),
                            gearing=FloatingRateCoupon_gearing(cflCouponFloat),
                            adjFix=FloatingRateCoupon_adjustedFixing(cflCouponFloat),
                            dtFix=Date_ISO(FloatingRateCoupon_fixingDate(cflCouponFloat)),
                            convexityAdj=FloatingRateCoupon_convexityAdjustment(cflCouponFloat)
  ))
}

print(df)

