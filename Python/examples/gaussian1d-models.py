
# Copyright (C) 2018 Angus Lee
#
# This file is part of QuantLib, a free-software/open-source library
# for financial quantitative analysts and developers - http://quantlib.org/
#
# QuantLib is free software: you can redistribute it and/or modify it under the
# terms of the QuantLib license.  You should have received a copy of the
# license along with this program; if not, please email
# <quantlib-dev@lists.sf.net>. The license is also available online at
# <http://quantlib.org/license.shtml>.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the license for more details.

import QuantLib as ql

def printBasket(basket):
    print ("%-20s %-20s %-20s %-20s %-20s %-20s" % ("Expiry", "Maturity", "Nominal", "Rate", "MarketVol", "Pay/Rec"))
    print ("==================================================================================================================")

    for i in range(0, len(basket)):
        expiryDate = ql.as_black_helper(basket[i]).swaptionExpiryDate()
        endDate = ql.as_black_helper(basket[i]).swaptionMaturityDate()
        nominal = ql.as_black_helper(basket[i]).swaptionNominal()
        vol     = ql.as_black_helper(basket[i]).volatility().value()
        rate    = ql.as_black_helper(basket[i]).swaptionStrike()
        print ("%-20s %-20s %-20f %-20f %-20f" % (str(expiryDate), str(endDate), nominal, rate, vol))

    print("==================================================================================================================")

def printModelCalibration(basket, volatility):
    print ("%-20s %-20s %-20s %-20s %-20s %-20s" % ("Expiry","Model sigma","ModelPrice","MarketPrice","Model impVol","Market impVol"))
    print ("=================================================================================================================")

    for i in range(0, len(basket)):
        expiryDate = ql.as_black_helper(basket[i]).swaptionExpiryDate()
        modelValue = ql.as_black_helper(basket[i]).modelValue()
        marketValue= ql.as_black_helper(basket[i]).marketValue()
        impVol     = ql.as_black_helper(basket[i]).impliedVolatility(modelValue, 1e-6, 1000, 0.0, 2.0)
        vol    = ql.as_black_helper(basket[i]).volatility().value()
        print ("%-20s %-20f %-20f %-20f %-20f %-20f" % (str(expiryDate), volatility[i], modelValue, marketValue, impVol, vol))

    print("==================================================================================================================")

print('This exercise tries to replicate the Quantlib C++ Gaussian1dModel example on how to use the GSR and Markov Functional model.')

refDate = ql.Date(30, 4, 2014)
ql.Settings.instance().setEvaluationDate(refDate)
forward6mQuote = ql.QuoteHandle(ql.SimpleQuote(0.025))
oisQuote = ql.QuoteHandle(ql.SimpleQuote(0.02))
volQuote = ql.QuoteHandle(ql.SimpleQuote(0.2))

print("\nThe evaluation date for this example is set to", refDate)

print('\nWe assume a multicurve setup, for simplicity with flat yield term structures. '
      '\nThe discounting curve is an Eonia curve at a level of %f '
      '\nand the forwarding curve is an Euribior 6m curve at a level of %f. ' % (oisQuote.value(), forward6mQuote.value()))

print("\nFor the volatility we assume a flat swaption volatility at %f. " % volQuote.value())

dc = ql.Actual365Fixed()
yts6m = ql.FlatForward(refDate, forward6mQuote, dc)
ytsOis= ql.FlatForward(refDate, oisQuote, dc)
yts6m.enableExtrapolation()
ytsOis.enableExtrapolation()
hyts6m = ql.RelinkableYieldTermStructureHandle(yts6m)
t0_curve = ql.YieldTermStructureHandle(yts6m)
t0_Ois = ql.YieldTermStructureHandle(ytsOis)
euribor6m = ql.Euribor6M(hyts6m)
swaptionVol = ql.ConstantSwaptionVolatility(0, ql.TARGET(), ql.ModifiedFollowing, volQuote, ql.Actual365Fixed())

effectiveDate = ql.TARGET().advance(refDate, ql.Period('2D'))
maturityDate = ql.TARGET().advance(effectiveDate, ql.Period('10Y'))

fixedSchedule = ql.Schedule(effectiveDate,
                            maturityDate,
                            ql.Period('1Y'),
                            ql.TARGET(),
                            ql.ModifiedFollowing,
                            ql.ModifiedFollowing,
                            ql.DateGeneration.Forward, False)

floatSchedule = ql.Schedule(effectiveDate,
                            maturityDate,
                            ql.Period('6M'),
                            ql.TARGET(),
                            ql.ModifiedFollowing,
                            ql.ModifiedFollowing,
                            ql.DateGeneration.Forward, False)

# Vector input for the NonstandardSwap obj
fixedNominal    = [1]*(len(fixedSchedule)-1)
floatingNominal = [1]*(len(floatSchedule)-1)
strike          = [0.04]*(len(fixedSchedule)-1)
gearing         = [1]*(len(floatSchedule)-1)
spread          = [0]*(len(floatSchedule)-1)

print("\nWe consider a standard 10y bermudan payer swaption "
      "\nwith yearly exercises at a strike of %f" % strike[0])

underlying = ql.NonstandardSwap(ql.VanillaSwap.Payer,
                            fixedNominal, floatingNominal, fixedSchedule, strike,
                            ql.Thirty360(), floatSchedule,
                            euribor6m, gearing, spread, ql.Actual360(), False, False, ql.ModifiedFollowing)

exerciseDates = [ql.TARGET().advance(x, -ql.Period('2D')) for x in fixedSchedule]
exerciseDates = exerciseDates[1:-1]
exercise = ql.BermudanExercise(exerciseDates)
swaption = ql.NonstandardSwaption(underlying,exercise,ql.Settlement.Physical)

stepDates = exerciseDates[:-1]
sigmas = [ql.QuoteHandle(ql.SimpleQuote(0.01)) for x in range(1, 10)]
reversion = [ql.QuoteHandle(ql.SimpleQuote(0.01))]

print("\nThe model is a one factor Hull White model with piecewise "
      "\nvolatility adapted to our exercise dates."
      "\nThe reversion is just kept constant at a level of %f" % reversion[0].value())

print("\nThe model's curve is set to the 6m forward curve. Note that "
      "\nthe model adapts automatically to other curves where appropriate "
      "\n(e.g. if an index requires a different forwarding curve) or "
      "\nwhere explicitly specified (e.g. in a swaption pricing engine)." )

gsr = ql.Gsr(t0_curve, stepDates, sigmas, reversion)
swaptionEngine = ql.Gaussian1dSwaptionEngine(gsr, 64, 7.0, True, False, t0_Ois)
nonstandardSwaptionEngine = ql.Gaussian1dNonstandardSwaptionEngine(gsr, 64, 7.0, True, False, ql.QuoteHandle(ql.SimpleQuote(0)), t0_Ois)

swaption.setPricingEngine(nonstandardSwaptionEngine)

swapBase = ql.EuriborSwapIsdaFixA(ql.Period('10Y'), t0_curve, t0_Ois)
basket = swaption.calibrationBasket(swapBase, swaptionVol, 'Naive')

for basket_i in basket:
    ql.as_black_helper(basket_i).setPricingEngine(swaptionEngine)

method = ql.LevenbergMarquardt()
ec = ql.EndCriteria(1000, 10, 1e-8, 1e-8, 1e-8)

gsr.calibrateVolatilitiesIterative(basket, method, ec)


print("\nThe engine can generate a calibration basket in two modes."
      "\nThe first one is called Naive and generates ATM swaptions adapted to"
      "\nthe exercise dates of the swaption and its maturity date"
      "\nThe resulting basket looks as follows:" )
printBasket(basket)
print("\nLet's calibrate our model to this basket. We use a specialized"
       "\ncalibration method calibrating the sigma function one by one to"
       "\nthe calibrating vanilla swaptions. The result of this is as follows:" )
printModelCalibration(basket, gsr.volatility())

npv = swaption.NPV()
print("Bermudan swaption NPV (ATM calibrated GSR) = %f" % npv)

# Recalibration to 4$ Strike swaption
print( "\nThere is another mode to generate a calibration basket called"
       "\nMaturityStrikeByDeltaGamma. This means that the maturity, "
       "\nthe strike and the nominal of the calibrating swaption are "
       "\ncomputed such that the npv and its first and second "
       "\nderivative with respect to the model's state variable) of"
       "\nthe exotics underlying match with the calibrating swaption's"
       "\nunderlying. Let's try this in our case." )

basket = swaption.calibrationBasket(swapBase, swaptionVol, 'MaturityStrikeByDeltaGamma')
printBasket(basket)

for basket_i in basket:
    ql.as_black_helper(basket_i).setPricingEngine(swaptionEngine)

print( "\nThe calibrated nominal is close to the exotics nominal."
       "\nThe expiries and maturity dates of the vanillas are the same"
       "\nas in the case above. The difference is the strike which"
       "\nis now equal to the exotics strike."
       "\nLet's see how this affects the exotics npv. The "
       "\nrecalibrated model is:" )

gsr.calibrateVolatilitiesIterative(basket, method, ec)
printModelCalibration(basket, gsr.volatility())
print("Bermudan swaption NPV (deal strike calibrated GSR) = %f" % swaption.NPV())

# Calibrated to amortizing nominial schedule
print("\nWe can do more complicated things, let's e.g. modify the"
       "\nnominal schedule to be linear amortizing and see what"
       "\nthe effect on the generated calibration basket is:")

for i in range(0,len(fixedSchedule)-1):
    tmp = 1 - i/ (len(fixedSchedule)-1)
    fixedNominal[i]        = tmp
    floatingNominal[i*2]   = tmp
    floatingNominal[i*2+1] = tmp

underlying2 = ql.NonstandardSwap(ql.VanillaSwap.Payer,
                            fixedNominal, floatingNominal, fixedSchedule, strike,
                            ql.Thirty360(), floatSchedule,
                            euribor6m, gearing, spread, ql.Actual360(), False, False, ql.ModifiedFollowing)

swaption2 = ql.NonstandardSwaption(underlying2,exercise,ql.Settlement.Physical)

swaption2.setPricingEngine(nonstandardSwaptionEngine)
basket = swaption2.calibrationBasket(swapBase, swaptionVol, 'MaturityStrikeByDeltaGamma')

printBasket(basket)

print("\nThe notional is weighted over the underlying exercised "
     "\ninto and the maturity is adjusted downwards. The rate"
     "\non the other hand is not affected." )

############## Pricing Bond features###################
print("\nYou can also price exotic bond's features. If you have e.g. a"
       "\nbermudan callable fixed bond you can set up the call right "
       "\nas a swaption to enter into a one leg swap with notional"
       "\nreimbursement at maturity."
       "\nThe exercise should then be written as a rebated exercise"
       "\npaying the notional in case of exercise."
       "\nThe calibration basket looks like this:")

fixedNominal2    = [1]*(len(fixedSchedule)-1)
floatingNominal2 = [0]*(len(floatSchedule)-1) #null the second leg

underlying3 = ql.NonstandardSwap(ql.VanillaSwap.Receiver,
                            fixedNominal2, floatingNominal2, fixedSchedule, strike,
                            ql.Thirty360(), floatSchedule,
                            euribor6m, gearing, spread, ql.Actual360(), False, True, ql.ModifiedFollowing)

rebateAmount = [-1 for x in range(0,len(exerciseDates))]
exercise2 = ql.RebatedExercise(exercise, rebateAmount, 2, ql.TARGET())
swaption3 = ql.NonstandardSwaption(underlying3,exercise2,ql.Settlement.Physical)

oas0 = ql.SimpleQuote(0)
oas100 = ql.SimpleQuote(0.01)
oas = ql.RelinkableQuoteHandle(oas0)

nonstandardSwaptionEngine2 = ql.Gaussian1dNonstandardSwaptionEngine(gsr, 64, 7.0, True, False, oas, t0_curve) # Change discounting to 6m

swaption3.setPricingEngine(nonstandardSwaptionEngine2)
basket = swaption3.calibrationBasket(swapBase, swaptionVol, 'MaturityStrikeByDeltaGamma')

printBasket(basket)

print( "\nNote that nominals are not exactly 1.0 here. This is"
        "\nbecause we do our bond discounting on 6m level while"
        "\nthe swaptions are still discounted on OIS level."
        "\n(You can try this by changing the OIS level to the "
        "\n6m level, which will produce nominals near 1.0)."
        "\nThe npv of the call right is (after recalibrating the model)")

for basket_i in basket:
    ql.as_black_helper(basket_i).setPricingEngine(swaptionEngine)

gsr.calibrateVolatilitiesIterative(basket, method, ec)

npv3 = swaption3.NPV()
print("Bond's bermudan call right npv = %f" % npv3)
print("\nUp to now, no credit spread is included in the pricing."
       "\nWe can do so by specifying an oas in the pricing engine."
       "\nLet's set the spread level to 100bp and regenerate"
       "\nthe calibration basket." )

oas.linkTo(oas100)
basket = swaption3.calibrationBasket(swapBase, swaptionVol, 'MaturityStrikeByDeltaGamma')
printBasket(basket)

print("The adjusted basket takes the credit spread into account"
      "\nThis is consistent to a hedge where you would have a"
      "\nmargin on the float leg around 100bp,too.")

for basket_i in basket:
    ql.as_black_helper(basket_i).setPricingEngine(swaptionEngine)

gsr.calibrateVolatilitiesIterative(basket, method, ec)

npv4 = swaption3.NPV()
print("Bond's bermudan call right npv (oas = 100bp) = %f" % npv4)

############### CMS10Y vs Euribor 6M swaption ###############
print("\nThe next instrument we look at is a CMS 10Y vs Euribor "
       "\n6M swaption. The maturity is again 10 years and the option"
       "\nis exercisable on a yearly basis" )

CMSNominal     = [1]*(len(fixedSchedule)-1)
CMSgearing     = [1]*(len(fixedSchedule)-1)
CMSspread      = [0]*(len(fixedSchedule)-1)
EuriborNominal = [1]*(len(floatSchedule)-1)
Euriborgearing = [1]*(len(floatSchedule)-1)
Euriborspread  = [0.001]*(len(floatSchedule)-1)
underlying4 = ql.FloatFloatSwap(ql.VanillaSwap.Payer,
                                CMSNominal, EuriborNominal,
                                fixedSchedule, swapBase, ql.Thirty360(),
                                floatSchedule, euribor6m, ql.Actual360(),
                                False, False, CMSgearing, CMSspread, [], [],
                                Euriborgearing, Euriborspread)

swaption4 = ql.FloatFloatSwaption(underlying4, exercise)
floatSwaptionEngine = ql.Gaussian1dFloatFloatSwaptionEngine(gsr, 64, 7.0, True, False, ql.QuoteHandle(ql.SimpleQuote(0)), t0_Ois, True)
swaption4.setPricingEngine(floatSwaptionEngine)

print("\nSince the underlying is quite exotic already, we start with"
      "\npricing this using the LinearTsrPricer for CMS coupon "
      "estimation" )

leg0 = underlying4.leg(0)
leg1 = underlying4.leg(1)
reversionQuote = ql.QuoteHandle(ql.SimpleQuote(0.01))
swaptionVolHandle = ql.SwaptionVolatilityStructureHandle(swaptionVol)
cmsPricer = ql.LinearTsrPricer(swaptionVolHandle, reversionQuote)
#cmsPricer = ql.CmsCouponPricer(ql.LinearTsrPricer(swaptionVolHandle, reversionQuote))
iborPricer = ql.BlackIborCouponPricer()

ql.setCouponPricer(leg0, cmsPricer)
ql.setCouponPricer(leg1, iborPricer)

swapPricer = ql.DiscountingSwapEngine(t0_Ois)
underlying4.setPricingEngine(swapPricer)

npv5 = underlying4.NPV()

print("Underlying CMS Swap NPV = %f" % npv5)
print("Underlying CMS Leg  NPV = %f" % underlying4.legNPV(0))
print("Underlying Euribor  NPV = %f" % underlying4.legNPV(1))

print("\nWe generate a naive calibration basket and calibrate "
      "\nthe GSR model to it:" )

# Recalibration to 4$ Strike swaption
basket = swaption4.calibrationBasket(swapBase, swaptionVol, 'Naive')

for basket_i in basket:
    ql.as_black_helper(basket_i).setPricingEngine(swaptionEngine)

gsr.calibrateVolatilitiesIterative(basket, method, ec)
printBasket(basket)
printModelCalibration(basket, gsr.volatility())

print("\nThe npv of the bermudan swaption is "
      "\nFloat swaption NPV (GSR) = %f" % swaption4.NPV())

print("\nIn this case it is also interesting to look at the "
                     "\nunderlying swap npv in the GSR model." )

print("\nFloat swap NPV (GSR) = %f" % swaption4.underlyingValue())

print("\nNot surprisingly, the underlying is priced differently"
      "\ncompared to the LinearTsrPricer, since a different"
      "\nsmile is implied by the GSR model."
      "\nThis is exactly where the Markov functional model"
      "\ncomes into play, because it can calibrate to any"
      "\ngiven underlying smile (as long as it is arbitrage"
      "\nfree). We try this now. Of course the usual use case"
      "\nis not to calibrate to a flat smile as in our simple"
      "\nexample, still it should be possible, of course..."
      )

markovStepDates = exerciseDates
cmsFixingDates = markovStepDates
markovSimgas = [0.01]* (len(markovStepDates)+1)
tenors = [ql.Period('10Y')]*len(cmsFixingDates)
markov = ql.MarkovFunctional(t0_curve, reversionQuote.value(), markovStepDates, markovSimgas, swaptionVolHandle,
                             cmsFixingDates, tenors, swapBase, 16)

swaptionEngineMarkov = ql.Gaussian1dSwaptionEngine(markov, 8, 5.0, True,
                                                   False, t0_Ois)

floatEngineMarkov = ql.Gaussian1dFloatFloatSwaptionEngine(markov, 16, 7.0, True, False, ql.QuoteHandle(ql.SimpleQuote(0)), t0_Ois, True)

swaption4.setPricingEngine(floatEngineMarkov)
npv7 = swaption4.NPV()

print("\nThe option npv is the markov model is:"
      "\nFloat swaption NPV (Markov) = %f" % npv7)
print("\nThis is not too far from the GSR price."
      "\nMore interesting is the question how well the Markov"
      "\nmodel did its job to match our input smile. For this"
      "\nwe look at the underlying npv under the Markov model"
      "\nFloat swap NPV (Markov) = %f" % swaption4.underlyingValue())
print("\nThis is closer to our terminal swap rate model price."
     "\nA perfect match is not expected anyway, because the"
     "\ndynamics of the underlying rate in the linear"
     "\nmodel is different from the Markov model, of"
     "\ncourse."
     "\nThe Markov model can not only calibrate to the"
     "\nunderlying smile, but has at the same time a"
     "\nsigma function (similar to the GSR model) which"
     "\ncan be used to calibrate to a second instrument"
     "\nset. We do this here to calibrate to our coterminal"
     "\nATM swaptions from above."
     "\nThis is a computationally demanding task, so"
     "\ndepending on your machine, this may take a"
     "\nwhile now...")

for basket_i in basket:
    ql.as_black_helper(basket_i).setPricingEngine(swaptionEngineMarkov)

markov.calibrate(basket, method, ec)
printModelCalibration(basket, markov.volatility())

print("\nNow let's have a look again at the underlying pricing."
      "\nIt shouldn't have changed much, because the underlying"
      "\nsmile is still matched." )

npv8 = swaption4.underlyingValue()
print("\nFloat swap NPV (Markov) = %f" % npv8)

print("\nThis is close to the previous value as expected."
      "\nAs a final remark we note that the calibration to"
      "\ncoterminal swaptions is not particularly reasonable"
      "\nhere, because the european call rights are not"
      "\nwell represented by these swaptions."
      "\nSecondly, our CMS swaption is sensitive to the"
      "\ncorrelation between the 10y swap rate and the"
      "\nEuribor 6M rate. Since the Markov model is one factor"
      "\nit will most probably underestimate the market value"
      "\nby construction."
      "\nThat was it. Thank you for running this demo. Bye."
      )
