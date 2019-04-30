"""
 Copyright (C) 2019 Klaus Spanderen

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
"""

import math
import unittest

from QuantLib import *

class SlvTest(unittest.TestCase):
    
    def setUp(self):
        self.todaysDate = Date(15,May,2019)
        Settings.instance().evaluationDate = self.todaysDate
        self.settlementDate = self.todaysDate + Period(2, Days)
        self.dc = Actual365Fixed()
        self.riskFreeRate = YieldTermStructureHandle(FlatForward(self.settlementDate, 0.05, self.dc))
        self.dividendYield = YieldTermStructureHandle(FlatForward(self.settlementDate, 0.025, self.dc))
        self.underlying = QuoteHandle(SimpleQuote(100.0))
    
    def tearDown(self):
        QuantLib.Settings.instance().setEvaluationDate(QuantLib.Date())


    def constVol(self, vol):
        return BlackVolTermStructureHandle(
            BlackConstantVol(self.settlementDate, TARGET(), vol, self.dc)) 

    def testSlvProcess(self):
        """ Testing HestonSLVProcess generation """
                
        hestonProcess = HestonProcess(
            self.riskFreeRate,
            self.riskFreeRate,
            self.underlying,
            0.1*0.1, 1.0, 0.25*0.25, 0.15, -0.75)
        
        localVol = LocalVolSurface(
            BlackVolTermStructureHandle(
                BlackConstantVol(self.settlementDate, TARGET(), 0.10, self.dc)),
            self.riskFreeRate,
            self.riskFreeRate,
            self.underlying)
            
        hestonSLVProcess = HestonSLVProcess(hestonProcess, localVol)

        
    def testSlvProcessAsBlackScholes(self):
        """ Testing HestonSLVProcess equal to Black-Scholes process """
        
        hestonProcess = HestonProcess(
            self.riskFreeRate,
            self.dividendYield,
            self.underlying,
            0.01, 1.0, 0.01, 1e-4, 0.0)
        
        localVol = LocalVolSurface(
            self.constVol(0.1),
            self.riskFreeRate,
            self.dividendYield,
            self.underlying)
            
        exercise = EuropeanExercise(self.todaysDate + Period(1, Years))
        payoff = PlainVanillaPayoff(Option.Call, self.underlying.value())                                    
        
        option = VanillaOption(payoff, exercise)
        
        hestonModel = HestonModel(hestonProcess)
        option.setPricingEngine(FdHestonVanillaEngine(hestonModel, 20, 100, 3))
        
        hestonNPV = option.NPV()
        
        option.setPricingEngine(
            AnalyticEuropeanEngine(BlackScholesMertonProcess(
                self.underlying, self.dividendYield, self.riskFreeRate, self.constVol(0.1)))
        )
                
        bsNPV = option.NPV()
        
        self.assertAlmostEqual(hestonNPV, bsNPV, 2, 
            msg="Unable to reproduce Heston vanilla option price with Black-Scholes process")
        
        option.setPricingEngine(
            FdHestonVanillaEngine(
                hestonModel, 20, 100, 3, 1, FdmSchemeDesc.Douglas(), 
                LocalVolSurface(
                    self.constVol(2.0),
                    self.riskFreeRate,
                    self.dividendYield,
                    self.underlying)))
        
        slvNPV = option.NPV()
        
        option.setPricingEngine(
            AnalyticEuropeanEngine(BlackScholesMertonProcess(
                self.underlying, self.dividendYield, self.riskFreeRate, self.constVol(0.20)))
        )
        
        bsNPV = option.NPV()
        
        self.assertAlmostEqual(slvNPV, bsNPV, 2, 
            msg="Unable to reproduce Heston plus constant local"
            "vol option price with Black-Scholes formula")


    def testSlvMonteCarloCalibration(self):
        """ Testing Monte-Carlo calibration of a HestonSLVProcess """
        
        localVol = LocalVolSurface(
            self.constVol(0.25),
            self.riskFreeRate,
            self.dividendYield,
            self.underlying)
        
        hestonProcess = HestonProcess(
            self.riskFreeRate,
            self.dividendYield,
            self.underlying,
            0.1*0.1, 5.0, 0.25*0.25, 0.25, -0.75)

        hestonModel = HestonModel(hestonProcess)

        exerciseDate = self.todaysDate + Period(1, Months)
        
        exercise = EuropeanExercise(exerciseDate)
        payoff = PlainVanillaPayoff(Option.Call, 1.1*self.underlying.value())                                    
        
        option = VanillaOption(payoff, exercise)
        
        option.setPricingEngine(
            AnalyticEuropeanEngine(
                BlackScholesMertonProcess(
                    self.underlying,
                    self.dividendYield,
                    self.riskFreeRate,
                    self.constVol(0.25))))
        
        bsNPV = option.NPV()
        
        option.setPricingEngine(
            FdHestonVanillaEngine(
                hestonModel, 25, 100, 50, 0, 
                FdmSchemeDesc.Hundsdorfer(), 
                HestonSLVMCModel(
                    localVol,
                    hestonModel,
                    MTBrownianGeneratorFactory(1234),
                    exerciseDate, 91).leverageFunction()))
                 
        slvNPV = option.NPV()

        self.assertAlmostEqual(slvNPV, bsNPV, 2, 
            msg="Unable to reproduce HestonSLV option price"
            " with Black-Scholes formula based on MC calibration.")
        
        fdmParams = HestonSLVFokkerPlanckFdmParams(
            201, 401, 200, 30, 2.0, 0, 2,
            0.1, 1e-4, 10000,
            1e-5, 1e-5, 0.0000025, 1.0, 0.1, 0.9, 1e-5,
            FdmHestonGreensFct.Gaussian,
            FdmSquareRootFwdOp.Log,
            FdmSchemeDesc.Hundsdorfer()
        )
        
        option.setPricingEngine(
            FdHestonVanillaEngine(
                hestonModel, 25, 100, 50, 0, 
                FdmSchemeDesc.Hundsdorfer(), 
                HestonSLVFDMModel(
                    localVol,
                    hestonModel,
                    exerciseDate,
                    fdmParams).leverageFunction()))
                 
        slvNPV = option.NPV()
        
        self.assertAlmostEqual(slvNPV, bsNPV, 2, 
            msg="Unable to reproduce HestonSLV option price"
            " with Black-Scholes formula based on FDM calibration.")
        
if __name__ == '__main__':
    print('testing QuantLib ' + QuantLib.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(SlvTest,'test'))
    unittest.TextTestRunner(verbosity=2).run(suite)
   
        
