/*
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
*/

#ifndef quantlib_slv_i
#define quantlib_slv_i

%include common.i
%include stl.i
%include vectors.i
%include volatilities.i
%include stochasticprocess.i
%include calibrationhelpers.i
%include fdm.i

%{
using QuantLib::HestonSLVProcess;
%}

%shared_ptr(HestonSLVProcess);
class HestonSLVProcess : public StochasticProcess {
  public:
    HestonSLVProcess(const boost::shared_ptr<HestonProcess>& hestonProcess,
                     const boost::shared_ptr<LocalVolTermStructure>& leverageFct);
};


%{
using QuantLib::HestonSLVMCModel;
using QuantLib::BrownianGeneratorFactory;
using QuantLib::MTBrownianGeneratorFactory;
using QuantLib::SobolBrownianGeneratorFactory;
%}

%shared_ptr(BrownianGeneratorFactory);
class BrownianGeneratorFactory {
  private:
    BrownianGeneratorFactory();
};


%shared_ptr(MTBrownianGeneratorFactory);
class MTBrownianGeneratorFactory : public BrownianGeneratorFactory {
  public:
    MTBrownianGeneratorFactory(unsigned long seed = 0);
};

class SobolBrownianGenerator {
  public:
    enum Ordering { Factors, Steps, Diagonal };
  private: 
    SobolBrownianGenerator();
};

%shared_ptr(SobolBrownianGeneratorFactory);
class SobolBrownianGeneratorFactory : public BrownianGeneratorFactory {
  public:
    SobolBrownianGeneratorFactory(
            SobolBrownianGenerator::Ordering ordering,
            unsigned long seed = 0,
            SobolRsg::DirectionIntegers directionIntegers = SobolRsg::Jaeckel);
};

class HestonSLVMCModel {
  public:
    %extend {
        HestonSLVMCModel(
           const boost::shared_ptr<LocalVolTermStructure>& localVol,
           const boost::shared_ptr<HestonModel>& model,
           const boost::shared_ptr<BrownianGeneratorFactory>& brownianGeneratorFactory,
           const Date& endDate,
           Size timeStepsPerYear = 365,
           Size nBins = 201,
           Size calibrationPaths = (1 << 15),
           const std::vector<Date>& mandatoryDates = std::vector<Date>()) {
            return new HestonSLVMCModel(
                Handle<LocalVolTermStructure>(localVol), Handle<HestonModel>(model),
                brownianGeneratorFactory, endDate, timeStepsPerYear,
                nBins, calibrationPaths, mandatoryDates);
        }
    }
    boost::shared_ptr<HestonProcess> hestonProcess() const;
    boost::shared_ptr<LocalVolTermStructure> localVol() const;
    boost::shared_ptr<LocalVolTermStructure> leverageFunction() const;
};


%{
using QuantLib::FdmHestonGreensFct;
using QuantLib::HestonSLVFDMModel;
using QuantLib::HestonSLVFokkerPlanckFdmParams;
%}

struct FdmHestonGreensFct {
    enum Algorithm { ZeroCorrelation, Gaussian, SemiAnalytical };
  private:
    FdmHestonGreensFct();
};

class HestonSLVFokkerPlanckFdmParams {
  public:
    %extend {
        HestonSLVFokkerPlanckFdmParams(
            Size xGrid, Size vGrid, 
            Size tMaxStepsPerYear, Size tMinStepsPerYear,
            Real tStepNumberDecay,
            Size nRannacherTimeSteps,
            Size predictionCorretionSteps,
            Real x0Density, Real localVolEpsProb,
            Size maxIntegrationIterations,
            Real vLowerEps, Real vUpperEps, Real vMin,
            Real v0Density, Real vLowerBoundDensity, Real vUpperBoundDensity,
            Real leverageFctPropEps,
            FdmHestonGreensFct::Algorithm greensAlgorithm,
            FdmSquareRootFwdOp::TransformationType trafoType,
            FdmSchemeDesc schemeDesc) {
            
                const HestonSLVFokkerPlanckFdmParams params = {
                    xGrid, vGrid,
                    tMaxStepsPerYear, tMinStepsPerYear,
                    tStepNumberDecay,
                    nRannacherTimeSteps,
                    predictionCorretionSteps,
                    x0Density,
                    localVolEpsProb,
                    maxIntegrationIterations,
                    vLowerEps, vUpperEps, vMin,
                    v0Density, vLowerBoundDensity, vUpperBoundDensity,
                    leverageFctPropEps,
                    greensAlgorithm,
                    trafoType,
                    schemeDesc };

                return new HestonSLVFokkerPlanckFdmParams(params);
        }
    }
};

class HestonSLVFDMModel {
  public:
    %extend {
        HestonSLVFDMModel(
            const boost::shared_ptr<LocalVolTermStructure>& localVol,
            const boost::shared_ptr<HestonModel>& model,
            const Date& endDate,
            const HestonSLVFokkerPlanckFdmParams& params,
            const bool logging = false,
            const std::vector<Date>& mandatoryDates = std::vector<Date>()) {
            return new HestonSLVFDMModel(
                Handle<LocalVolTermStructure>(localVol), Handle<HestonModel>(model),
                endDate, params, logging, mandatoryDates);
        }
    }
    boost::shared_ptr<HestonProcess> hestonProcess() const;
    boost::shared_ptr<LocalVolTermStructure> localVol() const;
    boost::shared_ptr<LocalVolTermStructure> leverageFunction() const;
};


#endif
