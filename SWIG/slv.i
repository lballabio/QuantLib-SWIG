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

%{
using QuantLib::HestonSLVProcess;
typedef boost::shared_ptr<StochasticProcess> HestonSLVProcessPtr;
%}

%rename(HestonSLVProcess) HestonSLVProcessPtr;
class HestonSLVProcessPtr : public boost::shared_ptr<StochasticProcess> {
  public:
    %extend {
      HestonSLVProcessPtr(const HestonProcessPtr& hestonProcess,
                          const boost::shared_ptr<LocalVolTermStructure>& leverageFct) {
        
        boost::shared_ptr<HestonProcess> hProcess =
            boost::dynamic_pointer_cast<HestonProcess>(hestonProcess);
        QL_REQUIRE(hProcess, "Heston process required");
        
        return new HestonSLVProcessPtr(
          new HestonSLVProcess(hProcess, leverageFct));
      }
    }
};


%{
using QuantLib::HestonProcess;
using QuantLib::HestonSLVMCModel;
using QuantLib::BrownianGeneratorFactory;
using QuantLib::MTBrownianGeneratorFactory;
using QuantLib::SobolBrownianGeneratorFactory;

typedef boost::shared_ptr<BrownianGeneratorFactory> MTBrownianGeneratorFactoryPtr;
typedef boost::shared_ptr<BrownianGeneratorFactory> SobolBrownianGeneratorFactoryPtr;
%}

%template(BrownianGeneratorFactory) boost::shared_ptr<BrownianGeneratorFactory>;

%rename(MTBrownianGeneratorFactory) MTBrownianGeneratorFactoryPtr;
class MTBrownianGeneratorFactoryPtr : public boost::shared_ptr<BrownianGeneratorFactory> {
  public:
    %extend {
        MTBrownianGeneratorFactoryPtr(unsigned long seed = 0) {
            return new MTBrownianGeneratorFactoryPtr(new MTBrownianGeneratorFactory(seed));
        }
    }
};

class SobolBrownianGenerator {
  public:
    enum Ordering { Factors, Steps, Diagonal };
    
  private: 
    SobolBrownianGenerator();
};

%rename(SobolBrownianGeneratorFactory) SobolBrownianGeneratorFactoryPtr;
class SobolBrownianGeneratorFactoryPtr : public boost::shared_ptr<BrownianGeneratorFactory> {
  public:
    %extend {
        SobolBrownianGeneratorFactoryPtr(
            SobolBrownianGenerator::Ordering ordering,
            unsigned long seed = 0,
            SobolRsg::DirectionIntegers directionIntegers = SobolRsg::Jaeckel) {
            
            return new SobolBrownianGeneratorFactoryPtr(
                new SobolBrownianGeneratorFactory(ordering, seed));
        }
    }
};

class HestonSLVMCModel {
  public:
    %extend {
        HestonSLVMCModel(
           const boost::shared_ptr<LocalVolTermStructure>& localVol,
           const HestonModelPtr& model,
           const boost::shared_ptr<BrownianGeneratorFactory>& brownianGeneratorFactory,
           const Date& endDate,
           Size timeStepsPerYear = 365,
           Size nBins = 201,
           Size calibrationPaths = (1 << 15),
           const std::vector<Date>& mandatoryDates = std::vector<Date>()) {
           
           boost::shared_ptr<HestonModel> hModel =
               boost::dynamic_pointer_cast<HestonModel>(model);
           QL_REQUIRE(hModel, "Heston model required");

           return new HestonSLVMCModel(
               Handle<LocalVolTermStructure>(localVol), Handle<HestonModel>(hModel), 
               brownianGeneratorFactory, endDate, timeStepsPerYear, nBins,
               calibrationPaths, mandatoryDates);
        }
    }    
    boost::shared_ptr<HestonProcess> hestonProcess() const;
    
    boost::shared_ptr<LocalVolTermStructure> localVol() const;
    
    boost::shared_ptr<LocalVolTermStructure> leverageFunction() const;
};

%{
using QuantLib::FdmSquareRootFwdOp;
using QuantLib::FdmHestonGreensFct;
using QuantLib::HestonSLVFDMModel;
using QuantLib::HestonSLVFokkerPlanckFdmParams;

typedef boost::shared_ptr<HestonSLVFokkerPlanckFdmParams> HestonSLVFokkerPlanckFdmParamsPtr;
%}

struct FdmSquareRootFwdOp {
    enum TransformationType { Plain, Power, Log };
  private:
    FdmSquareRootFwdOp();
};

struct FdmHestonGreensFct {
    enum Algorithm { ZeroCorrelation, Gaussian, SemiAnalytical };
  private:
    FdmHestonGreensFct();
};

%template(_HestonSLVFokkerPlanckFdmParams) boost::shared_ptr<HestonSLVFokkerPlanckFdmParams>;

%rename(HestonSLVFokkerPlanckFdmParams) HestonSLVFokkerPlanckFdmParamsPtr;
class HestonSLVFokkerPlanckFdmParamsPtr 
    : public boost::shared_ptr<HestonSLVFokkerPlanckFdmParams> {
  public:
    %extend {
        HestonSLVFokkerPlanckFdmParamsPtr(
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
                    
                return new HestonSLVFokkerPlanckFdmParamsPtr(
                    new HestonSLVFokkerPlanckFdmParams(params));
            }
            
        private:
            HestonSLVFokkerPlanckFdmParamsPtr cstr(
                const boost::shared_ptr<HestonSLVFokkerPlanckFdmParams>& params) : params_(params) {}
                
            boost::shared_ptr<HestonSLVFokkerPlanckFdmParams> params_;
    }
};

class HestonSLVFDMModel {
  public:
    %extend {
        HestonSLVFDMModel(
            const boost::shared_ptr<LocalVolTermStructure>& localVol,
            const HestonModelPtr& model,
            const Date& endDate,
            const HestonSLVFokkerPlanckFdmParamsPtr& params,
            const bool logging = false,
            const std::vector<Date>& mandatoryDates = std::vector<Date>()) {
            
            boost::shared_ptr<HestonModel> hModel =
                boost::dynamic_pointer_cast<HestonModel>(model);
            QL_REQUIRE(hModel, "Heston model required");

            
            return new HestonSLVFDMModel(
                Handle<LocalVolTermStructure>(localVol), Handle<HestonModel>(hModel), 
                endDate, *(params.get()), logging, mandatoryDates);            
        }
    }
    boost::shared_ptr<HestonProcess> hestonProcess() const;
    
    boost::shared_ptr<LocalVolTermStructure> localVol() const;
    
    boost::shared_ptr<LocalVolTermStructure> leverageFunction() const;
};


%{
using QuantLib::FdHestonBarrierEngine;
typedef boost::shared_ptr<PricingEngine> FdHestonBarrierEnginePtr;
%}

%rename(FdHestonBarrierEngine) FdHestonBarrierEnginePtr;
class FdHestonBarrierEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdHestonBarrierEnginePtr(
            const HestonModelPtr& model,
            Size tGrid = 100, Size xGrid = 100,
            Size vGrid = 50, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
            const boost::shared_ptr<LocalVolTermStructure>& leverageFct
                = boost::shared_ptr<LocalVolTermStructure>()) {
            
            boost::shared_ptr<HestonModel> hModel =
                 boost::dynamic_pointer_cast<HestonModel>(model);
            QL_REQUIRE(hModel, "Heston model required");
            return new FdHestonBarrierEnginePtr(
                new FdHestonBarrierEngine(hModel, tGrid, xGrid,
                    vGrid, dampingSteps, schemeDesc, leverageFct));
        }
    }
};


%{
using QuantLib::FdHestonDoubleBarrierEngine;
typedef boost::shared_ptr<PricingEngine> FdHestonDoubleBarrierEnginePtr;
%}

%rename(FdHestonDoubleBarrierEngine) FdHestonDoubleBarrierEnginePtr;
class FdHestonDoubleBarrierEnginePtr : public boost::shared_ptr<PricingEngine> {
  public:
    %extend {
        FdHestonDoubleBarrierEnginePtr(
            const HestonModelPtr& model,
            Size tGrid = 100, Size xGrid = 100,
            Size vGrid = 50, Size dampingSteps = 0,
            const FdmSchemeDesc& schemeDesc = FdmSchemeDesc::Hundsdorfer(),
            const boost::shared_ptr<LocalVolTermStructure>& leverageFct
                = boost::shared_ptr<LocalVolTermStructure>()) {
            
            boost::shared_ptr<HestonModel> hModel =
                 boost::dynamic_pointer_cast<HestonModel>(model);
            QL_REQUIRE(hModel, "Heston model required");
            return new FdHestonDoubleBarrierEnginePtr(
                new FdHestonDoubleBarrierEngine(hModel, tGrid, xGrid,
                    vGrid, dampingSteps, schemeDesc, leverageFct));
        }
    }
};

#endif
