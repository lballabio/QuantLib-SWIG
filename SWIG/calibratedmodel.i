/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2007, 2009 StatPro Italia srl
 Copyright (C) 2005 Dominic Thuillier
 Copyright (C) 2007 Luis Cota
 Copyright (C) 2016 Gouthaman Balaraman
 Copyright (C) 2016 Peter Caspers
 Copyright (C) 2018 Matthias Lungwitz

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

#ifndef quantlib_calibrated_model_i
#define quantlib_calibrated_model_i

%include termstructures.i
%include optimizers.i
%include linearalgebra.i
%include types.i
%include vectors.i

%{
using QuantLib::CalibrationHelper;
%}

// calibration helpers
%shared_ptr(CalibrationHelper)
class CalibrationHelper {
  public:
    Real calibrationError();
  private:
    CalibrationHelper();
};


// allow use of vectors of helpers
#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( ext::shared_ptr<CalibrationHelper> )
#endif
namespace std {
    %template(CalibrationHelperVector) vector<ext::shared_ptr<CalibrationHelper> >;
}

// the base class for calibrated models
%{
using QuantLib::CalibratedModel;
using QuantLib::TermStructureConsistentModel;
%}

%shared_ptr(CalibratedModel)
class CalibratedModel : public virtual Observable {
    #if defined(SWIGCSHARP)
    %rename("parameters") params;
    #endif
  public:
    Array params() const;
    virtual void calibrate(
        const std::vector<ext::shared_ptr<CalibrationHelper> >&,
        OptimizationMethod&, const EndCriteria &,
        const Constraint& constraint = Constraint(),
        const std::vector<Real>& weights = std::vector<Real>(),
        const std::vector<bool>& fixParameters = std::vector<bool>());

    void setParams(const Array& params);
    Real value(const Array& params,
               const std::vector<ext::shared_ptr<CalibrationHelper> >&);
    const ext::shared_ptr<Constraint>& constraint() const;
    EndCriteria::Type endCriteria() const;
    const Array& problemValues() const;
    Integer functionEvaluation() const;
  private:
    CalibratedModel();
};

%shared_ptr(TermStructureConsistentModel)
class TermStructureConsistentModel : public virtual Observable{
  public:
    const Handle<YieldTermStructure>& termStructure() const;
  private:
    TermStructureConsistentModel();
};

%template(CalibratedModelHandle) Handle<CalibratedModel>;
%template(RelinkableCalibratedModelHandle) RelinkableHandle<CalibratedModel>;



#endif
