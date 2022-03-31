/*
 Copyright (C) 2014 StatPro Italia srl

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

#ifndef quantlib_fitted_bond_i
#define quantlib_fitted_bond_i

%include termstructures.i
%include interpolation.i
%include ratehelpers.i

%{
using QuantLib::FittedBondDiscountCurve;

typedef ext::shared_ptr<YieldTermStructure> FittedBondDiscountCurvePtr;
typedef QuantLib::FittedBondDiscountCurve::FittingMethod FittingMethod;

std::vector<ext::shared_ptr<BondHelper> > convert_bond_helpers(
                 const std::vector<ext::shared_ptr<RateHelper> >& helpers) {
    std::vector<ext::shared_ptr<BondHelper> > result(helpers.size());
    for (Size i=0; i<helpers.size(); ++i)
        result[i] = ext::dynamic_pointer_cast<BondHelper>(helpers[i]);
    return result;
}
%}

%shared_ptr(FittingMethod)
class FittingMethod {
  public:
    virtual ~FittingMethod() = 0;
    Size size() const;
    Array solution() const;
    Integer numberOfIterations() const;
    Real minimumCostValue() const;
    bool constrainAtZero() const;
    Array weights() const;
};

%shared_ptr(FittedBondDiscountCurve);
class FittedBondDiscountCurve : public YieldTermStructure {
  public:
    FittedBondDiscountCurve(
                   Natural settlementDays,
                   const Calendar& calendar,
                   const std::vector<ext::shared_ptr<BondHelper> >& helpers,
                   const DayCounter& dayCounter,
                   const FittingMethod& fittingMethod,
                   Real accuracy = 1.0e-10,
                   Size maxEvaluations = 10000,
                   const Array& guess = Array(),
                   Real simplexLambda = 1.0);
    FittedBondDiscountCurve(
                   const Date &referenceDate,
                   const std::vector<ext::shared_ptr<BondHelper> >& helpers,
                   const DayCounter& dayCounter,
                   const FittingMethod& fittingMethod,
                   Real accuracy = 1.0e-10,
                   Size maxEvaluations = 10000,
                   const Array &guess = Array(),
                   Real simplexLambda = 1.0);
    const FittingMethod& fitResults() const;
};


%{
using QuantLib::ExponentialSplinesFitting;
using QuantLib::NelsonSiegelFitting;
using QuantLib::SvenssonFitting;
using QuantLib::CubicBSplinesFitting;
using QuantLib::SimplePolynomialFitting;
using QuantLib::SpreadFittingMethod;
%}

%shared_ptr(ExponentialSplinesFitting)
class ExponentialSplinesFitting : public FittingMethod {
  public:
    ExponentialSplinesFitting(bool constrainAtZero = true,
                              const Array& weights = Array(),
                              const Array& l2 = Array(),
                              Real minCutoffTime = 0.0,
                              Real maxCutoffTime = QL_MAX_REAL,
                              Size numCoeffs = 9,
                              Real fixedKappa = Null<Real>());
};

%shared_ptr(NelsonSiegelFitting)
class NelsonSiegelFitting : public FittingMethod {
  public:
    NelsonSiegelFitting(const Array& weights = Array());
};

%shared_ptr(SvenssonFitting)
class SvenssonFitting : public FittingMethod {
  public:
    SvenssonFitting(const Array& weights = Array());
};

%shared_ptr(CubicBSplinesFitting)
class CubicBSplinesFitting : public FittingMethod {
  public:
    CubicBSplinesFitting(const std::vector<Time>& knotVector,
                         bool constrainAtZero = true,
                         const Array& weights = Array());
    Real basisFunction(Integer i, Time t);
};

%shared_ptr(SimplePolynomialFitting)
class SimplePolynomialFitting : public FittingMethod {
  public:
    #if defined(SWIGJAVA)
    SimplePolynomialFitting(Natural degree);
    #else
    SimplePolynomialFitting(Natural degree,
                            bool constrainAtZero = true,
                            const Array& weights = Array());
    #endif
};

%shared_ptr(SpreadFittingMethod)
class SpreadFittingMethod : public FittingMethod {
  public:
    SpreadFittingMethod(const ext::shared_ptr<FittingMethod>& method,
                        Handle<YieldTermStructure> discountCurve,
                        Real minCutoffTime = 0.0,
                        Real maxCutoffTime = QL_MAX_REAL);
};


#endif
