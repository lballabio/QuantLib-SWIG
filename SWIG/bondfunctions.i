
/*
 Copyright (C) 2013 Simon Shakeshaft
 Copyright (C) 2016 Gouthaman Balaraman

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

#ifndef quantlib_bond_functions_i
#define quantlib_bond_functions_i

%include bonds.i
%include common.i
%include types.i
%include daycounters.i


%{
using QuantLib::BondFunctions;
%}

class BondFunctions {
    #if defined(SWIGPYTHON)
    %rename(bondYield) yield;
    #endif
  public:
    %extend {
        static Date startDate(const ext::shared_ptr<Bond>& bond) {
            return QuantLib::BondFunctions::startDate(
                    *(ext::dynamic_pointer_cast<Bond>(bond)));
        }
        static Date maturityDate(const ext::shared_ptr<Bond>& bond) {
            return QuantLib::BondFunctions::maturityDate(
                    *(ext::dynamic_pointer_cast<Bond>(bond)));
        }
        static bool isTradable(const ext::shared_ptr<Bond>& bond,
                               Date settlementDate = Date()) {
            return QuantLib::BondFunctions::isTradable(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    settlementDate);
        }
        static Date previousCashFlowDate(const ext::shared_ptr<Bond>& bond,
                                         Date refDate = Date()) {
            return QuantLib::BondFunctions::previousCashFlowDate(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    refDate);
        }
        static Date nextCashFlowDate(const ext::shared_ptr<Bond>& bond,
                                     Date refDate = Date()) {
            return QuantLib::BondFunctions::nextCashFlowDate(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    refDate);
        }
        static Real previousCashFlowAmount(const ext::shared_ptr<Bond>& bond,
                                           Date refDate = Date()) {
            return QuantLib::BondFunctions::previousCashFlowAmount(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    refDate);
        }
        static Real nextCashFlowAmount(const ext::shared_ptr<Bond>& bond,
                                       Date refDate = Date()) {
            return QuantLib::BondFunctions::nextCashFlowAmount(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    refDate);
        }
        static Rate previousCouponRate(const ext::shared_ptr<Bond>& bond,
                                       Date settlementDate = Date()) {
            return QuantLib::BondFunctions::previousCouponRate(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    settlementDate);
        }
        static Rate nextCouponRate(const ext::shared_ptr<Bond>& bond,
                                   Date settlementDate = Date()) {
            return QuantLib::BondFunctions::nextCouponRate(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    settlementDate);
        }
        static Date accrualStartDate(const ext::shared_ptr<Bond>& bond,
                                     Date settlementDate = Date()) {
            return QuantLib::BondFunctions::accrualStartDate(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    settlementDate);
        }
        static Date accrualEndDate(const ext::shared_ptr<Bond>& bond,
                                   Date settlementDate = Date()) {
            return QuantLib::BondFunctions::accrualEndDate(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    settlementDate);
        }
        static Time accrualPeriod(const ext::shared_ptr<Bond>& bond,
                                  Date settlementDate = Date()) {
            return QuantLib::BondFunctions::accrualPeriod(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    settlementDate);
        }
        static BigInteger accrualDays(const ext::shared_ptr<Bond>& bond,
                                      Date settlementDate = Date()) {
            return QuantLib::BondFunctions::accrualDays(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    settlementDate);
        }
        static Time accruedPeriod(const ext::shared_ptr<Bond>& bond,
                                  Date settlementDate = Date()) {
            return QuantLib::BondFunctions::accruedPeriod(
                *(ext::dynamic_pointer_cast<Bond>(bond)),
                settlementDate);
        }
        static BigInteger accruedDays(const ext::shared_ptr<Bond>& bond,
                                      Date settlementDate = Date()) {
            return QuantLib::BondFunctions::accruedDays(
                *(ext::dynamic_pointer_cast<Bond>(bond)),
                settlementDate);
        }
        static Real accruedAmount(const ext::shared_ptr<Bond>& bond,
                                  Date settlementDate = Date()){

            return QuantLib::BondFunctions::accruedAmount(
                *(ext::dynamic_pointer_cast<Bond>(bond)),
                settlementDate);
        }

        static Real cleanPrice(
                   const ext::shared_ptr<Bond>& bond,
                   const ext::shared_ptr<YieldTermStructure>& discountCurve,
                   Date settlementDate = Date()) {
            return QuantLib::BondFunctions::cleanPrice(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    *discountCurve,
                    settlementDate);
        }
        static Real bps(
                   const ext::shared_ptr<Bond>& bond,
                   const ext::shared_ptr<YieldTermStructure>& discountCurve,
                   Date settlementDate = Date()) {
            return QuantLib::BondFunctions::bps(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    *discountCurve,
                    settlementDate);
        }
        static Rate atmRate(
                   const ext::shared_ptr<Bond>& bond,
                   const ext::shared_ptr<YieldTermStructure>& discountCurve,
                   Date settlementDate = Date(),
                   Real cleanPrice = Null<Real>()) {
            return QuantLib::BondFunctions::atmRate(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    *discountCurve,
                    settlementDate,
                    cleanPrice);
        }
        static Real cleanPrice(const ext::shared_ptr<Bond>& bond,
                               const InterestRate& yield,
                               Date settlementDate = Date()) {
            return QuantLib::BondFunctions::cleanPrice(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    yield,
                    settlementDate);
        }
        static Real cleanPrice(const ext::shared_ptr<Bond>& bond,
                               Rate yield,
                               const DayCounter& dayCounter,
                               Compounding compounding,
                               Frequency frequency,
                               Date settlementDate = Date()) {
            return QuantLib::BondFunctions::cleanPrice(
                    *(ext::dynamic_pointer_cast<Bond>(bond)),
                    yield,
                    dayCounter,
                    compounding,
                    frequency,
                    settlementDate);
        }
        static Real bps(const ext::shared_ptr<Bond>& bond,
                        const InterestRate& yield,
                        Date settlementDate = Date()) {
            return QuantLib::BondFunctions::bps(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        yield,
                        settlementDate);
        }
        static Real bps(const ext::shared_ptr<Bond>& bond,
                        Rate yield,
                        const DayCounter& dayCounter,
                        Compounding compounding,
                        Frequency frequency,
                        Date settlementDate = Date()) {
            return QuantLib::BondFunctions::bps(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        yield,
                        dayCounter,
                        compounding,
                        frequency,
                        settlementDate);
        }
        static Rate yield(const ext::shared_ptr<Bond>& bond,
                          Real cleanPrice,
                          const DayCounter& dayCounter,
                          Compounding compounding,
                          Frequency frequency,
                          Date settlementDate = Date(),
                          Real accuracy = 1.0e-10,
                          Size maxIterations = 100,
                          Rate guess = 0.05) {
            return QuantLib::BondFunctions::yield(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        cleanPrice,
                        dayCounter,
                        compounding,
                        frequency,
                        settlementDate,
                        accuracy,
                        maxIterations,
                        guess);
        }

        %define DefineYieldFunctionSolver(SolverType)
        static Rate yield ## SolverType(SolverType solver,
                                         const ext::shared_ptr<Bond>& bond,
                                         Real cleanPrice,
                                         const DayCounter& dayCounter,
                                         Compounding compounding,
                                         Frequency frequency,
                                         Date settlementDate = Date(),
                                         Real accuracy = 1.0e-10,
                                         Rate guess = 0.05) {
            return QuantLib::BondFunctions::yield<SolverType>(
                        solver,
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        cleanPrice,
                        dayCounter,
                        compounding,
                        frequency,
                        settlementDate,
                        accuracy,
                        guess);
        }
        %enddef

        // See optimizers.i for solver definitions.
        DefineYieldFunctionSolver(Brent);
        DefineYieldFunctionSolver(Bisection);
        DefineYieldFunctionSolver(FalsePosition);
        DefineYieldFunctionSolver(Ridder);
        DefineYieldFunctionSolver(Secant);
        #if defined(SWIGPYTHON)
        DefineYieldFunctionSolver(Newton);
        DefineYieldFunctionSolver(NewtonSafe);
        #endif

        static Time duration(const ext::shared_ptr<Bond>& bond,
                             const InterestRate& yield,
                             Duration::Type type = Duration::Modified,
                             Date settlementDate = Date() ) {
            return QuantLib::BondFunctions::duration(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        yield,
                        type,
                        settlementDate);
        }
        static Time duration(const ext::shared_ptr<Bond>& bond,
                        Rate yield,
                        const DayCounter& dayCounter,
                        Compounding compounding,
                        Frequency frequency,
                        Duration::Type type = Duration::Modified,
                        Date settlementDate = Date()) {
            return QuantLib::BondFunctions::duration(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        yield,
                        dayCounter,
                        compounding,
                        frequency,
                        type,
                        settlementDate);
        }
        static Real convexity(const ext::shared_ptr<Bond>& bond,
                              const InterestRate& yield,
                              Date settlementDate = Date()) {
            return QuantLib::BondFunctions::convexity(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        yield,
                        settlementDate);
        }
        static Real convexity(const ext::shared_ptr<Bond>& bond,
                              Rate yield,
                              const DayCounter& dayCounter,
                              Compounding compounding,
                              Frequency frequency,
                              Date settlementDate = Date()) {
            return QuantLib::BondFunctions::convexity(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        yield,
                        dayCounter,
                        compounding,
                        frequency,
                        settlementDate);
        }
        static Real basisPointValue(const ext::shared_ptr<Bond>& bond,
                                    const InterestRate& yield,
                                    Date settlementDate = Date()) {
            return QuantLib::BondFunctions::basisPointValue(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        yield,
                        settlementDate);
        }
        static Real basisPointValue(const ext::shared_ptr<Bond>& bond,
                                    Rate yield,
                                    const DayCounter& dayCounter,
                                    Compounding compounding,
                                    Frequency frequency,
                                    Date settlementDate = Date()) {
            return QuantLib::BondFunctions::basisPointValue(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        yield,
                        dayCounter,
                        compounding,
                        frequency,
                        settlementDate);
        }
        static Real yieldValueBasisPoint(const ext::shared_ptr<Bond>& bond,
                                         const InterestRate& yield,
                                         Date settlementDate = Date()) {
            return QuantLib::BondFunctions::yieldValueBasisPoint(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        yield,
                        settlementDate);
        }
        static Real yieldValueBasisPoint(const ext::shared_ptr<Bond>& bond,
                                         Rate yield,
                                         const DayCounter& dayCounter,
                                         Compounding compounding,
                                         Frequency frequency,
                                         Date settlementDate = Date()) {
            return QuantLib::BondFunctions::yieldValueBasisPoint(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        yield,
                        dayCounter,
                        compounding,
                        frequency,
                        settlementDate);
        }
        static Spread zSpread(const ext::shared_ptr<Bond>& bond,
                              Real cleanPrice,
                              const ext::shared_ptr<YieldTermStructure>& discountCurve,
                              const DayCounter& dayCounter,
                              Compounding compounding,
                              Frequency frequency,
                              Date settlementDate = Date(),
                              Real accuracy = 1.0e-10,
                              Size maxIterations = 100,
                              Rate guess = 0.0){
            return QuantLib::BondFunctions::zSpread(
                        *(ext::dynamic_pointer_cast<Bond>(bond)),
                        cleanPrice,
                        discountCurve,
                        dayCounter,
                        compounding,
                        frequency,
                        settlementDate,
                        accuracy,
                        maxIterations,
                        guess);

        }

    }
};


#endif
