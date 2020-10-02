
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
    static Date startDate(const Bond& bond);
    static Date maturityDate(const Bond& bond);
    static bool isTradable(const Bond& bond,
                           Date settlementDate = Date());
    static Date previousCashFlowDate(const Bond& bond,
                                     Date refDate = Date());
    static Date nextCashFlowDate(const Bond& bond,
                                 Date refDate = Date());
    static Real previousCashFlowAmount(const Bond& bond,
                                       Date refDate = Date());
    static Real nextCashFlowAmount(const Bond& bond,
                                   Date refDate = Date());
    static Rate previousCouponRate(const Bond& bond,
                                   Date settlementDate = Date());
    static Rate nextCouponRate(const Bond& bond,
                               Date settlementDate = Date());
    static Date accrualStartDate(const Bond& bond,
                                 Date settlementDate = Date());
    static Date accrualEndDate(const Bond& bond,
                               Date settlementDate = Date());
    static Time accrualPeriod(const Bond& bond,
                              Date settlementDate = Date());
    static BigInteger accrualDays(const Bond& bond,
                                  Date settlementDate = Date());
    static Time accruedPeriod(const Bond& bond,
                              Date settlementDate = Date());
    static BigInteger accruedDays(const Bond& bond,
                                  Date settlementDate = Date());
    static Real accruedAmount(const Bond& bond,
                              Date settlementDate = Date());

    static Real cleanPrice(const Bond& bond,
                           const YieldTermStructure& discountCurve,
                           Date settlementDate = Date());
    static Real bps(const Bond& bond,
                    const YieldTermStructure& discountCurve,
                    Date settlementDate = Date());
    static Rate atmRate(const Bond& bond,
                        const YieldTermStructure& discountCurve,
                        Date settlementDate = Date(),
                        Real cleanPrice = Null<Real>());

    static Real cleanPrice(const Bond& bond,
                           const InterestRate& yield,
                           Date settlementDate = Date());
    static Real cleanPrice(const Bond& bond,
                           Rate yield,
                           const DayCounter& dayCounter,
                           Compounding compounding,
                           Frequency frequency,
                           Date settlementDate = Date());
    static Real bps(const Bond& bond,
                    const InterestRate& yield,
                    Date settlementDate = Date());
    static Real bps(const Bond& bond,
                    Rate yield,
                    const DayCounter& dayCounter,
                    Compounding compounding,
                    Frequency frequency,
                    Date settlementDate = Date());
    static Rate yield(const Bond& bond,
                      Real cleanPrice,
                      const DayCounter& dayCounter,
                      Compounding compounding,
                      Frequency frequency,
                      Date settlementDate = Date(),
                      Real accuracy = 1.0e-10,
                      Size maxIterations = 100,
                      Rate guess = 0.05);

    static Time duration(const Bond& bond,
                         const InterestRate& yield,
                         Duration::Type type = Duration::Modified,
                         Date settlementDate = Date());
    static Time duration(const Bond& bond,
                         Rate yield,
                         const DayCounter& dayCounter,
                         Compounding compounding,
                         Frequency frequency,
                         Duration::Type type = Duration::Modified,
                         Date settlementDate = Date());
    static Real convexity(const Bond& bond,
                          const InterestRate& yield,
                          Date settlementDate = Date());
    static Real convexity(const Bond& bond,
                          Rate yield,
                          const DayCounter& dayCounter,
                          Compounding compounding,
                          Frequency frequency,
                          Date settlementDate = Date());
    static Real basisPointValue(const Bond& bond,
                                const InterestRate& yield,
                                Date settlementDate = Date());
    static Real basisPointValue(const Bond& bond,
                                Rate yield,
                                const DayCounter& dayCounter,
                                Compounding compounding,
                                Frequency frequency,
                                Date settlementDate = Date());
    static Real yieldValueBasisPoint(const Bond& bond,
                                     const InterestRate& yield,
                                     Date settlementDate = Date());
    static Real yieldValueBasisPoint(const Bond& bond,
                                     Rate yield,
                                     const DayCounter& dayCounter,
                                     Compounding compounding,
                                     Frequency frequency,
                                     Date settlementDate = Date());
    static Spread zSpread(const Bond& bond,
                          Real cleanPrice,
                          const ext::shared_ptr<YieldTermStructure>& discountCurve,
                          const DayCounter& dayCounter,
                          Compounding compounding,
                          Frequency frequency,
                          Date settlementDate = Date(),
                          Real accuracy = 1.0e-10,
                          Size maxIterations = 100,
                          Rate guess = 0.0);

    %extend {

        %define DefineYieldFunctionSolver(SolverType)
        static Rate yield ## SolverType(SolverType solver,
                                         const Bond& bond,
                                         Real cleanPrice,
                                         const DayCounter& dayCounter,
                                         Compounding compounding,
                                         Frequency frequency,
                                         Date settlementDate = Date(),
                                         Real accuracy = 1.0e-10,
                                         Rate guess = 0.05) {
            return QuantLib::BondFunctions::yield<SolverType>(
                        solver,
                        bond,
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
    }
};


#endif
