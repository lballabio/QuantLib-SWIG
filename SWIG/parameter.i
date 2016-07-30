
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 StatPro Italia srl
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


#ifndef quantlib_parameter_i
#define quantlib_parameter_i

%include types.i
%include optimizers.i

%{
using QuantLib::Parameter;
using QuantLib::ConstantParameter;
using QuantLib::NullParameter;
using QuantLib::PiecewiseConstantParameter;
%}

class Parameter {
    public:
        Parameter();
        const Array& params() const;
        void setParam(Size i, Real x);
        bool testParams(const Array& params) const; 
        Size size() const ;
        Real operator()(Time t) const;
        const Constraint& constraint() const ;
};

//! Standard constant parameter \f$ a(t) = a \f$
class ConstantParameter : public Parameter {
    public: 
        ConstantParameter(const Constraint& constraint);
        ConstantParameter(Real value, const Constraint& constraint);

};

//! %Parameter which is always zero \f$ a(t) = 0 \f$
class NullParameter : public Parameter {
public:
    NullParameter();
};

//! Piecewise-constant parameter
/*! \f$ a(t) = a_i if t_{i-1} \geq t < t_i \f$.
    This kind of parameter is usually used to enhance the fitting of a
    model
*/
class PiecewiseConstantParameter : public Parameter {
public:
     PiecewiseConstantParameter(const std::vector<Time>& times,
                                const Constraint& constraint=QuantLib::NoConstraint());
    
};



#endif
