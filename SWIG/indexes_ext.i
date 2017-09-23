/*
 Copyright (C) 2017 Cheng Li
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

#ifndef quantlib_indexes_ext_i
#define quantlib_indexes_ext_i

%include indexes.i

export_xibor_instance(Shibor);

%{
    #include <qlext/indexes/repo.hpp>
    using QuantLib::Repo;
%}

export_xibor_instance(Repo);

#endif