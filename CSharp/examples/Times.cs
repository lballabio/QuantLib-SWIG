/*
 Copyright (C) 2021 Ralf Konrad Eckel

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

using System;
using System.Collections.Generic;
using System.Diagnostics;
using ql = QuantLib;

namespace TimesTest
{
    class Times
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main(string[] args)
        {
            var tenor03M = new ql.Period("03M");
            var tenor06M = new ql.Period("06M");
            var tenor01Y = new ql.Period("01Y");
            var tenor02Y = new ql.Period("02Y");

            var periods = new List<ql.Period>() { tenor01Y, tenor02Y, tenor06M, tenor03M };
            periods.Sort();

            Debug.Assert(periods[0] == tenor03M);
            Debug.Assert(periods[1] == tenor06M);
            Debug.Assert(periods[2] == tenor01Y);
            Debug.Assert(periods[3] == tenor02Y);

            var tenor12M = new ql.Period("12M");

            Debug.Assert(tenor03M.CompareTo(tenor12M) < 0);
            Debug.Assert(tenor06M.CompareTo(tenor12M) < 0);
            Debug.Assert(tenor01Y.CompareTo(tenor12M) == 0);
            Debug.Assert(tenor02Y.CompareTo(tenor12M) > 0);

            Debug.Assert(tenor03M != tenor12M);
            Debug.Assert(tenor06M != tenor12M);
            Debug.Assert(tenor01Y == tenor12M);
            Debug.Assert(tenor02Y != tenor12M);

            Debug.Assert(tenor01Y.ToString() == tenor12M.ToString());
            Debug.Assert(tenor01Y.GetHashCode() == tenor12M.GetHashCode());

            var tenor91D = new ql.Period("91D");
            Func<bool> compare91Dversus03MthrowsApplicationException = () =>
            {
                bool hasThrown = false;
                try
                {
                    tenor91D.CompareTo(tenor03M);
                }
                catch (System.ApplicationException)
                {
                    hasThrown = true;
                }
                return hasThrown;
            };

            Debug.Assert(compare91Dversus03MthrowsApplicationException());
        }
    }
}

