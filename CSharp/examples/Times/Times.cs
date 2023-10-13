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

#pragma warning disable CS1718 // Comparison made to same variable

using System;
using System.Collections.Generic;
using System.Diagnostics;
using QL = QuantLib;

namespace TimesTest
{
    class Times
    {
        class TestCaseException : Exception { }

        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main(string[] args)
        {
            try
            {
                DateTime startTime = DateTime.Now;

                RunTestCases();

                DateTime endTime = DateTime.Now;
                TimeSpan delta = endTime - startTime;
                Console.WriteLine();
                Console.WriteLine("Run completed in {0} s", delta.TotalSeconds);
                Console.WriteLine();
            }
            catch (TestCaseException exc)
            {
                Console.Error.WriteLine(exc);
                throw;
            }
        }

        private static void RunTestCases()
        {
            Action<string, IEnumerable<QL.Period>> writePeriods = (heading, periods2Write) =>
            {
                Console.Write($"  {heading}:  ");
                foreach (var item in periods2Write)
                {
                    var itemAsString = item != null ? item.ToString() : "null";
                    Console.Write($"{itemAsString}  ");
                }
                Console.WriteLine();
            };

            Action<bool> testCase = (testResult) =>
            {
                if (!testResult) throw new TestCaseException();
            };

            var tenorNull = null as QL.Period;
            var tenor91D = new QL.Period("91D");
            var tenor03M = new QL.Period("03M");
            var tenor06M = new QL.Period("06M");
            var tenor12M = new QL.Period("12M");
            var tenor01Y = new QL.Period("01Y");
            var tenor02Y = new QL.Period("02Y");

            var periods = new List<QL.Period>() { tenor01Y, tenorNull, tenor02Y, tenor06M, tenor03M };

            Console.WriteLine("Testing sorting of a list.");
            writePeriods("Before sorting", periods);

            periods.Sort();

            writePeriods(" After sorting", periods);

            testCase(periods[0] == tenorNull);
            testCase(periods[1] == tenor03M);
            testCase(periods[2] == tenor06M);
            testCase(periods[3] == tenor01Y);
            testCase(periods[4] == tenor02Y);


            #region test Period.CompareTo(Period)

            Console.WriteLine("test Period.CompareTo(Period)");

            testCase(tenor12M.CompareTo(tenorNull) > 0);
            testCase(tenor12M.CompareTo(tenor03M) > 0);
            testCase(tenor12M.CompareTo(tenor06M) > 0);
            testCase(tenor12M.CompareTo(tenor01Y) == 0);
            testCase(tenor01Y.CompareTo(tenor01Y) == 0);
            testCase(tenor12M.CompareTo(tenor02Y) < 0);

            #endregion

            #region test Period == Period

            Console.WriteLine("test Period == Period");

            testCase(tenorNull == null);
            testCase(null == tenorNull);
            testCase(!(tenorNull == tenor12M));
            testCase(!(tenor12M == null));
            testCase(tenor12M == tenor12M);
            testCase(tenor12M == tenor01Y);
            testCase(!(tenor12M == tenor06M));
            testCase(!(tenor06M == tenor12M));

            #endregion

            #region test Period != Period

            Console.WriteLine("test Period != Period");

            testCase(!(tenorNull != null));
            testCase(!(null != tenorNull));
            testCase(tenorNull != tenor12M);
            testCase(tenor12M != null);
            testCase(!(tenor12M != tenor12M));
            testCase(!(tenor12M != tenor01Y));
            testCase(tenor12M != tenor06M);

            #endregion

            #region test Period < Period

            Console.WriteLine("test Period < Period");

            testCase(!(tenorNull < null));
            testCase(!(null < tenorNull));
            testCase(tenorNull < tenor12M);
            testCase(!(tenor12M < null));
            testCase(!(tenor12M < tenor12M));
            testCase(!(tenor12M < tenor01Y));
            testCase(!(tenor12M < tenor06M));
            testCase(tenor06M < tenor12M);

            #endregion

            #region test Period <= Period

            Console.WriteLine("test Period <= Period");

            testCase(tenorNull <= null);
            testCase(null <= tenorNull);
            testCase(tenorNull <= tenor12M);
            testCase(!(tenor12M <= null));
            testCase(tenor12M <= tenor12M);
            testCase(tenor12M <= tenor01Y);
            testCase(!(tenor12M <= tenor06M));
            testCase(tenor06M <= tenor12M);

            #endregion

            #region test Period > Period

            Console.WriteLine("test Period > Period");

            testCase(!(tenorNull > null));
            testCase(!(null > tenorNull));
            testCase(!(tenorNull > tenor12M));
            testCase(tenor12M > null);
            testCase(!(tenor12M > tenor12M));
            testCase(!(tenor12M > tenor01Y));
            testCase(tenor12M > tenor06M);
            testCase(!(tenor06M > tenor12M));

            #endregion

            #region test Period >= Period

            Console.WriteLine("test Period >= Period");

            testCase(tenorNull >= null);
            testCase(null >= tenorNull);
            testCase(!(tenorNull >= tenor12M));
            testCase(tenor12M >= null);
            testCase(tenor12M >= tenor12M);
            testCase(tenor12M >= tenor01Y);
            testCase(tenor12M >= tenor06M);
            testCase(!(tenor06M >= tenor12M));

            #endregion

            Console.WriteLine("test Period.GetHashCode()");
            testCase(tenor01Y.GetHashCode() == tenor12M.GetHashCode());

            Console.WriteLine("test that uncomparable periods throw");
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

            testCase(compare91Dversus03MthrowsApplicationException());
        }
    }
}

#pragma warning restore CS1718 // Comparison made to same variable
