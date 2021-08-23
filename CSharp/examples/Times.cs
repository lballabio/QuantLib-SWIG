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
            var tenor01Y = new ql.Period("01Y");
            var tenor12M = new ql.Period("12M");

            Console.WriteLine($"{tenor01Y}, {tenor01Y.GetHashCode()}");
            Console.WriteLine($"{tenor12M}, {tenor12M.GetHashCode()}");

            var tenor03M = new ql.Period("03M");
            var tenor90D = new ql.Period("91D");

            Console.WriteLine(tenor03M.ToString());
            Console.WriteLine(tenor03M.__repr__());
            Console.WriteLine(tenor03M.__str__());

            tenor90D.CompareTo(tenor03M);
        }
    }
}

