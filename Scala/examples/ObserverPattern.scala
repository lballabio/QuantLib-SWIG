/*
 Copyright (C) 2018 Klaus Spanderen


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

package examples;

import org.quantlib.{Array => QArray, _}


/**
 * Check the thread- and gc safeness of QL's observer pattern.
 * This program is likely to fail if QL was compiled with 
 * QL_ENABLE_THREAD_SAFE_OBSERVER_PATTERN not being set.
 *
 * You need to run this using the Java Jar file and JNI library
  */

  
object ObserverPattern {
    def main(args: Array[String]) : Unit = {

        try {
            System.loadLibrary("QuantLibJNI");
        } 
        catch {
            case ex: UnsatisfiedLinkError => {
                  println("please check your LD_LIBRARY_PATH variable")
                throw ex
            }
        }

        
        val quote = new SimpleQuote(100d)
        (0 until 100).par.foreach{j =>
            (0 until 10000).par.foreach{i =>
                val underlying = new QuoteHandle(quote)
                quote.setValue(quote.value())
                if (i == 534) System.gc
            }
        }
    }
}

