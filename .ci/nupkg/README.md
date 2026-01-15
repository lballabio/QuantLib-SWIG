# QuantLib: the free/open-source library for quantitative finance

The QuantLib project (<https://www.quantlib.org>) is aimed at providing a comprehensive software framework for quantitative finance. QuantLib is a free/open-source library for modeling, trading, and risk management in real-life.

QuantLib is Non-Copylefted Free Software and OSI Certified Open Source Software.


## Note on multi-threading

The underlying C++ library is not thread-safe.  However, it provides some configuration options that make it possible to use these wrappers from C# using some limited multi-threading.  In particular:

- the globals in the library are per-thread; most notably, the evaluation date and the stored index fixings;
- there is code in place to prevent the garbage collector to interfere with notifications inside the library.


Any more complex attempt at multi-threading will probably fail.  We suggest to avoid sharing objects and state across threads; each thread should have its set of objects, evaluation date etc.
