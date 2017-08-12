
QuantLib-SWIG: language bindings for QuantLib
=============================================

QuantLib-SWIG provides the means to use QuantLib from a number of
languages including Python, Ruby, Perl, C# and Java.

The QuantLib project (<http://quantlib.org>) is aimed at providing a
comprehensive software framework for quantitative finance. QuantLib is
a free/open-source library for modeling, trading, and risk management
in real-life.

QuantLib is Non-Copylefted Free Software and OSI Certified Open Source
Software.


Download and usage
------------------

QuantLib-SWIG can be downloaded from <http://quantlib.org/download.shtml>.

On Linux/Unix, you can run:

    ./configure
    make
    make check
    sudo make install

to build, test and install al modules. 

If you want to use your *local* QuantLib version (e.g. if you want to test
your local dev or have no rights `make install` on your system) you need to 
extend your compiler flags to include the relevant path:

    export CPPFLAGS="-I/path/to/your/ql/include/dir -L/path/to/your/ql/lib/dir":$CPPFLAGS

and update your PATH environment to include the quantlib-config like:

    export PATH=/path/to/your/ql/dir:$PATH

If you're only interested in a
specific language, you can tell make to only work in its subdirectory,
as in:

    make -C Python

Alternatively, you can cd to a specific subdirectory and follow the
instructions in its README file. This is also the procedure for
Windows users.


Questions and feedback
----------------------

Bugs can be reported as a GitHub issue at
<https://github.com/lballabio/QuantLib-SWIG/issues>; if you have a
patch available, you can open a pull request instead (see
"Contributing" below).

You can also use the `quantlib-users` and `quantlib-dev` mailing lists
for feedback, questions, etc.  More information and instructions for
subscribing are at <http://quantlib.org/mailinglists.shtml>.


Contributing
------------

The easiest way to contribute is through pull requests on GitHub.  Get
a GitHub account if you don't have it already and clone the repository
at <https://github.com/lballabio/QuantLib-SWIG> with the "Fork" button
in the top right corner of the page. Check out your clone to your
machine, code away, push your changes to your clone and submit a pull
request; instructions are available at
<https://help.github.com/articles/fork-a-repo>.  (In case you need
them, more detailed instructions for creating pull requests are at
<https://help.github.com/articles/using-pull-requests>, and a basic
guide to GitHub is at
<https://guides.github.com/activities/hello-world/>.

It's likely that we won't merge your code right away, and we'll ask
for some changes instead. Don't be discouraged! That's normal; the
library is complex, and thus it might take some time to become
familiar with it and to use it in an idiomatic way.

We're looking forward to your contributions.

