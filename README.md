
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

If you want to use your *local* QuantLib dev version without installing it (i.e.
**no** `make install` after your `make` - e.g. if you want to test your 
local dev or have no rights to `make install` on your system) you can do the 
following (this needs to be done before `configure`):

- configure your QuantLib build to use your dev include and lib path like
`$./configure --prefix=/path/to/your/qldir/ 
--libdir=/path/to/your/qldir/ql/.libs 
--includedir=/path/to/your/qldir/`
since you are not installing, this will only provide the necessary path
information for the quantlib-config file.
- add *quantlib-config* to your path enviroment by extending it like
`$ export PATH:/path/to/your/qldir:$PATH`
- add search path to your swigged QuantLib shared object lib (.so extension),
so that it can find your dev libQuantLib.so.
`$ export CPPFLAGS=-Wl,-rpath,/path/to/your/qldir/ql/.libs`
This flag will include the absolute path to your dev lib.

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

