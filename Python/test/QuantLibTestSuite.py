"""
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2009 Joseph Malicki

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
"""

import os
import sys
import unittest

import QuantLib as ql


def test():
    print('testing QuantLib', ql.__version__)
    sys.argv[1:1] = ['discover', '-s', os.path.dirname(__file__)]
    unittest.main(module=None, verbosity=2)


if __name__ == '__main__':
    test()
