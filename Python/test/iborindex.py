"""
 Copyright (C) 2018 Wojciech Ślusarski

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

from QuantLib import *
import unittest


class IborIndexTest(unittest.TestCase):
    
    def setUp(self):
        self.euribor3m = Euribor3M()
        # values are not real due to copyrights of the fixing
        force_overwrite = True
        self.euribor3m.addFixing(Date(17, 7, 2018), -0.3, force_overwrite)
        self.euribor3m.addFixings([Date(12, 7, 2018), Date(13, 7, 2018)],
                                  [-0.3, - 0.3],
                                  force_overwrite)

    def testAddFixingFail(self):
        """Testing for RuntimeError while trying to overwrite fixing value"""

        with self.assertRaises(RuntimeError):
            # attempt to overwrite value that is already set at different level
            self.euribor3m.addFixing(Date(17, 7, 2018), -0.4)
            self.euribor3m.addFixings([Date(12, 7, 2018), Date(13, 7, 2018)],
                                      [-0.4, - 0.4])
                                      
    def testAddFixing(self):
        """Testing for overwriting fixing value"""

        force_overwrite = True
        try:
            # attempt to overwrite value that is already set at different level
            self.euribor3m.addFixing(Date(17, 7, 2018), -0.4,
                                     force_overwrite)
            self.euribor3m.addFixings([Date(12, 7, 2018), Date(13, 7, 2018)],
                                      [-0.4, - 0.4],
                                      force_overwrite)
        except RuntimeError as err:
            raise AssertionError("Failed to overwrite index fixixng "
                                 + "{}".format(err))
            
        

if __name__ == '__main__':
    print('testing QuantLib ' + QuantLib.__version__)
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(IborIndexTest, 'test'))
    unittest.TextTestRunner(verbosity=2).run(suite)
