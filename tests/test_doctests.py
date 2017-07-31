"""Runs doctests.

Written by Jesse Bloom."""


import os
import sys
import glob
import doctest
import unittest
import doctest
import pkgutil
import dms_tools2


def main():
    """Run the doctests."""

    sys.stderr.write('Running doctests...\n')
    failurestrings = []

    # test all modules with doctest
    for (importer, modname, ispkg) in pkgutil.iter_modules(dms_tools2.__path__):
        if (not ispkg) and modname[0] != '_':
            sys.stderr.write('\nTesting {0} with doctest... '.format(modname))
            module = __import__('dms_tools2.{0}'.format(modname), 
                    None, None, modname.split('.'))
            suite = doctest.DocTestSuite(module)
            del module
            result = unittest.TestResult()
            suite.run(result)
            if result.wasSuccessful():
                sys.stderr.write('all {0} tests were successful.\n'.format(
                        result.testsRun))
            else:
                sys.stderr.write('test FAILED!\n')
                for (testcase, failstring) in result.failures:
                    failurestrings.append(failstring)

    # print summary of failures
    if not failurestrings:
        sys.stderr.write('\nTesting complete. All tests passed.\n')
    else:
        sys.stderr.write('\nTesting complete. Failed the following tests:\n')
        for failstring in failurestrings:
            sys.stderr.write('\n*******\n{0}\n\n********\n'.format(failstring))


if __name__ == '__main__':
    main() # run the program
