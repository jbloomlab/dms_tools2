"""Runs all tests."""


import os
import sys
import glob
import unittest


def main():
    """Runs the tests."""

    failurestrings = []
    for test in glob.glob('test_*.py'):
        sys.stderr.write('\nRunning tests in {0}...\n'.format(test))
        test = os.path.splitext(test)[0]
        suite = unittest.TestLoader().loadTestsFromName(test)
        result = unittest.TestResult()
        suite.run(result)
        if result.wasSuccessful():
            sys.stderr.write('All tests were successful.\n')
        else:
            sys.stderr.write('Test(s) FAILED!\n')
            for (testcase, failstring) in result.failures + result.errors:
                failurestrings.append(failstring)

    if not failurestrings:
        sys.stderr.write('\nTesting complete. All passed successfully.\n')
    else:
        sys.stderr.write('\nTesting complete. Failed on the following:\n')
        for fstring in failurestrings:
            sys.stderr.write('\n*********\n{0}\n********\n'.format(fstring))


if __name__ == '__main__':
    main()
