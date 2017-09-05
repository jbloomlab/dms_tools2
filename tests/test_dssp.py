"""Tests `dms_tools2.dssp` module."""


import os
import unittest
import numpy
import pandas
import dms_tools2.dssp


class test_processDSSP(unittest.TestCase):
    """Tests `dms_tools2.dssp.processDSSP`."""

    def test_processDSSP(self):
        """Tests `dms_tools2.dssp.processDSSP`."""
        testdir = os.path.join(os.path.dirname(__file__),
                'dssp_input_files')
        dsspfile = os.path.join(testdir, '1RVX_trimer_sequentialnumbering.dssp')
        df = dms_tools2.dssp.processDSSP(dsspfile, 'A').sort_values('site')
        expected = (pandas.read_csv(os.path.join(testdir, 'expected_output.csv'))
                    .sort_values('site'))
        for c in ['site', 'ASA', 'RSA']:
            self.assertTrue(numpy.allclose(expected[c], df[c]), 
                    'mismatch for {0}'.format(c))
        self.assertTrue(all(expected['SS'].values == df['SS'].values))
        self.assertTrue(all(expected['SS_class'].values == 
                df['SS_class'].values))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
