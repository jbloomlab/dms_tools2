"""Tests ``dms2_logoplot`` with ``--prefs`` input.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import random
import numpy
import pandas



class test_logoplot_prefs(unittest.TestCase):
    """Runs ``dms2_logoplot`` with prefs as input."""

    def setUp(self):
        """Set up input data."""

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_logoplot_prefs_files/')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        self.indir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './logoplot_input_files/')

        self.name = 'prefs'

        # define output files
        self.outfiles = [os.path.join(self.testdir, 
                self.name + suffix) for suffix in
                ['.log', '_prefs.pdf']]

        for f in self.outfiles:
            if os.path.isfile(f):
                os.remove(f)

    def test_logoplot_prefs(self):
        """Run ``dms2_logoplot`` with prefs as input."""

        subprocess.check_call([
                'dms2_logoplot',
                '--name', self.name,
                '--prefs', os.path.join(self.indir, 'HA_prefs.csv'),
                '--outdir', self.testdir,
                '--overlay1', os.path.join(self.indir, 'continuous_overlay.csv'),
                        'cp', 'continuous-property',
                '--overlay2', os.path.join(self.indir, 'discrete_overlay.csv'),
                        'dp', 'discrete-property',
                '--stringency', '2.0',
                '--underlay', 'yes',
                ])

        for f in self.outfiles:
            self.assertTrue(os.path.isfile(f), "did not create {0}".format(f))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
