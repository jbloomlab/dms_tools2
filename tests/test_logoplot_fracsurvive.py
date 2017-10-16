"""Tests ``dms2_logoplot`` with ``--fracsurvive`` input.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import random
import numpy
import pandas


class test_logoplot_fracsurvive(unittest.TestCase):
    """Runs ``dms2_logoplot`` with fracsurvive as input."""

    OVERLAY = True

    def setUp(self):
        """Set up input data."""

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_logoplot_fracsurvive_files/')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        self.indir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './logoplot_input_files/')

        if self.OVERLAY:
            self.name = 'fracsurvive'
            self.overlay = [
                    '--overlay1',
                    os.path.join(self.indir, 'HA_fracsurvive.csv'),
                    'wildtype', 'wildtype',
                    '--overlay2',
                    os.path.join(self.indir, 'diffsel_overlay.csv'),
                    'esc', 'known escape mutation'
                    ]
        else:
            self.name = 'fracsurvive-no-overlay'
            self.overlay = []

        # define output files
        self.outfiles = [os.path.join(self.testdir, 
                self.name + suffix) for suffix in
                ['.log', '_fracsurvive.pdf']]

        for f in self.outfiles:
            if os.path.isfile(f):
                os.remove(f)


    def test_logoplot_fracsurvive(self):
        """Run ``dms2_logoplot`` with fracsurvive as input."""

        subprocess.check_call([
                'dms2_logoplot',
                '--name', self.name,
                '--fracsurvive', os.path.join(self.indir, 'HA_fracsurvive.csv'),
                '--outdir', self.testdir,
                '--sepline', 'no',
                '--underlay', 'yes',
                '--nperline', '120',
                '--fracsurvivemax', '5',
                '--scalebar', '1', 'fraction surviving = 1',
                ] + self.overlay)

        for f in self.outfiles:
            self.assertTrue(os.path.isfile(f), "did not create {0}".format(f))


class test_logoplot_fracsurvive_no_overlay(test_logoplot_fracsurvive):
    """Runs ``dms2_logoplot`` with fracsurvive as input and no overlay."""
    OVERLAY = False



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
