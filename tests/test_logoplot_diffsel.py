"""Tests ``dms2_logoplot`` with ``--diffsel`` input.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import random
import numpy
import pandas



class test_logoplot_diffsel(unittest.TestCase):
    """Runs ``dms2_logoplot`` with diffsel as input."""

    def setUp(self):
        """Set up input data."""

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_logoplot_diffsel_files/')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        self.indir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './logoplot_input_files/')

        self.name = 'diffsel'

        # define output files
        self.outfiles = [os.path.join(self.testdir, 
                self.name + suffix) for suffix in
                ['.log', '_diffsel.pdf']]

        for f in self.outfiles:
            if os.path.isfile(f):
                os.remove(f)

    def test_logoplot_diffsel(self):
        """Run ``dms2_logoplot`` with diffsel as input."""

        subprocess.check_call([
                'dms2_logoplot',
                '--name', self.name,
                '--diffsel', os.path.join(self.indir, 'HA_mutdiffsel.csv'),
                '--outdir', self.testdir,
                '--overlay1', os.path.join(self.indir, 'HA_mutdiffsel.csv'),
                        'wildtype', 'wildtype',
                '--overlay2', os.path.join(self.indir, 'diffsel_overlay.csv'),
                        'esc', 'known escape mutation',
                '--restrictdiffsel', 'positive',
                '--sepline', 'no',
                '--underlay', 'yes',
                ])

        for f in self.outfiles:
            self.assertTrue(os.path.isfile(f), "did not create {0}".format(f))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
