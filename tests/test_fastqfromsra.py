"""Tests `dms_tools2.sra.fastqFromSRA`.

Requires installation of ``fastq-dump`` and ``ascp``
for tests to run.
"""

import os
import time
import unittest
import pandas
import dms_tools2.sra


class test_fastqFromSRA(unittest.TestCase):
    """Tests `fastqFromSRA` with `aspera=None`.
    
    Before running, ensure that `FASTQ-DUMP` is a valid path to a
    recent version (>= 2.8.2) of the ``fastq-dump`` program."""
    ASPERA = None
    FASTQ_DUMP = 'fastq-dump' # test fails if this executable not installed
    NCPUS = 1

    def test_fastqFromSRA(self):
        """Tests `fastqFromSRA`."""
        fastqdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_fastqfromsra_files/')

        # small file that downloads quickly
        samples = pandas.DataFrame({'name':['temp'], 'run':['SRR5314806']})

        r1 = os.path.join(fastqdir, 'temp_R1.fastq.gz')
        r2 = os.path.join(fastqdir, 'temp_R2.fastq.gz')
        for f in [r1, r2]:
            if os.path.isfile(f):
                os.remove(f)

        starttime = time.time()
        dms_tools2.sra.fastqFromSRA(samples, 
                self.FASTQ_DUMP,
                fastqdir,
                aspera=self.ASPERA,
                ncpus=self.NCPUS,
                )
        tottime = time.time() - starttime
        print("Time with{0} aspera: {1}".format(
                {True:'', False:'out'}[bool(self.ASPERA)], tottime))

        self.assertTrue(os.path.isfile(r1))
        self.assertTrue(os.path.isfile(r2))
        self.assertTrue(os.path.abspath(r1) == os.path.abspath(
                os.path.join(fastqdir, samples['R1'][0])))
        self.assertTrue(os.path.abspath(r2) == os.path.abspath(
                os.path.join(fastqdir, samples['R2'][0])))

        # now make sure that we can run with overwrite False 
        # and not need to do anything
        dms_tools2.sra.fastqFromSRA(samples, 
                'non existent file',
                fastqdir,
                aspera=('non existent', 'non existent'))


class test_fastqFromSRA_aspera(test_fastqFromSRA):
    """Tests `fastqFromSRA` with `aspera` option.

    Before running, ensure that `ASPERA` gives a valid path
    to ``ascp`` and aspera key.
    """
    ASPERA = ('/app/aspera-connect/3.7.5/bin/ascp',
              '/app/aspera-connect/3.7.5/etc/asperaweb_id_dsa.openssh')


class test_fastqFromSRA_multi(test_fastqFromSRA):
    """Like `test_fastqFromSRA` but multiple CPUs."""
    NCPUS = 2


class test_fastqFromSRA_aspera_multi(test_fastqFromSRA_aspera):
    """Like `test_fastqFromSRA_aspera` but multiple CPUs."""
    NCPUS = 2


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
