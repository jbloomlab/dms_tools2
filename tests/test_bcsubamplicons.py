"""Tests ``dms2_bcsubamplicons``.

Written by Jesse Bloom."""


import sys
import re
import os
import unittest
import subprocess
import random
import gzip
import pandas
import dms_tools2
import dms_tools2.utils


def randSeq(n):
    """Returns a random nucleotide sequence of length `n`."""
    return ''.join([random.choice(dms_tools2.NTS) for i in range(n)])


def generateReadPair(refseq, alignspec, bc1, bc2, r1ext=0, r2ext=0,
        minq=15, r1lowq=0, r2lowq=0, r1extlowq=0, r2extlowq=0, 
        r1mut=0, r2mut=0, r1extmut=0, r2extmut=0, r1fail=False, 
        r2fail=False, rname='simulatedread'):
    """Generates a pair of FASTQ reads.

    Args:
        `refseq` (str)
            reference sequence
        `alignspec` (str)
            alignment specs, require even length subamplicons
        `bc1` (str)
            R1 barcode
        `bc2` (str)
            R2 barcode
        `r1ext` (int)
            How many codons does R1 extend past center of subamplicon?
        `r2ext` (int)
            How many codons does R2 extend past center of subamplicon?
        `minq` (int)
            Q scores < this are low quality
        `r1lowq` (int)
            How many low quality codons in R1 prior to amplicon center?
        `r2lowq` (int)
            How many low quality codons in R2 prior to amplicon center?
        `r1extlowq` (int)
            How many low quality codons in R1 after amplicon center?
        `r2extlowq` (int)
            How many low quality codons in R2 after amplicon center?
        `r1mut` (int)
            How many codon mutations in R1 before amplicon center?
        `r2mut` (int)
            How many codon mutations in R2 before amplicon center?
        `r1extmut` (int)
            How many codon mutations in R1 after amplicon center?
        `r2extmut` (int)
            How many codon mutations in R2 after amplicon center?
        `r1fail` (bool)
            Does R1 fail the quality filter?
        `r2fail` (bool) 
            Does R2 fail the quality filter?
        `rname` (str)
            Read name (first part, no spaces)

    Returns:
        The tuple `(r1, r2)` which are each a FASTQ block (4 lines)
        for that read.
    """
    (refseqstart, refseqend, r1start, r2start) = map(int, alignspec.split(','))
    refseqlen = refseqend - refseqstart + 1
    assert refseqlen % 3 == 0
    assert (refseqlen // 3) % 2 == 0, "not even number of codons in subamplicon"

    # first create "perfect" R1 and R2
    r1 = refseq[refseqstart - 1 : refseqlen // 2 + 3 * r1ext]
    r2 = refseq[refseqlen // 2 - 3 * r2ext : refseqend]

    # create full reads by adding barcodes and flanking regions
    r1 = bc1 + randSeq(r1start - len(bc1) - 1) + r1
    r2 = bc2 + randSeq(r2start - len(bc2) - 1) + \
            dms_tools2.utils.reverseComplement(r2)

    # create Q scores
    q1 = ''.join([chr(random.randint(minq, 40) + 33) for i in range(len(r1))])
    q2 = ''.join([chr(random.randint(minq, 40) + 33) for i in range(len(r2))])

    # create names
    assert not re.search('\s', rname)
    name1 = '@{0} 1:{1}:18:ATCACG'.format(rname, {False:'N', True:'Y'}[r1fail])
    name2 = '@{0} 2:{1}:18:ATCACG'.format(rname, {False:'N', True:'Y'}[r2fail])

    return ('{0}\n{1}\n+\n{2}\n'.format(name1, r1, q1),
            '{0}\n{1}\n+\n{2}\n'.format(name2, r2, q2))


class test_bcsubamplicons(unittest.TestCase):
    """Runs ``dms2_bcsubamplicons`` on test data.
    """

    def setUp(self):
        """Set up input data."""

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_bcsubamplicons_files/')

        # define alignment positions, barcode length, and refseq
        random.seed(1)
        self.alignspecs = ['1,30,15,14', '31,66,12,14']
        self.bclen = 8
        refseqlen = max([int(a.split(',')[1]) for a in self.alignspecs])
        refseq = ''.join([random.choice(dms_tools2.CODONS) for
                i in range(refseqlen)])
        self.refseqfile = os.path.join(self.testdir, 'refseq.fasta')
        with open(self.refseqfile, 'w') as f:
            f.write('>refseq\n{0}'.format(refseq))

        # now generate reads
        reads = []

        self.failfilter = 2
        for i in range(self.failfilter):
            r1fail = random.choice([True, False])
            r2fail = not r1fail
            reads.append(generateReadPair(refseq, 
                    random.choice(self.alignspecs),
                    randSeq(self.bclen), randSeq(self.bclen),
                    r1fail=r1fail, r2fail=r2fail))

        self.nreads = len(reads)

        # write reads files
        self.r1file = os.path.join(self.testdir, 'read_R1.fastq.gz')
        self.r2file = self.r1file.replace('_R1', '_R2')
        with gzip.open(self.r1file, 'wt') as f1, \
                gzip.open(self.r2file, 'wt') as f2:
            for (r1, r2) in reads:
                f1.write(r1)
                f2.write(r2)

        # defines output files
        self.name = 'test'
        files = ['readstats']
        self.outfiles = dict([(f, '{0}/{1}_{2}.csv'.format(
                self.testdir, self.name, f)) for f in files])

        for f in self.outfiles.values():
            if os.path.isfile(f):
                os.remove(f)


    def test_dms2_bcsubamplicons(self):
        """Runs ``dms2_bcsubamplicons`` on test data."""
        cmds = [
                'dms2_bcsubamplicons',
                '--name', self.name,
                '--refseq', self.refseqfile,
                '--alignspecs'] + self.alignspecs + [
                '--outdir', self.testdir,
                '--R1', self.r1file
               ] 
        sys.stderr.write('\nRunning the following command:\n{0}\n'.format(
                ' '.join(cmds)))
        subprocess.check_call(cmds)

        for f in self.outfiles.values():
            self.assertTrue(os.path.isfile(f), "Failed to create {0}".format(f))

        # check on read stats
        readstats = pandas.read_csv(self.outfiles['readstats'], 
                index_col='category')
        self.assertEqual(self.nreads, 
                readstats.at['total', 'number of reads'])
        self.assertEqual(self.failfilter, 
                readstats.at['fail filter', 'number of reads'])



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
