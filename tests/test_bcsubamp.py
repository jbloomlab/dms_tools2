"""Tests ``dms2_bcsubamp``.

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


def randSeq(n, num_N=0):
    """Returns a random nucleotide sequence of length `n`.
    
    This sequence contains `num_N` ``N`` nucleotides."""
    s = [random.choice(dms_tools2.NTS) for i in range(n)]
    for i in random.sample(range(n), num_N):
        s[i] = 'N'
    return ''.join(s)


def generateReadPair(refseq, alignspec, bc1, bc2, r1ext=0, r2ext=0,
        minq=15, r1lowq=0, r2lowq=0, r1bclowq=0, r2bclowq=0, 
        r1mut=0, r2mut=0, r1extmut=0, r2extmut=0,
        r1fail=False, r2fail=False, rname='simulatedread'):
    """Generates a pair of FASTQ reads.

    Args:
        `refseq` (str)
            reference sequence
        `alignspec` (str)
            alignment specs, require even length subamplicons
        `bc1`, `bc2` (str)
            R1 and R2 barcodes
        `r1ext`, `r2ext` (int)
            How many nucleotides R1 / R2 extend past center of subamplicon
        `minq` (int)
            Q scores < this are low quality
        `r1lowq`, `r2lowq` (int)
            How many low quality nucleotides in R1 / R2?
        `r1bclowq`, `r2bclowq` (int)
            How many low quality nucleotide in R1 / R2 barcode?
        `r1mut`, `r2mut` (int)
            How many nucleotide mutations in R1 / R2 before amplicon center?
        `r1extmut`, `r2extmut` (int)
            How many nucleotide mutations in R1 / R2 extensions.
        `r1fail`, `r2fail` (bool)
            Does R1 / R2 fail the quality filter?
        `rname` (str)
            Read name (first part, no spaces)

    Returns:
        The tuple `(r1, r2)` which are each a FASTQ block (4 lines)
        for that read.
    """
    (refseqstart, refseqend, r1start, r2start) = map(int, 
            alignspec.split(','))
    refseqlen = refseqend - refseqstart + 1
    assert refseqlen % 3 == 0
    assert (refseqlen // 3) % 2 == 0, "odd number of codons"

    # first create "perfect" R1 and R2, then add mutations
    r1 = list(refseq[refseqstart - 1 : 
            refseqstart - 1 + refseqlen // 2 + r1ext])
    r2 = list(refseq[refseqend - refseqlen // 2 - r2ext : refseqend])
    for i in random.sample(range(len(r1) - r1ext), r1mut):
        r1[i] = random.choice([nt for nt in dms_tools2.NTS if nt != r1[i]])
    for i in random.sample(range(len(r2) - r2ext), r2mut):
        r2[i] = random.choice([nt for nt in dms_tools2.NTS if nt != r2[i]])
    for i in random.sample(range(len(r1) - r1ext, len(r1)), r1extmut):
        r1[i] = random.choice([nt for nt in dms_tools2.NTS if nt != r1[i]])
    for i in random.sample(range(0, r2ext), r2extmut):
        r2[i] = random.choice([nt for nt in dms_tools2.NTS if nt != r2[i]])
    r1 = ''.join(r1)
    r2 = ''.join(r2)

    # create full reads by adding barcodes and flanking regions
    r1 = bc1 + randSeq(r1start - len(bc1) - 1) + r1
    r2 = bc2 + randSeq(r2start - len(bc2) - 1) + \
            dms_tools2.utils.reverseComplement(r2)

    # create Q scores
    assert minq > 2, "2 assumed low quality"
    q1 = [chr(random.randint(minq, 40) + 33) for i in range(len(r1))]
    for i in random.sample(range(len(bc1)), r1bclowq):
        q1[i] = chr(2 + 33)
    for i in random.sample(range(refseqlen // 2), r1lowq):
        q1[i + r1start - 1] = chr(2 + 33)
    q2 = [chr(random.randint(minq, 40) + 33) for i in range(len(r2))]
    for i in random.sample(range(len(bc2)), r2bclowq):
        q2[i] = chr(2 + 33)
    for i in random.sample(range(refseqlen // 2), r2lowq):
        q2[i + r2start - 1] = chr(2 + 33)
    q1 = ''.join(q1)
    q2 = ''.join(q2)

    # create names
    assert not re.search(r'\s', rname)
    name1 = '@{0} 1:{1}:18:ATCACG'.format(rname, {False:'N', True:'Y'}[r1fail])
    name2 = '@{0} 2:{1}:18:ATCACG'.format(rname, {False:'N', True:'Y'}[r2fail])

    return ('{0}\n{1}\n+\n{2}\n'.format(name1, r1, q1),
            '{0}\n{1}\n+\n{2}\n'.format(name2, r2, q2))



class test_bcsubamp(unittest.TestCase):
    """Runs ``dms2_bcsubamp`` on test data
    
    Use settings that tolerate modest numbers of mutations
    and low-quality nucleotides.
    """

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MAXMUTS = 2
    MINFRACCALL = 0.9
    MINCONCUR = 0.75
    NAME = 'test'
    BCLEN2 = None

    def setUp(self):
        """Set up input data."""

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_bcsubamp_files/')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        # define alignment positions, barcode length, and refseq
        random.seed(1)
        self.alignspecs = ['1,30,15,14', '31,66,12,14']
        self.bclen = 8
        if self.BCLEN2 is None:
            self.bclen2 = self.bclen
        else:
            self.bclen2 = self.BCLEN2
        totrefseqlen = max([int(a.split(',')[1]) for a in self.alignspecs])
        refseq = ''.join([random.choice(dms_tools2.NTS) for
                i in range(totrefseqlen)])
        self.refseqfile = '{0}/{1}_refseq.fasta'.format(self.testdir,
                self.NAME)
        with open(self.refseqfile, 'w') as f:
            f.write('>refseq\n{0}'.format(refseq))

        # now generate reads
        reads = []

        # some reads that fail chastity filter
        self.failfilter = 2
        for i in range(self.failfilter):
            r1fail = random.choice([True, False])
            r2fail = not r1fail
            reads.append(generateReadPair(refseq, 
                    random.choice(self.alignspecs),
                    randSeq(self.bclen), randSeq(self.bclen2),
                    r1fail=r1fail, r2fail=r2fail))

        # some reads with low quality barcodes
        self.lowQbc = 2
        reads.append(generateReadPair(refseq,
                random.choice(self.alignspecs),
                randSeq(self.bclen, num_N=1),
                randSeq(self.bclen2)))
        reads.append(generateReadPair(refseq,
                random.choice(self.alignspecs),
                randSeq(self.bclen),
                randSeq(self.bclen2),
                r1bclowq=1))
        if self.bclen2 > 0:
            self.lowQbc += 2
            reads.append(generateReadPair(refseq,
                    random.choice(self.alignspecs),
                    randSeq(self.bclen),
                    randSeq(self.bclen2, num_N=1)))
            reads.append(generateReadPair(refseq,
                    random.choice(self.alignspecs),
                    randSeq(self.bclen),
                    randSeq(self.bclen2),
                    r2bclowq=1))

        self.nbarcodes = 0
        self.nbarcodesaligned = 0
        self.nbarcodesunaligned = 0
        # create some barcodes with 1 to 4 reads each
        self.barcodes_with_nreads = {}
        for i in range(40):
            self.nbarcodes += 1
            nperbc = random.randint(1, 4)
            if nperbc > 1:
                self.nbarcodesaligned += 1
            if nperbc not in self.barcodes_with_nreads:
                self.barcodes_with_nreads[nperbc] = 1
            else:
                self.barcodes_with_nreads[nperbc] += 1
            bc1 = randSeq(self.bclen)
            bc2 = randSeq(self.bclen2)
            alignspec = random.choice(self.alignspecs)
            for i in range(nperbc):
                r1ext = random.randint(0, 2)
                r2ext = random.randint(0, 2)
                reads.append(generateReadPair(refseq,
                        alignspec, bc1, bc2,
                        r1ext=r1ext, r2ext=r2ext))

        # create a barcode with 3 of 4 reads concurring
        self.barcodes_with_nreads[4] += 1
        self.nbarcodes += 1
        if self.MINCONCUR <= 0.75:
            self.nbarcodesaligned += 1
        else:
            self.nbarcodesunaligned += 1
        bc1 = randSeq(self.bclen)
        bc2 = randSeq(self.bclen2)
        alignspec = random.choice(self.alignspecs)
        (refseqstart, refseqend, r1start, r2start) = map(
                int, alignspec.split(','))
        refseqlen = refseqend - refseqstart + 1
        for i in range(3):
            reads.append(generateReadPair(refseq,
                    alignspec, bc1, bc2))
        reads.append(generateReadPair(refseq, alignspec,
                bc1, bc2, r1mut=refseqlen // 2, r2mut=refseqlen // 2))

        # create a barcode with two mutations
        self.barcodes_with_nreads[2] += 1
        self.nbarcodes += 1
        if self.MAXMUTS < 2:
            self.nbarcodesunaligned += 1
        else:
            self.nbarcodesaligned += 1
        bc1 = randSeq(self.bclen)
        bc2 = randSeq(self.bclen2)
        alignspec = random.choice(self.alignspecs)
        readtup = generateReadPair(refseq, alignspec, bc1, bc2,
                r1mut=1, r2mut=1)
        reads += [readtup, readtup]

        # barcode with an uncalled nucleotide due to mutation in one read
        self.barcodes_with_nreads[2] += 1
        self.nbarcodes += 1
        alignspec = random.choice(self.alignspecs)
        if self.MINFRACCALL == 1:
            self.nbarcodesunaligned += 1
        else:
            self.nbarcodesaligned += 1
        bc1 = randSeq(self.bclen)
        bc2 = randSeq(self.bclen2)
        reads += [generateReadPair(refseq, alignspec, bc1, bc2),
                  generateReadPair(refseq, alignspec, bc1, bc2, r1mut=1)]

        # barcode with an uncalled nucleotide due to low quality
        self.barcodes_with_nreads[2] += 1
        self.nbarcodes += 1
        alignspec = random.choice(self.alignspecs)
        if self.MINFRACCALL == 1:
            self.nbarcodesunaligned += 1
        else:
            self.nbarcodesaligned += 1
        bc1 = randSeq(self.bclen)
        bc2 = randSeq(self.bclen2)
        reads += [generateReadPair(refseq, alignspec, bc1, bc2),
                  generateReadPair(refseq, alignspec, bc1, bc2, r1lowq=1)]

        # barcode with an uncalled nucleotide due to mutation in extension

        random.shuffle(reads)
        self.nreads = len(reads)

        # write reads files
        self.r1file = os.path.join(self.testdir, 
                '{0}_reads_R1_*.fastq.gz'.format(self.NAME))
        self.r2file = self.r1file.replace('_R1', '_R2')
        with gzip.open(self.r1file.replace('*', '1'), 'wt') as f1, \
                gzip.open(self.r2file.replace('*', '1'), 'wt') as f2:
            for (r1, r2) in reads[: 10]:
                f1.write(r1)
                f2.write(r2)
        with gzip.open(self.r1file.replace('*', '2'), 'wt') as f1, \
                gzip.open(self.r2file.replace('*', '2'), 'wt') as f2:
            for (r1, r2) in reads[10 : ]:
                f1.write(r1)
                f2.write(r2)

        # defines output files
        files = ['readstats', 'bcstats', 'readsperbc', 'codoncounts']
        self.outfiles = dict([(f, '{0}/{1}_{2}.csv'.format(
                self.testdir, self.NAME, f)) for f in files])
        for f in self.outfiles.values():
            if os.path.isfile(f):
                os.remove(f)


    def test_dms2_bcsubamp(self):
        """Runs ``dms2_bcsubamp`` on test data."""
        cmds = [
                'dms2_bcsubamp',
                '--name', self.NAME,
                '--refseq', self.refseqfile,
                '--alignspecs'] + self.alignspecs + [
                '--outdir', self.testdir,
                '--R1', os.path.basename(self.r1file),
                '--fastqdir', self.testdir,
                '--bclen', str(self.bclen),
                '--maxmuts', str(self.MAXMUTS),
                '--minfraccall', str(self.MINFRACCALL),
                '--minconcur', str(self.MINCONCUR),
                '--bcinfo'
               ]
        if self.BCLEN2 is not None:
            cmds += ['--bclen2', str(self.bclen2)]
        sys.stderr.write('\nRunning the following command:\n{0}\n'.format(
                ' '.join(cmds)))
        subprocess.check_call(cmds)

        for f in self.outfiles.values():
            self.assertTrue(os.path.isfile(f), "Failed to create {0}".format(f))

        # check on read stats
        readstats = pandas.read_csv(self.outfiles['readstats'])
        self.assertEqual(self.nreads, readstats.at[0, 'total'])
        self.assertEqual(self.failfilter, readstats.at[0, 'fail filter'])
        self.assertEqual(self.lowQbc, readstats.at[0, 'low Q barcode'])

        # check on reads per barcode
        readsperbcstats = pandas.read_csv(self.outfiles['readsperbc'],
                index_col=['number of reads'])
        for (nreads, nbcs) in self.barcodes_with_nreads.items():
            self.assertEqual(nbcs,
                    readsperbcstats.at[nreads, 'number of barcodes'])

        # check on barcode stats
        bcstats = pandas.read_csv(self.outfiles['bcstats'])
        self.assertEqual(self.nbarcodes, bcstats.at[0, 'total'])
        self.assertEqual(self.barcodes_with_nreads[1], 
                bcstats.at[0, 'too few reads'])
        self.assertEqual(self.nbarcodesaligned,
                bcstats.at[0, 'aligned'])
        self.assertEqual(self.nbarcodesunaligned,
                bcstats.at[0, 'not alignable'])


class test_bcsubamp_strictconcur(test_bcsubamp):
    """Tests ``dms2_bcsubamp`` with stricter ``--minconcur``."""
    MINCONCUR = 0.9
    NAME = 'test-minconcur'


class test_bcsubamp_bclen2_4(test_bcsubamp):
    """Tests ``dms2_bcsubamp`` with ``--bclen2 4``."""
    NAME = 'test-bclen2-4'
    BCLEN2 = 4


class test_bcsubamp_bclen2_0(test_bcsubamp):
    """Tests ``dms2_bcsubamp`` with ``--bclen2 0``."""
    NAME = 'test-bclen2-0'
    BCLEN2 = 0


class test_bcsubamp_strictmaxmuts(test_bcsubamp):
    """Tests ``dms2_bcsubamp`` with stricter ``--maxmuts``."""
    MAXMUTS = 0
    NAME = 'test-maxmuts'


class test_bcsubamp_strictminfraccall(test_bcsubamp):
    """Tests ``dms2_bcsubamp`` with stricter ``--minfraccall``."""
    MINFRACCALL = 1.0
    NAME = 'test-minfraccall'


class test_bcsubamp_trimreads(unittest.TestCase):
    """Tests trim reads feature of ``dms2_bcsubamp``."""

    def setUp(self):
        """Set up input data."""

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_bcsubamp_files/')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        self.name = 'test-trim'

        # define alignment positions, barcode length, and refseq
        random.seed(1)
        self.alignspecs = ['1,30,15,14', '31,66,12,14']
        self.bclen = 8
        totrefseqlen = max([int(a.split(',')[1]) for a in self.alignspecs])
        refseq = ''.join([random.choice(dms_tools2.NTS) for
                i in range(totrefseqlen)])
        self.refseqfile = '{0}/{1}_refseq.fasta'.format(self.testdir,
                self.name)
        with open(self.refseqfile, 'w') as f:
            f.write('>refseq\n{0}'.format(refseq))

        # generate reads for each alignspec that need trimming
        reads = []
        r1ext = 6
        r2ext = 9
        for alignspec in self.alignspecs:
            bc1 = randSeq(self.bclen)
            bc2 = randSeq(self.bclen)
            r = generateReadPair(refseq, alignspec, bc1, bc2,
                    r1ext=r1ext, r2ext=r2ext, r1extmut=2, r2extmut=3)
            reads += [r, r]

        # write reads files
        self.r1file = os.path.join(self.testdir, 
                '{0}_reads_R1.fastq.gz'.format(self.name))
        self.r2file = self.r1file.replace('_R1', '_R2')
        with gzip.open(self.r1file, 'wt') as f1, \
                gzip.open(self.r2file, 'wt') as f2:
            for (r1, r2) in reads:
                f1.write(r1)
                f2.write(r2)

    def test_dms2_bcsubamp(self):
        """Runs ``dms2_bcsubamp`` on test data +/- trimming."""
       
        fullr1trim = []
        fullr2trim = []
        for alignspec in self.alignspecs:
            (refseqstart, refseqend, r1start, r2start) = map(
                    int, alignspec.split(','))
            alignlen = refseqend - refseqstart + 1
            fullr1trim.append(str(r1start + alignlen // 2 - 1))
            fullr2trim.append(str(r2start + alignlen // 2 - 1))

        for (r1trim, r2trim, aligned, unaligned, desc) in [
                (['300'], ['300'], 0, len(self.alignspecs), 'notrim'),
                (fullr1trim, fullr2trim, len(self.alignspecs), 0, 'trim')
                ]:
            name = '{0}-{1}'.format(self.name, desc)
            cmds = [
                    'dms2_bcsubamp',
                    '--name', name,
                    '--refseq', self.refseqfile,
                    '--alignspecs'] + self.alignspecs + [
                    '--outdir', self.testdir,
                    '--R1', self.r1file,
                    '--maxmuts', '0',
                    '--minfraccall', '1.0',
                    '--bcinfo',
                    '--bcinfo_csv',
                    '--R1trim'] + r1trim + [
                    '--R2trim'] + r2trim
            sys.stderr.write('\nRunning:\n{0}\n'.format(
                    ' '.join(cmds)))
            subprocess.check_call(cmds)

            # check on barcode stats
            bcstatsfile = '{0}/{1}_bcstats.csv'.format(
                    self.testdir, name)
            bcstats = pandas.read_csv(bcstatsfile)
            self.assertEqual(max(aligned, unaligned),
                    bcstats.at[0, 'total'])
            self.assertEqual(aligned, bcstats.at[0, 'aligned'])
            self.assertEqual(unaligned, bcstats.at[0, 'not alignable'])

            # check barcode info is expected csv
            bcinfofile = '{0}/{1}_bcinfo.csv.gz'.format(self.testdir, name)
            self.assertTrue(os.path.isfile(bcinfofile))
            bcinfo = pandas.read_csv(bcinfofile)
            self.assertEqual(bcstats.at[0, 'total'], len(bcinfo))

            # check on counts at each site
            countsfile = '{0}/{1}_codoncounts.csv'.format(
                    self.testdir, name)
            counts = pandas.read_csv(countsfile).set_index('site')
            self.assertTrue((int(bool(aligned)) == 
                    counts[dms_tools2.CODONS].sum(axis=1)).all())


class test_bcsubamp_counts(unittest.TestCase):
    """Tests generation of counts file."""

    def setUp(self):
        """Set up input data."""

        random.seed(1)

        self.testdir = os.path.join(
                os.path.abspath(os.path.dirname(__file__)),
                './test_bcsubamp_files/')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        self.name = 'test-counts'

        # define alignment positions, barcode length, and refseq
        refseq = 'ATGGACTTCGGG'
        self.refseqfile = '{0}/{1}_refseq.fasta'.format(self.testdir,
                self.name)
        with open(self.refseqfile, 'w') as f:
            f.write('>refseq\n{0}'.format(refseq))

        self.alignspec = '1,12,9,9'
        bc1 = randSeq(8)
        bc2 = randSeq(8)
        r1s = [bc1 + 'ATGGACTNCGGG' + bc1,
               bc1 + 'ATGGACTTCGGG' + bc1,
               bc2 + 'GAGGACTTCGGG' + bc2,
               bc2 + 'GAGGACTTCGAG' + bc2,
              ]

        # generate reads for each alignspec that need trimming
        # write reads files
        self.r1file = os.path.join(self.testdir, 
                '{0}_reads_R1.fastq.gz'.format(self.name))
        self.r2file = self.r1file.replace('_R1', '_R2')
        with gzip.open(self.r1file, 'wt') as f1, \
                gzip.open(self.r2file, 'wt') as f2:
            for r1 in r1s:
                q = chr(73) * len(r1)
                f1.write('@name\n{0}\n+\n{1}\n'.format(r1, q))
                f2.write('@name\n{0}\n+\n{1}\n'.format(
                        dms_tools2.utils.reverseComplement(r1), q))

    def test_counts(self):
        """Make sure we generate the expected codon counts."""
        cmds = [
                'dms2_bcsubamp',
                '--name', self.name,
                '--refseq', self.refseqfile,
                '--alignspecs', self.alignspec,
                '--outdir', self.testdir,
                '--R1', self.r1file,
                '--minfraccall', '0.9',
               ] 
        sys.stderr.write('\nRunning the following command:\n{0}\n'.format(
                ' '.join(cmds)))
        subprocess.check_call(cmds)

        counts = pandas.read_csv('{0}/{1}_codoncounts.csv'.format(
                self.testdir, self.name)).set_index('site')
        self.assertEqual(counts.at[1, 'ATG'], 1)
        self.assertEqual(counts.at[1, 'GAG'], 1)
        self.assertEqual(counts.at[2, 'GAC'], 2)
        self.assertEqual(counts.at[3, 'TTC'], 1)
        self.assertEqual(counts.at[4, 'GGG'], 1)
        self.assertEqual(counts[dms_tools2.CODONS].sum(axis=1).sum(), 6)




if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
