"""Tests `dms_tools2.pacbio.CCS` fitlering and alignment.

In the process, also tests `dms_tools2.minimap2`.
"""

import os
import unittest
import collections
import random
import subprocess
import itertools

import numpy
import pandas
from pandas.testing import assert_frame_equal

import dms_tools2.pacbio
import dms_tools2.minimap2
from dms_tools2.pacbio import qvalsToAccuracy
from dms_tools2 import NTS


class test_pacbio_CCS_align(unittest.TestCase):
    """Tests `dms_tools2.pacbio.CCS` and related functions."""


    def setUp(self):
        """Create target and query, initialize `CCS` object"""

        cwd = os.path.abspath(os.path.dirname(__file__))

        self.testdir = os.path.join(cwd, 'test_pacbio_ccs_align_files')
        if not os.path.isdir(self.testdir):
            os.mkdir(self.testdir)

        # specify target sequence and flanking sequences
        random.seed(0)
        self.target = ''.join(random.choice(NTS) for _ in range(1000))
        self.targetfile = os.path.join(self.testdir, 'target.fasta')
        with open(self.targetfile, 'w') as f:
            f.write('>target\n{0}'.format(self.target))
        self.flank5 = ''.join(random.choice(NTS) for _ in range(20))
        self.flank3 = ''.join(random.choice(NTS) for _ in range(18))
        self.bclen = 12

        # create queries
        Query = collections.namedtuple('Query', ['name', 'barcoded',
                'barcode', 'aligned', 'cigar', 'seq', 'qvals', 'accuracy'])
        self.queries = []
        for i in range(1, 5000):
            name = 'query{0}'.format(i)
            random.seed(i)
            rand = random.random()
            if rand < 0.1:
                # should fail filtering and aligning
                barcoded = aligned = False
                barcode = cigar = ''
                seq = ''.join(random.choice(NTS) for _ in 
                        range(random.randint(300, 2 * len(self.target))))
            elif rand < 0.2:
                # should pass filtering, fail aligning
                barcoded = True
                barcode = ''.join(random.choice(NTS) for _ in range(self.bclen))
                aligned = False
                cigar = ''
                seq = (self.flank5 + 
                       ''.join(random.choice(NTS) for _ in
                               range(random.randint(300, 2 * len(self.target)))) +
                       barcode +
                       self.flank3
                       )
            else:
                # should pass filtering and aligning
                barcoded = aligned = True
                barcode = ''.join(random.choice(NTS) for _ in range(self.bclen))
                if random.random() < 0.4:
                    # mutations up to 5 in length, not too close to ends
                    mutlen = random.randint(1, 5)
                    mutloc = random.randint(12, len(self.target) - 20)
                    mutcigar = mut = ''
                    for i in range(mutloc, mutloc + mutlen):
                        mut += random.choice([nt for nt in NTS if
                                             nt != self.target[i]])
                        mutcigar += '*' + self.target[i].lower() + mut[-1].lower()
                    read = self.target[ : mutloc] + \
                           mut + \
                           self.target[mutloc + mutlen : ]
                    cigar = '=' + self.target[ : mutloc] + \
                            mutcigar + \
                            '=' + self.target[mutloc + mutlen : ]
                elif random.random() < 0.8:
                    del_len = random.randint(1, 200)
                    # deletion no closer than 90 nt to termini
                    mutloc = random.randint(90, len(self.target) - del_len - 90)
                    # different nts on sides of del to avoid ambiguous location
                    while self.target[mutloc - 1] == self.target[mutloc + del_len - 1]:
                        del_len += 1
                    read = self.target[ : mutloc] + self.target[mutloc + del_len : ]
                    cigar = '=' + self.target[ : mutloc] + \
                            '-' + self.target[mutloc : mutloc + del_len].lower() + \
                            '=' + self.target[mutloc + del_len : ]
                else:
                    # no mutations
                    read = self.target
                    cigar = '=' + self.target
                seq = (self.flank5 + 
                       read +
                       barcode +
                       self.flank3
                       )
            qvals = '?' * len(seq)
            self.queries.append(
                    Query(name=name, barcoded=barcoded, barcode=barcode,
                          aligned=aligned, cigar=cigar, seq=seq,
                          qvals=qvals,
                          accuracy=qvalsToAccuracy(qvals, encoding='sanger')))

        # create bamfile of queries
        sam_template = '{0[name]}\t4\t*\t0\t255\t*\t*\t0\t0\t{0[seq]}\t' +\
                       '{0[qvals]}\tnp:i:6\trq:f:{0[accuracy]}'
        samfile = os.path.join(self.testdir, 'queries.sam')
        with open(samfile, 'w') as f:
            for q in self.queries:
                f.write(sam_template.format(q._asdict()) + '\n')
        bamfile = os.path.join(self.testdir, 'queries.bam')
        _ = subprocess.check_call(['samtools', 'view', '-b', '-o',
                                   bamfile, samfile])

        # create CCS object for tests
        self.ccs = dms_tools2.pacbio.CCS('test', bamfile, reportfile=None)


    def test_filter_and_align(self):
        """Tests filtering and alignment on `CCS`."""
        # make sure all queries in `CCS` data frame
        self.assertCountEqual(self.ccs.df.name, [q.name for q in self.queries])

        # now filter and check that we get the right entries
        match_str = self.flank5 + \
                    '(?P<read>N+)' + \
                    '(?P<barcode>N{{{0}}})'.format(self.bclen) + \
                    self.flank3
        self.ccs.filterSeqs(match_str, 'barcoded')
        self.assertCountEqual(self.ccs.df.query('barcoded').barcode,
                              [q.barcode for q in self.queries if q.barcoded])

        # now align and check that we get the right entries
        mapper = dms_tools2.minimap2.Mapper(self.targetfile)
        self.ccs.align(mapper, 'read')
        self.assertCountEqual(self.ccs.df.query('aligned').name,
                             [q.name for q in self.queries if q.aligned])
        expected_cigars = dict((q.name, q.cigar) for q in self.queries)
        for row in self.ccs.df.query('aligned').itertuples():
            name = getattr(row, 'name')
            cigar = getattr(row, 'aligned_cigar')
            self.assertEqual(expected_cigars[name], cigar,
                    "\nexpected:\n{0}\nactual:\n{1}\nclip:{2}, {3}\nread:\n{4}"
                    .format(expected_cigars[name], cigar, 
                            getattr(row, 'aligned_clip_start'), 
                            getattr(row, 'aligned_clip_end'),
                            getattr(row, 'read')))

            



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
