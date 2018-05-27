"""Tests `dms_tools2.pacbio.CCS` fitlering and alignment.

In the process, also tests `dms_tools2.minimap2`.
"""

from pathlib import Path
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


class test_pacbio_CCS_align_1000(unittest.TestCase):
    """Tests `dms_tools2.pacbio.CCS` and related functions.
    
    Tests on simulated queries on target of length 1000."""

    TARGET_LEN = 1000

    def setUp(self):
        """Create target and query, initialize `CCS` object"""

        cwd = Path(__file__).absolute().parent

        self.testdir = (cwd.joinpath('test_pacbio_ccs_align_files')
                           .joinpath(str(self.TARGET_LEN))
                           )
        Path.mkdir(self.testdir, parents=True, exist_ok=True)

        # specify target sequence and flanking sequences
        random.seed(0)
        self.target = ''.join(random.choice(NTS)
                for _ in range(self.TARGET_LEN))
        self.targetfile = self.testdir.joinpath('target.fasta')
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
                        range(random.randint(300, 2 * self.TARGET_LEN)))
            elif rand < 0.2:
                # should pass filtering, fail aligning
                barcoded = True
                barcode = ''.join(random.choice(NTS) for _ in range(self.bclen))
                aligned = False
                cigar = ''
                seq = (self.flank5 + 
                       ''.join(random.choice(NTS) for _ in
                               range(random.randint(300, 2 * self.TARGET_LEN))) +
                       barcode +
                       self.flank3
                       )
            else:
                # should pass filtering and aligning
                barcoded = aligned = True
                barcode = ''.join(random.choice(NTS) for _ in range(self.bclen))
                rand = random.random()
                if rand < 0.3:
                    # mutations up to 5 in length, not too close to ends
                    mutlen = random.randint(1, 5)
                    mutloc = random.randint(12, self.TARGET_LEN - 20)
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
                elif rand < 0.6:
                    # random deletion
                    #del_len = random.randint(1, int(0.5 * self.TARGET_LEN))
                    del_len = random.randint(1, 10)
                    # deletion no closer than 25 nt to termini
                    delmargin = 25
                    del_len = min(del_len, self.TARGET_LEN - 2 * delmargin)
                    mutloc = random.randint(delmargin,
                            self.TARGET_LEN - del_len - delmargin)
                    read = self.target[ : mutloc] + self.target[mutloc + del_len : ]
                    cigar = ('=' + self.target[ : mutloc] + 
                             '-' + self.target[mutloc : mutloc + del_len].lower() +
                             '=' + self.target[mutloc + del_len : ])
                elif rand < 0.9:
                    # random insertion
                    ins_len = random.randint(1, 10)
                    # deletion no closer than 25 nt to termini
                    insmargin = 25
                    ins_len = min(ins_len, self.TARGET_LEN - 2 * insmargin)
                    mutloc = random.randint(insmargin,
                            self.TARGET_LEN - ins_len - insmargin)
                    ins = ''.join(random.choice(NTS) for _ in range(ins_len))
                    read = self.target[ : mutloc] + ins + self.target[mutloc : ]
                    cigar = ('=' + self.target[ : mutloc] +
                             '+' + ins.lower() + 
                             '=' + self.target[mutloc : ])
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
                    Query(name=name, 
                          barcoded=barcoded,
                          barcode=barcode,
                          aligned=aligned, 
                          cigar=dms_tools2.minimap2.shiftIndels(cigar),
                          seq=seq,
                          qvals=qvals,
                          accuracy=qvalsToAccuracy(qvals, encoding='sanger')))

        # create fasta file of queries
        with self.testdir.joinpath('queries.fasta').open('w') as f:
            f.write('\n'.join('>{0}\n{1}'.format(q.name, q.seq)
                    for q in self.queries))

        # create bamfile of queries
        sam_template = '{0[name]}\t4\t*\t0\t255\t*\t*\t0\t0\t{0[seq]}\t' +\
                       '{0[qvals]}\tnp:i:6\trq:f:{0[accuracy]}'
        samfile = self.testdir.joinpath('queries.sam')
        with open(samfile, 'w') as f:
            for q in self.queries:
                f.write(sam_template.format(q._asdict()) + '\n')
        bamfile = self.testdir.joinpath('queries.bam')
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
        mapper = dms_tools2.minimap2.Mapper(str(self.targetfile))
        self.ccs.align(mapper, 'read',
                paf_file=str(self.testdir.joinpath('alignment.paf')))
        self.assertEqual(len(self.ccs.df.query('aligned')),
                         len([q.name for q in self.queries if q.aligned]))
        self.assertCountEqual(self.ccs.df.query('aligned').name,
                             [q.name for q in self.queries if q.aligned])
        expected_cigars = dict((q.name, q.cigar) for q in self.queries)
        for row in self.ccs.df.query('aligned').itertuples():
            name = getattr(row, 'name')
            cigar = getattr(row, 'aligned_cigar')
            self.assertEqual(expected_cigars[name], cigar,
                    "\nquery: {0}\n\nexpected:\n{1}\n\nactual:\n{2}\n\n"
                    "clip:{3}, {4}\n\nread:\n{5}"
                    .format(name, expected_cigars[name], cigar, 
                            getattr(row, 'aligned_clip_start'), 
                            getattr(row, 'aligned_clip_end'),
                            getattr(row, 'read')))

            

#class test_pacbio_CCS_align_2500(test_pacbio_CCS_align_1000):
    """Tests `dms_tools2.pacbio.CCS` and related functions.
    
    Tests on simulated queries on target of length 2500."""

#    TARGET_LEN = 2500


#class test_pacbio_CCS_align_5000(test_pacbio_CCS_align_1000):
    """Tests `dms_tools2.pacbio.CCS` and related functions.
    
    Tests on simulated queries on target of length 5000."""

#    TARGET_LEN = 5000


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
