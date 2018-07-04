"""Tests `dms_tools2.pacbio.CCS` matching and alignment.

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
from pandas.testing import assert_frame_equal, assert_series_equal

import dms_tools2.pacbio
import dms_tools2.minimap2
from dms_tools2.pacbio import qvalsToAccuracy
from dms_tools2 import NTS


Query = collections.namedtuple('Query', ['name', 'barcoded',
        'barcode', 'aligned', 'cigar', 'seq', 'qvals', 'accuracy', 
        'trimmed', 'n_additional', 'target'])
Query.__doc__ = "Holds queries for simulated alignments."


def randSeq(seqlen):
    """Random nucleotide sequence of length `seqlen`."""
    return ''.join([random.choice(NTS) for _ in range(seqlen)])


class test_pacbio_CCS_align_short_codonDMS(unittest.TestCase):
    """Tests `dms_tools2.pacbio.CCS` and related functions.
    
    Data simulates codon-level DMS of a short target."""

    #: length of target sequence
    TARGET_LEN = 1000

    #: random number seed
    SEED = 1

    #: number of queries to simulate
    NQUERIES = 4000

    #: options to :py:mod:`dms_tools2.minimap2.Mapper`
    MAPPER_OPTIONS = dms_tools2.minimap2.OPTIONS_CODON_DMS

    #: deletion lengths range from 1 to this number.
    MAX_DEL_LEN = 10

    #: insertion lengths range from 1 to this number.
    MAX_INS_LEN = 10

    #: probability sequence has a deletion
    DEL_PROB = 0.3

    #: probability sequence has an insertion
    INS_PROB = 0.3

    #: minimum spacing between indels and other mutations
    INDEL_SPACING = 10

    #: number of mutations ranges from 0 to this number
    NMUTS = 4

    #: each mutation changes this many consecutive nucleotides
    MUTLEN = 3

    #: no mutations within this distance from termini
    MUT_BUFFER = 30

    def setUp(self):
        """Create target and query, initialize `CCS` object"""

        self.testdir = (Path(__file__).absolute().parent
                        .joinpath('test_pacbio_ccs_align_files')
                        .joinpath(self.__class__.__name__)
                        )
        Path.mkdir(self.testdir, parents=True, exist_ok=True)

        # target sequence
        random.seed(self.SEED)
        self.targets = {'target{0}'.format(i + 1):randSeq(self.TARGET_LEN)
                        for i in range(2)}
        self.targetfile = self.testdir.joinpath('target.fasta')
        with open(self.targetfile, 'w') as f:
            f.write('\n'.join('>{0}\n{1}'.format(*tup) for tup
                    in self.targets.items()))

        # flanking sequences and barcodes
        self.flank5 = randSeq(20)
        self.flank3 = randSeq(18)
        self.bclen = 12

        # create queries
        self.queries = []
        for iquery in range(self.NQUERIES):
            name = 'query{0}'.format(iquery + 1)
            rand = random.random()
            barcode = cigar = target = ''
            n_additional = -1
            trimmed = barcoded = aligned = False

            if rand < 0.1:
                # should fail matching and aligning
                if rand < 0.5:
                    seq = randSeq(random.randint(self.TARGET_LEN // 2,
                                                self.TARGET_LEN * 2))
                else:
                    # reverse complement won't match
                    seq = dms_tools2.utils.reverseComplement(
                            random.choice(list(self.targets.values())))
            elif rand < 0.2:
                # should pass matching, fail aligning
                barcoded = True
                barcode = randSeq(self.bclen)
                seq = (self.flank5 + 
                       randSeq(random.randint(self.TARGET_LEN // 2,
                                              self.TARGET_LEN * 2)) +
                       barcode +
                       self.flank3
                       )

            else:
                # should pass matching and aligning
                barcoded = aligned = True
                barcode = randSeq(self.bclen)

                # get sites eligible for mutating
                mutsites = list(range(self.MUT_BUFFER,
                        self.TARGET_LEN - self.MUT_BUFFER))

                deletions = []
                if random.random() < self.DEL_PROB:
                    del_len = random.randint(1, self.MAX_DEL_LEN)
                    max_i = max(mutsites) - del_len
                    del_start = random.choice([i for i in mutsites
                            if i < max_i])
                    deletions.append((del_start, del_len))
                    mutsites = [i for i in mutsites if
                            (i < del_start - self.INDEL_SPACING) or
                            (i > del_start + del_len + self.INDEL_SPACING)]

                insertions = []
                if random.random() < self.INS_PROB:
                    ins = randSeq(random.randint(1, self.MAX_INS_LEN))
                    ins_start = random.choice(mutsites)
                    insertions.append((ins_start, ins))
                    mutsites = [i for i in mutsites if
                            (i < ins_start - self.INDEL_SPACING) or
                            (i > ins_start + self.INDEL_SPACING)]

                mutations = []
                for imut in range(random.randint(0, self.NMUTS)):
                    i = random.choice(mutsites)
                    for j in range(i, i + self.MUTLEN):
                        if j in mutsites:
                            mutsites.remove(j)
                            mutations.append((j, random.choice(NTS)))
                (target, targetseq) = random.choice(list(self.targets.items()))
                n_additional = 0
                if random.random() < 0.2 and not deletions:
                    targetseq2 = random.choice(list(self.targets.values()))
                    targetseq = targetseq + targetseq2[ : self.TARGET_LEN // 2]
                    n_additional = 1
                    trimmed = True
                (read, cigar) = dms_tools2.minimap2.mutateSeq(
                        targetseq,
                        mutations, insertions, deletions)
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
                          accuracy=qvalsToAccuracy(qvals, encoding='sanger'),
                          trimmed=trimmed,
                          n_additional=n_additional,
                          target=target))

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


    def test_match_and_align(self):
        """Tests match and alignment on `CCS`."""
        # make sure all queries in `CCS` data frame
        self.assertCountEqual(self.ccs.df.name, [q.name for q in self.queries])

        # now match and check that we get the right entries
        match_str = self.flank5 + \
                    '(?P<read>N+)' + \
                    '(?P<barcode>N{{{0}}})'.format(self.bclen) + \
                    self.flank3
        df = dms_tools2.pacbio.matchSeqs(self.ccs.df,
                match_str, 'CCS', 'barcoded')
        self.assertCountEqual(df.query('barcoded').barcode,
                              [q.barcode for q in self.queries if q.barcoded])

        # now align and check that we get the right entries
        mapper = dms_tools2.minimap2.Mapper(str(self.targetfile),
                self.MAPPER_OPTIONS)
        df = dms_tools2.pacbio.alignSeqs(df,
                mapper, 'read', 'aligned',
                paf_file=str(self.testdir.joinpath('alignment.paf')))
        self.assertCountEqual(df.query('aligned').name,
                [q.name for q in self.queries if q.aligned])
        assert_series_equal(df.aligned_target,
                pandas.Series([q.target for q in self.queries]),
                check_names=False)
        self.assertListEqual(
                    (df.aligned_n_trimmed_query_start.replace(-1, 0) +
                    df.aligned_n_trimmed_query_end.replace(-1, 0) +
                    df.aligned_n_trimmed_target_start.replace(-1, 0) +
                    df.aligned_n_trimmed_target_end.replace(-1, 0)
                    ).astype('bool').tolist(),
                [q.trimmed for q in self.queries])
        assert_series_equal(df.aligned_n_additional,
                pandas.Series([q.n_additional for q in self.queries]),
                check_names=False)
        expected_cigars = dict((q.name, q.cigar) for q in self.queries)
        for row in df.query('aligned').itertuples():
            name = getattr(row, 'name')
            cigar = getattr(row, 'aligned_alignment').cigar_str
            if not getattr(row, 'aligned_n_additional'):
                self.assertEqual(expected_cigars[name], cigar)

        # now test the `matchAndAlignCCS` function
        df2 = dms_tools2.pacbio.matchAndAlignCCS(self.ccs, mapper,
                termini5=self.flank5, gene='N+', spacer=None, umi=None,
                barcode='N{{{0}}}'.format(self.bclen), termini3=self.flank3)
        assert_series_equal(df.aligned, df2.gene_aligned, check_names=False)
        assert_series_equal(df.aligned_target, df2.gene_aligned_target,
                check_names=False)
        assert_series_equal(df.aligned_alignment, df2.gene_aligned_alignment,
                check_names=False)

            
class test_pacbio_CCS_align_long_codonDMS(test_pacbio_CCS_align_short_codonDMS):
    """Tests `dms_tools2.pacbio.CCS` and related functions.
    
    Data simulates codon-level DMS of a long target."""

    #: length of target sequence
    TARGET_LEN = 4000


class test_pacbio_CCS_align_short_virus_w_del(
        test_pacbio_CCS_align_short_codonDMS):
    """Tests `dms_tools2.pacbio.CCS` and related functions.
    
    Data simulates short viral sequences with some long deletions."""

    #: options to :py:mod:`dms_tools2.minimap2.Mapper`
    MAPPER_OPTIONS = dms_tools2.minimap2.OPTIONS_VIRUS_W_DEL

    #: deletion lengths range from 1 to this number.
    MAX_DEL_LEN = 600

    #: insertion lengths range from 1 to this number.
    MAX_INS_LEN = 30

    #: minimum spacing between indels and other mutations
    INDEL_SPACING = 50

    #: each mutation changes this many consecutive nucleotides
    MUTLEN = 1

    #: no mutations within this distance from termini
    MUT_BUFFER = 80


class test_pacbio_CCS_align_long_virus_w_del(
        test_pacbio_CCS_align_short_virus_w_del):
    """Tests `dms_tools2.pacbio.CCS` and related functions.

    Data simulates long viral sequences with some long deletions."""

    #: length of target sequence
    TARGET_LEN = 4000

    #: deletion lengths range from 1 to this number.
    MAX_DEL_LEN = 3200



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
