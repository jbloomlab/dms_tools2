"""Tests `dms_tools2.minimap2.Mapper` with multiply mapping queries."""

from pathlib import Path
import unittest
import random

import dms_tools2.minimap2
from dms_tools2 import NTS


def randSeq(seqlen):
    """Random nucleotide sequence of length `seqlen`."""
    return ''.join([random.choice(NTS) for _ in range(seqlen)])


class test_minimap2_Mapper_multiple_alignments(unittest.TestCase):
    """Tests `dms_tools2.minimap2.Mapper` for multiple alignments."""

    def test_multiple_alignments(self):
        """Test on simulated queries with multiple alignments."""

        testdir = (Path(__file__).absolute().parent
                   .joinpath('test_minimap2_multi_alignments_files')
                   )
        Path.mkdir(testdir, parents=True, exist_ok=True)

        # make two targets
        random.seed(1)
        targets = {'target{0}'.format(i + 1):randSeq(500) for i in range(2)}
        targetfile = testdir.joinpath('target.fasta')
        with open(targetfile, 'w') as f:
            f.write('\n'.join('>{0}\n{1}'.format(*tup)
                    for tup in targets.items()))

        # create queries and expected targets
        query_targets = {}
        queries = {}

        queries['query1'] = targets['target1']
        query_targets['query1'] = ['target1']

        queries['query2'] = targets['target2'] + randSeq(500)
        query_targets['query2'] = ['target2']

        queries['query3'] = targets['target1'] + targets['target1'][20 : ]
        query_targets['query3'] = ['target1', 'target1']

        queries['query4'] = targets['target2'] + targets['target1'][100 : ]
        query_targets['query4'] = ['target2', 'target1']

        queries['query5'] = (targets['target1'][ : 200] + randSeq(100)
                + targets['target2'])
        query_targets['query5'] = ['target2', 'target1']

        queryfile = testdir.joinpath('queries.fasta')
        with open(queryfile, 'w') as f:
            f.write('\n'.join('>{0}\n{1}'.format(*tup)
                    for tup in queries.items()))

        # now map queries
        mapper = dms_tools2.minimap2.Mapper(str(targetfile),
                dms_tools2.minimap2.OPTIONS_VIRUS_W_DEL)
        results = mapper.map(str(queryfile))

        # now check mappings
        for (query, targets_expect) in query_targets.items():
            a = results[query]
            targets_actual = [a.target] + [a2.target for a2 in a.additional]
            self.assertEqual(targets_actual, targets_expect)




if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
