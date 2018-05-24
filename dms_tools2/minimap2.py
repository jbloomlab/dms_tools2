"""
===============
minimap2
===============

Runs `minimap2 <https://lh3.github.io/minimap2/>`_
aligner. Although this aligner has it's own Python AP
Python API (`mappy <https://github.com/lh3/minimap2/tree/master/python>`_),
it does not support all options.
Therefore, this module interacts directly with
the ``minimap2`` executable on the command line.
The module is tested with ``minimap2`` version 2.10.
"""


import os
import re
import io
import functools
import subprocess
import tempfile
import collections
import random

import Bio.SeqIO

from dms_tools2 import NTS


#: ``minimap2`` settings that work well for a query
#: being matched to a short target where the target
#: is known to be in the same orientation as the query,
#: and some of the mutations are codon mutations (this
#: last point requires care on opening small gaps).
ORIENTED_READ = ['-uf', 
                 '-k15',
                 '-w5',
                 '-g2000',
                 '-G200k',
                 '-A1',
                 '-B2',
                 '-O6,40',
                 '-E1,0',
                 '-z200',
                 '--secondary=no']

# namedtuple to hold alignments
Alignment = collections.namedtuple('Alignment',
        ['target', 'r_st', 'r_en', 'q_len', 'q_st',
         'q_en', 'strand', 'mapq', 'cigar_str'])
Alignment.__doc__ = "Alignment of a query to a target (reference)."
Alignment.target.__doc__ = "Target (reference) to which query was aligned."
Alignment.r_st.__doc__ = "Alignment start in target (0 based)."
Alignment.r_en.__doc__ = "Alignment end in target (0 based)."
Alignment.q_st.__doc__ = "Alignment start in query (0 based)."
Alignment.q_en.__doc__ = "Alignment end in query (0 based)."
Alignment.q_len.__doc__ = "Total length of query prior to any clipping."
Alignment.strand.__doc__ = "1 if aligns in forward polarity, -1 if in reverse."
Alignment.mapq.__doc__ = "Mapping quality."
Alignment.cigar_str.__doc__ = "CIGAR in `PAF long format <https://github.com/lh3/minimap2>`_"


class Mapper:
    """Class to run ``minimap2`` at command line and get results.

    Args:
        `target` (str)
            FASTA file with target (reference) to which we align
            reads.
        `prog` (str)
            Path to ``minimap2`` executable.
        `options` (list)
            Command line options to ``minimap2``.

    Attributes:
        `target` (str)
            Target (reference) set at initialization.
        `prog` (str)
            Path to ``minimap2`` set at initialization.
        `options` (list)
            Options to ``minimap2`` set at initialization.
        `version` (str)
            Version of ``minimap2``.

    Here is an example where we align a few reads to two target
    sequences.

    First, we generate a few target sequences:

    >>> targetlen = 200
    >>> random.seed(1)
    >>> targets = {}
    >>> for i in [1, 2]:
    ...     targets['target{0}'.format(i)] = ''.join(random.choice(NTS)
    ...             for _ in range(targetlen))

    Now we generate some queries. One is a random sequence that should not
    align, and the other two are substrings of the targets into which we
    have introduced a single mutation or indel. The names of the queries
    give their target, start in query, end in query, cigar string:

    >>> queries = {'randseq':''.join(random.choice(NTS) for _ in range(180))}
    >>> for qstart, qend, mut in [(0, 183, 'mut53'), (36, 194, 'del140')]:
    ...     target = random.choice(list(targets.keys()))
    ...     qseq = targets[target][qstart : qend]
    ...     mutsite = int(mut[3 : ])
    ...     if 'mut' in mut:
    ...         wt = qseq[mutsite]
    ...         mut = random.choice([nt for nt in NTS if nt != wt])
    ...         cigar = ('=' + qseq[ : mutsite] + '*' + wt.lower() +
    ...                  mut.lower() + '=' + qseq[mutsite + 1 : ])
    ...         qseq = qseq[ : mutsite] + mut + qseq[mutsite + 1 : ]
    ...     elif 'del' in mut:
    ...         cigar = ('=' + qseq[ : mutsite] + '-' +
    ...                  qseq[mutsite].lower() + '=' + qseq[mutsite + 1 : ])
    ...         qseq = qseq[ : mutsite] + qseq[mutsite + 1 : ]
    ...     queryname = '_'.join(map(str, [target, qstart, qend, cigar]))
    ...     queries[queryname] = qseq


    Now map the queries to the targets:

    >>> TempFile = functools.partial(tempfile.NamedTemporaryFile, mode='w')
    >>> with TempFile() as targetfile, TempFile() as queryfile:
    ...     _ = targetfile.write('\\n'.join('>{0}\\n{1}'.format(*tup)
    ...                          for tup in targets.items()))
    ...     targetfile.flush()
    ...     _ = queryfile.write('\\n'.join('>{0}\\n{1}'.format(*tup)
    ...                         for tup in queries.items()))
    ...     queryfile.flush()
    ...     mapper = Mapper(targetfile.name)
    ...     alignments = mapper.map(queryfile.name)

    Now make sure we find the expected alignments:

    >>> set(alignments.keys()) == set(q for q in queries if q != 'randseq')
    True
    >>> matched = []
    >>> for (query, a) in alignments.items():
    ...     expected = query.split('_')
    ...     matched.append(a.target == expected[0])
    ...     matched.append([a.r_st, a.r_en] == list(map(int, expected[1 : 3])))
    ...     matched.append([a.q_st, a.q_en] == [0, len(queries[query])])
    ...     matched.append(a.cigar_str == expected[3])
    ...     matched.append(a.strand == 1)
    >>> all(matched)
    True
    """

    def __init__(self, target, prog='minimap2', options=ORIENTED_READ):
        """See main :class:`Mapper` doc string."""
        try:
            version = subprocess.check_output([prog, '--version'])
        except:
            raise ValueError("Can't execute `prog` {0}".format(prog))
        self.version = version.strip().decode('utf-8')
        self.prog = prog
        self.options = options
        assert os.path.isfile(target), "no `target` {0}".format(target)
        self.target = target


    def map(self, queryfile, return_dict=True, outfile=None,
            unclip_matches=True):
        """Map query sequences to reference target.

        Aligns query sequences to `target`. Adds ``--c --cs=long``
        arguments to `options` to get a long CIGAR string, and
        returns the results as a dictionary and/or writes them
        to a PAF file.

        Will raise an error if there are multiple mappings for
        a query.

        Args:
            `queryfile` (str)
                FASTA file with query sequences to align.
                Headers should be unique.
            `return_dict` (bool)
                If `True`, return a dictionary keyed by
                each header in `query` and with values
                being the alignments. If you have a very
                large `query`, you might set this to `False`
                to avoid reading the alignments into memory.
            `outfile` (`None` or str)
                Name of output file containing alignment
                results in PAF format if a str is provided.
                Provide `None` if you just want to access
                results via `return_dict` and don't want to
                create a permanent output file.
            `unclip_matches` (bool)
                Pass alignments through :meth:`unclipMatches`.
                Only applies to alignments returned in dict,
                does not affect any results in `outfile`.

        Returns:
            Either a dict (if `return_dict` is `True`)
            or `None` (if `return_dict` is `False`). If returning
            a dictionary, the keys are the name each sequence
            in `queryfile`, and the values are alignments as 
            :class:`Alignment` objects.
        """
        assert os.path.isfile(queryfile), "no `queryfile` {0}".format(queryfile)

        assert '-a' not in self.options, \
                "output should be PAF format, not SAM"
        for arg in ['-c', '--cs=long']:
            if arg not in self.options:
                self.options.append(arg)

        if outfile is None:
            fout = tempfile.TemporaryFile('w+')
        else:
            fout = open(outfile, 'w+')

        stderr = tempfile.TemporaryFile()
        try:
            _ = subprocess.check_call(
                    [self.prog] + self.options + [self.target, queryfile],
                    stdout=fout, stderr=stderr)
            if return_dict:
                fout.seek(0)
                d = {}
                for query, alignment in parsePAF(fout):
                    assert query not in d, "duplicate query {0}".format(query)
                    d[query] = alignment
            else:
                d = None
        finally:
            fout.close()
            stderr.close()

        if unclip_matches:
            targetseqs = {seq.name:str(seq.seq) for seq in 
                          Bio.SeqIO.parse(self.target, 'fasta')}
            queryseqs = {seq.name:str(seq.seq) for seq in
                         Bio.SeqIO.parse(queryfile, 'fasta')}
            for query in list(d.keys()):
                a = d[query]
                assert a.q_len == len(queryseqs[query])
                a = unclipMatches(a, targetseqs[a.target], queryseqs[query])
                d[query] = a

        return d


def unclipMatches(a, target, query):
    """Undo any soft clipping of exact matches.

    For some unknown reason, ``minimap2`` sometimes soft clips
    queries even when the clipped regions exactly match the 
    target. This function takes such alignments and undoes 
    the clipping to the extent possible without introducing
    any mutations / gaps or otherwise affecting the alignment.

    Args:
        `a` (:class:`Alignment`)
            Alignment of `query` to `target` to un-clip
        `target` (str)
            Target sequence to which query is aligned.
        `query` (str)
            Query sequence.

    Returns:
        A new :class:`Alignment` with soft clipping of exact
        matches undone.

    Here are examples:

    Do nothing if no soft clipping.

    >>> target = 'ATGCAATGA'
    >>> query = 'TACAAT'
    >>> a = Alignment(mapq=60, strand=1, r_st=1, r_en=7,
    ...         target='target', q_st=0, q_en=6, q_len=6,
    ...         cigar_str='=T*gaCAAT')
    >>> unclipMatches(a, target, query) == a
    True

    Now some soft clipping that is undone at start:

    >>> query = 'AAGCAATA'
    >>> a = Alignment(mapq=60, strand=1, r_st=1, r_en=7,
    ...         target='target', q_st=1, q_en=7, q_len=8,
    ...         cigar_str='*taGCAAT')
    >>> a2 = unclipMatches(a, target, query)
    >>> a2.q_st
    0
    >>> a2.r_st
    0
    >>> a2.q_en == a.q_en
    True
    >>> a2.r_en == a.r_en
    True
    >>> a2.cigar_str
    '=A*taGCAAT'

    Some soft clipping that is undone at both ends:

    >>> query = 'ACTGCAATGA'
    >>> target = 'ATGCAATGA'
    >>> a = Alignment(mapq=60, strand=1, r_st=3, r_en=7,
    ...         target='target', q_st=4, q_en=8, q_len=10,
    ...         cigar_str='=CAAT')
    >>> a3 = unclipMatches(a, target, query)
    >>> a3.q_st
    2
    >>> a3.q_en
    10
    >>> a3.r_st
    1
    >>> a3.r_en
    9
    >>> a3.cigar_str
    '=TGCAATGA'
    """
    assert a.strand == 1, "not implemented for - strand"
    assert a.q_len == len(query), "wrong length query"

    # unclip from start
    unclip = 0
    while ((unclip < a.q_st) and (unclip < a.r_st) and
           (query[a.q_st - unclip - 1] == target[a.r_st - unclip - 1])):
        unclip += 1
    if unclip > 0:
        if a.cigar_str[0] == '=':
            trimfirst = 1
        else:
            trimfirst = 0
        new_cigar = '=' + query[a.q_st - unclip : a.q_st] + \
                        a.cigar_str[trimfirst : ]
        a = Alignment(
                mapq=a.mapq,
                strand=a.strand,
                r_st=a.r_st - unclip,
                r_en=a.r_en,
                target=a.target,
                q_st=a.q_st - unclip,
                q_en=a.q_en,
                q_len=a.q_len,
                cigar_str=new_cigar,
                )

    # unclip from end
    unclip = 0
    while ((unclip + a.q_en < a.q_len) and
           (unclip + a.r_en < len(target)) and
           (query[a.q_en + unclip] == target[a.r_en + unclip])):
        unclip += 1
    if unclip > 0:
        if a.cigar_str[-1].isupper():
            addchar = ''
        else:
            addchar = '='
        new_cigar = a.cigar_str + addchar + query[a.q_en : a.q_en + unclip]
        a = Alignment(
                mapq=a.mapq,
                strand=a.strand,
                r_st=a.r_st,
                r_en=a.r_en + unclip,
                target=a.target,
                q_st=a.q_st,
                q_en=a.q_en + unclip,
                q_len=a.q_len,
                cigar_str=new_cigar,
                )

    return a



def parsePAF(paf_file):
    """Parse ``*.paf`` file as created by ``minimap2``.

    The ``*.paf`` file is assumed to be created with
    the ``minimap2`` options ``-c --cs=long``, which
    creates long `cs` tags with the CIGAR string.

    PAF format is `described here <https://github.com/lh3/miniasm/blob/master/PAF.md>`_.
    PAF long CIGAR string format from ``minimap2`` is 
    `detailed here <https://github.com/lh3/minimap2>`_.

    Args:
        `paf_file` (str or iterator)
            If str, should be name of ``*.paf`` file.
            Otherwise should be an iterator that returns
            lines as would be read from a PAF file.

    Returns:
        A generator that yields query / alignments on each line in 
        `paf_file`. Returned as the 2-tuple `(query_name, a)`,
        where `query_name` is a str giving the name of the query
        sequence, and `a` is an :class:`Alignment`.

    Here is a short example:

    >>> paf_file = io.StringIO("queryname\\t10\\t0\\t10\\t+\\t"
    ...         "targetname\\t20\\t5\\t15\\t9\\t10\\t60\\t"
    ...         "cs:Z:=ATG*ga=GAACAT")
    >>> alignments = [tup for tup in parsePAF(paf_file)]
    >>> len(alignments)
    1
    >>> (queryname, alignment) = alignments[0]
    >>> queryname
    'queryname'
    >>> alignment.target
    'targetname'
    >>> (alignment.r_st, alignment.r_en)
    (5, 15)
    >>> (alignment.q_st, alignment.q_en)
    (0, 10)
    >>> alignment.strand
    1
    >>> alignment.mapq
    60
    >>> alignment.cigar_str
    '=ATG*ga=GAACAT'
    >>> alignment.q_len
    10
    """
    
    cigar_m = re.compile('cs:Z:(?P<cigar_str>'
            '(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)+)(?:\s+|$)')

    close_paf_file = False
    if isinstance(paf_file, str):
        assert os.path.isfile(paf_file), "no `paf_file` {0}".format(
                paf_file)
        paf_file = open(paf_file, 'r')
        close_paf_file = True

    elif not isinstance(paf_file, collections.Iterable):
        raise ValueError("`paf_file` must be file name or iterable")

    for line in paf_file:
        entries = line.split('\t', maxsplit=12)
        try:
            cigar_str = cigar_m.search(entries[12]).group('cigar_str')
        except:
            raise ValueError("Cannot match CIGAR:\n{0}".format(entries[12]))
        query_name = entries[0]
        a = Alignment(target=entries[5],
                      r_st=int(entries[7]),
                      r_en=int(entries[8]),
                      q_st=int(entries[2]),
                      q_en=int(entries[3]),
                      q_len=int(entries[1]),
                      mapq=int(entries[11]),
                      strand={'+':1, '-':-1}[entries[4]],
                      cigar_str=cigar_str)
        yield (query_name, a)

    if close_paf_file:
        paf_file.close()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
