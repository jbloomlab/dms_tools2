"""
===============
minimap2
===============

Runs `minimap2 <https://lh3.github.io/minimap2/>`_
aligner. 

This module is tested to work with the version of
``minimap2`` installed internally with `dms_tools2` (the
default when you initialize a :class:`Mapper` object
with `prog=None`). If you use that version of ``minimap2``,
you do not need to install ``minimap2`` separately.
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

#: `options` argument to :class:`Mapper` that works well
#: for codon-mutant libraries such as those created for
#: deep mutational scanning. Indels are highly penalized
#: as they are not expected, and settings are used that
#: effectively call strings of consecutive nucleotide
#: mutations as expected during codon mutagenesis. The
#: queries are assumed to be in the same orientation as
#: the target.
OPTIONS_CODON_DMS = ['--for-only',
                     '-A2',
                     '-B4',
                     '-O12',
                     '-E2',
                     '--secondary=no',
                     '--end-bonus=8',
                    ]

# namedtuple to hold alignments
Alignment = collections.namedtuple('Alignment',
        ['target', 'r_st', 'r_en', 'q_len', 'q_st',
         'q_en', 'strand', 'cigar_str'])
Alignment.__doc__ = "Alignment of a query to a target (reference)."
Alignment.target.__doc__ = "Target (reference) to which query was aligned."
Alignment.r_st.__doc__ = "Alignment start in target (0 based)."
Alignment.r_en.__doc__ = "Alignment end in target (0 based)."
Alignment.q_st.__doc__ = "Alignment start in query (0 based)."
Alignment.q_en.__doc__ = "Alignment end in query (0 based)."
Alignment.q_len.__doc__ = "Total length of query prior to any clipping."
Alignment.strand.__doc__ = "1 if aligns in forward polarity, -1 if in reverse."
Alignment.cigar_str.__doc__ = "CIGAR in `PAF long format <https://github.com/lh3/minimap2>`_"


class Mapper:
    """Class to run ``minimap2`` and get results.

    Args:
        `target` (str)
            FASTA file with target (reference) to which we align
            reads.
        `options` (list)
            Command line options to ``minimap2``. For 
            recommended options, for different situations, see:
                - :data:`OPTIONS_CODON_DMS`
        `prog` (str or `None`)
            Path to ``minimap2`` executable. `None` uses the
            version of ``minimap2`` installed internally with
            `dms_tools2`. This is recommended unless you have 
            some other preferred version.

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
    ...     mapper = Mapper(targetfile.name, OPTIONS_CODON_DMS)
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

    def __init__(self, target, options, prog=None):
        """See main :class:`Mapper` doc string."""
        if prog is None:
            # use default ``minimap2`` installed as package data
            prog = os.path.join(os.path.dirname(__file__),
                                'minimap2_prog')

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
            unclip_matches=True, join_gapped=True, shift_indels=True):
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
            `join_gapped` (bool)
                Pass alignments through :meth:`joinGappedAlignments`.
                Only applies to alignments returned in dict,
                does not affect any results in `outfile`.
            `shift_indels` (bool)
                Pass alignments through :meth:`shiftIndels`.
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
                dlist = collections.defaultdict(list)
                for query, alignment in parsePAF(fout):
                    dlist[query].append(alignment)
        finally:
            fout.close()
            stderr.close()

        if return_dict:
            if join_gapped or unclip_matches:
                targetseqs = {seq.name:str(seq.seq) for seq in 
                              Bio.SeqIO.parse(self.target, 'fasta')}
                queryseqs = {seq.name:str(seq.seq) for seq in
                             Bio.SeqIO.parse(queryfile, 'fasta')}
            d = {}
            for query in list(dlist.keys()):
                if len(dlist[query]) == 1:
                    d[query] = dlist[query][0]
                    del dlist[query]
                elif join_gapped:
                    a = joinGappedAlignments(dlist[query],
                            targetseqs[dlist[query][0].target],
                            queryseqs[query])
                    if a:
                        d[query] = a
                        del dlist[query]
            if dlist:
                raise ValueError("Multiple alignments for:\n\t{0}"
                        .format('\n\t'.join(dlist.keys())))
            if unclip_matches:
                for query in list(d.keys()):
                    a = d[query]
                    assert a.q_len == len(queryseqs[query])
                    a = unclipMatches(a, targetseqs[a.target], queryseqs[query])
                    d[query] = a
            if shift_indels:
                for query in list(d.keys()):
                    d[query] = d[query]._replace(cigar_str=
                            shiftIndels(d[query].cigar_str))
        else:
            d = None

        return d


def shiftIndels(cigar):
    """Shifts indels to consistent position.

    In some cases it is ambiguous where to place insertions /
    deletions in the CIGAR string. This function moves them
    to a consistent location (as far forward as possible).

    Args:
        `cigar` (str)
            PAF long CIGAR string, format is 
            `detailed here <https://github.com/lh3/minimap2>`_.

    Returns:
        A version of `cigar` with indels shifted as far
        forward as possible.

    >>> shiftIndels('=AAC-atagcc=GGG-ac=T')
    '=AA-catagc=CGGG-ac=T'

    >>> shiftIndels('=AAC-atagac=GGG-acg=AT')
    '=A-acatag=ACGG-gac=GAT'

    >>> shiftIndels('=TCC+c=TCAGA+aga=CT')
    '=T+c=CCTC+aga=AGACT'
    """
    # shift deletions
    indelmatch = re.compile('(?P<lead>=[A-Z]+)'
                          '(?P<indeltype>[\-\+])'
                          '(?P<indel>[a-z]+)'
                          '(?P<trail>=[A-Z]+)')
    i = 0
    m = indelmatch.search(cigar[i : ])
    while m:
        n = 0
        indel = m.group('indel').upper()
        while m.group('lead')[-n - 1 : ] == indel[-n - 1 : ]:
            n += 1
        if n > 0:
            if n == len(m.group('lead')) - 1:
                lead = '' # removed entire lead
            else:
                lead = m.group('lead')[ : -n]
            shiftseq = m.group('lead')[-n : ] # sequence to shift
            cigar = ''.join([
                    cigar[ : i + m.start('lead')], # sequence before match
                    lead, # remaining portion of lead
                    m.group('indeltype'),
                    shiftseq.lower(), m.group('indel')[ : -n], # new indel
                    '=', shiftseq, m.group('trail')[1 : ], # trail after indel
                    cigar[i + m.end('trail') : ] # sequence after match
                    ])
        else:
            i += m.start('trail')
        m = indelmatch.search(cigar[i : ])

    return cigar


def trimCigar(side, cigar):
    """Trims a nucleotide from CIGAR string.

    Currently just trims one site.

    Args:
        `side` (str)
            "start" trim from start, "end" to trim from end.
        `cigar` (str)
            PAF long CIGAR string, format is 
            `detailed here <https://github.com/lh3/minimap2>`_.

    Returns:
        A version of `cigar` with a single site trimmed
        from start or end.

    >>> trimCigar('start', '=ATG')
    '=TG'
    >>> trimCigar('end', '=ATG')
    '=AT'
    >>> trimCigar('start', '*ac=TG')
    '=TG'
    >>> trimCigar('end', '=AT*ag')
    '=AT'
    >>> trimCigar('start', '-aac=TG')
    '-ac=TG'
    >>> trimCigar('end', '=TG+aac')
    '=TG+aa'
    """
    if side == 'start':
        if re.match('=[A-Z]{2}', cigar):
            return '=' + cigar[2 : ]
        elif re.match('=[A-Z][\*\-\+]', cigar):
            return cigar[2 : ]
        elif re.match('\*[a-z]{2}', cigar):
            return cigar[3 : ]
        elif re.match('[\-\+][a-z]{2}', cigar):
            return cigar[0] + cigar[2 : ]
        elif re.match('[\-\+][a-z][^a-z]'):
            return cigar[2 : ]
        else:
            raise ValueError("Cannot match start of {0}".format(cigar))
    elif side == 'end':
        if re.search('[A-Z]{2}$', cigar):
            return cigar[ : -1]
        elif re.search('=[A-Z]$', cigar):
            return cigar[ : -2]
        elif re.search('\*[a-z]{2}$', cigar):
            return cigar[ : -3]
        elif re.search('[\-\+][a-z]$', cigar):
            return cigar[ : -2]
        elif re.search('[a-z]{2}$', cigar):
            return cigar[ : -1]
        else:
            raise ValueError("Cannot match end of {0}".format(cigar))
    else:
        raise ValueError("`side` must be 'start' or 'end', got {0}"
                         .format(side))


def joinGappedAlignments(alignments, target, query):
    """Join :class:`Alignment` objects of same query with long gaps.

    If a query aligns to a target with a very long gap,
    ``minimap2`` will return two entries in the PAF file,
    each of which can be captured as an :class:`Alignment`.
    This function joins them into a single :class:`Alignmnent`
    with a long gap. Generalizes to multiple alignments. Only
    joins when the alignments are separated by simple gaps, and
    the immediate flanking sequence aligns exactly.

    Args:
        `alignments` (list)
            List of :class:`Alignment` objects.
        `target` (str)
            Target sequence to which `query` is aligned.
        `query` (str)
            Query sequence that is aligned to `target`.

    Returns:
        If the alignments can be joined, return a single
        :class:`Alignment` object with the joined alignments.
        Otherwise return `None` if the alignments cannot be
        joined by this function, which may be the case if
        complex overlap with multiple types of mutations.

    Example of joining three alignments with simple gaps:

    >>> target = 'ATGCAGTCAGACATGA'
    >>> query =   'TGC  TCAG  ATG'.replace(' ', '')
    >>> a0 = Alignment(q_st=0, q_en=3, r_st=1, r_en=4, q_len=9,
    ...         cigar_str='=' + query[0 : 3], strand=1, target='target')
    >>> a1 = Alignment(q_st=3, q_en=7, r_st=6, r_en=10, q_len=9,
    ...         cigar_str='=' + query[3 : 7], strand=1, target='target')
    >>> a2 = Alignment(q_st=7, q_en=10, r_st=12, r_en=15, q_len=9,
    ...         cigar_str='=' + query[7 : 10], strand=1, target='target')
    >>> a = joinGappedAlignments([a1, a0, a2], target, query)
    >>> a.cigar_str
    '=TGC-ag=TCAG-ac=ATG'
    >>> a.q_st == 0
    True
    >>> a.q_en == len(query)
    True
    >>> a.r_st
    1
    >>> a.r_en == len(target) - 1
    True

    But we cannot join the two most distant alignments are they are not
    separated by a simple gap (there is also a gap in the query in this case):

    >>> joinGappedAlignments([a0, a2], target, query) is None
    True

    Now a more complex example where `a2` overlaps with `a1` near
    the gap, and so duplicated sequence needs to be trimmed:

    >>> a2_overlap = a2._replace(q_st=a2.q_st - 1, r_st=a2.r_st - 1,
    ...        cigar_str='=' + query[a2.q_st - 1 : a2.q_en]) 
    >>> a_overlap = joinGappedAlignments([a0, a2_overlap, a1], target, query)
    >>> a_overlap == a
    True

    Now even more complex example over overlap, where `a1` overlaps
    with `a0` via a mutation:

    >>> a0_overlap = a0._replace(q_en=a0.q_en + 1, r_en=a0.r_en + 1,
    ...         cigar_str=a0.cigar_str + '*' + target[a0.r_en].lower()
    ...                   + query[a0.q_en].lower())
    >>> a_overlap2 = joinGappedAlignments([a0_overlap, a2_overlap, a1],
    ...                                   target, query)
    >>> a_overlap2 == a
    True
    """
    assert (isinstance(alignments, collections.Iterable) and
            all(isinstance(a, Alignment) for a in alignments) and
            (len(alignments) >= 1)), \
            "`alignments` not non-empty list of `Alignment`s"

    for attr in ['q_len', 'strand', 'target']:
        if len(set(getattr(a, attr) for a in alignments)) != 1:
            raise ValueError("`alignments` don't all have same "
                             "value of `{0}`".format(attr))

    assert alignments[0].strand == 1, "only works for + strand"

    # sort by start in query
    alignments = [x[1] for x in sorted([(a.q_st, a) for a in alignments])]

    # join alignments
    while len(alignments) > 1:
        a0 = alignments[0]
        a1 = alignments[1]
        #if a0.r_en >= a1.r_st:
        #    return None # not a gap in the reference
        if a0.q_en >= a1.q_st:
            # try to resolve overlap between a1 and a0
            while a0.q_en > a1.q_st:
                noverlap = a0.q_en - a1.q_st
                if query[a0.q_en - noverlap] == target[a0.r_en - noverlap]:
                    # assign first nt of overlap to a0, not a1
                    a1 = a1._replace(
                            q_st=a1.q_st + 1,
                            r_st=a1.r_st + 1,
                            cigar_str=trimCigar('start', a1.cigar_str))
                elif (query[a1.q_st + noverlap - 1] ==
                        target[a1.r_st + noverlap - 1]):
                    # assign last nt of overlap to a1, not a0
                    a0 = a0._replace(
                            q_en=a0.q_en - 1,
                            r_en=a0.r_en - 1,
                            cigar_str=trimCigar('end', a0.cigar_str))
                else:
                    # no exact match in overlap region, suggesting a
                    # mutation complicating things... don't handle now
                    return None
            # now a simple join: a1 starts where a0 ends
            if target[a0.r_en : a1.r_st]:
                gap = '-' + target[a0.r_en : a1.r_st].lower()
                a1_cigar = a1.cigar_str
            else:
                gap = ''
                if a0.cigar_str[-1].isupper() and a1.cigar_str[0] == '=':
                    a1_cigar = a1.cigar_str[1 : ]
            cigar = a0.cigar_str + gap + a1_cigar
            a = Alignment(cigar_str=cigar,
                          q_st=a0.q_st,
                          q_en=a1.q_en,
                          r_st=a0.r_st,
                          r_en=a1.r_en,
                          q_len=a0.q_len,
                          strand=a0.strand,
                          target=a0.target
                          )
            alignments = [a] + alignments[2 : ]
        elif a0.q_en < a1.q_st:
            # gap in query between a0 and a1, might be caused
            # by a mutation near gap site or insertion in query
            return None 
        else:
            raise RuntimeError('should never get here')

    return alignments[0]            



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
    >>> a = Alignment(strand=1, r_st=1, r_en=7,
    ...         target='target', q_st=0, q_en=6, q_len=6,
    ...         cigar_str='=T*gaCAAT')
    >>> unclipMatches(a, target, query) == a
    True

    Now some soft clipping that is undone at start:

    >>> query = 'AAGCAATA'
    >>> a = Alignment(strand=1, r_st=1, r_en=7,
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
    >>> a = Alignment(strand=1, r_st=3, r_en=7,
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
                      strand={'+':1, '-':-1}[entries[4]],
                      cigar_str=cigar_str)
        yield (query_name, a)

    if close_paf_file:
        paf_file.close()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
