"""
===================
utils
===================

Miscellaneous utilities for ``dms_tools2``.
"""


import os
import sys
import time
import platform
import importlib
import logging
import tempfile
import six
import HTSeq
import dms_tools2


def sessionInfo():
    """Returns string with information about session / packages."""
    s = [
            'Version information:',
            '\tTime and date: {0}'.format(time.asctime()),
            '\tPlatform: {0}'.format(platform.platform()),
            '\tPython version: {0}'.format(
                    sys.version.replace('\n', ' ')),
            '\tdms_tools2 version: {0}'.format(dms_tools2.__version__),
            ]
    for modname in ['Bio', 'HTSeq', 'pandas']:
        try:
            v = importlib.import_module(modname).__version__
            s.append('\t{0} version: {1}'.format(modname, v))
        except ImportError:
            raise ImportError("Cannot import {0}".format(modname))
    return '\n'.join(s)


def initLogger(logfile, prog, args):
    """Initialize `logging.Logger` for scripts.

    Args:
        `logfile` (str)
            Name of file to which log is written.
        `prog` (str)
            Name of program for which we are logging.
        `args` (dict)
            Program arguments as arg / value pairs.

    Returns:
        An opened and initialized `logging.Logger`.
        Basic information about the program and args
        will have already been written to this logger.
    """
    if os.path.isfile(logfile):
        os.remove(logfile)
    logging.basicConfig(level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(prog)
    logfile_handler = logging.FileHandler(logfile)
    logger.addHandler(logfile_handler)
    formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s')
    logfile_handler.setFormatter(formatter)
    logger.info("Beginning execution of {0} in directory {1}\n".format(
            prog, os.getcwd()))
    logger.info("Progress is being logged to {0}".format(logfile))
    logger.info("{0}\n".format(sessionInfo()))
    logger.info("Parsed the following arguments:\n\t{0}\n".format(
            '\n\t'.join(['{0} = {1}'.format(arg, val) for (arg, val)
            in args.items()])))
    return logger


def iteratePairedFASTQ(r1files, r2files, r1trim=None, r2trim=None):
    """Iterates over FASTQ file pairs for paired-end sequencing reads.

    Args:
        `r1files` (list or str)
            Name of R1 FASTQ file or list of such files. Can optionally
            be gzipped.
        `r2files` (list or str)
            Like `r1files` but for R2 files.
        `r1trim` (int or `None`)
            If not `None`, trim `r1` and `q1` to be no longer than this.
        `r2trim` (int or `None`)
            Like `r1trim` but for R2.

    Returns:
        Each iteration returns `(name, r1, r2, q1, q2, fail)` where:

            - `name` is a string giving the read name

            - `r1` and `r2` are strings giving the reads

            - `q1` and `q2` are strings giving the PHRED Q scores

            - `fail` is `True` if either read failed Illumina chastity
              filter, `False` if both passed, `None` if info not present.

    We run a simple test by first writing an example FASTQ file and
    then testing on it.

    >>> n1_1 = '@DH1DQQN1:933:HMLH5BCXY:1:1101:2165:1984 1:N:0:CGATGT'
    >>> r1_1 = 'ATGCAATTG'
    >>> q1_1 = 'GGGGGIIII'
    >>> n2_1 = '@DH1DQQN1:933:HMLH5BCXY:1:1101:2165:1984 2:N:0:CGATGT'
    >>> r2_1 = 'CATGCATA'
    >>> q2_1 = 'G<GGGIII'
    >>> tf = tempfile.NamedTemporaryFile
    >>> with tf(mode='w') as r1file, tf(mode='w') as r2file:
    ...     dummyvar = r1file.write('\\n'.join([
    ...             n1_1, r1_1, '+', q1_1,
    ...             n1_1.replace(':N:', ':Y:'), r1_1, '+', q1_1,
    ...             n1_1.split()[0], r1_1, '+', q1_1,
    ...             ]))
    ...     r1file.flush()
    ...     dummyvar = r2file.write('\\n'.join([
    ...             n2_1, r2_1, '+', q2_1,
    ...             n2_1, r2_1, '+', q2_1,
    ...             n2_1, r2_1, '+', q2_1,
    ...             ]))
    ...     r2file.flush()
    ...     itr = iteratePairedFASTQ(r1file.name, r2file.name, r1trim=4, r2trim=5)
    ...     next(itr) == (n1_1.split()[0][1 : ], r1_1[ : 4], 
    ...             r2_1[ : 5], q1_1[ : 4], q2_1[ : 5], False)
    ...     next(itr) == (n1_1.split()[0][1 : ], r1_1[ : 4], 
    ...             r2_1[ : 5], q1_1[ : 4], q2_1[ : 5], True)
    ...     next(itr) == (n1_1.split()[0][1 : ], r1_1[ : 4], 
    ...             r2_1[ : 5], q1_1[ : 4], q2_1[ : 5], None)
    True
    True
    True

    """
    if isinstance(r1files, str):
        r1files = [r1files]
        assert isinstance(r2files, str)
        r2files = [r2files]
    assert len(r1files) == len(r2files) > 0
    assert isinstance(r1files, list) and isinstance(r2files, list)
    for (r1file, r2file) in zip(r1files, r2files):
        r1reader = HTSeq.FastqReader(r1file, raw_iterator=True)
        r2reader = HTSeq.FastqReader(r2file, raw_iterator=True)
        for ((r1, id1, q1, qs1), (r2, id2, q2, qs2)) in six.moves.zip(
                r1reader, r2reader):
            id1 = id1.split()
            id2 = id2.split()
            name1 = id1[0]
            name2 = id2[0]
            assert name1 == name2, "{0} vs {1}".format(name1, name2)
            # parse chastity filter assuming CASAVA 1.8 header
            fail = None
            try:
                f1 = id1[1][2]
                f2 = id2[1][2]
                if f1 == 'N' and f2 == 'N':
                    fail = False
                elif f1 in ['N', 'Y'] and f2 in ['N', 'Y']:
                    fail = True
            except IndexError:
                pass # header does not specify chastity filter
            if r1trim is not None:
                r1 = r1[ : r1trim]
                q1 = q1[ : r1trim]
            if r2trim is not None:
                r2 = r2[ : r2trim]
                q2 = q2[ : r2trim]
            yield (name1, r1, r2, q1, q2, fail)


def lowQtoN(r, q, minq):
    """Replaces low quality nucleotides with ``N`` characters.

    Args:
        `r` (str)
            A string representing a sequencing read.
        `q` (iterable)
            Iterable of same length as `r` holding Q scores.
        `minq` (int or length-one string)
            Replace all positions in `r` where `q` is < this.

    Returns:
        A version of `r` where all positions `i` where 
        `q[i] < minq` have been replaced with ``N``.

    >>> r = 'ATGCAT'
    >>> q = 'GB<.0+'
    >>> minq = '0'
    >>> qn = [38, 33, 27, 13, 15, 10]
    >>> minqn = 15
    >>> lowQtoN(r, q, minq) == lowQtoN(r, qn, minqn) == 'ATGNAN'
    True
    """
    assert len(r) == len(q)
    return ''.join([ri if qi >= minq else 'N'
            for (ri, qi) in zip(r, q)])


def buildReadConsensus(reads, minreads, minconcur):
    """Builds consensus sequence of some reads.

    You may want to pre-fill low-quality sites with ``N``
    using `lowQtoN`. An ``N`` is considered a non-called identity.

    Args:
        `reads` (list)
            List of reads as strings. If reads are not all same
            length, shorter ones are extended from 3' end with ``N``
            to match maximal length. 
        `minreads` (int)
            Only call consensus at a site if at least this many reads 
            have called identity.
        `minconcur` (float)
            Only call consensus at site if >= this fraction of called
            identities agree.

    Returns:
        A string giving the consensus sequence. Non-called 
        sites are returned as ``N```.

    >>> reads = ['ATGCAT',
    ...          'NTGNANA',
    ...          'ACGNNTAT',
    ...          'NTGNTA']
    >>> buildReadConsensus(reads, 2, 0.75) == 'ATGNNNAN'
    True
    >>> reads.append('CTGCATAT')
    >>> buildReadConsensus(reads, 2, 0.75) == 'NTGCATAT'
    True
    """
    maxlen = max(map(len, reads))
    consensus = []
    for i in range(maxlen):
        counts = {}
        for r in reads:
            if len(r) > i:
                x = r[i]
                if x != 'N':
                    if x in counts:
                        counts[x] += 1
                    else:
                        counts[x] = 1
        ntot = sum(counts.values())
        if ntot < minreads:
            consensus.append('N')
        else:
            (nmax, xmax) = sorted([(n, x) for (x, n) in counts.items()])[-1]
            if nmax / float(ntot) >= minconcur:
                consensus.append(xmax)
            else:
                consensus.append('N')
    return ''.join(consensus)


def reverseComplement(s):
    """Gets reverse complement of DNA sequence `s`.

    >>> s = 'ATGCAAN'
    >>> reverseComplement(s) == 'NTTGCAT'
    True
    """
    return ''.join(reversed([dms_tools2.NTCOMPLEMENT[nt] for nt in s]))


def alignSubamplicon(refseq, r1, r2, refseqstart, refseqend, maxmuts,
        maxN, chartype):
    """Try to align subamplicon to reference sequence at defined location.

    Tries to align reads `r1` and `r2` to `refseq` at location
    specified by `refseqstart` and `refseqend`. Determines how many
    sites of type `chartype` have mutations, and if <= `maxmuts` conside
    the subamplicon to align if fraction of ambiguous nucleotides <= `maxN`.
    In `r1` and `r2`, an ``N`` indicates a non-called ambiguous identity.
    If the reads disagree in a region of overlap that is set to ``N`` in
    the final subamplicon, but if one read has ``N`` and the other a called
    identity, then the called identity is used in the final subamplicon.

    Args:
        `refseq` (str)
            Sequence to which we align. if `chartype` is 'codon',
            must be a valid coding (length multiple of 3).
        `r1` (str)
            The forward sequence to align.
        `r2` (str)
            The reverse sequence to align. When reverse complemented,
            should read backwards in `refseq`.
        `refseqstart` (int)
            The nucleotide in `refseq` (1, 2, ... numbering) where the
            first nucleotide in `r1` aligns.
        `refseqend` (int)
            The nucleotide in `refseq` (1, 2, ... numbering) where the
            first nucleotide in `r2` aligns (note that `r2` then reads
            backwards towards the 5' end of `refseq`).
        `maxmuts` (int)
            Maximum number of mutations of character `chartype` that
            are allowed in the aligned subamplicons from the two reads.
        `maxN` (int or float)
            Maximum number of nucleotides for which we allow
            ambiguous (``N``) identities in final subamplicon.
        `chartype` (str)
            Character type for which we count mutations.
            Currently, the only allowable value is 'codon'.

    Returns:
        If reads align, return aligned subamplicon as string (of length
        `refseqend - refseqstart + 1`). Otherwise return `None`.

    >>> refseq = 'ATGGGGAAA'
    >>> s = alignSubamplicon(refseq, 'GGGGAA', 'TTTCCC', 3, 9, 1, 1, 'codon')
    >>> s == 'GGGGAAA' 
    True
    >>> s = alignSubamplicon(refseq, 'GGGGAA', 'TTTCCC', 1, 9, 1, 1, 'codon')
    >>> s == False
    True
    >>> s = alignSubamplicon(refseq, 'GGGGAT', 'TTTCCC', 3, 9, 1, 0, 'codon')
    >>> s == False
    True
    >>> s = alignSubamplicon(refseq, 'GGGGAT', 'TTTCCC', 3, 9, 1, 1, 'codon')
    >>> s == 'GGGGANA'
    True
    >>> s = alignSubamplicon(refseq, 'GGGGAT', 'TATCCC', 3, 9, 1, 0, 'codon')
    >>> s == 'GGGGATA'
    True
    >>> s = alignSubamplicon(refseq, 'GGGGAT', 'TATCCC', 3, 9, 0, 0, 'codon')
    >>> s == False
    True
    >>> s = alignSubamplicon(refseq, 'GGGNAA', 'TTTCCC', 3, 9, 0, 0, 'codon')
    >>> s == 'GGGGAAA'
    True
    >>> s = alignSubamplicon(refseq, 'GGGNAA', 'TTNCCC', 3, 9, 0, 0, 'codon')
    >>> s == 'GGGGAAA'
    True
    >>> s = alignSubamplicon(refseq, 'GTTTAA', 'TTTAAA', 3, 9, 1, 0, 'codon')
    >>> s == 'GTTTAAA' 
    True
    >>> s = alignSubamplicon(refseq, 'GGGGTA', 'TTACCC', 3, 9, 1, 0, 'codon')
    >>> s == 'GGGGTAA' 
    True
    >>> s = alignSubamplicon(refseq, 'GGGCTA', 'TTAGCC', 3, 9, 1, 0, 'codon')
    >>> s == False 
    True
    """
    assert chartype in ['codon'], "Invalid chartype"
    if chartype == 'codon':
        assert len(refseq) % 3 == 0, "refseq length not divisible by 3"

    r2 = reverseComplement(r2)

    len_subamplicon = refseqend - refseqstart + 1
    len_r1 = len(r1)
    len_subamplicon_minus_len_r2 = len_subamplicon - len(r2)
    subamplicon = []
    for i in range(len_subamplicon):
        if i < len_subamplicon_minus_len_r2: # site not in r2
            if i < len_r1: # site in r1
                subamplicon.append(r1[i])
            else: # site not in r1
                subamplicon.append('N')
        else: # site in r2
            if i < len_r1: # site in r1
                r1i = r1[i]
                r2i = r2[i - len_subamplicon_minus_len_r2]
                if r1i == r2i:
                    subamplicon.append(r1i)
                elif r1i == 'N':
                    subamplicon.append(r2i)
                elif r2i == 'N':
                    subamplicon.append(r1i)
                else:
                    subamplicon.append('N')
            else: # site not in r1
                subamplicon.append(r2[i - len_subamplicon_minus_len_r2])
    subamplicon = ''.join(subamplicon)

    if subamplicon.count('N') > maxN:
        return False

    if chartype == 'codon':
        if refseqstart % 3 == 1:
            startcodon = (refseqstart + 2) // 3
            codonshift = 0
        elif refseqstart % 3 == 2:
            startcodon = (refseqstart + 1) // 3 + 1
            codonshift = 2
        elif refseqstart % 3 == 0:
            startcodon = refseqstart // 3 + 1
            codonshift = 1
        nmuts = 0
        for icodon in range(startcodon, refseqend // 3 + 1):
            mutcodon = subamplicon[3 * (icodon - startcodon) + codonshift : 
                    3 * (icodon - startcodon) + 3 + codonshift]
            if ('N' not in mutcodon) and (mutcodon != 
                    refseq[3 * icodon - 3 : 3 * icodon]):
                nmuts += 1
                if nmuts > maxmuts:
                    return False
    else:
        raise ValueError("Invalid chartype")

    return subamplicon


def incrementCounts(refseqstart, subamplicon, chartype, counts):
    """Increment counts dict based on an aligned subamplicon.

    This is designed for keeping track of counts of different
    mutations / identities when aligning many subamplicons to
    a sequence.

    Any positions where `subamplicon` has an ``N`` are ignored,
    and not added to `counts`.

    Args:
        `refseqstart` (int)
            First nucleotide position in 1, 2, ... numbering 
            where `subamplicon` aligns.
        `subamplicon` (str)
            The subamplicon.
        `chartype` (str)
            Character type for which we are counting mutations.
            Currently, only allowable value is 'codon'.
        `counts` (dict)
            Stores counts of identities, and is incremented by
            this function. Is a dict keyed by every possible
            character (e.g., codon), with values lists with
            element `i` holding the counts for position `i`
            in 0, 1, ... numbering.

    Returns:
        On completion, `counts` has been incremented.

    >>> codonlen = 10
    >>> counts = dict([(codon, [0] * codonlen) for codon 
    ...         in dms_tools2.CODONS])
    >>> subamplicon1 = 'ATGGACTTTC'
    >>> incrementCounts(1, subamplicon1, 'codon', counts)
    >>> subamplicon2 = 'GGTCTTTCCCGGN'
    >>> incrementCounts(3, subamplicon2, 'codon', counts)
    >>> counts['ATG'][0] == 1
    True
    >>> counts['GAC'][1] == 1
    True
    >>> counts['GTC'][1] == 1
    True
    >>> counts['TTT'][2] == 2
    True
    >>> counts['CCC'][3] == 1
    True
    >>> sum([sum(c) for c in counts.values()]) == 6
    True
    """
    if chartype == 'codon':
        if refseqstart % 3 == 1:
            startcodon = (refseqstart + 2) // 3 - 1
            codonshift = 0
        elif refseqstart % 3 == 2:
            startcodon = (refseqstart + 1) // 3 
            codonshift = 2
        elif refseqstart % 3 == 0:
            startcodon = refseqstart // 3 
            codonshift = 1
    else:
        raise ValueError("Invalid chartype")

    shiftedsubamplicon = subamplicon[codonshift : ]
    for i in range(len(shiftedsubamplicon) // 3):
        codon = shiftedsubamplicon[3 * i : 3 * i + 3]
        if 'N' not in codon:
            counts[codon][startcodon + i] += 1


if __name__ == '__main__':
    import doctest
    doctest.testmod()
