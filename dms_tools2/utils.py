"""
===================
utils
===================

Miscellaneous utilities for ``dms_tools2``.
"""


import os
import math
import sys
import time
import platform
import importlib
import logging
import tempfile
import textwrap
import itertools
import collections
import random
import re

import pysam
import numpy
import scipy.misc
import scipy.special
import pandas
import gzip

import dms_tools2
from dms_tools2 import CODONS, CODON_TO_AA, AAS_WITHSTOP, AA_TO_CODONS, NTS
import dms_tools2._cutils


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
    for modname in ['Bio', 'pandas', 'numpy', 'IPython',
            'jupyter', 'matplotlib', 'plotnine', 'natsort', 'pystan',
            'scipy', 'seaborn', 'phydmslib', 'statsmodels', 'rpy2',
            'regex', 'umi_tools']:
        try:
            v = importlib.import_module(modname).__version__
            s.append('\t{0} version: {1}'.format(modname, v))
        except AttributeError:
            s.append('\t{0} version unknown'.format(modname))
        except ImportError:
            s.append("\t{0} cannot be imported".format(modname))
    return '\n'.join(s)


def initLogger(logfile, prog, args):
    """Initialize output logging for scripts.

    Args:
        `logfile` (str or `sys.stdout`)
            Name of file to which log is written, or 
            `sys.stdout` if you just want to write information
            to standard output.
        `prog` (str)
            Name of program for which we are logging.
        `args` (dict)
            Program arguments as arg / value pairs.

    Returns:
        If `logfile` is a string giving a file name, returns
        an opened and initialized `logging.Logger`. If `logfile`
        is `sys.stdout`, then writes information to `sys.stdout`.
        In either case, basic information is written about the program 
        and args.
    """
    if logfile == sys.stdout:
        logfile.write("Beginning execution of {0} in directory {1}\n\n".format(
                prog, os.getcwd()))
        logfile.write("{0}\n\n".format(sessionInfo()))
        logfile.write("Parsed the following arguments:\n\t{0}\n\n".format(
                '\n\t'.join(['{0} = {1}'.format(arg, val) for (arg, val)
                in args.items()])))
    else:
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
        try:
            logger.info("Beginning execution of {0} in directory {1}\n"
                .format(prog, os.getcwd()))
            logger.info("Progress is being logged to {0}".format(logfile))
            logger.info("{0}\n".format(sessionInfo()))
            logger.info("Parsed the following arguments:\n\t{0}\n".format(
                    '\n\t'.join(['{0} = {1}'.format(arg, val) for (arg, val)
                    in args.items()])))
        except:
            logger.exception("Error")
            raise

        return logger


def iteratePairedFASTQ(r1files, r2files, r1trim=None, r2trim=None):
    """Iterates over FASTQ files for single or paired-end sequencing.

    Args:
        `r1files` (list or str)
            Name of R1 FASTQ file or list of such files. Can optionally
            be gzipped.
        `r2files` (list or str or `None`)
            Like `r1files` but for R2 files, or `None` if no R2.
        `r1trim` (int or `None`)
            If not `None`, trim `r1` and `q1` to be no longer than this.
        `r2trim` (int or `None`)
            Like `r1trim` but for R2.

    Returns:
        Each iteration returns `(name, r1, r2, q1, q2, fail)` where:

            - `name` is a string giving the read name

            - `r1` and `r2` are strings giving the reads; `r2`
              is `None` if no R2.

            - `q1` and `q2` are strings giving the PHRED Q scores;
              `q2` is none if no R2.

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
    ...     _ = r1file.write('\\n'.join([
    ...             n1_1, r1_1, '+', q1_1,
    ...             n1_1.replace(':N:', ':Y:'), r1_1, '+', q1_1,
    ...             n1_1.split()[0], r1_1, '+', q1_1,
    ...             ]))
    ...     r1file.flush()
    ...     _ = r2file.write('\\n'.join([
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

    Now do the same test but for just R1:

    >>> with tf(mode='w') as r1file:
    ...     _ = r1file.write('\\n'.join([
    ...             n1_1, r1_1, '+', q1_1,
    ...             n1_1.replace(':N:', ':Y:'), r1_1, '+', q1_1,
    ...             n1_1.split()[0], r1_1, '+', q1_1,
    ...             ]))
    ...     r1file.flush()
    ...     itr_R1 = iteratePairedFASTQ(r1file.name, None, r1trim=4)
    ...     next(itr_R1) == (n1_1.split()[0][1 : ], r1_1[ : 4],
    ...             None, q1_1[ : 4], None, False)
    ...     next(itr_R1) == (n1_1.split()[0][1 : ], r1_1[ : 4],
    ...             None, q1_1[ : 4], None, True)
    ...     next(itr_R1) == (n1_1.split()[0][1 : ], r1_1[ : 4],
    ...             None, q1_1[ : 4], None, None)
    True
    True
    True

    """
    if isinstance(r1files, str):
        r1files = [r1files]
        if r2files is not None:
            r2files = [r2files]
    if not all(map(os.path.isfile, r1files)):
        raise ValueError('cannot find all `r1files`')
    if r2files is None:
        r2files = [None] * len(r1files)
    elif len(r1files) != len(r2files):
        raise ValueError('`r1files` and `r2files` differ in length')
    elif not all(map(os.path.isfile, r2files)):
        raise ValueError('cannot find all `r2files`')
    for (r1file, r2file) in zip(r1files, r2files):
        r1reader = pysam.FastxFile(r1file)
        if r2file is None:
            read_iterator = r1reader
        else:
            r2reader = pysam.FastxFile(r2file)
            read_iterator = zip(r1reader, r2reader)
        for tup in read_iterator:
            if r2file is None:
               a1 = tup
               r2 = q2 = None
            else:
                a1, a2 = tup
                r2 = a2.sequence
                q2 = a2.quality
                if a2.comment is not None:
                    id2 = f"{a2.name} {a2.comment}".split()
                else:
                    id2 = a2.name.split()
                name2 = id2[0]
            r1 = a1.sequence
            q1 = a1.quality
            if a1.comment is not None:
                id1 = f"{a1.name} {a1.comment}".split()
            else:
                id1 = a1.name.split()
            name1 = id1[0]
            if r2file is not None:
                # trims last two chars, need for SRA downloaded files
                if name1[-2 : ] == '.1' and name2[-2 : ] == '.2':
                    name1 = name1[ : -2]
                    name2 = name2[ : -2]
                if name1 != name2:
                    raise ValueError(f"name mismatch {name1} vs {name2}")
            # parse chastity filter assuming CASAVA 1.8 header
            try:
                f1 = id1[1][2]
                if r2file is None:
                    f2 = 'N'
                else:
                    f2 = id2[1][2]
                if f1 == 'N' and f2 == 'N':
                    fail = False
                elif f1 in ['N', 'Y'] and f2 in ['N', 'Y']:
                    fail = True
            except IndexError:
                fail = None # header does not specify chastity filter
            if r1trim is not None:
                r1 = r1[ : r1trim]
                q1 = q1[ : r1trim]
            if (r2trim is not None) and (r2file is not None):
                r2 = r2[ : r2trim]
                q2 = q2[ : r2trim]
            yield (name1, r1, r2, q1, q2, fail)


def lowQtoN(r, q, minq, use_cutils=True):
    """Replaces low quality nucleotides with ``N`` characters.

    Args:
        `r` (str)
            A string representing a sequencing read.
        `q` (str)
            String of same length as `r` holding Q scores
            in Sanger ASCII encoding.
        `minq` (length-one string)
            Replace all positions in `r` where `q` is < this.
        `use_cutils` (bool)
            Use the faster implementation in the `_cutils` module.

    Returns:
        A version of `r` where all positions `i` where 
        `q[i] < minq` have been replaced with ``N``.

    >>> r = 'ATGCAT'
    >>> q = 'GB<.0+'
    >>> minq = '0'
    >>> lowQtoN(r, q, minq) == 'ATGNAN'
    True
    """
    if use_cutils:
        return dms_tools2._cutils.lowQtoN(r, q, minq)
    assert len(r) == len(q)
    return ''.join([ri if qi >= minq else 'N'
            for (ri, qi) in zip(r, q)])


def buildReadConsensus(reads, minreads, minconcur, use_cutils=True):
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
        `use_cutils` (bool)
            Use the faster implementation in the `_cutils` module.

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
    if use_cutils:
        return dms_tools2._cutils.buildReadConsensus(reads, 
                minreads, minconcur)
    readlens = list(map(len, reads))
    maxlen = max(readlens)
    consensus = []
    for i in range(maxlen):
        counts = {}
        for (r, lenr) in zip(reads, readlens):
            if lenr > i:
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


def rarefactionCurve(barcodes, *, maxpoints=1e5, logspace=True):
    """Rarefaction curve from list of barcodes.

    Uses the analytical formula for the rarefaction curve defined
    `on Wikipedia <https://en.wikipedia.org/wiki/Rarefaction_(ecology)#Derivation>`_.

    Args:
        `barcodes` (list or pandas Series)
            Holds the list of unique barcodes for which we calculate
            the rarefaction curve. It is expected that some of these
            barcodes will be repeated multiple times in the list if
            the sampling is approaching saturation.
        `maxpoints` (int)
            Only calculate values at this many points. The benefit
            of this is that it can become very costly to calculate
            the curve at every point when there are many points.
        `logspace` (True)
            Logarithmically space the `maxpoints` points for
            the calculation. This will give better results if
            we are subsampling and the curve saturates. Only
            done if we have to subsample.

    Returns:
        The 2-tuple `(nreads, nbarcodes)`, where both `nreads` and
        `nbarcodes` are lists of the same length, and `nbarcodes[i]`
        is the expected number of barcodes observed when there are
        `nreads[i]` reads.

    Here we take a very small list and show that the results given
    by the function are equivalent to those obtained by random
    subsampling:

    >>> barcodes = ['A', 'A', 'A', 'A', 'G', 'G', 'C', 'T']
    >>> (nreads, nbarcodes) = rarefactionCurve(barcodes)
    >>> random.seed(1)
    >>> nrand = 100000
    >>> sim_equal_calc = []
    >>> for n in range(1, len(barcodes) + 1):
    ...     nbarcodes_sim = sum([len(set(random.sample(barcodes, n)))
    ...             for _ in range(nrand)]) / nrand
    ...     sim_equal_calc.append(numpy.allclose(nbarcodes_sim,
    ...             nbarcodes[nreads.index(n)], atol=1e-2))
    >>> all(sim_equal_calc)
    True
    """
    N = len(barcodes) # total number of items
    Ni = collections.Counter(barcodes)
    K = len(Ni)
    Mj = collections.Counter(Ni.values())

    Nk, num = map(numpy.array, zip(*Mj.items()))

    # use simplification that (N - Ni)Cr(n) / (N)Cr(n) =
    # [(N - Ni)! * (N - n)!] / [N! * (N - Ni - n)!]
    #
    # Also use fact that gamma(x + 1) = x!
    nbarcodes = []
    lnFactorial_N = scipy.special.gammaln(N + 1)
    if logspace and N > maxpoints:
        nreads = list(numpy.unique(numpy.logspace(
                math.log10(1), math.log10(N),
                num=int(min(N, maxpoints))).astype('int')))
    else:
        nreads = list(numpy.unique(numpy.linspace(
                1, N, num=min(N, maxpoints)).astype('int')))
    for n in nreads:
        lnFactorial_N_minus_n = scipy.special.gammaln(N - n + 1)
        i = numpy.nonzero(N - Nk - n >= 0) # indices where this is true
        nbarcodes.append(
                K - (num[i] * numpy.exp(
                            scipy.special.gammaln(N - Nk[i] + 1) +
                            lnFactorial_N_minus_n -
                            lnFactorial_N -
                            scipy.special.gammaln(N - Nk[i] - n + 1))
                    ).sum()
                )
    return (nreads, nbarcodes)


def reverseComplement(s, use_cutils=True):
    """Gets reverse complement of DNA sequence `s`.

    Args:
        `s` (str)
            Sequence to reverse complement.
        `use_cutils` (bool)
            Use the faster implementation in the `_cutils` module.

    Returns:
        Reverse complement of `s` as a str.

    >>> s = 'ATGCAAN'
    >>> reverseComplement(s) == 'NTTGCAT'
    True
    """
    if use_cutils:
        return dms_tools2._cutils.reverseComplement(s)
    return ''.join(reversed([dms_tools2.NTCOMPLEMENT[nt] for nt in s]))


def alignSubamplicon(refseq, r1, r2, refseqstart, refseqend, maxmuts,
        maxN, chartype, use_cutils=True):
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
        `maxmuts` (int or float)
            Maximum number of mutations of character `chartype` that
            are allowed in the aligned subamplicons from the two reads.
        `maxN` (int or float)
            Maximum number of nucleotides for which we allow
            ambiguous (``N``) identities in final subamplicon.
        `chartype` (str)
            Character type for which we count mutations.
            Currently, the only allowable value is 'codon'.
        `use_cutils` (bool)
            Use the faster implementation in the `_cutils` module.

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
    r2 = reverseComplement(r2)

    if use_cutils:
        return dms_tools2._cutils.alignSubamplicon(refseq, r1, r2, 
                refseqstart, refseqend, maxmuts, maxN, chartype)

    assert chartype in ['codon'], "Invalid chartype"
    if chartype == 'codon':
        assert len(refseq) % 3 == 0, "refseq length not divisible by 3"

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
    ...         in CODONS])
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


def codonToAACounts(counts):
    """Makes amino-acid counts `pandas.DataFrame` from codon counts.

    Args:
        `counts` (`pandas.DataFrame`)
            Columns are the string `site` `wildtype` and all codons
            in `CODONS`. Additional columns are allowed
            but ignored.

    Returns:
        `aacounts` (`pandas.DataFrame`)
            Columns are the string `site` and all amino acids
            in `AAS_WITHSTOP` with counts for each
            amino acid made by summing counts for encoding codons.

    >>> d = {'site':[1, 2], 'othercol':[0, 0], 'ATG':[105, 1],
    ...         'GGG':[3, 117], 'GGA':[2, 20], 'TGA':[0, 1],
    ...         'wildtype':['ATG', 'GGG']}
    >>> for codon in CODONS:
    ...     if codon not in d:
    ...         d[codon] = [0, 0]
    >>> counts = pandas.DataFrame(d)
    >>> aacounts = codonToAACounts(counts)
    >>> 'othercol' in aacounts.columns
    False
    >>> all(aacounts['site'] == [1, 2])
    True
    >>> all(aacounts['wildtype'] == ['M', 'G'])
    True
    >>> all(aacounts['M'] == [105, 1])
    True
    >>> all(aacounts['G'] == [5, 137])
    True
    >>> all(aacounts['*'] == [0, 1])
    True
    >>> all(aacounts['V'] == [0, 0])
    True
    """
    d = dict([(key, []) for key in ['site', 'wildtype'] + 
            AAS_WITHSTOP])
    for (i, row) in counts.iterrows():
        d['site'].append(row['site'])
        d['wildtype'].append(CODON_TO_AA[row['wildtype']])
        for aa in AAS_WITHSTOP:
            d[aa].append(0)
        for c in CODONS:
            d[CODON_TO_AA[c]][-1] += (row[c])
    return pandas.DataFrame(d)


def annotateCodonCounts(counts):
    """Gets annotated `pandas.DataFrame` from codon counts.

    Some of the programs (e.g., `dms2_bcsubamplicons`) create 
    ``*_codoncounts.csv`` files when run with ``--chartype codon``.
    These CSV files have columns indicating the `site` and `wildtype`
    codon, as well as a column for each codon giving the counts for that 
    codon. This function reads that file (or a `pandas.DataFrame` read
    from it) to return a `pandas.DataFrame` where a variety of additional
    useful annotations have been added.

    Args:
        `counts` (str)
            Name of existing codon counts CSV file, or `pandas.DataFrame`
            holding counts.

    Returns:
        `df` (`pandas.DataFrame`)
            The DataFrame with the information in `counts` plus
            the following added columns for each site:

                `ncounts` : number of counts at site

                `mutfreq` : mutation frequency at site

                `nstop` : number of stop-codon mutations

                `nsyn` : number of synonymous mutations

                `nnonsyn` : number of nonsynonymous mutations

                `n1nt` : number of 1-nucleotide codon mutations

                `n2nt` : number of 2-nucleotide codon mutations

                `n3nt` : number of 3-nucleotide codon mutations

                `AtoC`, `AtoG`, etc : number of each nucleotide mutation
                type among codon mutations with **one** nucleotide change.

                `mutfreq1nt`, `mutfreq2nt`, `mutfreq3nt` : frequency
                of 1-, 2-, and 3-nucleotide codon mutations at site.

    >>> d = {'site':[1, 2], 'wildtype':['ATG', 'GGG'], 'ATG':[105, 1],
    ...         'GGG':[3, 117], 'GGA':[2, 20], 'TGA':[0, 1]}
    >>> for codon in CODONS:
    ...     if codon not in d:
    ...         d[codon] = [0, 0]
    >>> counts = pandas.DataFrame(d)
    >>> with tempfile.NamedTemporaryFile(mode='w') as f:
    ...     counts.to_csv(f, index=False)
    ...     f.flush()
    ...     df = annotateCodonCounts(f.name)
    >>> all([all(df[col] == counts[col]) for col in counts.columns])
    True
    >>> all(df['ncounts'] == [110, 139])
    True
    >>> all(df['mutfreq'] == [5 / 110., 22 / 139.])
    True
    >>> all(df['nstop'] == [0, 1])
    True
    >>> all(df['nsyn'] == [0, 20])
    True
    >>> all(df['nnonsyn'] == [5, 1])
    True
    >>> all(df['n1nt'] == [0, 20])
    True
    >>> all(df['n2nt'] == [3, 2])
    True
    >>> all(df['n3nt'] == [2, 0])
    True
    >>> all(df['GtoA'] == [0, 20])
    True
    >>> all(df['AtoC'] == [0, 0])
    True
    >>> all(df['mutfreq1nt'] == [0, 20 / 139.])
    True
    >>> all(df['mutfreq3nt'] == [2 / 110., 0])
    True
    """
    if isinstance(counts, str):
        df = pandas.read_csv(counts)
    elif isinstance(counts, pandas.DataFrame):
        df = counts.copy()
    else:
        raise ValueError("invalid counts")
    assert set(CODONS) <= set(df.columns), \
            "Did not find counts for all codons".format(counts)

    df['ncounts'] = df[CODONS].sum(axis=1)

    df['mutfreq'] = (((df['ncounts'] - df.lookup(df['wildtype'].index,
            df['wildtype'].values)) / df['ncounts'].astype('float'))
            .fillna(0))

    ntchanges = ['{0}to{1}'.format(nt1, nt2) for nt1 in dms_tools2.NTS
            for nt2 in dms_tools2.NTS if nt1 != nt2]

    nstoplist = []
    nsynlist = []
    nnonsynlist = []
    nXntlists = dict([(n + 1, []) for n in range(3)])
    nntchangeslists = dict([(ntchange, []) for ntchange in ntchanges])
    for (i, row) in df.iterrows():
        nstop = nsyn = nnonsyn = 0
        nXnt = dict([(n + 1, 0) for n in range(3)])
        nntchanges = dict([(ntchange, 0) for ntchange in ntchanges])
        wt = row['wildtype']
        wtaa = CODON_TO_AA[wt]
        for c in CODONS:
            if c == wt:
                continue
            aa = CODON_TO_AA[c]
            if aa == '*':
                nstop += row[c]
            elif aa == wtaa:
                nsyn += row[c]
            else:
                nnonsyn += row[c]
            ntdiffs = ['{0}to{1}'.format(nt1, nt2) for (nt1, nt2) 
                    in zip(wt, c) if nt1 != nt2]
            nXnt[len(ntdiffs)] += row[c]
            if len(ntdiffs) == 1:
                nntchanges[ntdiffs[0]] += row[c]
        nstoplist.append(nstop)
        nsynlist.append(nsyn)
        nnonsynlist.append(nnonsyn)
        for n in range(3):
            nXntlists[n + 1].append(nXnt[n + 1])
        for ntchange in ntchanges:
            nntchangeslists[ntchange].append(nntchanges[ntchange])
    df = df.assign(nstop=nstoplist, nsyn=nsynlist, nnonsyn=nnonsynlist)
    df = df.assign(n1nt=nXntlists[1], n2nt=nXntlists[2], n3nt=nXntlists[3])
    for ntchange in ntchanges:
        df[ntchange] = nntchangeslists[ntchange]

    for nnt in range(3):
        df['mutfreq{0}nt'.format(nnt + 1)] = (df['n{0}nt'.format(nnt + 1)]
            / df['ncounts'].astype('float')).fillna(0)

    return df


def adjustErrorCounts(errcounts, counts, charlist, maxexcess):
    """Adjust error counts to not greatly exceed counts of interest.

    This function is useful when estimating preferences. Under the
    model, the error-control should not have a higher rate of error
    than the actual sample. However, this could happen if the experimental
    data don't fully meet the assumptions. So this function scales
    down the error counts in that case.

    Args:
        `errcounts` (pandas.DataFrame)
            Holds counts for error control.
        `counts` (pandas.DataFrame)
            Holds counts for which we are correcting errors.
        `charlist` (list)
            Characters for which we have counts.
        `maxexcess` (int)
            Only let error-control counts exceed actual by this much.

    Returns:
        A copy of `errcounts` except for any non-wildtype character,
        the maximum frequency of that character is adjusted to be
        at most the number predicted by the frequency in `counts`
        plus `maxexcess`.

    >>> counts = pandas.DataFrame({'site':[1], 'wildtype':['A'], 
    ...         'A':500, 'C':10, 'G':40, 'T':20})
    >>> errcounts = pandas.DataFrame({'site':[1], 'wildtype':['A'],
    ...         'A':250, 'C':1, 'G':30, 'T':10})
    >>> charlist = ['A', 'C', 'G', 'T']
    >>> errcounts = errcounts[['site', 'wildtype'] + charlist]
    >>> adj_errcounts = adjustErrorCounts(errcounts, counts, charlist, 1)
    >>> set(adj_errcounts.columns) == set(errcounts.columns)
    True
    >>> all(adj_errcounts['site'] == errcounts['site'])
    True
    >>> all(adj_errcounts['wildtype'] == errcounts['wildtype'])
    True
    >>> (adj_errcounts[adj_errcounts['site'] == 1][charlist].values[0]
    ...         == numpy.array([250, 1, 21, 10])).all()
    True
    """
    cols = counts.columns
    counts = counts.sort_values('site')
    errcounts = errcounts.sort_values('site')
    assert all(counts['site'] == errcounts['site'])
    assert all(counts['wildtype'] == errcounts['wildtype'])
    counts['total'] = counts[charlist].sum(axis=1).astype('float')
    errcounts['total'] = errcounts[charlist].sum(axis=1)
    maxallowed = (counts[charlist].div(counts['total'], axis=0).multiply(
            errcounts['total'], axis=0) + maxexcess).round().astype('int')
    adj_errcounts = errcounts[charlist].where(errcounts[charlist] < maxallowed, 
            maxallowed[charlist])
    for c in charlist:
        adj_errcounts[c] = adj_errcounts[c].where(counts['wildtype'] != c,
                errcounts[c])
    for col in cols:
        if col not in charlist:
            adj_errcounts[col] = counts[col]
    return adj_errcounts[cols]


def convertCountsFormat(oldfile, newfile, charlist):
    """Convert counts file from ``dms_tools`` to ``dms_tools2`` format.

    Args:
        `oldfile` (str)
            Name of counts file in the old ``dms_tools`` format:
            http://jbloomlab.github.io/dms_tools/fileformats.html
        `newfile` (str)
            Name of created counts file in the ``dms_tools2`` format:
            https://jbloomlab.github.io/dms_tools2/dms2_bcsubamp.html
        `charlist` (list)
            List of characters that we expect in the counts files.
            For instance, could be `CODONS`.
    """
    with open(oldfile) as f:
        header = f.readline()
    assert header[0] == '#'
    cols = header[1 : ].split()
    assert cols[0] == 'POSITION' and cols[1] == 'WT'
    cols = ['site', 'wildtype'] + cols[2 : ]
    assert set(charlist) == set(cols[2 : ])
    old = pandas.read_csv(oldfile, delim_whitespace=True, 
            names=cols, comment='#')
    old.to_csv(newfile, index=False)


def renumberSites(renumbfile, infiles, missing='error',
        outfiles=None, outprefix=None, outdir=None):
    """Renumber sites in CSV files.

    Switch numbering scheme in files with a column named `site`. 

    You must specify **exactly one** of `outfiles`,
    `outprefix`, and `outdir` as something other than `None`.

    Args:
        `renumbfile` (str)
            Name of existing CSV file with the re-numbering scheme.
            Should have columns with name `original` and `new`.
            Each entry in `original` should refer to a site in 
            the input files, and each entry in `new` should be
            the new number for this site. If an entry in `new`
            is `None` or `nan` then it is dropped from the newly
            numbered files regardless of `missing`.
        `infiles` (list)
            List of existing CSV files that we are re-numbering.
            Each file must have an entry of `site`.
        `missing` (str)
            How to handle sites in `infiles` but not `renumbfile`.
                - `error`: raise an error
                - `skip`: skip renumbering, leave with original number 
                - `drop`: drop any sites not in `renumbfile`
        `outfiles` (list)
            List of output files of the same length as `infiles`.
            The numbered version of `infiles` is named as the
            corresponding entry in `outfiles`.
        `outdir` (str)
            A directory name. The renumbered files have the same
            names as in `infile`, but are now placed in `outdir`.
        `outprefix` (str)
            The renumbered files have the same names and locations
            as `infiles`, but have the pre-pended filename extension
            `outprefix`.
    """
    assert os.path.isfile(renumbfile), "no renumbfile {0}".format(renumbfile)
    renumb = pandas.read_csv(renumbfile)
    assert {'original', 'new'} <=  set(renumb.columns), \
            "renumbfile lacks columns `original` and/or `new`"
    for col in ['original', 'new']:
        assert len(renumb[col]) == len(set(renumb[col])), \
                "duplicate sites for {0} in {1}".format(col, renumbfile)
        renumb[col] = renumb[col].astype('str')

    assert isinstance(infiles, list), "infiles is not a list"
    nin = len(infiles)
    infiles = [os.path.abspath(f) for f in infiles]
    assert len(set(infiles)) == nin, "duplicate files in `infiles`"

    if outfiles is not None:
        assert isinstance(outfiles, list), "`outfiles` not list"
        assert (outdir is None) and (outprefix is None), \
                "only specify one of `outfiles`, `outdir`, and `outprefix`"
        nout = len(outfiles)
        assert nout == nin, "`outfiles` and `infiles` different length"

    elif outdir is not None:
        assert isinstance(outdir, str), "`outdir` should be string"
        assert (outfiles is None) and (outprefix is None), \
                "only specify one of `outfiles`, `outdir`, and `outprefix`"
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        outfiles = [os.path.join(outdir, os.path.basename(f))
                for f in infiles]

    elif outprefix is not None:
        assert isinstance(outprefix, str), "`outdir` should be string"
        assert (outfiles is None) and (outdir is None), \
                "only specify one of `outfiles`, `outdir`, and `outprefix`"
        outfiles = [os.path.join(os.path.dirname(f), outprefix + 
                os.path.basename(f)) for f in infiles]

    else:
        raise ValueError("specify `outdir`, `outprefix`, `outfiles`")

    outfiles = [os.path.abspath(f) for f in outfiles]
    assert len(set(outfiles)) == len(outfiles), "duplicate files in `outfiles`"
    assert not set(outfiles).intersection(set(infiles)), \
            "some in and outfiles the same"

    for (fin, fout) in zip(infiles, outfiles):
        df_in = pandas.read_csv(fin)
        assert 'site' in df_in.columns, "no `site` column in {0}".format(fin)
        df_in['site'] = df_in['site'].astype('str')
        if missing == 'error':
            if set(df_in['site']) > set(renumb['original']):
                raise ValueError("`missing` is `error`, excess sites in {0}"
                    .format(fin))
        elif missing == 'skip':
            pass
        elif missing == 'drop':
            df_in = df_in[df_in['site'].isin(renumb['original'])]
        else:
            raise ValueError("invalid `missing` of {0}".format(missing))

        # can't just use replace below because of this bug:
        # https://github.com/pandas-dev/pandas/issues/16051
        unmappedsites = df_in[~df_in['site'].isin(renumb['original'])]['site']
        replacemap = dict(zip(
                renumb['original'].append(unmappedsites),
                renumb['new'].append(unmappedsites)))
        df_in['site'] = df_in['site'].map(replacemap)

        df_in = (df_in[df_in['site'].notnull()]
                      .query('site != "NaN"')
                      .query('site != "nan"')
                      .query('site != "None"')
                      )

        df_in.to_csv(fout, index=False)


def codonEvolAccessibility(seqs):
    """Accessibility of amino acids by nucleotide mutations.

    Args:
        `seqs` (str or list)
            A single coding sequence or a list of such sequences.

    Returns:
        A pandas DataFrame listing all sites in the sequence(s)
        numbered 1, 2, ..., with columns giving the accessibility
        of each amino acid by single nucleotide mutations.

    The accessibility of codon :math:`c` to amino-acid :math:`a`
    by single-nucleotide mutations is defined as the minimum
    number of nucleotide mutations needed to generate that
    amino-acid. 

    For a collection of sequences, we calculate the 
    accessibility as the weighted average of the accessibilities
    of all codons observed at that site in the collection of
    sequences.

    As an example, compute accessibility for one sequence:

    >>> s = "ATGGGA"
    >>> acc = codonEvolAccessibility(s)

    The returned pandas DataFrame `acc` is has a column named
    `site` plus columns for all amino acids:

    >>> all(acc.columns == ['site'] + AAS_WITHSTOP)
    True

    We look at entries for a few amino acids. At the first 
    site, the wildtype entry in the sequence `s` is the codon
    for *M* (``ATG``). So at this site, the distance to *M*
    is 0. The distance to *I* (which has codon ``ATA`` as a
    codon) is 1, and the distance to *W* (which has only ``TGG``
    as a codon) is 2.

    >>> acc[['site', 'G', 'I', 'M', 'W']]
       site    G    I    M    W
    0     1  2.0  1.0  0.0  2.0
    1     2  0.0  2.0  3.0  2.0

    If we pass the function a list of multiple sequences,
    then the accessibilities are averaged over the sequences:

    >>> acc2 = codonEvolAccessibility(['ATGGGA', 'ATAGGA'])
    >>> acc2[['site', 'G', 'I', 'M', 'W']]
       site    G    I    M    W
    0     1  2.0  0.5  0.5  2.5
    1     2  0.0  2.0  3.0  2.0
    """
    # get number of nucleotide diffs between all pairs of codons
    nt_diffs = dict([
            ((c1, c2), sum(1 for x1, x2 in zip(c1, c2) if x1 != x2))
            for c1, c2 in itertools.product(CODONS, repeat=2)])

    # get number of nucleotide diffs to nearest codon for amino acid
    aa_nt_diffs = {}
    for c in CODONS:
        for aa, othercs in AA_TO_CODONS.items():
            aa_nt_diffs[(c, aa)] = min([nt_diffs[(c, c2)] 
                    for c2 in othercs])

    # make sure seqs are of same valid length
    if isinstance(seqs, str):
        seqs = [seqs]
    assert len(seqs[0]) % 3 == 0, "seqs not of length divisible by 3"
    assert all([len(seqs[0]) == len(s) for s in seqs[1 : ]]), \
            "seqs not all of same length"

    # get nucleotide distances, summing for all sequences
    dists = collections.defaultdict(lambda: collections.defaultdict(float))
    for s in seqs:
        for r in range(len(s) // 3):
            c = s[3 * r : 3 * r + 3]
            assert c in CODONS, "invalid codon {0}".format(c)
            for aa in AAS_WITHSTOP:
                dists[r + 1][aa] += aa_nt_diffs[(c, aa)]

    return (pandas.DataFrame.from_dict(dists, orient='index')
            .rename_axis('site')
            [AAS_WITHSTOP]
            / len(seqs)).reset_index()


def sigFigStr(x, nsig):
    """Get str of `x` with `nsig` significant figures.

    >>> sigFigStr(11190, 2)
    '11000'
    >>> sigFigStr(117, 2)
    '120'
    >>> sigFigStr(6, 2)
    '6.0'
    >>> sigFigStr(0.213, 2)
    '0.21'
    >>> sigFigStr(0.007517, 3)
    '0.00752'
    """
    if x <= 0:
        raise ValueError('currently only handles numbers > 0')
    x = float(f"{{:.{nsig}g}}".format(x))
    if x >= 10**(nsig - 1):
        return '{:d}'.format(round(x))
    else:
        predecimal = math.floor(math.log10(x)) + 1
        postdecimal = nsig - predecimal
        assert postdecimal > 0, str(x)
        return f"{{:.{postdecimal}f}}".format(x)

def getSubstitutions(wildtype, mutant, amino_acid=False):
    """Get space delimited string of substitutions

    Args:
        `wildtype` (str):
             The wildtype sequence
        `mutant` (str):
             The mutant sequence
        `amino_acid` (bool)
             Specify whether the sequence is amino acid.
             Default is False
    Returns:
        A space delimited string of substitutions present in the
        mutant sequence

    >>> getSubstitutions('AGT', 'TGT')
    'A1T'
    >>> getSubstitutions('AAGTAACGA', 'ATCTAACGA')
    'A2T G3C'
    >>> getSubstitutions('TYARV', 'GYAGV', amino_acid=True)
    'T1G R4G'
    """
    if len(wildtype) != len(mutant):
        raise ValueError('wildtype and mutant must be same length')
    subs = []
    for site in range(len(wildtype)):
        wt = wildtype[site]
        mut = mutant[site]
        if amino_acid:
            if wt not in AAS_WITHSTOP:
                raise ValueError (f"Invalid wt residue {wt} at site {site+1}")
            if mut not in AAS_WITHSTOP:
                raise ValueError (f"Invalid mutant residue {mut} at site {site+1}")
        else:
            if wt not in NTS:
                raise ValueError (f"Invalid wt nucleotide {wt} at site {site+1}")
            if mut not in NTS:
                raise ValueError (f"Invalid mutant nucleotide {mut} at site {site+1}")
        if wt!=mut:
            pos = str(site + 1)
            subs.append(f"{wt}{pos}{mut}")
    subs = ' '.join(subs)

    return subs


def codon_to_nt_counts(codoncounts):
    """Convert codon counts file to nucleotide counts.

    Args:
        `codoncounts` (str or pandas.DataFrame)
            Codon counts in format produced by ``dms2_bcsubamp``,
            either as CSV file or data frame holding CSV.

    Returns:
        pandas.DataFrame with nucleotide counts.

    Example:

    >>> with tempfile.NamedTemporaryFile('w') as f:
    ...     _ = f.write(textwrap.dedent('''
    ...             site,wildtype,AAA,AAC,AAG,AAT,ACA,ACC,ACG,ACT,AGA,AGC,AGG,AGT,ATA,ATC,ATG,ATT,CAA,CAC,CAG,CAT,CCA,CCC,CCG,CCT,CGA,CGC,CGG,CGT,CTA,CTC,CTG,CTT,GAA,GAC,GAG,GAT,GCA,GCC,GCG,GCT,GGA,GGC,GGG,GGT,GTA,GTC,GTG,GTT,TAA,TAC,TAG,TAT,TCA,TCC,TCG,TCT,TGA,TGC,TGG,TGT,TTA,TTC,TTG,TTT
    ...             1,ATG,0,0,0,0,0,0,2,0,0,0,0,0,8,0,333985,14,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    ...             2,AAG,16,20,333132,41,13,12,27,14,8,6,67,8,9,13,29,9,10,11,12,8,10,15,15,11,6,9,3,7,8,10,17,4,3,7,49,7,9,14,9,4,10,7,7,7,9,11,11,5,14,14,11,6,13,16,15,14,9,9,15,8,9,11,8,15
    ...             3,GCA,2,3,8,3,34,11,7,6,7,6,9,8,4,3,5,0,6,14,10,12,6,8,7,10,5,11,7,6,6,1,3,12,19,6,11,9,333250,10,6,9,15,3,5,5,37,9,9,7,8,4,8,3,23,5,7,8,6,11,7,10,7,9,3,6 
    ...             '''.strip()))
    ...     f.flush()
    ...     nt_counts = codon_to_nt_counts(f.name)
    >>> nt_counts
       site wildtype       A       C       G       T
    0     1        A  334009       0       6       0
    1     2        T       0       2       0  334013
    2     3        G       8       0  333993      14
    3     4        A  333424     156     169     187
    4     5        A  333361     211     186     178
    5     6        G     156     185  333427     168
    6     7        G     116     124  333410     125
    7     8        C     126  333407     121     121
    8     9        A  333435     114     112     114

    """
    if not isinstance(codoncounts, pandas.DataFrame):
        codoncounts = pandas.read_csv(codoncounts)

    if codoncounts['site'].dtype != int:
        raise ValueError('`site` column in `codoncounts` must be integer')

    nt_counts = []
    for i_nt in [0, 1, 2]:
        nt_counts.append(
                codoncounts
                .melt(id_vars=['site', 'wildtype'],
                      var_name='codon',
                      value_name='count',
                      )
                .assign(
                    site=lambda x: 3 * (x['site'] - 1) + i_nt + 1,
                    wildtype=lambda x: x['wildtype'].str[i_nt],
                    nucleotide=lambda x: x['codon'].str[i_nt],
                    )
                .groupby(['site', 'wildtype', 'nucleotide'])
                .aggregate({'count': 'sum'})
                .reset_index()
                .pivot_table(values='count',
                             columns='nucleotide',
                             index=['site', 'wildtype'])
                .reset_index()
                )

    nt_counts = (pandas.concat(nt_counts)
                 .sort_values('site')
                 .reset_index(drop=True)
                 )
    del nt_counts.columns.name
    return nt_counts


def barcodeInfoToCodonVariantTable(samples, geneseq, path=None):
    """Convert barcode info files into a CodonVariantTable

    Convert barcode info files output from `dms2_bcsubamp` into a
    `CodonVariantTable`. Barcode info files contain reads and barcodes from
    barcoded subamplicon sequencing, described 
    `here <https://jbloomlab.github.io/dms_tools2/bcsubamp.html>`_.
    This function takes consensus reads retained by `dms2_bcsubamp`, 
    gives each unique sequence a numerical barcode (since the barcodes from 
    `dms2_bcsubamp` could come from the same variant), and counts the number
    of retained consensus reads corresponding to each sequence. Then, a
    `CodonVariantTable` is made using the sequences and their numerical 
    barcodes, and counts are added based on the number of retained consensus
    reads of those sequences. Therefore, the `CodonVariantTable` will only 
    contain one 'variant' for each unique sequence with the total count for all
    the unbarcoded variants in the experiment which had the same sequence.

    Args:
        `samples` (dict):
            Dictionary with libraries as keys and lists of info file prefixes
            (file names without the '_bcinfo.txt.gz') for files corresponding
            to those libraries as values.
            
            Example: {'library-1':['condition-1-library-1'],
                      'library-2':['condition-1-library-2']}
        `geneseq` (str):
            The wildtype gene sequence
        `path` (str)
            Directory in which barcode info files are located

    Returns:
        A `dms_variants.codonvarianttable.CodonVariantTable` with 'counts'
        generated from the barcode info files
    """

    # Set up re matchers for looking at lines
    matcher = re.compile(r'(?P<linetype>^.*\:) '
                         r'(?P<contents>.*$)')

    alt_matcher = re.compile(r'(?P<linetype>^R\d READS:$)')

    read_matcher = re.compile(r'(?P<read>^[ATGCN\s]*$)')

    # Create a dictionary to contain dictionaries of each library's barcodes
    libraries = {}

    # Initialize lists for making the codonvarianttable
    barcodes = []
    subs = []
    variant_call_support = []
    library_list = []

    # For each library, go through each sample file and collect data
    for library in samples.keys():
        # Initialize dictionary to contain this library's reads and barcodes
        barcode_dictionary = {}
        # Start a barcode count for this library
        cur_barcode = 1
        # For each barcode info file corresponding to a sample in this library
        for sample in samples[library]:

            # Set initial conditions
            take_next = False
            description_skipped = False

            # Find the file
            f = f"{sample}_bcinfo.txt.gz"
            if path:
                file_path = os.path.join(os.path.abspath(path), f)
            else:
                file_path = f

            # Open the file and loop through it to find retained consensus
            # reads and give them each a new barcode
            with gzip.open(file_path, 'r') as f:
                # Make sure the first line looks like it is supposed to
                firstline = f.readline()
                firstline = firstline.decode()
                first_match = matcher.match(firstline)
                if first_match.group('linetype') != 'BARCODE:':
                    raise ValueError(f"Unexpected first line {firstline}: may be "
                    "unexpected file type")
                else:
                    previous_line = first_match

                # Go through the lines, making they are in the expected order
                for line in f:
                    line = line.decode()
                    line_match = matcher.match(line)
                    if not line_match:
                        line_match = alt_matcher.match(line)
                    if not line_match:
                        read_match = read_matcher.match(line)
                        if not read_match:
                            raise ValueError(f"Unable to recognize line {line}")
                        else:
                            line_is_read = True
                            previous_linetype = previous_line.group('linetype')
                            if previous_linetype != 'R1 READS:' and \
                               previous_linetype != 'R2 READS:':
                               raise ValueError(f"Unexpected line {line}")
                    else:
                        line_is_read = False
                    if previous_line.group('linetype') == 'BARCODE:':
                        if line_match.group('linetype') != 'RETAINED:':
                            raise ValueError(f"Unexpected line {line}")
                        # Decide whether to retain the next consensus or not
                        else:
                            if line_match.group('contents') == 'False':
                                retain = False
                            elif line_match.group('contents') == 'True':
                                retain = True
                            else:
                                raise ValueError(f"Unexpected line {line}")
                    elif previous_line.group('linetype') == 'RETAINED:':
                        if line_match.group('linetype') != 'DESCRIPTION:':
                            raise ValueError(f"Unexpected line {line}")
                    elif previous_line.group('linetype') == 'DESCRIPTION:':
                        if line_match.group('linetype') != 'CONSENSUS:':
                            raise ValueError(f"Unexpected line {line}")
                        # Make sure we know whether to retain or not
                        elif not isinstance(retain, bool):
                            raise ValueError(
                            f"Unclear whether to retain {line_match.group('contents')}"
                            )
                        elif retain:
                            read = line_match.group('contents')
                            # Add the read to the dictionary if not in it
                            # Also give it a barcode
                            if 'N' not in read:
                                if read not in barcode_dictionary:
                                    # Create the sequence in the dictionary
                                    barcode_dictionary[read] = {}
                                    # Give it an initial count of 1 for this sample
                                    barcode_dictionary[read][sample] = 1
                                    # Give it the next barcode
                                    barcode_dictionary[read]['barcode'] = cur_barcode
                                    # Save values for making CodonVariantTable
                                    barcodes.append(cur_barcode)
                                    subs.append(getSubstitutions(geneseq, read))
                                    variant_call_support.append(1)
                                    library_list.append(library)
                                    # Advance current barcode
                                    cur_barcode += 1
                                else:
                                    # Add a counter for the sample if sequence
                                    # not seen for this sample yet
                                    if sample not in barcode_dictionary[read]:
                                        barcode_dictionary[read][sample] = 1
                                    else:
                                        # Add another count to this read for
                                        # this sample
                                        barcode_dictionary[read][sample] += 1
                        # Set retain to None
                        retain = None
                    elif previous_line.group('linetype') == 'CONSENSUS:':
                        if line_match.group('linetype') != 'R1 READS:':
                            raise ValueError(f"Unexpected line {line}")
                    elif previous_line.group('linetype') == 'R1 READS:':
                        if not line_is_read:
                            if line_match.group('linetype') != 'R2 READS:':
                                raise ValueError(f"Unexpected line {line}")
                    elif previous_line.group('linetype') == 'R2 READS:':
                        if not line_is_read:
                            if line_match.group('linetype') != 'BARCODE:':
                                raise ValueError(f"Unexpected line {line}")
                    # Save this line as the previous line if it is not a read
                    if not line_is_read:
                        previous_line = line_match
        # After going through each file for a library, save its dictionary with
        # reads and barcodes
        libraries[library] = barcode_dictionary

    # Make the dataframe for creating the codonvarianttable
    df =  {'barcode':barcodes,
            'substitutions':subs,
            'library':library_list,
            'variant_call_support':variant_call_support,
           }
    df = pandas.DataFrame(df)

    # Make the codonvarianttable
    with tempfile.NamedTemporaryFile(mode='w') as f:
        df.to_csv(f, index=False)
        f.flush()
        variants = dms_variants.codonvarianttable.CodonVariantTable(
                    barcode_variant_file=f.name,
                    geneseq=geneseq)

    # Make the counts dataframe:
    # Initialize list of dataframes
    dfs = []

    # Loop through each library and produce count dataframes for each sample
    for library in libraries:
        barcode_dictionary = libraries[library]
        for sample in samples[library]:
            barcodes_list = []
            counts_list = []
            sample_list = []
            library_list = []
            # Get counts for this sample
            for sequence in barcode_dictionary.keys():
                if sample not in barcode_dictionary[sequence].keys():
                    counts_list.append(0)
                else:
                    counts_list.append(barcode_dictionary[sequence][sample])
                barcodes_list.append(barcode_dictionary[sequence]['barcode'])
                sample_list.append(sample)
                library_list.append(library)
            # Make a dataframe for this sample
            data = {'barcode':barcodes_list,
                    'count':counts_list,
                    'sample':sample_list,
                    'library':library_list,
                   }
            data = pandas.DataFrame(data)
            # Append it to the list of dataframes
            dfs.append(data)
        # Concatenate the list of dataframes into a counts dataframe
        barcode_counts = pandas.concat(dfs)
    # Add the counts for each sample to the codonvarianttable
    for library in libraries:
        for sample in samples[library]:
            icounts = barcode_counts.query('library == @library & sample == @sample')
            icounts = icounts[['barcode', 'count']]
            variants.addSampleCounts(library, sample, icounts)

    return(variants)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
