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
import pandas
import HTSeq
import dms_tools2
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
    for modname in ['six', 'Bio', 'HTSeq', 'pandas', 'numpy', 'IPython',
            'plotnine']:
        try:
            v = importlib.import_module(modname).__version__
            s.append('\t{0} version: {1}'.format(modname, v))
        except ImportError:
            raise ImportError("Cannot import {0}".format(modname))
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
            # next check trims last two chars, need for SRA downloaded files
            if name1[-2 : ] == '.1' and name2[-2 : ] == '.2':
                name1 = name1[ : -2]
                name2 = name2[ : -2]
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


def annotateCodonCounts(countsfile):
    """Gets annotated `pandas` DataFrame from ``*_codoncounts.csv`` files.

    Some of the programs (e.g., `dms2_bcsubamplicons`) create 
    ``*_codoncounts.csv`` files when run with ``--chartype codon``.
    These CSV files have columns indicating the `site` and `wildtype`
    codon, as well as a column for each codon giving the counts for that 
    codon. This function reads that file to return a `pandas` DataFrame
    where a variety of additional useful annotations have been added
    for each site.

    Args:
        `countsfile` (str)
            Name of existing codon counts file.

    Returns:
        `df` (`pandas.DataFrame`)
            The DataFrame with the information in `countsfile` plus
            the following added columns for each site:

                `ncounts` : number of counts at site

                `mutfreq` : mutation frequency at site

                `nstop` : number of stop-codon mutations

                `nsyn` : number of synonymous mutations

                `nnonsyn` : number of nonsynonymous mutations

                `n1nt` : number of 1-nucleotide codon mutations

                `n2nt` : number of 2-nucleotide codon mutations
                
                `n3nt` : number of 3-nucleotide codon mutations

    >>> d = {'site':[1, 2], 'wildtype':['ATG', 'GGG'], 'ATG':[105, 1],
    ...         'GGG':[3, 117], 'GGA':[2, 20], 'TGA':[0, 1]}
    >>> for codon in dms_tools2.CODONS:
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
    """
    df = pandas.read_csv(countsfile)
    assert set(dms_tools2.CODONS) <= set(df.columns), \
            "{0} does not have counts for all codons".format(countsfile)

    df['ncounts'] = df[dms_tools2.CODONS].sum(axis=1)

    df['mutfreq'] = (df['ncounts'] - df.lookup(df['wildtype'].index,
            df['wildtype'].values)) / df['ncounts'].astype('float')

    nstoplist = []
    nsynlist = []
    nnonsynlist = []
    nXntlists = dict([(n + 1, []) for n in range(3)])
    for (i, row) in df.iterrows():
        nstop = nsyn = nnonsyn = 0
        nXnt = dict([(n + 1, 0) for n in range(3)])
        wt = row['wildtype']
        wtaa = dms_tools2.CODON_TO_AA[wt]
        for c in dms_tools2.CODONS:
            if c == wt:
                continue
            aa = dms_tools2.CODON_TO_AA[c]
            if aa == '*':
                nstop += row[c]
            elif aa == wtaa:
                nsyn += row[c]
            else:
                nnonsyn += row[c]
            nXnt[sum([nt1 != nt2 for (nt1, nt2) in zip(wt, c)])] += row[c]
        nstoplist.append(nstop)
        nsynlist.append(nsyn)
        nnonsynlist.append(nnonsyn)
        for n in range(3):
            nXntlists[n + 1].append(nXnt[n + 1])
    df = df.assign(nstop=nstoplist, nsyn=nsynlist, nnonsyn=nnonsynlist)
    df = df.assign(n1nt=nXntlists[1], n2nt=nXntlists[2], n3nt=nXntlists[3])

    return df


if __name__ == '__main__':
    import doctest
    doctest.testmod()
