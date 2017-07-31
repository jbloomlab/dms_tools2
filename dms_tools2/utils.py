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
    for modname in ['Bio', 'HTSeq']:
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
    logger.info(sessionInfo())
    logger.info("Parsed the following arguments:\n\t{0}".format(
            '\n\t'.join(['{0} = {1}'.format(arg, val) for (arg, val)
            in args.items()])))
    return logger


def iteratePairedFASTQ(r1files, r2files):
    """Iterates over FASTQ file pairs for paired-end sequencing reads.

    Args:
        `r1files` (list or str)
            Name of R1 FASTQ file or list of such files. Can optionally
            be gzipped.
        `r2files` (list or str)
            Like `r1files` but for R2 files.

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
    ...     itr = iteratePairedFASTQ(r1file.name, r2file.name)
    ...     next(itr) == (n1_1.split()[0][1 : ], r1_1, r2_1, q1_1, q2_1, False)
    ...     next(itr) == (n1_1.split()[0][1 : ], r1_1, r2_1, q1_1, q2_1, True)
    ...     next(itr) == (n1_1.split()[0][1 : ], r1_1, r2_1, q1_1, q2_1, None)
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
            yield (name1, r1, r2, q1, q2, fail)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
