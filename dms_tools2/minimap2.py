"""
===============
minimap2
===============

Runs `minimap2 <https://lh3.github.io/minimap2/>`_
aligner. Although this aligner comes with it's own
Python API (`mappy`), it does not support all option.

Therefore, this module interacts directly with
the ``minimap2`` executable on the command line.
"""


import os
import re
import subprocess
import tempfile
import collections


#: ``minimap2`` settings that work well for a query
#: being matched to a short target where the target
#: is known to be in the same orientation at the query.
ORIENTED_READ = ['-uf', 
                 '-k15',
                 '-w5',
                 '-g2000',
                 '-G200k',
                 '-A1',
                 '-B2',
                 '-O2,32',
                 '-E1,0',
                 '-z200',
                 '--secondary=no']



class MapperCmdLine:
    """Class to run ``minimap2`` at command line and get results.

    Known to work with ``minimap2`` version 2.10.

    Args:
        `prog` (str)
            Path to ``minimap2`` executable.
        `options` (list)
            Command line options to ``minimap2``.

    Attributes:
        `prog` (str)
            Path to ``minimap2`` set at initialization.
        `options` (list)
            Options to ``minimap2`` set at initialization.
        `version` (str)
            Version of ``minimap2``.
    """

    def __init__(self, prog='minimap2', options=ORIENTED_READ):
        """See main :class:`MapperCmdLine` doc string."""
        try:
            version = subprocess.check_output([prog, '--version'])
        except:
            raise ValueError("Can't execute `prog` {0}".format(prog))
        self.version = version.strip().decode('utf-8')
        self.prog = prog
        self.options = options


    def align(self, target, query, return_dict=True, outfile=None):
        """Align sequences.

        Aligns query sequences to a target. Adds ``--c --cs=long``
        arguments to `options` to get a long CIGAR string, and
        returns the results as a dictionary and/or writes them
        to a PAF file.

        Args:
            `target` (str)
                FASTA file with target to which we align.
            `query` (str)
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

        Returns:
            Either a dict (if `return_dict` is `True`)
            or `None` (if `return_dict` is `False`). If returning
            a dictionary, the keys are the header for each sequence
            in `query`, and the values are alignments in the format
            returned by :meth:`parsePAF`.
        """
        assert os.path.isfile(target), "no `target` {0}".format(target)
        assert os.path.isfile(query), "no `query` {0}".format(query)

        assert '-a' not in self.options, \
                "output should be PAF format, not SAM"
        for arg in ['-c', '--cs=long']:
            if arg not in self.options:
                self.options.append(arg)

        if outfile is None:
            fout = tempfile.TemporaryFile()
        else:
            fout = open(outfile, 'w')
        stderr = tempfile.TemporaryFile()
        try:
            _ = subprocess.check_call(
                    [self.prog] + args + [target, query],
                    stdout=fout, stderr=stderr)
        finally:
            fout.close()
            stderr.close()


def parsePAF(paf_file):
    """Parse ``*.paf`` file as created by ``minimap2``.

    The ``*.paf`` file is assumed to be created with
    the ``minimap2`` options ``-c --cs=long``, which
    creates long `cs` tags with the CIGAR string.

    PAF format is `described here <https://github.com/lh3/miniasm/blob/master/PAF.md>`_.
    PAF long CIGAR string format from ``minimap2`` is 
    `described here <https://github.com/lh3/minimap2>`_.

    Args:
        `paf_file` (str)
            Name of ``*.paf`` file.

    Returns:
        A generator that yields alignments on each line in 
        `paf_file`. Each alignment is returned as the 2-tuple
        `(query_name, a)`, where `query_name` is a str giving
        the name of the query sequence, and `a` is a named tuple
        with the following attributes:
            `ctg`: name of reference to which query is mapped
            `r_st`, `r_en`: start / end in reference (0 based)
            `q_st`, `q_en`: start / end in query (0 based)
            `strand`: 1 for forward, -1 for reverse
            `mapq`: mapping quality
            `cigar_str`: CIGAR string
    """
    assert os.path.isfile(paf_file), "no `paf_file` {0}".format(paf_file)
    
    Alignment = collections.namedtuple('Alignment', ['ctg', 'r_st',
            'r_en', 'q_st', 'q_en', 'strand', 'cigar_str'])
    cigar_m = re.compile('cs:Z:(?P<cigar_str>'
            '[:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+]+)(?:\s+|$)')

    with open(paf_file) as f:
        for line in f:
            entries = line.split('\t', maxsplit=12)
            try:
                cigar_str = cigar_m.search(entries[12]).group('cigar_str')
            except:
                raise ValueError("Cannot match CIGAR:\n{0}".format(entries[12]))
            query_name = entries[0]
            a = Alignment(ctg=entries[5],
                          r_st=int(entries[7]),
                          r_en=int(entries[8]),
                          q_st=int(entries[2]),
                          q_en=int(entries[3]),
                          strand={'+':1, '-':-1}[entries[4]],
                          cigar_str=cigar_str)
            yield (query_name, a)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
