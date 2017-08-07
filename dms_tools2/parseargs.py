"""
===================
parseargs
===================

Argument parsing for the executable scripts in ``dms_tools2``.
"""


import argparse
import dms_tools2


def parentParser():
    """Returns parent parser with some common options added.

    Returns:
        `argparse.ArgumentParser` with the following arguments
        already added: 

            - ``--use_existing``

            - ``-v / --version``
    """
    parser = argparse.ArgumentParser(add_help=False)

    parser.set_defaults(use_existing=False)
    parser.add_argument('--use_existing', dest='use_existing',
            action='store_true', help=('Do nothing if files with '
            'names of all expected output already exist.'))

    parser.add_argument('-v', '--version', action='version', 
            version='%(prog)s {0}'.format(dms_tools2.__version__))

    return parser


def parserDescription(description):
    """Augments description with program information.

    Args:
        `description` (str)
            Description of program

    Returns:
        A string with `description` augmented with information
        on the `dms_tools2` package / version.
    """
    return ("{0} Part of {1} (version {2}) written by {3}. "
            "See {4} for documentation.".format(
            description, dms_tools2.__name__,
            dms_tools2.__version__, dms_tools2.__author__,
            dms_tools2.__url__))


def bcsubampParser():
    """Returns `argparse.ArgumentParser` for ``dms2_bcsubamp``."""
    parser = argparse.ArgumentParser(
            description=parserDescription(
                'Align barcoded subamplicons and count mutations.'),
            parents=[parentParser()],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--name', required=True, 
            help='Sample name used for output files.')

    parser.add_argument('--refseq', required=True, 
            help='Align subamplicons to gene in this FASTA file.')

    parser.add_argument('--alignspecs', required=True, nargs='+',
            help=("Subamplicon alignment positions as "
            "'REFSEQSTART,REFSEQEND,R1START,R2START'. "
            "REFSEQSTART is nt (1, 2, ... numbering) in "
            "'refseq' where nt R1START in R1 aligns. "
            "REFSEQEND is nt in 'refseq' where nt R2START "
            "in R2 aligns.'"))

    parser.add_argument('--R1', required=True, nargs='+',
            help="Read 1 (R1) FASTQ files, can be gzipped. "
            "See also '--fastqdir'.")

    parser.add_argument('--bclen', type=int, default=8,
            help='Length of NNN... barcode at start of each read.')

    parser.add_argument('--outdir',  
            help='Output files to this directory (create if needed).')

    parser.add_argument('--fastqdir',
            help='R1 and R2 files in this directory.')

    parser.add_argument('--R2', nargs='+', help=("Read 2 (R2) FASTQ "
            "files assumed to have same names as R1 but with "
            "'_R1' replaced by '_R2'. If that is not case, provide "
            "names here."))

    parser.add_argument('--R1trim', type=int, default=300, nargs='+',
        help=("Trim R1 from 3' end to this length. One value for "
        "all reads or values for each subamplicon in 'alignspecs'."))

    parser.add_argument('--R2trim', type=int, default=300, nargs='+',
        help="Like '--R1trim', but for R2.")

    parser.add_argument('--chartype', default='codon', choices=['codon'],
            help='Character type for which we count mutations.')

    parser.add_argument('--maxmuts', type=int, default=4, 
            help=("Max allowed mismatches in alignment of subamplicon; "
            "mismatches counted in terms of character '--chartype'."))

    parser.add_argument('--minq', type=int, default=15,
            help="Only call nucleotides with Q score >= this.")

    parser.add_argument ('--minreads', type=int, default=2, 
            help=("Require this many reads in a barcode to agree to "
            "call consensus nucleotide identity."))

    parser.add_argument('--minfraccall', type=float, default=0.9, 
            help=("Retain only barcodes where trimmed consensus "
            "sequence for each read has >= this frac sites called."))

    parser.add_argument('--minconcur', default=0.75, type=float,
            help=('Only call consensus identity for barcode when >= '
            'this fraction of reads concur.'))

    parser.add_argument('--purgeread', type=float, default=0,
            help=("Randomly purge read pairs with this probability "
            "to subsample data."))

    parser.add_argument('--purgebc', type=float, default=0, 
            help=("Randomly purge barcodes with this probability to "
            "subsample data."))

    parser.set_defaults(bcinfo=False)
    parser.add_argument('--bcinfo', dest='bcinfo', action='store_true', 
            help=("Create file with suffix 'bcinfo.txt.gz' with info "
            "about each barcode."))

    parser.add_argument('--seed', type=int, default=1, 
            help="Random number seed for '--purgeread' and '--purgebc'.")
    
    return parser



if __name__ == '__main__':
    import doctest
    doctest.testmod()
