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

            - ``--outdir``

            - ``--ncpus``

            - ``--use_existing``

            - ``-v / --version``
    """
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument('--outdir',  
            help='Output files to this directory (create if needed).')

    parser.add_argument('--ncpus', type=int, default=-1, 
            help="Number of CPUs to use, -1 is all available.")

    parser.add_argument('--use_existing', choices=['yes', 'no'],
            default='no', help=('If files with names of expected '
            'output already exist, do not re-run.'))

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
    return ("{0} Part of `{1} <{4}>`_ (version {2}) written by {3}."
            .format(description, dms_tools2.__name__,
            dms_tools2.__version__, dms_tools2.__author__,
            dms_tools2.__url__))


def bcsubampParentParser():
    """Parent parser for ``dms2_bcsubamp`` / ``dms2_batch_bcsubamp``."""
    parser = argparse.ArgumentParser(
            parents=[parentParser()],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--refseq', required=True, 
            help='Align subamplicons to gene in this FASTA file.')

    parser.add_argument('--alignspecs', required=True, nargs='+',
            help=("Subamplicon alignment positions as "
            "'REFSEQSTART,REFSEQEND,R1START,R2START'. "
            "REFSEQSTART is nt (1, 2, ... numbering) in "
            "'refseq' where nt R1START in R1 aligns. "
            "REFSEQEND is nt in 'refseq' where nt R2START "
            "in R2 aligns.'"))

    parser.add_argument('--bclen', type=int, default=8,
            help='Length of NNN... barcode at start of each read.')

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

    parser.add_argument('--minfraccall', type=float, default=0.95, 
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

    return parser


def bcsubampParser():
    """Returns `argparse.ArgumentParser` for ``dms2_bcsubamp``."""
    parser = argparse.ArgumentParser(
            description=parserDescription(
                'Align barcoded subamplicons and count mutations.'),
            parents=[bcsubampParentParser()],
            conflict_handler='resolve',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--name', required=True, 
            help='Sample name used for output files.')

    parser.add_argument('--R1', required=True, nargs='+',
            help="Read 1 (R1) FASTQ files, can be gzipped. "
            "See also '--fastqdir'.")

    return parser


def batch_bcsubampParser():
    """Returns `argparse.ArgumentParser` for ``dms2_batch_bcsubamp``."""
    parser = argparse.ArgumentParser(
            description=parserDescription(
                'Perform many runs of ``dms2_bcsubamp`` and plot results.'),
            parents=[bcsubampParentParser()],
            conflict_handler='resolve',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--batchfile', help="CSV file specifying each "
            "``dms2_bcsubamp`` run. Must have these columns: "
            "`name`, `R1`", required=True)

    parser.add_argument('--summaryprefix', required=True,
            help="Prefix of output summary plots.")

    return parser


def prefsParentParser():
    """Parent parser for ``dms2_prefs`` / ``dms2_batch_prefs``."""
    parser = argparse.ArgumentParser(
            parents=[parentParser()],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--method', default='bayesian', 
            choices=['ratio', 'bayesian'], help="Method to "
            "estimate preferences: normalized enrichment ratios "
            "or Bayesian inference.")

    parser.add_argument('--indir', help="Input counts files in this "
            "directory.")

    parser.add_argument('--chartype', default='codon_to_aa',
            choices=['codon_to_aa'], help="Characters for which "
            "preferences are estimated. `codon_to_aa` = amino acids "
            "from codon counts.")

    parser.add_argument('--excludestop', default='yes', choices=['yes', 'no'],
            help="Exclude stop codons as a possible amino acid?")

    parser.add_argument('--conc', nargs=3, default=[1, 1, 1],
            type=float, metavar=('Cprefs', 'Cmut', 'Cerr'),
            help="Concentration parameters for priors for "
            "``--method bayesian``. Priors are over preferences, "
            "mutagenesis rate, and error rate(s).")

    parser.add_argument('--pseudocount', default=1,
            help="Pseudocount used with ``--method ratio``.")

    return parser


def prefsParser():
    """Returns `argparse.ArgumentParser` for ``dms2_prefs``."""
    parser = argparse.ArgumentParser(
            description=parserDescription(
                'Estimate preferences from mutation counts.'),
            parents=[prefsParentParser()],
            conflict_handler='resolve',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--pre', required=True, 
            help='Pre-selection counts file or prefix used when creating '
            'this file.')

    parser.add_argument('--post', required=True, 
            help='Like ``--pre`` but for post-selection counts.')

    parser.add_argument('--name', required=True, 
            help='Name used for output files.')

    parser.add_argument('--err', nargs=2, metavar=('ERRPRE', 'ERRPOST'),
            help="Like ``--pre`` but for counts for error control(s) for "
            "``--pre`` and ``--post``. Specify same file twice for same "
            "control for both.")

    return parser


def batch_prefsParser():
    """Returns `argparse.ArgumentParser` for ``dms2_batch_prefs``."""
    parser = argparse.ArgumentParser(
            description=parserDescription(
                'Perform many runs of ``dms2_prefs`` and summarize results.'),
            parents=[prefsParentParser()],
            conflict_handler='resolve',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--batchfile', help="CSV file specifying each "
            "``dms2_prefs`` run. Must have these columns: "
            "`name`, `pre`, `post`. Can also have these columns: "
            "`err` or `errpre` and `errpost`.", required=True)

    parser.add_argument('--summaryprefix', required=True,
            help="Prefix of output summary files and plots.")

    return parser


def logoplotParser():
    """Returns `argparse.ArgumentParser` for ``dms2_logoplot``."""
    parser = argparse.ArgumentParser(
            description=parserDescription(
                'Create logo plot visualization.'),
            parents=[parentParser()],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prefs', help='CSV file of amino-acid preferences.')
    group.add_argument('--diffsel', help='CSV file of amino-acid '
            'differential selection.')

    parser.add_argument('--name', required=True, 
            help='Name used for output files.')

    parser.add_argument('--nperline', help='Number of sites per line.',
            type=int, default=70)

    parser.add_argument('--numberevery', type=int, default=10,
            help='Number sites at this interval.')

    parser.add_argument('--excludestop', choices=['yes', 'no'],
            default='no', help='Exclude stop codons a possible amino '
            'acid?')

    parser.add_argument('--stringency', type=float, default=1,
            help='Stringency parameter to re-scale prefs.')

    parser.add_argument('--sortsites', choices=['yes', 'no'],
            default='yes', help='Sort sites from first to last '
            'before plotting.')

    parser.add_argument('--mapmetric', default='functionalgroup', choices=[
            'kd', 'mw', 'charge', 'functionalgroup'], help='Color amino acids '
            'by Kyte-Doolittle hydrophobicity, molecular weight, charge, '
            'or functional group.')

    parser.add_argument('--colormap', default='jet', help="`matplotlib "
            "color map <http://matplotlib.org/users/colormaps.html>`_ for"
            " amino acids when `--mapmetric` is 'kd' or 'mw'.")

    parser.add_argument('--overlaycolormap', default='jet', help="`matplotlib "
            "color map <http://matplotlib.org/users/colormaps.html>`_ for"
            " overlay bars (also consider 'YlOrRd').")

    parser.add_argument('--letterheight', type=int, default=1,
            help="Relative height of letter stacks in logo plot.")

    parser.add_argument('--sepline', default='yes', choices=['yes', 'no'],
            help="Separate positive and negative diffsel with black line?")

    return parser



if __name__ == '__main__':
    import doctest
    doctest.testmod()
