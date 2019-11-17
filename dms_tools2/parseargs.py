"""
===================
parseargs
===================

Argument parsing for the executable scripts in ``dms_tools2``.
"""


import re
import argparse
import dms_tools2


def checkName(name, nametype):
    """Check if `name` is an allowable name.

    Allowed names can contain most characters but **not**
    LaTex special characters.

    Args:
        `name` (str)
            Name to check
        `nametype` (str)
            If we print an exception what do we call the parameter
            that failed? For instance, `name` or `group`.

    Returns:
        Returns `True` if `name` is an allowable name.
        Otherwise raises a `ValueError` explaining why the
        `name` is invalid.

    >>> checkName('sample-1', 'name')
    True

    >>> checkName('sample 1', 'name')
    True

    >>> checkName('PGT151 - 5 nM', 'name')
    True

    >>> checkName('sample_1', 'name')
    Traceback (most recent call last):
        ...
    ValueError: name sample_1 contains following illegal characters: _

    >>> checkName('sample_1', 'group')
    Traceback (most recent call last):
        ...
    ValueError: group sample_1 contains following illegal characters: _
    """
    if not name or name.isspace():
        raise ValueError("{0} is all whitespace".format(nametype))
    illegal_chars = [c for c in name if 
            re.search(r'^[a-zA-Z0-9\- \.]$', c) is None]
    if illegal_chars:
        raise ValueError("{0} {1} contains following illegal characters: "
                "{2}".format(nametype, name, ', '.join(illegal_chars)))
    return True

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
            help='Length of NNN... barcode at start of each read. '
                 'Assumed to be same for R1 and R2, use `--bclen2` '
                 'if this is not the case.')

    parser.add_argument('--fastqdir',
            help='R1 and R2 files in this directory.')

    parser.add_argument('--R2', nargs='+', help=("Read 2 (R2) FASTQ "
            "files assumed to have same names as R1 but with "
            "'_R1' replaced by '_R2'. If that is not case, provide "
            "names here."))

    parser.add_argument('--R1trim', type=int, nargs='+',
        help=("Trim R1 from 3' end to this length. One value for all "
        "reads or values for each subamplicon in ``--alignspecs``."))

    parser.add_argument('--R2trim', type=int, nargs='+',
        help="Like '--R1trim', but for R2.")

    parser.add_argument('--bclen2', type=int, help='If R1 and R2 have '
            'different length barcodes, use `--bclen` for R1 length '
            'and `--bclen2` for R2 length.')

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

    parser.add_argument('--sitemask', help='Use to only consider '
            'mutations at a subset of sites. Should be a CSV file '
            'with column named `site` listing all sites to include.')

    parser.add_argument('--purgeread', type=float, default=0,
            help=("Randomly purge read pairs with this probability "
            "to subsample data."))

    parser.add_argument('--purgebc', type=float, default=0, 
            help=("Randomly purge barcodes with this probability to "
            "subsample data."))

    parser.set_defaults(bcinfo=False, bcinfo_csv=False)
    parser.add_argument('--bcinfo', dest='bcinfo', action='store_true',
            help=("Create file with suffix 'bcinfo.txt.gz' with info "
            "about each barcode."))
    parser.add_argument('--bcinfo_csv', dest='bcinfo_csv',
            action='store_true', help=("Store 'bcinfo' file as a csv "
            "with the suffix 'bcinfo.csv.gz'. Only has an effect if "
            "`--bcinfo` is used."))

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
            "`name`, `R1`. Can optionally have columns `R1trim` and "
            "`R2trim` with spaces delimiting subamplicon-specific trimming. "
            "If `R1trim` / `R2trim` in batch file, do **not** "
            "also give values for ``--R1trim`` and ``--R2trim``. "
            "Other columns are ignored, so other "
            "``dms2_bcsubamp`` args should be passed as separate "
            "command line args rather than in ``--batchfile``.", 
            required=True)

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
            "`err` or `errpre` and `errpost`. Other columns are ignored, "
            "so other ``dms2_prefs`` args should be passed as separate "
            "command line args rather than in ``--batchfile``.",
            required=True)

    parser.add_argument('--summaryprefix', required=True,
            help="Prefix of output summary files and plots.")

    parser.set_defaults(no_corr=False)
    parser.add_argument('--no_corr', dest='no_corr', action='store_true',
            help='Do not create correlation plot.')

    parser.set_defaults(no_avg=False)
    parser.add_argument('--no_avg', dest='no_avg', action='store_true',
            help='Do not create average prefs CSV.')

    return parser


def diffselParentParser():
    """Parent parser for ``dms2_diffsel`` / ``dms2_batch_diffsel``."""
    parser = argparse.ArgumentParser(
            parents=[parentParser()],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--indir', help="Input counts files in this "
            "directory.")

    parser.add_argument('--chartype', default='codon_to_aa',
            choices=['codon_to_aa'], help="Characters for which "
            "differential selection is estimated. `codon_to_aa` = amino "
            "acids from codon counts.")

    parser.add_argument('--excludestop', default='yes', choices=['yes', 'no'],
            help="Exclude stop codons as a possible amino acid?")

    parser.add_argument('--pseudocount', default=5, type=float,
            help="Pseudocount added to each count for sample with smaller "
            "depth; pseudocount for other sample scaled by relative depth.")

    parser.add_argument('--mincount', default=0, type=float,
            help="Report as `NaN` the diffsel of mutations for which both "
            "selected and mock-selected samples have < this many counts.")

    return parser


def diffselParser():
    """Returns `argparse.ArgumentParser` for ``dms2_diffsel``."""
    parser = argparse.ArgumentParser(
            description=parserDescription(
                'Estimate differential selection from mutation counts.'),
            parents=[diffselParentParser()],
            conflict_handler='resolve',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--name', required=True, 
            help='Name used for output files.')

    parser.add_argument('--sel', required=True, help="Post-selection "
            "counts file or prefix used when creating this file.")

    parser.add_argument('--mock', required=True, help="Like ``--sel``, "
            "but for mock-selection counts.")

    parser.add_argument('--err', help="Like ``--sel`` but for "
            "error-control to correct mutation counts.")

    return parser


def batch_diffselParser():
    """Returns `argparse.ArgumentParser` for ``dms2_batch_diffsel``."""
    parser = argparse.ArgumentParser(
            description=parserDescription(
                'Perform many runs of ``dms2_diffsel`` and summarize results.'),
            parents=[diffselParentParser()],
            conflict_handler='resolve',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--batchfile', required=True, help='CSV file '
            'specifying each ``dms2_diffsel`` run. Must have these '
            'columns: `name`, `sel`, `mock`. Can also have these: '
            '`err`, `group`, `grouplabel`. If `group` is used, '
            'samples are grouped in summary plots labeled by `group`, '
            'or by `grouplabel` if specified. Other columns are ignored, '
            'so other ``dms2_diffsel`` args should be passed as separate '
            'command line args rather than in ``--batchfile``.')

    parser.add_argument('--summaryprefix', required=True,
            help='Prefix of output summary files and plots.')

    return parser


def fracsurviveParser():
    """Returns `argparse.ArgumentParser` for ``dms2_fracsurvive``."""
    parser = argparse.ArgumentParser(
            description=parserDescription(
                'Estimate fraction surviving for each mutation.'),
            parents=[fracsurviveParentParser()],
            conflict_handler='resolve',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--name', required=True, 
            help='Name used for output files.')

    parser.add_argument('--sel', required=True, help="Post-selection "
            "counts file or prefix used when creating this file.")

    parser.add_argument('--mock', required=True, help="Like ``--sel``, "
            "but for mock-selection counts.")

    parser.add_argument('--libfracsurvive', required=True, type=float,
            help='Overall fraction of total library surviving selection '
            'versus mock condition. Should be between 0 and 1.')

    parser.add_argument('--err', help="Like ``--sel`` but for "
            "error-control to correct mutation counts.")

    return parser


def batch_fracsurviveParser():
    """Returns `argparse.ArgumentParser` for ``dms2_batch_fracsurvive``"""
    parser = argparse.ArgumentParser(
            description=parserDescription(
                'Perform runs of ``dms2_fracsurvive`` and summarize results.'),
            parents=[fracsurviveParentParser()],
            conflict_handler='resolve',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--batchfile', required=True, help='CSV file '
            'specifying each ``dms2_fracsurvive`` run. Must have these '
            'columns: `name`, `sel`, `mock`, `libfracsurvive`. Can also ' 
            'have these `err`, `group`, `grouplabel`. If `group` is used, '
            'samples are grouped in summary plots labeled by `group`, or by '
            '`grouplabel` if it is specified. Other columns are ignored, so '
            'other ``dms2_fracsurvive`` args should be passed as separate '
            'command line args rather than in ``--batchfile``.')

    parser.add_argument('--summaryprefix', required=True,
            help='Prefix of output summary files and plots.')

    return parser


def fracsurviveParentParser():
    """Parent parser ``dms2_fracsurvive`` / ``dms2_batch_fracsurvive``"""
    parser = argparse.ArgumentParser(
            parents=[parentParser()],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--indir', help="Input counts files in this "
            "directory.")

    parser.add_argument('--chartype', default='codon_to_aa',
            choices=['codon_to_aa'], help="Characters for which "
            "fraction surviving selection is estimated. `codon_to_aa` ="
            " amino acids from codon counts.")

    parser.add_argument('--aboveavg', default='no', choices=['yes', 'no'],
            help="Report fracsurvive **above** the library average "
            "rather than direct fracsurvive values.")

    parser.add_argument('--excludestop', default='yes', choices=['yes', 'no'],
            help="Exclude stop codons as a possible amino acid?")

    parser.add_argument('--pseudocount', default=5, type=float,
            help="Pseudocount added to each count for sample with smaller "
            "depth; pseudocount for other sample scaled by relative depth.")

    parser.add_argument('--mincount', default=0, type=float,
            help="Report as `NaN` the fracsurvive of mutations for which "
            "both selected and mock-selected samples have < this many counts.")

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
    group.add_argument('--fracsurvive', help='CSV file of amino-acid '
            'fraction surviving.')
    group.add_argument('--diffprefs', help='CSV file of differences in '
            'amino-acid preferences.')
    group.add_argument('--muteffects', help='CSV file of amino-acid '
            'mutational effects.')

    parser.add_argument('--name', required=True, 
            help='Name used for output files.')

    parser.add_argument('--nperline', help='Number of sites per line.',
            type=int, default=70)

    parser.add_argument('--numberevery', type=int, default=10,
            help='Number sites at this interval.')

    parser.add_argument('--excludestop', choices=['yes', 'no'],
            default='no', help='Exclude stop codons as possible amino '
            'acid?')

    parser.add_argument('--stringency', type=float, default=1,
            help='Stringency parameter to re-scale prefs.')

    parser.add_argument('--restrictdiffsel', default='all',
            choices=['all', 'positive', 'negative'], 
            help='Plot all diffsel, or only positive or negative.')

    parser.add_argument('--diffselrange', nargs=2, type=float,
            metavar=('MINDIFFSEL', 'MAXDIFFSEL'), 
            help='Specify a fixed range for `diffsel`. Otherwise '
            'determined from data range.')

    parser.add_argument('--muteffectrange', nargs=2, type=float,
            metavar=('MINMUTEFFECT', 'MAXMUTEFFECT'),
            help='Specify a fixed range for `muteffects`. Otherwise '
            'determined from data range.')

    parser.add_argument('--fracsurvivemax', type=float,
            help='Specify maximum value for `fracsurvive`. '
            'Otherwise determined from data range.')

    parser.add_argument('--sortsites', choices=['yes', 'no'],
            default='yes', help='Sort sites from first to last '
            'before plotting.')

    parser.add_argument('--mapmetric', default='functionalgroup', choices=[
            'kd', 'mw', 'charge', 'functionalgroup', 'singlecolor'],
            help='Color amino acids '
            'by Kyte-Doolittle hydrophobicity, molecular weight, charge, '
            'or functional group.')

    parser.add_argument('--colormap', default='jet', help="`matplotlib "
            "color map <http://matplotlib.org/users/colormaps.html>`_ for"
            " amino acids when `--mapmetric` is 'kd' or 'mw'; name of "
            "single color when it is 'singlecolor'.")

    parser.add_argument('--overlay1', nargs=3,
            metavar=('FILE', 'SHORTNAME', 'LONGNAME'),
            help="Color bar above logo plot to denote per-residue "
            "property. FILE is CSV format with column names `site` "
            "and SHORTNAME. SHORTNAME is <= 5 character property name. "
            "LONGNAME is longer name for legend. Sites not in FILE "
            "are colored white. To show wildtype identity, make "
            "SHORTNAME and LONGNAME both `wildtype` and have this "
            "column in FILE give 1-letter wildtype amino-acid code. "
            "To show ``omegabysite.txt`` file from ``phydms``, "
            "give that file and set both SHORTNAME and LONGNAME "
            "to `omegabysite`.")

    parser.add_argument('--overlay2', default=None, nargs=3,
            metavar=('FILE', 'SHORTNAME', 'LONGNAME'),
            help='Second overlay color bar.')

    parser.add_argument('--overlay3', default=None, nargs=3,
            metavar=('FILE', 'SHORTNAME', 'LONGNAME'),
            help='Third overlay color bar.')

    parser.add_argument('--underlay', default='no', choices=['yes', 'no'],
            help='Plot underlay rather than overlay bars.')

    parser.add_argument('--scalebar', default=None, nargs=2,
            metavar=('BARHEIGHT', 'LABEL'), help='Plot a scale bar '
            'indicating BARHEIGHT with LABEL. Only for `diffsel`, '
            '`fracsurvive`, and `muteffects`.')

    parser.add_argument('--overlaycolormap', default='jet', help="`matplotlib "
            "color map <http://matplotlib.org/users/colormaps.html>`_ for"
            " overlay bars (e.g., 'jet' or 'YlOrRd').")

    parser.add_argument('--letterheight', type=float, default=1,
            help="Relative height of letter stacks in logo plot.")

    parser.add_argument('--ignore_extracols', default='no',
            choices=['yes', 'no'], help='Ignore extra columns in data')

    parser.add_argument('--sepline', default='yes', choices=['yes', 'no'],
            help="Separate positive and negative diffsel with black line?")

    return parser



if __name__ == '__main__':
    import doctest
    doctest.testmod()
