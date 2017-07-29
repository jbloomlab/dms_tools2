"""
===========================
``dms_tools2`` package
===========================
This package is for analyzing deep mutational scanning (DMS) data.

See http://jbloomlab.github.io/dms_tools2 for documentation.

"""

# import package-level metadata
from ._metadata import __version__
from ._metadata import __author__
from ._metadata import __author_email__
from ._metadata import __url__
from ._metadata import __acknowledgments__

# import all modules as here
import os.path 
import glob
__all__ = [os.path.basename(f)[ : -3] for f in 
        glob.glob(os.path.join(os.path.dirname(__file__), "*.py"))
        if os.path.isfile(f) and not f.endswith('__init__.py')]


# define constants related to nucleotides / amino acids / codons
import Bio.Alphabet.IUPAC
import Bio.Seq

AAS = sorted([aa.upper() for aa in 
        Bio.Alphabet.IUPAC.IUPACProtein.letters])

AAS_WITHSTOP = AAS + ['*']

NTS = sorted([nt.upper() for nt in 
        Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters])

CODONS = ['{0}{1}{2}'.format(nt1, nt2, nt3) for nt1 in NTS 
        for nt2 in NTS for nt3 in NTS]

CODON_TO_AA = dict([(codon, str(Bio.Seq.Seq(codon).translate())) 
        for codon in CODONS])

AA_TO_CODONS = {}
for aa in AAS_WITHSTOP:
    AA_TO_CODONS[aa] = [codon for codon in CODONS if 
            CODON_TO_AA[codon] == aa]

# following lines needed because list comprehension variables remain
# in Python2 but not Python3, and we want to be compatible with both
for var in ['f', 'aa', 'nt', 'nt1', 'nt2', 'nt3', 'codon']:
    del locals()[var]
