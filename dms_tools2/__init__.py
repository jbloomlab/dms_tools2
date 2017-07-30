"""
======================
dms_tools2
======================

Python package for analyzing deep mutational scanning (DMS) data.

See http://jbloomlab.github.io/dms_tools2 for documentation.

Package-level constants
-------------------------
The following constants are defined for the package:

`AAS` : alphabetized list of all 20 amino acids

`AAS_WITHSTOP` : like `AAS` but includes stop codon (`*`) at end

`NTS` : alphabetized list of all nucleotides

`CODONS` : alphabetized list of all codons

`CODON_TO_AA` : dict translating codons to amino acids

`AA_TO_CODONS` : dict back-translating amino acid to list of codons.

"""

# import package-level metadata
from ._metadata import __version__
from ._metadata import __author__
from ._metadata import __author_email__
from ._metadata import __url__

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
for var in ['aa', 'nt', 'nt1', 'nt2', 'nt3', 'codon']:
    if var in locals():
        del locals()[var]
