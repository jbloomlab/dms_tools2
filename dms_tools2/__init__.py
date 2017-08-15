"""
======================
dms_tools2
======================

Python package for analyzing deep mutational scanning (DMS) data.

See http://jbloomlab.github.io/dms_tools2 for documentation.
"""

# import package-level metadata
from ._metadata import __version__
from ._metadata import __author__
from ._metadata import __author_email__
from ._metadata import __url__

# define constants related to nucleotides / amino acids / codons
import Bio.Alphabet.IUPAC
import Bio.Seq

#: alphabetized list of all 20 amino acids
AAS = sorted([_aa.upper() for _aa in 
        Bio.Alphabet.IUPAC.IUPACProtein.letters])

#: like `AAS` but includes stop codon (``*``) at end
AAS_WITHSTOP = AAS + ['*']

#: alphabetized list of all nucleotides
NTS = sorted([_nt.upper() for _nt in 
        Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters])

#: dict mapping each nucleotide to its complement.
NTCOMPLEMENT = dict([(_nt, str(Bio.Seq.Seq(_nt).reverse_complement()))
        for _nt in NTS] + [('N', 'N')])

#: alphabetized list of all codons
CODONS = ['{0}{1}{2}'.format(_nt1, _nt2, _nt3) for _nt1 in NTS 
        for _nt2 in NTS for _nt3 in NTS]

#: dict translating codons to amino acids
CODON_TO_AA = dict([(_codon, str(Bio.Seq.Seq(_codon).translate())) 
        for _codon in CODONS])

#: dict back-translating amino acid to list of codons.
AA_TO_CODONS = {}
for _aa in AAS_WITHSTOP:
    AA_TO_CODONS[_aa] = [_codon for _codon in CODONS if 
            CODON_TO_AA[_codon] == _aa]

# following lines needed because list comprehension variables remain
# in Python2 but not Python3, and we want to be compatible with both
for var in ['_aa', '_nt', '_nt1', '_nt2', '_nt3', '_codon']:
    if var in locals():
        del locals()[var]
