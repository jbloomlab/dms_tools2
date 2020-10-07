"""
======================
dms_tools2
======================

Python package for analyzing deep mutational scanning (DMS) data.

See http://jbloomlab.github.io/dms_tools2 for documentation.
"""

import os
import glob


# import package-level metadata
from ._metadata import __version__
from ._metadata import __author__
from ._metadata import __author_email__
from ._metadata import __url__

# define constants related to nucleotides / amino acids / codons
import Bio.Seq

#: alphabetized list of all 20 amino acids
AAS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
       'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

#: like `AAS` but includes stop codon (``*``) at end
AAS_WITHSTOP = AAS + ['*']

#: alphabetized list of all nucleotides
NTS = ['A', 'C', 'G', 'T']

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
AA_TO_CODONS = {
        _aa:[_c for _c in CODONS if CODON_TO_AA[_c] == _aa]
        for _aa in AAS_WITHSTOP}

#: dict keyed by nucleotide code, values `re` matches for code.
NT_TO_REGEXP = dict(
        map(lambda tup: (tup[0], tup[1]) if len(tup[1]) == 1 else
                        (tup[0], '[' + ''.join(sorted(tup[1])) + ']'), 
        Bio.Data.IUPACData.ambiguous_dna_values.items()))

# import all modules as here: https://stackoverflow.com/a/1057534
# but do not import `rplot` as it is optional.
__all__ = [os.path.basename(f)[ : -3] for f in
           glob.glob(os.path.dirname(__file__) + '/*.py') if
           os.path.isfile(f) and not os.path.basename(f).startswith('_')
           and not f.endswith('rplot.py')]
from . import *
