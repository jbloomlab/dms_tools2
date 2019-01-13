"""
=================
codonvarianttable
=================

This module defines :class:`CodonVariantTable` objects
for storing and handling codon variants of a gene. They are
useful when you are performing deep mutational scanning
on barcoded codon-mutant libraries of a gene.

.. contents:: Table of Contents
   :local:
   :depth: 2

Example on simulated data
--------------------------

.. code-links::

Here is an example on simulated data:

Import Python modules / packages. For ggplot2-like plotting,
we use `plotnine <https://plotnine.readthedocs.io>`_.

.. nbplot::

   >>> import collections
   >>> import os
   >>> import math
   >>> import tempfile
   >>> import scipy
   >>> import pandas
   >>> from plotnine import *
   >>> from dms_tools2 import CODONS, CODON_TO_AA, NTS, AAS_WITHSTOP
   >>> from dms_tools2.plot import COLOR_BLIND_PALETTE_GRAY as PALETTE
   >>> NONSTOP_CODONS = [c for c in CODONS if CODON_TO_AA[c] != '*']
   >>> from dms_tools2.codonvarianttable import CodonVariantTable

Seed random number generator for reproducible output:

.. nbplot::

   >>> scipy.random.seed(1)

Set pandas display options to show large chunks of Data Frames
in this example:

.. nbplot::

    >>> pandas.set_option('display.max_columns', 10)
    >>> pandas.set_option('display.width', 500)

Initialize a CodonVariantTable
++++++++++++++++++++++++++++++++++++++++++++

Define a wildtype gene sequence of 100 codons. Our variants
will all be of this sequence.

.. nbplot::

    >>> genelength = 100
    >>> geneseq = ''.join(scipy.random.choice(NONSTOP_CODONS)
    ...                   for _ in range(genelength))
    >>> geneseq
    'GCCGGTATAAGAAGCAGTACCATTAAACAAAACATAACTGTCTTTACGCGCTCACCAGCCCAGCCAAGTTTTGGGCTACTCATGTCAACACCTCCTGGCTATTGGCTGGAACCGATCGGCAGCACTCCGTTAAACAAATTTCACAGACGAATCTCCGTTGGGTTAAGACTGACTAATACGCCCTATTTAAATACACGATATGGTATACGGTCGCAAGTCTCCGGCCAGATTAAAACACGCGTTGAGCCTACTCGGCGCGGACCGAGCAATGCTCCTGCACGTGCCTTACATGCGAGAGAA'

To initialize a :class:`CodonVariantTable`, we need a CSV
file giving the barcodes and nucleotide substitutions in
the gene (in 1, 2, ... numbering) for each variant. You also
need to specify which library each variant belongs to (you will
have at least one library, but hopefully you performed experimental
replicates), as well as the variant call support. This is how
confident you are that the variant is called correctly, and might
be the number of PacBio CCSs that support that variant. Here we
simulate the data for two libraries and put that simulated data
into a pandas Data Frame.

We will simulate two libraries, each with 10,000 variants
variants, and with barcodes of length 16. We will have a Poisson
distributed number of codon mutations per variant with a mean of
1.3. We simulate a pandas Data Frame with such data. So the only
purpose the large code block below is to simulate the type of
data that we might get for a real library in order to initialize
a plausible :class:`CodonVariantTable`:

.. nbplot::

   >>> libs = ['lib_1', 'lib_2']
   >>> variants_per_lib = 10000
   >>> bclen = 16
   >>> barcode_variant_dict = collections.OrderedDict([
   ...     ('library', []), ('barcode', []), ('substitutions', []),
   ...     ('variant_call_support', [])])
   >>> for lib in libs:
   ...     for ivariant in range(variants_per_lib):
   ...         barcode = ''.join(scipy.random.choice(NTS)
   ...                           for _ in range(bclen))
   ...         while barcode in barcode_variant_dict['barcode']:
   ...             barcode = ''.join(scipy.random.choice(NTS)
   ...                               for _ in range(bclen))
   ...         substitutions = []
   ...         ncodonmuts = scipy.random.poisson(1.3)
   ...         for icodon in scipy.random.randint(1, genelength + 1, ncodonmuts):
   ...             wtcodon = geneseq[3 * (icodon - 1) : 3 * icodon]
   ...             mutcodon = scipy.random.choice([c for c in CODONS if c != wtcodon])
   ...             for i_nt, (wt_nt, mut_nt) in enumerate(zip(wtcodon, mutcodon)):
   ...                 if wt_nt != mut_nt:
   ...                     igene = 3 * (icodon - 1) + i_nt + 1 # nucleotide in gene
   ...                     substitutions.append(f'{wt_nt}{igene}{mut_nt}')
   ...         barcode_variant_dict['library'].append(lib)  
   ...         barcode_variant_dict['barcode'].append(barcode)  
   ...         barcode_variant_dict['substitutions'].append(' '.join(substitutions))
   ...         barcode_variant_dict['variant_call_support'].append(scipy.random.randint(1, 4))
   >>> barcode_variants = pandas.DataFrame(barcode_variant_dict)

Here are the first and last few lines of the Data Frame. As you can
see, it gives the nucleotide mutations (in 1, 2, ... numbering)
for each barcode:

.. nbplot::

   >>> barcode_variants.head(n=4)
     library           barcode                substitutions  variant_call_support
   0   lib_1  GGTTTTCTTAGGGAAT  A232T A233C A234C A19C C21G                     3
   1   lib_1  GATCTTCATATAGATT            C247G C248A T249C                     2
   2   lib_1  GTACTCAGAAAACGAC                                                  2
   3   lib_1  GAATGGGACGGGGTTA                        G129A                     3
   >>> barcode_variants.tail(n=4)
         library           barcode                                 substitutions  variant_call_support
   19996   lib_2  ATGTACTGTGATTTAT              G118C G119A C120G T43G T44A T45C                     1
   19997   lib_2  GCTCATTATGCCGGCC  G1A C2G C3A C226A A227G G228T A25C A26T A27G                     2
   19998   lib_2  GACAATCTGCAACCGT                         C263G G264C A37C T39C                     1
   19999   lib_2  TCATGTGCTGAGACGA                             A265G G266A C267A                     2

We now write this Data Frame to a CSV file and use it to initialize
a :class:`CodonVariantTable`. Obviously in a real experiment you
would have determined the CSV file giving your barcode-variant
identities experimentally, and would be passing that experimentally
determined data in CSV format:

.. nbplot::

   >>> with tempfile.NamedTemporaryFile(mode='w') as f:
   ...     barcode_variants.to_csv(f, index=False)
   ...     f.flush()
   ...     variants = CodonVariantTable(
   ...                 barcode_variant_file=f.name,
   ...                 geneseq=geneseq)

Now we have a :class:`CodonVariantTable` with the data on our barcodes
and variants. We can get basic information such as:

The codon sites:

.. nbplot::

   >>> variants.sites
   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]

The codon identities at each site:

.. nbplot::

   >>> variants.codons
   OrderedDict([(1, 'GCC'), (2, 'GGT'), (3, 'ATA'), (4, 'AGA'), (5, 'AGC'), (6, 'AGT'), (7, 'ACC'), (8, 'ATT'), (9, 'AAA'), (10, 'CAA'), (11, 'AAC'), (12, 'ATA'), (13, 'ACT'), (14, 'GTC'), (15, 'TTT'), (16, 'ACG'), (17, 'CGC'), (18, 'TCA'), (19, 'CCA'), (20, 'GCC'), (21, 'CAG'), (22, 'CCA'), (23, 'AGT'), (24, 'TTT'), (25, 'GGG'), (26, 'CTA'), (27, 'CTC'), (28, 'ATG'), (29, 'TCA'), (30, 'ACA'), (31, 'CCT'), (32, 'CCT'), (33, 'GGC'), (34, 'TAT'), (35, 'TGG'), (36, 'CTG'), (37, 'GAA'), (38, 'CCG'), (39, 'ATC'), (40, 'GGC'), (41, 'AGC'), (42, 'ACT'), (43, 'CCG'), (44, 'TTA'), (45, 'AAC'), (46, 'AAA'), (47, 'TTT'), (48, 'CAC'), (49, 'AGA'), (50, 'CGA'), (51, 'ATC'), (52, 'TCC'), (53, 'GTT'), (54, 'GGG'), (55, 'TTA'), (56, 'AGA'), (57, 'CTG'), (58, 'ACT'), (59, 'AAT'), (60, 'ACG'), (61, 'CCC'), (62, 'TAT'), (63, 'TTA'), (64, 'AAT'), (65, 'ACA'), (66, 'CGA'), (67, 'TAT'), (68, 'GGT'), (69, 'ATA'), (70, 'CGG'), (71, 'TCG'), (72, 'CAA'), (73, 'GTC'), (74, 'TCC'), (75, 'GGC'), (76, 'CAG'), (77, 'ATT'), (78, 'AAA'), (79, 'ACA'), (80, 'CGC'), (81, 'GTT'), (82, 'GAG'), (83, 'CCT'), (84, 'ACT'), (85, 'CGG'), (86, 'CGC'), (87, 'GGA'), (88, 'CCG'), (89, 'AGC'), (90, 'AAT'), (91, 'GCT'), (92, 'CCT'), (93, 'GCA'), (94, 'CGT'), (95, 'GCC'), (96, 'TTA'), (97, 'CAT'), (98, 'GCG'), (99, 'AGA'), (100, 'GAA')])

The amino-acid identities at each site:

.. nbplot::

   >>> variants.aas
   OrderedDict([(1, 'A'), (2, 'G'), (3, 'I'), (4, 'R'), (5, 'S'), (6, 'S'), (7, 'T'), (8, 'I'), (9, 'K'), (10, 'Q'), (11, 'N'), (12, 'I'), (13, 'T'), (14, 'V'), (15, 'F'), (16, 'T'), (17, 'R'), (18, 'S'), (19, 'P'), (20, 'A'), (21, 'Q'), (22, 'P'), (23, 'S'), (24, 'F'), (25, 'G'), (26, 'L'), (27, 'L'), (28, 'M'), (29, 'S'), (30, 'T'), (31, 'P'), (32, 'P'), (33, 'G'), (34, 'Y'), (35, 'W'), (36, 'L'), (37, 'E'), (38, 'P'), (39, 'I'), (40, 'G'), (41, 'S'), (42, 'T'), (43, 'P'), (44, 'L'), (45, 'N'), (46, 'K'), (47, 'F'), (48, 'H'), (49, 'R'), (50, 'R'), (51, 'I'), (52, 'S'), (53, 'V'), (54, 'G'), (55, 'L'), (56, 'R'), (57, 'L'), (58, 'T'), (59, 'N'), (60, 'T'), (61, 'P'), (62, 'Y'), (63, 'L'), (64, 'N'), (65, 'T'), (66, 'R'), (67, 'Y'), (68, 'G'), (69, 'I'), (70, 'R'), (71, 'S'), (72, 'Q'), (73, 'V'), (74, 'S'), (75, 'G'), (76, 'Q'), (77, 'I'), (78, 'K'), (79, 'T'), (80, 'R'), (81, 'V'), (82, 'E'), (83, 'P'), (84, 'T'), (85, 'R'), (86, 'R'), (87, 'G'), (88, 'P'), (89, 'S'), (90, 'N'), (91, 'A'), (92, 'P'), (93, 'A'), (94, 'R'), (95, 'A'), (96, 'L'), (97, 'H'), (98, 'A'), (99, 'R'), (100, 'E')])

The different libraries for which we have barcodes:

.. nbplot::

   >>> variants.libraries
   ['lib_1', 'lib_2']

We can also get a Data Frame that includes the information we
passed about the variants along with additional columns containing
amino-acid mutations and mutation counts via the
`barcode_variant_df` attribute of the :class:`CodonVariantTable`:

.. nbplot::

   >>> variants.barcode_variant_df.head(n=4)
     library           barcode  variant_call_support codon_substitutions aa_substitutions  n_codon_substitutions  n_aa_substitutions
   0   lib_1  AAAAAAACGTCTTGAC                     2            GTT81ACA             V81T                      1                   1
   1   lib_1  AAAAAAATGCTAGAAC                     2            ATC51CCC             I51P                      1                   1
   2   lib_1  AAAAACGTTTTAAGTC                     1            TAT67TAC                                       1                   0
   3   lib_1  AAAAACTAAAAGTTAA                     3                                                           0                   0

The :class:`CodonVariantTable` has a large number of methods
to get or plot information about the variants. Eventually,
we will add the counts of different variants for specific
samples using :meth:`CodonVariantTable.addSampleCounts`,
and then we can get information for those samples. But
for now, we will use some of the methods setting the
`samples` parameter to `None`. When we do that, we simply get
information about the properties of the variants that
are represented in the table.

:meth:`CodonVariantTable.n_variants_df` gives us
information about the number of variants:

.. nbplot::

   >>> variants.n_variants_df(samples=None)
            library             sample  count
   0          lib_1  barcoded variants  10000
   1          lib_2  barcoded variants  10000
   2  all libraries  barcoded variants  20000

You can use :meth:`CodonVariantTable.valid_barcodes`
to get the set of valid barcodes for each library. Here
we just show the first 10 for `lib_1`:

.. nbplot::

   >>> sorted(variants.valid_barcodes('lib_1'))[ : 10]
   ['AAAAAAACGTCTTGAC', 'AAAAAAATGCTAGAAC', 'AAAAACGTTTTAAGTC', 'AAAAACTAAAAGTTAA', 'AAAAAGAAGCTCACTC', 'AAAAAGAATTCTTGTA', 'AAAAAGCAGGCCAGAA', 'AAAAAGGGTGTAAGGA', 'AAAAATAAAAGAGTAA', 'AAAAATTAAGCCGATG']

Plot information about variants
++++++++++++++++++++++++++++++++++++++++++++

The :class:`CodonVariantTable` also has methods to plot
data about the variants. Since we have not yet added samples,
to the :class:`CodonVariantTable`, we call these methods
with `samples=None` to just get information about the
variants in the table.

Note that the plotting methods all return
`plotnine.ggplot <https://plotnine.readthedocs.io/en/stable/generated/plotnine.ggplot.html>`_
objects.

:meth:`CodonVariantTable.plotNumMutsHistogram` plots the
number of mutations per variant. Here we do that for
amino-acid mutations:

.. nbplot::

   >>> p = variants.plotNumMutsHistogram('aa', samples=None)
   >>> _ = p.draw()

Most of the plotting methods can also plot values
for codon mutations, such as here:

.. nbplot::

   >>> p = variants.plotNumMutsHistogram('codon', samples=None)
   >>> _ = p.draw()

:meth:`CodonVariantTable.plotNumCodonMutsByType` shows
the number of each type of codon mutation. It can be
called for just single mutants (and wildtype) variants:

.. nbplot::

   >>> p = variants.plotNumCodonMutsByType('single', samples=None)
   >>> _ = p.draw()

Or for variants with any number of mutations:

.. nbplot::

   >>> p = variants.plotNumCodonMutsByType('all', samples=None)
   >>> _ = p.draw()


:meth:`CodonVariantTable.plotCumulMutCoverage` shows the 
fraction of mutations seen <= some number of times:

.. nbplot::

   >>> p = variants.plotCumulMutCoverage('single', 'aa', samples=None)
   >>> _ = p.draw()

To get the numerical information about cumulative coverage
plotted above, use :meth:`CodonVariantTable.mutCounts`:

.. nbplot::

   >>> variants.mutCounts('single', 'aa', samples=None).head(n=4)
     library             sample mutation  count  mutation_type  site
   0   lib_1  barcoded variants     P88R     11  nonsynonymous    88
   1   lib_1  barcoded variants     E82S      9  nonsynonymous    82
   2   lib_1  barcoded variants     G40R      8  nonsynonymous    40
   3   lib_1  barcoded variants     K46R      8  nonsynonymous    46

:meth:`CodonVariantTable.plotMutFreqs` shows the frequency of
mutations along the gene:

.. nbplot::

   >>> p = variants.plotMutFreqs('single', 'aa', samples=None)
   >>> _ = p.draw()

:meth:`CodonVariantTable.plotMutHeatmap` shows the frequency of
mutations as a heat map:

.. nbplot::
 
   >>> p = variants.plotMutHeatmap('single', 'aa', samples=None)
   >>> _ = p.draw()

Simulate variant counts and add to table
++++++++++++++++++++++++++++++++++++++++++++

API documentation
--------------------
"""

import re
import os
import collections
import itertools
import tempfile

import numpy
import scipy
import pandas
import Bio.SeqUtils.ProtParamData

# use plotnine for plotting
from plotnine import *

from dms_tools2.plot import COLOR_BLIND_PALETTE_GRAY, latexSciNot
from dms_tools2 import CODON_TO_AA, CODONS, AAS_WITHSTOP, AA_TO_CODONS


class CodonVariantTable:
    """Associates barcodes with codon mutants of gene.

    Args:
        `barcode_variant_file` (str)
            CSV file giving barcodes and variants. Must have
            columns named "library", "barcode", "substitutions",
            (nucleotide mutations in 1, ... numbering in a format
            like "G301A A302T G856C"), and "variant_call_support"
            (sequences supporting barcode-variant call). Any
            additional columns are removed unless they are specified
            in `extra_cols`.
        `geneseq` (str)
            Sequence of wildtype protein-coding gene.
        `substitutions_are_codon` (bool)
            If `True`, then the "substitutions" column in
            `barcode_variant_file` gives the substitutions
            as codon rather than nucleotide mutations (e.g.,
            "ATG1ATA GTA5CCC" for substitutions at codons
            1 and 5 in 1, 2, ... numbering).
        `extra_cols` (list)
            Additional columns in `barcode_variant_file` to
            retain when creating `barcode_variant_df` and
            `variant_count_df` attributes.

    Attributes:
        `geneseq` (str)
            `geneseq` passed at initialization.
        `sites` (list)
            List of all codon sites in 1, 2, ... numbering.
        `codons` (OrderedDict)
            `codons[r]` is wildtype codon at site `r`, ordered
            by sites.
        `aas` (OrderedDict)
            `aas[r]` is wildtype amino acid at site `r`,
            ordered by sites.
        `libraries` (list)
            List of libraries in `barcode_variantfile`.
        `barcode_variant_df` (pandas DataFrame)
            Info about codon mutations parsed from `barcode_variantfile`.
        `variant_count_df` (pandas DataFrame or None)
            Initially `None`, but after data are added with
            :class:`CodonVariantTable.addSampleCounts`, holds
            counts of each variant for each sample. Differs from
            `barcode_variant_df` in that the former just holds
            initial barcode-variant data, whereas `variant_count_df`
            is updated with variant counts for samples.

    Here is an example.

    First, initialize a :class:`CodonVariantTable`:

    >>> geneseq = 'ATGGGATGA'
    >>> variantfile = '_variantfile.csv'
    >>> with open(variantfile, 'w') as f:
    ...     _ = f.write(
    ...           'library,barcode,substitutions,variant_call_support\\n'
    ...           'lib_1,AAC,,2\\n'
    ...           'lib_1,GAT,G4C A6C,1\\n'
    ...           'lib_2,AAC,T2A G4T,2\\n'
    ...           'lib_2,CAT,A6C,3'
    ...           )
    >>> variants = CodonVariantTable(
    ...             barcode_variant_file=variantfile,
    ...             geneseq=geneseq
    ...             )
    >>> os.remove(variantfile)

    Check attributes of the :class:`CodonVariantTable`:

    >>> variants.sites
    [1, 2, 3]
    >>> variants.codons == {1:'ATG', 2:'GGA', 3:'TGA'}
    True
    >>> variants.aas == {1:'M', 2:'G', 3:'*'}
    True
    >>> variants.libraries
    ['lib_1', 'lib_2']
    >>> variants.valid_barcodes('lib_1') == {'AAC', 'GAT'}
    True
    >>> variants.valid_barcodes('lib_2') == {'AAC', 'CAT'}
    True
    >>> pandas.set_option('display.max_columns', 10)
    >>> pandas.set_option('display.width', 500)
    >>> variants.barcode_variant_df
      library barcode  variant_call_support codon_substitutions aa_substitutions  n_codon_substitutions  n_aa_substitutions
    0   lib_1     AAC                     2                                                           0                   0
    1   lib_1     GAT                     1             GGA2CGC              G2R                      1                   1
    2   lib_2     AAC                     2     ATG1AAG GGA2TGA          M1K G2*                      2                   2
    3   lib_2     CAT                     3             GGA2GGC                                       1                   0

    We can also look at the number of variants; we get this
    by calling :class:`CodonVariantTable.n_variants_df` with
    `samples=None` since we don't have samples, and just want
    stats across barcoded variants:

    >>> variants.n_variants_df(samples=None)
             library             sample  count
    0          lib_1  barcoded variants      2
    1          lib_2  barcoded variants      2
    2  all libraries  barcoded variants      4

    We can also see how these numbers change if we require a
    variant call support of at least 2:

    >>> variants.n_variants_df(samples=None, min_support=2)
             library             sample  count
    0          lib_1  barcoded variants      1
    1          lib_2  barcoded variants      2
    2  all libraries  barcoded variants      3

    If we want to combine the data for the two libraries, we can use
    :class:`CodonVariantTable.addMergedLibraries`, which creates a
    new combined library called "all libraries":

    >>> variants.addMergedLibraries(variants.barcode_variant_df)
             library    barcode  variant_call_support codon_substitutions aa_substitutions  n_codon_substitutions  n_aa_substitutions
    0          lib_1        AAC                     2                                                           0                   0
    1          lib_1        GAT                     1             GGA2CGC              G2R                      1                   1
    2          lib_2        AAC                     2     ATG1AAG GGA2TGA          M1K G2*                      2                   2
    3          lib_2        CAT                     3             GGA2GGC                                       1                   0
    4  all libraries  lib_1-AAC                     2                                                           0                   0
    5  all libraries  lib_1-GAT                     1             GGA2CGC              G2R                      1                   1
    6  all libraries  lib_2-AAC                     2     ATG1AAG GGA2TGA          M1K G2*                      2                   2
    7  all libraries  lib_2-CAT                     3             GGA2GGC                                       1                   0

    Note however that :class:`CodonVariantTable.addMergedLibraries`
    doesn't do anything if there is only one library:

    >>> variants.addMergedLibraries(variants.barcode_variant_df
    ...                             .query('library == "lib_1"'))
      library barcode  variant_call_support codon_substitutions aa_substitutions  n_codon_substitutions  n_aa_substitutions
    0   lib_1     AAC                     2                                                           0                   0
    1   lib_1     GAT                     1             GGA2CGC              G2R                      1                   1

    Count number of barcoded variants with each mutation:

    >>> variants.mutCounts('all', 'aa', samples=None)[ : 2]
      library             sample mutation  count  mutation_type  site
    0   lib_1  barcoded variants      G2R      1  nonsynonymous     2
    1   lib_1  barcoded variants      *3A      0  nonsynonymous     3

    We can do the same for codon mutations (here for only a
    single library), first for all variants:

    >>> variants.mutCounts('all', 'codon', samples=None,
    ...         libraries=['lib_2'])[ : 4]
      library             sample mutation  count  mutation_type  site
    0   lib_2  barcoded variants  ATG1AAG      1  nonsynonymous     1
    1   lib_2  barcoded variants  GGA2GGC      1     synonymous     2
    2   lib_2  barcoded variants  GGA2TGA      1           stop     2
    3   lib_2  barcoded variants  ATG1AAA      0  nonsynonymous     1

    Then for just single-mutant variants:

    >>> variants.mutCounts('single', 'codon', samples=None,
    ...         libraries=['lib_2'])[ : 3]
      library             sample mutation  count  mutation_type  site
    0   lib_2  barcoded variants  GGA2GGC      1     synonymous     2
    1   lib_2  barcoded variants  ATG1AAA      0  nonsynonymous     1
    2   lib_2  barcoded variants  ATG1AAC      0  nonsynonymous     1

    Initially we haven't added any barcode count information
    for any samples:

    >>> all(variants.samples(lib) == [] for lib in variants.libraries)
    True
    >>> variants.variant_count_df is None
    True

    Now we add barcode count information for sample "input"
    from library 1 using :class:`CodonVariantTable.addSampleCounts`:

    >>> counts_lib1_input = pandas.DataFrame(
    ...         {'barcode':['AAC', 'GAT'],
    ...          'count'  :[  253,  1101]})
    >>> variants.addSampleCounts('lib_1', 'input', counts_lib1_input)
    >>> variants.variant_count_df
      barcode  count library sample  variant_call_support codon_substitutions aa_substitutions  n_codon_substitutions  n_aa_substitutions
    0     GAT   1101   lib_1  input                     1             GGA2CGC              G2R                      1                   1
    1     AAC    253   lib_1  input                     2                                                           0                   0

    We get an error if we try to add these same data again,
    as they are already added for that sample to that library:

    >>> variants.addSampleCounts('lib_1', 'input', counts_lib1_input)
    Traceback (most recent call last):
    ...
    ValueError: `library` lib_1 already has `sample` input

    But, we can add barcode counts for another
    sample (named "selected" in this case) to library 1:

    >>> counts_lib1_selected = pandas.DataFrame(
    ...         {'barcode':['AAC', 'GAT'],
    ...          'count'  :[  513,  401]})
    >>> variants.addSampleCounts('lib_1', 'selected', counts_lib1_selected)

    As well as barcode counts for the same two samples
    ("input" and "selected") to our other library (library 2):

    >>> counts_lib2_input = pandas.DataFrame(
    ...         {'barcode':['AAC', 'CAT'],
    ...          'count'  :[ 1253,  923]})
    >>> variants.addSampleCounts('lib_2', 'input', counts_lib2_input)
    >>> counts_lib2_selected = pandas.DataFrame(
    ...         {'barcode':['AAC', 'CAT'],
    ...          'count'  :[  113,  1200]})
    >>> variants.addSampleCounts('lib_2', 'selected', counts_lib2_selected)
    >>> variants.variant_count_df
      barcode  count library    sample  variant_call_support codon_substitutions aa_substitutions  n_codon_substitutions  n_aa_substitutions
    0     GAT   1101   lib_1     input                     1             GGA2CGC              G2R                      1                   1
    1     AAC    253   lib_1     input                     2                                                           0                   0
    2     AAC    513   lib_1  selected                     2                                                           0                   0
    3     GAT    401   lib_1  selected                     1             GGA2CGC              G2R                      1                   1
    4     AAC   1253   lib_2     input                     2     ATG1AAG GGA2TGA          M1K G2*                      2                   2
    5     CAT    923   lib_2     input                     3             GGA2GGC                                       1                   0
    6     CAT   1200   lib_2  selected                     3             GGA2GGC                                       1                   0
    7     AAC    113   lib_2  selected                     2     ATG1AAG GGA2TGA          M1K G2*                      2                   2

    We can also use :class:`CodonVariantTable.mutCounts`
    to look at total counts of each mutation:

    >>> variants.mutCounts('all', 'aa')[ : 2]
      library sample mutation  count  mutation_type  site
    0   lib_1  input      G2R   1101  nonsynonymous     2
    1   lib_1  input      *3A      0  nonsynonymous     3

    >>> variants.mutCounts('all', 'aa', libraries=['lib_2'])[ : 3]
      library sample mutation  count  mutation_type  site
    0   lib_2  input      G2*   1253           stop     2
    1   lib_2  input      M1K   1253  nonsynonymous     1
    2   lib_2  input      *3A      0  nonsynonymous     3

    We can use :class:`CodonVariantTable.writeCodonCounts` to
    write codon count files. First for only **single** mutants:

    >>> with tempfile.TemporaryDirectory() as tmpdir:
    ...     countfiles = variants.writeCodonCounts("single", outdir=tmpdir)
    ...     lib1_input = pandas.read_csv(f'{tmpdir}/lib_1_input_codoncounts.csv')
    ...     all_sel = pandas.read_csv(f'{tmpdir}/all-libraries_selected_codoncounts.csv')

    Make sure we created expected countfiles:

    >>> countfiles.assign(countfile=lambda x: x.countfile.apply(os.path.basename))
             library    sample                               countfile
    0          lib_1     input             lib_1_input_codoncounts.csv
    1          lib_1  selected          lib_1_selected_codoncounts.csv
    2          lib_2     input             lib_2_input_codoncounts.csv
    3          lib_2  selected          lib_2_selected_codoncounts.csv
    4  all-libraries     input     all-libraries_input_codoncounts.csv
    5  all-libraries  selected  all-libraries_selected_codoncounts.csv

    Check for expected values in a few of these counts files, only
    showing columns with non-zero entries:

    >>> lib1_input.iloc[:, (lib1_input != 0).any(axis='rows').values]
       site wildtype  ATG   CGC  GGA  TGA
    0     1      ATG  253     0    0    0
    1     2      GGA    0  1101  253    0
    2     3      TGA    0     0    0  253
    >>> all_sel.iloc[:, (all_sel != 0).any(axis='rows').values]
       site wildtype  ATG  CGC  GGA   GGC  TGA
    0     1      ATG  513    0    0     0    0
    1     2      GGA    0  401  513  1200    0
    2     3      TGA    0    0    0     0  513

    Now write codon counts files for **all** mutants:

    >>> with tempfile.TemporaryDirectory() as tmpdir:
    ...     _ = variants.writeCodonCounts("all", outdir=tmpdir)
    ...     lib1_input_all = pandas.read_csv(f'{tmpdir}/lib_1_input_codoncounts.csv')
    ...     all_sel_all = pandas.read_csv(f'{tmpdir}/all-libraries_selected_codoncounts.csv')
    >>> lib1_input_all.iloc[:, (lib1_input_all != 0).any(axis='rows').values]
       site wildtype   ATG   CGC  GGA   TGA
    0     1      ATG  1354     0    0     0
    1     2      GGA     0  1101  253     0
    2     3      TGA     0     0    0  1354
    >>> all_sel_all.iloc[:, (all_sel_all != 0).any(axis='rows').values]
       site wildtype  AAG   ATG  CGC  GGA   GGC   TGA
    0     1      ATG  113  2114    0    0     0     0
    1     2      GGA    0     0  401  513  1200   113
    2     3      TGA    0     0    0    0     0  2227

    We can also initialize a :class:`CodonVariantTable` from the
    `variant_count_df` if we have written that to a CSV file.
    We do this using :meth:`CodonVariantTable.from_variant_count_df`.
    The example below shows how this newly initialized variant table
    is equal to the original one used to write the CSV file:

    >>> with tempfile.NamedTemporaryFile(mode='w') as f:
    ...     variants.variant_count_df.to_csv(f, index=False)
    ...     f.flush()
    ...     variants_eq = CodonVariantTable.from_variant_count_df(
    ...                     variant_count_df_file=f.name,
    ...                     geneseq=geneseq)
    >>> variants == variants_eq
    True

    Of course, they initializd variant table is **not** equal
    to original one if we don't write the full `variant_count_df`
    to the CSV file:

    >>> with tempfile.NamedTemporaryFile(mode='w') as f:
    ...     (variants
    ...         .variant_count_df.query('sample == "input"')
    ...         .to_csv(f, index=False)
    ...         )
    ...     f.flush()
    ...     variants_ne = CodonVariantTable.from_variant_count_df(
    ...                     variant_count_df_file=f.name,
    ...                     geneseq=geneseq)
    >>> variants == variants_ne
    False
    """

    def __eq__(self, other):
        # following here: https://stackoverflow.com/a/390640
        if type(other) is not type(self):
            return False
        elif self.__dict__.keys() != other.__dict__.keys():
            return False
        else:
            for key, val in self.__dict__.items():
                val2 = getattr(other, key)
                if isinstance(val, pandas.DataFrame):
                    if not val.equals(val2):
                        return False
                else:
                    if val != val2:
                        return False
            return True


    @classmethod
    def from_variant_count_df(cls, *, variant_count_df_file, geneseq,
            drop_all_libs=True):
        """:class:`CodonVariantTable` from CSV of `variant_count_df`.

        Use this method when you have written a CSV file of the
        `variant_count_df` attribute of a :class:`CodonVariantTable`,
        and now wish to re-initialize that :class:`CodonVariantTable`.

        Args:
            `variant_count_df_file` (str)
                Name of CSV file containing the `variant_count_df`.
                Must have following columns: "barcode", "library",
                "variant_call_support", "codon_substitutions",
                "sample", and "count".
            `geneseq` (str)
                Sequence of wildtype protein-coding gene.
            `drop_all_libs` (bool)
                If there is a library named "all libraries",
                drop it from the list of libraries in the created
                :class:`CodonVariantTable` as it probably was
                added by :meth:`CodonVariantTable.addMergedLibraries`
                and duplicates information for the individual libraries.

        Returns:
            The :class:`CodonVariantTable` used to write
            `variant_count_df_file`.
        """
        df = pandas.read_csv(variant_count_df_file)

        req_cols = ['barcode', 'library', 'variant_call_support',
                    'codon_substitutions', 'sample', 'count']
        if not (set(req_cols) < set(df.columns)):
            raise ValueError(f"{variant_count_df} lacks required "
                             f"columns {req_cols}")
        else:
            df = df[req_cols]

        if drop_all_libs:
            dropcol = "all libraries"
            if dropcol in df['library'].unique():
                df = df.query('library != @dropcol')

        with tempfile.NamedTemporaryFile(mode='w') as f:
            (df
                .drop(columns=['sample', 'count'])
                .rename(columns={'codon_substitutions':'substitutions'})
                .drop_duplicates()
                .to_csv(f, index=False)
                )
            f.flush()
            cvt = cls(barcode_variant_file=f.name,
                      geneseq=geneseq,
                      substitutions_are_codon=True)

        for sample in df['sample'].unique():
            for lib in cvt.libraries:
                idf = df.query('sample == @sample & library == @lib')
                if len(idf):
                    cvt.addSampleCounts(lib,
                                        sample,
                                        idf[['barcode', 'count']]
                                        )

        return cvt


    def __init__(self, *, barcode_variant_file, geneseq,
                 substitutions_are_codon=False, extra_cols=[]):
        """See main class doc string."""

        self.geneseq = geneseq.upper()
        if not re.match('^[ATGC]+$', self.geneseq):
            raise ValueError(f"invalid nucleotides in {self.geneseq}")
        if ((len(geneseq) % 3) != 0) or len(geneseq) == 0:
            raise ValueError(f"`geneseq` of invalid length {len(self.geneseq)}")
        self.sites = list(range(1, len(self.geneseq) // 3 + 1))
        self.codons = collections.OrderedDict([
                (r, self.geneseq[3 * (r - 1) : 3 * r]) for r in self.sites])
        self.aas = collections.OrderedDict([
                (r, CODON_TO_AA[codon]) for r, codon in self.codons.items()])

        df = pandas.read_csv(barcode_variant_file)
        required_cols = ['library', 'barcode',
                         'substitutions', 'variant_call_support']
        if not set(df.columns).issuperset(set(required_cols)):
            raise ValueError("`variantfile` does not have "
                             f"required columns {required_cols}")
        if extra_cols and not set(df.columns).issuperset(set(extra_cols)):
            raise ValueError("`variantfile` does not have "
                             f"`extra_cols` {extra_cols}")
        df = df[required_cols + extra_cols]

        self.libraries = sorted(df.library.unique().tolist())
        self._valid_barcodes = {}
        for lib in self.libraries:
            barcodes = df.query('library == @lib').barcode
            if len(set(barcodes)) != len(barcodes):
                raise ValueError(f"duplicated barcodes for {lib}")
            self._valid_barcodes[lib] = set(barcodes)

        self._samples = {lib:[] for lib in self.libraries}
        self.variant_count_df = None

        if substitutions_are_codon:
            codonSubsFunc = lambda x: x.upper()
        else:
            codonSubsFunc = self._ntToCodonMuts

        self.barcode_variant_df = (
                df
                # info about codon and amino-acid substitutions
                .assign(codon_substitutions=
                            lambda x: x.substitutions
                                       .fillna('')
                                       .apply(codonSubsFunc),
                        aa_substitutions=
                            lambda x: x.codon_substitutions
                                       .apply(self.codonToAAMuts),
                        n_codon_substitutions=
                            lambda x: x.codon_substitutions
                                       .str
                                       .split()
                                       .apply(len),
                        n_aa_substitutions=
                            lambda x: x.aa_substitutions
                                       .str
                                       .split()
                                       .apply(len)
                        )
                # we no longer need initial `substitutions` column
                .drop('substitutions', axis='columns')
                # sort to ensure consistent order
                .assign(library=lambda x:
                                pandas.Categorical(
                                    x.library,
                                    categories=self.libraries,
                                    ordered=True
                                    )
                        )
                .sort_values(['library', 'barcode'])
                .reset_index(drop=True)
                )

        # check validity of codon substitutions given `geneseq`
        for codonmut in itertools.chain.from_iterable(
                        self.barcode_variant_df
                          .codon_substitutions.str.split()):
            m = re.match('^(?P<wt>[ATGC]{3})'
                          '(?P<r>\d+)'
                          '(?P<mut>[ATGC]{3})$',
                         codonmut)
            if m is None:
                raise ValueError(f"invalid mutation {codonmut}")
            wt = m.group('wt')
            r = int(m.group('r'))
            mut = m.group('mut')
            if r not in self.sites:
                raise ValueError(f"invalid site {r} in codon mutation {codonmut}")
            if self.codons[r] != wt:
                raise ValueError(f"Wrong wildtype codon in {codonmut}. "
                                 f"Expected wildtype of {codons[r]}.")
            if wt == mut:
                raise ValueError(f"invalid mutation {codonmut}")

        # define some colors for plotting
        self._mutation_type_colors = {
                'nonsynonymous':COLOR_BLIND_PALETTE_GRAY[1],
                'synonymous':COLOR_BLIND_PALETTE_GRAY[2],
                'stop':COLOR_BLIND_PALETTE_GRAY[3]
                }


    def samples(self, library):
        """List of all samples for `library`.

        Args:
            `library` (str)
                Valid `library` for the :class:`CodonVariantTable`.

        Returns:
            List of all samples for which barcode counts have
            been added.
        """
        try:
            return self._samples[library]
        except KeyError:
            raise ValueError(f"invalid `library` {library}")


    def addSampleCounts(self, library, sample, barcodecounts):
        """Add variant counts for a sample to `variant_count_df`.

        Args:
            `library` (str)
                Valid `library` for the :class:`CodonVariantTable`.
            `sample` (str)
                Sample name, must **not** already be in
                :class:`CodonVariantTable.samples` for `library`.
            `barcodecounts` (pandas DataFrame)
                Gives counts for each variant by barcode. Must
                have columns named "barcode" and "count". The
                "barcode" column must contain all the barcodes
                in :class:`CodonVariantTable.valid_barcodes` for
                `library`. Such data frames are returned
                by :class:`dms_tools2.barcodes.IlluminaBarcodeParser.parse`.
        """
        if library not in self.libraries:
            raise ValueError(f"invalid library {library}")

        if sample in self.samples(library):
            raise ValueError(f"`library` {library} already "
                             f"has `sample` {sample}")

        req_cols = ['barcode', 'count']
        if not set(barcodecounts.columns).issuperset(set(req_cols)):
            raise ValueError(f"`barcodecounts` lacks columns {req_cols}")
        if len(barcodecounts) != len(set(barcodecounts.barcode.unique())):
            raise ValueError("`barcodecounts` has non-unique barcodes")
        if set(barcodecounts.barcode.unique()) != self.valid_barcodes(library):
            raise ValueError("barcodes in `barcodecounts` do not match "
                             f"those expected for `library` {library}")

        self._samples[library].append(sample)

        df = (barcodecounts
              [req_cols]
              .assign(library=library, sample=sample)
              .merge(self.barcode_variant_df,
                     how='inner',
                     on=['library', 'barcode'],
                     sort=False,
                     validate='one_to_one')
              )

        if self.variant_count_df is None:
            self.variant_count_df = df
        else:
            self.variant_count_df = pandas.concat(
                              [self.variant_count_df, df],
                              axis='index',
                              ignore_index=True,
                              sort=False
                              )

        # samples in order added after ordering by library, getting
        # unique ones as here: https://stackoverflow.com/a/39835527
        unique_samples = list(collections.OrderedDict.fromkeys(
                itertools.chain.from_iterable(
                    [self.samples(lib) for lib in self.libraries])
                ))

        # make library and sample categorical and sort
        self.variant_count_df = (
                self.variant_count_df
                .assign(library=lambda x:
                                pandas.Categorical(
                                    x.library,
                                    categories=self.libraries,
                                    ordered=True
                                    ),
                        sample=lambda x:
                               pandas.Categorical(
                                    x['sample'],
                                    categories=unique_samples,
                                    ordered=True
                                    )
                         )
                .sort_values(['library', 'sample', 'count'],
                             ascending=[True, True, False])
                .reset_index(drop=True)
                )


    def valid_barcodes(self, library):
        """Set of valid barcodes for `library`."""
        if library not in self.libraries:
            raise ValueError(f"invalid `library` {library}; "
                             f"valid libraries are {self.libraries}")
        else:
            return self._valid_barcodes[library]


    def n_variants_df(self, *, libraries='all', samples='all',
                      min_support=1):
        """Number of variants per library / sample.

        Args:
            Same meaning as for
            :class:`CodonVariantTable.plotNumMutsHistogram`.

        Returns:
            DataFrame giving number of variants per library /
            sample.
        """
        df, nlibraries, nsamples = self._getPlotData(libraries,
                                                     samples,
                                                     min_support)

        return (df
                .groupby(['library', 'sample'])
                .aggregate({'count':'sum'})
                .reset_index()
                )


    def mutCounts(self, variant_type, mut_type, *,
            libraries='all', samples='all', min_support=1):
        """Get counts of each individual mutation.

        Args:
            `variant_type` ("single" or "all")
                Include just single mutants, or all mutants?
            Other args:
                Same meaning as for
                :class:`CodonVariantTable.plotNumMutsHistogram`

        Returns:
            A tidy data frame with columns named "library",
            "sample", "mutation", "count", "mutation_type",
            and "site".
        """
        df, nlibraries, nsamples = self._getPlotData(libraries,
                                                     samples,
                                                     min_support)

        samplelist = df['sample'].unique().tolist()
        librarylist = df['library'].unique().tolist()

        if mut_type == 'codon':
            wts = self.codons
            chars = CODONS
            mutation_types = ['nonsynonymous', 'synonymous', 'stop']
        elif mut_type == 'aa':
            wts = self.aas
            chars = AAS_WITHSTOP
            mutation_types = ['nonsynonymous', 'stop']
        else:
            raise ValueError(f"invalid mut_type {mut_type}")

        # data frame listing all mutations with count 0
        mut_list = []
        for r, wt in wts.items():
            for mut in chars:
                if mut != wt:
                    mut_list.append(f'{wt}{r}{mut}')
        all_muts = pandas.concat([
                    pandas.DataFrame({'mutation':mut_list,
                                     'library':library,
                                     'sample':sample,
                                     'count':0})
                    for library, sample in
                    itertools.product(librarylist, samplelist)])

        if variant_type == 'single':
            df = df.query(f'n_{mut_type}_substitutions == 1')
        elif variant_type == 'all':
            df = df.query(f'n_{mut_type}_substitutions >= 1')
        else:
            raise ValueError(f"invalid variant_type {variant_type}")

        def _classify_mutation(mut_str):
            if mut_type == 'aa':
                m = re.match('^(?P<wt>[A-Z\*])\d+(?P<mut>[A-Z\*])$',
                             mut_str)
                wt_aa = m.group('wt')
                mut_aa = m.group('mut')
            else:
                m = re.match('^(?P<wt>[ACTG]{3})\d+(?P<mut>[ACTG]{3})$',
                             mut_str)
                wt_aa = CODON_TO_AA[m.group('wt')]
                mut_aa = CODON_TO_AA[m.group('mut')]
            if wt_aa == mut_aa:
                return 'synonymous'
            elif mut_aa == '*':
                return 'stop'
            else:
                return 'nonsynonymous'

        def _get_site(mut_str):
            if mut_type == 'aa':
                m = re.match('^[A-Z\*](?P<site>\d+)[A-Z\*]$', mut_str)
            else:
                m = re.match('^[ACTG]{3}(?P<site>\d+)[ACTG]{3}$', mut_str)
            site = int(m.group('site'))
            assert site in self.sites
            return site

        df = (df
              .rename(columns=
                  {f"{mut_type}_substitutions":"mutation"}
                  )
              [['library', 'sample', 'mutation', 'count']]
              .pipe(tidy_split, column='mutation')
              .merge(all_muts, how='outer')
              .groupby(['library', 'sample', 'mutation'])
              .aggregate({'count':'sum'})
              .reset_index()
              .assign(
                library=lambda x:
                         pandas.Categorical(
                          x['library'],
                          librarylist,
                          ordered=True),
                sample=lambda x:
                         pandas.Categorical(
                          x['sample'],
                          samplelist,
                          ordered=True),
                mutation_type=lambda x:
                         pandas.Categorical(
                          x['mutation'].apply(_classify_mutation),
                          mutation_types,
                          ordered=True),
                site=lambda x: x['mutation'].apply(_get_site),
                )
              .sort_values(
                ['library', 'sample', 'count', 'mutation'],
                ascending=[True, True, False, True])
              .reset_index(drop=True)
              )

        return df


    def plotMutHeatmap(self, variant_type, mut_type, *,
            count_or_frequency='frequency',
            libraries='all', samples='all', plotfile=None,
            orientation='h', widthscale=1, heightscale=1,
            min_support=1):
        """Heatmap of mutation counts or frequencies.

        Args:
            `count_or_frequency` ("count" or "frequency")
                Plot mutation counts or frequencies?
            All other args
                Same meaning as for
                :meth:`CodonVariantTable.plotCumulMutCoverage`.

        Returns:
            A `plotnine <https://plotnine.readthedocs.io/en/stable/>`_
            plot; can be displayed in a Jupyter notebook with `p.draw()`.
        """

        df = self.mutCounts(variant_type, mut_type, samples=samples,
                            libraries=libraries, min_support=min_support)

        n_variants = (self.n_variants_df(libraries=libraries,
                                         samples=samples,
                                         min_support=min_support)
                      .rename(columns={'count':'nseqs'})
                      )

        # order amino acids by Kyte-Doolittle hydrophobicity,
        # codons by the amino acid they encode
        aa_order = [tup[0] for tup in sorted(
                    Bio.SeqUtils.ProtParamData.kd.items(),
                    key=lambda tup: tup[1])] + ['*']
        codon_order = list(itertools.chain.from_iterable(
                        [AA_TO_CODONS[aa] for aa in aa_order]))

        df = (df
              [['library', 'sample', 'mutation', 'site', 'count']]
              .merge(n_variants, on=['library', 'sample'])
              .assign(frequency=lambda x: x['count'] / x['nseqs'],
                      mut_char=lambda x:
                        pandas.Categorical(
                         x.mutation.str.extract(
                            '^[A-Z\*]+\d+(?P<mut_char>[A-Z\*]+)$')
                            .mut_char,
                         {'aa':aa_order, 'codon':codon_order}[mut_type],
                         ordered=True)
                      )
              )

        if count_or_frequency not in {'count', 'frequency'}:
            raise ValueError(f"invalid count_or_frequency "
                             f"{count_or_frequency}")

        nlibraries = len(df['library'].unique())
        nsamples = len(df['sample'].unique())

        if mut_type == 'codon':
            height_per = 5.5
            mut_desc = 'codon'
        elif mut_type == 'aa':
            height_per = 1.7
            mut_desc = 'amino acid'
        else:
            raise ValueError(f"invalid mut_type {mut_type}")

        if orientation == 'h':
            facet_str = 'sample ~ library'
            width = widthscale * (1.6 + 3.5 * nlibraries)
            height = heightscale * (0.8 + height_per * nsamples)
        elif orientation == 'v':
            facet_str = 'library ~ sample'
            width = widthscale * (1.6 + 3.5 * nsamples)
            height = heightscale * (0.8 + height_per * nlibraries)
        else:
            raise ValueError(f"invalid `orientation` {orientation}")

        p = (ggplot(df, aes('site', 'mut_char',
                            fill=count_or_frequency)) +
             geom_tile() +
             facet_grid(facet_str) +
             theme(figure_size=(width, height),
                   legend_key=element_blank(),
                   axis_text_y=element_text(size=6)
                   ) +
             scale_x_continuous(
                name=f'{mut_desc} site',
                limits=(min(self.sites) - 1, max(self.sites) + 1),
                expand=(0, 0)
                ) +
             ylab(mut_desc) +
             scale_fill_cmap('gnuplot')
             )


        if plotfile:
            p.save(plotfile, height=height, width=width,
                   verbose=False, limitsize=False)

        return p


    def plotMutFreqs(self, variant_type, mut_type, *,
            libraries='all', samples='all', plotfile=None,
            orientation='h', widthscale=1, heightscale=1,
            min_support=1):
        """Mutation frequency along length of gene.

        Args:
            Args have same meaning as for
            :meth:`CodonVariantTable.plotCumulMutCoverage`.

        Returns:
            A `plotnine <https://plotnine.readthedocs.io/en/stable/>`_
            plot; can be displayed in a Jupyter notebook with `p.draw()`.
        """

        df = self.mutCounts(variant_type, mut_type, samples=samples,
                            libraries=libraries, min_support=min_support)

        n_variants = (self.n_variants_df(libraries=libraries,
                                         samples=samples,
                                         min_support=min_support)
                      .rename(columns={'count':'nseqs'})
                      )

        df = (df
              .groupby(['library', 'sample', 'mutation_type', 'site'])
              .aggregate({'count':'sum'})
              .reset_index()
              .merge(n_variants, on=['library', 'sample'])
              .assign(freq=lambda x: x['count'] / x['nseqs'])
              )

        nlibraries = len(df['library'].unique())
        nsamples = len(df['sample'].unique())

        if orientation == 'h':
            facet_str = 'sample ~ library'
            width = widthscale * (1.6 + 1.8 * nlibraries)
            height = heightscale * (0.8 + 1 * nsamples)
        elif orientation == 'v':
            facet_str = 'library ~ sample'
            width = widthscale * (1.6 + 1.8 * nsamples)
            height = heightscale * (0.8 + 1 * nlibraries)
        else:
            raise ValueError(f"invalid `orientation` {orientation}")

        if mut_type == 'aa':
            mut_desc = 'amino-acid'
        else:
            mut_desc = mut_type

        if height < 3:
            ylabel = (f'{mut_desc} mutation\nfrequency '
                      f'({variant_type} mutants)')
        else:
            ylabel = (f'{mut_desc} mutation frequency '
                      f'({variant_type} mutants)')

        p = (ggplot(df, aes('site', 'freq', color='mutation_type')) +
             geom_step() +
             scale_color_manual(
                [self._mutation_type_colors[m] for m in
                 df.mutation_type.unique().sort_values().tolist()],
                name='mutation type'
                ) +
             scale_x_continuous(
                name=f'{mut_desc} site',
                limits=(min(self.sites), max(self.sites))
                ) +
             ylab(ylabel) +
             facet_grid(facet_str) +
             theme(figure_size=(width, height),
                   legend_key=element_blank(),
                   legend_text=element_text(size=11)
                   )
             )

        if plotfile:
            p.save(plotfile, height=height, width=width, verbose=False)

        return p


    def plotCumulMutCoverage(self, variant_type, mut_type, *,
            libraries='all', samples='all', plotfile=None,
            orientation='h', widthscale=1, heightscale=1,
            min_support=1, max_count=None):
        """Fraction of mutation seen <= some number of times.

        Args:
            `variant_type` ("single" or "all")
                Include just single mutants, or all mutants?
                Mutations are counted relative to `mut_type`.
            `max_count` (`None` or int)
                Plot cumulative fraction plot out to this
                number of observations of mutation. If `None`,
                a reasonable value is automatically determined.
            Other args:
                Same meaning as for
                :class:`CodonVariantTable.plotNumMutsHistogram`

        Returns:
            A `plotnine <https://plotnine.readthedocs.io/en/stable/>`_
            plot; can be displayed in a Jupyter notebook with `p.draw()`.
        """

        df = self.mutCounts(variant_type, mut_type, samples=samples,
                            libraries=libraries, min_support=min_support)

        # add one to counts to plot fraction found < this many
        # as stat_ecdf by default does <=
        df = df.assign(count=lambda x: x['count'] + 1)

        if max_count is None:
            max_count = df['count'].quantile(0.75)

        nlibraries = len(df['library'].unique())
        nsamples = len(df['sample'].unique())

        if orientation == 'h':
            facet_str = 'sample ~ library'
            width = widthscale * (1.6 + 1.3 * nlibraries)
            height = heightscale * (1 + 1.2 * nsamples)
        elif orientation == 'v':
            facet_str = 'library ~ sample'
            width = widthscale * (1.6 + 1.3 * nsamples)
            height = heightscale * (1 + 1.2 * nlibraries)
        else:
            raise ValueError(f"invalid `orientation` {orientation}")

        if width > 4:
            xlabel = f'counts among {variant_type} mutants'
        else:
            xlabel = f'counts among\n{variant_type} mutants'

        mut_desc = {'aa':'amino-acid', 'codon':'codon'}[mut_type]
        if height > 3:
            ylabel = f'frac {mut_desc} mutations found < this many times'
        else:
            ylabel = f'frac {mut_desc} mutations\nfound < this many times'

        p = (ggplot(df, aes('count', color='mutation_type')) +
             stat_ecdf(geom='step', size=0.75) +
             coord_cartesian(xlim=(0, max_count)) +
             scale_color_manual(
                [self._mutation_type_colors[m] for m in
                 df.mutation_type.unique().sort_values().tolist()],
                name='mutation type'
                ) +
             xlab(xlabel) +
             ylab(ylabel) +
             facet_grid(facet_str) +
             theme(figure_size=(width, height),
                   legend_key=element_blank(),
                   legend_text=element_text(size=11),
                   axis_text_x=element_text(angle=90),
                   )
             )

        if plotfile:
            p.save(plotfile, height=height, width=width, verbose=False)

        return p


    def plotNumCodonMutsByType(self, variant_type, *,
            libraries='all', samples='all', plotfile=None,
            orientation='h', widthscale=1, heightscale=1,
            min_support=1):
        """Nonsynonymous, synonymous, stop mutations per variant.

        Args:
            `variant_type` ("single" or "all")
                Include just single-codon mutants and wildtype,
                or include all mutants.
            Other args:
                Same meaning as for
                :class:`CodonVariantTable.plotNumMutsHistogram`

        Returns:
            A `plotnine <https://plotnine.readthedocs.io/en/stable/>`_
            plot; can be displayed in a Jupyter notebook with `p.draw()`.
        """
        df, nlibraries, nsamples = self._getPlotData(libraries,
                                                     samples,
                                                     min_support)

        if variant_type == 'single':
            df = df.query('n_codon_substitutions <= 1')
        elif variant_type != 'all':
            raise ValueError(f"invalid variant_type {variant_type}")

        if orientation == 'h':
            facet_str = 'sample ~ library'
            width = widthscale * (1 + 1.4 * nlibraries)
            height = heightscale * (1 + 1.3 * nsamples)
        elif orientation == 'v':
            facet_str = 'library ~ sample'
            width = widthscale * (1 + 1.4 * nsamples)
            height = heightscale * (1 + 1.3 * nlibraries)
        else:
            raise ValueError(f"invalid `orientation` {orientation}")

        if height > 3:
            ylabel = f'mutations per variant ({variant_type} mutants)'
        else:
            ylabel = f'mutations per variant\n({variant_type} mutants)'

        codon_mut_types = ['nonsynonymous', 'synonymous', 'stop']

        # mutations from stop to another amino-acid counted as nonsyn
        df = (df
              .assign(
                synonymous=lambda x: x.n_codon_substitutions -
                                     x.n_aa_substitutions,
                stop=lambda x: x.aa_substitutions.str
                               .findall('[A-Z]\d+\*').apply(len),
                nonsynonymous=lambda x: x.n_codon_substitutions -
                                        x.synonymous - x.stop
                )
              .melt(id_vars=['library', 'sample', 'count'],
                    value_vars=codon_mut_types,
                    var_name='mutation_type',
                    value_name='num_muts')
              .assign(
                  mutation_type=lambda x:
                                pandas.Categorical(
                                 x.mutation_type,
                                 categories=codon_mut_types,
                                 ordered=True),
                  num_muts_count=lambda x: x.num_muts * x['count']
                  )
              .groupby(['library', 'sample', 'mutation_type'])
              .aggregate({'num_muts_count':'sum', 'count':'sum'})
              .reset_index()
              .assign(number=lambda x: x.num_muts_count / x['count'])
              )

        p = (ggplot(df, aes('mutation_type', 'number',
                            fill='mutation_type', label='number')) +
             geom_bar(stat='identity') +
             geom_text(size=8, va='bottom', format_string='{0:.3f}') +
             facet_grid(facet_str) +
             scale_y_continuous(name=ylabel,
                                expand=(0.03, 0, 0.12, 0)) +
             scale_fill_manual(
                [self._mutation_type_colors[m] for m in
                 df.mutation_type.unique().sort_values().tolist()]
                ) +
             theme(figure_size=(width, height),
                   axis_title_x=element_blank(),
                   axis_text_x=element_text(angle=90, size=11),
                   legend_position='none')
             )

        if plotfile:
            p.save(plotfile, height=height, width=width, verbose=False)

        return p


    def plotNumMutsHistogram(self, mut_type, *,
            libraries='all', samples='all', plotfile=None,
            orientation='h', widthscale=1, heightscale=1,
            min_support=1, max_muts=None):
        """Plot histograms of number of mutations per variant.

        Args:
            `mut_type` (str)
                Type of mutation: "codon" or "aa".
            `libraries` ("all", "all_only", or list)
                Set to "all" to include all libraries including
                a merge of the libraries, "all_only" to only
                include the merge of the libraries, or a list
                of libraries.
            `samples` (the str "all", `None`, or list)
                Set to "all" to include all samples, `None` to
                just count each barcoded variant once, or specify
                a list of samples.
            `plotfile` (`None` or str)
                Name of file to which we save plot.
            `orientation` (the str 'h' or 'v')
                Do we facet libraries horizontally or vertically?
            `widthscale` (float or int)
                Expand width of plot by this much.
            `heightscale` (float or int)
                Expand height of plot by this much.
            `min_support` (int)
                Only include variants with at least this large
                of a variant call support.
            `max_muts` (int or `None`)
                In histogram, group together all variants with
                >= this many mutations; set to `None` for no
                cutoff.

        Returns:
            `plotnine <https://plotnine.readthedocs.io/en/stable/>`_
            plot; can be displayed in a Jupyter notebook with `p.draw()`.
        """
        df, nlibraries, nsamples = self._getPlotData(libraries,
                                                     samples,
                                                     min_support)

        if mut_type == 'aa':
            mut_col = 'n_aa_substitutions'
            xlabel = 'amino-acid mutations'
        elif mut_type == 'codon':
            mut_col = 'n_codon_substitutions'
            xlabel = 'codon mutations'
        else:
            raise ValueError(f"invalid mut_type {mut_type}")

        if orientation == 'h':
            facet_str = 'sample ~ library'
            width = widthscale * (1 + 1.5 * nlibraries)
            height = heightscale * (0.6 + 1.5 * nsamples)
        elif orientation == 'v':
            facet_str = 'library ~ sample'
            width = widthscale * (1 + 1.5 * nsamples)
            height = heightscale * (0.6 + 1.5 * nlibraries)
        else:
            raise ValueError(f"invalid `orientation` {orientation}")

        df[mut_col] = numpy.clip(df[mut_col], None, max_muts)

        df = (df
              .groupby(['library', 'sample', mut_col])
              .aggregate({'count':'sum'})
              .reset_index()
              )

        p = (ggplot(df, aes(mut_col, 'count')) +
             geom_bar(stat='identity') +
             facet_grid(facet_str) +
             xlab(xlabel) +
             scale_y_continuous(labels=latexSciNot) +
             theme(figure_size=(width, height))
             )

        if plotfile:
            p.save(plotfile, height=height, width=width, verbose=False)

        return p


    def writeCodonCounts(self, single_or_all, *, outdir=None):
        """Writes codon counts files for all libraries and samples.

        Files are written in
        `format <https://jbloomlab.github.io/dms_tools2/dms2_bcsubamp.html#id12>`_
        produced by
        `dms2_bcsubamp <https://jbloomlab.github.io/dms_tools2/dms2_bcsubamp.html>`_.

        Args:
            `single_or_all` ("single" or "all")
                If "single", then counts are just from single-codon
                mutants and wildtype, and we count all wildtype codons
                and just mutated codon for single-codon mutants. If
                "all", we count all codons for all variants at all
                sites. This is appropriate if enrichment of each mutation
                is supposed to represent its effect for "single", and if
                enrichment of mutation is supposed to represent its
                average effect across genetic backgrounds in the library
                for "all" provided mutations are Poisson distributed.
            `outdir` (`None` or str)
                Name of directory into which we write counts,
                created if it does not exist. Use `None` to
                write to current directory.

        Returns:
            A pandas data frame with the columns "library", "sample",
            ("all-libraries" is one library) and "countfile".
            The "countfile" column gives the name of the created
            file, which is ``<library>_<sample>_codoncounts.csv``.
        """

        codonmatch = re.compile(f'^(?P<wt>{"|".join(CODONS)})'
                                 '(?P<r>\d+)'
                                f'(?P<mut>{"|".join(CODONS)})$'
                                )
        def _parseCodonMut(mutstr):
            m = codonmatch.match(mutstr)
            return (m.group('wt'), int(m.group('r')), m.group('mut'))

        if self.variant_count_df is None:
            raise ValueError("no samples with counts")

        if single_or_all not in {'single', 'all'}:
            raise ValueError(f"invalid `single_or_all` {single_or_all}")

        if outdir is not None:
            os.makedirs(outdir, exist_ok=True)
        else:
            outdir = ''

        df = self.addMergedLibraries(self.variant_count_df,
                                     all_lib='all-libraries')

        countfiles = []
        liblist = []
        samplelist = []

        for lib, sample in itertools.product(
                            df['library'].unique().tolist(),
                            df['sample'].unique().tolist()
                            ):

            countfile = os.path.join(outdir,
                            f'{lib}_{sample}_codoncounts.csv')
            countfiles.append(countfile)
            liblist.append(lib)
            samplelist.append(sample)

            i_df = df.query('library == @lib & sample == @sample')

            codoncounts = {codon:[0] * len(self.sites) for codon
                           in CODONS}

            if single_or_all == 'single':
                n_wt = (i_df
                        .query('n_codon_substitutions == 0')
                        ['count']
                        .sum()
                        )
                for isite, site in enumerate(self.sites):
                    codoncounts[self.codons[site]][isite] += n_wt
                for mut, count in (i_df
                                   .query('n_codon_substitutions == 1')
                                   [['codon_substitutions', 'count']]
                                   .itertuples(index=False, name=None)
                                   ):
                    wtcodon, r, mutcodon = _parseCodonMut(mut)
                    codoncounts[mutcodon][r - 1] += count

            elif single_or_all == 'all':
                n_wt = i_df['count'].sum()
                for isite, site in enumerate(self.sites):
                    codoncounts[self.codons[site]][isite] += n_wt
                for muts, count in (i_df
                                   .query('n_codon_substitutions > 0')
                                   [['codon_substitutions', 'count']]
                                   .itertuples(index=False, name=None)
                                   ):
                    for mut in muts.split():
                        wtcodon, r, mutcodon = _parseCodonMut(mut)
                        codoncounts[mutcodon][r - 1] += count
                        codoncounts[wtcodon][r - 1] -= count

            else:
                raise ValueError(f"invalid `single_or_all` {single_or_all}")

            counts_df = pandas.DataFrame(collections.OrderedDict(
                         [('site', self.sites),
                          ('wildtype', [self.codons[r] for r in self.sites])] +
                         [(codon, codoncounts[codon]) for codon in CODONS]
                         ))
            counts_df.to_csv(countfile, index=False)

        assert all(map(os.path.isfile, countfiles))

        return pandas.DataFrame({'library':liblist,
                                 'sample':samplelist,
                                 'countfile':countfiles})


    @staticmethod
    def addMergedLibraries(df, *, all_lib='all libraries'):
        """Add data to `df` for all libraries merged.

        Args:
            `df` (pandas DataFrame)
                DataFrame that includes columns named
                "library" and "barcode".
            `all_lib` (str)
                Name given to library that is merge of all
                other libraries.

        Returns:
            If `df` only has data for one library, just returns
            `df`. Otherwise returns a copy of `df` that has a
            new library with the name given by `all_lib`
            and contains the data for all individual libraries,
            with the "barcode" column giving the original
            library name followed by a hyphen and the barcode.
        """
        libs = df.library.unique().tolist()

        if len(libs) <= 1:
            return df

        if all_lib in libs:
            raise ValueError(f"library {all_lib} already exists")

        df = (pandas.concat([df,
                             df.assign(
                                barcode=lambda x:
                                    x.library.str
                                     .cat(x.barcode, sep='-'),
                                library=all_lib)
                             ],
                            axis='index',
                            ignore_index=True,
                            sort=False)
              .assign(library=lambda x:
                              pandas.Categorical(
                               x['library'],
                               categories=libs + [all_lib],
                               ordered=True)
                      )
              )

        return df


    def _getPlotData(self, libraries, samples, min_support):
        """Gets data to plot from library and sample filters.

        Args:
            `libraries`, `samples`, `min_support` have meaning as
            for :class:`CodonVariantTable.plotNumMutsHistogram`.

        Returns:
            The 3-tuple `(df, nlibraries, nsamples)` where:

                - `df`: DataFrame with data to plot.

                - `nlibraries`: number of libraries being plotted.

                - `nsamples`: number of samples being plotted.
        """
        if samples is None:
            df = (self.barcode_variant_df
                  .assign(sample='barcoded variants')
                  .assign(count=1)
                  )
        elif samples == 'all':
            if self.variant_count_df is None:
                raise ValueError('no samples have been added')
            df = self.variant_count_df
        elif isinstance(samples, list):
            all_samples = set(itertools.chain.from_iterable(
                    self.samples(lib) for lib in self.libraries))
            if not all_samples.issuperset(set(samples)):
                raise ValueError(f"invalid sample(s) in {samples}")
            if len(samples) != len(set(samples)):
                raise ValueError(f"duplicate samples in {samples}")
            df = (df
                  .query('sample in @samples')
                  .assign(sample=lambda x:
                                 pandas.Categorical(
                                   x['sample'],
                                   categories=samples,
                                   ordered=True)
                          )
                  )
        else:
            raise ValueError(f"invalid `samples` {samples}")

        df = df.query('variant_call_support >= @min_support')

        if not len(df):
            raise ValueError(f"no samples {samples}")
        else:
            nsamples = len(df['sample'].unique())

        if libraries == 'all':
            df = self.addMergedLibraries(df)
        elif libraries == 'all_only':
            df = (self.addMergedLibraries(df)
                  .query('library == "all libraries"')
                  )
        elif isinstance(libraries, list):
            if not set(self.libraries).issuperset(set(libraries)):
                raise ValueError(f"invalid library in {libraries}")
            if len(libraries) != len(set(libraries)):
                raise ValueError(f"duplicate library in {libraries}")
            df = (df
                  .query('library in @libraries')
                  .assign(library=lambda x:
                                  pandas.Categorical(
                                   x['library'],
                                   categories=libraries,
                                   ordered=True)
                          )
                  )
        else:
            raise ValueError(f"invalid `libraries` {libraries}")
        if not len(df):
            raise ValueError(f"no libraries {libraries}")
        else:
            nlibraries = len(df['library'].unique())

        return (df, nlibraries, nsamples)


    @classmethod
    def codonToAAMuts(self, codon_mut_str):
        """Converts string of codon mutations to amino-acid mutations.

        Args:
            `codon_mut_str` (str)
                Codon mutations, delimited by a space and in
                1, 2, ... numbering.

        Returns:
            String with amino acid mutations in 1, 2, ... numbering.

        >>> CodonVariantTable.codonToAAMuts('ATG1GTG GGA2GGC TGA3AGA')
        'M1V *3R'
        """
        aa_muts = {}
        for mut in codon_mut_str.upper().split():
            m = re.match('^(?P<wt>[ATGC]{3})(?P<r>\d+)(?P<mut>[ATGC]{3})$',
                         mut)
            if not m:
                raise ValueError(f"invalid mutation {mut} in {codon_mut_str}")
            r = int(m.group('r'))
            if r in aa_muts:
                raise ValueError(f"duplicate codon mutation for {r}")
            wt_codon = m.group('wt')
            mut_codon = m.group('mut')
            if wt_codon == mut_codon:
                raise ValueError(f"invalid mutation {mut}")
            wt_aa = CODON_TO_AA[wt_codon]
            mut_aa = CODON_TO_AA[mut_codon]
            if wt_aa != mut_aa:
                aa_muts[r] = f"{wt_aa}{r}{mut_aa}"

        return ' '.join([mut_str for r, mut_str in sorted(aa_muts.items())])


    def _ntToCodonMuts(self, nt_mut_str):
        """Converts string of nucleotide mutations to codon mutations.

        Args:
            `nt_mut_str` (str)
                Nucleotide mutations, delimited by a space and in
                1, 2, ... numbering.

        Returns:
            String with codon mutations in 1, 2, ... numbering of
            codon sites.

        >>> geneseq = 'ATGGGATGA'
        >>> with tempfile.NamedTemporaryFile(mode='w') as f:
        ...     _ = f.write('library,barcode,substitutions,variant_call_support')
        ...     f.flush()
        ...     variants = CodonVariantTable(
        ...                 barcode_variant_file=f.name,
        ...                 geneseq=geneseq
        ...                 )
        >>> variants._ntToCodonMuts('A1G G4C A6T')
        'ATG1GTG GGA2CGT'
        >>> variants._ntToCodonMuts('A1G G4C G6T')
        Traceback (most recent call last):
        ...
        ValueError: nucleotide 6 should be A not G
        """
        mut_codons = collections.defaultdict(set)
        for mut in nt_mut_str.upper().split():
            m = re.match('^(?P<wt>[ATCG])(?P<i>\d+)(?P<mut>[ATCG])$', mut)
            if not m:
                raise ValueError(f"invalid mutation {mut}")
            wt_nt = m.group('wt')
            i = int(m.group('i'))
            mut_nt = m.group('mut')
            if wt_nt == mut_nt:
                raise ValueError(f"invalid mutation {mut}")
            if i > len(self.geneseq):
                raise ValueError(f"invalid nucleotide site {i}")
            if self.geneseq[i - 1] != wt_nt:
                raise ValueError(f"nucleotide {i} should be "
                                 f"{self.geneseq[i - 1]} not {wt_nt}")
            icodon = (i - 1) // 3 + 1
            i_nt = (i - 1) % 3
            assert self.codons[icodon][i_nt] == wt_nt
            if i_nt in mut_codons[icodon]:
                raise ValueError(f"duplicate mutations {i_nt} in {icodon}")
            mut_codons[icodon].add((i_nt, mut_nt))

        codon_mut_list = []
        for r, r_muts in sorted(mut_codons.items()):
            wt_codon = self.codons[r]
            mut_codon = list(wt_codon)
            for i, mut_nt in r_muts:
                mut_codon[i] = mut_nt
            codon_mut_list.append(f"{wt_codon}{r}{''.join(mut_codon)}")

        return ' '.join(codon_mut_list)



def simulateSampleCounts(*,
                         variants,
                         phenotype_func,
                         variant_error_rate,
                         pre_sample,
                         post_samples,
                         seed=None):
    """Simulate pre- and post-selection variant counts.

    Simulate variant counts for experiment where barcodes
    are sequenced pre- and post-selection.

    Args:
        `variants` (:class:`CodonVariantTable`)
           Holds variants used in simulation.
        `phenotype_func` (function)
            Takes a row from `variants.barcode_variant_df`
            and returns the phenotype. Typically this is
            calculated from the "aa_substitutions"
            or "codon_substitutions" column. The phenotype
            is a number >= 0, and represents the expected
            enrichment of the variant relative to wildtype
            post selection (values > 1 indicate beneficial).
        `variant_error_rate` (float)
            Rate at which variants in `variants` are
            mis-called. Provide the probability that a
            variant has a spuriously called (or missing)
            codon mutation: with this probability, each
            variant then has a random codon mutation added
            or taken away before being passed to
            `phenotype_func`. Designed to simulate the fact
            that the method used to link barcodes to variants
            is probably not perfect.
        `pre_sample` (pandas DataFrame or dict)
            The counts of each variant pre-selection. You can
            specify counts or have them simulated:

                - To specify counts, provide a data frame with
                  columns "library", "barcode", and "count"
                  where "count" is the pre-selection counts. 

                - To simulate, provide dict with keys "total_count"
                  and "uniformity". For each library, we simulate
                  pre-selection counts as a draw of "total_count"
                  counts from a multinomial parameterized by
                  pre-selection frequences drawn from a Dirichlet
                  distribution with a concentration parameter of
                  of "uniformity" (larger uniformity values mean
                  more even libraries; 5 is reasonable value).

        `post_samples` (dict)
            Nested dict indicating post-selection samples.
            Keyed by name of each post-selection sample, with
            value being another dict. The post-selection counts
            are drawn from a multinomial parameterized
            by the pre-selection frequencies (frequencies calculated from
            counts if `pre_sample` gives counts, Dirichlet frequencies
            if `pre_sample` is simulated). Add noise by the
            following parameters provided in the dict for each sample:

                - "total_count": total overall counts per library.

                - "noise": add additional noise to selection by
                  multiplying phenotype times a random variable
                  drawn from a normal distribution with mean 1
                  and this standard deviation (truncated at lower
                  end to zero). Set noise to 0 for no noise.

                - "bottleneck": put the pre-selection frequencies
                  through a bottleneck of this size, then re-calcuate
                  initial frequencies that selection acts upon.

        `seed` (None or int)
            If not `None`, random number seed set before executing
            function (sets `scipy.random.seed`).

    Returns:
        A pandas DataFrame with the following columns:

            - "library"

            - "barcode"

            - "sample"

            - "count"

        The first two columns indicate the library and barcode
        for each variant as in the `barcode_variant_df` attribute
        of `variants`, the "sample" and "count" columns give counts
        for "pre-selection" and each post-selection sample in
        `post_samples`.
    """
    if seed is not None:
        scipy.random.seed(seed)

    #-----------------------------------------
    # internal function
    def _add_variant_errors(codon_substitutions):
        """Add errors to variant according to `variant_error_rate`."""
        if scipy.random.random() < variant_error_rate:
            muts = codon_substitutions.split()
            if len(muts) == 0 or scipy.random.random() < 0.5:
                # add mutation
                mutatedsites = set(map(
                        int,
                        [re.match('^[ATGC]{3}(?P<site>\d+)[ATGC]{3}$',
                                  mut).group('site')
                         for mut in muts]
                        ))
                unmutatedsites = list(set(variants.sites) - mutatedsites)
                if not unmutatedsites:
                    raise RuntimeError("variant already has all mutations")
                errorsite = scipy.random.choice(unmutatedsites)
                wtcodon = variants.codons[errorsite]
                mutcodon = scipy.random.choice(list(set(CODONS) - {wtcodon}))
                muts.append(f'{wtcodon}{errorsite}{mutcodon}')
                return ' '.join(muts)
            else:
                # remove mutation 
                muts = muts.pop(scipy.random.randint(0, len(muts)))
                return muts
        else:
            return codon_substitutions
    #-----------------------------------------

    barcode_variant_df = (
        variants.barcode_variant_df
        [['library', 'barcode', 'codon_substitutions']]
        .assign(
            codon_substitutions=lambda x: x.codon_substitutions
                                .apply(_add_variant_errors),
            aa_substitutions=lambda x: x.codon_substitutions
                             .apply(CodonVariantTable.codonToAAMuts),
            phenotype=lambda x: x.apply(phenotype_func, axis=1)
            )
        [['library', 'barcode', 'phenotype']]
        )

    libraries = variants.libraries

    if isinstance(pre_sample, pandas.DataFrame):
        # pre-sample counts specified
        req_cols = ['library', 'barcode', 'count']
        if not set(req_cols).issubset(set(pre_sample.columns)):
            raise ValueError(f"pre_sample lacks cols {req_cols}:"
                             f"\n{pre_sample}")
        cols = ['library', 'barcode']
        if (pre_sample[cols].sort_values(cols).reset_index() !=
            barcode_variant_df[cols].sort_values(cols).reset_index()
            ).any():
            raise ValueError("pre_sample DataFrame lacks required "
                             "library and barcode columns")
        barcode_variant_df = (
                barcode_variant_df
                .merge(pre_sample[req_cols], on=['library', 'barcode'])
                .rename(columns={'count':'pre-selection'})
                )
        # "true" pre-selection freqs are just input counts
        nperlib = (barcode_variant_df
                   .groupby('library')
                   .rename('total_count')
                   .reset_index()
                   )
        barcode_variant_df = (
                barcode_variant_df
                .merge(nperlib, on='library')
                .assign(pre_freqs=lambda x: x['pre-selection'] /
                                            x.total_count)
                .drop(columns='total_count')
                )

    elif isinstance(pre_sample, dict):
        pre_req_keys = {'uniformity', 'total_count'}
        if set(pre_sample.keys()) != pre_req_keys:
            raise ValueError(f"pre_sample lacks required keys {pre_req_keys}")

        pre_df_list = []
        for lib in libraries:
            df = (
                barcode_variant_df
                .query('library == @lib')
                .assign(
                    pre_freq=lambda x: scipy.random.dirichlet(
                             pre_sample['uniformity'] *
                             scipy.ones(len(x))),
                    count=lambda x: scipy.random.multinomial(
                            pre_sample['total_count'], x.pre_freq),
                    sample='pre-selection'
                    )
                )
            pre_df_list.append(df)
        barcode_variant_df = pandas.concat(pre_df_list)

    else:
        raise ValueError("pre_sample not DataFrame / dict: "
                         f"{pre_sample}")

    cols = ['library', 'barcode', 'sample', 'count',
            'pre_freq', 'phenotype']
    assert set(barcode_variant_df.columns) == set(cols), (
            f"cols = {set(cols)}\nbarcode_variant_df.columns = "
            f"{set(barcode_variant_df.columns)}")
    for col in cols:
        if col in post_samples:
            raise ValueError(f"post_samples can't have key {col}; "
                             "choose another sample name")

    df_list = [barcode_variant_df[cols[ : 4]]]

    post_req_keys = {'bottleneck', 'noise', 'total_count'}
    for lib, (sample, sample_dict) in itertools.product(
            libraries, post_samples.items()):

        if set(sample_dict.keys()) != post_req_keys:
            raise ValueError(f"post_samples {sample} lacks keys {post_req_keys}")

        lib_df = (
            barcode_variant_df.query('library == @lib')
            .assign(
                sample=sample,
                # simulated pre-selection freqs after bottleneck
                bottleneck_freq=lambda x: scipy.random.multinomial(
                        sample_dict['bottleneck'], x.pre_freq) /
                        sample_dict['bottleneck'],
                # post-selection freqs with noise
                noise=scipy.clip(scipy.random.normal(1, sample_dict['noise']),
                                 0, None),
                post_freq_nonorm=lambda x:
                            x.bottleneck_freq * x.phenotype * x.noise,
                post_freq=lambda x: x.post_freq_nonorm /
                            x.post_freq_nonorm.sum(),
                # post-selection counts simulated from frequencies
                count=lambda x: scipy.random.multinomial(
                            sample_dict['total_count'], x.post_freq)
                )
            .rename(columns={'post_counts':sample})
            [['library', 'barcode', 'sample', 'count']]
            )

        df_list.append(lib_df)

    return pandas.concat(df_list)


def tidy_split(df, column, sep=' ', keep=False):
    """
    Split values of a column and expand so new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Taken from https://stackoverflow.com/a/39946744

    Args:
        df : pandas.DataFrame
            dataframe with the column to split and expand
        column : str
            the column to split and expand
        sep : str
            the string used to split the column's values
        keep : bool
            whether to retain the presplit value as it's own row

    Returns:
        pandas.DataFrame
            Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    df = df.dropna(subset=[column])
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df


if __name__ == '__main__':
    import doctest
    doctest.testmod()
