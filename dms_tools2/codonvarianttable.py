"""
=================
codonvarianttable
=================

This module defines :class:`CodonVariantTable` objects
for storing and handling codon variants of a gene. They are
useful when you are performing deep mutational scanning
on barcoded codon-mutant libraries of a gene.
"""

import re
import os
import math
import collections
import itertools
import tempfile
import random

import scipy
import pandas as pd
import Bio.SeqUtils.ProtParamData
import gpmap
import numpy as np

# use plotnine for plotting
from plotnine import *

from dms_tools2.plot import latexSciNot
from dms_tools2 import (CODON_TO_AA,
                        CODONS,
                        AAS_WITHSTOP,
                        AA_TO_CODONS,
                        NTS
                        )

#: `color-blind safe palette <http://bconnelly.net/2013/10/creating-colorblind-friendly-figures/>`_
CBPALETTE = ["#999999", "#E69F00", "#56B4E9", "#009E73",
             "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]


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
    >>> pd.set_option('display.max_columns', 15)
    >>> pd.set_option('display.width', 500)
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

    Here is the number of variants if we look at just single amino-acid
    mutants (and wildtype):

    >>> variants.n_variants_df(samples=None, variant_type='single', mut_type='aa')
             library             sample  count
    0          lib_1  barcoded variants      2
    1          lib_2  barcoded variants      1
    2  all libraries  barcoded variants      3

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

    >>> counts_lib1_input = pd.DataFrame(
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

    >>> counts_lib1_selected = pd.DataFrame(
    ...         {'barcode':['AAC', 'GAT'],
    ...          'count'  :[  513,  401]})
    >>> variants.addSampleCounts('lib_1', 'selected', counts_lib1_selected)

    As well as barcode counts for the same two samples
    ("input" and "selected") to our other library (library 2):

    >>> counts_lib2_input = pd.DataFrame(
    ...         {'barcode':['AAC', 'CAT'],
    ...          'count'  :[ 1253,  923]})
    >>> variants.addSampleCounts('lib_2', 'input', counts_lib2_input)
    >>> counts_lib2_selected = pd.DataFrame(
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
    ...     countfiles = variants.writeCodonCounts("single", outdir=tmpdir, include_all_libs=True)
    ...     lib1_input = pd.read_csv(f'{tmpdir}/lib_1_input_codoncounts.csv')
    ...     all_sel = pd.read_csv(f'{tmpdir}/all-libraries_selected_codoncounts.csv')

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
    ...     _ = variants.writeCodonCounts("all", outdir=tmpdir, include_all_libs=True)
    ...     lib1_input_all = pd.read_csv(f'{tmpdir}/lib_1_input_codoncounts.csv')
    ...     all_sel_all = pd.read_csv(f'{tmpdir}/all-libraries_selected_codoncounts.csv')
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

    Of course, the initialized variant table is **not** equal
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

    We can use :meth:`CodonVariantTable.func_scores` to compute
    the functional effects of mutations. We cannot use this method
    with default options as we have no wildtype counts (needed for
    normalization) for `lib_2` so we get an error:

    >>> variants.func_scores('input')
    Traceback (most recent call last):
    ...
    ValueError: no wildtype counts:
      library    sample  count
    0   lib_1     input    253
    1   lib_1  selected    513
    2   lib_2     input      0
    3   lib_2  selected      0

    However, we can use the method with the `permit_zero_wt` option:

    >>> variants.func_scores('input', permit_zero_wt=True)
      library pre_sample post_sample barcode  func_score  func_score_var  pre_count  post_count  pre_count_wt  post_count_wt  pseudocount codon_substitutions  n_codon_substitutions aa_substitutions  n_aa_substitutions
    0   lib_1      input    selected     AAC    0.000000        0.024528        253         513           253            513          0.5                                          0                                    0
    1   lib_1      input    selected     GAT   -2.474376        0.019337       1101         401           253            513          0.5             GGA2CGC                      1              G2R                   1
    2   lib_2      input    selected     AAC   -3.465198        8.345474       1253         113             0              0          0.5     ATG1AAG GGA2TGA                      2          M1K G2*                   2
    3   lib_2      input    selected     CAT    0.378452        8.329463        923        1200             0              0          0.5             GGA2GGC                      1                                    0
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
                if isinstance(val, pd.DataFrame):
                    if not val.equals(val2):
                        return False
                else:
                    if val != val2:
                        return False
            return True

    @classmethod
    def from_simulation(cls, *, geneseq, bclen, library_specs,
                        seed=1, variant_call_support=1):
        """Simulate :class:`CodonVariantTable` variants.

        Use this method to simulate the variants in a
        :class:`CodonVariantTable`. Note that this only
        simulates the variants, not counts for samples.
        To add those, then use :func:`simulateSampleCounts`
        and :meth:`CodonVariantTable.addSampleCounts`.

        Args:
            `geneseq` (str)
                Sequence of wildtype protein-coding gene.
            `bclen` (int)
                Length of the barcodes; must enable complexity
                at least 10-fold greater than max number of variants.
            `library_specs` (dict)
                Specifications for each simulated library. Keys
                are 'avgmuts' and 'nvariants', and values are
                average-codon mutations per variant and number of
                variants. Mutations per variant are Poisson distributed.
            `seed` (int or `None`)
                Random number seed or `None` to set no seed.
            `variant_call_support` (int or 2-tuple)
                If an integer, all variant call supports are set to this.
                If a 2-tuple, variant call support is drawn as a random
                integer uniformly from this range (inclusive).

        Returns:
            The simulated :class:`CodonVariantTable`.
        """
        if seed is not None:
            scipy.random.seed(seed)
            random.seed(seed)

        if len(library_specs) < 1:
            raise ValueError('empty `library_specs`')

        if isinstance(variant_call_support, int):
            variant_call_support = tuple([variant_call_support] * 2)

        if len(geneseq) % 3 != 0:
            raise ValueError('length of `geneseq` not multiple of 3')
        genelength = len(geneseq) // 3

        barcode_variant_dict = collections.defaultdict(list)
        for lib, specs_dict in library_specs.items():

            nvariants = specs_dict['nvariants']
            avgmuts = specs_dict['avgmuts']
            if 10 * nvariants > (len(NTS))**bclen:  # safety factor 10
                raise ValueError('barcode too short for nvariants')
            existing_barcodes = set([])

            for ivariant in range(nvariants):

                barcode = ''.join(random.choices(NTS, k=bclen))
                while barcode in existing_barcodes:
                    barcode = ''.join(random.choices(NTS, k=bclen))
                existing_barcodes.add(barcode)

                support = random.randint(*variant_call_support)

                # get mutations
                substitutions = []
                nmuts = scipy.random.poisson(avgmuts)
                for icodon in random.sample(range(1, genelength + 1), nmuts):
                    wtcodon = geneseq[3 * (icodon - 1): 3 * icodon]
                    mutcodon = random.choice([c for c in CODONS
                                              if c != wtcodon])
                    for i_nt, (wt_nt, mut_nt) in enumerate(zip(wtcodon,
                                                               mutcodon)):
                        if wt_nt != mut_nt:
                            igene = 3 * (icodon - 1) + i_nt + 1
                            substitutions.append(f'{wt_nt}{igene}{mut_nt}')
                substitutions = ' '.join(substitutions)

                barcode_variant_dict['barcode'].append(barcode)
                barcode_variant_dict['substitutions'].append(substitutions)
                barcode_variant_dict['library'].append(lib)
                barcode_variant_dict['variant_call_support'].append(support)

        barcode_variants = pd.DataFrame(barcode_variant_dict)

        with tempfile.NamedTemporaryFile(mode='w') as f:
            barcode_variants.to_csv(f, index=False)
            f.flush()
            cvt = cls(barcode_variant_file=f.name, geneseq=geneseq)

        return cvt

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
        df = pd.read_csv(variant_count_df_file)

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

        df = pd.read_csv(barcode_variant_file)
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
            codonSubsFunc = self._sortCodonMuts
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
                                pd.Categorical(
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
                'nonsynonymous':CBPALETTE[1],
                'synonymous':CBPALETTE[2],
                'stop':CBPALETTE[3]
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
            self.variant_count_df = pd.concat(
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
                                pd.Categorical(
                                    x.library,
                                    categories=self.libraries,
                                    ordered=True
                                    ),
                        sample=lambda x:
                               pd.Categorical(
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


    def func_scores(self, preselection, *,
                pseudocount=0.5, by="barcode",
                combine_libs=False, syn_as_wt=False, logbase=2,
                permit_zero_wt=False, permit_self_comp=False):
        """Get data frame with functional scores for variants.

        The functional score is calculated from the change in counts
        for a variant pre- and post-selection using the formula in
        `Otwinoski et al (2018) <https://doi.org/10.1073/pnas.1804015115>`_.
        Specifically, if :math:`n^v_{pre}` and :math:`n^v_{post}` are
        the counts of variant :math:`v` pre- and post-selection, and
        if :math:`n^{wt}_{pre}` and :math:`n^{wt}_{post}` are
        the summed counts of **all** wildtype variants pre- and post-
        selection, then the functional score of the variant is
        :math:`f_v = \log_b\left(\\frac{n^v_{post} / n^{wt}_{post}}{n^v_{pre} / n^{wt}_{pre}}\\right)`
        and the variance due to Poisson sampling statistics is
        :math:`\\frac{1}{\left(\ln b\\right)^2}\left(\sigma^2_v = \\frac{1}{n^v_{post}} + \\frac{1}{n^{wt}_{post}} + \\frac{1}{n^v_{pre}} + \\frac{1}{n^{wt}_{pre}}\\right)`
        where :math:`b` is logarithm base (see `logbase` argument).
        For both calculations, a pseudocount (see `pseudocount` argument)
        is added to each count first. The wildtype counts are computed
        across all **fully wildtype** variants (see `syn_as_wt` for
        how this is defined).

        Args:
            `preselection` (str or dict)
                The pre-selection sample. If the same for all post-
                selection samples, then provide the name of the
                pre-selection sample. If it differs among post-selection
                samples, then provide a dict keyed by each post-selection
                sample with the pre-selection sample being the value.
            `pseudocount` (float)
                The pseudocount added to each count.
            `by` (str)
                Compute effects for each "barcode", set of
                "aa_substitutions", or set of "codon_substitutions".
                In the latter two cases, all barcodes with each
                set of substitutions are combined (see `combine_libs`).
                Note that if you use "aa_substitutions" then it may
                be more sensible to set `syn_as_wt` to `True`.
            `syn_as_wt` (bool)
                In the formula for the functional scores, consider
                variants with only synonymous mutations when determining
                wildtype counts? If `False`, only variants with **no**
                mutations of any type contribute to wildtype counts.
            `combine_libs` (bool)
                If `by` is "aa_substitutions" or "codon_substitutions",
                do we combine across libraries as well as barcodes?
            `logbase` (float)
                Base for logarithm when calculating functional score.
            `permit_zero_wt` (bool)
                If the wildtype counts are zero for any sample, do
                we raise an error or permit the calculation to proceed
                just using the pseudocount?
            `permit_self_comp` (bool)
                Permit comparisons where the same sample is used as
                the pre- and post-selection?

        Returns:
            A pandas Data Frame with the following columns:
              - "library": the library ("all libraries" if `combine_libs`)
              - "pre_sample": the pre-selection sample
              - "post_sample": the post-selection sample
              - the value corresponding to the grouping we are using
                to compute effects (the value of `by`)
              - "func_score": the functional score
              - "func_score_var": variance on the functional score
              - "pre_count": pre-selection counts
              - "post_count: post-selection counts
              - "pre_count_wt": pre-selection counts for all wildtype
              - "post_count_wt": post-selection counts for all wildtype
              - "pseudocount": the pseudocount value
              - as many of "aa_substitutions", "n_aa_substitutions",
                "codon_substitutions", and "n_codon_substitutions"
                as can be retained given the value of `by`.
        """
        ordered_samples = self.variant_count_df['sample'].unique()
        if isinstance(preselection, str):
            # make `preselection` into dict
            preselection = {s:preselection for s in ordered_samples
                            if s != preselection or permit_self_comp}
        elif not isinstance(preselection, dict):
            raise ValueError('`preselection` not str or dict')
        if not permit_self_comp:
            if any(pre == post for pre, post in preselection.items()):
                raise ValueError('`permit_self_comp` is False but there'
                                 ' are identical pre and post samples')

        # all samples of interest
        samples = set(preselection.keys()).union(set(preselection.values()))
        if not samples.issubset(set(ordered_samples)):
            extra_samples = samples - set(ordered_samples)
            raise ValueError(f"invalid samples: {extra_samples}")

        # get data frame with samples of interest
        if self.variant_count_df is None:
            raise ValueError('no sample variant counts have been added')
        df = (self.variant_count_df
              .query('sample in @samples')
              )

        if combine_libs and (len(self.libraries) > 1):
            if any(s not in self.samples(lib) for lib in self.libraries
                                              for s in samples):
                raise ValueError('cannot use `combine_libs`, not every '
                                 f"library has every sample: {samples}")
            df = (self.addMergedLibraries(df)
                  .query('library == "all libraries"')
                  )

        # get wildtype counts for each sample and library
        if syn_as_wt:
            wt_col = 'n_aa_substitutions'
        else:
            wt_col = 'n_codon_substitutions'
        wt_counts = (
                df
                .assign(count=lambda x: x['count'] * (0 == x[wt_col])
                                        .astype('int'))
                .groupby(['library', 'sample'])
                .aggregate({'count':'sum'})
                .reset_index()
                )
        if (wt_counts['count'] <= 0).any() and not permit_zero_wt:
            raise ValueError(f"no wildtype counts:\n{wt_counts}")

        # sum counts in groups specified by `by`
        group_cols = ['codon_substitutions', 'n_codon_substitutions',
                      'aa_substitutions', 'n_aa_substitutions']
        if by in {'aa_substitutions', 'codon_substitutions'}:
            group_cols = group_cols[group_cols.index(by) + 1 : ]
        elif by != 'barcode':
            raise ValueError(f"invalid `by` of {by}")
        df = (df
              .groupby(['library', 'sample', by] + group_cols)
              .aggregate({'count':'sum'})
              .reset_index()
              )

        # get data frame with pre- and post-selection samples / counts
        df_func_scores = []
        for post_sample in ordered_samples:
            if post_sample not in preselection:
                continue
            pre_sample = preselection[post_sample]
            sample_dfs = []
            for stype, s in [('pre', pre_sample), ('post', post_sample)]:
                sample_dfs.append(
                        df
                        .query('sample == @s')
                        .rename(columns={'count':f"{stype}_count"})
                        .merge(wt_counts
                               .rename(columns={'count':f"{stype}_count_wt"}),
                               how='inner', validate='many_to_one'
                               )
                        .rename(columns={'sample':f"{stype}_sample"})
                        )
            df_func_scores.append(
                    pd.merge(sample_dfs[0], sample_dfs[1],
                             how='inner', validate='1:1')
                    )
        df_func_scores = pd.concat(df_func_scores,
                                   ignore_index=True, sort=False)

        # check pseudocount
        if pseudocount < 0:
            raise ValueError(f"`pseudocount` is < 0: {pseudocount}")
        elif (pseudocount == 0) and any((df_func_scores[c] <= 0).any()
                                    for c in ['pre_count', 'post_count',
                                    'pre_count_wt', 'post_count_wt']):
            raise ValueError('some counts are zero, you must use '
                             '`pseudocount` > 0')

        # calculate functional score and variance
        df_func_scores = (
                df_func_scores
                .assign(
                    pseudocount=pseudocount,
                    func_score=lambda x:
                               scipy.log(
                                    ((x.post_count + x.pseudocount) /
                                     (x.post_count_wt + x.pseudocount)) /
                                    ((x.pre_count + x.pseudocount) /
                                     (x.pre_count_wt + x.pseudocount))
                               ) / scipy.log(logbase),
                    func_score_var=lambda x:
                               (1 / (x.post_count + x.pseudocount) +
                                1 / (x.post_count_wt + x.pseudocount) +
                                1 / (x.pre_count + x.pseudocount) +
                                1 / (x.pre_count_wt + x.pseudocount)
                               ) / (scipy.log(logbase)**2)
                    )
                # set column order in data frame
                [['library', 'pre_sample', 'post_sample', by,
                  'func_score', 'func_score_var', 'pre_count',
                  'post_count', 'pre_count_wt', 'post_count_wt',
                  'pseudocount'] + group_cols]
                )

        return df_func_scores


    def n_variants_df(self, *, libraries='all', samples='all',
                      min_support=1, variant_type='all',
                      mut_type=None):
        """Number of variants per library / sample.

        Args:
            `variant_type` ("single" or "all")
                Include just variants with <= 1 `mut_type` mutation.
            `mut_type` ("aa", "codon", or `None`)
                If `variant_type` is single, set to "codon" or "aa"
                to indicate how to filter for single mutants.
            All other args:
                Same meaning as for
                :class:`CodonVariantTable.plotNumMutsHistogram`.

        Returns:
            DataFrame giving number of variants per library /
            sample.
        """
        df, nlibraries, nsamples = self._getPlotData(libraries,
                                                     samples,
                                                     min_support)

        if variant_type == 'single':
            if mut_type in {'aa', 'codon'}:
                df = df.query(f"n_{mut_type}_substitutions <= 1")
            else:
                raise ValueError('`mut_type` must be "aa" or "single"')
        elif variant_type != 'all':
            raise ValueError(f"invalid `variant_type` {variant_type}")

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
        all_muts = pd.concat([
                    pd.DataFrame({'mutation':mut_list,
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
                         pd.Categorical(
                          x['library'],
                          librarylist,
                          ordered=True),
                sample=lambda x:
                         pd.Categorical(
                          x['sample'],
                          samplelist,
                          ordered=True),
                mutation_type=lambda x:
                         pd.Categorical(
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
            A `plotnine <https://plotnine.readthedocs.io>`_
            plot; can be displayed in a Jupyter notebook with `p.draw()`.
        """

        df = self.mutCounts(variant_type, mut_type, samples=samples,
                            libraries=libraries, min_support=min_support)

        n_variants = (self.n_variants_df(libraries=libraries,
                                         samples=samples,
                                         min_support=min_support,
                                         variant_type=variant_type,
                                         mut_type=mut_type)
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
                        pd.Categorical(
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
            A `plotnine <https://plotnine.readthedocs.io>`_
            plot; can be displayed in a Jupyter notebook with `p.draw()`.
        """

        df = self.mutCounts(variant_type, mut_type, samples=samples,
                            libraries=libraries, min_support=min_support)

        n_variants = (self.n_variants_df(libraries=libraries,
                                         samples=samples,
                                         min_support=min_support,
                                         variant_type=variant_type,
                                         mut_type=mut_type)
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


    def plotCumulVariantCounts(self, *, variant_type='all',
            libraries='all', samples='all', plotfile=None,
            orientation='h', widthscale=1, heightscale=1,
            min_support=1, mut_type='aa', tot_variants_hline=True):
        """Plot number variants with >= that each number of counts.

        Args:
            `variant_type` ("single" or "all")
                Include just variants with <= one mutation of type,
                `mut_type` or all mutants?
            `tot_variants_hline` (bool)
                Include a dotted horizontal line indicating the total
                number of variants (useful since ones not observed
                at all will not apparent on plot).
            Other args:
                Same meaning as for
                :meth:`CodonVariantTable.plotNumMutsHistogram`.

        Returns:
            `plotnine <https://plotnine.readthedocs.io>`_
            plot; can be displayed in a Jupyter notebook with `p.draw()`.
        """
        df, nlibraries, nsamples = self._getPlotData(libraries,
                                                     samples,
                                                     min_support)

        if variant_type == 'single':
            if mut_type == 'aa':
                mutstr = 'amino acid'
            elif mut_type == 'codon':
                mutstr = mut_type
            else:
                raise ValueError(f"invalid `mut_type` {mut_type}")
            ylabel = f"single {mutstr} variants with >= this many counts"
            df = df.query(f"n_{mut_type}_substitutions <= 1")
        elif variant_type == 'all':
            ylabel = 'variants with >= this many counts'
        else:
            raise ValueError(f"invalid `variant_type` {variant_type}")

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

        df = (getCumulVariantsByCount(df, group_cols=['library', 'sample'])
              .query('count > 0')
              )

        p = (ggplot(df, aes('count', 'nvariants')) +
             geom_step() +
             facet_grid(facet_str) +
             xlab('number of counts') +
             ylab(ylabel) +
             scale_x_log10(labels=latexSciNot) +
             scale_y_continuous(labels=latexSciNot) +
             theme(figure_size=(width, height))
             )

        if tot_variants_hline:
            p = p + geom_hline(aes(yintercept='total_variants'),
                               linetype='dashed', color=CBPALETTE[1])

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
            A `plotnine <https://plotnine.readthedocs.io>`_
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
            A `plotnine <https://plotnine.readthedocs.io>`_
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
                                pd.Categorical(
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
            `plotnine <https://plotnine.readthedocs.io>`_
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

        df[mut_col] = scipy.clip(df[mut_col], None, max_muts)

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


    def writeCodonCounts(self, single_or_all, *,
                         outdir=None, include_all_libs=False):
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
            `include_all_libs` (bool)
                Include data for a library (named "all-libraries")
                that has the summmed data for all individual libraries
                if there are multiple libraries?

        Returns:
            Pandas data frame with columns "library", "sample",
            and "countfile". The "countfile" columns gives name of the
            createad CSV file, ``<library>_<sample>_codoncounts.csv``.
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

        if include_all_libs:
            df = self.addMergedLibraries(self.variant_count_df,
                                         all_lib='all-libraries')
        else:
            df = self.variant_count_df

        countfiles = []
        liblist = []
        samplelist = []

        for lib, sample in itertools.product(
                            df['library'].unique().tolist(),
                            df['sample'].unique().tolist()
                            ):

            i_df = df.query('library == @lib & sample == @sample')
            if len(i_df) == 0:
                continue # no data for this library and sample

            countfile = os.path.join(outdir,
                            f'{lib}_{sample}_codoncounts.csv')
            countfiles.append(countfile)
            liblist.append(lib)
            samplelist.append(sample)

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

            counts_df = pd.DataFrame(collections.OrderedDict(
                         [('site', self.sites),
                          ('wildtype', [self.codons[r] for r in self.sites])] +
                         [(codon, codoncounts[codon]) for codon in CODONS]
                         ))
            counts_df.to_csv(countfile, index=False)

        assert all(map(os.path.isfile, countfiles))

        return pd.DataFrame({'library':liblist,
                             'sample':samplelist,
                             'countfile':countfiles})

    @staticmethod
    def classifyVariants(df,
                         *,
                         variant_class_col='variant_class',
                         max_aa=2):
        """Classifies codon variants in `df`.

        Args:
            `df` (pandas DataFrame)
                Must have columns named 'aa_substitutions',
                'n_aa_substitutions', and 'n_codon_substitutions'.
                For instance, a data frame of this type can
                be obtained via the `variant_count_df` or
                `barcode_variant_df` of a :class:`CodonVariantTable`,
                or via :meth:`CodonVariantTable.func_scores`.
            `variant_class_col` (str)
                Name of column added to `df` that contains
                variant classification. Overwritten if already exists.
            `max_aa` (int)
                When classifying variants, group all with >=
                this many amino-acid mutations.

        Returns:
            A copy of `df` with the column specified by
            `variant_class_col` classifying variants as:

              - 'wildtype': no codon mutations

              - 'synonymous': only synonymous codon mutations

              - 'stop': at least one stop-codon mutation

              - '{n_aa} nonsynonymous' where `n_aa` is the
                number of amino-acid mutations, or is '>{max_aa}'
                if there are more than `max_aa` such mutations.

        >>> df = pd.DataFrame.from_records(
        ...         [('AAA', '', 0, 0),
        ...          ('AAG', '', 0, 1),
        ...          ('ATA', 'M1* G5K', 2, 3),
        ...          ('GAA', 'G5H', 1, 2),
        ...          ('CTT', 'M1C G5C', 2, 3),
        ...          ('CTT', 'M1A L3T G5C', 3, 3),
        ...          ],
        ...         columns=['barcode', 'aa_substitutions',
        ...                  'n_aa_substitutions', 'n_codon_substitutions']
        ...         )
        >>> CodonVariantTable.classifyVariants(df)
          barcode aa_substitutions  n_aa_substitutions  n_codon_substitutions      variant_class
        0     AAA                                    0                      0           wildtype
        1     AAG                                    0                      1         synonymous
        2     ATA          M1* G5K                   2                      3               stop
        3     GAA              G5H                   1                      2    1 nonsynonymous
        4     CTT          M1C G5C                   2                      3  >=2 nonsynonymous
        5     CTT      M1A L3T G5C                   3                      3  >=2 nonsynonymous
        """
        req_cols = ['aa_substitutions', 'n_aa_substitutions',
                    'n_codon_substitutions']
        if not (set(req_cols) <= set(df.columns)):
            raise ValueError(f"`df` does not have columns {req_cols}")

        def _classify_func(row):
            if row['n_codon_substitutions'] == 0:
                return 'wildtype'
            elif row['n_aa_substitutions'] == 0:
                return 'synonymous'
            elif '*' in row['aa_substitutions']:
                return 'stop'
            elif row['n_aa_substitutions'] < max_aa:
                return f"{row['n_aa_substitutions']} nonsynonymous"
            else:
                return f">={max_aa} nonsynonymous"

        return df.assign(**{variant_class_col: lambda x:
                                               x.apply(_classify_func,
                                                       axis=1)
                            })


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

        df = (pd.concat([df,
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
                              pd.Categorical(
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
            df = self.variant_count_df.query('sample in @samples')
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
            df = df.query('library in @libraries')
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


    def _sortCodonMuts(self, mut_str):
        """Sort space-delimited codon mutations and make uppercase.

        >>> geneseq = 'ATGGGATGA'
        >>> with tempfile.NamedTemporaryFile(mode='w') as f:
        ...     _ = f.write('library,barcode,substitutions,variant_call_support')
        ...     f.flush()
        ...     variants = CodonVariantTable(
        ...                 barcode_variant_file=f.name,
        ...                 geneseq=geneseq
        ...                 )
        >>> variants._sortCodonMuts('GGA2CGT ATG1GTG')
        'ATG1GTG GGA2CGT'
        """
        muts = {}
        for mut in mut_str.upper().split():
            m = re.match('^(?P<wt>[ATCG]{3})(?P<r>\d+)(?P<mut>[ACTG]{3})$',
                         mut)
            if not m:
                raise ValueError(f"invalid codon mutation {mut}")
            wt_codon = m.group('wt')
            r = int(m.group('r'))
            mut_codon = m.group('mut')
            if wt_codon == mut_codon:
                raise ValueError(f"invalid codon mutation {mut}")
            if r not in self.sites:
                raise ValueError(f"invalid site in codon mutation {mut}")
            if wt_codon != self.codons[r]:
                raise ValueError(f"invalid wt in codon mutation {mut}")
            if r in muts:
                raise ValueError(f"duplicate mutation at codon {mut}")
            muts[r] = mut
        return ' '.join(mut for r, mut in sorted(muts.items()))


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
        >>> variants._ntToCodonMuts('G4C A6T A1G')
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
            if i > len(self.geneseq) or i < 1:
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
                         pre_sample_name='pre-selection',
                         seed=1):
    """Simulate pre- and post-selection variant counts.

    Simulate variant counts for experiment where barcodes
    are sequenced pre- and post-selection.

    Args:
        `variants` (:class:`CodonVariantTable`)
           Holds variants used in simulation.
        `phenotype_func` (function)
            Takes a row from `variants.barcode_variant_df`
            and returns the phenotype. Typically this is
            calculated from the "aa_substitutions" or
            "codon_substitutions" column. The phenotype is a
            number >= 0, and represents the expected enrichment of
            the variant relative to wildtype post-selection (values
            > 1 indicate beneficial). For instance, you could pass
            :meth:`PhenotypeSimulator.observedPhenotype`.
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
                  initial frequencies that selection acts upon. Set
                  to `None` for no bottleneck.

        `pre_sample_name` (str)
            Name used for the pre-selection sample.
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
        for each sample.
    """
    if seed is not None:
        scipy.random.seed(seed)

    if pre_sample_name in post_samples:
        raise ValueError('`pre_sample_name` is in `post_samples`')

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
                unmutatedsites = [r for r in variants.sites
                                  if r not in mutatedsites]
                if not unmutatedsites:
                    raise RuntimeError("variant already has all mutations")
                errorsite = scipy.random.choice(unmutatedsites)
                wtcodon = variants.codons[errorsite]
                mutcodon = scipy.random.choice([c for c in CODONS
                                                if c != wtcodon])
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

    if isinstance(pre_sample, pd.DataFrame):
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
                .rename(columns={'count':pre_sample_name})
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
                .assign(pre_freqs=lambda x: x[pre_sample_name] /
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
                    sample=pre_sample_name
                    )
                )
            pre_df_list.append(df)
        barcode_variant_df = pd.concat(pre_df_list)

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

    def _bottleneck_freqs(pre_freq, bottleneck):
        if bottleneck is None:
            return pre_freq
        else:
            return scipy.random.multinomial(bottleneck, pre_freq) / bottleneck

    post_req_keys = {'bottleneck', 'noise', 'total_count'}
    for lib, (sample, sample_dict) in itertools.product(
            libraries, sorted(post_samples.items())):

        if set(sample_dict.keys()) != post_req_keys:
            raise ValueError(f"post_samples {sample} lacks keys {post_req_keys}")

        lib_df = (
            barcode_variant_df.query('library == @lib')
            .assign(
                sample=sample,
                # simulated pre-selection freqs after bottleneck
                bottleneck_freq=lambda x:
                                _bottleneck_freqs(x.pre_freq,
                                                  sample_dict['bottleneck']),
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

    return pd.concat(df_list)


class PhenotypeSimulator:
    """Simulates phenotypes of variants under plausible model.

    Mutational effects on latent phenotype are simulated to follow
    compound normal distribution; latent phenotype maps to observed
    phenotype via sigmoid. This distinction between latent and
    observed phenotype parallel the "global epistasis" models of
    `Otwinoski et al <https://doi.org/10.1073/pnas.1804015115>`_ and
    `Sailer and Harms <http://www.genetics.org/content/205/3/1079>`_.

    Initialize with a codon sequence and random number seed.
    The effects of mutations on the latent phenotype are then drawn
    from a compound normal distribution biased to negative values.

    To calculate latent and observed phenotypes, pass rows of a
    :meth:`CodonVariantTable.barcode_variant_df` to
    :meth:`PhenotypeSimulator.latentPhenotype` or
    :meth:`PhenotypeSimulator.observedPhenotype`.

    Args:
        `geneseq` (str)
            Codon sequence of wild-type gene.
        `seed` (int)
            Random number seed.
        `wt_latent` (float)
            Latent phenotype of wildtype.
        `norm_weights` (list of tuples)
            Specifies compound normal distribution of mutational
            effects on latent phenotype. Each tuple is
            `(weight, mean, sd)`, giving weight, mean, and standard
            deviation of each Gaussian in compound normal.
        `stop_effect` (float)
            Effect of stop codon at any position.

    Attributes:
        `wt_latent` (float)
            Value for wildtype latent phenotype.
        `muteffects` (dict)
            Effects on latent phenotype of each amino-acid mutation.
    """

    def __init__(self, geneseq, *, seed=1, wt_latent=1,
                 norm_weights=[(0.4, -0.5, 1), (0.6, -5, 2.5)],
                 stop_effect=-10):
        """See main class docstring for how to initialize."""
        self.wt_latent = wt_latent

        # simulate muteffects from compound normal distribution
        self.muteffects = {}
        scipy.random.seed(seed)
        weights, means, sds = zip(*norm_weights)
        cumweights = scipy.cumsum(weights)
        for icodon in range(len(geneseq) // 3):
            wt_aa = CODON_TO_AA[geneseq[3 * icodon : 3 * icodon + 3]]
            for mut_aa in AAS_WITHSTOP:
                if mut_aa != wt_aa:
                    if mut_aa == '*':
                        muteffect = stop_effect
                    else:
                        # choose Gaussian from compound normal
                        i = scipy.argmin(cumweights < scipy.random.rand())
                        # draw mutational effect from chosen Gaussian
                        muteffect = scipy.random.normal(means[i], sds[i])
                    self.muteffects[f'{wt_aa}{icodon + 1}{mut_aa}'] = muteffect

    def latentPhenotype(self, v):
        """Returns latent phenotype of a variant.

        Args:
            `v` (dict or row of pandas DataFrame)
                Must have key 'aa_substitutions' that gives
                space de-limited list of amino-acid mutations.
        """
        return self.wt_latent + sum([self.muteffects[m] for m in
                                     v['aa_substitutions'].split()])

    def observedPhenotype(self, v):
        """Like `latentPhenotype` but returns observed phenotype."""
        return self.latentToObservedPhenotype(self.latentPhenotype(v))

    @staticmethod
    def latentToObservedPhenotype(latent):
        """Returns observed phenotype from latent phenotype."""
        return 1 / (1 + math.exp(-latent - 3))

    def plotLatentVersusObservedPhenotype(self, *,
            latent_min=-15, latent_max=5, npoints=200):
        """Plots observed phenotype as function of latent phenotype.

        Plot includes a vertical line at wildtype latent phenotype.

        Args:
            `latent_min` (float)
                Smallest value of latent phenotype on plot.
            `latent_max` (float)
                Largest value of latent phenotype on plot.
            `npoints` (int)
                Plot a line fit to this many points.

        Returns:
            A `plotnine <https://plotnine.readthedocs.io>`_
            plot; can be displayed in a Jupyter notebook with `p.draw()`.
        """
        latent = scipy.linspace(latent_min, latent_max, npoints)
        p = (ggplot(pd.DataFrame(dict(latent=latent))
                       .assign(observed=lambda x: x.latent.apply(
                                        self.latentToObservedPhenotype)),
                    aes('latent', 'observed')
                    ) +
             geom_line() +
             geom_vline(xintercept=self.wt_latent, color=CBPALETTE[1],
                        linetype='dashed') +
            theme(figure_size=(3.5, 2.5)) +
            xlab('latent phenotype') +
            ylab('observed phenotype')
            )
        return p

    def plotMutsHistogram(self, latent_or_observed, *,
            mutant_order=1, bins=30):
        """Plots distribution of phenotype for all mutants.

        Plot includes a vertical line at wildtype phenotype.

        Args:
            `latent_or_observed` ("latent" or "observed")
                Which type of phenotype to plot.
            `mutant_order` (int)
                Plot mutations of this order. Currently only works
                for 1 (single mutants).
            `bins` (int)
                Number of bins in histogram.

        Returns:
            A `plotnine <https://plotnine.readthedocs.io>`_
            plot; can be displayed in a Jupyter notebook with `p.draw()`.
        """
        if mutant_order != 1:
            raise ValueError('only implemented for `mutant_order` of 1')
        if latent_or_observed == 'latent':
            phenoFunc = self.latentPhenotype
        elif latent_or_observed == 'observed':
            phenoFunc = self.observedPhenotype
        else:
            raise ValueError('invalid value of `latent_or_observed`')
        phenotypes = [phenoFunc({'aa_substitutions':m}) for m in
                      self.muteffects.keys()]
        p = (ggplot(pd.DataFrame({'phenotype':phenotypes}),
                    aes('phenotype')) +
             geom_histogram(bins=bins) +
             theme(figure_size=(3.5, 2.5)) +
             ylab(f"number of {mutant_order}-mutants") +
             xlab(f"{latent_or_observed} phenotype") +
             geom_vline(xintercept=phenoFunc({'aa_substitutions':''}),
                        color=CBPALETTE[1], linetype='dashed')
             )
        return p



def tidy_split(df, column, sep=' ', keep=False):
    """
    Split values of a column and expand so new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Taken from https://stackoverflow.com/a/39946744

    Args:
        df : pandas DataFrame
            dataframe with the column to split and expand
        column : str
            the column to split and expand
        sep : str
            the string used to split the column's values
        keep : bool
            whether to retain the presplit value as it's own row

    Returns:
        pandas DataFrame
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


def rarefyBarcodes(barcodecounts, *,
                   barcodecol='barcode', countcol='count',
                   maxpoints=1e5, logspace=True):
    """Rarefaction curve of barcode observations.

    Uses analytical formula for rarefaction defined
    `here <https://en.wikipedia.org/wiki/Rarefaction_(ecology)#Derivation>`_.

    Args:
        `barcodecounts` (pandas DataFrame)
            Has columns with names matching `barcodecol` and `countcol`.
        `barcodecol` (str)
            Name of column listing all unique barcodes.
        `countcol` (str)
            Name of column with observed counts of each barcode.
        `maxpoints` (int)
            Only calculate rarefaction curve at this many points.
            The benefit is that it can become very costly to calculate
            the curve for many points.
        `logspace` (bool)
            Logarithmically space the `maxpoint` points. If False,
            space them linearly.

    Returns:
        A pandas DataFrame with the columns `ncounts` and
        `nbarcodes`, giving the number of unique barcodes
        observed for each total number observed counts.

    Here is a small example:

    >>> barcodecounts = pd.DataFrame({'barcode':['A', 'G', 'C', 'T'],
    ...                               'count':  [  4,   2,   1,   0]})
    >>> rarefaction_curve = rarefyBarcodes(barcodecounts)
    >>> rarefaction_curve
       ncounts  nbarcodes
    0        1   1.000000
    1        2   1.666667
    2        3   2.114286
    3        4   2.428571
    4        5   2.666667
    5        6   2.857143
    6        7   3.000000

    Verify that the result from this example matches what is
    obtained by random sampling:

    >>> random.seed(1)
    >>> barcodelist = []
    >>> for tup in barcodecounts.itertuples(index=False):
    ...     barcodelist += [tup.barcode] * tup.count
    >>> nrand = 10000
    >>> ncounts = list(range(1, barcodecounts['count'].sum() + 1))
    >>> nbarcodes = []
    >>> for ncount in ncounts:
    ...     nbarcodes.append(sum(len(set(random.sample(barcodelist, ncount)))
    ...                      for _ in range(nrand)) / nrand)
    >>> sim_rarefaction_curve = pd.DataFrame(dict(ncounts=ncounts,
    ...                                           nbarcodes=nbarcodes))
    >>> scipy.allclose(rarefaction_curve, sim_rarefaction_curve, atol=1e-2)
    True
    """
    if len(barcodecounts) != len(barcodecounts[barcodecol].unique()):
        raise ValueError('non-unique barcodes in `barcodecounts`')

    # follow nomenclature at
    # https://en.wikipedia.org/wiki/Rarefaction_(ecology)#Derivation
    Ni = barcodecounts.set_index(barcodecol)[countcol].to_dict()
    N = sum(Ni.values())
    K = len(barcodecounts)
    Mj = collections.Counter(Ni.values())

    Nk, num = map(scipy.array, zip(*Mj.items()))

    # use simplification that (N - Ni)Cr(n) / (N)Cr(n) =
    # [(N - Ni)! * (N - n)!] / [N! * (N - Ni - n)!]
    #
    # Also use fact that gamma(x + 1) = x!
    nbarcodes = []
    lnFactorial_N = scipy.special.gammaln(N + 1)
    if logspace and N > maxpoints:
        ncounts = list(scipy.unique(scipy.logspace(
                       math.log10(1), math.log10(N),
                       num=min(N, maxpoints)).astype('int')))
    else:
        ncounts = list(scipy.unique(scipy.linspace(
                       1, N, num=min(N, maxpoints)).astype('int')))
    for n in ncounts:
        lnFactorial_N_minus_n = scipy.special.gammaln(N - n + 1)
        i = scipy.nonzero(N - Nk - n >= 0) # indices where this is true
        nbarcodes.append(
                K - (num[i] * scipy.exp(
                            scipy.special.gammaln(N - Nk[i] + 1) +
                            lnFactorial_N_minus_n -
                            lnFactorial_N -
                            scipy.special.gammaln(N - Nk[i] - n + 1))
                    ).sum()
                )
    return pd.DataFrame(dict(ncounts=ncounts, nbarcodes=nbarcodes))


def getCumulVariantsByCount(df, *, group_cols=None,
                            group_cols_as_str=False):
    """Get number of variants with observed >= each number of times.

    Args:
        `df` (pandas DataFrame)
            A row for each variant, and a column named "count" giving
            how many times variant observed. Can have other columns.
        `group_cols` (None or list)
            Group by these columns and analyze each group separately.
        `group_cols_as_str` (bool)
            Explicitly convert any `group_cols` columns to str type.
            For some reason, this is needed if you are calling in ``R``
            using `reticulate <https://rstudio.github.io/reticulate/>`_.

    Returns:
        A pandas Data Frame with columns named "count", "nvariants",
        and "total_variants" (plus any columns in `group_cols`).
        For each value of "count", "nvariants" gives the number of
        variants with >= that many counts. The "total_variants"
        column gives the total number of variants.

    Here is an example. First, create input Data Frame:

    >>> df = pd.DataFrame(dict(
    ...          sample= ['a', 'a', 'b', 'b', 'a', 'a'],
    ...          count=  [  9,   0,   1,   4,   3,   3]))

    Now apply function, first **not** grouping by sample. Note
    how the sample column is therefore simply ignored and dropped:

    >>> getCumulVariantsByCount(df)
       count  nvariants  total_variants
    0      9          1               6
    1      4          2               6
    2      3          4               6
    3      1          5               6
    4      0          6               6

    Now apply function, grouping by sample. Note how each sample
    is now analyzed separately:

    >>> getCumulVariantsByCount(df, group_cols=['sample'])
      sample  count  nvariants  total_variants
    0      a      9          1               4
    1      a      3          3               4
    2      a      0          4               4
    3      b      4          1               2
    4      b      1          2               2
    """
    if 'count' not in df.columns:
        raise ValueError('df does not have a column named "count"')

    if not group_cols:
        drop_group_cols = True
        group_cols=['dummy_col']
        df[group_cols[0]] = '_'
    else:
        drop_group_cols = False
        if isinstance(group_cols, str):
            group_cols = [group_cols]

    if not set(group_cols).issubset(set(df.columns)):
        raise ValueError(f"invalid `group_cols` {group_cols}")

    df = (
        df

        # get number of variants with each count
        .assign(nvariants=1)
        .groupby(group_cols + ['count'])
        .aggregate({'nvariants':'count'})
        .reset_index()
        .sort_values(group_cols + ['count'],
                     ascending=[True] * len(group_cols) + [False])
        .reset_index(drop=True)

        # make nvariants cumulative number with <= number of counts
        .assign(nvariants=lambda x: x.groupby(group_cols)
                                    ['nvariants']
                                    .cumsum())

        # add new column that is total number of variants
        .assign(total_variants=lambda x: x.groupby(group_cols)
                                         ['nvariants']
                                         .transform('max'))
        )

    if drop_group_cols:
        df = df.drop(group_cols, axis='columns')
    elif group_cols_as_str:
        for col in group_cols:
            df[col] = df[col].astype('str')

    return df


def codonSubsToSeq(wildtype, codon_subs, return_aa=False, aa_subs=False):
    """Convert codon substitutions to sequence.

    Args:
        `wildtype` (str)
            The wildtype sequence
        `codon_subs` (str)
            String of space delimited codon substitutions, in the format:
            OldCodonSiteNewCodon
        `return_aa` (bool)
            Specify whether to return sequence as nucleotide or amino acid.
            Default is nucleotide.
        `aa_subs` (bool)
            Specify whether the substitutions are in amino acid form
            rather than codon. Default is codon. `return_aa` must be True
            in order for `aa_subs` to be True, since there are numerous
            possible nucleotide sequences for an amino acid sequence.

    Returns:
        A str of the sequence with all the codon substitutions

    Here is an example:

    >>> codonSubsToSeq('ATGGAACAA', '')
    'ATGGAACAA'
    >>> codonSubsToSeq('ATGGAACAA', 'GAA2CAG')
    'ATGCAGCAA'
    >>> codonSubsToSeq('ATGGAACAA', 'ATG1GGG GAA2CAG')
    'GGGCAGCAA'
    """
    # Make sure you are not trying to convert amino acids to codons
    if aa_subs:
        if not return_aa==True:
            raise ValueError('Cannot return nucleotide sequence using aa subs')
    # Make sure the wildtype sequence is divisible into codons
    if len(wildtype) % 3  != 0:
        raise ValueError('`wildtype` not divisible by 3')

    codon_list = [wildtype[3 * r : 3 * r + 3] for
                  r in range(len(wildtype) // 3)]

    # Create parser for the codon substitutions
    if aa_subs:
        submatcher = re.compile('^(?P<wt>.)'
                                '(?P<site>\d+)'
                                '(?P<mut>.)$')
    if not aa_subs:
        submatcher = re.compile('^(?P<wt>[ATGC]{3})'
                                '(?P<site>\d+)'
                                '(?P<mut>[ATGC]{3})$')
    for codon in codon_subs.split():
        m = submatcher.match(codon)
        if not m:
            raise ValueError(f"Invalid codon substitution {codon}")
        site = int(m.group('site')) - 1
        if not (0 <= site < len(codon_list)):
            raise ValueError('Codon site out of bounds')
        if not aa_subs:
            if m.group('wt') != codon_list[site]:
                raise ValueError(f"Invalid wildtype codon in {codon}")
        else:
            if m.group('wt') != CODON_TO_AA[codon_list[site]]:
                raise ValueError(f"Invalid wildtype aa in {codon}")
        if m.group('wt') == m.group('mut'):
            raise ValueError(f"wildtype and mutant the same in {codon}")
        if not aa_subs:
            if m.group('mut') not in CODON_TO_AA:
                raise ValueError(f'Invalid mutant codon in {codon}')
        else:
            if m.group('mut') not in AAS_WITHSTOP:
                raise ValueError(f'Invalid mutant aa in {codon}')
        # Change to the mutant codon
        sub = m.group('mut')
        if aa_subs:
            sub = AA_TO_CODONS[sub][0]
            codon_list[site] = sub
        else:
            codon_list[site] = m.group('mut')

    # Return the sequence
    if not return_aa:
        return ''.join(codon_list)
    else:
        return ''.join(CODON_TO_AA[codon] for codon in codon_list)


def func_score_to_gpm(func_scores_df, wildtype, metric='func_score'):
    """Generate a genotype phenotype map from a functional score dataframe.

    Args:
        `func_scores_df` (functional score dataframe)
            A functional score dataframe (from the
            :meth:`CodonVariantTable.func_scores` method), narrowed down to one
            post sample condition, typically with something like:
            `func_scores_df.query('library == @library & sample == @sample')`
        `wildtype` (str)
            A string containing the wildtype sequence
        `metric` (str)
            A string specifying which metric to use as a phenotype

    Returns:
         A genotype phenotype map object from the Harms lab
         `gpmap package <https://github.com/harmslab/gpmap>`_
         to be used as in input into epistasis models in the epistasis package.
    """
    # Put the phenotypes into a list
    phenotypes = func_scores_df[metric].tolist()

    # Get standard deviations for each phenotype
    var = func_scores_df['func_score_var'].tolist()
    stdev = np.sqrt(var)

    # Get codon substitutions in a list
    substitutions = func_scores_df['aa_substitutions'].tolist()

    # Get a list of genotypes
    genotypes = []
    for subs in substitutions:
        genotype = codonSubsToSeq(wildtype, subs, return_aa=True, aa_subs=True)
        genotypes.append(genotype)

    # Get the wildtype amino acid sequence
    wildtype = codonSubsToSeq(wildtype, '', return_aa=True)

    # Generate the genotype phenotype map
    gpm = gpmap.GenotypePhenotypeMap(wildtype=wildtype, genotypes=genotypes,
                               phenotypes=phenotypes, stdeviations=stdev)

    return gpm




if __name__ == '__main__':
    import doctest
    doctest.testmod()
