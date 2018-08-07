"""
===============
seqnumbering
===============

Deals with sequence numbering conversions.
"""


import collections
import tempfile
import re

import Bio.SeqIO
import Bio.SeqUtils

from dms_tools2 import NTS


class TranscriptConverter:
    """Convert sites in mRNA transcript to chromosome or CDS sites.

    In mRNA sequencing, we identify mutations with respect
    to their position in the mRNA. But we may want map them
    to the corresponding numbers in the entire chromosome 
    (gene segment in case of a segmented virus, or entire
    viral genome in the case of a non-segmented virus), or
    to CDSs encoded on that chromosome.

    For all site numbers used as input and output for this
    class, numbering is assumed to be 1, 2, .... Note that
    this is different than the 0, 1, ... numbering used
    by Python for strings.

    Args:
        `genbankfile` (str)
            Genbank file with one or more loci, each of which
            should be a separate chromosome. The features
            of relevance are called `mRNA` and `CDS`.
            Each of these features should have a qualifier
            called `label` that gives its name. 
        `ignore_other_features` (bool)
            If `genbankfile` contains features not `mRNA`
            or `CDS`, ignore them or raise error?
        `to_upper` (bool)
            Convert all sequences to upper case letters?

    Attributes:
        `chromosomes` (dict)
            Keyed by chromosome names, values are
            `Bio.SeqRecord.SeqRecord` for chromosomes.
        `mRNAs` (dict)
            Keyed by mRNA names, values are
            `Bio.SeqFeature.SeqFeature` for mRNA.
        `CDSs` (dict)
            Keyed by CDS names, values are
            `Bio.SeqFeature.SeqFeature` for CDS.
        `mRNA_chromosome` (dict)
            Keyed by mRNA name, values is name of chromosome
            containing mRNA.
        `chromosome_CDSs` (dict)
            Keyed by chromosome names, values are names of
            CDSs on chromosome.

    Specify example `genbankfile` contents with required information.
    In this file, the chromosome is `fluNS`, and it encodes two mRNAs
    and two CDSs (for `fluNS1` and `fluNS2`):

    >>> genbank_text = '''
    ... LOCUS       fluNS                    890 bp    DNA              UNK 01-JAN-1980
    ... FEATURES             Location/Qualifiers
    ...      mRNA            2..868
    ...                      /label="fluNS1"
    ...      mRNA            join(2..56,529..868)
    ...                      /label="fluNS2"
    ...      CDS             27..719
    ...                      /label="fluNS1"
    ...      CDS             join(27..56,529..864)
    ...                      /label="fluNS2"
    ... ORIGIN
    ...         1 agcaaaagca gggtgacaaa gacataatgg atccaaacac tgtgtcaagc tttcaggtag
    ...        61 attgctttct ttggcatgtc cgcaaaagag ttgcagacca agaactaggt gatgccccat
    ...       121 tccttgatcg gcttcgccga gatcagaagt ccctaagagg aagaggcagc actcttggtc
    ...       181 tggacatcga aacagccacc cgtgctggaa agcaaatagt ggagcggatt ctgaaggaag
    ...       241 aatctgatga ggcactcaaa atgaccatgg cctctgtacc tgcatcgcgc tacctaactg
    ...       301 acatgactct tgaggaaatg tcaaggcact ggttcatgct catgcccaag cagaaagtgg
    ...       361 caggccctct ttgtatcaga atggaccagg cgatcatgga taagaacatc atactgaaag
    ...       421 cgaacttcag tgtgattttt gaccggctgg agactctaat attactaagg gccttcaccg
    ...       481 aagaggggac aattgttggc gaaatttcac cactgccctc tcttccagga catactgatg
    ...       541 aggatgtcaa aaatgcagtt ggggtcctca tcggaggact tgaatggaat aataacacag
    ...       601 ttcgagtctc tgaaactcta cagagattcg cttggagaag cagtaatgag aatgggagac
    ...       661 ctccactcac tccaaaacag aaacggaaaa tggcgggaac aattaggtca gaagtttgaa
    ...       721 gaaataaggt ggttgattga agaagtgaga cacagactga agataacaga gaatagtttt
    ...       781 gagcaaataa catttatgca agccttacaa ctattgcttg aagtggagca agagataaga
    ...       841 actttctcgt ttcagcttat ttaataataa aaaacaccct tgtttctact
    ... //
    ... '''

    Now initialize a :class:`TranscriptConverter`:

    >>> with tempfile.NamedTemporaryFile(mode='r+') as genbankfile:
    ...     _ = genbankfile.write(genbank_text)
    ...     genbankfile.flush()
    ...     _ = genbankfile.seek(0)
    ...     converter = TranscriptConverter(genbankfile)

    Confirm resulting :class:`TranscriptConverter` contains one
    chromosome (`fluNS`) with the expected to mRNAs and CDSs:

    >>> list(converter.chromosomes.keys())
    ['fluNS']
    >>> sorted(converter.CDSs.keys())
    ['fluNS1', 'fluNS2']
    >>> sorted(converter.mRNAs.keys())
    ['fluNS1', 'fluNS2']
    >>> converter.mRNA_chromosome['fluNS1']
    'fluNS'
    >>> converter.mRNA_chromosome['fluNS2']
    'fluNS'
    >>> converter.chromosome_CDSs['fluNS']
    ['fluNS1', 'fluNS2']

    Get site in chromosome (`fluNS`) that corresponds to a
    position in the `fluNS1` mRNA, and then do the same for
    the `fluNS2` mRNA, using
    :class:`TranscriptConverter.i_mRNAtoChromosome`. Then
    check nucleotide identities with
    :class:`TranscriptConverter.ntIdentity`:

    >>> converter.i_mRNAtoChromosome('fluNS1', 60)
    61
    >>> converter.ntIdentity('fluNS', 61)
    'A'
    >>> converter.i_mRNAtoChromosome('fluNS2', 58)
    531
    >>> converter.ntIdentity('fluNS', 531)
    'C'

    Do the same for substrings of `fluNS1` and `fluNS2`:

    >>> converter.i_mRNAtoChromosome('fluNS1', 2, mRNAfragment='GATTGCTTTCT')
    61
    >>> converter.i_mRNAtoChromosome('fluNS2', 9, mRNAfragment='TTTCAGGACATACTGATG')
    531

    Get amino-acid substitutions caused by chromosome point mutation.
    We can specify point mutation as string or tuple, and get
    output with 3 or 1-letter amino-acid codes:

    >>> converter.aaSubstitutions('fluNS', 'A61T')
    'fluNS1-Asp12Val'
    >>> converter.aaSubstitutions('fluNS', ('A', 61, 'T'))
    'fluNS1-Asp12Val'
    >>> converter.aaSubstitutions('fluNS', 'A61T', aa_3letter=False)
    'fluNS1-D12V'

    Now look at a point mutation that affects `fluNS1` and `fluNS2`, but
    only causes an amino-acid substitution in the former (is synonymous
    in the latter):

    >>> converter.aaSubstitutions('fluNS', 'C531T')
    'fluNS1-His169Tyr'

    Now mutation that causes amino-acid substitutions in
    `fluNS1` and `fluNS2`:

    >>> converter.aaSubstitutions('fluNS', 'A532T')
    'fluNS1-His169Leu_fluNS2-Ile12Leu'

    """

    def __init__(self, genbankfile, *,
            ignore_other_features=False, to_upper=True):
        """See main class docstring."""

        self.to_upper = to_upper

        self.chromosomes = {c.name:c for c in
                Bio.SeqIO.parse(genbankfile, 'genbank')}

        self.CDSs = {}
        self.mRNAs = {}
        self.mRNA_chromosome = {}
        self.chromosome_CDSs = collections.defaultdict(list)
        for c in self.chromosomes.values():
            if to_upper:
                c.seq = c.seq.upper()
            for feature in c.features:
                if feature.type in {'mRNA', 'CDS'}:
                    name = feature.qualifiers['label']
                    if len(name) != 1:
                        raise ValueError("multiple 'label' qualifiers")
                    feature.name = name[0]
                    if feature.type == 'mRNA':
                        d = self.mRNAs
                        self.mRNA_chromosome[feature.name] = c.name
                    elif feature.type == 'CDS':
                        d = self.CDSs
                        self.chromosome_CDSs[c.name].append(feature.name)
                    else:
                        raise RuntimeError('should not get here')
                    if feature.name in d:
                        raise ValueError("duplicate {0} {1}".format(
                                feature.type, feature.name))
                    d[feature.name] = feature
                elif not ignore_other_features:
                    raise ValueError("unrecognized feature type {0}"
                            .format(feature.type))


    def i_mRNAtoChromosome(self, mRNA, i, *, mRNAfragment=None):
        """Convert site number in mRNA to number in chromosome.

        Args:
            `mRNA` (str)
                Name of a valid mRNA.
            `i` (int)
                Site number in mRNA.
            `mRNAfragment` (str or `None`)
                Substring sequence of `mRNA`. In this case, `i`
                is taken as the site in this substring of the
                `mRNA`. Useful because sometimes mutations may
                be called in a fragment of the full mRNA. Set
                to `None` if `i` is site in full mRNA.

        Returns:
            Site number in chromosome that contains `mRNA`.
        """
        try:
            chromosome = self.chromosomes[self.mRNA_chromosome[mRNA]]
            mRNA_feature = self.mRNAs[mRNA]
        except KeyError:
            raise ValueError("invalid `mRNA` {0}".format(mRNA))

        mRNA_seq = mRNA_feature.extract(chromosome).seq

        if mRNAfragment:
            if self.to_upper:
                mRNAfragment = mRNAfragment.upper()
            n = mRNA_seq.count(mRNAfragment)
            if n == 1:
                i += mRNA_seq.find(mRNAfragment)
            elif n > 1:
                raise ValueError("`mRNA` {0} does contains {1} "
                        "copies of `mRNAfragment`".format(mRNA, n))
            else:
                raise ValueError("`mRNA` {0} does not contain "
                        "specified `mRNAfragment`".format(mRNA))

        mRNA_sites = list(mRNA_feature)
        if not (1 <= i <= len(mRNA_sites)):
            raise ValueError("`i` of {0} not valid site in `mRNA` {1}"
                    .format(i, mRNA))

        return mRNA_sites[i - 1] + 1


    def ntIdentity(self, chromosome, i):
        """Gets identity at site in chromosome.

        Args:
            `chromosome` (str)
                Name of chromosome.
            `i` (int)
                Site in chromosome.

        Returns:
            Nucleotide in `chromosome` at site `i`.
        """
        try:
            chromosome_seq = self.chromosomes[chromosome]
        except KeyError:
            raise ValueError("`chromosome` {0} does not exist"
                    .format(chromosome))

        if not (1 <= i <= len(chromosome_seq)):
            raise ValueError("`i` of {0} not in `chromosome` {1}"
                    .format(i, chromosome))

        return str(chromosome_seq[i - 1])


    def aaSubstitutions(self, chromosome, mutation, *,
            aa_3letter=True):
        """Amino-acid substitutions from point mutation to chromosome.

        Gets all amino-acid substitutions in CDSs caused by a
        point mutation to a chromosome.

        Args:
            `chromosome` (str)
                Name of chromosome.
            `mutation` (str or 3-tuple)
                Point mutation in chromosome based numbering. A str
                like "A50T", or tuple `(wt, i, mut)` where `wt` is
                wildtype nt, `i` is site, and `mut` is mutant nt.
            `aa_3letter` (bool)
                Use 3-letter rather than 1-letter amino-acid codes?

        Returns:
            A string giving the amino-acid mutations:
            
                - If no mutations, returns empty str.
                
                - If one mutation, will be str like 'fluNS1-Asp12Val'.

                - If several mutations, will be str like
                  'fluNS1-His169Leu_fluNS2-Ile12Leu'.
        """
        try:
            cds_names = self.chromosome_CDSs[chromosome]
            chromosome_seq = self.chromosomes[chromosome].seq
        except KeyError:
            raise ValueError("invalid `chromosome` {0}".format(chromosome))

        if isinstance(mutation, str):
            m = re.match('^(?P<wt>[{0}])(?P<i>\d+)(?P<mut>[{0}])$'
                    .format(''.join(NTS)), mutation)
            if not m:
                raise ValueError("cannot parse mutation {0}"
                        .format(mutation))
            mutation = (m.group('wt'), int(m.group('i')), m.group('mut'))

        if (len(mutation) != 3) or (mutation[0] not in NTS
                    ) or (mutation[2] not in NTS):
            raise ValueError("invalid mutation {0}".format(mutation))

        (wt, i, mut) = mutation

        if not (1 <= i <= len(chromosome_seq)):
            raise ValueError("`i` of {0} is not valid site in "
                    "`chromosome` {1}".format(i, chromosome))

        if wt != self.ntIdentity(chromosome, i):
            raise ValueError("Invalid wildtype nucleotide at site "
                    "{0} of {1}. Specified {2}, should be {3}".format(
                    i, chromosome, i, self.ntIdentity(chromosome, i)))
        
        if wt == mut:
            raise ValueError("wildtype and mutant nt are the same")

        cds_w_mut = [cds for cds in [self.CDSs[cds_name]
                for cds_name in cds_names] if (i - 1) in cds]

        mutstring = []
        mut_chromosome = None
        for cds_name in cds_names:
            cds = self.CDSs[cds_name]
            if (i - 1) in cds:
                if mut_chromosome is None:
                    mut_chromosome = chromosome_seq.tomutable()
                    mut_chromosome[i - 1] = mut
                wtprot = cds.extract(chromosome_seq).translate()
                mutprot = cds.extract(mut_chromosome).translate()
                aa_mut = [(wt, r0 + 1, mut) for r0, (wt, mut) in
                        enumerate(zip(wtprot, mutprot)) if wt != mut]
                if len(aa_mut) > 1:
                    raise ValueError(">1 amino-acid mutation")
                elif aa_mut:
                    wt_aa, r, mut_aa = aa_mut[0]
                    if aa_3letter:
                        wt_aa = Bio.SeqUtils.seq3(wt_aa)
                        mut_aa = Bio.SeqUtils.seq3(mut_aa)
                    mutstring.append('{0}-{1}{2}{3}'.format(cds_name,
                            wt_aa, r, mut_aa))

        return '_'.join(mutstring)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
