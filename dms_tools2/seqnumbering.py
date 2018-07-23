"""
===============
seqnumbering
===============

Deals with sequence numbering conversions.
"""


import collections
import tempfile

import Bio.SeqIO


class transcriptConverter:
    """Convert sites in mRNA transcript to chromosome or CDS sites.

    In mRNA sequencing, we identify mutations with respect
    to their position in the mRNA. But we may want map them
    to the corresponding numbers in the entire chromosome 
    (gene segment in case of a segmented virus, or entire
    viral genome in the case of a non-segmented virus), or
    to CDSs encoded on that chromosome. 

    Args:
        `genbankfile` (str)
            Genbank file with one or more loci, each of which
            should be a separate chromosome. The features
            of relevance are called `mRNA` and `CDS`.
            Each of these features should have a qualifier
            called `name` that gives its name. 
        `indexstart` (int)
            Number assigned to first nucleotide in feature.
            Although Python uses 0-based indexing, you may
            want to set to 1 as sequences are normally
            numbered starting at 1.
        `ignore_other_features` (bool)
            If `genbankfile` contains features not `mRNA`
            or `CDS`, ignore them or raise error?

    Attributes:
        `indexstart` (int)
            Value set at initialization.
        `chromosomes` (list)
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
    ...                      /name="fluNS1"
    ...      mRNA            join(2..56,529..868)
    ...                      /name="fluNS2"
    ...      CDS             27..719
    ...                      /name="fluNS1"
    ...      CDS             join(27..56,529..864)
    ...                      /name="fluNS2"
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

    Now initialize a :class:`transcriptConverter`:

    >>> with tempfile.NamedTemporaryFile(mode='r+') as genbankfile:
    ...     _ = genbankfile.write(genbank_text)
    ...     genbankfile.flush()
    ...     _ = genbankfile.seek(0)
    ...     converter = transcriptConverter(genbankfile)

    Confirm resulting :class:`transcriptConverter` contains one
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
    

    """

    def __init__(self, genbankfile, *, indexstart=1,
            ignore_other_features=False):
        """See main class docstring."""

        self.indexstart = indexstart

        self.chromosomes = {c.name:c for c in
                Bio.SeqIO.parse(genbankfile, 'genbank')}

        self.CDSs = {}
        self.mRNAs = {}
        self.mRNA_chromosome = {}
        self.chromosome_CDSs = collections.defaultdict(list)
        for c in self.chromosomes.values():
            for feature in c.features:
                if feature.type in {'mRNA', 'CDS'}:
                    name = feature.qualifiers['name']
                    if len(name) != 1:
                        raise ValueError("multiple names")
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
                            .feature.type)

    def i_mRNAtoChromosome(self, mRNA, i, mRNAfragment=None):
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
        """
        raise RuntimeError('not yet implemented')



if __name__ == '__main__':
    import doctest
    doctest.testmod()
