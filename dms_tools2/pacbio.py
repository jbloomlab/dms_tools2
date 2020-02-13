"""
===================
pacbio
===================

Tools for processing PacBio sequencing data.
"""


import os
import gzip
import re
import io
import math
import subprocess
import collections
import tempfile
import numbers

import regex
import numpy
import pandas
import pysam
import Bio.SeqFeature

# import dms_tools2.plot to set plotting contexts / themes
import dms_tools2
from dms_tools2.plot import COLOR_BLIND_PALETTE
from dms_tools2.plot import COLOR_BLIND_PALETTE_GRAY
import matplotlib.pyplot as plt
from plotnine import *


class CCS:
    """Class to handle results of ``ccs``.

    Holds results of PacBio ``ccs``.
    Has been tested on output of ``ccs`` version 3.0.0.

    This class reads all data into memory, and so you
    may need a lot of RAM if `ccsfile` is large.

    Args:
        `samplename` (str)
            Sample or sequencing run
        `ccsfile` (str)
            File created by ``ccs`` that holds the CCSs. The
            ``ccs`` program outputs BAM files. However, you
            can also pass FASTQ files generated from these
            BAM files using ``samtools bam2fq -T np,rq <bamfile>``
            (note that ``-T np,rq`` flag which is needed to
            preserve the number of passes and accuracy flags).
            The file format is determined from the file extension,
            and can be ``*.bam``, ``*.fastq``, ``*.fq``,
            ``*.fastq.gz``, or ``*.fq.gz``.
        `reportfile` (str or `None`)
            Report file created by ``ccs``, or
            `None` if you have no reports.

    Attributes:
        `samplename` (str)
            Name set at initialization
        `ccsfile` (str)
            ``ccs`` BAM file set at initialization
        `reportfile` (str or `None`)
            ``ccs`` report file set at initialization
        `zmw_report` (pandas.DataFrame or `None`):
            ZMW stats in `reportfile`, or `None` if no
            `reportfile`. Columns are *status*, *number*,
            *percent*, and *fraction*.
        `subread_report` (pandas.DataFrame or `None`)
            Like `zmw_report` but for subreads.
        `df` (pandas.DataFrame)
            The CCSs in `ccsfile`. Each row is a different CCS
            On creation, there will be the following columns (you
            can modify to add more): 

              - "name": the name of the CCS
              - "samplename": the sample as set via `samplename`
              - "CCS": the circular consensus sequence
              - "CCS_qvals": the Q-values as a numpy array
              - "passes": the number of passes of the CCS
              - "CCS_accuracy": the accuracy of the CCS
              - "CCS_length": the length of the CCS

    Here is an example.

    First, define the sequences, quality scores,
    and names for 3 example sequences. The names indicate
    the barcodes, the accuracy of the barcode, and the polarity.
    Two of the sequences have the desired termini and
    a barcode. The other does not. Note that the second
    sequence has an extra nucleotide at each end, this
    will turn out to be fine with the `match_str` we write.
    The second sequence is also reverse complemented:

    >>> termini5 = 'ACG'
    >>> termini3 = 'CTT'
    >>> ccs_seqs = [
    ...         {'name':'barcoded_TTC_0.999_plus',
    ...          'seq':termini5 + 'TTC' + 'ACG' + termini3,
    ...          'qvals':'?' * 12,
    ...         },
    ...         {'name':'barcoded_AGA_0.995_minus',
    ...          'seq':dms_tools2.utils.reverseComplement(
    ...                'T' + termini5 + 'AGA' + 'GCA' + termini3 + 'A'),
    ...          'qvals':''.join(reversed('?' * 4 + '5?9' + '?' * 7)),
    ...         },
    ...         {'name':'invalid',
    ...          'seq':'GGG' + 'CAT' + 'GCA' + termini3,
    ...          'qvals':'?' * 12,
    ...         }
    ...         ]
    >>> for iccs in ccs_seqs:
    ...     iccs['accuracy'] = qvalsToAccuracy(iccs['qvals'], encoding='sanger')

    Now place these in a block of text that meets the
    `CCS SAM specification <https://github.com/PacificBiosciences/unanimity/blob/develop/doc/PBCCS.md>`_:

    >>> sam_template = '\\t'.join([
    ...        '{0[name]}',
    ...        '4', '*', '0', '255', '*', '*', '0', '0',
    ...        '{0[seq]}',
    ...        '{0[qvals]}',
    ...        'np:i:6',
    ...        'rq:f:{0[accuracy]}',
    ...        ])
    >>> samtext = '\\n'.join([sam_template.format(iccs) for
    ...                      iccs in ccs_seqs])

    Create small SAM file with these sequences, then
    convert to BAM file used to initialize a :class:`CCS`
    (note this requires ``samtools`` to be installed):

    >>> samfile = '_temp.sam'
    >>> bamfile = '_temp.bam'
    >>> with open(samfile, 'w') as f:
    ...     _ = f.write(samtext)
    >>> _ = subprocess.check_call(['samtools', 'view',
    ...         '-b', '-o', bamfile, samfile])
    >>> ccs = CCS('test', bamfile, None)
    >>> os.remove(samfile)

    We also sometimes create the BAM files created by PacBio
    ``ccs`` to FASTQ. Do that using ``samtools bam2fq -T np,rq``
    to keep flags with number of passes and overall read quality:

    >>> fastq_data = subprocess.check_output(
    ...         ['samtools', 'bam2fq', '-T', 'np,rq', bamfile])

    Show how the resulting FASTQ data keeps the *np* and *rq* tags:

    >>> print(fastq_data.decode('utf-8').strip().replace('\\t', ' '))
    @barcoded_TTC_0.999_plus np:i:6 rq:f:0.999
    ACGTTCACGCTT
    +
    ????????????
    @barcoded_AGA_0.995_minus np:i:6 rq:f:0.998144
    TAAGTGCTCTCGTA
    +
    ???????9?5????
    @invalid np:i:6 rq:f:0.999
    GGGCATGCACTT
    +
    ????????????

    Write the FASTQ to a file, and check that :class:`CCS`
    initialized from the FASTQ is the same as one from the BAM:

    >>> fastqfile = '_temp.fastq'
    >>> gzfastqfile = '_temp.fastq.gz'
    >>> with open(fastqfile, 'wb') as f:
    ...     _ = f.write(fastq_data)
    >>> with gzip.open(gzfastqfile, 'wb') as f:
    ...     _ = f.write(fastq_data)
    >>> ccs_fastq = CCS('test', fastqfile, None)
    >>> ccs_gzfastq = CCS('test', gzfastqfile, None)
    >>> pandas.testing.assert_frame_equal(ccs_fastq.df, ccs.df)
    >>> pandas.testing.assert_frame_equal(ccs_gzfastq.df, ccs.df)
    >>> os.remove(fastqfile)
    >>> os.remove(gzfastqfile)
    >>> os.remove(bamfile)

    Check `ccs.df` has correct names, samplename, CCS sequences,
    and columns:

    >>> set(ccs.df.name) == {s['name'] for s in ccs_seqs}
    True
    >>> all(ccs.df.samplename == 'test')
    True
    >>> set(ccs.df.CCS) == {s['seq'] for s in ccs_seqs}
    True
    >>> set(ccs.df.columns) == {'CCS', 'CCS_qvals', 'name',
    ...         'passes', 'CCS_accuracy', 'CCS_length', 'samplename'}
    True

    Use :meth:`matchSeqs` to match sequences with expected termini
    and define barcodes and reads in these:

    >>> match_str = (termini5 + '(?P<barcode>N{3})' +
    ...         '(?P<read>N+)' + termini3)
    >>> ccs.df = matchSeqs(ccs.df, match_str, 'CCS', 'barcoded')

    This matching adds new columns to the new `ccs.df`:

    >>> set(ccs.df.columns) >= {'barcode', 'barcode_qvals',
    ...         'barcode_accuracy', 'read', 'read_qvals',
    ...         'read_accuracy', 'barcoded', 'barcoded_polarity'}
    True

    Now make sure `df` indicates that the correct sequences
    are barcoded, and that they have the correct barcodes:

    >>> bc_names = sorted([s['name'] for s in ccs_seqs if
    ...         'barcoded' in s['name']])
    >>> ccs.df = ccs.df.sort_values('barcode')
    >>> (ccs.df.query('barcoded').name == bc_names).all()
    True
    >>> barcodes = [x.split('_')[1] for x in bc_names]
    >>> (ccs.df.query('barcoded').barcode == barcodes).all()
    True
    >>> (ccs.df.query('not barcoded').barcode == ['']).all()
    True
    >>> barcode_accuracies = [float(x.split('_')[2]) for x in bc_names]
    >>> numpy.allclose(ccs.df.query('barcoded').barcode_accuracy,
    ...     barcode_accuracies, atol=1e-4)
    True
    >>> numpy.allclose(ccs.df.query('barcoded').barcode_accuracy,
    ...         [qvalsToAccuracy(qvals) for qvals in
    ...         ccs.df.query('barcoded').barcode_qvals])
    True
    >>> numpy.allclose(ccs.df.query('not barcoded').barcode_accuracy,
    ...     -1, atol=1e-4)
    True
    >>> barcoded_polarity = [{'plus':1, 'minus':-1}[x.split('_')[3]]
    ...         for x in bc_names]
    >>> (ccs.df.query('barcoded').barcoded_polarity == barcoded_polarity).all()
    True

    """

    def __init__(self, samplename, ccsfile, reportfile):
        """See main class doc string."""
        self.samplename = samplename

        assert os.path.isfile(ccsfile), f"can't find {ccsfile}"
        self.ccsfile = ccsfile

        self.reportfile = reportfile
        if self.reportfile is None:
            self.zmw_report = None
            self.subread_report = None
        else:
            assert os.path.isfile(reportfile), \
                    "can't find {0}".format(reportfile)
            # set `zmw_report` and `subread_report`
            self._parse_report()

        self._build_df_from_ccsfile()


    def __eq__(self, other):
        return self.__dict__ == other.__dict__


    def _parse_report(self):
        """Set `zmw_report` and `subread_report` using `reportfile`."""
        # match reports made by ccs 3.0.0
        reportmatch = regex.compile('^ZMW Yield\n(?P<zmw>(.+\n)+)\n\n'
                            'Subread Yield\n(?P<subread>(.+\n)+)$')

        with open(self.reportfile) as f:
            report = f.read()
        m = reportmatch.search(report)
        assert m, "Cannot match {0}\n\n{1}".format(
                self.reportfile, report)

        for read_type in ['zmw', 'subread']:
            df = (pandas.read_csv(
                        io.StringIO(m.group(read_type)),
                        names=['status', 'number', 'percent']
                        )
                  .assign(fraction=lambda x: 
                        x.percent.str.slice(None, -1)
                        .astype('float') / 100)
                  )
            setattr(self, read_type + '_report', df)


    def _build_df_from_ccsfile(self):
        """Builds `df` from `ccsfile`."""
        # read into dictionary
        d = collections.defaultdict(list)

        # get file type by extensions
        base, ext = [s.lower() for s in os.path.splitext(self.ccsfile)]
        if ext in {'.gz', '.gzip'}:
            gzipped = True
            ext = os.path.splitext(base)[1].lower()
        else:
            gzipped = False

        # extract data based on file extension
        if ext == '.bam':
            if gzipped:
                raise ValueError("Cannot handle gzipped BAM")
            for s in pysam.AlignmentFile(self.ccsfile, 'rb',
                                         check_sq=False):
                d['CCS'].append(s.query_sequence)
                d['CCS_qvals'].append(numpy.asarray(s.query_qualities,
                                                    dtype='int'))
                d['name'].append(s.query_name)
                d['passes'].append(s.get_tag('np'))
                d['CCS_accuracy'].append(s.get_tag('rq'))
                d['CCS_length'].append(s.query_length)
                d['samplename'].append(self.samplename)

        elif ext in {'.fq', '.fastq'}:
            headmatch = re.compile(r'^(?P<name>\S+)\s+'
                                   r'np:i:(?P<passes>\d+)\s+'
                                   r'rq:f:(?P<accuracy>\d+\.{0,1}\d*)')
            for a in pysam.FastxFile(self.ccsfile):
                if a.comment is not None:
                    head = f"{a.name} {a.comment}"
                else:
                    head = a.name
                m = headmatch.match(head)
                if not m:
                    raise ValueError(f"could not match {head}")
                d['CCS'].append(a.sequence)
                qvals = numpy.array([ord(qi) - 33 for qi in a.quality],
                                    dtype='int')
                d['CCS_qvals'].append(qvals)
                d['name'].append(m.group('name'))
                d['passes'].append(int(m.group('passes')))
                d['CCS_accuracy'].append(float(m.group('accuracy')))
                d['CCS_length'].append(len(a.sequence))
                d['samplename'].append(self.samplename)

        else:
            raise ValueError(f"invalid file extension {ext}")

        # create data frame
        self.df = pandas.DataFrame(d)

        # some checks on `df`
        assert self.df.name.size == self.df.name.unique().size,\
                "non-unique names for {0}".format(self.name)
        assert (self.df.CCS_length == self.df.CCS.apply(len)).all(),\
                "CCS not correct length"
        assert (self.df.CCS_length == self.df.CCS_qvals.apply(len)).all(),\
                "qvals not correct length"


TerminiVariantTag = collections.namedtuple(
        'TerminiVariantTag', ['termini', 'site', 'nucleotides'])
TerminiVariantTag.__doc__ = "Variant tag at termini."
TerminiVariantTag.termini.__doc__ = \
        "Location of tag: `termini5` or `termini3`."
TerminiVariantTag.site.__doc__ = \
        "Site of tag in termini (0, 1, ... numbering)."
TerminiVariantTag.nucleotides.__doc__ = \
        "A dict keyed variant nucleotides, values variant name."


class TerminiVariantTagCaller:
    """Call variant tags at termini of CCSs.

    Args:
        `features` (list)
            List of BioPython `SeqFeature` objects. Any
            features with a type attribute of `variant_tag`
            are taken to specify variant tags. These should
            consist of a single nucleotide, and have qualifiers
            that give the nucleotide for each variant. The
            features list should also have features with
            type attributes `termini5` and `termini3` used
            to determine which termin each tag falls in.
        `variants` (list)
            List of variant names, must have nucleotide for
            variant specified in qualifier for each variant tag.
        `trim_termini` (int)
            The amount trimmed from the 5' termini of `termini5`
            and the 3' termini of `termini3` when these are
            passed to :class:`TerminiVariantTagCaller.call`.
    
    Attributes:
        `variant_tags` (list)
            List of :class:`TerminiVariantTag` objects.
        `variants` (list)
            List of variant names set on initialization.
        `trim_termini` (int)
            Value set as argument on initialization.
        `termini5` (Bio.SeqFeature.SeqFeature)
            The 5' termini in `features`
        `termini3` (Bio.SeqFeature.SeqFeature)
            The 3' termini in `features`

    Here is an example. First, create list of features that has a
    single-nucleotide variant tag for each of two possible variants
    ('variant_1' and 'variant_2') in each termini:

    >>> SeqFeature = Bio.SeqFeature.SeqFeature
    >>> FeatureLocation = Bio.SeqFeature.FeatureLocation
    >>> features = [
    ...     SeqFeature(type='termini5', location=FeatureLocation(0, 147)),
    ...     SeqFeature(type='termini3', location=FeatureLocation(1303, 1342)),
    ...     SeqFeature(type='variant_tag', location=FeatureLocation(32, 33),
    ...         qualifiers={'variant_1':['A'], 'variant_2':['G']}),
    ...     SeqFeature(type='variant_tag', location=FeatureLocation(1310, 1311),
    ...         qualifiers={'variant_1':['T'], 'variant_2':['C']})
    ...     ]

    Now initialize the :class:`TerminiVariantTagCaller`:

    >>> caller = TerminiVariantTagCaller(features, trim_termini=4)
    >>> caller.variants
    ['variant_1', 'variant_2']
    >>> int(caller.termini5.location.start)
    0
    >>> int(caller.termini5.location.end)
    147
    >>> int(caller.termini3.location.start)
    1303
    >>> int(caller.termini3.location.end)
    1342
    >>> len(caller.variant_tags)
    2
    >>> caller.variant_tags[0].termini
    'termini5'
    >>> caller.variant_tags[0].site
    32
    >>> caller.variant_tags[0].nucleotides == {'A':'variant_1', 'G':'variant_2'}
    True
    >>> caller.variant_tags[1].termini
    'termini3'
    >>> caller.variant_tags[1].site
    7
    >>> caller.variant_tags[1].nucleotides == {'T':'variant_1', 'C':'variant_2'}
    True

    Do some example variant calling:

    >>> caller.call({'termini5':'GGCGTCACACTTTGCTATGCCATAGCATATTTATCC',
    ...              'termini3':'AGATCGGTAGAGCGTCGTGTAGGGAAAGAGTGTGG'})
    'variant_1'
    >>> caller.call({'termini5':'GGCGTCACACTTTGCTATGCCATAGCATGTTTATCC',
    ...              'termini3':'AGATCGGCAGAGCGTCGTGTAGGGAAAGAGTGTGG'})
    'variant_2'
    >>> caller.call({'termini5':'GGCGTCACACTTTGCTATGCCATAGCATGTTTATCC',
    ...              'termini3':'AGATCGGTAGAGCGTCGTGTAGGGAAAGAGTGTGG'})
    'mixed'
    >>> caller.call({'termini5':'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC',
    ...              'termini3':'AGATCGGTAGAGCGTCGTGTAGGGAAAGAGTGTGG'})
    'invalid'
    >>> caller.call({'termini5':'GGC', 'termini3':'AGAT'})
    'unknown'
    >>> caller.call({})
    'unknown'
    """

    def __init__(self, features, *, variants=['variant_1', 'variant_2'],
            trim_termini):
        """See main class docs."""
        features_dict = collections.defaultdict(list)
        for feature in features:
            features_dict[feature.type].append(feature)
        for termini in ['termini5', 'termini3']:
            if len(features_dict[termini]) != 1:
                raise ValueError(f"Failed to find exactly one {termini}")
            else:
                setattr(self, termini, features_dict[termini][0])

        if len(features_dict['variant_tag']) < 1:
            raise ValueError("no `variant_tag`s specified")

        self.trim_termini = trim_termini
        if trim_termini < 0:
            raise ValueError("trim_termini must be >= 0")

        self.variants = variants
        if len(self.variants) < 1:
            raise ValueError("no variants specified")

        self.variant_tags = []
        for variant_feature in features_dict['variant_tag']:
            if len(variant_feature) != 1:
                raise ValueError(f"variant not length 1: {variant_feature}")
            termini = [f for f in [self.termini5, self.termini3] if
                    variant_feature.location.start in f]
            if len(termini) != 1:
                raise ValueError("variant tag not in exactly one termini")
            self.variant_tags.append(
                    TerminiVariantTag(
                        termini=termini[0].type,
                        site=variant_feature.location.start - 
                            termini[0].location.start,
                        nucleotides={variant_feature.qualifiers[v][0]:v
                            for v in self.variants})
                    )

    def checkTagNucleotides(self, amplicon):
        """Check amplicon carrying tags has right ambiguous nucleotides.

        Arguments:
            `amplicon` (BioPython `SeqRecord`)
                The full amplicon that contains the termini
                and variant tags

        This method checks that the `amplicon` correctly has
        IUPAC ambiguous nucleotides that cover the possible
        diversity at the site of each variant tag. If so,
        it does nothing and returns `None`. If not, it raises
        a `ValueError`.
        """
        for variant_tag in self.variant_tags:
            termini = {'termini5':self.termini5,
                       'termini3':self.termini3}[variant_tag.termini]
            terminiseq = str(termini.location.extract(amplicon).seq)
            nt = terminiseq[variant_tag.site]
            if not all(re.match(dms_tools2.NT_TO_REGEXP[nt], variant_nt)
                    for variant_nt in variant_tag.nucleotides.keys()):
                raise ValueError(f"Nucleotide {nt} invalid for {variant_tag}")

    def call(self, termini_seqs):
        """Call variant identity.

        Args:
            `termini_seqs` (dict, namedtuple, pandas row)
                Some object that has attributes that can be
                accessed as `termini5` and `termini3`. These
                termini are assumed to have the amount specified
                by :class:`TerminiVariantTagCaller.trim_termini`
                trimmed from the 5' termini of `termini5` and
                the 3' termini of `termini3`.

        Returns:
            A str that can be any of the following:

                - If all tag sites in termini match the same variant,
                  return the name of that variant.
                - If different tag sites match different variants,
                  return "mixed".
                - If any tag sites have a nucleotide that matches no
                  known variants, return "invalid".
                - If `termini_seqs` lacks a termini or has a termini
                  that is too short to contain the tag sites, return
                  "unknown".
        """
        if not ('termini5' in termini_seqs and 'termini3' in termini_seqs):
            return "unknown"

        variants = []
        for variant_tag in self.variant_tags:
            i = variant_tag.site
            if variant_tag.termini == 'termini5':
                i -= self.trim_termini
            if len(termini_seqs[variant_tag.termini]) <= i:
                return 'unknown'
            nt = termini_seqs[variant_tag.termini][i]
            if nt in variant_tag.nucleotides:
                variants.append(variant_tag.nucleotides[nt])
            else:
                return 'invalid'
        if len(set(variants)) == 1:
            return variants[0]
        else:
            return 'mixed'


def matchAndAlignCCS(ccslist, mapper, *,
        termini5, gene, spacer, umi, barcode, termini3,
        termini5_fuzziness=0, gene_fuzziness=0,
        spacer_fuzziness=0, umi_fuzziness=0,
        barcode_fuzziness=0, termini3_fuzziness=0,
        targetvariants=None, mutationcaller=None,
        terminiVariantTagCaller=None,
        tagged_termini_remove_indels=True,
        rc_barcode_umi=True):
    """Identify CCSs that match pattern and align them.

    This is a convenience function that runs :meth:`matchSeqs`
    and :meth:`alignSeqs` for a common use case. It takes one
    or more :class:`CCS` objects, looks for CCS sequences in them
    that match a specific pattern, and aligns them to targets. It
    returns a pandas data frame with all the results. The CCS
    sequences are assumed to be molecules that have the following
    structure, although potentially in either orientation::

        5'-...-termini5-gene-spacer-umi-barcode-termini3-...-3'

    As indicated by the ``...``, there can be sequence before and
    after our expected pattern that we ignore. The gene element
    is the aligned to the targets. The full CCS is also aligned
    in the absence of the pattern matching.

    Args:
        `ccslist` (:class:`CCS` object or list of them)
            Analyze the CCS's in the `df` attributes. If there are
            multiple :class:`CCS` objectes, they are concatenated.
            However, they must have the same columns.
        `mapper` (:py:mod:`dms_tools2.minimap2.Mapper`)
            Mapper used to perform alignments.
        `termini5` (str or `None`)
            Expected sequence at 5' end as str that can be compiled
            to `regex` object. Passed through :meth:`re_expandIUPAC`.
            For instance, make it 'ATG|CTG' if the sequence might
            start with either `ATG` or `CTG`. Set to `None` if
            no expected 5' termini.
        `gene` (str)
            Like `termini5` but gives the gene to match. For instance,
            'N+' if the gene can be arbitrary sequence and length.
        `spacer` (str or `None`)
            Like `termini5`, but for the spacer after `gene`.
        `umi` (str or `None`)
            Like `termini5`, but for UMI.
        `barcode` (str or `None`)
            Like `termini5`, but for barcode. For instance, 'N{10}'
            if 10-nucleotide barcode.
        `termini3` (str or `None`)
            Like `termini5`, but for termini3.
        `termini5_fuzziness`, ..., `termini3_fuzziness` (int)
            The matching for the sequence patterns uses `regex`,
            which enables fuzzy matching. Set `termini5_fuzziness`
            to enable a specific number of differences (can be
            insertion, deletion, or mismatch) when matching
            `termini5`. Likewise for `gene_fuzziness`, etc.
            Note that the fuzzy matching uses the *BESTMATCH*
            flag to try to find the best fuzzy match.
            Note also that you can **not** both use fuzzy
            matching charcters in the strings to match (e.g.,
            `termini5` and set fuzziness to a value > 0:
            choose one or the other way to specify fuzzy matches.
        `targetvariants` (:class:`dms_tools2.minimap2.TargetVariants`)
            Call target variants. See docs for same argument to
            :meth:`alignSeqs`.
        `mutationcaller` (:class:`dms_tools2.minimap2.MutationCaller`)
            Call mutations. See docs for same argument to :meth:`alignSeqs`.
        `terminiVariantTagCaller` (:class:`TerminiVariantTagCaller`)
            Call variants in termini.
        `tagged_termini_remove_indels` (bool)
            If `terminiVariantTagCaller` is being used and this,
            is `True`, then use `remove_indels` flag when calling
            `matchSeqs` for the termini. This is useful if
            using fuzzy matching from the termini, as it aids
            in the calling of tags as it doesn't cause indels
            to misplace the tag. Has no meaning if
            `terminiVariantTagCaller` is not being used.
        `rc_barcode_umi` (bool)
            Do we reverse complement the `barcode` and `UMI` in the
            returned data frame relative to the orientation of
            the gene. Typically this is desirable because actual
            barcode sequencing goes in the reverse direction of the
            gene.

    Returns:
        A pandas dataframe that will have all columns already in the
        `df` attribute of the input :class:`CCS` objects with the
        following columns added:

        - `barcoded`: `True` if CCS matches full expected pattern,
          `False` otherwise.

        - `barcoded_polarity`: 1 of the match is in the polarity of
          the CCS, -1 if to the reverse complement, 0 if no match.

        - Columns named `termini5`, `gene`, `spacer`, `UMI`,
          `barcode`, and `termini3` (except if any of these elements
          are `None`). If `barcoded` is `True` for that CCS, these
          columns give the sequence for that element. If it is `False`,
          they are empty strings. There are likewise columns with
          these same names suffixed with "_accuracy" that give the CCS
          accuracy for that element, and columns suffixed with "_qvals"
          that give the quality scores for the elements.

        - For each of `termini5`, `spacer`, and `termini3` that are
          not `None`, a column named `has_termini5`, etc that
          indicates if that element is matched in isolate even if
          the full pattern is not matched.

        - `gene_aligned` is True if the CCS matches the expected
          pattern (is `barcoded`), and `gene` can further be
          aligned using `mapper`. It is `False` otherwise.

        - `gene_aligned_alignment`, `gene_aligned_target`,
          `gene_aligned_n_trimmed_query_start`,
          `gene_aligned_n_trimmed_query_end`,
          `gene_aligned_n_trimmed_target_start`,
          `gene_aligned_n_trimmed_target_end`,
          `gene_aligned_n_additional`, and
          `gene_aligned_n_additional_difftarget` give the
          :py:mod:`dms_tools2.minimap2.Alignment`, the alignment
          target, number of nucleotides trimmed from ends of
          the query gene or target, the number
          of additional alignments if `gene_aligned`,
          and the number of additional alignments to different
          targets (see `target_isoforms` attribute of
          :py:mod:`dms_tools2.minimap2.Mapper`). If
          the gene is not aligned, these are `None`,
          empty strings, or -1.

        - If `targetvariants` is not `None`, column named
          `gene_aligned_target_variant` giving target variant
          returned by :class:`dms_tools2.minimap2.TargtVariants.call`.

        - If `mutationcaller` is not `None`, column named
          `gene_aligned_mutations` giving the
          :class:`dms_tools2.minimap2.Mutations` object returned
          by :class:`dms_tools2.minimap2.MutationCaller.call`,
          or `None` if there is no alignment.

        - If `terminiVariantTagCaller` is not `None`, column
          named `termini_variant` giving the termini variant
          returned by :class:`TerminiVariantTagCaller.call`,
          or the str "unknown" if both termini are not matched.

        - `CCS_aligned` is `True` if the CCS can be aligned
          using `mapper` even if a gene cannot be matched,
          and `False` otherwise. `CCS_aligned_alignment`
          and `CCS_aligned_target` give the
          :py:mod:`dms_tools2.minimap2.Alignment` (or `None`)
          and the target (or empty string).
    """
    if isinstance(ccslist, collections.Iterable):
        col_list = [ccs.df.columns for ccs in ccslist]
        assert all([col_list[0].equals(col) for col in col_list]),\
                "the CCS.df's in `ccslist` don't have same columns"
        df = pandas.concat([ccs.df for ccs in ccslist])
    else:
        df = ccslist.df

    # internal function:
    def _align_CCS_both_orientations(df, mapper):
        """Try align CCS both ways, adds columns.
          `CCS_aligned`, `CCS_aligned_alignment`, and
        `CCS_aligned_target`."""
        df_bi = (df.pipe(dms_tools2.pacbio.alignSeqs,
                         mapper=mapper,
                         query_col='CCS',
                         aligned_col='CCS_for_aligned')
                   .assign(CCS_rev=lambda x: x.CCS.map(
                           dms_tools2.utils.reverseComplement))
                   .pipe(dms_tools2.pacbio.alignSeqs,
                         mapper=mapper,
                         query_col='CCS_rev',
                         aligned_col='CCS_rev_aligned')
                   )
        return (df.assign(CCS_aligned=df_bi.CCS_for_aligned |
                          df_bi.CCS_rev_aligned)
                .assign(CCS_aligned_alignment=
                        df_bi.CCS_for_aligned_alignment.where(
                        df_bi.CCS_for_aligned,
                        df_bi.CCS_rev_aligned_alignment))
                .assign(CCS_aligned_target=lambda x:
                        x.CCS_aligned_alignment.map(
                        lambda x: x.target if x is not None else ''))
                )

    # build match_str
    match_str = collections.OrderedDict()

    fuzz = {'termini5':termini5_fuzziness,
            'gene':gene_fuzziness,
            'spacer':spacer_fuzziness,
            'UMI':umi_fuzziness,
            'barcode':barcode_fuzziness,
            'termini3':termini3_fuzziness}
    has_fuzz = any(f > 0 for f in fuzz.values())

    seqs = {'termini5':termini5,
            'gene':gene,
            'spacer':spacer,
            'UMI':umi,
            'barcode':barcode,
            'termini3':termini3}

    for s in ['termini5', 'gene', 'spacer', 'UMI', 'barcode', 'termini3']:
        if seqs[s] is not None:
            if has_fuzz:
                match_str[s] = f"(?P<{s}>{seqs[s]}){{e<={fuzz[s]}}}"
                if '{' in seqs[s] or '}' in seqs[s]:
                    raise ValueError('Using fuzziness and fuzzy match in'
                        f" {s}:\nfuzziness = {fuzz[s]}\nseq = {seqs[s]}")
            else:
                match_str[s] = f"(?P<{s}>{seqs[s]})"
        else:
            match_str[s] = None

    if tagged_termini_remove_indels and (
            terminiVariantTagCaller is not None):
        remove_indels = ['termini5', 'termini3']
    else:
        remove_indels = []

    # now create df
    df = (
        df

        # match barcoded sequences
        .pipe(dms_tools2.pacbio.matchSeqs,
              match_str=''.join(m for m in match_str.values()
                    if m is not None),
              col_to_match='CCS',
              match_col='barcoded',
              remove_indels=remove_indels)
    
        # look for just termini or spacer
        .pipe(dms_tools2.pacbio.matchSeqs, 
              match_str=match_str['termini5'],
              col_to_match='CCS',
              match_col='has_termini5',
              add_polarity=False,
              add_group_cols=False)
        .pipe(dms_tools2.pacbio.matchSeqs, 
              match_str=match_str['termini3'],
              col_to_match='CCS',
              match_col='has_termini3',
              add_polarity=False,
              add_group_cols=False)
        .pipe(dms_tools2.pacbio.matchSeqs, 
              match_str=match_str['spacer'],
              col_to_match='CCS',
              match_col='has_spacer',
              add_polarity=False,
              add_group_cols=False)
    
        # see if gene aligns in correct orientation
        .pipe(dms_tools2.pacbio.alignSeqs,
              mapper=mapper,
              query_col='gene',
              aligned_col='gene_aligned',
              targetvariants=targetvariants,
              mutationcaller=mutationcaller)
    
        # look for any alignment of CCS, take best in either orientation
        .pipe(_align_CCS_both_orientations,
              mapper=mapper)
        )

    if terminiVariantTagCaller is not None:
        df = df.assign(termini_variant=lambda x: x.apply(
                terminiVariantTagCaller.call, axis=1))

    # reverse complement barcode and UMI
    if rc_barcode_umi:
        if barcode is not None:
            df.barcode = df.barcode.map(dms_tools2.utils.reverseComplement)

        if umi is not None:
            df.UMI = df.UMI.map(dms_tools2.utils.reverseComplement)

    return df


def matchSeqs(df, match_str, col_to_match, match_col, *,
        add_polarity=True, add_group_cols=True, remove_indels=[],
        add_accuracy=True, add_qvals=True,
        expandIUPAC=True, overwrite=False):
    """Identify sequences in a dataframe that match a specific pattern.

    Args:
        `df` (pandas DataFrame)
            Data frame with column holding sequences to match.
        `match_str` (str)
            A string that can be passed to `regex.compile` that gives
            the pattern that we are looking for, with target 
            subsequences as named groups. See also the `expandIUPAC`
            parameter, which simplifies writing `match_str`.
            If `None` we just return `df`. Note that we use
            `regex` rather than `re`, so fuzzy matching is
            enabled. Note that the matching uses the *BESTMATCH*
            flag to find the best match.
        `col_to_match` (str)
            Name of column in `df` that contains the sequences
            to match.
        `match_col` (str)
            Name of column added to `df`. Elements of columns are
            `True` if `col_to_match` matches `match_str` for that
            row, and `False` otherwise.
        `add_polarity` (bool)
            Add a column specifying the polarity of the match?
        `add_group_cols` (bool)
            Add columns with the sequence of every group in
            `match_str`?
        `remove_indels` (list)
            Only meaningful if `match_str` specifies to allow
            fuzzy matching for a group, and `add_group_cols`
            is `True`. Then for each named group in `match_str`,
            in sequence for that group that is added to the
            returned `df`, indicate indels by adding a `-`
            gap character for deletions, and removing the
            inserted nucleotide called by regex if there
            is an insertion.
        `add_accuracy` (bool)
            For each group in the match, add a column giving
            the accuracy of that group's sequence? Only used
            if `add_group_cols` is `True`.
        `add_qvals` (bool)
            For each group in the match, add a column giving
            the Q values for that group's sequence? Only used if
            `add_group_cols` is `True`.
        `expandIUPAC` (bool)
            Use `IUPAC code <https://en.wikipedia.org/wiki/Nucleic_acid_notation>`_
            to expand ambiguous nucleotides (e.g., "N") by passing
            `match_str` through the :meth:`re_expandIUPAC` function.
        `overwrite` (bool)
            If `True`, we overwrite any existing columns to
            be created that already exist. If `False`, raise
            an error if any of the columns already exist.

    Returns:
        A **copy** of `df` with new columns added. The exact columns
        to add are specified by the calling arguments. Specifically:

            - We always add a column with the name given by `match_col`
              that is `True` if there was a match and `False` otherwise.

            - If `add_polarity` is `True`, add a column that is
              `match_col` suffixed by "_polarity" which is 1 if
              the match is directly to the sequence in `col_to_match`,
              and -1 if it is to the reverse complement of this sequence.
              The value is 0 if there is no match.

            - If `add_group_cols` is `True`, then for each group
              in `match_str` specified using the `re` group naming
              syntax, add a column with that group name that
              gives the sequence matching that group. These
              sequences are empty strings if there is no match.
              These added sequences are in the polarity of the
              match, so if the sequence in `match_col` has
              to be reverse complemented for a match, then these
              sequences will be the reverse complement that matches.
              Additionally, when `add_group_cols` is True:

                - If `add_accuracy` is `True`, we also add a column
                  suffixed by "_accuracy" that gives the
                  accuracy of that group as computed from the Q-values.
                  The value -1 if there is match for that row. Adding
                  accuracy requires a colum in `df` with the name
                  given by `match_col` suffixed by "_qvals."

                - If `add_qvals` is `True`, we also add a column 
                  suffixed by "_qvals" that gives the Q-values
                  for that sequence. Adding these Q-values requires
                  that there by a column in `df` with the name given by
                  `match_col` suffixed by "_qvals". The Q-values are
                  in the form of a numpy array, or an empty numpy array
                  if there is no match for that row.
              
    See docs for :class:`CCS` for example uses of this function.

    Here is a short example that uses the fuzzy matching of
    the `regex` model for the polyA tail:

    >>> gene = 'ATGGCT'
    >>> polyA = 'AAAACAAAA'
    >>> df_in = pandas.DataFrame({'CCS':[gene + polyA]})
    >>> match_str = '(?P<gene>N+)(?P<polyA>AA(A{5,}){e<=1}AA)'
    >>> df = matchSeqs(df_in, match_str, 'CCS', 'matched',
    ...         add_accuracy=False, add_qvals=False)
    >>> expected = df.assign(gene=gene, polyA=polyA,
    ...         matched=True, matched_polarity=1)
    >>> (df.sort_index(axis=1) == expected.sort_index(axis=1)).all().all()
    True

    Here is a short example with fuzzy matching that uses the
    `remove_indels` option.
    First, do not remove the indels:

    >>> termini5 = 'ACAT'
    >>> termini3 = 'ATAC'
    >>> match_str2 = '^(?P<termini5>AAT){e<=1}(?P<gene>ATGGCT){e<=1}(?P<termini3>ATGAC){e<=1}$'
    >>> df_in2 = pandas.DataFrame({'CCS':[termini5 + gene + termini3]})
    >>> df2 = matchSeqs(df_in2, match_str2, 'CCS', 'matched',
    ...         add_accuracy=False, add_qvals=False)
    >>> df2.gene.values[0] == gene
    True
    >>> df2.termini5.values[0] == termini5
    True
    >>> df2.termini3.values[0] == termini3
    True

    Now remove the indels in just *termini3*:

    >>> df2_rm = matchSeqs(df_in2, match_str2, 'CCS', 'matched',
    ...         remove_indels=['termini3'],
    ...         add_accuracy=False, add_qvals=False)
    >>> df2_rm.gene.values[0] == gene
    True
    >>> df2_rm.termini5.values[0] == termini5
    True
    >>> df2_rm.termini3.values[0] == termini3
    False

    Now remove indels in **both** termini:

    >>> df2_rm2 = matchSeqs(df_in2, match_str2, 'CCS', 'matched',
    ...         remove_indels=['termini5', 'termini3'],
    ...         add_accuracy=False, add_qvals=False)
    >>> df2_rm2.gene.values[0] == gene
    True
    >>> df2_rm2.termini5.values[0] == termini5
    False
    >>> df2_rm2.termini5.values[0]
    'AAT'
    >>> df2_rm2.termini3.values[0] == termini3
    False
    """

    if match_str is None:
        return df

    assert col_to_match in df.columns, \
            "`df` lacks `col_to_match` column {0}".format(col_to_match)

    if expandIUPAC:
        match_str = re_expandIUPAC(match_str)
    matcher = regex.compile(match_str, flags=regex.BESTMATCH)

    newcols = [match_col]
    if add_polarity:
        polarity_col = match_col + '_polarity'
        newcols.append(polarity_col)

    if add_group_cols:
        groupnames = list(matcher.groupindex.keys())
        if len(set(groupnames)) != len(groupnames):
            raise ValueError("duplicate group names in {0}"
                             .format(match_str))
        newcols += groupnames
        if add_accuracy:
            newcols += [g + '_accuracy' for g in groupnames]
        if add_qvals:
            newcols += [g + '_qvals' for g in groupnames]
        if add_accuracy or add_qvals:
            match_qvals_col = col_to_match + '_qvals'
            if match_qvals_col not in df.columns:
                raise ValueError("To use `add_accuracy` or "
                        "`add_qvals`, you need a column in `df` "
                        "named {0}".format(match_qvals_col))
        if remove_indels:
            if set(remove_indels) > set(remove_indels):
                raise ValueError("`remove_indels` specifies "
                        "unknown group(s)")
        else:
            remove_indels = []
    else:
        groupnames = []
        if remove_indels:
            raise ValueError("can't use `remove_indels` without "
                             "using `add_group_cols`")

    # make sure created columns don't already exist
    dup_cols = set(newcols).intersection(set(df.columns))
    if not overwrite and dup_cols:
        raise ValueError("`df` already contains some of the "
                "columns that we are supposed to add:\n{0}"
                .format(dup_cols))

    # look for matches for each row
    match_d = {c:[] for c in newcols}
    for tup in df.itertuples():
        s = getattr(tup, col_to_match)
        m = matcher.search(s)
        if add_group_cols and (add_accuracy or add_qvals):
            qs = getattr(tup, match_qvals_col)
        if m:
            polarity = 1
        else:
            m = matcher.search(dms_tools2.utils.reverseComplement(s))
            polarity = -1
            if add_group_cols and (add_accuracy or add_qvals):
                qs = numpy.flip(qs, axis=0)
        if m:
            match_d[match_col].append(True)
            if add_polarity:
                match_d[polarity_col].append(polarity)
            ins_sites = m.fuzzy_changes[1]
            del_sites = m.fuzzy_changes[2]
            for g in groupnames:
                g_start = m.start(g)
                g_end = m.end(g)
                if g in remove_indels:
                    g_ins_sites = [i for i in ins_sites
                                  if g_start <= i < g_end]
                    g_del_sites = [i for i in del_sites
                                   if g_start <= i < g_end]
                else:
                    g_ins_sites = []
                    g_del_sites = []
                if add_qvals:
                    g_qs = qs[g_start : g_end]
                    g_qs_list = []
                if g_ins_sites or g_del_sites:
                    g_seq = []
                    for i, x in enumerate(m.group(g)):
                        if i + g_start in g_ins_sites:
                            pass
                        elif i + g_start in g_del_sites:
                            g_seq.append(x + '-')
                            if add_qvals:
                                g_qs_list.append(g_qs[i])
                                g_qs_list.append(numpy.nan)
                        else:
                            g_seq.append(x)
                            if add_qvals:
                                g_qs_list.append(g_qs[i])
                    g_seq = ''.join(g_seq)
                    if add_qvals:
                        g_qs = numpy.array(g_qs_list)
                else:
                    g_seq = m.group(g)
                match_d[g].append(g_seq)
                if add_qvals:
                    match_d[g + '_qvals'].append(g_qs)
                if add_accuracy:
                    match_d[g + '_accuracy'].append(qvalsToAccuracy(
                            qs[g_start : g_end]))
        else:
            match_d[match_col].append(False)
            if add_polarity:
                match_d[polarity_col].append(0)
            for g in groupnames:
                match_d[g].append('')
                if add_qvals:
                    match_d[g + '_qvals'].append(numpy.array([], dtype='int'))
                if add_accuracy:
                    match_d[g + '_accuracy'].append(-1)

    # set index to make sure matches `df`
    indexname = df.index.name
    assert indexname not in match_d
    match_d[indexname] = df.index.tolist()
    if (not overwrite) and dup_cols:
        raise ValueError("overwriting columns")
    return pandas.concat(
            [df.drop(dup_cols, axis=1),
                pandas.DataFrame(match_d).set_index(indexname),
            ],
            axis=1)


def alignSeqs(df, mapper, query_col, aligned_col, *,
        add_alignment=True, add_target=True,
        add_n_trimmed=True, add_n_additional=True,
        add_n_additional_difftarget=True, targetvariants=None,
        mutationcaller=None, overwrite=True, paf_file=None):
    """Align sequences in a dataframe to target sequence(s).

    Arguments:
        `df` (pandas DataFrame)
            Data frame in which one column holds sequences to match.
            There also must be a column named "name" with unique names.
        `mapper` (:py:mod:`dms_tools2.minimap2.Mapper`)
            Align using the :py:mod:`dms_tools2.minimap2.Mapper.map`
            function of `mapper`. Target sequence(s) to which
            we align are specified when initializing `mapper`.
        `query_col` (str)
            Name of column in `df` with query sequences to align.
            If we are to use Q-values, there must also be a column
            with this name suffixed by "_qvals".
        `aligned_col` (str)
            Name of column added to `df`. Elements of column are
            `True` if `query_col` aligns, and `False` otherwise.
        `add_alignment` (bool)
            Add column with the :py:mod:`dms_tools2.minimap2.Alignment`.
        `add_target` (bool)
            Add column giving target (reference) to which sequence
            aligns.
        `add_n_trimmed` (bool)
            Add columns giving number of nucleotides trimmed from
            ends of both the query and target in the alignment.
        `add_n_additional` (bool)
            Add column specifying the number of additional
            alignments.
        `targetvariants` (:class:`dms_tools2.minimap2.TargetVariants`)
            Call target variants of aligned genes using the `call`
            function of this object. Note that this also adjusts
            the returned alignments / CIGAR if a variant is called.
            If the `variantsites_min_acc` attribute is not `None`,
            then `df` must have a column with the name of `query_col`
            suffixed by '_qvals' that gives the Q-values to compute
            accuracies.
        `mutationcaller` (:class:`dms_tools2.minimap2.MutationCaller`)
            Call mutations of aligned genes using the `call` function
            of this object. Note that any target variant mutations are
            handled first and then removed and not called here.
        `add_n_additional_difftarget` (bool)
            Add columns specifying number of additional alignments
            to a target other than the one in the primary alignment.
        `overwrite` (bool)
            If `True`, we overwrite any existing columns to
            be created that already exist. If `False`, raise
            an error if any of the columns already exist.
        `paf_file` (`None` or str)
            If a str, is the name of the PAF file created
            by `mapper` (see `outfile` argument of
            :py:mod:`dms_tools2.minimap2.Mapper.map`) Otherwise
            this file is not saved.

    Returns:
        A **copy** of `df` with new columns added. The exact
        columns to add are specified by the calling arguments.
        Specifically:

            - We always add a column with the name given by
              `aligned_col` that is `True` if there was an
              alignment and `False` otherwise.
              
            - If `add_alignment` is `True`, add column named
              `aligned_col` suffixed by "_alignment" that gives
              the alignment as a :py:mod:`dms_tools2.minimap2.Alignment`
              object, or `None` if there is no alignment. Note that
              if there are multiple alignments, then this is the
              "best" alignment, and the remaining alignments are in
              the :py:mod:`dms_tools2.minimap2.Alignment.additional`
              attribute.

            - If `add_target` is `True`, add column named
              `aligned_col` suffixed by "_target" that gives
              the target to which the sequence aligns in the
              "best" alignment, or an empty string if no alignment.

            - If `add_n_trimmed` is `True`, add column named
              `aligned_col` suffixed by "_n_trimmed_query_start",
              "_n_trimmed_query_end", "_n_trimmed_target_start",
              and "_n_trimmed_target_end" that give the number
              of nucleotides trimmed from the query and target
              in the "best" alignment. Are all zero if the
              zero if the alignment is end-to-end. Are -1 if no
              alignment.

            - If `add_n_additional` is `True`, add column
              named `aligned_col` suffixed by "_n_additional" that
              gives the number of additional alignments (in
              :py:mod:`dms_tools2.minimap2.Alignment.additional`),
              or -1 if there is no alignment.

            - If `add_n_additional_difftarget` is `True`, add column
              named `aligned_col` suffixed by "_n_additional_difftarget"
              that gives the number of additional alignments to
              **different** targets that are not isoforms, or -1
              if if there is no alignment. See the `target_isoforms`
              attribute of :py:mod:`dms_tools2.minimap2.Mapper`.

            - If `targetvariants` is not `None`, add a column
              named `aligned_col` suffixed by "_target_variant"
              that has the values returned for that alignment by
              :class:`dms_tools2.minimap2.TargetVariants.call`, or
              an empty string if no alignment.

            - If `mutationcaller` is not `None`, column named
              `aligned_col` suffixed by "_mutations" giving the
              :class:`dms_tools2.minimap2.Mutations` object returned
              by :class:`dms_tools2.minimap2.MutationCaller.call`,
              or `None` if there is no alignment.
    """
    assert query_col in df.columns, "no `query_col` {0}".format(query_col)

    newcols = [aligned_col]
    if add_alignment:
        alignment_col = aligned_col + '_alignment'
        newcols.append(alignment_col)
    if add_target:
        target_col = aligned_col + '_target'
        newcols.append(target_col)
    if add_n_trimmed:
        n_trimmed_prefix = aligned_col + '_n_trimmed_'
        for suffix in ['query_start', 'query_end',
                'target_start', 'target_end']:
            newcols.append(n_trimmed_prefix + suffix)
    if add_n_additional:
        n_additional_col = aligned_col + '_n_additional'
        newcols.append(n_additional_col)
    if add_n_additional_difftarget:
        n_additional_difftarget_col = (
                aligned_col + '_n_additional_difftarget')
        newcols.append(n_additional_difftarget_col)
    qvals_col = query_col + '_qvals'
    if qvals_col in df.columns:
        qvals = pandas.Series(df[qvals_col].values,
                              index=df.name).to_dict()
    else:
        qvals = collections.defaultdict(lambda: math.nan)
    if targetvariants is not None:
        targetvariant_col = aligned_col + '_target_variant'
        newcols.append(targetvariant_col)
        if targetvariants.variantsites_min_acc is not None:
            if qvals_col not in df.columns:
                raise ValueError("Cannot use `variantsites_min_acc` "
                        "of `targetvariants` as there is not a column "
                        "in `df` named {0}".format(qvals_col))
    if mutationcaller is not None:
        mutations_col = aligned_col + '_mutations'
        newcols.append(mutations_col)

    assert len(newcols) == len(set(newcols))

    dup_cols = set(newcols).intersection(set(df.columns))
    if (not overwrite) and dup_cols:
        raise ValueError("`df` already contains these columns:\n{0}"
                         .format(dup_cols))

    # perform the mapping
    assert len(df.name) == len(df.name.unique()), \
            "`name` in `df` not unique"
    with tempfile.NamedTemporaryFile(mode='w') as queryfile:
        queryfile.write('\n'.join([
                        '>{0}\n{1}'.format(*tup) for tup in
                        df.query('{0} != ""'.format(query_col))
                            [['name', query_col]]
                            .itertuples(index=False, name=None)
                        ]))
        map_dict = mapper.map(queryfile.name, outfile=paf_file)

    align_d = {c:[] for c in newcols}
    for name in df.name:
        if name in map_dict:
            a = map_dict[name]
            assert a.strand == 1, "method does not handle - polarity"
            if targetvariants:
                (variant, a) = targetvariants.call(a, qvals[name])
                align_d[targetvariant_col].append(variant)
            if mutationcaller:
                align_d[mutations_col].append(mutationcaller.call(a,
                        qvals[name]))
            align_d[aligned_col].append(True)
            if add_alignment:
                align_d[alignment_col].append(a)
            if add_target:
                align_d[target_col].append(a.target)
            if add_n_trimmed:
                align_d[n_trimmed_prefix + 'query_start'].append(
                        a.q_st)
                align_d[n_trimmed_prefix + 'query_end'].append(
                        a.q_len - a.q_en)
                align_d[n_trimmed_prefix + 'target_start'].append(
                        a.r_st)
                align_d[n_trimmed_prefix + 'target_end'].append(
                        a.r_len - a.r_en)
            if add_n_additional:
                align_d[n_additional_col].append(len(a.additional))
            if add_n_additional_difftarget:
                align_d[n_additional_difftarget_col].append(
                        len([a2.target for a2 in a.additional if
                        a2.target not in mapper.target_isoforms[a.target]]))

        else:
            align_d[aligned_col].append(False)
            if add_alignment:
                align_d[alignment_col].append(None)
            if add_target:
                align_d[target_col].append('')
            if add_n_trimmed:
                for suffix in ['query_start', 'query_end',
                        'target_start', 'target_end']:
                    align_d[n_trimmed_prefix + suffix].append(-1)
            if add_n_additional:
                align_d[n_additional_col].append(-1)
            if add_n_additional_difftarget:
                align_d[n_additional_difftarget_col].append(-1)
            if targetvariants:
                align_d[targetvariant_col].append('')
            if mutationcaller:
                align_d[mutations_col].append(None)

    # set index to make sure matches `df`
    index_name = df.index.name
    assert index_name not in align_d
    align_d[index_name] = df.index.tolist()
    if (not overwrite) and dup_cols:
        raise ValueError("overwriting columns")
    return pandas.concat(
            [df.drop(dup_cols, axis=1),
                pandas.DataFrame(align_d).set_index(index_name),
            ],
            axis=1)


def qvalsToAccuracy(qvals, encoding='numbers', no_avg=False):
    r"""Converts set of quality scores into average accuracy.

    Args:
        `qvals` (numpy array or number or str)
            List of Q-values, assumed to be Phred scores.
            For how they are encoded, see `encoding`.
        `encoding` (str)
            If it is "numbers" then `qvals` should be a
            numpy array giving the Q-values, or a number
            with one Q-value. If it is "sanger", then `qvals`
            is a string, with the score being the ASCII value
            minus 33.
        `no_avg` (bool)
            Compute the accuracies of individual Q-values
            rather than the average of the array or list.

    Returns:
        A number giving the average accuracy, or 
        `nan` if `qvals` is empty.

    Note that the probability :math:`p` of an error at a
    given site is related to the Q-value :math:`Q` by
    :math:`Q = -10 \log_{10} p`.

    >>> qvals = numpy.array([13, 77, 93])
    >>> round(qvalsToAccuracy(qvals), 3) == 0.983
    True
    >>> round(qvalsToAccuracy(qvals[1 : ]), 3) == 1
    True
    >>> qvalsToAccuracy(numpy.array([]))
    nan

    >>> qvals_str = '.n~'
    >>> round(qvalsToAccuracy(qvals_str, encoding='sanger'), 3) == 0.983
    True

    >>> round(qvalsToAccuracy(15), 3) == 0.968
    True

    >>> [round(a, 5) for a in qvalsToAccuracy(qvals, no_avg=True)] == [0.94988, 1, 1]
    True
    """
    if encoding == 'numbers':
        if isinstance(qvals, numbers.Number):
            qvals = numpy.array([qvals])
            no_avg = False
        elif isinstance(qvals, list):
            qvals = numpy.array(qvals)

    if qvals is None or len(qvals) == 0:
        return math.nan

    if encoding == 'numbers':
        pass
    elif encoding == 'sanger':
        qvals = numpy.array([ord(q) - 33 for q in qvals])
    else:
        raise RuntimeError("invalid `encoding`: {0}".format(encoding))

    if no_avg:
        return 1 - 10**(qvals / -10)
    else:
        return (1 - 10**(qvals / -10)).sum() / len(qvals)


def summarizeCCSreports(ccslist, report_type, plotfile,
                        plotminfrac=0.005):
    """Summarize and plot `CCS` reports.

    Args:
        `ccslist` (`CCS` object or list of them)
            `CCS` objects to summarize
        `report_type` (str "zmw" or "subread")
            Which type of report to summarize
        `plotfile` (str or `None`)
            Name of created bar plot, or `None`
            if you want to return the created plot.
        `plotminfrac` (float)
            Only plot status categories with >=
            this fraction in at least one `CCS`

    Returns:

        - If `plotfile` is a str, returns a pandas DataFrame
          aggregating the reports and creates `plotfile`.

        - If `plotfile` is `None`, returns the 2-tuple
          containing the data frame and the plot.
    """
    if isinstance(ccslist, CCS):
        ccslist = [ccslist]
    assert all([isinstance(ccs, CCS) for ccs in ccslist]), \
            "`ccslist` not a list of `CCS` objects"

    assert report_type in ['zmw', 'subread']
    report = report_type + '_report'

    df = (pandas.concat([getattr(ccs, report).assign(sample=ccs.samplename)
                for ccs in ccslist])
          .sort_values(['sample', 'number'], ascending=False)
          [['sample', 'status', 'number', 'fraction']]
          )

    # version of df that only has categories with `plotminfrac`
    plot_df = (df.assign(maxfrac=lambda x: x.groupby('status')
                         .fraction.transform('max'))
                 .query('maxfrac >= @plotminfrac')
                 )
    nstatus = len(plot_df.status.unique())

    p = (ggplot(plot_df) +
            geom_col(aes(x='sample', y='number', fill='status'),
                     position='stack') +
            theme(axis_text_x=element_text(angle=90, vjust=1,
                  hjust=0.5)) +
            ylab({'zmw':'ZMWs', 'subread':'subreads'}[report_type])
            )
    
    if nstatus <= len(COLOR_BLIND_PALETTE):
        p = p + scale_fill_manual(list(reversed(
                COLOR_BLIND_PALETTE[ : nstatus])))

    if plotfile is None:
        return (df, p)
    else:
        p.save(plotfile, 
               height=3,
               width=(2 + 0.3 * len(ccslist)),
               verbose=False)
        plt.close()
        return df

def re_expandIUPAC(re_str):
    """Expand IUPAC ambiguous nucleotide codes in `re` search string.

    Simplifies writing `re` search strings that include ambiguous
    nucleotide codes.

    Args:
        `re_str` (str)
            String appropriate to be passed to `regex.compile`.

    Returns:
        A version of `re_str` where any characters not in the group
        names that correspond to upper-case ambiguous nucleotide codes
        are expanded according to their definitions in the
        `IUPAC code <https://en.wikipedia.org/wiki/Nucleic_acid_notation>`_.

    >>> re_str = '^(?P<termini5>ATG)(?P<cDNA>N+)A+(?P<barcode>N{4})$'
    >>> re_expandIUPAC(re_str)
    '^(?P<termini5>ATG)(?P<cDNA>[ACGT]+)A+(?P<barcode>[ACGT]{4})$'
    """
    # We simply do a simple replacement on all characters not in group
    # names. So first we must find group names:
    groupname_indices = set([])
    groupname_matcher = regex.compile(r'\(\?P<[^>]*>')
    for m in groupname_matcher.finditer(re_str):
        for i in range(m.start(), m.end()):
            groupname_indices.add(i)
    
    # now replace ambiguous characters
    new_re_str = []
    for i, c in enumerate(re_str):
        if (i not in groupname_indices) and c in dms_tools2.NT_TO_REGEXP:
            new_re_str.append(dms_tools2.NT_TO_REGEXP[c])
        else:
            new_re_str.append(c)

    return ''.join(new_re_str)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
