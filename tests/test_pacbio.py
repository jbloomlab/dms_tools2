"""Tests core functions of `dms_tools2.pacbio`.

Based off some small slice of the data from the 
analysis at https://github.com/jbloomlab/IFNsorted_flu_single_cell"""


import os
import unittest

import pandas
from pandas.testing import assert_frame_equal

import dms_tools2.utils
import dms_tools2.pacbio
import dms_tools2.minimap2
import dms_tools2.seqnumbering


class test_pacbio(unittest.TestCase):
    """Tests `dms_tools2.pacbio.matchAndAlignCCS`"""

    def test_matchAndAlignCCS(self):
        """Tests `dms_tools2.pacbio.matchAndAlignCCS`."""

        cwd = os.path.abspath(os.path.dirname(__file__))
        indir = f'{cwd}/pacbio_files'

        ccslist = [
                dms_tools2.pacbio.CCS(
                    name,
                    f'{indir}/{name}.bam',
                    None)
                for name in ['CCSs-1', 'CCSs-2']]

        mapper = dms_tools2.minimap2.Mapper(
                f'{indir}/targets.fasta',
                dms_tools2.minimap2.OPTIONS_VIRUS_W_DEL,
                target_isoforms={'fluM1':['fluM2'], 'fluM2':['fluM1'],
                                 'fluNS1':['fluNS2'], 'fluNS2':['fluNS1']}
                )

        targetvariants = dms_tools2.minimap2.TargetVariants(
                {'wildtype':f'{indir}/targets.fasta',
                 'synonymously barcoded':f'{indir}/targets-double-syn.fasta'},
                mapper,
                variantsites_min_acc=0.99
                )

        transcriptconverter = dms_tools2.seqnumbering.TranscriptConverter(
                f'{indir}/flu-wsn.gb',
                ignore_other_features=True
                )

        mutationcaller = dms_tools2.minimap2.MutationCaller(
                mapper,
                transcriptconverter=transcriptconverter,
                query_softclip=10,
                target_clip=20
                )

        primer3 = 'CTACACGACGCTCTTCCGATCT'
        trim_termini = 5
        primer5_mix = {'fluPB2':'GCGAAAGCAGGTCAATTATATTCAATATGGAAAG',
                       'fluPB1':'GCGAAAGCAGGCAAACCATTTG',
                       'fluPA' :'GCGAAAGCAGGTACTGATTCAAAATGG',
                       'fluHA' :'GCAAAAGCAGGGGAAAATAAAAACAACC',
                       'fluNP' :'GCAAAAGCAGGGTAGATAATCACTCAC',
                       'fluNA' :'GCGAAAGCAGGAGTTTAAATGAATCCAAAC',
                       'fluM'  :'GCAAAAGCAGGTAGATATTGAAAGATGAGTC',
                       'fluNS' :'GCAAAAGCAGGGTGACAAAGACATAATG',
                       }

        df_ccs = dms_tools2.pacbio.matchAndAlignCCS(
                ccslist=ccslist,
                mapper=mapper,
                termini5='|'.join([s[trim_termini : ] for s in primer5_mix.values()]),
                gene='N+B',
                spacer='AAA(A{19,}){e<=2}AAA',
                umi='N{10}',
                barcode='N{16}',
                termini3=dms_tools2.utils.reverseComplement(primer3)[ : -trim_termini],
                targetvariants=targetvariants,
                mutationcaller=mutationcaller
                ).rename(columns={'has_spacer':'has_polyA'})

        # test alignment stats
        align_stats_aligned = (
                df_ccs
                .assign(n=1)
                .groupby(
                    ['barcoded', 'barcoded_polarity', 'has_termini3',
                     'has_termini5', 'has_polyA', 'gene_aligned',
                     'CCS_aligned'])
                .aggregate({'n':'sum'})
                .reset_index()
                )
        expected = pandas.read_csv(f'{indir}/align_stats_aligned.csv',
                keep_default_na=False)
        assert_frame_equal(align_stats_aligned, expected)
        align_stats_gene = (
                df_ccs
                .query('gene_aligned')
                .assign(n=1)
                .groupby('gene_aligned_target')
                .aggregate({'n':'sum'})
                .reset_index()
                )
        expected = pandas.read_csv(f'{indir}/align_stats_gene.csv',
                keep_default_na=False)
        assert_frame_equal(align_stats_gene, expected)

        # check mutations and alignments on aligned reads
        alignment_info = (
                df_ccs
                .query('gene_aligned')
                .assign(
                    cigar=lambda x:
                        x.gene_aligned_alignment
                         .apply(lambda a: a.cigar_str),
                    substitutions=lambda x:
                        x.gene_aligned_mutations
                         .apply(lambda m: ' '.join(m.substitutions())),
                    insertions=lambda x:
                        x.gene_aligned_mutations
                         .apply(lambda m: ' '.join(m.insertions())),
                    deletions=lambda x:
                        x.gene_aligned_mutations
                         .apply(lambda m: ' '.join(m.deletions()))
                    )
                [['samplename',
                  'name',
                  'gene_aligned_target',
                  'gene_aligned_n_trimmed_query_start',
                  'gene_aligned_n_trimmed_query_end',
                  'gene_aligned_n_trimmed_target_start',
                  'gene_aligned_n_trimmed_target_end',
                  'gene_aligned_n_additional',
                  'gene_aligned_n_additional_difftarget',
                  'gene_aligned_target_variant',
                  'cigar',
                  'substitutions',
                  'insertions',
                  'deletions']]
                .reset_index(drop=True)
                )
        expected_alignment_info = pandas.read_csv(
                f'{indir}/alignment_info.csv',
                keep_default_na=False)
        assert_frame_equal(alignment_info, expected_alignment_info)



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
