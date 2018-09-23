"""Profiles `dms_tools2.pacbio.matchAndAlignCCS`."""


import os
import cProfile
import pstats

import pandas

import dms_tools2.utils
import dms_tools2.pacbio
import dms_tools2.minimap2
import dms_tools2.seqnumbering


def run_matchAndAlignCCS():
    """Runs `dms_tools2.pacbio.matchAndAlignCCS`."""

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


if __name__ == '__main__':
    statsfile = 'pacbio_pstats'
    cProfile.run('run_matchAndAlignCCS()', statsfile)
    p = pstats.Stats(statsfile)
    for t in ['cumtime', 'tottime']:
        print(f'\n{t}:')
        p.strip_dirs().sort_stats(t).print_stats(10)
