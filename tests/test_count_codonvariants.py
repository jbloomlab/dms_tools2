"""Tests counting codon variants from `CodonVariantTable`.

Runs on a snippet of real data for RecA.

This test doesn't test for correct output, but rather
simply that output is the same as a version of the
code that we are confident is correct. Therefore,
if later changes break this code, this test will identify
that fact.

Written by Jesse Bloom."""

import os
import re
import collections
import itertools
import unittest

import numpy
import pandas
import pandas.testing
import Bio.SeqIO

import dms_tools2


class test_count_codonvariants(unittest.TestCase):
    
    def test_count_codonvariants(self):

        indir = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                             'count_codonvariant_files/')
        outdir = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                              'test_count_codonvariants_files/')
        os.makedirs(outdir, exist_ok=True)

        barcode_variant_file = os.path.join(indir, 'barcode_variant_table.csv')

        ampliconfile = os.path.join(indir, 'PacBio_amplicon.gb')
        amplicon = Bio.SeqIO.read(ampliconfile, 'genbank')

        gene = [f for f in amplicon.features if f.type == 'gene']
        assert len(gene) == 1, "Failed to find exactly one gene feature"
        geneseq = str(gene[0].location.extract(amplicon).seq)


        variants = dms_tools2.barcodes.CodonVariantTable(
                        barcode_variant_file=barcode_variant_file,
                        geneseq=geneseq)

        n_variants_file = os.path.join(indir, 'n_variants.csv')
        pandas.testing.assert_frame_equal(
                variants.n_variants_df(samples=None)
                        .assign(library=lambda x: x['library'].astype('str'),
                                sample=lambda x: x['sample'].astype('str')
                                ),
                pandas.read_csv(n_variants_file),
                check_dtype=False
                )

        for mut_type in ['aa', 'codon']:
            p = variants.plotNumMutsHistogram(mut_type, samples=None, max_muts=7)

        for variant_type in ['all', 'single']:
            p = variants.plotNumCodonMutsByType(variant_type, samples=None,)

        p = variants.plotNumCodonMutsByType('all', samples=None, min_support=2)


        for variant_type, mut_type in itertools.product(['all', 'single'],
                                                        ['aa', 'codon']):
            p = variants.plotCumulMutCoverage(variant_type,
                                              mut_type,
                                              samples=None)

        for variant_type, mut_type in itertools.product(['all', 'single'],
                                                        ['aa', 'codon']):
            p = variants.plotMutFreqs(variant_type, mut_type, samples=None)

        for mut_type in ['aa', 'codon']:
            p = variants.plotMutHeatmap('all', mut_type,
                                samples=None, libraries='all_only',
                                widthscale=2)


        samples = (pandas.DataFrame({
                        'library':['library-1', 'library-2', 'library-1'],
                        'sample':['plasmid', 'plasmid', 'uninduced'],
                        })
                   .assign(
                        run='run-1',
                        upstream='TTTTTAAGTTGTAAGGATATGCCATTCTAGA',
                        downstream='',
                        R1file=lambda x: indir + x['library'] + '_' +
                                         x['sample'] + '_R1.fastq',
                        R2file=lambda x: indir + x['library'] + '_' +
                                         x['sample'] + '_R2.fastq',
                        )     
                   ) 

        barcode = [f for f in amplicon.features if f.type == 'barcode']
        assert len(barcode) == 1, "Failed to find exactly one barcode feature"
        bclen = len(barcode[0])

        fates = []
        for (lib, sample), runs in samples.groupby(['library', 'sample']):
    
            # read barcodes for all runs for library / sample
            barcodes = []
            for run_tup in runs.itertuples(index=False):
                parser = dms_tools2.barcodes.IlluminaBarcodeParser(
                            bclen=bclen,
                            upstream=run_tup.upstream,
                            downstream=run_tup.downstream,
                            valid_barcodes=variants.valid_barcodes(lib),
                            upstream_mismatch=3,
                            downstream_mismatch=3)
                run_barcodes, run_fates = parser.parse(
                                    r1files=run_tup.R1file,
                                    r2files=None)
                barcodes.append(run_barcodes)
                fates.append(run_fates.assign(library=lib,
                                              sample=sample,
                                              run=run_tup.run
                                              ))
    
            # combine barcodes read for each run
            barcodes = (pandas.concat(barcodes,
                              ignore_index=True,
                              sort=False)
                         .groupby('barcode')
                         .aggregate({'count':'sum'})
                         .reset_index()
                         )
    
            # add barcode counts to variant table
            variants.addSampleCounts(lib, sample, barcodes)
    
        # concatenate read fates into one data frame
        fates = (pandas.concat(fates, ignore_index=True, sort=False)
                  .assign(count=lambda x: x['count'].fillna(0).astype('int'))
                  )

        fatesfile = os.path.join(indir, 'fates.csv')
        pandas.testing.assert_frame_equal(
                fates,
                pandas.read_csv(fatesfile)
                )

        libs_to_analyze = ['library-1']

        for mut_type in ['aa', 'codon']:
            p = variants.plotNumMutsHistogram(mut_type,
                        libraries=libs_to_analyze, max_muts=7,
                        orientation='v')

        for variant_type in ['all', 'single']:
            p = variants.plotNumCodonMutsByType(variant_type,
                                        libraries=libs_to_analyze,
                                        orientation='v')

        for variant_type, mut_type in itertools.product(['all', 'single'],
                                                        ['aa', 'codon']):
            p = variants.plotCumulMutCoverage(variant_type, mut_type,
                                      libraries=libs_to_analyze,
                                      orientation='v')

        for variant_type, mut_type in itertools.product(['all', 'single'],
                                                        ['aa', 'codon']):
            p = variants.plotMutFreqs(variant_type, mut_type,
                                libraries=libs_to_analyze,
                                orientation='v')

        for mut_type in ['aa', 'codon']:
            p = variants.plotMutHeatmap('all', mut_type,
                                    libraries=libs_to_analyze,
                                    widthscale=2)

        variant_count_file = os.path.join(indir, 'variant_counts.csv')
        pandas.testing.assert_frame_equal(
            variants.variant_count_df
                    .assign(library=lambda x: x['library'].astype(str),
                            sample=lambda x: x['sample'].astype(str)
                            ),
            pandas.read_csv(variant_count_file).fillna('')
            )

        countfiles = []
        for mut_type in ['single', 'all']:
            countsdir = os.path.join(outdir, f'{mut_type}_mutant_codoncounts')
            mut_type_countfiles = variants.writeCodonCounts(mut_type,
                                    outdir=countsdir)
            countfiles.append(mut_type_countfiles
                              .assign(mutant_type=f'{mut_type} mutants')
                              )


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
