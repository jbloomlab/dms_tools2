"""Tests building barcode variant tables.

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

import dms_tools2


class test_build_barcode_variant_table(unittest.TestCase):
    
    def test_build_barcode_variant_table(self):

        indir = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            'build_barcode_variant_table_files/')
        outdir = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                              'test_build_barcode_variant_table/')
        os.makedirs(outdir, exist_ok=True)

        lookup_table = os.path.join(outdir, 'barcode_variant_table.csv')

        processed_file = os.path.join(indir, 'processed_CCSs.csv')

        df_processed = (
            pandas.read_csv(processed_file)
            .assign(mutations=lambda x: x.mutations.apply(
                dms_tools2.minimap2.Mutations.from_str))
            )
    
    
        min_seq_acc = 0.999

        df_processed = (
            df_processed
            .assign(accurate=lambda x: 
                (x.barcode_accuracy >= min_seq_acc) & 
                (x.gene_accuracy >= min_seq_acc))
            )


        mut_acc_df = (
            df_processed
            .query('accurate')
            .assign(insertion=lambda x: x.mutations.apply(
                        dms_tools2.minimap2.Mutations.insertions,
                                returnval='accuracy'),
                    deletion=lambda x: x.mutations.apply(
                        dms_tools2.minimap2.Mutations.deletions,
                                returnval='accuracy'),
                    substitution=lambda x: x.mutations.apply(
                        dms_tools2.minimap2.Mutations.substitutions,
                                returnval='accuracy'))
            .melt(value_vars=['insertion', 'deletion', 'substitution'],
                  var_name='mutation_type', value_name='accuracy')
            .groupby('mutation_type')
            .accuracy
            .apply(itertools.chain.from_iterable)
            .apply(list)
            .apply(pandas.Series)
            .stack()
            .rename('accuracy')
            .reset_index()
            .assign(error_rate=lambda x: numpy.clip(1 - x.accuracy, 1e-6, None))
            )

        mut_acc_df_file = os.path.join(indir, 'mut_acc_df.csv')
        pandas.testing.assert_frame_equal(
                mut_acc_df,
                pandas.read_csv(mut_acc_df_file)
                )

        min_mut_acc = 0.999

        df_processed = (
            df_processed
            .assign(
                substitutions=lambda x: x.mutations.apply(
                    dms_tools2.minimap2.Mutations.substitutions,
                    min_acc=min_mut_acc, min_acc_filter_nan=False),
                insertions=lambda x: x.mutations.apply(
                    dms_tools2.minimap2.Mutations.insertions,
                    min_acc=min_mut_acc, min_acc_filter_nan=False),
                deletions=lambda x: x.mutations.apply(
                    dms_tools2.minimap2.Mutations.deletions,
                    min_acc=min_mut_acc, min_acc_filter_nan=False),
                substitutions_str=lambda x: x.substitutions.str.join(', '),
                deletions_str=lambda x: x.deletions.str.join(', '),
                insertions_str=lambda x: x.insertions.str.join(', '))
            )

        df_processed = (
            df_processed
            .assign(has_indel=lambda x: 
                (x.insertions.apply(len) + x.deletions.apply(len)).astype('bool'))
            )

        consensus, dropped = dms_tools2.barcodes.simpleConsensus(
                df_processed.query('accurate'), library_col='library')

        consensus = (
            consensus
            .assign(
                n_substitutions=lambda x: x.substitutions.apply(len),
                n_insertions=lambda x: x.insertions.apply(len),
                n_deletions=lambda x: x.deletions.apply(len),
                has_indel=lambda x: (x.insertions + x.deletions).astype('bool')
                )
            )

        barcode_lookup = (
            consensus
            .query('not has_indel')
            .assign(substitutions=lambda x: x.substitutions.str.join(' '))
            [['library', 'barcode', 'substitutions', 'variant_call_support']]
            )

        barcode_lookup_file = os.path.join(indir, 'barcode_lookup.csv')
        pandas.testing.assert_frame_equal(
                barcode_lookup.reset_index(drop=True),
                pandas.read_csv(barcode_lookup_file).fillna('')
                )


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
