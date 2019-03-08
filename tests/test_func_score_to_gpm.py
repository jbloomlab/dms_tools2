"""Tests func_score_to_gpm function in `codonvarianttable.py`.
 
Checks if a previously produced gpm matches a gpm produced
using the same code. Will only check if the code produces 
the same gpm, not if it is working 'correctly'"""

import os
import unittest
import tempfile

import pandas as pd

import dms_tools2.codonvarianttable
import gpmap

class test_func_score_to_gpm(unittest.TestCase):
    """Tests func_score_to_gpm"""
    
    def test_produced_gpm(self):
        """Tests func_score_to_gpm"""
        
        self.indir = os.path.join(
                    os.path.abspath(os.path.dirname(__file__)),
                    'func_score_to_gpm_input_files/')

        geneseq = 'ATGAGGTGA'
        barcode_variants = {}
        barcode_variants['library'] = ['library-1', 'library-1', 'library-1']
        barcode_variants['barcode'] = ['AAA', 'ATT', 'AAT']
        barcode_variants['substitutions'] = ['', 'G3A', 'A9C']
        barcode_variants['variant_call_support'] = [1, 1, 1]
        barcode_variants = pd.DataFrame(barcode_variants)
        with tempfile.NamedTemporaryFile(mode='w') as f:
            barcode_variants.to_csv(f, index=False)
            f.flush()
            variants = dms_tools2.codonvarianttable.CodonVariantTable(
                        barcode_variant_file=f.name,
                        geneseq=geneseq)
        phenosimulator = dms_tools2.codonvarianttable.PhenotypeSimulator(geneseq)
        counts_df = []
        post_samples = {'post_selection':{
                        'total_count':50,
                        'noise':0,
                        'bottleneck':50,
                         }}
        post_to_pre = {}
        post_to_pre['post_selection'] = 'pre_selection'
        pre_sample_name = 'pre_selection'
        counts_df.append(
                dms_tools2.codonvarianttable.simulateSampleCounts(
                        variants=variants,
                        phenotype_func=phenosimulator.observedPhenotype,
                        variant_error_rate=0,
                        pre_sample={'total_count':50,
                                    'uniformity':10},
                        pre_sample_name=pre_sample_name,
                        post_samples=post_samples
                        )
                 )
        counts_df = pd.concat(counts_df)
        pd.DataFrame.from_records(list(post_to_pre.items()),
                                columns=['post-selection', 'pre-selection'])
        for sample in counts_df['sample'].unique():
            icounts = counts_df.query('sample == @sample')
            variants.addSampleCounts('library-1', sample, icounts)
        func_scores = variants.func_scores(post_to_pre)
        library = 'library-1'
        for sample in func_scores.post_sample.unique():
            gpm = dms_tools2.codonvarianttable.func_score_to_gpm(
                    func_scores.query(
                        'library == @library & post_sample == @sample'
                        ), 
                    variants.geneseq)
        prior_gpm = gpmap.GenotypePhenotypeMap.read_csv(
                        os.path.join(
                            self.indir,
                            'func_score_to_gpm_test.csv'
                            ),
                        wildtype = 'MR*'
                        )
        self.assertTrue(gpm.data.round(decimals=10).equals(
                            prior_gpm.data.round(decimals=10)
                            )
                        )

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
