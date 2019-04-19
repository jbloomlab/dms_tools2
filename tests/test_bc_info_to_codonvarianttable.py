import unittest
import os

import pandas as pd

from dms_tools2.utils import bc_info_to_codonvarianttable

class test_bc_info_to_codonvarianttable(unittest.TestCase):
    """Tests bc_info_to_codonvarianttable
    
    Only test if the function is writing an output identical to a 
    previous output, not that the function is working as expected
    in every way.
    """

    def test_produced_codonvarianttable(self):
        """Tests bc_info_to_codonvarianttable"""
        
        self.indir = os.path.join(
                    os.path.abspath(os.path.dirname(__file__)),
                    'bc_info_to_codonvarianttable_input_files/')

        samples = {'library-1':['test']}
        geneseq = 'ATGTCTAAGAAACCAGGAGGGCCCGGCAAAAGCCGGGCTGTCAATATGCTAAAACGCGGAATGCCCCGCGTGTTGTCCTTAATTGGACTGAAGAGGGCTATGCTGAGCCTGATCGACGGTAGGGGGCCAATACGGTTTGTGTTGGCTCTCTTGGCGTTTTTTAGGTTCACGGCAATTGCTCCGACCCGGGCAGTGCTGGATCGATGGAGAAGTGTGAACAAACAAACAGCGATGAAACACCTCCTGAGTTTCAAGAAGGAACTAGGAACCTTGACCAGCGCTATCAACCGGCGGAGTTCAAAACAGAAG'
        variants = bc_info_to_codonvarianttable(samples,
                                                geneseq,
                                                path=self.indir
                                               )
        variants.writeCodonCounts(single_or_all='all', outdir=self.indir)
        test = pd.read_csv(os.path.join(
                               self.indir,
                               'library-1_test_codoncounts.csv'
                               )
                          )
        previous = pd.read_csv(os.path.join(
                                  self.indir,
                                  'previous_codoncounts.csv'
                                  )
                              )
        
        self.assertTrue(test.equals(previous))
        
if __name__ == '__main__':
    runner = unitteset.TextTestRunner()
    unittest.main(testRunner=runner)

