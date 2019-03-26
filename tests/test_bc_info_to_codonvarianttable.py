import unittest

import pandas as pd

from dms_tools2.codonvarianttable import bc_info_to_codonvarianttable

class test_bc_info_to_codonvarianttable(unittest.TestCase):
    """Tests bc_info_to_codonvarianttable"""

    def test_produced_codonvarianttable(self):
        """Tests bc_info_to_codonvarianttable"""

        samples = ['test']
        geneseq = 'ATGTCTAAGAAACCAGGAGGGCCCGGCAAAAGCCGGGCTGTCAATATGCTAAAACGCGGAATGCCCCGCGTGTTGTCCTTAATTGGACTGAAGAGGGCTATGCTGAGCCTGATCGACGGTAGGGGGCCAATACGGTTTGTGTTGGCTCTCTTGGCGTTTTTTAGGTTCACGGCAATTGCTCCGACCCGGGCAGTGCTGGATCGATGGAGAAGTGTGAACAAACAAACAGCGATGAAACACCTCCTGAGTTTCAAGAAGGAACTAGGAACCTTGACCAGCGCTATCAACCGGCGGAGTTCAAAACAGAAG'
        variants = bc_info_to_codonvarianttable(samples,
                                                geneseq,
                                                path='bc_info_to_codonvarianttable_input_files/',
                                               )
        variants.writeCodonCounts(single_or_all='all', outdir='bc_info_to_codonvarianttable_input_files/')
        test = pd.read_csv('bc_info_to_codonvarianttable_input_files/library-1_test_codoncounts.csv')
        previous = pd.read_csv('bc_info_to_codonvarianttable_input_files/previous_codoncounts.csv')
        
        self.assertTrue(test.equals(previous))
        
if __name__ == '__main__':
    runner = unitteset.TextTestRunner()
    unittest.main(testRunner=runner)

