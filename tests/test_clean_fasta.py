#!/usr/bin/env python
# File created on 8 Sep 2010

__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

from os.path import basename
from os import remove

from cogent3.util.unit_test import TestCase, main
from primerprospector.old_cogent import get_tmp_filename
from cogent3.util.misc import remove_files

from primerprospector.clean_fasta import filter_fasta




class FilterFastaFileTests(TestCase):
    """ Unit tests for the filter_fasta.py module """
    
    def setUp(self):
        # create the temporary input files that will be used with the
        # filter_fasta function
        self.input_aligned_fasta1 = input_aligned_fasta1
        
        self.input_aligned_fasta1_fp = get_tmp_filename(\
         prefix = 'input_aligned_fasta1',
         suffix = '.fasta')
        seq_file = open(self.input_aligned_fasta1_fp, 'w')
        seq_file.write(self.input_aligned_fasta1)
        seq_file.close()
        
        self.input_fasta2 = input_fasta2
        
        self.input_fasta2_fp = get_tmp_filename(\
         prefix = 'input_fasta2',
         suffix = '.fasta')
        seq_file = open(self.input_fasta2_fp, 'w')
        seq_file.write(self.input_fasta2)
        seq_file.close()
        
        # Get output filenames, are based on input tmp file names
        self.expected_output_fasta1 = '/tmp/' +\
         basename(self.input_aligned_fasta1_fp).split('.')[0] +\
         '_filtered.fasta'
        self.expected_output_fasta2 = '/tmp/' +\
         basename(self.input_fasta2_fp).split('.')[0] +\
         '_filtered.fasta'
         

        self.expected_output_combined1 = expected_output_combined1
        self.expected_output_combined2 = expected_output_combined2
        self.expected_output_disabled_filtering = input_aligned_fasta1
        
        self.expected_output_fasta1_lines = expected_output_fasta1_lines
        
        self._files_to_remove =\
         [self.input_aligned_fasta1_fp, self.input_fasta2_fp]
                 
    def tearDown(self):
        remove_files(self._files_to_remove)
        
    def test_filter_fasta_single_file(self):
        """ Properly filters given fasta file """
        
        self._files_to_remove.append(self.expected_output_fasta1)
        
        
        filter_fasta(self.input_aligned_fasta1_fp, output_dir = '/tmp/')
        
        actual_result = []
        
        f = open(self.expected_output_fasta1, "U")
        for line in f:
            actual_result.append(line.strip())
            
        self.assertEqual(actual_result, self.expected_output_fasta1_lines)
        
    def test_filter_fasta_multi_files(self):
        """ Properly filters given fasta files """
        
        self._files_to_remove.append(self.expected_output_fasta1)
        self._files_to_remove.append(self.expected_output_fasta2)
        
        filter_fasta(self.input_aligned_fasta1_fp + ":" + self.input_fasta2_fp,\
         output_dir = '/tmp/')
        
        actual_result = []
        
        f = open(self.expected_output_fasta1, "U")
        for line in f:
            actual_result.append(line.strip())
            
        self.assertEqual(actual_result, self.expected_output_combined1)
        
        actual_result = []
        
        f = open(self.expected_output_fasta2, "U")
        for line in f:
            actual_result.append(line.strip())
            
        self.assertEqual(actual_result, self.expected_output_combined2)
        
    def test_filter_fasta_disabled_filtering(self):
        """ Skips filtering options if flagged """
        
        
        self._files_to_remove.append(self.expected_output_fasta1)
        
        
        filter_fasta(self.input_aligned_fasta1_fp, output_dir = '/tmp/',
         gap_chars = False, space_chars = False, convert_uracil = False)
        
        actual_result = []
        
        f = open(self.input_aligned_fasta1_fp, "U")
        for line in f:
            actual_result.append(line.strip())
            
        self.assertEqual("".join(actual_result),
         self.expected_output_disabled_filtering.replace('\n', ''))
        
        
        
        
        

        
# Place large strings at end of file for better readability

# Simple aligned fasta with both forms of gap chars, U char.
input_aligned_fasta1 = "\n".join(""">seq1
----  --AATC-GGU--...
>seq2
---C--AATC-G-U--...""".split('\n'))

# expected filter_fasta result of the input_aligned_fasta1 file.
expected_output_fasta1_lines = """>seq1
AATCGGT
>seq2
CAATCGT""".split('\n')

input_fasta2 = "\n".join(""">seq3
AATCTTTACG
>seq4
---CT---CG""".split('\n'))

expected_output_combined1 = """>seq1
AATCGGT
>seq2
CAATCGT""".split('\n')

expected_output_combined2 = """>seq3
AATCTTTACG
>seq4
CTCG""".split('\n')
        
        
if __name__ == "__main__":
    main()
