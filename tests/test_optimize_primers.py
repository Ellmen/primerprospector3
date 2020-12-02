#!/usr/bin/env python
from __future__ import division

__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

from os.path import basename, isdir, exists
from shutil import rmtree

from cogent3.util.unit_test import TestCase, main
from primerprospector.old_cogent import get_tmp_filename
from cogent3.util.misc import remove_files, create_dir, get_random_directory_name
from cogent3 import DNA

from primerprospector.optimize_primers import generate_primers_pos_data,\
 filter_uneven_sequence_hits, reverse_sequence_hits, get_base_freqs,\
 get_primer_base_freqs_fp




class GenerateLinkersTests(TestCase):
    """unit tests for generate_linkers """
    
    def setUp(self):
        """ Generate temporary files for use in unit testing"""
        
        self.main_output_dir = get_random_directory_name(prefix = '/tmp/')
        self.main_output_dir += '/'
        
        create_dir(self.main_output_dir)
        
        self.small_forward_hits_file = small_forward_hits_file
        
        self.input_hits_f_fp = \
         self.main_output_dir +'100f_hits.txt'
        hits_file = open(self.input_hits_f_fp, 'w')
        hits_file.write(self.small_forward_hits_file)
        hits_file.close()
        
        self.small_reverse_hits_file = small_reverse_hits_file
        
        self.input_hits_r_fp = \
         self.main_output_dir + '206r_hits.txt'
        hits_file = open(self.input_hits_r_fp, 'w')
        hits_file.write(self.small_reverse_hits_file)
        hits_file.close()
        
        self._files_to_remove = [self.input_hits_f_fp, self.input_hits_r_fp]
        
        self.expected_f_default_results = expected_f_default_results
        
        self.expected_r_results = expected_r_results
        
    def tearDown(self):
        remove_files(self._files_to_remove)
        if exists(self.main_output_dir):
            rmtree(self.main_output_dir)
            
    def test_generate_primers_pos_data_default_settings(self):
        """ Overall module functionality test with default settings """
        
        # Test using default settings
        
        generate_primers_pos_data(self.input_hits_f_fp,
         output_dir = self.main_output_dir)
         
        actual_results_fp = self.main_output_dir +\
         "100f_hits_base_frequencies.txt"
         
        actual_results_f = open(actual_results_fp, "U")
        
        actual_results = actual_results_f.read()
        
        self.assertEqual(actual_results, self.expected_f_default_results)
        
        
    def test_generate_primers_pos_data_altered_scoring(self):
        """ Overall module functionality handles different settings """
        
        # Test with reverse primer hits, altered score settings
        
        generate_primers_pos_data(self.input_hits_r_fp,
         output_dir = self.main_output_dir, score_type = "overall_mismatches",
         score_threshold = 3)
         
        actual_results_fp = self.main_output_dir +\
         "206r_hits_base_frequencies.txt"
         
        actual_results_f = open(actual_results_fp, "U")
        
        actual_results = actual_results_f.read()
        
        self.assertEqual(actual_results, self.expected_r_results)
        
    def test_filter_uneven_sequence_hits(self):
        """ Removes sequences that do not match the primer length """
        
        
        primer = "ACGACG"
                
        # Should not remove any sequences as they all are equal in len to primer
        seqs = ['ACCACC', 'ATTACG', 'ACAGAG']
        
        expected_seqs = seqs
        
        actual_seqs = filter_uneven_sequence_hits(seqs, primer)
        
        self.assertEqual(actual_seqs, expected_seqs)
        
        # Should remove the two sequences that are not the same size
        seqs = ['CCACC', 'AATTACG', 'ACAGAG']
        
        expected_seqs = ['ACAGAG']
        
        actual_seqs = filter_uneven_sequence_hits(seqs, primer)
        
        self.assertEqual(actual_seqs, expected_seqs)
        
        
        
    def test_reverse_sequence_hits(self):
        """ Reverse complements sequences correctly """
        
        seqs = ['CCACC', 'AATTACG', 'ACAGAG']
        
        expected_seqs = ['GGTGG', 'CGTAATT', 'CTCTGT']
        
        actual_seqs = reverse_sequence_hits(seqs)
        
        self.assertEqual(actual_seqs, expected_seqs)
        
    def test_get_base_freqs(self):
        """ Properly sorts base frequencies to a list of dicts """
        
        primer = "ACGAGC"
        
        seqs = ['ACCACC', 'ATTACG', 'ACAGAG', 'ACAGAG']
        
        expected_base_freqs = [{'A': 1.0, 'C': 0.0, 'T': 0.0, 'G': 0.0},
                               {'A': 0.0, 'C': 0.75, 'T': 0.25, 'G': 0.0},
                               {'A': 0.5, 'C': 0.25, 'T': 0.25, 'G': 0.0},
                               {'A': 0.5, 'C': 0.0, 'T': 0.0, 'G': 0.5},
                               {'A': 0.5, 'C': 0.5, 'T': 0.0, 'G': 0.0},
                               {'A': 0.0, 'C': 0.25, 'T': 0.0, 'G': 0.75}]
        
        actual_base_freqs = get_base_freqs(seqs, primer)
        
        self.assertEqual(actual_base_freqs, expected_base_freqs)
        
        
    def test_get_primer_base_freqs_fp(self):
        """ Returns output filepath correctly """
        
        output_dir = "/test/"
        
        hits_fp = "../100f_hits.txt"
        
        expected_output_fp = "/test/100f_hits_base_frequencies.txt"
        
        actual_fp = get_primer_base_freqs_fp(output_dir, hits_fp)
        
        self.assertEqual(actual_fp, expected_output_fp)
        

    
# Large strings placed at the ends for better readability

small_forward_hits_file = """# Primer: 100f 5'-GTGCC-3'
# Input fasta file: archaeal_v15.fasta
# Parameters
# 3' length: 5
# non 3' mismatch penalty: 0.40 per mismatch
# 3' mismatch penalty: 1.00 per mismatch
# last base mismatch penalty: 3.00
# non 3' gap penalty: 1.00 per gap
# 3' gap penalty: 3.00 per gap
# Note - seq hit and primer hit are the best local pairwise alignment results for a given sequence and primer pair.  A gap in seq hit represents a deletion in the sequence, whereas a gap in the primer hit signifies an insertion in the target sequence.
#
# seq ID, seq hit, primer hit, hit start position, non 3' mismatches, 3' mismatches (except last base), last base mismatch, non 3' gaps, 3' gaps, overall weighted score, hits sequence end 
AB329763,GTGCC,GTGCC,5,1,2,False,0,0,2.4,False
AJ937875,GTGTG,GTGCC,5,1,0,False,0,0,0.4,False
DQ243733,GTGTG,GTGCC,4,1,0,False,0,0,0.4,False
DQ278134,GTGTC,GTGCC,5,2,0,False,0,0,0.4,False
DQ243735,GTGTC,GTGCC,5,1,3,False,0,0,0.4,False
AB213058,GTGTC,GT-CC,5,1,0,False,0,0,0.4,False"""

small_reverse_hits_file = """# Primer: 206r 5'-GGACT-3'
# Input fasta file: archaeal_v15.fasta
# Parameters
# 3' length: 5
# non 3' mismatch penalty: 0.40 per mismatch
# 3' mismatch penalty: 1.00 per mismatch
# last base mismatch penalty: 3.00
# non 3' gap penalty: 1.00 per gap
# 3' gap penalty: 3.00 per gap
# Note - seq hit and primer hit are the best local pairwise alignment results for a given sequence and primer pair.  A gap in seq hit represents a deletion in the sequence, whereas a gap in the primer hit signifies an insertion in the target sequence.
#
# seq ID, seq hit, primer hit, hit start position, non 3' mismatches, 3' mismatches (except last base), last base mismatch, non 3' gaps, 3' gaps, overall weighted score, hits sequence end 
AB329763,ATCGG,ATCAG,10,0,3,False,0,0,2.4,False
AJ937875,CACAG,ATCAG,10,0,0,False,0,0,0.0,False
DQ243733,ATCAG,ATCAG,10,3,1,False,0,0,1.2,False
DQ278134,ATCAG,ATCAG,10,0,0,False,0,0,0.0,False
DQ243735,ATCAG,ATCAG,10,0,0,False,0,0,0.0,False
AB213058,ATCAG,ATCAG,10,0,0,False,0,0,0.0,False"""


expected_f_default_results = """# Base frequency report for optimizing primers
# This file is tab separated for easy importation into Excel or other spreadsheets
# Degenerate DNA codes (listed here for convenience): R=AG, Y=CT, M=AC, K=GT, W=AT, S=CG, B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT
# The primer is listed in 5' to 3' orientation.
# Hits file used to generate base frequency data: 100f_hits.txt
# Score type used: weighted_score
# Score threshold: 2
# Primer sequence: GTGCC
#
Primer\tG\tT\tG\tC\tC
Base
A\t0.00\t0.00\t0.00\t0.00\t0.00
T\t0.00\t1.00\t0.00\t1.00\t0.00
C\t0.00\t0.00\t0.00\t0.00\t0.60
G\t1.00\t0.00\t1.00\t0.00\t0.40
"""

expected_r_results = """# Base frequency report for optimizing primers
# This file is tab separated for easy importation into Excel or other spreadsheets
# Degenerate DNA codes (listed here for convenience): R=AG, Y=CT, M=AC, K=GT, W=AT, S=CG, B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT
# The primer is listed in 5' to 3' orientation.
# Hits file used to generate base frequency data: 206r_hits.txt
# Score type used: overall_mismatches
# Score threshold: 3
# Primer sequence: GGACT
#
Primer\tG\tG\tA\tC\tT
Base
A\t0.00\t0.00\t0.00\t0.80\t0.00
T\t0.00\t0.80\t0.00\t0.20\t0.80
C\t1.00\t0.20\t0.00\t0.00\t0.00
G\t0.00\t0.00\t1.00\t0.00\t0.20
"""


if __name__ == "__main__":
    main()
