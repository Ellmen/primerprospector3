#!/usr/bin/env python
#created 11-23-09
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
from primerprospector.cogentutil.misc import remove_files, create_dir, get_random_directory_name
from cogent3 import DNA

from primerprospector.generate_linkers import calc_base_freqs,\
 generate_linkers, get_linker_site_seqs, get_linkers, get_linkers_fp,\
 sort_linker_freqs
from primerprospector.parse import build_fasta_data



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
        
        self.small_fasta_file = small_fasta_file
        
        self.input_fasta_fp = get_tmp_filename(\
         prefix = "test_seqs_", suffix = ".fasta")
        fasta_file = open(self.input_fasta_fp, "w")
        fasta_file.write(self.small_fasta_file)
        fasta_file.close()
        
        self.expected_f_linkers_weighted_score =\
         expected_f_linkers_weighted_score
        self.expected_r_linkers_weighted_score =\
         expected_r_linkers_weighted_score
        self.expected_f_linkers_weighted_score_three_bases =\
         expected_f_linkers_weighted_score_three_bases
        self.expected_f_linkers_overall_mismatches =\
         expected_f_linkers_overall_mismatches
        self.expected_f_linkers_tp_mismatches =\
         expected_f_linkers_tp_mismatches

        
        self._files_to_remove = [self.input_hits_f_fp, self.input_hits_r_fp,
         self.input_fasta_fp]
        
    def tearDown(self):
        remove_files(self._files_to_remove)
        if exists(self.main_output_dir):
            rmtree(self.main_output_dir)
        
        
    def test_calc_base_freqs(self):
        """ Returns counts for nts at each linker position correctly """
        
        linker_len = 2
        linker_site_seqs = ["AT", "AT", "AG", "CT", "YG"]
        
        # Should not count degenerate nucleotides towards frequencies
        expected_base_freqs = [{'A':0.60, 'T':0, 'C':0.20, 'G':0}, 
         {'A':0, 'T':0.6, 'C':0, 'G':0.4}]
         
        actual_base_freqs = calc_base_freqs(linker_site_seqs, linker_len)
        
        self.assertEqual(actual_base_freqs, expected_base_freqs)
        
    def test_generate_linkers_forward(self):
        """ Handles forward linkers correctly """
        
        # Test forward hits file
        hits = self.input_hits_f_fp
        fasta_data = build_fasta_data([self.input_fasta_fp])
        linker_len = 2
        output_dir = self.main_output_dir
        score_type = "weighted_score"
        score_threshold = 1
        
        generate_linkers(hits, fasta_data, linker_len, output_dir, score_type,
         score_threshold)
         
        result_fp = output_dir + '100f_hits_suggested_linkers.txt'
        
        result_f = open(result_fp, "U")
        
        actual_result = "\n".join(line.strip() for line in result_f)
        
        self.assertEqual(actual_result, self.expected_f_linkers_weighted_score)
        
        
    def test_generate_linkers_reverse(self):
        """ Handles reverse linkers correctly """
        
        # Test reverse hits file
        hits = self.input_hits_r_fp
        fasta_data = build_fasta_data([self.input_fasta_fp])
        linker_len = 2
        output_dir = self.main_output_dir
        score_type = "weighted_score"
        score_threshold = 1
        
        generate_linkers(hits, fasta_data, linker_len, output_dir, score_type,
         score_threshold)
         
        result_fp = output_dir + '206r_hits_suggested_linkers.txt'
        
        result_f = open(result_fp, "U")
        
        actual_result = "\n".join(line.strip() for line in result_f)
        
        self.assertEqual(actual_result, self.expected_r_linkers_weighted_score)
        
    def test_generate_linkers_forward_three_bases(self):
        """ Handles forward linkers correctly with 3 base pair linker"""
        
        # Test forward hits file
        hits = self.input_hits_f_fp
        fasta_data = build_fasta_data([self.input_fasta_fp])
        linker_len = 3
        output_dir = self.main_output_dir
        score_type = "weighted_score"
        score_threshold = 1
        
        generate_linkers(hits, fasta_data, linker_len, output_dir, score_type,
         score_threshold)
         
        result_fp = output_dir + '100f_hits_suggested_linkers.txt'
        
        result_f = open(result_fp, "U")
        
        actual_result = "\n".join(line.strip() for line in result_f)
        
        self.assertEqual(actual_result, \
         self.expected_f_linkers_weighted_score_three_bases)
         
    def test_generate_linkers_forward_overall_mismatches(self):
        """ Handles forward linkers correctly """
        
        # Test forward hits file
        hits = self.input_hits_f_fp
        fasta_data = build_fasta_data([self.input_fasta_fp])
        linker_len = 2
        output_dir = self.main_output_dir
        score_type = "overall_mismatches"
        score_threshold = 1
        
        generate_linkers(hits, fasta_data, linker_len, output_dir, score_type,
         score_threshold)
         
        result_fp = output_dir + '100f_hits_suggested_linkers.txt'
        
        result_f = open(result_fp, "U")
        
        actual_result = "\n".join(line.strip() for line in result_f)
        
        self.assertEqual(actual_result, 
         self.expected_f_linkers_overall_mismatches)
        
    def test_generate_linkers_forward_tp_mismatches(self):
        """ Handles forward linkers correctly """
        
        # Test forward hits file
        hits = self.input_hits_f_fp
        fasta_data = build_fasta_data([self.input_fasta_fp])
        linker_len = 2
        output_dir = self.main_output_dir
        score_type = "tp_mismatches"
        score_threshold = 2
        
        generate_linkers(hits, fasta_data, linker_len, output_dir, score_type,
         score_threshold)
         
        result_fp = output_dir + '100f_hits_suggested_linkers.txt'
        
        result_f = open(result_fp, "U")
        
        actual_result = "\n".join(line.strip() for line in result_f)
        
        self.assertEqual(actual_result, self.expected_f_linkers_tp_mismatches)
        
    def test_get_linker_site_seqs(self):
        """ Slices out linker region sequences correctly """
        
        hits_data = ["AB329763,GTGCC,GTGCC,5,1,2,False,0,0,2.4,False",
                     "AJ937875,GTGTC,GTGCC,5,1,0,False,0,0,0.4,False",
                     "DQ243733,GTGTC,GTGCC,4,1,0,False,0,0,0.4,False"]
        fasta_data = {"AB329763":"ATAATCGATTAACCT", 
                      "AJ937875":"CCATGATCCGAAACTAC",
                      "DQ243733":"CTATCGAAACCTAGCG"}
        primer_direction = "f"
        linker_len = 2
        
        expected_linker_site_seqs = ["AT", "TG", "AT"]
        
        actual_linker_site_seqs = get_linker_site_seqs(hits_data, fasta_data,
         primer_direction, linker_len)
         
        self.assertEqual(actual_linker_site_seqs, expected_linker_site_seqs)
        
        # Test for correct reverse linkers
        
        primer_direction = "r"
        
        expected_linker_site_seqs = ["TT", "TT", "GG"]
        
        actual_linker_site_seqs = get_linker_site_seqs(hits_data, fasta_data,
         primer_direction, linker_len)
         
        self.assertEqual(actual_linker_site_seqs, expected_linker_site_seqs)
        
    def test_get_linkers_single_file(self):
        """ Gets linkers handles single file correctly """
        
        # Test with a single hits filepath
        hits_fps = self.input_hits_f_fp
        fasta_fps = self.input_fasta_fp
        linker_len = 2
        all_files = False
        output_dir = self.main_output_dir
        score_type = "overall_mismatches"
        score_threshold = 1
        
        get_linkers(hits_fps, fasta_fps, linker_len, all_files, output_dir,
         score_type, score_threshold)
         
        result_fp = output_dir + '100f_hits_suggested_linkers.txt'
        
        result_f = open(result_fp, "U")
        
        actual_result = "\n".join(line.strip() for line in result_f)
        
        self.assertEqual(actual_result, 
         self.expected_f_linkers_overall_mismatches)
        
    def test_get_linkers_multiple_files(self):
        """ Gets linkers handles multiple files correctly """
        
        # Test for multiple input hits files
        # Point to directory where hits files are
        hits_fps = self.main_output_dir


        fasta_fps = self.input_fasta_fp
        linker_len = 2
        all_files = True
        output_dir = self.main_output_dir
        score_type = "weighted_score"
        score_threshold = 1
        
        get_linkers(hits_fps, fasta_fps, linker_len, all_files, output_dir,
         score_type, score_threshold)
        
        # Test forward hits file result
        result_fp = output_dir + '100f_hits_suggested_linkers.txt'
        
        result_f = open(result_fp, "U")
        
        actual_result = "\n".join(line.strip() for line in result_f)
        
        self.assertEqual(actual_result, self.expected_f_linkers_weighted_score)
        
        # Test reverse hits file result
        result_fp = output_dir + '206r_hits_suggested_linkers.txt'
        
        result_f = open(result_fp, "U")
        
        actual_result = "\n".join(line.strip() for line in result_f)
        
        self.assertEqual(actual_result, self.expected_r_linkers_weighted_score)
        
    def test_get_linkers_fp(self):
        """ Returns output filenames for suggested linkers correctly """
        
        # No file IO, just string manipulation
        output_dir = "../linkers"
        hits_fp = "hits_dir/107f_hits.txt"
        
        expected_linkers_fp = "../linkers/107f_hits_suggested_linkers.txt"
        
        actual_linkers_fp = get_linkers_fp(output_dir, hits_fp)
        
        self.assertEqual(actual_linkers_fp, expected_linkers_fp)
        
    def test_sort_linker_freqs(self):
        """ Generates string of linker base freqs, best/worst linkers """
        
        base_freqs = [{'A':0.60, 'T':0, 'C':0.20, 'G':0}, 
                      {'A':0, 'T':0.6, 'C':0, 'G':0.4}]
                      
        expected_linker_summary = 'Base position 0\nA: 0.600 T: 0.000 C: 0.200 G: 0.000 \nBase position 1\nA: 0.000 T: 0.600 C: 0.000 G: 0.400 \n'
        expected_best_seq = "TA"
        expected_worst_seq = "AT"
        
        actual_linker_summary, actual_best_seq, actual_worst_seq =\
         sort_linker_freqs(base_freqs)
         
        self.assertEqual(actual_linker_summary, expected_linker_summary)
        self.assertEqual(actual_best_seq, expected_best_seq)
        self.assertEqual(actual_worst_seq, expected_worst_seq)


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
AJ937875,GTGTC,GTGCC,5,1,0,False,0,0,0.4,False
DQ243733,GTGTC,GTGCC,4,1,0,False,0,0,0.4,False
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
AB329763,ATCGG,ATCAG,10,1,2,False,0,0,2.4,False
AJ937875,ATCAG,ATCAG,10,0,0,False,0,0,0.0,False
DQ243733,ATCAG,ATCAG,10,3,0,False,0,0,1.2,False
DQ278134,ATCAG,ATCAG,10,0,0,False,0,0,0.0,False
DQ243735,ATCAG,ATCAG,10,0,0,False,0,0,0.0,False
AB213058,ATCAG,ATCAG,10,0,0,False,0,0,0.0,False"""

small_fasta_file = """>AB329763 1 1516 Archaea/Euryarchaeota/Halobacteriales/uncultured
AGACCATTGCTATCCGAGTTCTA
>AJ937875 1 1523 Archaea/Crenarchaeota/uncultured/uncultured
TCTTTGTTGATCCTGCCGGACCCGA
>DQ243733 1 1537 Archaea/Crenarchaeota/uncultured
AATTCTGCCGAGGGCGTGG
>DQ278134 1 1647 Archaea/Crenarchaeota/uncultured/uncultured
ATCCGCATGCCCAGGGGACGTGGA
>DQ243735 1 1589 Archaea/Crenarchaeota/uncultured
GGGGAAATTCCCCGGGGCGTTTTCCGGTTGA
>AB213058 1 1549 Archaea/Crenarchaeota/uncultured/uncultured
ATCGAGATCCTGCCGGACCC
"""

expected_f_linkers_weighted_score = """# Summary data for suggested linkers
# Note-degenerate bases are ignored for linker results, as are linkers that
# exceed the length of the input fasta sequence.  Base position starts at the
# 5' position of the linker, so position 0 would be the 5' most base in the
# linker, position 1 would be the next base, and so on.
Hits file used to generate linkers: 100f_hits.txt
Score type used: weighted_score
Score threshold: 1

Base position 0
A: 0.000 T: 0.500 C: 0.250 G: 0.250
Base position 1
A: 0.250 T: 0.500 C: 0.000 G: 0.250

Suggested linker sequence:
AC
Worst linker sequence:
TT"""

expected_f_linkers_overall_mismatches = """# Summary data for suggested linkers
# Note-degenerate bases are ignored for linker results, as are linkers that
# exceed the length of the input fasta sequence.  Base position starts at the
# 5' position of the linker, so position 0 would be the 5' most base in the
# linker, position 1 would be the next base, and so on.
Hits file used to generate linkers: 100f_hits.txt
Score type used: overall_mismatches
Score threshold: 1

Base position 0
A: 0.000 T: 1.000 C: 0.000 G: 0.000
Base position 1
A: 0.000 T: 1.000 C: 0.000 G: 0.000

Suggested linker sequence:
AA
Worst linker sequence:
TT"""

expected_f_linkers_tp_mismatches = """# Summary data for suggested linkers
# Note-degenerate bases are ignored for linker results, as are linkers that
# exceed the length of the input fasta sequence.  Base position starts at the
# 5' position of the linker, so position 0 would be the 5' most base in the
# linker, position 1 would be the next base, and so on.
Hits file used to generate linkers: 100f_hits.txt
Score type used: tp_mismatches
Score threshold: 2

Base position 0
A: 0.000 T: 0.500 C: 0.500 G: 0.000
Base position 1
A: 0.000 T: 0.500 C: 0.250 G: 0.250

Suggested linker sequence:
AA
Worst linker sequence:
CT"""

expected_f_linkers_weighted_score_three_bases = """# Summary data for suggested linkers
# Note-degenerate bases are ignored for linker results, as are linkers that
# exceed the length of the input fasta sequence.  Base position starts at the
# 5' position of the linker, so position 0 would be the 5' most base in the
# linker, position 1 would be the next base, and so on.
Hits file used to generate linkers: 100f_hits.txt
Score type used: weighted_score
Score threshold: 1

Base position 0
A: 0.250 T: 0.250 C: 0.250 G: 0.250
Base position 1
A: 0.000 T: 0.500 C: 0.250 G: 0.250
Base position 2
A: 0.250 T: 0.500 C: 0.000 G: 0.250

Suggested linker sequence:
AAC
Worst linker sequence:
ATT"""

expected_r_linkers_weighted_score = """# Summary data for suggested linkers
# Note-degenerate bases are ignored for linker results, as are linkers that
# exceed the length of the input fasta sequence.  Base position starts at the
# 5' position of the linker, so position 0 would be the 5' most base in the
# linker, position 1 would be the next base, and so on.
Hits file used to generate linkers: 206r_hits.txt
Score type used: weighted_score
Score threshold: 1

Base position 0
A: 0.000 T: 0.250 C: 0.500 G: 0.250
Base position 1
A: 0.000 T: 0.000 C: 0.750 G: 0.250

Suggested linker sequence:
AA
Worst linker sequence:
CC"""

if __name__ == "__main__":
    main()
