#!/usr/bin/env python
#created 7-28-09

from __future__ import division

__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

from os.path import isdir, exists
from os import remove, makedirs
from shutil import rmtree

from cogent3.core.alphabet import AlphabetError
from cogent3 import DNA, LoadSeqs
from numpy import array, bitwise_or
from cogent3.util.unit_test import TestCase, main
from primerprospector.old_cogent import get_tmp_filename
from cogent3.util.misc import remove_files, get_random_directory_name

from primerprospector.sort_denovo_primers import Primer, KnownPrimer,\
 convert_to_numeric, build_primers, find_overlap, find_3prime_match,\
 find_matches, find_overlap, parse_primer_data, get_known_primers,\
 sort_primers_data, convert_to_DNA, build_known_primers, analyze_primers,\
 check_for_std_alignment




 

class SortDenovoPrimersTests(TestCase):
    """ Unit tests for the sort_denovo_primers.py module """
    def setUp(self):
        """ Initializes variables for unit tests"""
        
        self.sample_hits_data = sample_hits_data
        self.known_primers_no_overlap_data = known_primers_no_overlap_data
        self.known_primers_data = known_primers_data
        self.bad_known_primers_data = bad_known_primers_data
    
    def test_primer_filter_base_freq(self):
        """ method returns proper data"""
        # Should error with invalid data
        self.assertRaises(TypeError, Primer.filter_base_freq)
        self.assertRaises(TypeError, Primer.filter_base_freq, "string")
        self.assertRaises(TypeError, Primer.filter_base_freq, "string",
         0.1)
        bad_list=['a', 'b', 'c']
        self.assertRaises(TypeError, Primer.filter_base_freq,
         bad_list, 0.1)
        # The default threshold to delete bases with low occurrances
        # is 0.10.  There should be no deletions with the frequences
        # in the following instance of primer.
        expected_result=[{'G':1.0}, {'A':0.2, 'G':0.8}]
        small_seqs=LoadSeqs(data=['>Seq1', 'GA', '>Seq2', 'GG',
         '>Seq3', 'GG', '>Seq4', 'GG', '>Seq5', 'GG'], moltype=DNA)
        small_primer_instance=Primer("GA,51,2041,9,90.0%,10.0%,102,137",
         small_seqs, small_seqs)
        self.assertEqual(small_primer_instance.filtered_f_base_freq,
         expected_result)
        # By default, this method should delete the 'A' dictionary
        # item that is occurring at 0.20 frequency with a
        # variable_pos_freq parameter of 0.25
        expected_result=[{'G':1.0}, {'G':0.8}]
        small_seqs=LoadSeqs(data=['>Seq1', 'GA', '>Seq2', 'GG',
         '>Seq3', 'GG', '>Seq4', 'GG', '>Seq5', 'GG'], moltype=DNA)
        small_primer_instance=Primer("GA,51,2041,9,90.0%,10.0%,102,137",
         small_seqs, small_seqs, 0.25)
        self.assertEqual(small_primer_instance.filtered_f_base_freq,
         expected_result)
        
    def test_primer_get_filtered_degenerate_seq(self):
        """ Properly filters degenerate sequences """
        # method raises errors with bad data 
        self.assertRaises(TypeError, Primer.get_filtered_degenerate_seq)
        self.assertRaises(TypeError, Primer.get_filtered_degenerate_seq,
         "string")
        bad_list=['X', 'O']
        self.assertRaises(TypeError, Primer.get_filtered_degenerate_seq,
         bad_list)
        # method returns proper results with valid data
        # With the default retention of characters, should get a
        # degenerate resulting string for the last 3 bp
        expected_result="GRKS"
        small_seqs=LoadSeqs(data=['>Seq1', 'GATC', '>Seq2', 'GGGG',
         '>Seq3', 'GGGG', '>Seq4', 'GGGG', '>Seq5', 'GGGG'], moltype=DNA)
        small_primer_instance=Primer("GATC,51,2041,9,90.0%,10.0%,102,137",
         small_seqs, small_seqs)
        self.assertEqual(small_primer_instance.f_degenerate_seq,
         expected_result)
        # With a more stringent threshold for variable position 
        # retention, should get no degeneracy in result.
        expected_result="GGGG"
        small_seqs=LoadSeqs(data=['>Seq1', 'GATC', '>Seq2', 'GGGG',
         '>Seq3', 'GGGG', '>Seq4', 'GGGG', '>Seq5', 'GGGG'], moltype=DNA)
        small_primer_instance=Primer("GATC,51,2041,9,90.0%,10.0%,102,137",
         small_seqs, small_seqs, 0.25)
        self.assertEqual(small_primer_instance.f_degenerate_seq,
         expected_result)
        # Test for all degenerate codes working properly
        # Identity test
        expected_result="ATTCGNRYMKWSBDHV-A"
        small_seqs=LoadSeqs(data=['>Seq1', 'ATUCGNRYMKWSBDHV--',
         '>Seq2', 'ATUCGNRYMKWSBDHV-A',], moltype=DNA)
        small_primer_instance=Primer("GATC,51,2041,9,90.0%,10.0%,102,137",\
         small_seqs, small_seqs)
        self.assertEqual(small_primer_instance.f_degenerate_seq,\
         expected_result)
        # Complementation test
        expected_result="NNNNNNNN"
        small_seqs=LoadSeqs(data=['>Seq1', 'ACGTUWMR', '>Seq2', 'BDHVVSKY',],
         moltype=DNA)
        small_primer_instance=Primer("GATC,51,2041,9,90.0%,10.0%,102,137",
         small_seqs, small_seqs)
        self.assertEqual(small_primer_instance.f_degenerate_seq,
         expected_result)
         
    def test_primer_orientations(self):
        """ Generates numeric versions of primer in all 4 orientations """
        seq_all_bases = LoadSeqs(data=['>Seq1', 'ATUCGNRYMKWSBDHV-'],
         moltype=DNA)
        test_primer = Primer("GATC,51,2041,9,90.0%,10.0%,102,137",
         seq_all_bases, seq_all_bases)
        expected_f_seq = \
         [1, 2, 2, 4, 8, 15, 9, 6, 5, 10, 3, 12, 14, 11, 7, 13, 0]
        expected_f_complement = \
         [2, 1, 1, 8, 4, 15, 6, 9, 10, 5, 3, 12, 13, 7, 11, 14, 0]
        expected_f_reverse = \
         [0, 13, 7, 11, 14, 12, 3, 10, 5, 6, 9, 15, 8, 4, 2, 2, 1]
        expected_f_rc = \
         [0, 14, 11, 7, 13, 12, 3, 5, 10, 9, 6, 15, 4, 8, 1, 1, 2]
        expected_r_seq =\
         [0, 14, 11, 7, 13, 12, 3, 5, 10, 9, 6, 15, 4, 8, 1, 1, 2]
        expected_r_complement =\
         [0, 13, 7, 11, 14, 12, 3, 10, 5, 6, 9, 15, 8, 4, 2, 2, 1]
        expected_r_reverse =\
         [2, 1, 1, 8, 4, 15, 6, 9, 10, 5, 3, 12, 13, 7, 11, 14, 0]
        expected_r_rc =\
         [1, 2, 2, 4, 8, 15, 9, 6, 5, 10, 3, 12, 14, 11, 7, 13, 0]
        # Test to make sure all orientations are generated properly
        self.assertEqual(test_primer.f_seq, expected_f_seq)
        self.assertEqual(test_primer.f_complement, expected_f_complement)
        self.assertEqual(test_primer.f_reverse, expected_f_reverse)
        self.assertEqual(test_primer.f_rc, expected_f_rc)
        self.assertEqual(test_primer.r_seq, expected_r_seq)
        self.assertEqual(test_primer.r_complement, expected_r_complement)
        self.assertEqual(test_primer.r_reverse, expected_r_reverse)
        self.assertEqual(test_primer.r_rc, expected_r_rc)
        
    def test_primer_name(self):
        """ Names primer properly according to supplied indices """
        
        small_seqs = LoadSeqs(data=['>Seq1','ATUCGNRYMKWSBDHV--',
         '>Seq2', 'ATUCGNRYMKWSBDHV-A',], moltype=DNA)
        test_primer=Primer("GATC,51,2041,9,90.0%,10.0%",\
         small_seqs, small_seqs)
        # As no standard aligned index was supplied, should name primers
        # according to the first unaligned index, 51
        expected_primer_f_name = "51f"
        expected_primer_r_name = "51r"
        self.assertEqual(test_primer.primer_f_name, expected_primer_f_name)
        self.assertEqual(test_primer.primer_r_name, expected_primer_r_name)
        
        # Same data, but with supplied standard alignment values
        small_seqs = LoadSeqs(data=['>Seq1','ATUCGNRYMKWSBDHV--',\
         '>Seq2','ATUCGNRYMKWSBDHV-A',], moltype=DNA)
        test_primer=Primer("GATC,51,2041,9,90.0%,10.0%,102,137",\
         small_seqs, small_seqs)
        expected_primer_f_name = "102f"
        expected_primer_r_name = "137r"
        self.assertEqual(test_primer.primer_f_name, expected_primer_f_name)
        self.assertEqual(test_primer.primer_r_name, expected_primer_r_name)

    def test_convert_to_DNA(self):
        """ properly converts numeric codes to DNA """
        bad_seq=[2, 1, 'A', 4]
        self.assertRaises(KeyError, convert_to_DNA,bad_seq)
        expected_result="ATTCGNRYMKWSBDHV"
        numeric_list=[1,2,2,4,8,15,9,6,5,10,3,12,14,11,7,13]
        self.assertEqual(convert_to_DNA(numeric_list), expected_result)
        
    def test_convert_to_numeric(self):
        """ convert_to_numeric raises errors if invalid data passed """
        bad_primer="AATXOB"
        self.assertRaises(KeyError, convert_to_numeric, bad_primer)
        """ convert_to_numeric converts IUPAC characters to numeric list.
        These lists are used for bitwise comparisons """
        sample_seq="ATUCGNRYMKWSBDHV"
        expected_list=[1,2,2,4,8,15,9,6,5,10,3,12,14,11,7,13]
        self.assertEqual(convert_to_numeric(sample_seq), expected_list)
        
    def test_build_primers(self):
        """ Properly builds primers from list of lines of data """
        
        # See sample_hits_data at the end of the file
        test_variable_pos_freq = 0.20
        test_primer_name = "New_Primer"
        test_cmp_truncate_len = 5
        
        # first forward primer
        expected_f_primer1_result = ["New_Primer103f", "RSMDDMMBGCTCWGTAACAC",
         "GGCRDACGGCTCAGTAACAC", "GGCATACGGCTCAGTAACAC", 90.0, 10.0,
         [1, 1, 4, 1, 4]]

        # second reverse primer
        expected_r_primer2_result = ["New_Primer137r", "GKYADRTTRHCYRCGTGTTA",
         "GGYAKRTTRHCYACGTGTTA", "GGCAGGTTAACTACGTGTTA", 90.0, 10.0,
         [2, 8, 2, 2, 1]]

        test_primers = build_primers(self.sample_hits_data,
         test_variable_pos_freq, test_primer_name, test_cmp_truncate_len)

        f_primer_result = [test_primers[0].primer_f_name,
         test_primers[0].f_IUPAC, test_primers[0].f_degenerate_seq,
         test_primers[0].f_consensus, test_primers[0].sensitivity,
         test_primers[0].specificity, test_primers[0].f_trunc_seq]
        r_primer_result = [test_primers[1].primer_r_name,
         test_primers[1].r_IUPAC, test_primers[1].r_degenerate_seq,
         test_primers[1].r_consensus, test_primers[1].sensitivity,
         test_primers[1].specificity, test_primers[1].r_trunc_seq]
         
        self.assertEqual(f_primer_result, expected_f_primer1_result)
        self.assertEqual(r_primer_result, expected_r_primer2_result)
        
        
    def test_find_3prime_match(self):
        """ Properly finds 3' match between denovo and known primers """
        
        test_variable_pos_freq = 0.20
        test_primer_name = ""
        test_cmp_truncate_len = 5
        test_primers = build_primers(self.sample_hits_data,
         test_variable_pos_freq, test_primer_name, test_cmp_truncate_len)
        
        # Partial overlap, no 3' match
        test_known_primer_f = "GGCATACGGCTCAGT"
        # overlap and 3' match
        test_known_primer_r = "GTTAACTACGTGTTA"
        
        trunc_len = 5
        test_known_primer = KnownPrimer(test_known_primer_f, "overlap_f",
         trunc_len)
        for p in test_primers:
            find_3prime_match(p, test_known_primer, trunc_len)
            # f_unique value should remain True as no 3' matches are found
            self.assertEqual(p.f_unique, True)
        
        test_known_primer = KnownPrimer(test_known_primer_r, "match_r",
         trunc_len)
        # First reverse primer does not have a perfect 3' match, but the 
        # second reverse primer does
        for p in test_primers:
            find_3prime_match(p, test_known_primer, trunc_len)
            
        self.assertEqual(test_primers[0].r_unique, True)
        self.assertEqual(test_primers[1].r_unique, False)
        
        expected_3prime_match =\
         ['137r\tGGYAKRTTRHCYACGTGTTA\tmatch_r\tGTTAACTACGTGTTA\n']
        self.assertEqual(test_primers[1].r_3prime_match, expected_3prime_match)

        
    def test_find_matches(self):
        """ Properly finds overlap, 3' matches between primers """
        
        test_variable_pos_freq = 0.20
        test_primer_name = ""
        test_cmp_truncate_len = 5
        test_match_len = 10
        test_primers = build_primers(self.sample_hits_data,
         test_variable_pos_freq, test_primer_name, test_cmp_truncate_len)
         
         
        test_primers = find_matches(test_primers, self.known_primers_data,
         test_match_len, test_cmp_truncate_len)
         
        # Should get partial overlap for all primers, and a 3' match for the 
        # second reverse primer
        expected_overlap_results = \
         [False, ['103f\tGGCRDACGGCTCAGTAACAC\ttest_known_primer_f\tGGCATACGGCTCAGT\tGGCATACGGC\n'],
         [], False,
         ['138r\tGGGYAKRTTRHCYACGTGTT\ttest_known_primer_r\tGTTAACTACGTGTTA\tGTTAACTACG\n'],
         [], False,
         ['102f\tHGGCRDACGGCTCAGTAACA\ttest_known_primer_f\tGGCATACGGCTCAGT\tGGCATACGGC\n'],
         [], False,
         ['137r\tGGYAKRTTRHCYACGTGTTA\ttest_known_primer_r\tGTTAACTACGTGTTA\tGTTAACTACG\n'],
         ['137r\tGGYAKRTTRHCYACGTGTTA\ttest_known_primer_r\tGTTAACTACGTGTTA\n']]
        actual_overlap_results = \
         [test_primers[0].f_unique, test_primers[0].f_partial_overlap,
          test_primers[0].f_3prime_match, test_primers[0].r_unique,
          test_primers[0].r_partial_overlap, test_primers[0].r_3prime_match,
          test_primers[1].f_unique, test_primers[1].f_partial_overlap,
          test_primers[1].f_3prime_match, test_primers[1].r_unique,
          test_primers[1].r_partial_overlap, test_primers[1].r_3prime_match]
        self.assertEqual(actual_overlap_results, expected_overlap_results)
        
        # Test again, should get no overlap with known primers that do not
        # share sequences
        test_primers = build_primers(self.sample_hits_data,
         test_variable_pos_freq, test_primer_name, test_cmp_truncate_len)
        
        test_primers = find_matches(test_primers,
         self.known_primers_no_overlap_data, test_match_len,
         test_cmp_truncate_len)
        
        expected_overlap_results = \
         [True, [], [], True, [], [], True, [], [], True, [], []]
        actual_overlap_results = \
         [test_primers[0].f_unique, test_primers[0].f_partial_overlap,
          test_primers[0].f_3prime_match, test_primers[0].r_unique,
          test_primers[0].r_partial_overlap, test_primers[0].r_3prime_match,
          test_primers[1].f_unique, test_primers[1].f_partial_overlap,
          test_primers[1].f_3prime_match, test_primers[1].r_unique,
          test_primers[1].r_partial_overlap, test_primers[1].r_3prime_match]
        self.assertEqual(actual_overlap_results, expected_overlap_results)

        
    def test_find_overlap(self):
        """ Finds overlap between two primers"""
        trunc_len=5
        primer1=KnownPrimer("AATCGACG", "primer1", trunc_len)
        primer2=KnownPrimer("AATCGACG", "primer2", trunc_len)
        match_len=6
        expected_result="AATCGA"
        self.assertEqual(find_overlap(match_len,primer1.numeric_seq,
         primer2.numeric_seq),expected_result)
        # Change the first base of one of the primers, should shift result
        # one base to the right
        primer2=KnownPrimer("CATCGACG","primer2",trunc_len)
        expected_result="ATCGAC"
        self.assertEqual(find_overlap(match_len,primer1.numeric_seq,
         primer2.numeric_seq),expected_result)
        # Primers that don't have 6 (match_len value) consecutive matching
        # base pairs should return False result
        primer2=KnownPrimer("CATCCACG","primer2",trunc_len)
        expected_result=False
        self.assertEqual(find_overlap(match_len,primer1.numeric_seq,
         primer2.numeric_seq),expected_result)
        # Test for degeneracy values
        primer2=KnownPrimer("RWNSVDMG","primer2",trunc_len)
        expected_result="AATCGA"
        self.assertEqual(find_overlap(match_len,primer1.numeric_seq,
         primer2.numeric_seq),expected_result)
        
    def test_build_known_primers(self):
        """ Properly builds list of KnownPrimer objects for known primers """
        
        test_cmp_truncate_len = 5
        
        test_known_primers = get_known_primers(self.known_primers_data)
        test_known_primers = build_known_primers(test_known_primers,
         test_cmp_truncate_len)
        expected_primer_name = '1f'
        expected_primer_seq = "ATTCCGGTTGATCCTGC"
        expected_primer_trunc_seq = [4, 4, 2, 8, 4]
        expected_primer_num_seq = \
         [1, 2, 2, 4, 4, 8, 8, 2, 2, 8, 1, 2, 4, 4, 2, 8, 4]
        
        self.assertEqual(test_known_primers[0].name, expected_primer_name)
        self.assertEqual(test_known_primers[0].seq, expected_primer_seq)
        self.assertEqual(test_known_primers[0].trunc_seq,
         expected_primer_trunc_seq)
        self.assertEqual(test_known_primers[0].numeric_seq,
         expected_primer_num_seq)
         
        # Should raise error if incorrect format of input data
        self.assertRaises(ValueError, get_known_primers,
         self.bad_known_primers_data)
         
class OverallInputOutputTests(TestCase):
    """ Tests overall input/output functionality of module """
    
    def setUp(self):
        # create the temporary input files that will be used with the
        # search_sequences function
        
        self.expected_formatted_primers_default =\
         expected_formatted_primers_default
        self.expected_primers_details_default =\
         expected_primers_details_default
        self.expected_formatted_primers_45perc_var_pos_freq =\
         expected_formatted_primers_45perc_var_pos_freq
        self.expected_primers_details_45perc_var_pos_freq =\
         expected_primers_details_45perc_var_pos_freq
        self.expected_formatted_primers_sensitivity_sort =\
         expected_formatted_primers_sensitivity_sort
        self.expected_formatted_primers_specificity_sort =\
         expected_formatted_primers_specificity_sort
        self.expected_formatted_primers_degeneracy_sort =\
         expected_formatted_primers_degeneracy_sort
        self.expected_primers_overlap = expected_primers_overlap
        self.expected_formatted_primers_with_name = \
         expected_formatted_primers_with_name
        self.expected_primers_overlap_match_len_5 =\
         expected_primers_overlap_match_len_5
        self.expected_primers_overlap_trunc_len_1 =\
         expected_primers_overlap_trunc_len_1
        self.expected_formatted_primers_std_index =\
         expected_formatted_primers_std_index
        self.expected_primers_details_default_std_index =\
         expected_primers_details_default_std_index
        self.expected_amplicon_pairs =\
         expected_amplicon_pairs
                 
        self.sample_primer_hits_file = get_tmp_filename(\
         prefix='sample_hits_',
         suffix='.txt')
        seq_file = open(self.sample_primer_hits_file, 'w')
        seq_file.write(sample_primer_hits_file)
        seq_file.close()        
        
        self.sample_primer_hits_file_std_alignment = get_tmp_filename(\
         prefix='sample_hits_std_alignment_',
         suffix='.txt')
        seq_file = open(self.sample_primer_hits_file_std_alignment, 'w')
        seq_file.write(sample_primer_hits_file_std_alignment)
        seq_file.close()
        
        self.sample_known_primers_file = get_tmp_filename(\
         prefix='known_primers_',
         suffix='.txt')
        seq_file = open(self.sample_known_primers_file, 'w')
        seq_file.write(sample_known_primers_file)
        seq_file.close()
        
        
        self.output_dir = get_random_directory_name(\
         prefix='/tmp/sort_denovo_primers')
        self.output_dir += '/'
         
        if not exists(self.output_dir):
            make_dirs(self.output_dir)
        
        self._files_to_remove =\
         [self.sample_primer_hits_file, 
          self.sample_primer_hits_file_std_alignment, 
          self.sample_known_primers_file]
        
    def tearDown(self):
        remove_files(self._files_to_remove)
        if exists(self.output_dir):
            rmtree(self.output_dir)
        
    def test_analyze_primers_default_settings(self):
        """ Test output with default settings """
        
            
        formatted_primers_output = self.output_dir + "formatted_primers.txt"
        primers_details_output = self.output_dir + "primers_details.txt"
        
        self._files_to_remove.append(formatted_primers_output)
        self._files_to_remove.append(primers_details_output)
        
        # Use the default values
        analyze_primers(hits_file = self.sample_primer_hits_file,
                    output_dir = self.output_dir,
                    verbose = False,
                    variable_pos_freq = 0.20,
                    sort_method = 'S', 
                    known_primers_filepath = None,
                    primer_name = "",
                    match_len = 10, 
                    cmp_truncate_len = 10,
                    amplicon_len = False)
                    
        formatted_primers = open(formatted_primers_output, "U")
        primers_details = open(primers_details_output, "U")
        
        formatted_primers_result = [l.strip('\n') for l in formatted_primers]
        primers_details_result = [l.strip('\n') for l in primers_details]
            
        self.assertEqual(formatted_primers_result,
         self.expected_formatted_primers_default)
        self.assertEqual(primers_details_result,
         self.expected_primers_details_default)
                    
        
        
        
        
    def test_analyze_primers_variable_pos_freq(self):
        """ Test output with altered degeneracy filter """
        
        formatted_primers_output = self.output_dir + "formatted_primers.txt"
        primers_details_output = self.output_dir + "primers_details.txt"
        
        self._files_to_remove.append(formatted_primers_output)
        self._files_to_remove.append(primers_details_output)
        
        analyze_primers(hits_file = self.sample_primer_hits_file,
                    output_dir = self.output_dir,
                    verbose = False,
                    # Set a higher threshold, will lower degeneracy
                    variable_pos_freq = 0.45,
                    sort_method = 'S', 
                    known_primers_filepath = None,
                    primer_name = "",
                    match_len = 10, 
                    cmp_truncate_len = 10,
                    amplicon_len = False)
                    
        formatted_primers = open(formatted_primers_output, "U")
        primers_details = open(primers_details_output, "U")
        
        formatted_primers_result = [l.strip('\n') for l in formatted_primers]
        primers_details_result = [l.strip('\n') for l in primers_details]
            
        self.assertEqual(formatted_primers_result,
         self.expected_formatted_primers_45perc_var_pos_freq)
        self.assertEqual(primers_details_result,
         self.expected_primers_details_45perc_var_pos_freq)
                    
        
        
    def test_analyze_primers_sort_method(self):
        """ Test analyze primers with different sorting methods """
        
        formatted_primers_output = self.output_dir + "formatted_primers.txt"
        primers_details_output = self.output_dir + "primers_details.txt"
        
        self._files_to_remove.append(formatted_primers_output)
        self._files_to_remove.append(primers_details_output)
        
        # Start with default, sensitivity sorting
        analyze_primers(hits_file = self.sample_primer_hits_file,
                    output_dir = self.output_dir,
                    verbose = False,
                    variable_pos_freq = 0.20,
                    sort_method = 'S', 
                    known_primers_filepath = None,
                    primer_name = "",
                    match_len = 10, 
                    cmp_truncate_len = 10,
                    amplicon_len = False)
                    
        formatted_primers = open(formatted_primers_output, "U")
        
        formatted_primers_result = [l.strip('\n') for l in formatted_primers]
            
        self.assertEqual(formatted_primers_result,
         self.expected_formatted_primers_sensitivity_sort)
         
        # Specificity sorting test
        analyze_primers(hits_file = self.sample_primer_hits_file,
                    output_dir = self.output_dir,
                    verbose = False,
                    variable_pos_freq = 0.20,
                    sort_method = 'P', 
                    known_primers_filepath = None,
                    primer_name = "",
                    match_len = 10, 
                    cmp_truncate_len = 10,
                    amplicon_len = False)
                    
        formatted_primers = open(formatted_primers_output, "U")
        
        formatted_primers_result = [l.strip('\n') for l in formatted_primers]
            
        self.assertEqual(formatted_primers_result,
         self.expected_formatted_primers_specificity_sort)
         
        # Degeneracy sorting test
        analyze_primers(hits_file = self.sample_primer_hits_file,
                    output_dir = self.output_dir,
                    verbose = False,
                    variable_pos_freq = 0.20,
                    sort_method = 'O', 
                    known_primers_filepath = None,
                    primer_name = "",
                    match_len = 10, 
                    cmp_truncate_len = 10,
                    amplicon_len = False)
                    
        formatted_primers = open(formatted_primers_output, "U")
        
        formatted_primers_result = [l.strip('\n') for l in formatted_primers]
            
        self.assertEqual(formatted_primers_result,
         self.expected_formatted_primers_degeneracy_sort)

        
        
    def test_analyzie_primers_known_primers(self):
        """ Test analyze primers with comparison to known primers """
        
        formatted_primers_output = self.output_dir + "formatted_primers.txt"
        primers_details_output = self.output_dir + "primers_details.txt"
        primers_overlap_output = self.output_dir + "primers_overlap.txt"
        
        self._files_to_remove.append(formatted_primers_output)
        self._files_to_remove.append(primers_details_output)
        self._files_to_remove.append(primers_overlap_output)
        
        # Use the default values
        analyze_primers(hits_file = self.sample_primer_hits_file,
                    output_dir = self.output_dir,
                    verbose = False,
                    variable_pos_freq = 0.20,
                    sort_method = 'S', 
                    known_primers_filepath = self.sample_known_primers_file,
                    primer_name = "",
                    match_len = 10, 
                    cmp_truncate_len = 10,
                    amplicon_len = False)
                    
        primers_overlap = open(primers_overlap_output, "U")
        
        primers_overlap_result = [l.strip('\n') for l in primers_overlap]
            
        self.assertEqual(primers_overlap_result,
         self.expected_primers_overlap)

                    
        
        
    def test_analyze_primers_primer_name(self):
        """ Test analyze primers with added primer name to output """
        
        formatted_primers_output = self.output_dir + "formatted_primers.txt"
        primers_details_output = self.output_dir + "primers_details.txt"
        
        self._files_to_remove.append(formatted_primers_output)
        self._files_to_remove.append(primers_details_output)
        
        # Use the default values, with added primer prefix name
        analyze_primers(hits_file = self.sample_primer_hits_file,
                    output_dir = self.output_dir,
                    verbose = False,
                    variable_pos_freq = 0.20,
                    sort_method = 'S', 
                    known_primers_filepath = None,
                    primer_name = "CaptHammer",
                    match_len = 10, 
                    cmp_truncate_len = 10,
                    amplicon_len = False)
                    
        formatted_primers = open(formatted_primers_output, "U")
        
        formatted_primers_result = [l.strip('\n') for l in formatted_primers]
            
        self.assertEqual(formatted_primers_result,
         self.expected_formatted_primers_with_name)

        
    def test_analyze_primers_match_len(self):
        """ Test analyze primers to altered match len to known primers """
        
        formatted_primers_output = self.output_dir + "formatted_primers.txt"
        primers_details_output = self.output_dir + "primers_details.txt"
        primers_overlap_output = self.output_dir + "primers_overlap.txt"
        
        self._files_to_remove.append(formatted_primers_output)
        self._files_to_remove.append(primers_details_output)
        self._files_to_remove.append(primers_overlap_output)
        
        # Change match len to a smaller value, will result in more
        # overlap hits
        analyze_primers(hits_file = self.sample_primer_hits_file,
                    output_dir = self.output_dir,
                    verbose = False,
                    variable_pos_freq = 0.20,
                    sort_method = 'S', 
                    known_primers_filepath = self.sample_known_primers_file,
                    primer_name = "",
                    match_len = 5, 
                    cmp_truncate_len = 10,
                    amplicon_len = False)
                    
        primers_overlap = open(primers_overlap_output, "U")
        
        primers_overlap_result = [l.strip('\n') for l in primers_overlap]
            
        self.assertEqual(primers_overlap_result,
         self.expected_primers_overlap_match_len_5)
        
    def test_analyze_primers_cmp_truncate_len(self):
        """ Test analyze primers with altered match to 3' of known primers """
        
        formatted_primers_output = self.output_dir + "formatted_primers.txt"
        primers_details_output = self.output_dir + "primers_details.txt"
        primers_overlap_output = self.output_dir + "primers_overlap.txt"
        
        self._files_to_remove.append(formatted_primers_output)
        self._files_to_remove.append(primers_details_output)
        self._files_to_remove.append(primers_overlap_output)
        
        # Change trunc len to 1, will result in many more primers being
        # 'complete' matches to known primers
        analyze_primers(hits_file = self.sample_primer_hits_file,
                    output_dir = self.output_dir,
                    verbose = False,
                    variable_pos_freq = 0.20,
                    sort_method = 'S', 
                    known_primers_filepath = self.sample_known_primers_file,
                    primer_name = "",
                    match_len = 10, 
                    cmp_truncate_len = 1,
                    amplicon_len = False)
                    
        primers_overlap = open(primers_overlap_output, "U")
        
        primers_overlap_result = [l.strip('\n') for l in primers_overlap]
            
        self.assertEqual(primers_overlap_result,
         self.expected_primers_overlap_trunc_len_1)
         
    def test_analyze_primers_amplicon_len(self):
        """ Test output with amplicon_len settings """
        
            
        formatted_primers_output = self.output_dir + "formatted_primers.txt"
        primers_details_output = self.output_dir + "primers_details.txt"
        amplicon_details_output = self.output_dir + "amplicon_len_pairs.txt"
        
        self._files_to_remove.append(formatted_primers_output)
        self._files_to_remove.append(primers_details_output)
        self._files_to_remove.append(amplicon_details_output)
        
        # Use the default values, but test values with a standard alignment
        # value
        analyze_primers(hits_file = self.sample_primer_hits_file_std_alignment,
                    output_dir = self.output_dir,
                    verbose = False,
                    variable_pos_freq = 0.20,
                    sort_method = 'S', 
                    known_primers_filepath = None,
                    primer_name = "",
                    match_len = 10, 
                    cmp_truncate_len = 10,
                    amplicon_len = "200:600")
                    
        formatted_primers = open(formatted_primers_output, "U")
        primers_details = open(primers_details_output, "U")
        amplicon_details = open(amplicon_details_output, "U")
        
        formatted_primers_result = [l.strip('\n') for l in formatted_primers]
        primers_details_result = [l.strip('\n') for l in primers_details]
        amplicon_details_result = [l.strip('\n') for l in amplicon_details]
        
                    
        # Should only find one pair of primers with the estimated amplicon
        # in the range 200:600
        self.assertEqual(formatted_primers_result,
         self.expected_formatted_primers_std_index)
        self.assertEqual(primers_details_result,
         self.expected_primers_details_default_std_index)
        self.assertEqual(amplicon_details_result,
         self.expected_amplicon_pairs)
         
    def test_check_for_std_alignment(self):
        """ Raises error if standard alignment not in input data """
        
        hits_data_no_alignment = "\n".join(sample_primer_hits_file)
        
        self.assertRaises(IndexError, check_for_std_alignment, \
         hits_data_no_alignment)
        
        

        
        

# Start sample data, list format

sample_hits_data = ['C,TAACA,51,2041,9,90.0%,10.0%,102,137', 'M,2041,TGGCATACCGCTCTGTAACA,TAACACGTAGATAACCTACC,AY800210', 'M,2041,TGGCGTACGGCTCAGTAACA,TAACACGTGGATAACTTACC,EU883771', 'M,2041,TGGCATACGGCTCAGTAACA,TAACACGTAGTCAACATGCC,EF503699', 'M,2041,GGCATACAGGCTCAGTAACA,TAACACGTAGTCAACATGCC,EF503697', 'M,2041,AGGCGGACGGCTCAGTAACA,TAACACGCAGTCAACCTAAC,AJ535128', 'M,2041,TAGCATACGGCTCAGTAACA,TAACACGTAGTCAACCTGCC,EF021349', 'M,2041,CGGCGAACGGCTCAGTAACA,TAACACGTGGGTAATCTGCC,AJ969787', 'M,2041,CGGCATACGGCTCAGTAACA,TAACACGTGGGTAACCTGCC,AY360086', 'M,2041,AGGCGGACTGCTCAGTAACA,TAACACGTGGGTAATCTACC,DQ641863', 'C,AACAC,52,2043,9,90.0%,10.0%,103,138', 'M,2043,GGCATACCGCTCTGTAACAC,AACACGTAGATAACCTACCC,AY800210', 'M,2043,GGCGTACGGCTCAGTAACAC,AACACGTGGATAACTTACCC,EU883771', 'M,2043,GGCATACGGCTCAGTAACAC,AACACGTAGTCAACATGCCC,EF503699', 'M,2043,GCATACAGGCTCAGTAACAC,AACACGTAGTCAACATGCCC,EF503697', 'M,2043,GGCGGACGGCTCAGTAACAC,AACACGCAGTCAACCTAACT,AJ535128', 'M,2043,AGCATACGGCTCAGTAACAC,AACACGTAGTCAACCTGCCC,EF021349', 'M,2043,GGCGAACGGCTCAGTAACAC,AACACGTGGGTAATCTGCCC,AJ969787', 'M,2043,GGCATACGGCTCAGTAACAC,AACACGTGGGTAACCTGCCC,AY360086', 'M,2043,GGCGGACTGCTCAGTAACAC,AACACGTGGGTAATCTACCC,DQ641863']
expected_sample_f_primers = {'2043': ['AACAC,52,2043,9,90.0%,10.0%,103,138', 'GGCATACCGCTCTGTAACAC', 'GGCGTACGGCTCAGTAACAC', 'GGCATACGGCTCAGTAACAC', 'GCATACAGGCTCAGTAACAC', 'GGCGGACGGCTCAGTAACAC', 'AGCATACGGCTCAGTAACAC', 'GGCGAACGGCTCAGTAACAC', 'GGCATACGGCTCAGTAACAC', 'GGCGGACTGCTCAGTAACAC'], '2041': ['TAACA,51,2041,9,90.0%,10.0%,102,137', 'TGGCATACCGCTCTGTAACA', 'TGGCGTACGGCTCAGTAACA', 'TGGCATACGGCTCAGTAACA', 'GGCATACAGGCTCAGTAACA', 'AGGCGGACGGCTCAGTAACA', 'TAGCATACGGCTCAGTAACA', 'CGGCGAACGGCTCAGTAACA', 'CGGCATACGGCTCAGTAACA', 'AGGCGGACTGCTCAGTAACA']}
expected_sample_r_primers = {'2043': ['AACAC,52,2043,9,90.0%,10.0%,103,138', 'AACACGTAGATAACCTACCC', 'AACACGTGGATAACTTACCC', 'AACACGTAGTCAACATGCCC', 'AACACGTAGTCAACATGCCC', 'AACACGCAGTCAACCTAACT', 'AACACGTAGTCAACCTGCCC', 'AACACGTGGGTAATCTGCCC', 'AACACGTGGGTAACCTGCCC', 'AACACGTGGGTAATCTACCC'], '2041': ['TAACA,51,2041,9,90.0%,10.0%,102,137', 'TAACACGTAGATAACCTACC', 'TAACACGTGGATAACTTACC', 'TAACACGTAGTCAACATGCC', 'TAACACGTAGTCAACATGCC', 'TAACACGCAGTCAACCTAAC', 'TAACACGTAGTCAACCTGCC', 'TAACACGTGGGTAATCTGCC', 'TAACACGTGGGTAACCTGCC', 'TAACACGTGGGTAATCTACC']}
# Bad because first field should be 'C' or 'M'
bad_sample_hits_data = ['X,TAACA,51,2041,9,90.0%,10.0%,102,137', 'M,2041,TGGCATACCGCTCTGTAACA,TAACACGTAGATAACCTACC,AY800210', 'M,2041,TGGCGTACGGCTCAGTAACA,TAACACGTGGATAACTTACC,EU883771', 'M,2041,TGGCATACGGCTCAGTAACA,TAACACGTAGTCAACATGCC,EF503699', 'M,2041,GGCATACAGGCTCAGTAACA,TAACACGTAGTCAACATGCC,EF503697', 'M,2041,AGGCGGACGGCTCAGTAACA,TAACACGCAGTCAACCTAAC,AJ535128', 'M,2041,TAGCATACGGCTCAGTAACA,TAACACGTAGTCAACCTGCC,EF021349', 'M,2041,CGGCGAACGGCTCAGTAACA,TAACACGTGGGTAATCTGCC,AJ969787', 'M,2041,CGGCATACGGCTCAGTAACA,TAACACGTGGGTAACCTGCC,AY360086', 'M,2041,AGGCGGACTGCTCAGTAACA,TAACACGTGGGTAATCTACC,DQ641863', 'C,AACAC,52,2043,9,90.0%,10.0%,103,138', 'M,2043,GGCATACCGCTCTGTAACAC,AACACGTAGATAACCTACCC,AY800210', 'M,2043,GGCGTACGGCTCAGTAACAC,AACACGTGGATAACTTACCC,EU883771', 'M,2043,GGCATACGGCTCAGTAACAC,AACACGTAGTCAACATGCCC,EF503699', 'M,2043,GCATACAGGCTCAGTAACAC,AACACGTAGTCAACATGCCC,EF503697', 'M,2043,GGCGGACGGCTCAGTAACAC,AACACGCAGTCAACCTAACT,AJ535128', 'M,2043,AGCATACGGCTCAGTAACAC,AACACGTAGTCAACCTGCCC,EF021349', 'M,2043,GGCGAACGGCTCAGTAACAC,AACACGTGGGTAATCTGCCC,AJ969787', 'M,2043,GGCATACGGCTCAGTAACAC,AACACGTGGGTAACCTGCCC,AY360086', 'M,2043,GGCGGACTGCTCAGTAACAC,AACACGTGGGTAATCTACCC,DQ641863']

known_primers_data = ['# Current Univeral Primers',
"# primer_id <tab> primer sequence (5'->3')  <tab> primer citation'",
"1f\tATTCCGGTTGATCCTGC\tBaker and Cowan, 2004",
"21f\tTTCCGGTTGATCCYGCCGGA\tEllis 2008, DeLong 1992 (also listed as F1 in Baker and Cowan, 2004)",
"27f\tAGAGTTTGATCMTGGCTCAG",
"# Partial overlap, no 3' match",
"test_known_primer_f\tGGCATACGGCTCAGT",
"# overlap and 3' match",
"test_known_primer_r\tGTTAACTACGTGTTA"]

known_primers_no_overlap_data = ['# Current Univeral Primers',
"# primer_id <tab> primer sequence (5'->3')  <tab> primer citation'",
"1f\tATTCCGGTTGATCCTGC\tBaker and Cowan, 2004",
"21f\tTTCCGGTTGATCCYGCCGGA\tEllis 2008, DeLong 1992 (also listed as F1 in Baker and Cowan, 2004)",
"27f\tAGAGTTTGATCMTGGCTCAG"]

bad_known_primers_data = ['# Current Univeral Primers',
"# primer_id <tab> primer sequence (5'->3')  <tab> primer citation'",
"1f\tATTCCGGTTGATCCTGC\tBaker and Cowan, 2004",
"21f\tTTCCGGTTGATCCYGCCGGA\tEllis 2008, DeLong 1992 (also listed as F1 in Baker and Cowan, 2004)",
"27f\tAGAGTTTGATCMTGGCTCAG",
"# Bad primer format, does not end with 'r' or 'f'",
"Bad_Primer_format\tGGCATACGGCTCAGT"]

# Start sample data for input files

sample_primer_hits_file = """# Matches that are found in at least 60.00% target sequences and specific at or below the 1.00% threshold 
# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches
C,TTGGG,253,14978,6,60.0%,0.0%
M,14978,CGTTATCCGGATTTATTGGG,TTGGGCGTAAAGAGCTCGTA,AB240712
M,14978,CGTTGTCCGGAATTATTGGG,TTGGGCGTAAAGAGCGCGTA,EU909933
M,14978,GTTATCCCGAATTTATTGGG,TTGGGGTGTAAAAGGTTGCG,AJ863512
M,14978,GCGTTATCGGAATTATTGGG,TTGGGCGTAAGAGTGCGTAG,FJ032552
M,14978,CGTTATCCGGAATTATTGGG,TTGGGCGTAAAGCGCGCGCA,EU693511
M,14978,TCTTATTCGGATTTATTGGG,TTGGGCGTAAAGCGTCTGCA,AF523993
C,AGCCG,392,13138,8,80.0%,0.0%
M,13138,CAACTACGTGCCAGCAGCCG,AGCCGCGGTGACACGTAGGC,AB240712
M,13138,ATCTCTGTGCCCAGCAGCCG,AGCCGCGGTTAATACAGAGG,EU434431
M,13138,TAACTACGTGCCAGCAGCCG,AGCCGCGGTAATACGTANGA,EU909933
M,13138,TAACTACGTGCCAGCAGCCG,AGCCGCGTATACGTAGGGGG,FJ032552
M,13138,TAACTACGTGCCAGCAGCCG,AGCCGCGGTAATACGTAGGT,EU693511
M,13138,TAAATTTGTGCCAGCAGCCG,AGCCGCGGTAATACAAATGA,AF523993
M,13138,TAACTCCGTGCCAGCAGCCG,AGCCGCGGTAATACGGAGGG,AY664363
M,13138,TAACTTCGTGCCAGCAGCCG,AGCCGCGGTAAGACGAAGGT,DQ310742
C,TTATT,250,14972,6,60.0%,0.0%
M,14972,AAGCGTTATCCGGATTTATT,TTATTGGGCGTAAAGAGCTC,AB240712
M,14972,AAGCGTTGTCCGGAATTATT,TTATTGGGCGTAAAGAGCGC,EU909933
M,14972,AGCGTTATCCCGAATTTATT,TTATTGGGGTGTAAAAGGTT,AJ863512
M,14972,CAAGCGTTATCGGAATTATT,TTATTGGGCGTAAGAGTGCG,FJ032552
M,14972,AAGCGTTATCCGGAATTATT,TTATTGGGCGTAAAGCGCGC,EU693511
M,14972,AAGTCTTATTCGGATTTATT,TTATTGGGCGTAAAGCGTCT,AF523993"""

sample_primer_hits_file_std_alignment = """# Matches that are found in at least 60.00% target sequences and specific at or below the 1.00% threshold 
# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches,5' standard index, 3' standard index
C,TTGGG,253,14978,6,60.0%,0.0%,548,583
M,14978,CGTTATCCGGATTTATTGGG,TTGGGCGTAAAGAGCTCGTA,AB240712
M,14978,CGTTGTCCGGAATTATTGGG,TTGGGCGTAAAGAGCGCGTA,EU909933
M,14978,GTTATCCCGAATTTATTGGG,TTGGGGTGTAAAAGGTTGCG,AJ863512
M,14978,GCGTTATCGGAATTATTGGG,TTGGGCGTAAGAGTGCGTAG,FJ032552
M,14978,CGTTATCCGGAATTATTGGG,TTGGGCGTAAAGCGCGCGCA,EU693511
M,14978,TCTTATTCGGATTTATTGGG,TTGGGCGTAAAGCGTCTGCA,AF523993
C,AGCCG,392,13138,8,80.0%,0.0%,507,542
M,13138,CAACTACGTGCCAGCAGCCG,AGCCGCGGTGACACGTAGGC,AB240712
M,13138,ATCTCTGTGCCCAGCAGCCG,AGCCGCGGTTAATACAGAGG,EU434431
M,13138,TAACTACGTGCCAGCAGCCG,AGCCGCGGTAATACGTANGA,EU909933
M,13138,TAACTACGTGCCAGCAGCCG,AGCCGCGTATACGTAGGGGG,FJ032552
M,13138,TAACTACGTGCCAGCAGCCG,AGCCGCGGTAATACGTAGGT,EU693511
M,13138,TAAATTTGTGCCAGCAGCCG,AGCCGCGGTAATACAAATGA,AF523993
M,13138,TAACTCCGTGCCAGCAGCCG,AGCCGCGGTAATACGGAGGG,AY664363
M,13138,TAACTTCGTGCCAGCAGCCG,AGCCGCGGTAAGACGAAGGT,DQ310742
C,TTATT,250,14972,6,60.0%,0.0%,745,780
M,14972,AAGCGTTATCCGGATTTATT,TTATTGGGCGTAAAGAGCTC,AB240712
M,14972,AAGCGTTGTCCGGAATTATT,TTATTGGGCGTAAAGAGCGC,EU909933
M,14972,AGCGTTATCCCGAATTTATT,TTATTGGGGTGTAAAAGGTT,AJ863512
M,14972,CAAGCGTTATCGGAATTATT,TTATTGGGCGTAAGAGTGCG,FJ032552
M,14972,AAGCGTTATCCGGAATTATT,TTATTGGGCGTAAAGCGCGC,EU693511
M,14972,AAGTCTTATTCGGATTTATT,TTATTGGGCGTAAAGCGTCT,AF523993"""

sample_known_primers_file = """# Current Univeral Primers
# primer_id <tab> primer sequence (5'->3')  <tab> primer citation
1f	ATTCCGGTTGATCCTGC	 Baker and Cowan, 2004
21f	TTCCGGTTGATCCYGCCGGA	 Ellis 2008, DeLong 1992 (also listed as F1 in Baker and Cowan, 2004)
27f	AGAGTTTGATCMTGGCTCAG	
Samplef	GGATTTATTGGG"""


# Start expected output
expected_formatted_primers_default = """# primer_id <tab> primer sequence (5'->3')
392f	TAACTWCGTGCCAGCAGCCG
253f	SSTTWTYCGGAWTTATTGGG
250f	AAGSSTTWTYCGGAWTTATT
392r	HCCYHYGTRTWACCGCGGCT
253r	YRCRVSCKYTTTACGCCCAA
250r	RVSCKYTTTACGCCCAATAA""".split('\n')

expected_formatted_primers_std_index = """# primer_id <tab> primer sequence (5'->3')
507f\tTAACTWCGTGCCAGCAGCCG
548f\tSSTTWTYCGGAWTTATTGGG
745f\tAAGSSTTWTYCGGAWTTATT
542r\tHCCYHYGTRTWACCGCGGCT
583r\tYRCRVSCKYTTTACGCCCAA
780r\tRVSCKYTTTACGCCCAATAA""".split('\n')

expected_formatted_primers_45perc_var_pos_freq = """# primer_id <tab> primer sequence (5'->3')
392f	TAACTACGTGCCAGCAGCCG
253f	CGTTATCCGGAWTTATTGGG
250f	AAGCGTTATCCGGAWTTATT
392r	CCCTWCGTATTACCGCGGCT
253r	TGCGVGCTCTTTACGCCCAA
250r	GVGCTCTTTACGCCCAATAA""".split('\n')

expected_formatted_primers_degeneracy_sort = """# primer_id <tab> primer sequence (5'->3')
392f	TAACTWCGTGCCAGCAGCCG
253f	SSTTWTYCGGAWTTATTGGG
250f	AAGSSTTWTYCGGAWTTATT
250r	RVSCKYTTTACGCCCAATAA
392r	HCCYHYGTRTWACCGCGGCT
253r	YRCRVSCKYTTTACGCCCAA""".split('\n')

expected_formatted_primers_specificity_sort = """# primer_id <tab> primer sequence (5'->3')
253f	SSTTWTYCGGAWTTATTGGG
392f	TAACTWCGTGCCAGCAGCCG
250f	AAGSSTTWTYCGGAWTTATT
253r	YRCRVSCKYTTTACGCCCAA
392r	HCCYHYGTRTWACCGCGGCT
250r	RVSCKYTTTACGCCCAATAA""".split('\n')

expected_formatted_primers_sensitivity_sort = """# primer_id <tab> primer sequence (5'->3')
392f	TAACTWCGTGCCAGCAGCCG
253f	SSTTWTYCGGAWTTATTGGG
250f	AAGSSTTWTYCGGAWTTATT
392r	HCCYHYGTRTWACCGCGGCT
253r	YRCRVSCKYTTTACGCCCAA
250r	RVSCKYTTTACGCCCAATAA""".split('\n')

expected_formatted_primers_with_name = """# primer_id <tab> primer sequence (5'->3')
CaptHammer392f	TAACTWCGTGCCAGCAGCCG
CaptHammer253f	SSTTWTYCGGAWTTATTGGG
CaptHammer250f	AAGSSTTWTYCGGAWTTATT
CaptHammer392r	HCCYHYGTRTWACCGCGGCT
CaptHammer253r	YRCRVSCKYTTTACGCCCAA
CaptHammer250r	RVSCKYTTTACGCCCAATAA""".split('\n')

expected_primers_details_default = """# Primer Report
# This file contains details about the prospective primers given in the formatted_primers.txt file.
# Conserved_Xmer is the conserved site found with the sequence_searcher module that serves as the primer 3' region.
# Unaligned_Index is the first position the X-mer was found when running the sequence_searcher module after degapping.
# Aligned_Index is the first position the X-mer was found when running the sequence_searcher module.
# Number_Hits is the number of sequences that had a perfect match to the X-mer.
# Percent_Match is the percentage of total sequences with a perfect match.
# Nonspecific_Match is the percentage of sequences in exclusion files that had perfect matches.  0.0% if no files were specified.
# 5'_Standard_Index gives a corrected index based on an alignment given with the -a parameter with the sequence_searcher module.  Will be empty if no standard alignment file was specified.
# 3'_Standard_Index gives a corrected index based on an alignment given with the -a parameter with the sequence_searcher module.
# F_Primer_IUPAC_Sequence shows the fully degenerate version of the prospective forward primer.  If any gaps are present, this base will be given as a '.' character.
# F_Primer_Filtered_Sequence shows a filtered version of the degenerate sequence where bases have to be present more than the specified frequency (see -p parameter) to contribute to degeneracy.
# F_Primer_Majority_Consensus gives the most frequently occurring bases.
# R_Primer data are the same as the F_Primer, but for the reverse primer, and have been reverse complemented.
# Start primer data:
# Conserved_Xmer, Unaligned_Index, Aligned_Index, Number_Hits, Percent_Match, Nonspecific_Match, 5'_Standard_Index, 3'_Standard_Index, F_Primer_IUPAC_Sequence, F_Primer_Filtered_Sequence, F_Primer_Majority_Consensus, R_Primer_IUPAC_Sequence, R_Primer_Filtered_Sequence, R_Primer_Majority_Consensus
TTGGG,253,14978,6,60.0%,0.0%,,,BBKWDHYCGRAWTTATTGGG,SSTTWTYCGGAWTTATTGGG,CGTTATCCGGAATTATTGGG,YDMVVVMBYYTWMMSCCCAA,YRCRVSCKYTTTACGCCCAA,TGCGGGCTCTTTACGCCCAA
AGCCG,392,13138,8,80.0%,0.0%,,,HWMHYHBKKSCCAGCAGCCG,TAACTWCGTGCCAGCAGCCG,TAACTACGTGCCAGCAGCCG,NCNYHBDHNTHWMCGCGGCT,HCCYHYGTRTWACCGCGGCT,CCCTTCGTATTACCGCGGCT
TTATT,250,14972,6,60.0%,0.0%,,,MRVBBKWDHYCGRAWTTATT,AAGSSTTWTYCGGAWTTATT,AAGCGTTATCCGGAATTATT,VVVMBYYTWMMSCCCAATAA,RVSCKYTTTACGCCCAATAA,GGGCTCTTTACGCCCAATAA""".split('\n')

expected_primers_details_default_std_index = """# Primer Report
# This file contains details about the prospective primers given in the formatted_primers.txt file.
# Conserved_Xmer is the conserved site found with the sequence_searcher module that serves as the primer 3' region.
# Unaligned_Index is the first position the X-mer was found when running the sequence_searcher module after degapping.
# Aligned_Index is the first position the X-mer was found when running the sequence_searcher module.
# Number_Hits is the number of sequences that had a perfect match to the X-mer.
# Percent_Match is the percentage of total sequences with a perfect match.
# Nonspecific_Match is the percentage of sequences in exclusion files that had perfect matches.  0.0% if no files were specified.
# 5'_Standard_Index gives a corrected index based on an alignment given with the -a parameter with the sequence_searcher module.  Will be empty if no standard alignment file was specified.
# 3'_Standard_Index gives a corrected index based on an alignment given with the -a parameter with the sequence_searcher module.
# F_Primer_IUPAC_Sequence shows the fully degenerate version of the prospective forward primer.  If any gaps are present, this base will be given as a '.' character.
# F_Primer_Filtered_Sequence shows a filtered version of the degenerate sequence where bases have to be present more than the specified frequency (see -p parameter) to contribute to degeneracy.
# F_Primer_Majority_Consensus gives the most frequently occurring bases.
# R_Primer data are the same as the F_Primer, but for the reverse primer, and have been reverse complemented.
# Start primer data:
# Conserved_Xmer, Unaligned_Index, Aligned_Index, Number_Hits, Percent_Match, Nonspecific_Match, 5'_Standard_Index, 3'_Standard_Index, F_Primer_IUPAC_Sequence, F_Primer_Filtered_Sequence, F_Primer_Majority_Consensus, R_Primer_IUPAC_Sequence, R_Primer_Filtered_Sequence, R_Primer_Majority_Consensus
TTGGG,253,14978,6,60.0%,0.0%,548,583,BBKWDHYCGRAWTTATTGGG,SSTTWTYCGGAWTTATTGGG,CGTTATCCGGAATTATTGGG,YDMVVVMBYYTWMMSCCCAA,YRCRVSCKYTTTACGCCCAA,TGCGGGCTCTTTACGCCCAA
AGCCG,392,13138,8,80.0%,0.0%,507,542,HWMHYHBKKSCCAGCAGCCG,TAACTWCGTGCCAGCAGCCG,TAACTACGTGCCAGCAGCCG,NCNYHBDHNTHWMCGCGGCT,HCCYHYGTRTWACCGCGGCT,CCCTTCGTATTACCGCGGCT
TTATT,250,14972,6,60.0%,0.0%,745,780,MRVBBKWDHYCGRAWTTATT,AAGSSTTWTYCGGAWTTATT,AAGCGTTATCCGGAATTATT,VVVMBYYTWMMSCCCAATAA,RVSCKYTTTACGCCCAATAA,GGGCTCTTTACGCCCAATAA""".split('\n')


expected_primers_details_45perc_var_pos_freq = """# Primer Report
# This file contains details about the prospective primers given in the formatted_primers.txt file.
# Conserved_Xmer is the conserved site found with the sequence_searcher module that serves as the primer 3' region.
# Unaligned_Index is the first position the X-mer was found when running the sequence_searcher module after degapping.
# Aligned_Index is the first position the X-mer was found when running the sequence_searcher module.
# Number_Hits is the number of sequences that had a perfect match to the X-mer.
# Percent_Match is the percentage of total sequences with a perfect match.
# Nonspecific_Match is the percentage of sequences in exclusion files that had perfect matches.  0.0% if no files were specified.
# 5'_Standard_Index gives a corrected index based on an alignment given with the -a parameter with the sequence_searcher module.  Will be empty if no standard alignment file was specified.
# 3'_Standard_Index gives a corrected index based on an alignment given with the -a parameter with the sequence_searcher module.
# F_Primer_IUPAC_Sequence shows the fully degenerate version of the prospective forward primer.  If any gaps are present, this base will be given as a '.' character.
# F_Primer_Filtered_Sequence shows a filtered version of the degenerate sequence where bases have to be present more than the specified frequency (see -p parameter) to contribute to degeneracy.
# F_Primer_Majority_Consensus gives the most frequently occurring bases.
# R_Primer data are the same as the F_Primer, but for the reverse primer, and have been reverse complemented.
# Start primer data:
# Conserved_Xmer, Unaligned_Index, Aligned_Index, Number_Hits, Percent_Match, Nonspecific_Match, 5'_Standard_Index, 3'_Standard_Index, F_Primer_IUPAC_Sequence, F_Primer_Filtered_Sequence, F_Primer_Majority_Consensus, R_Primer_IUPAC_Sequence, R_Primer_Filtered_Sequence, R_Primer_Majority_Consensus
TTGGG,253,14978,6,60.0%,0.0%,,,BBKWDHYCGRAWTTATTGGG,CGTTATCCGGAWTTATTGGG,CGTTATCCGGAATTATTGGG,YDMVVVMBYYTWMMSCCCAA,TGCGVGCTCTTTACGCCCAA,TGCGGGCTCTTTACGCCCAA
AGCCG,392,13138,8,80.0%,0.0%,,,HWMHYHBKKSCCAGCAGCCG,TAACTACGTGCCAGCAGCCG,TAACTACGTGCCAGCAGCCG,NCNYHBDHNTHWMCGCGGCT,CCCTWCGTATTACCGCGGCT,CCCTTCGTATTACCGCGGCT
TTATT,250,14972,6,60.0%,0.0%,,,MRVBBKWDHYCGRAWTTATT,AAGCGTTATCCGGAWTTATT,AAGCGTTATCCGGAATTATT,VVVMBYYTWMMSCCCAATAA,GVGCTCTTTACGCCCAATAA,GGGCTCTTTACGCCCAATAA""".split('\n')

expected_primers_overlap = """#Unique primers that do not overlap or share 3' ends with known primers
#primer name<tab>primer sequence
392f	TAACTWCGTGCCAGCAGCCG
250f	AAGSSTTWTYCGGAWTTATT
253r	YRCRVSCKYTTTACGCCCAA
392r	HCCYHYGTRTWACCGCGGCT
250r	RVSCKYTTTACGCCCAATAA
#Primers overlapping with known primers.
#primer name<tab>primer sequence<tab>known primer name<tab>known primer sequence<tab>overlapping sequence
253f	SSTTWTYCGGAWTTATTGGG	Samplef	GGATTTATTGGG	GGATTTATTG
#Primers sharing 3' ends with known primers.
#primer name<tab>primer sequence<tab>known primer name<tab>known primer sequence
253f	SSTTWTYCGGAWTTATTGGG	Samplef	GGATTTATTGGG""".split('\n')

expected_primers_overlap_match_len_5 = """#Unique primers that do not overlap or share 3' ends with known primers
#primer name<tab>primer sequence
#Primers overlapping with known primers.
#primer name<tab>primer sequence<tab>known primer name<tab>known primer sequence<tab>overlapping sequence
253f	SSTTWTYCGGAWTTATTGGG	1f	ATTCCGGTTGATCCTGC	TTCCG
253f	SSTTWTYCGGAWTTATTGGG	21f	TTCCGGTTGATCCYGCCGGA	TTCCG
253f	SSTTWTYCGGAWTTATTGGG	Samplef	GGATTTATTGGG	TTATT
392f	TAACTWCGTGCCAGCAGCCG	1f	ATTCCGGTTGATCCTGC	TTGAT
392f	TAACTWCGTGCCAGCAGCCG	21f	TTCCGGTTGATCCYGCCGGA	TTGAT
392f	TAACTWCGTGCCAGCAGCCG	27f	AGAGTTTGATCMTGGCTCAG	TTGAT
250f	AAGSSTTWTYCGGAWTTATT	1f	ATTCCGGTTGATCCTGC	TTCCG
250f	AAGSSTTWTYCGGAWTTATT	21f	TTCCGGTTGATCCYGCCGGA	TTCCG
250f	AAGSSTTWTYCGGAWTTATT	Samplef	GGATTTATTGGG	TTATT
253r	YRCRVSCKYTTTACGCCCAA	1f	ATTCCGGTTGATCCTGC	TTCCG
253r	YRCRVSCKYTTTACGCCCAA	21f	TTCCGGTTGATCCYGCCGGA	CCCGC
253r	YRCRVSCKYTTTACGCCCAA	27f	AGAGTTTGATCMTGGCTCAG	GGCTC
253r	YRCRVSCKYTTTACGCCCAA	Samplef	GGATTTATTGGG	TTGGG
392r	HCCYHYGTRTWACCGCGGCT	1f	ATTCCGGTTGATCCTGC	ATCCT
392r	HCCYHYGTRTWACCGCGGCT	21f	TTCCGGTTGATCCYGCCGGA	CGCCG
392r	HCCYHYGTRTWACCGCGGCT	27f	AGAGTTTGATCMTGGCTCAG	ATGGC
392r	HCCYHYGTRTWACCGCGGCT	Samplef	GGATTTATTGGG	TATTG
250r	RVSCKYTTTACGCCCAATAA	1f	ATTCCGGTTGATCCTGC	TTCCG
250r	RVSCKYTTTACGCCCAATAA	21f	TTCCGGTTGATCCYGCCGGA	CCCGC
250r	RVSCKYTTTACGCCCAATAA	27f	AGAGTTTGATCMTGGCTCAG	GGCTC
250r	RVSCKYTTTACGCCCAATAA	Samplef	GGATTTATTGGG	TTATT
#Primers sharing 3' ends with known primers.
#primer name<tab>primer sequence<tab>known primer name<tab>known primer sequence
253f	SSTTWTYCGGAWTTATTGGG	Samplef	GGATTTATTGGG""".split('\n')

expected_primers_overlap_trunc_len_1 = """#Unique primers that do not overlap or share 3' ends with known primers
#primer name<tab>primer sequence
250f	AAGSSTTWTYCGGAWTTATT
253r	YRCRVSCKYTTTACGCCCAA
392r	HCCYHYGTRTWACCGCGGCT
250r	RVSCKYTTTACGCCCAATAA
#Primers overlapping with known primers.
#primer name<tab>primer sequence<tab>known primer name<tab>known primer sequence<tab>overlapping sequence
253f	SSTTWTYCGGAWTTATTGGG	Samplef	GGATTTATTGGG	GGATTTATTG
#Primers sharing 3' ends with known primers.
#primer name<tab>primer sequence<tab>known primer name<tab>known primer sequence
253f	SSTTWTYCGGAWTTATTGGG	27f	AGAGTTTGATCMTGGCTCAG
253f	SSTTWTYCGGAWTTATTGGG	Samplef	GGATTTATTGGG
392f	TAACTWCGTGCCAGCAGCCG	27f	AGAGTTTGATCMTGGCTCAG
392f	TAACTWCGTGCCAGCAGCCG	Samplef	GGATTTATTGGG""".split('\n')

expected_amplicon_pairs = """# Primer pairs within estimated amplicon size
# Min size: 200
# Max size: 600

507f	TAACTWCGTGCCAGCAGCCG
780r	RVSCKYTTTACGCCCAATAA
Estimated amplicon size: 233""".split('\n')

if __name__ == "__main__":
    main()
