#!/usr/bin/env python
#created 11-10-09

__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

from os.path import exists
from shutil import rmtree

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, get_random_directory_name
from cogent.app.util import get_tmp_filename

from primerprospector.get_amplicons_and_reads import generate_paired_amplicons,\
 generate_unidirectional_amplicons, get_output_name_single_primer,\
 get_output_name_primer_pair, \
 generate_amplicons_from_hits_files, get_output_name_reads, generate_reads,\
 get_hits_files, get_amplicons_and_reads


class GetAmpliconsAndReadsTests(TestCase):
    """unit tests for get_amplicons_and_reads """
    
    def setUp(self):
        """ Contains tmp files setup for unit tests"""
        
        self.small_forward_hits_file = small_forward_hits_file
        self.small_reverse_hits_file = small_reverse_hits_file
        self.expected_amplicons_score_three = expected_amplicons_score_three
        self.fasta_seqs_real = fasta_seqs_real
        self.sample_forward_primer_hits_real = sample_forward_primer_hits_real
        self.sample_reverse_primer_hits_real = sample_reverse_primer_hits_real
        self.forward_primer_hits_not_matching_labels_real =\
         forward_primer_hits_not_matching_labels_real
         
        self.small_forward_hits = small_forward_hits
        self.small_reverse_hits = small_reverse_hits
        self.small_fasta_data = small_fasta_data
        self.expected_amplicons_default_setting =\
         expected_amplicons_default_setting
        self.expected_amplicons_score_three = expected_amplicons_score_three
        self.expected_uni_f_amplicons_default_setting =\
         expected_uni_f_amplicons_default_setting
        self.expected_uni_f_amplicons_lenient_setting =\
         expected_uni_f_amplicons_lenient_setting
        self.expected_uni_r_amplicons_default_setting =\
         expected_uni_r_amplicons_default_setting
        self.expected_reads_f_4 = expected_reads_f_4
        self.expected_reads_r_4 = expected_reads_r_4
        self.expected_amplicons_default_setting_real_seqs =\
         expected_amplicons_default_setting_real_seqs
        self.expected_reads_default_setting_real_seqs =\
         expected_reads_default_setting_real_seqs
        self.expected_amplicons_score_two_real_seqs =\
         expected_amplicons_score_two_real_seqs
        self.expected_reads_score_two_real_seqs =\
         expected_reads_score_two_real_seqs
        self.expected_amplicons_reverse_primer_only_real_seqs =\
         expected_amplicons_reverse_primer_only_real_seqs
        self.expected_reads_reverse_primer_only_real_seqs =\
         expected_reads_reverse_primer_only_real_seqs
        self.expected_amplicons_forward_primer_only_real_seqs =\
         expected_amplicons_forward_primer_only_real_seqs
        self.expected_reads_forward_primer_only_real_seqs =\
         expected_reads_forward_primer_only_real_seqs
        
        self.input_forward_hits1_fp = get_tmp_filename(\
         prefix = '100f_', suffix = '.txt')
        hits_file = open(self.input_forward_hits1_fp, 'w')
        hits_file.write(self.small_forward_hits_file)
        hits_file.close()
        
        self.input_reverse_hits1_fp = get_tmp_filename(\
         prefix = '300r_', suffix = '.txt')
        hits_file = open(self.input_reverse_hits1_fp, 'w')
        hits_file.write(self.small_reverse_hits_file)
        hits_file.close()
        
        self.amplicons_fp = get_tmp_filename(\
         prefix = '100f_300r_amplicons', suffix = '.fasta')
        amplicon_file = open(self.amplicons_fp, 'w')
        amplicon_file.write('\n'.join(self.expected_amplicons_score_three))
        amplicon_file.close()
        
        self.real_fasta_fp = get_tmp_filename(\
         prefix = 'real_fasta_seqs_', suffix = '.fasta')
        real_fasta_f = open(self.real_fasta_fp, 'w')
        real_fasta_f.write(self.fasta_seqs_real)
        real_fasta_f.close()
        
        self.real_forward_hits_fp = get_tmp_filename(\
         prefix = '515f_hits', suffix = '.txt')
        real_forward_hits_f = open(self.real_forward_hits_fp, 'w')
        real_forward_hits_f.write(self.sample_forward_primer_hits_real)
        real_forward_hits_f.close()
        
        self.real_reverse_hits_fp = get_tmp_filename(\
         prefix = '806r_hits', suffix = '.txt')
        real_reverse_hits_f = open(self.real_reverse_hits_fp, 'w')
        real_reverse_hits_f.write(self.sample_reverse_primer_hits_real)
        real_reverse_hits_f.close()
        
        self.real_forward_hits_non_matching_labels_fp = get_tmp_filename(\
         prefix = '515f_hits', suffix = '.txt')
        real_forward_hits_f = \
         open(self.real_forward_hits_non_matching_labels_fp, 'w')
        real_forward_hits_f.write(self.forward_primer_hits_not_matching_labels_real)
        real_forward_hits_f.close()
        
        self.output_dir = get_random_directory_name(prefix = "/tmp/amplicons_",
         suffix = "_and_reads")
        self.output_dir += '/'
        
        self._files_to_remove =\
         [self.input_forward_hits1_fp, self.input_reverse_hits1_fp,
          self.amplicons_fp, self.real_fasta_fp, self.real_forward_hits_fp,
          self.real_reverse_hits_fp, 
          self.real_forward_hits_non_matching_labels_fp]
         
    def tearDown(self):
        remove_files(self._files_to_remove)
        if exists(self.output_dir):
            rmtree(self.output_dir)
        
    def test_generate_paired_amplicons_weighted_score_default(self):
        """ Creates list of amplicons from primer pair correctly """
        
        primer_hits_data = [self.small_forward_hits, self.small_reverse_hits]
        score_threshold = 1.0
        min_seq_len = 10
        fasta_data = self.small_fasta_data
        score_type = 'weighted_score'
        expected_amplicons = self.expected_amplicons_default_setting
        

        actual_amplicons = generate_paired_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_paired_amplicons_weighted_score_amplicons(self):
        """ Does not include amplicons that are too short """
        
        primer_hits_data = [self.small_forward_hits, self.small_reverse_hits]
        score_threshold = 1.0
        min_seq_len = 12
        fasta_data = self.small_fasta_data
        score_type = 'weighted_score'
        
        # Result should be empty since all amplicons are too short
        expected_amplicons = []
        

        actual_amplicons = generate_paired_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_paired_amplicons_weighted_score_lenient(self):
        """ Creates list of amplicons from primer pair correctly """
        
        primer_hits_data = [self.small_forward_hits, self.small_reverse_hits]
        score_threshold = 3.0
        min_seq_len = 10
        fasta_data = self.small_fasta_data
        score_type = 'weighted_score'
        expected_amplicons = self.expected_amplicons_score_three
        

        actual_amplicons = generate_paired_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_paired_amplicons_tp_mm_default(self):
        """ Creates list of amplicons from primer pair correctly """
        
        primer_hits_data = [self.small_forward_hits, self.small_reverse_hits]
        score_threshold = 1.0
        min_seq_len = 10
        fasta_data = self.small_fasta_data
        score_type = 'tp_mismatches'
        # One hit has 2 3' mismatches, so only should get one amplicon result
        expected_amplicons = self.expected_amplicons_default_setting
        

        actual_amplicons = generate_paired_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_paired_amplicons_tp_mm_lenient(self):
        """ Creates list of amplicons from primer pair correctly """
        
        primer_hits_data = [self.small_forward_hits, self.small_reverse_hits]
        score_threshold = 3.0
        min_seq_len = 10
        fasta_data = self.small_fasta_data
        score_type = 'tp_mismatches'
        # One hit has 2 3' mismatches, with score of 3.0 should get both
        expected_amplicons = self.expected_amplicons_score_three
        

        actual_amplicons = generate_paired_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_paired_amplicons_overall_mm_default(self):
        """ Creates list of amplicons from primer pair correctly """
        
        primer_hits_data = [self.small_forward_hits, self.small_reverse_hits]
        score_threshold = 1.0
        min_seq_len = 10
        fasta_data = self.small_fasta_data
        score_type = 'overall_mismatches'
        # One hit has 3 total mismatches, so only should get one amplicon result
        expected_amplicons = self.expected_amplicons_default_setting
        

        actual_amplicons = generate_paired_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_paired_amplicons_overall_mm_lenient(self):
        """ Creates list of amplicons from primer pair correctly """
        
        primer_hits_data = [self.small_forward_hits, self.small_reverse_hits]
        score_threshold = 3.0
        min_seq_len = 10
        fasta_data = self.small_fasta_data
        score_type = 'overall_mismatches'
        # Should get 2 amplicons with more lenient score
        expected_amplicons = self.expected_amplicons_score_three
        

        actual_amplicons = generate_paired_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_unidirectional_amplicons_forward_default(self):
        """ Creates list of amplicons from single primer correctly """
        
        primer_hits_data = self.small_forward_hits
        score_threshold = 1.0
        min_seq_len = 10
        fasta_data = self.small_fasta_data
        score_type = 'weighted_score'
        expected_amplicons = self.expected_uni_f_amplicons_default_setting
        primer_direction = "f"
        
        # Only should return a single amplicon
        actual_amplicons = generate_unidirectional_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type, primer_direction)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_unidirectional_amplicons_forward_lenient(self):
        """ Creates list of amplicons from single primer correctly """
        
        primer_hits_data = self.small_forward_hits
        score_threshold = 3.0
        min_seq_len = 10
        fasta_data = self.small_fasta_data
        score_type = 'weighted_score'
        expected_amplicons = self.expected_uni_f_amplicons_lenient_setting
        primer_direction = "f"
        
        # Should get 2 amplicons with lenient setting
        actual_amplicons = generate_unidirectional_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type, primer_direction)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_unidirectional_amplicons_forward_long_amplicons(self):
        """ Creates list of amplicons from single primer correctly """
        
        primer_hits_data = self.small_forward_hits
        score_threshold = 1.0
        min_seq_len = 40
        fasta_data = self.small_fasta_data
        score_type = 'weighted_score'
        expected_amplicons = []
        primer_direction = "f"    
        
        # Should not return anything because amplicons are too short
        actual_amplicons = generate_unidirectional_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type, primer_direction)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_unidirectional_amplicons_forward_tp_mismatches(self):
        """ Creates list of amplicons from single primer correctly """
        
        primer_hits_data = self.small_forward_hits
        score_threshold = 1.0
        min_seq_len = 10
        fasta_data = self.small_fasta_data
        score_type = 'tp_mismatches'
        expected_amplicons = self.expected_uni_f_amplicons_default_setting
        primer_direction = "f"
        
        # Only should return a single amplicon
        actual_amplicons = generate_unidirectional_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type, primer_direction)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_unidirectional_amplicons_forward_overall_mismatches(self):
        """ Creates list of amplicons from single primer correctly """
        
        primer_hits_data = self.small_forward_hits
        score_threshold = 1.0
        min_seq_len = 10
        fasta_data = self.small_fasta_data
        score_type = 'overall_mismatches'
        expected_amplicons = self.expected_uni_f_amplicons_default_setting
        primer_direction = "f"
        
        # Only should return a single amplicon
        actual_amplicons = generate_unidirectional_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type, primer_direction)
         
        self.assertEqual(actual_amplicons, expected_amplicons)
        
    def test_generate_unidirectional_amplicons_reverse_default(self):
        """ Creates list of amplicons from single primer correctly """
        
        primer_hits_data = self.small_reverse_hits
        score_threshold = 1.0
        min_seq_len = 5
        fasta_data = self.small_fasta_data
        score_type = 'weighted_score'
        expected_amplicons = self.expected_uni_r_amplicons_default_setting
        primer_direction = "r"
        
        # Only should return a single amplicon
        actual_amplicons = generate_unidirectional_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type, primer_direction)
         
        self.assertEqual(actual_amplicons, expected_amplicons)

    def test_get_output_name_single_primer(self):
        """ Generates unidirectional amplicon filename correctly """
        
        # Test with current dir as output directory, and hits file in current
        # directory.  No file IO, just string manipulation
        output_dir = "."
        hits_filepath = "623f_hits.txt"
        
        expected_outfile_name = "./623f_amplicons.fasta"
        
        actual_outfile_name = get_output_name_single_primer(hits_filepath,
         output_dir)
         
        self.assertEqual(actual_outfile_name, expected_outfile_name)
        
        # Test for relative paths result
        
        output_dir = "/home/user/Desktop/primers/"
        hits_filepath = "../436r_hits.txt"
        
        expected_outfile_name = \
         "/home/user/Desktop/primers/436r_amplicons.fasta"
         
        actual_outfile_name = get_output_name_single_primer(hits_filepath,
         output_dir)
         
        self.assertEqual(actual_outfile_name, expected_outfile_name)
        
        
    def test_get_output_name_primer_pair(self):
        """ Generates primer pair amplicons filename correctly """
        
        # Test for case of all paths in current directory
        # No file IO, just string manipulation
        output_dir = "."
        hits_filepath_f = "623f_hits.txt"
        hits_filepath_r = "1029r_hits.txt"
        hits_filepaths = (hits_filepath_f, hits_filepath_r)
        
        expected_outfile_name = "./623f_1029r_amplicons.fasta"
        
        actual_outfile_name = get_output_name_primer_pair(hits_filepaths,
         output_dir)
         
        self.assertEqual(actual_outfile_name, expected_outfile_name)
        
        # Test for relative paths result
        
        output_dir = "/home/user/Desktop/primers/"
        hits_filepath_f = "../436f_hits.txt"
        hits_filepath_r = "/test/1244r_hits.txt"
        hits_filepaths = (hits_filepath_f, hits_filepath_r)
        
        expected_outfile_name = \
         "/home/user/Desktop/primers/436f_1244r_amplicons.fasta"
         
        actual_outfile_name = get_output_name_primer_pair(hits_filepaths,
         output_dir)
         
        self.assertEqual(actual_outfile_name, expected_outfile_name)
        


        
    def test_generate_amplicons_from_hits_files_pair(self):
        """ Tests pair or unidirectional primers, writes amplicons """
        
        hits_files = [self.input_forward_hits1_fp, self.input_reverse_hits1_fp]
        fasta_data = self.small_fasta_data
        
        score_type = "weighted_score"
        score_threshold = 1.0
        min_seq_len = 5
        
        # Expected file to be created with these input data
        expected_outfile_path = self.output_dir + '100f_300r_amplicons.fasta'
        self._files_to_remove.append(expected_outfile_path)
        
        generate_amplicons_from_hits_files(hits_files, fasta_data, 
         self.output_dir, score_type, score_threshold, min_seq_len)
         
        # Read actual results
        actual_results_f = open(expected_outfile_path, "U")
        actual_results = [line.strip() for line in actual_results_f]
        
        self.assertEqual(actual_results, 
         self.expected_amplicons_default_setting)
        
    def test_generate_amplicons_from_hits_files_pair_lenient(self):
        """ Tests pair or unidirectional primers, writes amplicons """
        
        hits_files = [self.input_forward_hits1_fp, self.input_reverse_hits1_fp]
        fasta_data = self.small_fasta_data

        score_type = "weighted_score"
        # Set more lenient score, should get 2 amplicon result
        score_threshold = 3.0
        min_seq_len = 5
        
        # Expected file to be created with these input data
        expected_outfile_path = self.output_dir + '100f_300r_amplicons.fasta'
        self._files_to_remove.append(expected_outfile_path)
        
        generate_amplicons_from_hits_files(hits_files, fasta_data, 
         self.output_dir, score_type, score_threshold, min_seq_len)
         
        # Read actual results
        actual_results_f = open(expected_outfile_path, "U")
        actual_results = [line.strip() for line in actual_results_f]
        
        self.assertEqual(actual_results, self.expected_amplicons_score_three)
        
    def test_generate_amplicons_from_hits_files_pair_total_mm(self):
        """ Tests pair or unidirectional primers, writes amplicons """
        
        hits_files = [self.input_forward_hits1_fp, self.input_reverse_hits1_fp]
        fasta_data = self.small_fasta_data

        # With overall_mm and score of 1.0, should have a single amplicon
        score_type = "overall_mismatches"
        score_threshold = 1.0
        min_seq_len = 5
        
        # Expected file to be created with these input data
        expected_outfile_path = self.output_dir + '100f_300r_amplicons.fasta'
        self._files_to_remove.append(expected_outfile_path)
        
        generate_amplicons_from_hits_files(hits_files, fasta_data, 
         self.output_dir, score_type, score_threshold, min_seq_len)
         
        # Read actual results
        actual_results_f = open(expected_outfile_path, "U")
        actual_results = [line.strip() for line in actual_results_f]
        
        self.assertEqual(actual_results,
         self.expected_amplicons_default_setting)
        
    def test_generate_amplicons_from_hits_files_pair_tp_mm(self):
        """ Tests pair or unidirectional primers, writes amplicons """
        
        hits_files = [self.input_forward_hits1_fp, self.input_reverse_hits1_fp]
        fasta_data = self.small_fasta_data

        # With tp_mismatches and score of 3.0, should have both amplicons
        score_type = "tp_mismatches"
        score_threshold = 3.0
        min_seq_len = 5
        
        # Expected file to be created with these input data
        expected_outfile_path = self.output_dir + '100f_300r_amplicons.fasta'
        self._files_to_remove.append(expected_outfile_path)
        
        generate_amplicons_from_hits_files(hits_files, fasta_data, 
         self.output_dir, score_type, score_threshold, min_seq_len)
         
        # Read actual results
        actual_results_f = open(expected_outfile_path, "U")
        actual_results = [line.strip() for line in actual_results_f]
        
        self.assertEqual(actual_results, self.expected_amplicons_score_three)
        
    def test_generate_amplicons_from_hits_files_unidirection(self):
        """ Tests pair or unidirectional primers, writes amplicons """
        
        hits_files = [self.input_forward_hits1_fp]
        fasta_data = self.small_fasta_data

        score_type = "weighted_score"
        score_threshold = 1.0
        min_seq_len = 5
        
        # Expected file to be created with these input data
        expected_outfile_path = self.output_dir + '100f_amplicons.fasta'
        self._files_to_remove.append(expected_outfile_path)
        
        generate_amplicons_from_hits_files(hits_files, fasta_data, 
         self.output_dir, score_type, score_threshold, min_seq_len)
         
        # Read actual results
        actual_results_f = open(expected_outfile_path, "U")
        actual_results = [line.strip() for line in actual_results_f]
        
        self.assertEqual(actual_results, 
         self.expected_uni_f_amplicons_default_setting)
         
    def test_generate_amplicons_from_hits_files_unidirection_lenient(self):
        """ Tests pair or unidirectional primers, writes amplicons """
        
        hits_files = [self.input_forward_hits1_fp]
        fasta_data = self.small_fasta_data

        score_type = "weighted_score"
        # Should generate 2 amplicons with score set at 3.0
        score_threshold = 3.0
        min_seq_len = 5
        
        # Expected file to be created with these input data
        expected_outfile_path = self.output_dir + '100f_amplicons.fasta'
        self._files_to_remove.append(expected_outfile_path)
        
        generate_amplicons_from_hits_files(hits_files, fasta_data, 
         self.output_dir, score_type, score_threshold, min_seq_len)
         
        # Read actual results
        actual_results_f = open(expected_outfile_path, "U")
        actual_results = [line.strip() for line in actual_results_f]
        
        self.assertEqual(actual_results, 
         self.expected_uni_f_amplicons_lenient_setting)
         
    def test_generate_amplicons_from_hits_files_unidirection_r(self):
        """ Tests pair or unidirectional primers, writes amplicons """
        
        hits_files = [self.input_reverse_hits1_fp]
        fasta_data = self.small_fasta_data

        score_type = "weighted_score"
        score_threshold = 1.0
        min_seq_len = 5
        
        # Expected file to be created with these input data
        expected_outfile_path = self.output_dir + '300r_amplicons.fasta'
        self._files_to_remove.append(expected_outfile_path)
        
        generate_amplicons_from_hits_files(hits_files, fasta_data,
         self.output_dir, score_type, score_threshold, min_seq_len)
         
        # Read actual results
        actual_results_f = open(expected_outfile_path, "U")
        actual_results = [line.strip() for line in actual_results_f]
        
        self.assertEqual(actual_results, 
         self.expected_uni_r_amplicons_default_setting)
        
        
    def test_get_output_name_reads(self):
        """ Generates reads output filenames correctly """
        
                          
        # Test case where all files are in current directory
        amplicons_fp = "450f_899r_amplicons.fasta"
        read_direction = "f"
        read_len = 250
        output_dir = "."
        
        expected_reads_output_name = ["./450f_899r_f_250_reads.fasta"]
        
        actual_reads_output_name = get_output_name_reads(amplicons_fp,
         read_direction, read_len, output_dir)
         
        self.assertEqual(actual_reads_output_name, expected_reads_output_name)
        
        # Test reverse reads
        amplicons_fp = "450f_899r_amplicons.fasta"
        read_direction = "r"
        read_len = 250
        output_dir = "."
        
        expected_reads_output_name = ["./450f_899r_r_250_reads.fasta"]
        
        actual_reads_output_name = get_output_name_reads(amplicons_fp,
         read_direction, read_len, output_dir)
         
        self.assertEqual(actual_reads_output_name, expected_reads_output_name)
        
        # Test paired end reads
        amplicons_fp = "450f_899r_amplicons.fasta"
        read_direction = "p"
        read_len = 250
        output_dir = "."
        
        expected_reads_output_name = ["./450f_899r_f_250_reads.fasta",
         "./450f_899r_r_250_reads.fasta"]
        
        actual_reads_output_name = get_output_name_reads(amplicons_fp,
         read_direction, read_len, output_dir)
         
        self.assertEqual(actual_reads_output_name, expected_reads_output_name)
        
    def test_generate_reads_forward(self):
        """ Creates, writes reads file(s) correctly """
                   
        # Test forward reads
        amplicons_fp = self.amplicons_fp
        read_direction = 'f'
        read_len = 4

        
        reads_fp = self.output_dir + '100f_300r_f_4_reads.fasta'
        self._files_to_remove.append(reads_fp)
        
        generate_reads(amplicons_fp, read_direction, read_len, self.output_dir)
        
        expected_reads = self.expected_reads_f_4
        
        actual_reads_f = open(reads_fp, "U")
        actual_reads = [line.strip() for line in actual_reads_f]
        
        self.assertEqual(actual_reads, expected_reads)
        
    def test_generate_reads_reverse(self):
        """ Creates, writes reads file(s) correctly """
                   
        # Test reverse reads
        amplicons_fp = self.amplicons_fp
        read_direction = 'r'
        read_len = 4

        
        reads_fp = self.output_dir + '100f_300r_r_4_reads.fasta'
        self._files_to_remove.append(reads_fp)
        
        generate_reads(amplicons_fp, read_direction, read_len, self.output_dir)
        
        expected_reads = self.expected_reads_r_4
        
        actual_reads_f = open(reads_fp, "U")
        actual_reads = [line.strip() for line in actual_reads_f]
        
        self.assertEqual(actual_reads, expected_reads)
        
    def test_generate_reads_paired_end(self):
        """ Creates, writes reads file(s) correctly """
                   
        # Test paired end reads
        amplicons_fp = self.amplicons_fp
        read_direction = 'p'
        read_len = 4

        
        reads_f_fp = self.output_dir + '100f_300r_f_4_reads.fasta'
        self._files_to_remove.append(reads_f_fp)
        reads_r_fp = self.output_dir + '100f_300r_r_4_reads.fasta'
        self._files_to_remove.append(reads_r_fp)
        
        generate_reads(amplicons_fp, read_direction, read_len, self.output_dir)
        
        expected_reads_f = self.expected_reads_f_4
        
        actual_reads_f = open(reads_f_fp, "U")
        actual_reads = [line.strip() for line in actual_reads_f]
        
        self.assertEqual(actual_reads, expected_reads_f)
        
        expected_reads_r = self.expected_reads_r_4
        
        actual_reads_r = open(reads_r_fp, "U")
        actual_reads = [line.strip() for line in actual_reads_r]
        
        self.assertEqual(actual_reads, expected_reads_r)
        
    def test_generate_reads_hits_seq_end(self):
        
        # Should get complete amplicon if read_len longer than actual seq size
        amplicons_fp = self.amplicons_fp
        read_direction = 'f'
        read_len = 30

        
        reads_fp = self.output_dir + '100f_300r_f_30_reads.fasta'
        self._files_to_remove.append(reads_fp)
        
        generate_reads(amplicons_fp, read_direction, read_len, self.output_dir)
        
        # should match input amplicons
        expected_reads = self.expected_amplicons_score_three
        
        actual_reads_f = open(reads_fp, "U")
        actual_reads = [line.strip() for line in actual_reads_f]
        
        self.assertEqual(actual_reads, expected_reads)
        
        
        
    def test_get_hits_files(self):
        """ Generates list of hits files, test for validity """
        
        
        # Pass a single hits filepath
        expected_hits_fp = [self.input_forward_hits1_fp]
        
        actual_hits_fp = get_hits_files(self.input_forward_hits1_fp)
        
        self.assertEqual(actual_hits_fp, expected_hits_fp)
        
        # pass pair of filepaths
        
        expected_hits_fp = [self.input_forward_hits1_fp, 
         self.input_reverse_hits1_fp]
         
        actual_hits_fp = get_hits_files(self.input_forward_hits1_fp +\
         ":" + self.input_reverse_hits1_fp)
        
        self.assertEqual(actual_hits_fp, expected_hits_fp)
        
        # Raise error if more than 2 files passed
        
        three_hits_fp = self.input_forward_hits1_fp + ":" + \
         self.input_forward_hits1_fp + ":" + self.input_forward_hits1_fp
         
        self.assertRaises(ValueError, get_hits_files, three_hits_fp)
        
        # Raise error if cannot open one of the filepaths
        
        bad_filepath = self.input_forward_hits1_fp + ":" + "fake_file_path.txt"
        
        self.assertRaises(IOError, get_hits_files, bad_filepath)
        
    def test_get_amplicons_and_reads_default_settings(self):
        """ Test for overall program functionality """
        
                            
        # Using real data for overall tests
        
        # Test default settings
        primer_hits = self.real_forward_hits_fp + ":" +\
         self.real_reverse_hits_fp
        fasta_fps = self.real_fasta_fp

        score_type = 'weighted_score'
        score_threshold = 1.0
        min_seq_len = 100
        read_direction = 'r'
        read_len = 250
        
        output_amplicons_fp = self.output_dir + '515f_806r_amplicons.fasta'
        output_reads_fp = self.output_dir + '515f_806r_r_250_reads.fasta'
        self._files_to_remove.append(output_amplicons_fp)
        self._files_to_remove.append(output_reads_fp)
        
        get_amplicons_and_reads(primer_hits, fasta_fps, self.output_dir, 
         score_type, score_threshold, min_seq_len, read_direction, read_len)
         
        expected_amplicons = self.expected_amplicons_default_setting_real_seqs
        
        actual_amplicons_f = open(output_amplicons_fp, "U")
        actual_amplicons = [line.strip() for line in actual_amplicons_f]
        
        self.assertEqual(actual_amplicons, expected_amplicons)
        
        expected_reads = self.expected_reads_default_setting_real_seqs
        
        actual_reads_f = open(output_reads_fp, "U")
        actual_reads = [line.strip() for line in actual_reads_f]
        
        self.assertEqual(actual_reads, expected_reads)
        
    def test_get_amplicons_and_reads_lenient_settings(self):
        """ Test for overall program functionality """
        
                            
        # Using real data for overall tests
        
        # score 2.0 setting, should get 3 amplicons/reads
        primer_hits = self.real_forward_hits_fp + ":" +\
         self.real_reverse_hits_fp
        fasta_fps = self.real_fasta_fp

        score_type = 'weighted_score'
        score_threshold = 2.0
        min_seq_len = 100
        read_direction = 'r'
        read_len = 250
        
        output_amplicons_fp = self.output_dir + '515f_806r_amplicons.fasta'
        output_reads_fp = self.output_dir + '515f_806r_r_250_reads.fasta'
        self._files_to_remove.append(output_amplicons_fp)
        self._files_to_remove.append(output_reads_fp)
        
        get_amplicons_and_reads(primer_hits, fasta_fps, self.output_dir,
         score_type, score_threshold, min_seq_len, read_direction, read_len)
         
        expected_amplicons = self.expected_amplicons_score_two_real_seqs
        
        actual_amplicons_f = open(output_amplicons_fp, "U")
        actual_amplicons = [line.strip() for line in actual_amplicons_f]
        
        self.assertEqual(actual_amplicons, expected_amplicons)
        
        expected_reads = self.expected_reads_score_two_real_seqs
        
        actual_reads_f = open(output_reads_fp, "U")
        actual_reads = [line.strip() for line in actual_reads_f]
        
        self.assertEqual(actual_reads, expected_reads)
        
    def test_get_amplicons_and_reads_single_primer(self):
        """ Test for overall program functionality """
        
                            
        # Using real data for overall tests
        
        # Use forward primer only and forward reads, should give 4 amplicons
        # and reads with default score threshold of 1.0
        primer_hits = self.real_forward_hits_fp
        fasta_fps = self.real_fasta_fp

        score_type = 'weighted_score'
        score_threshold = 1.0
        min_seq_len = 100
        read_direction = 'f'
        read_len = 250
        
        output_amplicons_fp = self.output_dir + '515f_amplicons.fasta'
        output_reads_fp = self.output_dir + '515f_f_250_reads.fasta'
        self._files_to_remove.append(output_amplicons_fp)
        self._files_to_remove.append(output_reads_fp)
        
        get_amplicons_and_reads(primer_hits, fasta_fps, self.output_dir,
         score_type, score_threshold, min_seq_len, read_direction, read_len)
         
        expected_amplicons =\
         self.expected_amplicons_forward_primer_only_real_seqs
        
        actual_amplicons_f = open(output_amplicons_fp, "U")
        actual_amplicons = [line.strip() for line in actual_amplicons_f]
        
        self.assertEqual(actual_amplicons, expected_amplicons)
        
        expected_reads =\
         self.expected_reads_forward_primer_only_real_seqs
        
        actual_reads_f = open(output_reads_fp, "U")
        actual_reads = [line.strip() for line in actual_reads_f]
        
        self.assertEqual(actual_reads, expected_reads)
        
    def test_get_amplicons_and_reads_reverse_primer(self):
        """ Test for overall program functionality """
        
                            
        # Using real data for overall tests
        
        # Reverse primer should only generate 2 amplicons/reads
        primer_hits = self.real_reverse_hits_fp
        fasta_fps = self.real_fasta_fp

        score_type = 'weighted_score'
        score_threshold = 1.0
        min_seq_len = 100
        read_direction = 'r'
        read_len = 250
        
        output_amplicons_fp = self.output_dir + '806r_amplicons.fasta'
        output_reads_fp = self.output_dir + '806r_r_250_reads.fasta'
        self._files_to_remove.append(output_amplicons_fp)
        self._files_to_remove.append(output_reads_fp)
        
        get_amplicons_and_reads(primer_hits, fasta_fps, self.output_dir,
         score_type, score_threshold, min_seq_len, read_direction, read_len)
         
        expected_amplicons =\
         self.expected_amplicons_reverse_primer_only_real_seqs
        
        actual_amplicons_f = open(output_amplicons_fp, "U")
        actual_amplicons = [line.strip() for line in actual_amplicons_f]
        
        self.assertEqual(actual_amplicons, expected_amplicons)
        
        expected_reads = self.expected_reads_reverse_primer_only_real_seqs
        
        actual_reads_f = open(output_reads_fp, "U")
        actual_reads = [line.strip() for line in actual_reads_f]
        
        self.assertEqual(actual_reads, expected_reads)
        
    def test_get_amplicons_and_reads_error_non_matching_hits(self):
        """ Test for overall program functionality """
        
                            
        # Using real data for overall tests
        
        # Should raise an error if hits files do not have matching sequence IDs
        primer_hits = self.real_forward_hits_non_matching_labels_fp + ":" +\
         self.real_reverse_hits_fp
        fasta_fps = self.real_fasta_fp

        score_type = 'weighted_score'
        score_threshold = 1.0
        min_seq_len = 100
        read_direction = 'r'
        read_len = 250
        
        
        self.assertRaises(ValueError, get_amplicons_and_reads, primer_hits, 
         fasta_fps, self.output_dir, score_type, score_threshold, min_seq_len, 
         read_direction, read_len)
         

        

# Large strings at the end for better readability
small_forward_hits = """ABCDEFG,AACCTAGACT,ACGGGTACMG,10,1,2,False,0,0,2.4,False
1234567,AACCTAGACT,ACGGGTACMG,5,0,0,False,0,0,0.0,False""".split('\n')

small_forward_hits_file = """# Primer: 233f 5'-ACGGGTACMG-3'
# Input fasta file: test_seqs.fasta
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
ABCDEFG,AACCTAGACT,ACGGGTACMG,10,1,2,False,0,0,2.4,False
1234567,AACCTAGACT,ACGGGTACMG,5,0,0,False,0,0,0.0,False"""

small_reverse_hits = """ABCDEFG,AACCTAGACT,ACGGGTACMG,30,0,1,False,0,0,1.0,False
1234567,AACCTAGACT,ACGGGTACMG,25,0,0,False,0,0,0.0,False""".split('\n')

small_reverse_hits_file = """# Primer: 899r 5'-ACGGGTACMG-3'
# Input fasta file: test_seqs.fasta
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
ABCDEFG,AACCTAGACT,ACGGGTACMG,30,0,1,False,0,0,1.0,False
1234567,AACCTAGACT,ACGGGTACMG,25,0,0,False,0,0,0.0,False"""

small_fasta_seqs = """>ABCDEFG
TACCAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCG
>1234567
TACCAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCG"""

small_fasta_data = {"ABCDEFG":"TACCAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCG",
"1234567":"TACCAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCG"}

expected_amplicons_default_setting = """>1234567
TGGTTGGGAC""".split('\n')

expected_amplicons_score_three = """>ABCDEFG
GGGACATTTA
>1234567
TGGTTGGGAC""".split('\n')

expected_uni_f_amplicons_default_setting = """>1234567
TGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCG""".split('\n')

expected_uni_f_amplicons_lenient_setting = """>ABCDEFG
GGGACATTTATTGGGCCTAAAGCATTCGTAGCCG
>1234567
TGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCG""".split('\n')

expected_uni_r_amplicons_default_setting = """>ABCDEFG
TACCAGCTCCCCGAGTGGTTGGGACATTTA
>1234567
TACCAGCTCCCCGAGTGGTTGGGAC""".split('\n')

expected_reads_f_4 = """>ABCDEFG
GGGA
>1234567
TGGT""".split('\n')

expected_reads_r_4 = """>ABCDEFG
TTTA
>1234567
GGAC""".split('\n')


sample_reverse_primer_hits_real = """# Primer: 806r 5'-GGACTACVSGGGTATCTAAT-3'
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
AB329763,ACTAAAAACCCCGGTAGTCC,ATTAGATACCCSBGTAGTCC,823,1,2,False,0,0,2.4,False
AJ937875,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,748,0,0,False,0,0,0.0,False
DQ243733,ATTAGATACCCGGGTAGGTA,ATTAGATACCCSBGTAGTCC,686,3,0,False,0,0,1.2,False
EF533958,ATTAGATACC-GGGTAGTCC,ATTAGATACCCSBGTAGTCC,724,0,0,False,1,0,1.0,False
DQ278152,ATTAGATACCCGGGTAGTCT,ATTAGATACCCSBGTAGTCC,548,0,0,True,0,0,3.0,False"""

sample_forward_primer_hits_real = """# Primer: 515f 5'-GTGCCAGCMGCCGCGGTAA-3'
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
AB329763,GTGCCAGCCGCCACGGGTA,GTGCCAGCMGCCGCGGTAA,483,1,2,False,0,0,2.4,False
AJ937875,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,476,1,0,False,0,0,0.4,False
DQ243733,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,414,1,0,False,0,0,0.4,False
EF533958,GTGCCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,451,0,0,False,0,0,0.0,False
DQ278152,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,276,1,0,False,0,0,0.4,False"""

forward_primer_hits_not_matching_labels_real = """# Primer: 515f 5'-GTGCCAGCMGCCGCGGTAA-3'
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
XXXXXXXX,GTGCCAGCCGCCACGGGTA,GTGCCAGCMGCCGCGGTAA,483,1,2,False,0,0,2.4,False
AJ937875,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,476,1,0,False,0,0,0.4,False
DQ243733,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,414,1,0,False,0,0,0.4,False
EF533958,GTGCCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,451,0,0,False,0,0,0.0,False
DQ278152,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,276,1,0,False,0,0,0.4,False"""

# Sequences taken from the Silva 97 dataset
fasta_seqs_real = """>AB329763 1 1516 Archaea/Euryarchaeota/Halobacteriales/uncultured
AGACCATTGCTATCCGAGTTCTACTAAGCCATGCGAGCTGAGAGGTGAAAGACCTCGGCGAATTGCTCCGTAATACGTAGTCAACTTGCCCTCAAGTGGGAGATAATCTCGGGAAACTGAGGCTAATATCCCATAGTTTTACATTACTGGAATGTGTGTAAAGCGAAAGCCACGGCGCTTGAGGATAGGACTGCGCCTGATTAGGCTGTTGTTGGTGTAAAGATACCACAATATGGCTTTATAGAACCAGCAAACCGATATGGAGGTAAGGGTATTCTTTCCTCCACGGCAAAATCAGTACGGGCTATGGAAGTAGGAGCCCGGAGATGGATTCTGAGACATGAATCCAGGTACTACGGTACGCAGCAGGCGAGAAACCTTTGCAATCGGTTAACGCCGGACAGGGGAACTTGGAGTGTTCTAGTTTTACTAGAACTTTTGCTCACTGAAAACGGGTGAGTGAATAAGGGCTGGCCAAGACGGGTGCCAGCCGCCACGGGTAATACCCGCAGCCCAAGTGATGGTCACGATTATTGGGTCTAAAGCGTCCGTAGCCGGTCTAACAAGTTCCTGGTGAAATCTTACAGCAAACTGTAAGGCTTGCTGGGGATACTGTTAGACTTGAGACCGGGAGAGGACAGAGGTACTCGTAGGGTAGGGGTTAAATCCTATAATCCTACGGGGACCACCTGTGGCGAAAGCGTCTGTCTAGAACGGGTCTGACGGTGAGGGACGAAAGCTAGGGGAGCGATCGGGATTAGATATGGAGAAGGAAAAAGTAATTCAATAAACGAAAGTTTAGAGAATCACCTATTCCAGACCCACTAAAAACCCCGGTAGTCCTAGCTGTAAACATTGCCCGCTTGGTGTTGCAGACCTCTTGGGGGTTTGCAGTGCCGGAGCGTAGGTGTTAAGCGGGCCACCTGGGGAGTACAGTCGCAAGGCCGAAACTTAAGGGAATTGGCGGGGGAGCACACAAGGGGTGGGCAGTGCGGTTCAATTAGATTCAACGCTGGAAATCTTACCAGGGGCGACAGCAGGATGAGGGTCAGTCTGAAGGGCTTACCAGACAAGCTGAGAGGTGGTGCATGGCCATCGTCAGCTCGTACCGTGAGGCGTGCCGTTAAGTCGGTTAACGAGCGAGACCCACATTTCATGTTGCAACTATACTTTCCGAAGTATAGGCACTCATGAGAGACCGCTGGTGATAAACCAGAGGAAGGTGTGGGCGACGGTAGGTCTGTATGCCCCGGATCTCCTGGGCTACACGCGCTGCACAATGCGTGCCACAATGGGAAGCAACTCCGAGAGGAGAAGCGAATCCCCTAAAAGCACGCTTAGTTCAGATTGAGGGTTGCAACTCACCCTCATGAAGCCGGAATCCCTAGTAGGCGAATGTCACTAAGTAGTAGACTGCTAATGAATTACATGCACCACGCTTGATTAAAGATGATAGACAGAGGGGAAAGAATCAGTCTGCTATCGTGAAATCGTTCGCCGAATACGTCCCTGCTCC
>EF533958 1 1551 Archaea/Euryarchaeota/Halobacteriales/Halobacteriales/Halogeometricum
TCCGGTTGATCCTGCCGGAGGCCATTGCTATCGGAGTCCGACTTAGCCATGCTAGTTGCGCGAGTTTAGACTCGCAGCAGATAGCTCAGTAACACGTGGTCAAGCTACCCTGCAGACACGGACAACCTCGGGAAACTGAGGCTAATCCGCGATATCGATCCCACGCTGGAACAGCCGGGATCAGCAAACGCTCCGGCGCTGCAGGATGCGGCTGCGGCCGATTAGGTAGACGGTGGGGTAACGGCCCACCGTGCCGATAATCGGTACGGGTTGTGAGAGCAAGAGCCCGGAGACGGGATCTGAGACAAGATTCCGGGCCCTACGGGGCGCAGCAGGCGCGAAACCTTTACACTGCACGCAAGTGCGATAAGGGGACTCCGAGTGCGGGGGCATATAGTCCTCGCTTTTGAGAACCGTAAGGCGGTTCTCGAATAAGAGCTGGGCAAGACCGGTGCCAGCCGCCGCGGTAATACCGGCAGCTCAAGTGATGGCCAATCTTATTGGGCCTAAAGCGTCCGTAGCTGGCCGTGAAAGTTCGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGAAAACTTCACGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCTGGGGTAGGAGTGAAATCCCGTAATCCTGGACGGACCGCCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACAGTGAGGGACGAAAGCTAGGGTCTCGAACCGGATTAGATACCCGGGTAGTCCTAGCCGTAAACAATGTTCACTAGGTGTGGCACAGGCTACGAGCCTGTGCTGTGCCGTAGGGAAGCCGAGAAGTGAACCGCCTGGGAAGTACGTCCGCAAGGATGAAACTTAAAGGAATTGGCGGGGGAGCACTACAACCGGAGGAGCCTGCGGTTTAATTGGACTCAACGCCGGACATCTCACCAGCATCGACTGTAGTAATGACGACCAGGTTGATGACCTTGTCCGAGCTTCAGAGAGGAGGTGCATGGCCGCCGTCAGCTCGTACCGTGAGGCGTCCTGTTAAGTCAGGCAACGAGCGAGACCCGCATCCTTACTTGCCAGCAGCACTGCGAAGTGGCTGGGGACAGTAGGGAGACCGCCGTGGCCAACACGGAGGAAGGAACGGGCAACGGTAGGTCAGTATGCCCCGAATGAGCTGGGCTACACGCGGGCTACAATGGCCAGAACAATGGGTTGCAACCCCGAAAGGGAACGCTAATCTCCGAAATCTGGTCGTAGTTCGGATTGAGGGCTGAAACTCGCCCTCATGAAGCTGGATTCGGTAGTAATCGCATTTCAGAAGAGTGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCAAAGCACCCGAGTGAGGTCCGGATGAGGCCATCGCAAGATGGTCGAATCTGGGCTTCGCAAGGGGGCTTAAGTCGTAACAAGGTAGCCGTAGGGGAATCTGCGGCTGGATCACCTCCTGATGGTCGAATCTGGGCTTCGCAAGGGGGCTTAAGTCGTAACAAGGTAGCCGTAGGGGAATCTGCGGCTGGATCACCTCCT
>DQ278152 1 1571 Archaea/Crenarchaeota/uncultured/uncultured
TGCGCGAGCTTCTTTTGTGCATCCGCTCTAGGATGGGACTGCGGCCGATCAGGCTGTTGGTGGGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCCTTGAGAGAGGGGGCCCGGAGATGGACACTGAGAGAAGGGTCCAGGTCCTACGGGGCGCAGCAGGCGCGAAACCTTCGCAATGCGCGAAAGCGTGACGAGGTTAATCCGAGTGGTCCCCGCTGAGGGGATCTTTTCTTGGCTCTAAAACGGCCAGGGAATAAGGGGAGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCTCCTCGAGTGGTCGGGGCGATTATTGGGCCTAAAGCATCCGTAGCCGGTTCCGTAGGTCTTCTGTTAAATCCAACGGCTTAACCGTTGGCCTGCAGAAGATACCGCGGGACTAGGAGGCGGGAGAGGTGGACGGTACTCCACGTGTAGGGGTAAAATCCTTTGATCCGTGGAAGACCACCAGTGGCGAAGGCGGTCCACCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAACCGGATTAGATACCCGGGTAGTCCCAGCCGTAAACGATGCGAGCTAGATGATCCAATCGCAAATCGCGATTGGAGTGTCGCAGGGAAGCCGTTAAGCTCGCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAAGGGGTGAAGCCTGCGGTTCAATTGGAGTCAACGCCAGGAAACTTACCAGGAGAGACAGCAGGTTGAGGGTCAAGCTGAAGACTTTACCTGACGAGCTGAGAGGAGGTGCATGGCCGTCGCCAGCTCGTGCCGTGAGGCGTCCTGTTAAGTCAGGCAACGAGCGAGACCCCTATTGGTAATTGCTAACAGTCCAGAGATGGATTGAGGCTAGTTACCGAGACTGCCGCTGCTAAAGCGGAGGAAGGAGGGGGCCACGGCAGGTCAGTATGCCCCGAAACTCCTGGGCCACACGCGGGCTGCAATGGTAGGGACAATTGGTTCCGACTTCGAAAGAAGGAGGCAATCCGCAAACCCTACCCCAGTTGTGATTGAGGATTGAAACCCGTCCTCATGAATATGGAATCCCTAGTAACCGCGTGTCATCACCGCGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCGCTTCACCCGAGTTTGGCTAGGACGAGGTGACCTCTAATTGGGGTTATCGAATCTATGTTAGGCAAGGGGGGAGAAGCCGTAACAAGGTGGCTGTAGGGGAACCTGCAGCCGGATCACCTCCTTAGTACATACGGACGAACAACAGACGTACGTACTGAAATCGGGGGTTGGAATACGCACTGCAATGGAGTGAAAAGCAAATGCCAGTGTCAGACGTTGACAGACAGATGAAGATGCACTGTCTCACGTGAGTGAGATAGGAAGGGCTAGCAGGGAACGGAAGGGCCTTCGCGGGCGCCCCTGCGATGTTGCCATGCATAGGGTTCGCGCAAGCATCCGCCGGGGATGTTAATCATCCGACACAACGACGCTGTATGGTGAATGGCTTGGCTTGGT
>AJ937875 1 1523 Archaea/Crenarchaeota/uncultured/uncultured
TTCCGGTTGATCCTGCCGGACCCGACTGCTATCGGAGTGGGGCTAAGCCATGCGAGTCTTACGTCTCTTGACTGCAAGAGGCGTGGCGACCGGCTCAGTAACACGTAGCCAACCTAACCTAGAGACGTGGACAACCCCGGGAAACTGGGATTAATCCACGATAGATCAAGGCTCCTGGAATGGGCTTTGGTCCAAAAGGATGGATAAGCATGTTTTGTCCATCTGCTCTAGGATGGGGCTGCGGCCGATCAGGTTGTTGATGGGGTAATGGCCCACCAAGCCTATAACCGGTACGGGCTTTGAGAGAAGGAGCCCGGAGATGGACACTGAGATAAGGGTCCAGGCCCTACGGGGCGCAGCAGGCGCGAAACCTTCGCAATGCGCGAAAGCGTGACGAGGTTACTTCGAGTGCCACTCGCTGAGGGTGGCTTTTCCCAAATCTAAAAAGCTTGGGGAATAAGGGGAGGGTAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCGGTCTTGTAGGTCTCCTGTTAAACTCAACGGCTCAACCGTTGAACTGCAGGAGATACCGCAGGACTAGGAGGCGGGAGAGGCGGACGGTACTCCATGGGTAGCGGCAAAATGTTCTGATCCATGGAAGACCACCAGTGGCGTAGGCGGTCCGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAACCGGATTAGATACCCGGGTAGTCCCAGCCGTAAACGATGCGAGCTAGGTGATGGAATGGCTGCGTGCCATTCCAGTGCCGCAGGGAAGCCGTTAAGCTCGCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAAGGGGTGAAGCCTGCGGTTTAATTGGAGTCAACGCCGGGAACCTTACCCAGGGAGACAGCAGGATGAAGGTCAAGCTAAAGACTTTACCAGACAAGCTGAGAGGTGGTGCATGGCCGTCGCCAGCTCGTGCCGTGAGGCGTCCTGTTAAGTCAGGCAACGAGCGAGACCCTTGCCGCTAGTTGCTATCTCCATCCGAAAGGATCGAGAGGCTAATTAGCGGGACCGCTGCCGCTAAGGCGGAGGAAGGTAGGGGCCACGGCAGGTCAGTATGCCCCGAAACCCTGGGGCCACACGCGGGCTGCAATGGTAGGGACAATGGGTTTCGACCCCGAGAGGGGGAGACAATCCCTAAACCCTACCCCAGTTGTGACTGAGGGCTGCAACCCGCCCTCATGAATCTGGAATCCCTAGTAACCGCGTGTCACCATCGCGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCGCTTCACCCGAGTTGGGTTTGGGTGAGGTGGTATCATATTGGTACCATCGAACTCGAATTCGGCAAGGGGGGAGAAGTCGTAACAAGGTAACCAATCACTAGTGCGGCCGCCTGCAGGTCGACCATATGGGAGAGTCCCACGCGTGA
>DQ243733 1 1537 Archaea/Crenarchaeota/uncultured
CGTCCCTCTGCCGAGGGCGTGGCGCACGGCTGATAATACGCGGCTAACCTGCCCTCGGGAGGGGGGTAACCCCGGGAAACTGGGACTAATCCCCCCATAGACGGGGAGGCCTGGAAGGGTTCCCCGTCGAAAAGGGCGCTTAGGGGTCATCGCTACGCGCTCGCCCGAGGATGGGGCCGCGTCCCATCAGGTAGATGGCGGGGTAAAGGCCCGCCATGCCTATAACGGGTAGGGGCCGTGAGAGCGGGAGCCCCCAGATGGGCACTGAGACAAGGGCCCAGGCCCTACGGGGCGCACCAGGCGGGAAACCTCCGCAATGCGGGAAACCGTGACGGGGTCACCCCGAGTGCCCTCCTTTGGGAGGGCTTTTCCCCGCTGTAAACAGGCGGGGGCAATAAGCGGGGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCCCCGCGAGTGGTCGGGACGGTTATTGGGCCTAAAGCGCCCGTAGCCGGCCTAGCAAGTCTCCGCGGAAATCCCCGGGCTCAACCCGGGGGCATGCGGAGATACTGCTAGGCTAGGGGGCGGGAGAGGCCGGGGGTACTCCCGGGGTAGGGGCGAAATCCTATAATCCCGGGAGGACCACCAGTGGCGAAGGCGCCCGGCTGGAACGCGCCCGACGGTGAGGGGCGAAAGCCGGGGGAGCGAACCGGATTAGATACCCGGGTAGGTATGAGCCGCGGGATAAAGTGTTCTCTATCTAAGTGCAAAGATGGAGCCTACTGTACAGATATCAAAGGTTCGCATACACTTACAGATGTAAACCCGAACAGATGCGGGGAACCTATGAAAACCCTAAAGGGTGAGTGTTCCCTTTAGGGCGGGGAGAAGGTCATCCGCGGCCATACCCCGGAGAGTCCCGGCTGTAAACGATGCGGGCTAGGTGTCGGACGGGCTTAGAGCCCGTCCGGTGCCGCAGGGAAGCCGTTAAGCCCGCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAAGGGGTGGAGCCTGCGGCTTAATTGGAGTCAACGCCAGGAAACTTACCGGGGGCGACAGCAGGATGACGGCCAGGTTAAAGGCCTTGCCCGACGCGCTGAGAGGAGGTGCACGGCCGTCGCCAGCTCGTGCCGTGAGGTGTCCGGTTAAGTCCGGCAACGAGCGAGACCCCCGCCCCTAGTTGCTACCCGGTCCTTCGGGGCCGGGGCACACTAGGGGGACTGCCGGCGTTAAGCCGGAGGAAGGTGGGGGCCACGGCAGGTCAGCATAGGGGCCTATCCCGAGGGAGGGCCCCGCTGAAGCCCCGAAACCCCCGGGCTGCACGCGGGCTACAATGGCCGGGACAGCGGGATCCGACCCCGAAAGGGGGAGGCAATCCCTTAAACCCGGCCTCAGTTGGGATCGAGGGCTGTAACTCGCCCTCGTGAACGCGGAATCCCTAGTAACCGCGCGTCATCATCGCGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCAAGGGC
>AB293215 1 1537 Archaea/Euryarchaeota/Thermoplasmatales
TTCTGGTTGATCCTGCCGGAGGCCACTGCTATCCGGCTCCGACTAAGCCATGCGAGTCTGGGGCCGGGGAAACCCGGCACCGGCGGACGGCTCAGTAACACGTGGCCAACCTACCCTCGGGTGGGGGATAACCCCGGGAAAGCCGGTAGGGGCGCGGGAGAGCACTTGCATTTTACCCCCTTGGGCCGAGCAAAGGATCATGGACCCTCTGGCCCAGGGGCACATAGGGCAGCCATCAGGGCTGCTTAGGCTTCGGGGCGTTGCCCCATAGCGGCAGTCACTCTCAATGGCGCACAGGTCCATAGGCCCCTCCGGCCCCGAAACTGGGGCTAATCCCCCATAGGGGAGGGATTGCTGGAAGGCGCCCTCCCCGAAAGCGCGAGCGCCCGAGGATGGGGCTGCGGCCTATCAGGTAGTTGTGGGGGTAACGGCCCCACAAGCCTGTGACGGGTACGGGCCGTGGGAGCGGGAGCCCGGAGATGGCCCCTGAGACACGGGGCCAGGCCCTACGGGGCGCAGCAGGCGCGAAACCTGCGCAATGCGCGCAAGCGCGACGCGGGGAGCCGGAGTGCCCGTGCTTTGCACGGGCTTTTCCCGTGTCTAAAAAGCACGGGGAATAAGGGCTGGGCAAGGCCGGTGGCAGCCGCCGCGGTAATACCGGCAGCCCGAGTGGTGGCCGCTTCTATTGAGCCTAAAGCGTCCGTAGCCGGCCCCGTAAATCCCTGGGTAAATCGGCCCGCTCAACGGGCCGCTCCCCGGGGAGACTGCGGGGCTTGGGACCGGGAGAGGCCGAGGGTACTCCCGGGGTAGGGGTAAAATCCTATAATCCCGGGGGGACCACCTGTGGGGAAGCCGCTCGGCTGGAACGGCTCCGACGGTGAGGGACGAAGGCCAGGGGATCAAACCGGATTAGATACCCGGGTAGTCCTGGCTGTAAACGCTGCCCGCTTGGTGTTGCCCGCCTCTCGGGGGCGGGCAGTGCCGTAGGGAAGCCGTTAAGCGGGCCGCTTGGGGAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGGGCACTGCAACGGGAGGATGCGTGCGGTTTAATTGGATTCAACGCCGGAAACCTCACCGGGGGAGACCGGCACATGAGGGTCAGGCTGAAGACCTTACCCGATAGCCGGAGAGGTGGTGCATGGCGGTCGTCAGTTCGTACCGTGGGGCGTCCTCTTAAGTGAGATAACGAACGAGACCCCCGCCCCCAGTTGCCACCCGGATCTCCCGGTCCGGGGGGAAACTGGGGGGACTGCCGGCGCTAAGCCGGAGGAAGGTGGGGGCAACGGTACGTCCGTATGCCCCGAATCCCCCGGGCTACACGCGCGCTACAAAGGGCGGGACAATGGGTGTGGCGACCCCGAAAGGGGGAGCTGATCCCGAAACCCGCCCGTAGTTCGGATCGAGGGCTGTAACTCGCCCTCGTGAAGGTGGATTCCGTAGTAATCGCGGGTCAACATCCCGCGGTGAATACGCCCCTGCCCCTTGCACACACTGCCCGTC
>DQ278134 1 1647 Archaea/Crenarchaeota/uncultured/uncultured
AGTCAACATGCCCAGGGGACGTGGATAACCTCGGGAAACTGAGAATAAACCATGATAGGTCATGATTTCTGGAATGGATCATGGCTCAAATCTACATGGCCCCTGGATTGGACTGCGGCCGATCAGGCAGTTGGTGAGGTAAAGGCCCACCAAACCTATAACCGGTACGGGCTCTGAGAGGAGAAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCGCGAAACCTCTGCAATAGGCGAAAGCCTGACAGGGTTACTCTGAGTGATTTCCGCTAAGGGGATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTGGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCGGTTCTGCAAGTCTTCCGTTAAATCCAGCCGCTCAACGGATGGGCCGCGGAGGATACTACAGGACTAGGAGGCGGGAGAGGCAAGCGGTACTCAGTGGGTAGGGGTAAAATCCTTTGATCCATTGAAGACCACCAGTGGCGAAGGCGGCTTGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCAAACCGGATTAGATACCCGGGTAGTCCCAGCTGTAAACGATGCAGACTCAGTGATGGGCTAGCCTTGTGCCAGCCCAGTGCCGCAGGGAAGCCGTTAAGTCTGCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAAGGGGTGAAGCCTGCGGTTCAATTGGAGTCAACGCCGGGAATCTTACCGGGGGCGACAGCAGAGTGAAGGTCAAGCCGAAGACTTTACCAGACAAGCTGAGAGGAGGTGCATGGCCGTCGCCAGCTCGTGCCGTGAGGTGTCCTGTTAAGTCAGGTAACGAGCGAGATCCCTACCTCTAGTTGCTACTGCTATTCTCAGGAGTAGCAGAGCTAATTAGAGGGACCGCCGCCGCTGAGGCGGAGGAAGGAGGGGGCTACGGCAGGTCAGTATGCCCCGAAACCCTCGGGCCACACGCGGGCTGCAATGGTAAGGACAATGGGTTCCTATTCCGAAAGGAGGAGGCAATCCCTAAACCTTACCCCAGTTATGATTGAGGGCTGAAACTCGCCCTCATGAATATGGAATCCCTAGTAACCGCGTGACATCATCGCGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCGCTTCATCGAAGTTGGTTCTTGGCGAGGTTGTGTCTAATTGGCACTATCGAACCTGGGGTCAGCAACGAGGGAGAAGTCGTAACAAGGTGGCCGTAGGGGAACCTGCGGCCGGATCACCTCCTTCGATTAAGGCCGGCAACTTAAGTGTAAAGGGTCATGTGCCATTGTACGTACGAATGCAATTCTCTGCAGGATATGCAGGGAATGAAGGACCGAGCAAACGCAAGTTTGCAATGACCGTCATGTATAGGATAAATCTGCATCTCTCAAACGGTTGTTAGACATAATAGTCTGACAGCAATGATAGATGAAATGTGGCTTTACAATGTTCACGTAAAGTTTCAATGCTTTAGACTACGCTGGGAAAACTGGTGCAAAGACGCCGGTTGGTGGATG"""

expected_amplicons_default_setting_real_seqs = """>AJ937875
TACCAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCGGTCTTGTAGGTCTCCTGTTAAACTCAACGGCTCAACCGTTGAACTGCAGGAGATACCGCAGGACTAGGAGGCGGGAGAGGCGGACGGTACTCCATGGGTAGCGGCAAAATGTTCTGATCCATGGAAGACCACCAGTGGCGTAGGCGGTCCGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAACCGG
>EF533958
TACCGGCAGCTCAAGTGATGGCCAATCTTATTGGGCCTAAAGCGTCCGTAGCTGGCCGTGAAAGTTCGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGAAAACTTCACGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCTGGGGTAGGAGTGAAATCCCGTAATCCTGGACGGACCGCCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACAGTGAGGGACGAAAGCTAGGGTCTCGAACCGG""".split('\n')

expected_reads_default_setting_real_seqs = """>AJ937875
CAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCGGTCTTGTAGGTCTCCTGTTAAACTCAACGGCTCAACCGTTGAACTGCAGGAGATACCGCAGGACTAGGAGGCGGGAGAGGCGGACGGTACTCCATGGGTAGCGGCAAAATGTTCTGATCCATGGAAGACCACCAGTGGCGTAGGCGGTCCGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAACCGG
>EF533958
GGCAGCTCAAGTGATGGCCAATCTTATTGGGCCTAAAGCGTCCGTAGCTGGCCGTGAAAGTTCGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGAAAACTTCACGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCTGGGGTAGGAGTGAAATCCCGTAATCCTGGACGGACCGCCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACAGTGAGGGACGAAAGCTAGGGTCTCGAACCGG""".split('\n')


expected_amplicons_score_two_real_seqs = """>AJ937875
TACCAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCGGTCTTGTAGGTCTCCTGTTAAACTCAACGGCTCAACCGTTGAACTGCAGGAGATACCGCAGGACTAGGAGGCGGGAGAGGCGGACGGTACTCCATGGGTAGCGGCAAAATGTTCTGATCCATGGAAGACCACCAGTGGCGTAGGCGGTCCGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAACCGG
>DQ243733
TACCAGCCCCGCGAGTGGTCGGGACGGTTATTGGGCCTAAAGCGCCCGTAGCCGGCCTAGCAAGTCTCCGCGGAAATCCCCGGGCTCAACCCGGGGGCATGCGGAGATACTGCTAGGCTAGGGGGCGGGAGAGGCCGGGGGTACTCCCGGGGTAGGGGCGAAATCCTATAATCCCGGGAGGACCACCAGTGGCGAAGGCGCCCGGCTGGAACGCGCCCGACGGTGAGGGGCGAAAGCCGGGGGAGCGAACCGG
>EF533958
TACCGGCAGCTCAAGTGATGGCCAATCTTATTGGGCCTAAAGCGTCCGTAGCTGGCCGTGAAAGTTCGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGAAAACTTCACGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCTGGGGTAGGAGTGAAATCCCGTAATCCTGGACGGACCGCCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACAGTGAGGGACGAAAGCTAGGGTCTCGAACCGG""".split('\n')


expected_reads_score_two_real_seqs = """>AJ937875
CAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCGGTCTTGTAGGTCTCCTGTTAAACTCAACGGCTCAACCGTTGAACTGCAGGAGATACCGCAGGACTAGGAGGCGGGAGAGGCGGACGGTACTCCATGGGTAGCGGCAAAATGTTCTGATCCATGGAAGACCACCAGTGGCGTAGGCGGTCCGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAACCGG
>DQ243733
CAGCCCCGCGAGTGGTCGGGACGGTTATTGGGCCTAAAGCGCCCGTAGCCGGCCTAGCAAGTCTCCGCGGAAATCCCCGGGCTCAACCCGGGGGCATGCGGAGATACTGCTAGGCTAGGGGGCGGGAGAGGCCGGGGGTACTCCCGGGGTAGGGGCGAAATCCTATAATCCCGGGAGGACCACCAGTGGCGAAGGCGCCCGGCTGGAACGCGCCCGACGGTGAGGGGCGAAAGCCGGGGGAGCGAACCGG
>EF533958
GGCAGCTCAAGTGATGGCCAATCTTATTGGGCCTAAAGCGTCCGTAGCTGGCCGTGAAAGTTCGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGAAAACTTCACGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCTGGGGTAGGAGTGAAATCCCGTAATCCTGGACGGACCGCCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACAGTGAGGGACGAAAGCTAGGGTCTCGAACCGG""".split('\n')

expected_amplicons_forward_primer_only_real_seqs = """>AJ937875
TACCAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCGGTCTTGTAGGTCTCCTGTTAAACTCAACGGCTCAACCGTTGAACTGCAGGAGATACCGCAGGACTAGGAGGCGGGAGAGGCGGACGGTACTCCATGGGTAGCGGCAAAATGTTCTGATCCATGGAAGACCACCAGTGGCGTAGGCGGTCCGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAACCGGATTAGATACCCGGGTAGTCCCAGCCGTAAACGATGCGAGCTAGGTGATGGAATGGCTGCGTGCCATTCCAGTGCCGCAGGGAAGCCGTTAAGCTCGCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAAGGGGTGAAGCCTGCGGTTTAATTGGAGTCAACGCCGGGAACCTTACCCAGGGAGACAGCAGGATGAAGGTCAAGCTAAAGACTTTACCAGACAAGCTGAGAGGTGGTGCATGGCCGTCGCCAGCTCGTGCCGTGAGGCGTCCTGTTAAGTCAGGCAACGAGCGAGACCCTTGCCGCTAGTTGCTATCTCCATCCGAAAGGATCGAGAGGCTAATTAGCGGGACCGCTGCCGCTAAGGCGGAGGAAGGTAGGGGCCACGGCAGGTCAGTATGCCCCGAAACCCTGGGGCCACACGCGGGCTGCAATGGTAGGGACAATGGGTTTCGACCCCGAGAGGGGGAGACAATCCCTAAACCCTACCCCAGTTGTGACTGAGGGCTGCAACCCGCCCTCATGAATCTGGAATCCCTAGTAACCGCGTGTCACCATCGCGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCGCTTCACCCGAGTTGGGTTTGGGTGAGGTGGTATCATATTGGTACCATCGAACTCGAATTCGGCAAGGGGGGAGAAGTCGTAACAAGGTAACCAATCACTAGTGCGGCCGCCTGCAGGTCGACCATATGGGAGAGTCCCACGCGTGA
>DQ243733
TACCAGCCCCGCGAGTGGTCGGGACGGTTATTGGGCCTAAAGCGCCCGTAGCCGGCCTAGCAAGTCTCCGCGGAAATCCCCGGGCTCAACCCGGGGGCATGCGGAGATACTGCTAGGCTAGGGGGCGGGAGAGGCCGGGGGTACTCCCGGGGTAGGGGCGAAATCCTATAATCCCGGGAGGACCACCAGTGGCGAAGGCGCCCGGCTGGAACGCGCCCGACGGTGAGGGGCGAAAGCCGGGGGAGCGAACCGGATTAGATACCCGGGTAGGTATGAGCCGCGGGATAAAGTGTTCTCTATCTAAGTGCAAAGATGGAGCCTACTGTACAGATATCAAAGGTTCGCATACACTTACAGATGTAAACCCGAACAGATGCGGGGAACCTATGAAAACCCTAAAGGGTGAGTGTTCCCTTTAGGGCGGGGAGAAGGTCATCCGCGGCCATACCCCGGAGAGTCCCGGCTGTAAACGATGCGGGCTAGGTGTCGGACGGGCTTAGAGCCCGTCCGGTGCCGCAGGGAAGCCGTTAAGCCCGCCGCCTGGGGAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAAGGGGTGGAGCCTGCGGCTTAATTGGAGTCAACGCCAGGAAACTTACCGGGGGCGACAGCAGGATGACGGCCAGGTTAAAGGCCTTGCCCGACGCGCTGAGAGGAGGTGCACGGCCGTCGCCAGCTCGTGCCGTGAGGTGTCCGGTTAAGTCCGGCAACGAGCGAGACCCCCGCCCCTAGTTGCTACCCGGTCCTTCGGGGCCGGGGCACACTAGGGGGACTGCCGGCGTTAAGCCGGAGGAAGGTGGGGGCCACGGCAGGTCAGCATAGGGGCCTATCCCGAGGGAGGGCCCCGCTGAAGCCCCGAAACCCCCGGGCTGCACGCGGGCTACAATGGCCGGGACAGCGGGATCCGACCCCGAAAGGGGGAGGCAATCCCTTAAACCCGGCCTCAGTTGGGATCGAGGGCTGTAACTCGCCCTCGTGAACGCGGAATCCCTAGTAACCGCGCGTCATCATCGCGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCAAGGGC
>EF533958
TACCGGCAGCTCAAGTGATGGCCAATCTTATTGGGCCTAAAGCGTCCGTAGCTGGCCGTGAAAGTTCGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGAAAACTTCACGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCTGGGGTAGGAGTGAAATCCCGTAATCCTGGACGGACCGCCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACAGTGAGGGACGAAAGCTAGGGTCTCGAACCGGATTAGATACCCGGGTAGTCCTAGCCGTAAACAATGTTCACTAGGTGTGGCACAGGCTACGAGCCTGTGCTGTGCCGTAGGGAAGCCGAGAAGTGAACCGCCTGGGAAGTACGTCCGCAAGGATGAAACTTAAAGGAATTGGCGGGGGAGCACTACAACCGGAGGAGCCTGCGGTTTAATTGGACTCAACGCCGGACATCTCACCAGCATCGACTGTAGTAATGACGACCAGGTTGATGACCTTGTCCGAGCTTCAGAGAGGAGGTGCATGGCCGCCGTCAGCTCGTACCGTGAGGCGTCCTGTTAAGTCAGGCAACGAGCGAGACCCGCATCCTTACTTGCCAGCAGCACTGCGAAGTGGCTGGGGACAGTAGGGAGACCGCCGTGGCCAACACGGAGGAAGGAACGGGCAACGGTAGGTCAGTATGCCCCGAATGAGCTGGGCTACACGCGGGCTACAATGGCCAGAACAATGGGTTGCAACCCCGAAAGGGAACGCTAATCTCCGAAATCTGGTCGTAGTTCGGATTGAGGGCTGAAACTCGCCCTCATGAAGCTGGATTCGGTAGTAATCGCATTTCAGAAGAGTGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCAAAGCACCCGAGTGAGGTCCGGATGAGGCCATCGCAAGATGGTCGAATCTGGGCTTCGCAAGGGGGCTTAAGTCGTAACAAGGTAGCCGTAGGGGAATCTGCGGCTGGATCACCTCCTGATGGTCGAATCTGGGCTTCGCAAGGGGGCTTAAGTCGTAACAAGGTAGCCGTAGGGGAATCTGCGGCTGGATCACCTCCT
>DQ278152
TACCAGCTCCTCGAGTGGTCGGGGCGATTATTGGGCCTAAAGCATCCGTAGCCGGTTCCGTAGGTCTTCTGTTAAATCCAACGGCTTAACCGTTGGCCTGCAGAAGATACCGCGGGACTAGGAGGCGGGAGAGGTGGACGGTACTCCACGTGTAGGGGTAAAATCCTTTGATCCGTGGAAGACCACCAGTGGCGAAGGCGGTCCACCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAACCGGATTAGATACCCGGGTAGTCCCAGCCGTAAACGATGCGAGCTAGATGATCCAATCGCAAATCGCGATTGGAGTGTCGCAGGGAAGCCGTTAAGCTCGCCGCCTGGGGAGTACGGTCGCAAGACTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAAGGGGTGAAGCCTGCGGTTCAATTGGAGTCAACGCCAGGAAACTTACCAGGAGAGACAGCAGGTTGAGGGTCAAGCTGAAGACTTTACCTGACGAGCTGAGAGGAGGTGCATGGCCGTCGCCAGCTCGTGCCGTGAGGCGTCCTGTTAAGTCAGGCAACGAGCGAGACCCCTATTGGTAATTGCTAACAGTCCAGAGATGGATTGAGGCTAGTTACCGAGACTGCCGCTGCTAAAGCGGAGGAAGGAGGGGGCCACGGCAGGTCAGTATGCCCCGAAACTCCTGGGCCACACGCGGGCTGCAATGGTAGGGACAATTGGTTCCGACTTCGAAAGAAGGAGGCAATCCGCAAACCCTACCCCAGTTGTGATTGAGGATTGAAACCCGTCCTCATGAATATGGAATCCCTAGTAACCGCGTGTCATCACCGCGCGGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCGCTTCACCCGAGTTTGGCTAGGACGAGGTGACCTCTAATTGGGGTTATCGAATCTATGTTAGGCAAGGGGGGAGAAGCCGTAACAAGGTGGCTGTAGGGGAACCTGCAGCCGGATCACCTCCTTAGTACATACGGACGAACAACAGACGTACGTACTGAAATCGGGGGTTGGAATACGCACTGCAATGGAGTGAAAAGCAAATGCCAGTGTCAGACGTTGACAGACAGATGAAGATGCACTGTCTCACGTGAGTGAGATAGGAAGGGCTAGCAGGGAACGGAAGGGCCTTCGCGGGCGCCCCTGCGATGTTGCCATGCATAGGGTTCGCGCAAGCATCCGCCGGGGATGTTAATCATCCGACACAACGACGCTGTATGGTGAATGGCTTGGCTTGGT""".split('\n')

expected_reads_forward_primer_only_real_seqs = """>AJ937875
TACCAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCGGTCTTGTAGGTCTCCTGTTAAACTCAACGGCTCAACCGTTGAACTGCAGGAGATACCGCAGGACTAGGAGGCGGGAGAGGCGGACGGTACTCCATGGGTAGCGGCAAAATGTTCTGATCCATGGAAGACCACCAGTGGCGTAGGCGGTCCGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAAC
>DQ243733
TACCAGCCCCGCGAGTGGTCGGGACGGTTATTGGGCCTAAAGCGCCCGTAGCCGGCCTAGCAAGTCTCCGCGGAAATCCCCGGGCTCAACCCGGGGGCATGCGGAGATACTGCTAGGCTAGGGGGCGGGAGAGGCCGGGGGTACTCCCGGGGTAGGGGCGAAATCCTATAATCCCGGGAGGACCACCAGTGGCGAAGGCGCCCGGCTGGAACGCGCCCGACGGTGAGGGGCGAAAGCCGGGGGAGCGAAC
>EF533958
TACCGGCAGCTCAAGTGATGGCCAATCTTATTGGGCCTAAAGCGTCCGTAGCTGGCCGTGAAAGTTCGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGAAAACTTCACGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCTGGGGTAGGAGTGAAATCCCGTAATCCTGGACGGACCGCCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACAGTGAGGGACGAAAGCTAGGGTCTCGAA
>DQ278152
TACCAGCTCCTCGAGTGGTCGGGGCGATTATTGGGCCTAAAGCATCCGTAGCCGGTTCCGTAGGTCTTCTGTTAAATCCAACGGCTTAACCGTTGGCCTGCAGAAGATACCGCGGGACTAGGAGGCGGGAGAGGTGGACGGTACTCCACGTGTAGGGGTAAAATCCTTTGATCCGTGGAAGACCACCAGTGGCGAAGGCGGTCCACCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAAC""".split('\n')

expected_amplicons_reverse_primer_only_real_seqs = """>AJ937875
TTCCGGTTGATCCTGCCGGACCCGACTGCTATCGGAGTGGGGCTAAGCCATGCGAGTCTTACGTCTCTTGACTGCAAGAGGCGTGGCGACCGGCTCAGTAACACGTAGCCAACCTAACCTAGAGACGTGGACAACCCCGGGAAACTGGGATTAATCCACGATAGATCAAGGCTCCTGGAATGGGCTTTGGTCCAAAAGGATGGATAAGCATGTTTTGTCCATCTGCTCTAGGATGGGGCTGCGGCCGATCAGGTTGTTGATGGGGTAATGGCCCACCAAGCCTATAACCGGTACGGGCTTTGAGAGAAGGAGCCCGGAGATGGACACTGAGATAAGGGTCCAGGCCCTACGGGGCGCAGCAGGCGCGAAACCTTCGCAATGCGCGAAAGCGTGACGAGGTTACTTCGAGTGCCACTCGCTGAGGGTGGCTTTTCCCAAATCTAAAAAGCTTGGGGAATAAGGGGAGGGTAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCGGTCTTGTAGGTCTCCTGTTAAACTCAACGGCTCAACCGTTGAACTGCAGGAGATACCGCAGGACTAGGAGGCGGGAGAGGCGGACGGTACTCCATGGGTAGCGGCAAAATGTTCTGATCCATGGAAGACCACCAGTGGCGTAGGCGGTCCGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAACCGG
>EF533958
TCCGGTTGATCCTGCCGGAGGCCATTGCTATCGGAGTCCGACTTAGCCATGCTAGTTGCGCGAGTTTAGACTCGCAGCAGATAGCTCAGTAACACGTGGTCAAGCTACCCTGCAGACACGGACAACCTCGGGAAACTGAGGCTAATCCGCGATATCGATCCCACGCTGGAACAGCCGGGATCAGCAAACGCTCCGGCGCTGCAGGATGCGGCTGCGGCCGATTAGGTAGACGGTGGGGTAACGGCCCACCGTGCCGATAATCGGTACGGGTTGTGAGAGCAAGAGCCCGGAGACGGGATCTGAGACAAGATTCCGGGCCCTACGGGGCGCAGCAGGCGCGAAACCTTTACACTGCACGCAAGTGCGATAAGGGGACTCCGAGTGCGGGGGCATATAGTCCTCGCTTTTGAGAACCGTAAGGCGGTTCTCGAATAAGAGCTGGGCAAGACCGGTGCCAGCCGCCGCGGTAATACCGGCAGCTCAAGTGATGGCCAATCTTATTGGGCCTAAAGCGTCCGTAGCTGGCCGTGAAAGTTCGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGAAAACTTCACGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCTGGGGTAGGAGTGAAATCCCGTAATCCTGGACGGACCGCCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACAGTGAGGGACGAAAGCTAGGGTCTCGAACCGG""".split('\n')

expected_reads_reverse_primer_only_real_seqs = """>AJ937875
CAGCTCCCCGAGTGGTTGGGACATTTATTGGGCCTAAAGCATTCGTAGCCGGTCTTGTAGGTCTCCTGTTAAACTCAACGGCTCAACCGTTGAACTGCAGGAGATACCGCAGGACTAGGAGGCGGGAGAGGCGGACGGTACTCCATGGGTAGCGGCAAAATGTTCTGATCCATGGAAGACCACCAGTGGCGTAGGCGGTCCGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCGAACCGG
>EF533958
GGCAGCTCAAGTGATGGCCAATCTTATTGGGCCTAAAGCGTCCGTAGCTGGCCGTGAAAGTTCGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGAAAACTTCACGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCTGGGGTAGGAGTGAAATCCCGTAATCCTGGACGGACCGCCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACAGTGAGGGACGAAAGCTAGGGTCTCGAACCGG""".split('\n')





if __name__ == "__main__":
    main()
