#!/usr/bin/env python
# File created on 26 Apr 2010
from __future__ import division

__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"
 
from os.path import basename
from os import remove, rmdir

from cogent.util.unit_test import TestCase, main
from cogent.app.util import get_tmp_filename
from cogent.util.misc import remove_files, create_dir, get_random_directory_name

from primerprospector.parse import parse_primer_data, get_known_primers,\
 get_primers_data, parse_formatted_primers_data, get_fasta_filepaths,\
 build_fasta_data, get_primer_hits_data_pair, get_primer_hits_data,\
 parse_taxa_mapping_file, get_amplicons_filepaths, get_hits_files, \
 parse_hits_data, get_barcodes, get_hits_field, get_primer_seq



class ParseTests(TestCase):
    """ Unit tests for the parse.py module """
    
    def setUp(self):
        # create the temporary input files that will be used with the
        # get_fasta_filepath function
        self.input_fasta1_fp = get_tmp_filename(\
         prefix = 'input_fasta1',
         suffix = '.fasta')
        seq_file = open(self.input_fasta1_fp, 'w')
        seq_file.write(sample_fasta_file1)
        seq_file.close()
        
        self.input_fasta2_fp = get_tmp_filename(\
         prefix = 'input_fasta2',
         suffix = '.fasta')
        seq_file = open(self.input_fasta2_fp, 'w')
        seq_file.write(sample_fasta_file2)
        seq_file.close()
        
        self.bad_fasta1_fp = get_tmp_filename(\
         prefix = 'bad_fasta1',
         suffix = '.fasta')
        seq_file = open(self.bad_fasta1_fp, 'w')
        seq_file.write(bad_fasta_file1)
        seq_file.close()
        
        self.bad_fasta2_fp = get_tmp_filename(\
         prefix = 'bad_fasta2',
         suffix = '.fasta')
        seq_file = open(self.bad_fasta2_fp, 'w')
        seq_file.write(bad_fasta_file2)
        seq_file.close()
        
        self.sample_hits_data = sample_hits_data
        self.sample_hits_file = sample_hits_file
        self.expected_sample_f_primers = expected_sample_f_primers
        self.expected_sample_r_primers = expected_sample_r_primers
        self.bad_sample_hits_data = bad_sample_hits_data
        self.known_primers_data = known_primers_data
        self.known_primers_data_bad = known_primers_data_bad
        self.expected_parsed_hits_file = expected_parsed_hits_file
        self.sample_primers_file = sample_primers_file
        self.bad_primers_file = bad_primers_file
        self.expected_forward_reverse_hits_data =\
         expected_forward_reverse_hits_data
        self.forward_hits_file = forward_hits_file
        self.reverse_hits_file = reverse_hits_file
        self.forward_hits_file_diff_len = forward_hits_file_diff_len
        self.reverse_hits_file_diff_ids = reverse_hits_file_diff_ids
        self.expected_data_forward_hits = expected_data_forward_hits
        
        self.sample_taxa_mapping_file = sample_taxa_mapping_file
        self.expected_taxa_map = expected_taxa_map
        self.bad_taxa_mapping_file_not_tabbed = bad_taxa_mapping_file_not_tabbed
        self.bad_taxa_mapping_file_extra_tab = bad_taxa_mapping_file_extra_tab
        
        self.sample_fasta_file1 = sample_fasta_file1
        self.sample_fasta_file2 = sample_fasta_file2
        
        self._files_to_remove =\
         [self.input_fasta1_fp, self.input_fasta2_fp,
          self.bad_fasta1_fp, self.bad_fasta2_fp]
          

                 
    def tearDown(self):
        remove_files(self._files_to_remove)


    def test_get_hits_field(self):
        """ Returns list of items from a particular field of hits files """
        
        sample_hits_lines = ['A01,ACAAGACG,AYAAGACG,10,4.0',
                             'b01,ACTAGACG,AYAAGACG,11,4.0',
                             'd01,ACTAAACG,AYAAGACG,11,4.0']
                             
        expected_fields = ['ACAAGACG', 'ACTAGACG', 'ACTAAACG']
        
        actual_fields = get_hits_field(sample_hits_lines, 1)
        
        self.assertEqual(actual_fields, expected_fields)
        
        # Should raise error with invalid index given
        
        self.assertRaises(IndexError, get_hits_field, sample_hits_lines, 5)
        
    def test_get_primer_seq(self):
        """ Returns primer sequence from hits file """
        
        expected_primer_seq = 'GGACTACVSGGGTATCTAAT'
        
        hits_reverse_f = get_tmp_filename(prefix = 'reverse_', 
         suffix = '_hits.txt')
        hits_f = open(hits_reverse_f, "w")
        self._files_to_remove.append(hits_reverse_f)
        hits_f.write(self.reverse_hits_file)
        hits_f.close()
        
        actual_primer_seq = get_primer_seq(hits_reverse_f)
        
        self.assertEqual(actual_primer_seq, expected_primer_seq)

    def test_parse_primer_data(self):
        """ Properly parses primer data from given list of lines """
        
        # See sample data, expected data at the end of the file
        # Should parse into proper dictionary form
        test_f_primers, test_r_primers =\
         parse_primer_data(self.sample_hits_data)
        
        self.assertEqual(test_f_primers, self.expected_sample_f_primers)
        self.assertEqual(test_r_primers, self.expected_sample_r_primers)
        
        # Should raise error if hits file not formatted properly
        self.assertRaises(ValueError, parse_primer_data,
         self.bad_sample_hits_data)
        
    def test_get_known_primers(self):
        """ Correctly parses known primers data """
        
        # Should build list of tuples, with seq, primer name
        expected_known_primers = \
         [("ATTCCGGTTGATCCTGC","1f"),
          ("TTCCGGTTGATCCYGCCGGA","21f"),
          ("AGAGTTTGATCMTGGCTCAG","27f")]
         
        result_known_primers = get_known_primers(self.known_primers_data)
        
        self.assertEqual(result_known_primers, expected_known_primers)
        
        #Should raise ValueError if file format incorrect
        
        self.assertRaises(ValueError, get_known_primers,
         self.known_primers_data_bad)
        
    def test_get_primers_data(self):
        """ Correctly parses primer output from generate_primers_denovo.py """
        
        # Should ignore comment lines, empty lines
        actual_parsed_hits_file = get_primers_data(self.sample_hits_file)
        
        self.assertEqual(actual_parsed_hits_file,
         self.expected_parsed_hits_file)
        
    def test_parse_formatted_primers_data(self):
        """ Correctly parses formatted primers file """
        
        expected_result = [("123f","AACGATT"),("428r","ATCCAGGRC")]
        
        actual_result = parse_formatted_primers_data(self.sample_primers_file)
        
        self.assertEqual(actual_result, expected_result)
        
        # Raises error with incorrect format
        self.assertRaises(IndexError, \
         parse_formatted_primers_data, self.bad_primers_file)
        
    def test_get_fasta_filepaths(self):
        """ Should return list of filepaths, test validity of fasta data """
        
        # Test with a single file
        expected_result = [self.input_fasta1_fp]
        
        actual_result = get_fasta_filepaths(self.input_fasta1_fp)
        
        self.assertEqual(actual_result, expected_result)
        
        # Test with multiple files
        
        expected_result = [self.input_fasta1_fp, self.input_fasta2_fp]
        
        actual_result = get_fasta_filepaths(self.input_fasta1_fp + ":" + \
         self.input_fasta2_fp)
         
        self.assertEqual(actual_result, expected_result)
        
        # Raises error with invalid fasta file (gaps or U character)
        
        # first file has U character
        
        self.assertRaises(ValueError, get_fasta_filepaths, self.bad_fasta1_fp)
        
        # second file has gap characters
        
        self.assertRaises(ValueError, get_fasta_filepaths, self.bad_fasta2_fp)
        
    def test_build_fasta_data(self):
        """ Builds fasta dictionary from fasta files passed correctly"""
        
        # Should build dictionary from passed fasta filepaths
        expected_fasta_dict = {'seq1':'AATCGGT', 'seq2':'CAATCGT', 
         'seq3':'AATCTTTACG', 'seq4':'CTCG'}
         
        actual_fasta_dict = build_fasta_data([self.input_fasta1_fp, 
         self.input_fasta2_fp])
         
        self.assertEqual(actual_fasta_dict, expected_fasta_dict)
         
        
        
        
    def test_get_primer_hits_data_pair(self):
        """ Builds 2d list of data for primer hits pair correctly """
        
        hits_forward_f = get_tmp_filename(prefix = 'forward_hits', 
         suffix = '.txt')
        hits_f = open(hits_forward_f, "w")
        self._files_to_remove.append(hits_forward_f)
        hits_f.write(self.forward_hits_file)
        hits_f.close()
        
        hits_reverse_f = get_tmp_filename(prefix = 'reverse_hits', 
         suffix = '.txt')
        hits_f = open(hits_reverse_f, "w")
        self._files_to_remove.append(hits_reverse_f)
        hits_f.write(self.reverse_hits_file)
        hits_f.close()
        
        expected_hits_data = self.expected_forward_reverse_hits_data
        actual_hits_data = get_primer_hits_data_pair([hits_forward_f, 
         hits_reverse_f])
         
        self.assertEqual(actual_hits_data, expected_hits_data)
        
    def test_get_primer_hits_data_pair_len_diff(self):
        """ get_primer_hits raises errors if hits data incongruous """
        
        hits_forward_f = get_tmp_filename(prefix = 'forward_hits', 
         suffix = '.txt')
        hits_f = open(hits_forward_f, "w")
        self._files_to_remove.append(hits_forward_f)
        hits_f.write(self.forward_hits_file_diff_len)
        hits_f.close()
        
        hits_reverse_f = get_tmp_filename(prefix = 'reverse_hits', 
         suffix = '.txt')
        hits_f = open(hits_reverse_f, "w")
        self._files_to_remove.append(hits_reverse_f)
        hits_f.write(self.reverse_hits_file)
        hits_f.close()
        
        # should raise error due to differing amounts of hits data
        self.assertRaises(ValueError, get_primer_hits_data_pair, 
         [hits_forward_f, hits_reverse_f])
         
    def test_get_primer_hits_data_pair_diff_ids(self):
        """ get_primer_hits raises errors if hits data incongruous """
        
        hits_forward_f = get_tmp_filename(prefix = 'forward_hits', 
         suffix = '.txt')
        hits_f = open(hits_forward_f, "w")
        self._files_to_remove.append(hits_forward_f)
        hits_f.write(self.forward_hits_file)
        hits_f.close()
        
        hits_reverse_f = get_tmp_filename(prefix = 'reverse_hits', 
         suffix = '.txt')
        hits_f = open(hits_reverse_f, "w")
        self._files_to_remove.append(hits_reverse_f)
        hits_f.write(self.reverse_hits_file_diff_ids)
        hits_f.close()
        
        # should raise error due to IDs not matching
        self.assertRaises(ValueError, get_primer_hits_data_pair, 
         [hits_forward_f, hits_reverse_f])
        
    def test_get_primer_hits_data(self):
        """ Builds list of data from single primer hits file correctly """
        
        hits_forward_f = get_tmp_filename(prefix = 'forward_hits', 
         suffix = '.txt')
        hits_f = open(hits_forward_f, "w")
        self._files_to_remove.append(hits_forward_f)
        hits_f.write(self.forward_hits_file)
        hits_f.close()
        
        expected_data = self.expected_data_forward_hits
        
        actual_data = get_primer_hits_data(hits_forward_f)
        
        self.assertEqual(actual_data, expected_data)
        
    def test_parse_taxa_mapping_file(self):
        """ Properly builds dictionary mapping seq ID to taxa """
        
        # Test for proper dictionary format, comment/blank line skipping
        expected_result = self.expected_taxa_map
        
        actual_taxa_map = parse_taxa_mapping_file(self.sample_taxa_mapping_file)
        
        self.assertEqual(actual_taxa_map, expected_result)
        
        # Should raise error if non-comment line has no tab
        
        self.assertRaises(ValueError, parse_taxa_mapping_file, 
         self.bad_taxa_mapping_file_not_tabbed)
         
        # Should raise error if extra tab in mapping file
        
        self.assertRaises(ValueError, parse_taxa_mapping_file,
         self.bad_taxa_mapping_file_extra_tab)
         
    def test_get_amplicons_filepaths(self):
        """ Builds list of filepaths for amplicons file correctly """
        
        amp_fp1 = get_tmp_filename(prefix = 'amp_1_', 
         suffix = '_amplicons.fasta')
        self._files_to_remove.append(amp_fp1)
        amp_f1 = open(amp_fp1, 'w')
        amp_f1.write(self.sample_fasta_file1)
        amp_f1.close()
        
        
        
        amp_fp2 = get_tmp_filename(prefix = 'amp_2_', 
         suffix = '_amplicons.fasta')
        self._files_to_remove.append(amp_fp2)
        amp_f2 = open(amp_fp2, 'w')
        amp_f2.write(self.sample_fasta_file2)
        amp_f2.close()
        
        # First test for a single filepath return
        expected_amp_fp = [amp_fp1]
        
        actual_amp_fp = get_amplicons_filepaths(amp_fp1, all_files=False)
        
        self.assertEqual(actual_amp_fp, expected_amp_fp)
        
        # Test for all amplicon filepaths
        
        expected_amp_fp = [amp_fp1, amp_fp2]
        
        actual_amp_fp = get_amplicons_filepaths('/tmp/', all_files=True)
        
        # Test to make sure all filepaths present, may be in different order
        for fp in actual_amp_fp:
            self.assertTrue(fp in expected_amp_fp)
            
        # Should raise IOError for invalid filepath
        
        bad_filepath = get_tmp_filename(prefix = 'amp_bad_', 
         suffix = '_amplicons.fasta')
        
        self.assertRaises(IOError, get_amplicons_filepaths, bad_filepath)
        
    def test_get_hits_files(self):
        """ Returns list of hits filepaths correctly """
        
        # Should return single element list if single file specified
        
        expected_hits_fps = ['/test/23f_hits.txt']
        
        actual_hits_fps = get_hits_files('/test/23f_hits.txt', all_files=False)
        
        self.assertEqual(actual_hits_fps, expected_hits_fps)
        
        # Handle 2 filepaths separated by a colon
        expected_hits_fps = ['test_1f_hits.txt', '../test_300r_hits.txt']
        
        actual_hits_fps =\
         get_hits_files('test_1f_hits.txt:../test_300r_hits.txt', 
         all_files=False)
         
        self.assertEqual(actual_hits_fps, expected_hits_fps)
        
        
        hits_forward_f = get_tmp_filename(prefix = 'forward_', 
         suffix = '_hits.txt')
        hits_f = open(hits_forward_f, "w")
        self._files_to_remove.append(hits_forward_f)
        hits_f.write(self.forward_hits_file)
        hits_f.close()
        
        hits_reverse_f = get_tmp_filename(prefix = 'reverse_', 
         suffix = '_hits.txt')
        hits_f = open(hits_reverse_f, "w")
        self._files_to_remove.append(hits_reverse_f)
        hits_f.write(self.reverse_hits_file)
        hits_f.close()
        
        expected_hits_fps = [hits_forward_f, hits_reverse_f]
        
        actual_hits_fps = get_hits_files('/tmp/', all_files = True)
        
        
        # Test for list of exactly 2 elements
        self.assertTrue(len(actual_hits_fps) == 2)
        
        for hits_fps in actual_hits_fps:
            self.assertTrue(hits_fps in expected_hits_fps)
            
            
    def test_parse_hits_data(self):
        """ Returns list of lines according to score parameters correctly """
        
        test_hits_data = self.expected_data_forward_hits
        
        score_type = "weighted_score"
        score_threshold = 1.0
        
        expected_result = ['AJ937875,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,476,1,0,False,0,0,0.4,False']
        
        actual_result = parse_hits_data(test_hits_data, 
         score_threshold, score_type)
         
        self.assertEqual(actual_result, expected_result)
        
        score_type = "tp_mismatches"
        score_threshold = 0
        
        expected_result = ['AJ937875,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,476,1,0,False,0,0,0.4,False']
        
        actual_result = parse_hits_data(test_hits_data, 
         score_threshold, score_type)
         
        self.assertEqual(actual_result, expected_result)
        
        score_type = "overall_mismatches"
        score_threshold = 3
        
        expected_result = ['AB329763,GTGCCAGCCGCCACGGGTA,GTGCCAGCMGCCGCGGTAA,483,1,2,False,0,0,2.4,False\n', 'AJ937875,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,476,1,0,False,0,0,0.4,False']
        
        actual_result = parse_hits_data(test_hits_data, 
         score_threshold, score_type)
         
        self.assertEqual(actual_result, expected_result)
        
        score_type = "overall_mismatches"
        score_threshold = 2
        
        expected_result = ['AJ937875,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,476,1,0,False,0,0,0.4,False']
        
        actual_result = parse_hits_data(test_hits_data, 
         score_threshold, score_type)
         
        self.assertEqual(actual_result, expected_result)
        
    def test_get_barcodes(self):
        """ Parses barcodes file correctly """
        
        expected_barcodes = ['AACATAGCCAG', 'AGGACAGATAA', 'ACAGGATAAGA']
        
        actual_barcodes = get_barcodes(sample_barcodes_file)
        
        self.assertEqual(actual_barcodes, expected_barcodes)
        
        
        
        
        
        
        
        
sample_hits_data = ['C,TAACA,51,2041,9,90.0%,10.0%,102,137', 'M,2041,TGGCATACCGCTCTGTAACA,TAACACGTAGATAACCTACC,AY800210', 'M,2041,TGGCGTACGGCTCAGTAACA,TAACACGTGGATAACTTACC,EU883771', 'M,2041,TGGCATACGGCTCAGTAACA,TAACACGTAGTCAACATGCC,EF503699', 'M,2041,GGCATACAGGCTCAGTAACA,TAACACGTAGTCAACATGCC,EF503697', 'M,2041,AGGCGGACGGCTCAGTAACA,TAACACGCAGTCAACCTAAC,AJ535128', 'M,2041,TAGCATACGGCTCAGTAACA,TAACACGTAGTCAACCTGCC,EF021349', 'M,2041,CGGCGAACGGCTCAGTAACA,TAACACGTGGGTAATCTGCC,AJ969787', 'M,2041,CGGCATACGGCTCAGTAACA,TAACACGTGGGTAACCTGCC,AY360086', 'M,2041,AGGCGGACTGCTCAGTAACA,TAACACGTGGGTAATCTACC,DQ641863', 'C,AACAC,52,2043,9,90.0%,10.0%,103,138', 'M,2043,GGCATACCGCTCTGTAACAC,AACACGTAGATAACCTACCC,AY800210', 'M,2043,GGCGTACGGCTCAGTAACAC,AACACGTGGATAACTTACCC,EU883771', 'M,2043,GGCATACGGCTCAGTAACAC,AACACGTAGTCAACATGCCC,EF503699', 'M,2043,GCATACAGGCTCAGTAACAC,AACACGTAGTCAACATGCCC,EF503697', 'M,2043,GGCGGACGGCTCAGTAACAC,AACACGCAGTCAACCTAACT,AJ535128', 'M,2043,AGCATACGGCTCAGTAACAC,AACACGTAGTCAACCTGCCC,EF021349', 'M,2043,GGCGAACGGCTCAGTAACAC,AACACGTGGGTAATCTGCCC,AJ969787', 'M,2043,GGCATACGGCTCAGTAACAC,AACACGTGGGTAACCTGCCC,AY360086', 'M,2043,GGCGGACTGCTCAGTAACAC,AACACGTGGGTAATCTACCC,DQ641863']
expected_sample_f_primers = {'2043': ['AACAC,52,2043,9,90.0%,10.0%,103,138', 'GGCATACCGCTCTGTAACAC', 'GGCGTACGGCTCAGTAACAC', 'GGCATACGGCTCAGTAACAC', 'GCATACAGGCTCAGTAACAC', 'GGCGGACGGCTCAGTAACAC', 'AGCATACGGCTCAGTAACAC', 'GGCGAACGGCTCAGTAACAC', 'GGCATACGGCTCAGTAACAC', 'GGCGGACTGCTCAGTAACAC'], '2041': ['TAACA,51,2041,9,90.0%,10.0%,102,137', 'TGGCATACCGCTCTGTAACA', 'TGGCGTACGGCTCAGTAACA', 'TGGCATACGGCTCAGTAACA', 'GGCATACAGGCTCAGTAACA', 'AGGCGGACGGCTCAGTAACA', 'TAGCATACGGCTCAGTAACA', 'CGGCGAACGGCTCAGTAACA', 'CGGCATACGGCTCAGTAACA', 'AGGCGGACTGCTCAGTAACA']}
expected_sample_r_primers = {'2043': ['AACAC,52,2043,9,90.0%,10.0%,103,138', 'AACACGTAGATAACCTACCC', 'AACACGTGGATAACTTACCC', 'AACACGTAGTCAACATGCCC', 'AACACGTAGTCAACATGCCC', 'AACACGCAGTCAACCTAACT', 'AACACGTAGTCAACCTGCCC', 'AACACGTGGGTAATCTGCCC', 'AACACGTGGGTAACCTGCCC', 'AACACGTGGGTAATCTACCC'], '2041': ['TAACA,51,2041,9,90.0%,10.0%,102,137', 'TAACACGTAGATAACCTACC', 'TAACACGTGGATAACTTACC', 'TAACACGTAGTCAACATGCC', 'TAACACGTAGTCAACATGCC', 'TAACACGCAGTCAACCTAAC', 'TAACACGTAGTCAACCTGCC', 'TAACACGTGGGTAATCTGCC', 'TAACACGTGGGTAACCTGCC', 'TAACACGTGGGTAATCTACC']}
# Bad because first field should be 'C' or 'M'
bad_sample_hits_data = ['X,TAACA,51,2041,9,90.0%,10.0%,102,137', 'M,2041,TGGCATACCGCTCTGTAACA,TAACACGTAGATAACCTACC,AY800210', 'M,2041,TGGCGTACGGCTCAGTAACA,TAACACGTGGATAACTTACC,EU883771', 'M,2041,TGGCATACGGCTCAGTAACA,TAACACGTAGTCAACATGCC,EF503699', 'M,2041,GGCATACAGGCTCAGTAACA,TAACACGTAGTCAACATGCC,EF503697', 'M,2041,AGGCGGACGGCTCAGTAACA,TAACACGCAGTCAACCTAAC,AJ535128', 'M,2041,TAGCATACGGCTCAGTAACA,TAACACGTAGTCAACCTGCC,EF021349', 'M,2041,CGGCGAACGGCTCAGTAACA,TAACACGTGGGTAATCTGCC,AJ969787', 'M,2041,CGGCATACGGCTCAGTAACA,TAACACGTGGGTAACCTGCC,AY360086', 'M,2041,AGGCGGACTGCTCAGTAACA,TAACACGTGGGTAATCTACC,DQ641863', 'C,AACAC,52,2043,9,90.0%,10.0%,103,138', 'M,2043,GGCATACCGCTCTGTAACAC,AACACGTAGATAACCTACCC,AY800210', 'M,2043,GGCGTACGGCTCAGTAACAC,AACACGTGGATAACTTACCC,EU883771', 'M,2043,GGCATACGGCTCAGTAACAC,AACACGTAGTCAACATGCCC,EF503699', 'M,2043,GCATACAGGCTCAGTAACAC,AACACGTAGTCAACATGCCC,EF503697', 'M,2043,GGCGGACGGCTCAGTAACAC,AACACGCAGTCAACCTAACT,AJ535128', 'M,2043,AGCATACGGCTCAGTAACAC,AACACGTAGTCAACCTGCCC,EF021349', 'M,2043,GGCGAACGGCTCAGTAACAC,AACACGTGGGTAATCTGCCC,AJ969787', 'M,2043,GGCATACGGCTCAGTAACAC,AACACGTGGGTAACCTGCCC,AY360086', 'M,2043,GGCGGACTGCTCAGTAACAC,AACACGTGGGTAATCTACCC,DQ641863']

known_primers_data = ['# Current Univeral Primers',
"# primer_id <tab> primer sequence (5'->3')  <tab> primer citation'",
"1f\tATTCCGGTTGATCCTGC\tBaker and Cowan, 2004",
"21f\tTTCCGGTTGATCCYGCCGGA\tEllis 2008, DeLong 1992 (also listed as F1 in Baker and Cowan, 2004)",
"27f\tAGAGTTTGATCMTGGCTCAG"]

known_primers_data_bad = ['# Current Univeral Primers',
"# primer_id <tab> primer sequence (5'->3')  <tab> primer citation'",
"1f\tATTCCGGTTGATCCTGC\tBaker and Cowan, 2004",
"21f\tTTCCGGTTGATCCYGCCGGA\tEllis 2008, DeLong 1992 (also listed as F1 in Baker and Cowan, 2004)",
"27f\tAGAGTTTGATCMTGGCTCAG",
"323bad\tACCGATTCARG"]

sample_hits_file = ["# Matches that are found in at least 95.00% target sequences",
"# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches,5' standard index, 3' standard index",
"",
"C,AAACT,135,3105,10,100.0%,0.0%,144,179",
"M,3105,GGGGAAACTCCCGGGAAACT,AAACTGGGCCTAATCCCCGA,AY800210",
"M,3105,TGGGATAACTCTGGGAAACT,AAACTGGGGATAATACTGGA,EU883771",
"C,GTAGC,514,15657,10,100.0%,0.0%,565,600",
"M,15657,GGGCCTAAAGCGTCCGTAGC,GTAGCCGGGCGTGCAAGTCA,AY800210",
"M,15657,GGGCTTAAAGCGTTCGTAGC,GTAGCTTGATTTTTAAGTCT,EU883771",
"# comment line"]

expected_parsed_hits_file = ["C,AAACT,135,3105,10,100.0%,0.0%,144,179",
"M,3105,GGGGAAACTCCCGGGAAACT,AAACTGGGCCTAATCCCCGA,AY800210",
"M,3105,TGGGATAACTCTGGGAAACT,AAACTGGGGATAATACTGGA,EU883771",
"C,GTAGC,514,15657,10,100.0%,0.0%,565,600",
"M,15657,GGGCCTAAAGCGTCCGTAGC,GTAGCCGGGCGTGCAAGTCA,AY800210",
"M,15657,GGGCTTAAAGCGTTCGTAGC,GTAGCTTGATTTTTAAGTCT,EU883771"]

sample_primers_file = ["#Primer\tSequence\tCitation",
"123f\tAACGATT\t",
"428r\tATCCAGGRC\t"]

bad_primers_file = ["# Bad format ",
"123f AACGATT",
"428r\tATCCAGGRC\t"]

sample_fasta_file1 = """>seq1
AATCGGT
>seq2
CAATCGT"""

sample_fasta_file2 = """>seq3
AATCTTTACG
>seq4
CTCG"""

bad_fasta_file1 = """>seq5
AAAUAACG"""

bad_fasta_file2 = """>seq6
..---AACGTGC--.."""

forward_hits_file = """# Primer: 515f 5'-GTGCCAGCMGCCGCGGTAA-3'
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
AJ937875,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,476,1,0,False,0,0,0.4,False"""

reverse_hits_file = """# Primer: 806r 5'-GGACTACVSGGGTATCTAAT-3'
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
AJ937875,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,748,0,0,False,0,0,0.0,False"""

expected_forward_reverse_hits_data = [['AB329763,GTGCCAGCCGCCACGGGTA,GTGCCAGCMGCCGCGGTAA,483,1,2,False,0,0,2.4,False\n', 'AJ937875,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,476,1,0,False,0,0,0.4,False'], ['AB329763,ACTAAAAACCCCGGTAGTCC,ATTAGATACCCSBGTAGTCC,823,1,2,False,0,0,2.4,False\n', 'AJ937875,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,748,0,0,False,0,0,0.0,False']]

forward_hits_file_diff_len = """# Primer: 515f 5'-GTGCCAGCMGCCGCGGTAA-3'
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
DQ243733,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,414,1,0,False,0,0,0.4,False"""

reverse_hits_file_diff_ids = """# Primer: 806r 5'-GGACTACVSGGGTATCTAAT-3'
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
XXXXXXXX,ACTAAAAACCCCGGTAGTCC,ATTAGATACCCSBGTAGTCC,823,1,2,False,0,0,2.4,False
AJ937875,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,748,0,0,False,0,0,0.0,False"""

expected_data_forward_hits = ['AB329763,GTGCCAGCCGCCACGGGTA,GTGCCAGCMGCCGCGGTAA,483,1,2,False,0,0,2.4,False\n', 'AJ937875,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,476,1,0,False,0,0,0.4,False']

sample_taxa_mapping_file = """AB175604\tArchaea;Euryarchaeota;Halobacteriales;uncultured
AB177021\tArchaea;Crenarchaeota;uncultured;uncultured
AB302039\tArchaea;uncultured
# Bacteria

AJ318104\tBacteria;Actinobacteria;Microsphaera
EU801979\tBacteria;Actinobacteria;Sporichtya;uncultured
AJ441233\tBacteria;uncultured OP8""".split('\n')

expected_taxa_map = {'EU801979': 'Bacteria;Actinobacteria;Sporichtya;uncultured', 'AJ441233': 'Bacteria;uncultured OP8', 'AB302039': 'Archaea;uncultured', 'AJ318104': 'Bacteria;Actinobacteria;Microsphaera', 'AB175604': 'Archaea;Euryarchaeota;Halobacteriales;uncultured', 'AB177021': 'Archaea;Crenarchaeota;uncultured;uncultured'}

bad_taxa_mapping_file_not_tabbed = """AB175604\tArchaea;Euryarchaeota;Halobacteriales;uncultured
AB177021\tArchaea;Crenarchaeota;uncultured;uncultured
AB302039 Archaea;uncultured
# Bacteria

AJ318104\tBacteria;Actinobacteria;Microsphaera
EU801979\tBacteria;Actinobacteria;Sporichtya;uncultured
AJ441233\tBacteria;uncultured OP8""".split('\n')

bad_taxa_mapping_file_extra_tab = """AB175604\tArchaea;Euryarchaeota;Halobacteriales;uncultured
AB177021\tArchaea;Crenarchaeota;uncultured;uncultured
AB302039\t\tArchaea;uncultured
# Bacteria

AJ318104\tBacteria;Actinobacteria;Microsphaera
EU801979\tBacteria;Actinobacteria;Sporichtya;uncultured
AJ441233\tBacteria;uncultured OP8""".split('\n')


sample_barcodes_file = """# Test barcodes
AACATAGCCAG
AGGACAGATAA

ACAGGATAAGA
""".split('\n')
if __name__ == "__main__":
    main()
