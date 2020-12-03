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
 
from os.path import basename, isfile, isdir, exists
from os import remove
from shutil import rmtree


from cogent3.util.unit_test import TestCase, main
from primerprospector.old_cogent import get_tmp_filename
from primerprospector.cogentutil.misc import remove_files, create_dir
from primerprospector.cogentutil.misc import get_random_directory_name

from primerprospector.amplicons_histograms import get_amplicons_lens,\
 generate_histogram_combined, get_histogram_name, sort_combined_amplicons,\
 sort_amplicons, generate_histograms_domains, generate_histograms,\
 max_value_from_bins
 


class AmpliconsHistogramsTests(TestCase):
    """ Unit tests for the amplicons_histograms.py module """
    
    def setUp(self):
        # create the temporary input files that will be used
        
        self._files_to_remove = []

        
        self.taxa_mapping_file = taxa_mapping_file
        self.sample_bacterial_amplicons = sample_bacterial_amplicons
        self.amplicons_file = amplicons_file
        
        self.output_dir_overall_test =\
         get_random_directory_name(prefix = '/tmp/')
        self.output_dir_overall_test += '/'
        
    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if exists(self.output_dir_overall_test):
            rmtree(self.output_dir_overall_test)


    def test_get_amplicons_lens(self):
        """ Returns dictionary of seq ID: seq lens properly """
        
        expected_amplicons_lens =\
         {'AB240712': 33, 'EU434431': 31, 'EF446229': 37}
        
        actual_amplicons_lens =\
         get_amplicons_lens(self.sample_bacterial_amplicons)
        
        self.assertEqual(actual_amplicons_lens, expected_amplicons_lens)
        
        
    def test_generate_histogram_combined(self):
        """ Runs without errors (cannot test graphic output) """
        
        amplicons_data = { 100: 33, 99: 31, 101: 37, 102: 5}
        
        output_dir = '/tmp/'
        
        amp_f = get_tmp_filename(prefix='/test/234f_987r_amplicons_',
         suffix='.fasta')
         
        expected_outfile_path = '/tmp/' +\
         basename(amp_f).split('.')[0] + ".ps"
         
        self._files_to_remove.append(expected_outfile_path)
         
        generate_histogram_combined(amplicons_data, amp_f, output_dir)
        
        # Test for creation of outfile path
        self.assertTrue(isfile(expected_outfile_path))
        
    def test_max_value_from_bins(self):
        """ properly returns largest value from a list of lists """
        
        bins = [[1, 2, 3], [2, 3, 4], [5, 1, 2]]
        
        expected_max = 5
        
        actual_max = max_value_from_bins(bins)
        
        self.assertEqual(actual_max, expected_max)
        
        # Should also handle empty bins
        bins = [[1, 2, 3], [], [1, 3, 3]]
        
        expected_max = 3
        
        actual_max = max_value_from_bins(bins)
        
        self.assertEqual(actual_max, expected_max)
        
        
    def test_get_histogram_name(self):
        """ Properly generates histogram output file name """
        
        amplicons_fp = "test/723f_809r_amplicons.fasta"
        output_dir = "../"
        
        expected_histogram_fp = "..//723f_809r_amplicons.ps"
        
        actual_histogram_fp = get_histogram_name(amplicons_fp, output_dir)
        
        self.assertEqual(actual_histogram_fp, expected_histogram_fp)
        
    def test_sort_combined_amplicons(self):
        """ Properly bins combined amplicon lens into dict """
        
        raw_amplicons_data = {'seq1':10, 'seq2':11, 'seq3':12, 'seq4':11,
         'seq5':13}
         
        expected_amplicons_data = {10:1, 11:2, 12:1, 13:1}
        
        actual_amplicons_data = sort_combined_amplicons(raw_amplicons_data)
        
        self.assertEqual(actual_amplicons_data, expected_amplicons_data)
        
    def test_sort_amplicons(self):
        """ Properly bins amplicon sizes according to domain """
        
        raw_amplicons_data = {'seq_arch1':10, 'seq_arch2':10, 'seq_bact1':7,
         'seq_unknown':11, 'seq_bact2':8, 'seq_euk1':12}
         
        taxa_mapping = { 'seq_arch1':'archaea;crenarchaota;uncultured',
                         'seq_arch2':'archaea;euryarchaota',
                         'seq_bact1':'bacteria;alphaproteobacteria',
                         'seq_bact2':'bacteria;firmicutes;clostridia',
                         'seq_euk1':'eukarya;bacillorhiaphyta',
                         'seq_unknown':'root'}
                         
        expected_arch = { 10:2 }
        expected_bact = { 7:1, 8:1 }
        expected_euk = { 12:1 }
        expected_other = { 11:1 }
        
        actual_arch, actual_bact, actual_euk, actual_other =\
         sort_amplicons(raw_amplicons_data, taxa_mapping)
         
        self.assertEqual(actual_arch, expected_arch)
        self.assertEqual(actual_bact, expected_bact)
        self.assertEqual(actual_euk, expected_euk)
        self.assertEqual(actual_other, expected_other)
        
        # Should raise error if sequence not found in taxa mapping data
        raw_amplicons_data['bad_seq'] = 10
        
        self.assertRaises(KeyError, sort_amplicons, raw_amplicons_data,
         taxa_mapping)
        
    def test_generate_histograms_domains(self):
        """ Runs without errors (cannot test graphic output) """
        
        amplicons_data = [{ 100: 33, 99: 31, 101: 37, 102: 5},
                          { 50: 25, 25: 15},
                          { 25: 15, 16: 11},
                          {}]
        
        output_dir = '/tmp/'
        
        amp_f = get_tmp_filename(prefix='/test/234f_987r_amplicons_',
         suffix='.fasta')
         
        expected_outfile_path = '/tmp/' +\
         basename(amp_f).split('.')[0] + ".ps"
         
        self._files_to_remove.append(expected_outfile_path)
         
        generate_histograms_domains(amplicons_data, amp_f, output_dir)
        
        # Test for creation of outfile path
        self.assertTrue(isfile(expected_outfile_path))
        
    def test_generate_histograms_no_taxa_mapping(self):
        """ Main function, test for lack of errors, creation of files """
        
        
        
        create_dir(self.output_dir_overall_test)
        
        amp_fp = self.output_dir_overall_test + '234f_987r_amplicons.fasta'
         
         
        expected_outfile_path = self.output_dir_overall_test +\
         basename(amp_fp).split('.')[0] + ".ps"
         
        amp_f = open(amp_fp, "w")
        amp_f.write(self.amplicons_file)
        amp_f.close()
         
        self._files_to_remove.append(amp_fp)
        
        amplicons_filepath = amp_fp
        output_dir = self.output_dir_overall_test
        all_files = False 
        taxa_map_filepath = None
        
        # Should be called without errors
        generate_histograms(amplicons_filepath, output_dir, all_files,
         taxa_map_filepath)
         
        # Test for creation of histogram file
        
        self.assertTrue(isfile(expected_outfile_path))
        
    def test_generate_histograms_no_taxa_mapping_all_paths(self):
        """ Main function, test for lack of errors, creation of files """
        
        
        
        create_dir(self.output_dir_overall_test)
        
        amp_fp = self.output_dir_overall_test + '234f_987r_amplicons.fasta'
         
         
        expected_outfile_path = self.output_dir_overall_test +\
         basename(amp_fp).split('.')[0] + ".ps"
         
        amp_f = open(amp_fp, "w")
        amp_f.write(self.amplicons_file)
        amp_f.close()
         
        self._files_to_remove.append(amp_fp)
        
        amplicons_filepath = self.output_dir_overall_test
        output_dir = self.output_dir_overall_test
        all_files = True
        taxa_map_filepath = None
        
        # Should be called without errors
        generate_histograms(amplicons_filepath, output_dir, all_files,
         taxa_map_filepath)
         
        # Test for creation of histogram file
        
        self.assertTrue(isfile(expected_outfile_path))
        
    def test_generate_histograms_with_taxa_mapping(self):
        """ Main function, test for lack of errors, creation of files """
        
        
        
        create_dir(self.output_dir_overall_test)
        
        amp_fp = self.output_dir_overall_test + '234f_987r_amplicons.fasta'
         
         
        expected_outfile_path = self.output_dir_overall_test +\
         basename(amp_fp).split('.')[0] + ".ps"
         
        amp_f = open(amp_fp, "w")
        amp_f.write(self.amplicons_file)
        amp_f.close()
         
        self._files_to_remove.append(amp_fp)
        
        taxa_fp = self.output_dir_overall_test + 'test_taxa_mapping.txt'
         
         
         
        taxa_f = open(taxa_fp, "w")
        taxa_f.write(self.taxa_mapping_file)
        taxa_f.close()
         
        self._files_to_remove.append(taxa_fp)
        
        amplicons_filepath = amp_fp
        output_dir = self.output_dir_overall_test
        all_files = False
        taxa_map_filepath = taxa_fp
        
        # Should be called without errors
        generate_histograms(amplicons_filepath, output_dir, all_files,
         taxa_map_filepath)
         
        # Test for creation of histogram file
        
        self.assertTrue(isfile(expected_outfile_path))
    
    
# Placing large strings at the end for better readability
# 3 bacterial sequences from the Silva database.
sample_bacterial_amplicons = """>AB240712
AGAGTTTGATCCTGGCTCAGGGTGAACGCTGGC
>EF446229
TACGGTTACCCTTGTTACGACTTAAAGCCCTAACTTT
>EU434431
AATCTGCACACACCTTCTCAGACTTCCTACG""".split('\n')


amplicons_file = """>seq1
AGAGTTTGATCCTGGCTCAGGGTGAACGCTGGC
>seq2
TACGGTTACCCTTGTTACGACTTAAAGCCCTAACTTT
>seq3
AATCTGCACACACCTTCTCAGACTTCCTACG"""

taxa_mapping_file = """# Test mapping file
seq1\tarchaea;crenarchaota
seq2\tbacteria;firmicutes
seq3\teukarya;chordata"""


    
if __name__ == "__main__":
    main()
