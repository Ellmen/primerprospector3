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

from os.path import basename, isdir, isfile, exists
from shutil import rmtree

from cogent3.util.unit_test import TestCase, main
from primerprospector.old_cogent import get_tmp_filename
from primerprospector.cogentutil.misc import remove_files, create_dir, get_random_directory_name

from primerprospector.taxa_coverage import get_coverage_data_primer_pair,\
 get_coverage_data_single_primer, get_coverage_filepath,\
 get_coverage_graph_root_filepath, get_graph_fp, get_log_filepath,\
 get_output_dir, get_primer_pairs, graph_taxa_coverage, write_log_file,\
 write_primer_coverage, write_primer_coverage_graph, get_curr_taxa




class TaxaCoverageTests(TestCase):
    """unit tests for taxa_coverage """
    
    def setUp(self):
        """ Generate temporary files for use in unit testing"""
        
        self.expected_combined_coverage = expected_combined_coverage
        self.expected_combined_log = expected_combined_log
        self.expected_500r_coverage = expected_500r_coverage
        self.expected_500r_log = expected_500r_log
        self.expected_100f_coverage = expected_100f_coverage
        self.expected_100f_log = expected_100f_log
        self.pdf_filepaths = pdf_filepaths
        self.taxa_mapping_data = taxa_mapping_data
        self.taxa_mapping_file = taxa_mapping_file
        self.reverse_hits_data_taxa_unmapped = reverse_hits_data_taxa_unmapped
        self.reverse_hits_data_not_matched = reverse_hits_data_not_matched
        self.reverse_hits_data_tp_mm = reverse_hits_data_tp_mm
        self.reverse_hits_data = reverse_hits_data
        self.reverse_hits_file = reverse_hits_file
        self.forward_hits_data_taxa_unmapped = forward_hits_data_taxa_unmapped
        self.forward_hits_data = forward_hits_data
        self.forward_hits_file = forward_hits_file
        
        self.output_dir = get_random_directory_name(prefix = '/tmp/')
        self.output_dir += '/'
        
        create_dir(self.output_dir)
        
        self.forward_hits = self.output_dir + "/100f_hits.txt"
        f_hits = open(self.forward_hits, "w")
        f_hits.write(self.forward_hits_file)
        f_hits.close()
        
        self.reverse_hits = self.output_dir + "/500r_hits.txt"
        r_hits = open(self.reverse_hits, "w")
        r_hits.write(self.reverse_hits_file)
        r_hits.close()
        
        self.taxa_mapping = self.output_dir + "/taxa_mapping.txt"
        taxa_mapping = open(self.taxa_mapping, "w")
        taxa_mapping.write(self.taxa_mapping_file)
        taxa_mapping.close()
        
        
        self._files_to_remove = [self.forward_hits, self.reverse_hits,
         self.taxa_mapping]
        
    def tearDown(self):
        remove_files(self._files_to_remove)
        if exists(self.output_dir ):
            rmtree(self.output_dir )
            
    def test_get_coverage_data_primer_pair_weighted_score(self):
        """ Properly gets sequence coverage, count for primer pair hits data """
        
        
        
        pair_hits_data = [self.forward_hits_data, self.reverse_hits_data]
        taxa_map = self.taxa_mapping_data
        taxa_depth = 3
        score_type = 'weighted_score'
        score_threshold = 1
        hits_f_fp = '100f_hits.txt'
        hits_r_fp = '../500r_hits.txt'
        taxa_fp = '/home/BigErn/taxa_data/taxa_mapping.txt'
        
        '''These are the taxa that correspond to the sequence IDs given
         in the hits data:
        AB302407  Bacteria;Firmicutes;uncultured
        EF142855  Bacteria;Firmicutes;Clostridia
        DQ278118  Bacteria;Firmicutes;Clostridia
        DQ278125 Archaea;Euryarchaeota;Thermoplasmatales;uncultured
        
        The results in this case should be 1 Archaea domain (and 1 each of its
        descending taxa), and 3 bacteria.  On the second level of taxa, all
        3 bacteria should be under firmicutes.  The third level should have
        2 clostridia, and one uncultured (labeled as its parent taxa + 
        _Uncultured to distinguish it from other uncultured taxa at that level).
        
        All of the sequences should be considered passing, as neither of the
        two hits weighted scores are greater than 1, the threshold'''
        
        expected_coverage = [{('Archaea', 'Archaea'): [1, 1],
         ('Bacteria', 'Bacteria'): [3, 3]},
         {('Bacteria', 'Firmicutes'): [3, 3],
         ('Archaea', 'Euryarchaeota'): [1, 1]},
         {('Bacteria', 'Firmicutes_Uncultured'): [1, 1],
         ('Bacteria', 'Clostridia'): [2, 2],
         ('Archaea', 'Thermoplasmatales'): [1, 1]}]
        
        # Log file should give base names of all input files, scoring 
        # parameters, and unknown taxa (in this case 0)
        expected_log =\
         ['# Taxonomic coverage report\n# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.\n',
          '# Hits files analyzed: 100f_hits.txt,500r_hits.txt\n', 
          '# Taxonomy mapping file used: taxa_mapping.txt\n', 
          '# Scoring method used: weighted_score\n', 
          '# Score threshold: 1.00\n', 
          '# Total sequences in hits files: 4\n', 
          '# Total sequences with no corresponding taxonomy in taxonomy mapping file: 0\n', 
          '# Sequence IDs lacking taxonomy mapping:\n', 
          '']
        
        actual_coverage, actual_log =\
         get_coverage_data_primer_pair(pair_hits_data, taxa_map, taxa_depth,
         score_type, score_threshold, hits_f_fp, hits_r_fp, taxa_fp)
         
        self.assertEqual(actual_coverage, expected_coverage)
        self.assertEqual(actual_log, expected_log)
        
    def test_get_coverage_data_primer_pair_weighted_score_inc_threshold(self):
        """ Properly gets sequence coverage, count for primer pair hits data """
        
        pair_hits_data = [self.forward_hits_data, self.reverse_hits_data]
        taxa_map = self.taxa_mapping_data
        taxa_depth = 3
        score_type = 'weighted_score'
        score_threshold = 0.5
        hits_f_fp = '100f_hits.txt'
        hits_r_fp = '../500r_hits.txt'
        taxa_fp = '/home/BigErn/taxa_data/taxa_mapping.txt'
        
        ''' As above test, but sequence ID EF142855 has a reverse primer 
        weighted score of 0.8, so it should be considered not passing '''
        
        expected_coverage = [{('Archaea', 'Archaea'): [1, 1],
         ('Bacteria', 'Bacteria'): [3, 2]},
         {('Bacteria', 'Firmicutes'): [3, 2],
         ('Archaea', 'Euryarchaeota'): [1, 1]},
         {('Bacteria', 'Firmicutes_Uncultured'): [1, 1],
         ('Bacteria', 'Clostridia'): [2, 1],
         ('Archaea', 'Thermoplasmatales'): [1, 1]}]
        
        # Log file should give base names of all input files, scoring 
        # parameters, and unknown taxa (in this case 0)
        expected_log =\
         ['# Taxonomic coverage report\n# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.\n',
          '# Hits files analyzed: 100f_hits.txt,500r_hits.txt\n', 
          '# Taxonomy mapping file used: taxa_mapping.txt\n', 
          '# Scoring method used: weighted_score\n', 
          '# Score threshold: 0.50\n', 
          '# Total sequences in hits files: 4\n', 
          '# Total sequences with no corresponding taxonomy in taxonomy mapping file: 0\n', 
          '# Sequence IDs lacking taxonomy mapping:\n', 
          '']
        
        actual_coverage, actual_log =\
         get_coverage_data_primer_pair(pair_hits_data, taxa_map, taxa_depth,
         score_type, score_threshold, hits_f_fp, hits_r_fp, taxa_fp)
         
        self.assertEqual(actual_coverage, expected_coverage)
        self.assertEqual(actual_log, expected_log)
        
    def test_get_coverage_data_primer_pair_overall_mismatches(self):
        """ Properly gets sequence coverage, count for primer pair hits data """
        
        pair_hits_data = [self.forward_hits_data, self.reverse_hits_data]
        taxa_map = self.taxa_mapping_data
        taxa_depth = 3
        score_type = 'overall_mismatches'
        score_threshold = 1
        hits_f_fp = '100f_hits.txt'
        hits_r_fp = '../500r_hits.txt'
        taxa_fp = '/home/BigErn/taxa_data/taxa_mapping.txt'
        
        ''' As above test, but sequence ID EF142855 has 2 non three primer  
        mismatches, so it should be considered not passing '''
        
        expected_coverage = [{('Archaea', 'Archaea'): [1, 1],
         ('Bacteria', 'Bacteria'): [3, 2]},
         {('Bacteria', 'Firmicutes'): [3, 2],
         ('Archaea', 'Euryarchaeota'): [1, 1]},
         {('Bacteria', 'Firmicutes_Uncultured'): [1, 1],
         ('Bacteria', 'Clostridia'): [2, 1],
         ('Archaea', 'Thermoplasmatales'): [1, 1]}]
        
        # Log file should give base names of all input files, scoring 
        # parameters, and unknown taxa (in this case 0)
        expected_log =\
         ['# Taxonomic coverage report\n# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.\n',
          '# Hits files analyzed: 100f_hits.txt,500r_hits.txt\n', 
          '# Taxonomy mapping file used: taxa_mapping.txt\n', 
          '# Scoring method used: overall_mismatches\n', 
          '# Score threshold: 1.00\n', 
          '# Total sequences in hits files: 4\n', 
          '# Total sequences with no corresponding taxonomy in taxonomy mapping file: 0\n', 
          '# Sequence IDs lacking taxonomy mapping:\n', 
          '']
        
        actual_coverage, actual_log =\
         get_coverage_data_primer_pair(pair_hits_data, taxa_map, taxa_depth,
         score_type, score_threshold, hits_f_fp, hits_r_fp, taxa_fp)
         
        self.assertEqual(actual_coverage, expected_coverage)
        self.assertEqual(actual_log, expected_log)
        
    def test_get_coverage_data_primer_pair_tp_mismatches(self):
        """ Properly gets sequence coverage, count for primer pair hits data """
        
        pair_hits_data = [self.forward_hits_data, self.reverse_hits_data]
        taxa_map = self.taxa_mapping_data
        taxa_depth = 3
        score_type = 'tp_mismatches'
        score_threshold = 0
        hits_f_fp = '100f_hits.txt'
        hits_r_fp = '../500r_hits.txt'
        taxa_fp = '/home/BigErn/taxa_data/taxa_mapping.txt'
        
        ''' All sequences have zero three prime mismatches, so all should
        pass '''
        
        expected_coverage = [{('Archaea', 'Archaea'): [1, 1],
         ('Bacteria', 'Bacteria'): [3, 3]},
         {('Bacteria', 'Firmicutes'): [3, 3],
         ('Archaea', 'Euryarchaeota'): [1, 1]},
         {('Bacteria', 'Firmicutes_Uncultured'): [1, 1],
         ('Bacteria', 'Clostridia'): [2, 2],
         ('Archaea', 'Thermoplasmatales'): [1, 1]}]
        
        # Log file should give base names of all input files, scoring 
        # parameters, and unknown taxa (in this case 0)
        expected_log =\
         ['# Taxonomic coverage report\n# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.\n',
          '# Hits files analyzed: 100f_hits.txt,500r_hits.txt\n', 
          '# Taxonomy mapping file used: taxa_mapping.txt\n', 
          '# Scoring method used: tp_mismatches\n', 
          '# Score threshold: 0.00\n', 
          '# Total sequences in hits files: 4\n', 
          '# Total sequences with no corresponding taxonomy in taxonomy mapping file: 0\n', 
          '# Sequence IDs lacking taxonomy mapping:\n', 
          '']
        
        actual_coverage, actual_log =\
         get_coverage_data_primer_pair(pair_hits_data, taxa_map, taxa_depth,
         score_type, score_threshold, hits_f_fp, hits_r_fp, taxa_fp)
         
        self.assertEqual(actual_coverage, expected_coverage)
        self.assertEqual(actual_log, expected_log)
        
    def test_get_coverage_data_primer_pair_raises_error(self):
        """ Raises error if forward and reverse sequence IDs do not match """
        
        pair_hits_data = [self.forward_hits_data,
         self.reverse_hits_data_not_matched]
        taxa_map = self.taxa_mapping_data
        taxa_depth = 3
        score_type = 'tp_mismatches'
        score_threshold = 0
        hits_f_fp = '100f_hits.txt'
        hits_r_fp = '../500r_hits.txt'
        taxa_fp = '/home/BigErn/taxa_data/taxa_mapping.txt'
        
        
        self.assertRaises(ValueError, get_coverage_data_primer_pair,
         pair_hits_data, taxa_map, taxa_depth, score_type, score_threshold, 
         hits_f_fp, hits_r_fp, taxa_fp)
    
    def test_get_coverage_data_primer_pair_handles_unassigned(self):
        """ Properly handles cases of sequence IDs not in taxa mapping """
        
        
        
        pair_hits_data = [self.forward_hits_data_taxa_unmapped,
         self.reverse_hits_data_taxa_unmapped]
        taxa_map = self.taxa_mapping_data
        taxa_depth = 3
        score_type = 'weighted_score'
        score_threshold = 1
        hits_f_fp = '100f_hits.txt'
        hits_r_fp = '../500r_hits.txt'
        taxa_fp = '/home/BigErn/taxa_data/taxa_mapping.txt'
        
        '''These are the taxa that correspond to the sequence IDs given
         in the hits data:
        EF142855  Bacteria;Firmicutes;Clostridia
        DQ278118  Bacteria;Firmicutes;Clostridia
        DQ278125 Archaea;Euryarchaeota;Thermoplasmatales;uncultured
        
        XXXX is another sequence in both forward and reverse data, is not in
        the taxa mapping dictionary.
        
        The results in this case should be 1 Archaea domain (and 1 each of its
        descending taxa), 2 bacteria, and 1 Unknown.  On the second level of 
        taxa, both bacteria should be under firmicutes.  The third level should 
        have 2 clostridia.
        
        All of the sequences should be considered passing, as neither of the
        two hits weighted scores are greater than 1, the threshold'''
        
        expected_coverage = [{('Archaea', 'Archaea'): [1, 1],
         ('Unknown', 'Unknown'): [1, 1],
         ('Bacteria', 'Bacteria'): [2, 2]},
         {('Bacteria', 'Firmicutes'): [2, 2],
         ('Archaea', 'Euryarchaeota'): [1, 1]},
         {('Bacteria', 'Clostridia'): [2, 2],
         ('Archaea', 'Thermoplasmatales'): [1, 1]}]
        
        # Log file should give base names of all input files, scoring 
        # parameters, and unknown taxa (in this case 1)
        expected_log =\
         ['# Taxonomic coverage report\n# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.\n',
          '# Hits files analyzed: 100f_hits.txt,500r_hits.txt\n', 
          '# Taxonomy mapping file used: taxa_mapping.txt\n', 
          '# Scoring method used: weighted_score\n', 
          '# Score threshold: 1.00\n', 
          '# Total sequences in hits files: 4\n', 
          '# Total sequences with no corresponding taxonomy in taxonomy mapping file: 1\n', 
          '# Sequence IDs lacking taxonomy mapping:\n', 
          'XXXX\n']
        
        actual_coverage, actual_log =\
         get_coverage_data_primer_pair(pair_hits_data, taxa_map, taxa_depth,
         score_type, score_threshold, hits_f_fp, hits_r_fp, taxa_fp)
         
        self.assertEqual(actual_coverage, expected_coverage)
        self.assertEqual(actual_log, expected_log)


    def test_get_coverage_data_single_primer_weighted_score(self):
        """ Properly gets sequence coverage, count for individual primer """
        
        hits_data = self.reverse_hits_data
        taxa_map = self.taxa_mapping_data
        taxa_depth = 3
        score_type = 'weighted_score'
        score_threshold = 1
        hits_fp = '../500r_hits.txt'
        taxa_fp = '/home/BigErn/taxa_data/taxa_mapping.txt'
        
        '''These are the taxa that correspond to the sequence IDs given
         in the hits data:
        AB302407  Bacteria;Firmicutes;uncultured
        EF142855  Bacteria;Firmicutes;Clostridia
        DQ278118  Bacteria;Firmicutes;Clostridia
        DQ278125 Archaea;Euryarchaeota;Thermoplasmatales;uncultured
        
        The results in this case should be 1 Archaea domain (and 1 each of its
        descending taxa), and 3 bacteria.  On the second level of taxa, all
        3 bacteria should be under firmicutes.  The third level should have
        2 clostridia, and one uncultured (labeled as its parent taxa + 
        _Uncultured to distinguish it from other uncultured taxa at that level).
        
        All of the sequences should be considered passing, as neither of the
        two hits weighted scores are greater than 1, the threshold'''
        
        expected_coverage = [{('Archaea', 'Archaea'): [1, 1],
         ('Bacteria', 'Bacteria'): [3, 3]},
         {('Bacteria', 'Firmicutes'): [3, 3],
         ('Archaea', 'Euryarchaeota'): [1, 1]},
         {('Bacteria', 'Firmicutes_Uncultured'): [1, 1],
         ('Bacteria', 'Clostridia'): [2, 2],
         ('Archaea', 'Thermoplasmatales'): [1, 1]}]
        
        # Log file should give base names of all input files, scoring 
        # parameters, and unknown taxa (in this case 0)
        expected_log =\
         ['# Taxonomic coverage report\n# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.\n',
          '# Hits file analyzed: 500r_hits.txt\n', 
          '# Taxonomy mapping file used: taxa_mapping.txt\n', 
          '# Scoring method used: weighted_score\n', 
          '# Score threshold: 1.00\n', 
          '# Total sequences in hits file: 4\n', 
          '# Total sequences with no corresponding taxonomy in taxonomy mapping file: 0\n', 
          '# Sequence IDs lacking taxonomy mapping:\n', 
          '']
        
        actual_coverage, actual_log =\
         get_coverage_data_single_primer(hits_data, taxa_map, taxa_depth,
         score_type, score_threshold, hits_fp, taxa_fp)
         
        self.assertEqual(actual_coverage, expected_coverage)
        self.assertEqual(actual_log, expected_log)
        
    def test_get_coverage_data_single_primer_weighted_score_inc_score(self):
        """ Properly gets sequence coverage, count for individual primer """
        
        hits_data = self.reverse_hits_data
        taxa_map = self.taxa_mapping_data
        taxa_depth = 3
        score_type = 'weighted_score'
        score_threshold = 0.5
        hits_fp = '../500r_hits.txt'
        taxa_fp = '/home/BigErn/taxa_data/taxa_mapping.txt'
        
        '''As above, but with increased stringency for weighted score 
        threshold.
        
        EF142855 has a weighted score of 0.8, so should not pass'''
        
        expected_coverage = [{('Archaea', 'Archaea'): [1, 1],
         ('Bacteria', 'Bacteria'): [3, 2]},
         {('Bacteria', 'Firmicutes'): [3, 2],
         ('Archaea', 'Euryarchaeota'): [1, 1]},
         {('Bacteria', 'Firmicutes_Uncultured'): [1, 1],
         ('Bacteria', 'Clostridia'): [2, 1],
         ('Archaea', 'Thermoplasmatales'): [1, 1]}]
        
        # Log file should give base names of all input files, scoring 
        # parameters, and unknown taxa (in this case 0)
        expected_log =\
         ['# Taxonomic coverage report\n# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.\n',
          '# Hits file analyzed: 500r_hits.txt\n', 
          '# Taxonomy mapping file used: taxa_mapping.txt\n', 
          '# Scoring method used: weighted_score\n', 
          '# Score threshold: 0.50\n', 
          '# Total sequences in hits file: 4\n', 
          '# Total sequences with no corresponding taxonomy in taxonomy mapping file: 0\n', 
          '# Sequence IDs lacking taxonomy mapping:\n', 
          '']
        
        actual_coverage, actual_log =\
         get_coverage_data_single_primer(hits_data, taxa_map, taxa_depth,
         score_type, score_threshold, hits_fp, taxa_fp)
         
        self.assertEqual(actual_coverage, expected_coverage)
        self.assertEqual(actual_log, expected_log)
        
    def test_get_coverage_data_single_primer_overall_mismatches(self):
        """ Properly gets sequence coverage, count for individual primer """
        
        hits_data = self.reverse_hits_data
        taxa_map = self.taxa_mapping_data
        taxa_depth = 3
        score_type = 'overall_mismatches'
        score_threshold = 1.0
        hits_fp = '../500r_hits.txt'
        taxa_fp = '/home/BigErn/taxa_data/taxa_mapping.txt'
        
        '''As above, but with increased stringency for weighted score 
        threshold.
        
        EF142855 has 2 non three primer mismatches, so should not pass'''
        
        expected_coverage = [{('Archaea', 'Archaea'): [1, 1],
         ('Bacteria', 'Bacteria'): [3, 2]},
         {('Bacteria', 'Firmicutes'): [3, 2],
         ('Archaea', 'Euryarchaeota'): [1, 1]},
         {('Bacteria', 'Firmicutes_Uncultured'): [1, 1],
         ('Bacteria', 'Clostridia'): [2, 1],
         ('Archaea', 'Thermoplasmatales'): [1, 1]}]
        
        # Log file should give base names of all input files, scoring 
        # parameters, and unknown taxa (in this case 0)
        expected_log =\
         ['# Taxonomic coverage report\n# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.\n',
          '# Hits file analyzed: 500r_hits.txt\n', 
          '# Taxonomy mapping file used: taxa_mapping.txt\n', 
          '# Scoring method used: overall_mismatches\n', 
          '# Score threshold: 1.00\n', 
          '# Total sequences in hits file: 4\n', 
          '# Total sequences with no corresponding taxonomy in taxonomy mapping file: 0\n', 
          '# Sequence IDs lacking taxonomy mapping:\n', 
          '']
        
        actual_coverage, actual_log =\
         get_coverage_data_single_primer(hits_data, taxa_map, taxa_depth,
         score_type, score_threshold, hits_fp, taxa_fp)
         
        self.assertEqual(actual_coverage, expected_coverage)
        self.assertEqual(actual_log, expected_log)
        
    def test_get_coverage_data_single_primer_tp_mismatches(self):
        """ Properly gets sequence coverage, count for individual primer """
        
        hits_data = self.reverse_hits_data_tp_mm
        taxa_map = self.taxa_mapping_data
        taxa_depth = 3
        score_type = 'tp_mismatches'
        score_threshold = 1.0
        hits_fp = '../500r_hits.txt'
        taxa_fp = '/home/BigErn/taxa_data/taxa_mapping.txt'
        
        '''As above, but with increased stringency for weighted score 
        threshold.
        
        DQ278118, DQ278125 both have more than 1 three prime mismatch, so 
        should not be counted as passing'''
        
        expected_coverage = [{('Archaea', 'Archaea'): [1, 0],
         ('Bacteria', 'Bacteria'): [3, 2]},
         {('Bacteria', 'Firmicutes'): [3, 2],
         ('Archaea', 'Euryarchaeota'): [1, 0]},
         {('Bacteria', 'Firmicutes_Uncultured'): [1, 1],
         ('Bacteria', 'Clostridia'): [2, 1],
         ('Archaea', 'Thermoplasmatales'): [1, 0]}]
        
        # Log file should give base names of all input files, scoring 
        # parameters, and unknown taxa (in this case 0)
        expected_log =\
         ['# Taxonomic coverage report\n# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.\n',
          '# Hits file analyzed: 500r_hits.txt\n', 
          '# Taxonomy mapping file used: taxa_mapping.txt\n', 
          '# Scoring method used: tp_mismatches\n', 
          '# Score threshold: 1.00\n', 
          '# Total sequences in hits file: 4\n', 
          '# Total sequences with no corresponding taxonomy in taxonomy mapping file: 0\n', 
          '# Sequence IDs lacking taxonomy mapping:\n', 
          '']
        
        actual_coverage, actual_log =\
         get_coverage_data_single_primer(hits_data, taxa_map, taxa_depth,
         score_type, score_threshold, hits_fp, taxa_fp)
         
        self.assertEqual(actual_coverage, expected_coverage)
        self.assertEqual(actual_log, expected_log)
        
        
    def test_get_coverage_data_single_primer_missing_taxa(self):
        """ Properly gets sequence coverage, count for individual primer """
        
        hits_data = self.reverse_hits_data_taxa_unmapped
        taxa_map = self.taxa_mapping_data
        taxa_depth = 3
        score_type = 'weighted_score'
        score_threshold = 1
        hits_fp = '../500r_hits.txt'
        taxa_fp = '/home/BigErn/taxa_data/taxa_mapping.txt'
        
        '''Should get one sequence, XXXX as Unknown'''
        
        expected_coverage = [{('Archaea', 'Archaea'): [1, 1],
         ('Unknown', 'Unknown'): [1, 1],
         ('Bacteria', 'Bacteria'): [2, 2]},
         {('Bacteria', 'Firmicutes'): [2, 2],
         ('Archaea', 'Euryarchaeota'): [1, 1]},
         {('Bacteria', 'Clostridia'): [2, 2],
         ('Archaea', 'Thermoplasmatales'): [1, 1]}]
        
        # Log file should give base names of all input files, scoring 
        # parameters, and unknown taxa (in this case 0)
        expected_log =\
         ['# Taxonomic coverage report\n# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.\n',
          '# Hits file analyzed: 500r_hits.txt\n', 
          '# Taxonomy mapping file used: taxa_mapping.txt\n', 
          '# Scoring method used: weighted_score\n', 
          '# Score threshold: 1.00\n', 
          '# Total sequences in hits file: 4\n', 
          '# Total sequences with no corresponding taxonomy in taxonomy mapping file: 1\n', 
          '# Sequence IDs lacking taxonomy mapping:\n', 
          'XXXX\n']
        
        actual_coverage, actual_log =\
         get_coverage_data_single_primer(hits_data, taxa_map, taxa_depth,
         score_type, score_threshold, hits_fp, taxa_fp)
         
        self.assertEqual(actual_coverage, expected_coverage)
        self.assertEqual(actual_log, expected_log)
        
        
    def test_get_coverage_filepath(self):
        """ Properly generates filepath for text output file """
        
        # Test single primer output directory name generation
        curr_primer_output_dir = '/primer_test/200f/'
        hits = '../200f_hits.txt'
        
        expected_outf = '/primer_test/200f/200f_hits_coverage.txt'
        
        actual_outf = get_coverage_filepath(curr_primer_output_dir, hits)
         
        self.assertEqual(actual_outf, expected_outf)
        
        # Test primer pair
        hits = ['../200f_hits.txt', 'reverse_primer/500r_hits.txt']
        
        expected_outf =\
         '/primer_test/200f/200f_hits_500r_hits_coverage.txt'
         
        actual_outf = get_coverage_filepath(curr_primer_output_dir, hits,
         primer_pair = True)
         
        self.assertEqual(actual_outf, expected_outf)
        
        
    def test_get_coverage_graph_root_filepath(self):
        """ Properly gets output root filepath, title for graph at taxa level"""
        
        # Test single primer result
        curr_primer_output_dir = '/primer_test/200f/'
        hits = '../200f_hits.txt'
        
        expected_graph_fp = '/primer_test/200f/200f_hits_coverage'
        expected_graph_title = "Predicted Taxonomic Coverage\n200f_hits\n"

        actual_graph_fp, actual_graph_title =\
         get_coverage_graph_root_filepath(curr_primer_output_dir, hits)
         
        self.assertEqual(actual_graph_fp, expected_graph_fp)
        self.assertEqual(actual_graph_title, expected_graph_title)
        
        # Test primer pair
        curr_primer_output_dir = '/primer_test/200f/'
        hits = ['../200f_hits.txt', 'reverse_primer/500r_hits.txt']
        
        expected_graph_fp = '/primer_test/200f/200f_hits_500r_hits_coverage'
        expected_graph_title =\
         "Predicted Taxonomic Coverage\n200f_hits_500r_hits\n"

        actual_graph_fp, actual_graph_title =\
         get_coverage_graph_root_filepath(curr_primer_output_dir, hits,
         primer_pair = True)
         
        self.assertEqual(actual_graph_fp, expected_graph_fp)
        self.assertEqual(actual_graph_title, expected_graph_title)
        
        
    def test_get_graph_fp(self):
        """ Properly generates graph filepath for given taxa level """
        
        # Test at taxa level 0
        primer_graph_filepath = 'graphs/200f_hits'
        curr_taxa_level = 0
        
        expected_graph_fp = "graphs/200f_hits_taxonomy_level_0.pdf"
        
        actual_graph_fp = get_graph_fp(primer_graph_filepath, curr_taxa_level)
        
        self.assertEqual(actual_graph_fp, expected_graph_fp)
        
        # Test for lower taxa level, domain name
        rimer_graph_filepath = 'graphs/200f_hits'
        curr_taxa_level = 2
        domain = "Eukarya"
        
        expected_graph_fp =\
         "graphs/200f_hits_Eukarya_taxonomy_level_2.pdf"
        
        actual_graph_fp = get_graph_fp(primer_graph_filepath, curr_taxa_level,
         domain)
        
        self.assertEqual(actual_graph_fp, expected_graph_fp)
        
        
        
    def test_get_log_filepath(self):
        """ Properly generates log filepath """
        
                     
        # Test single primer
        curr_primer_output_dir = 'primer_test/'
        hits = '200f_hits.txt'
        
        expected_log_fp = 'primer_test/200f_hits.log'
        
        actual_log_fp = get_log_filepath(curr_primer_output_dir, hits)
        
        self.assertEqual(actual_log_fp, expected_log_fp)
        
        # Test primer pair
        curr_primer_output_dir = 'primer_test/'
        hits = ['200f_hits.txt', '500r_hits.txt']
        
        expected_log_fp = 'primer_test/200f_hits_500r_hits.log'
        
        actual_log_fp = get_log_filepath(curr_primer_output_dir, hits,
         primer_pair = True)
         
        self.assertEqual(actual_log_fp, expected_log_fp)
        
        
        
    def test_get_output_dir(self):
        """ Properly generates output dir"""
        

                   
        # Test output dir for single primer
        output_dir = "."
        hits = "../200f_hits.txt"
        
        expected_output_dir = "./200f_hits_primer_coverage/"
        
        actual_output_dir = get_output_dir(output_dir, hits)
        
        self.assertEqual(actual_output_dir, expected_output_dir)
        
        # Test output dir for primer pair
        output_dir = "/home/BigErn/primers"
        hits = ["../200f_hits.txt", "500r_hits.txt"]
        
        expected_output_dir =\
         "/home/BigErn/primers/200f_hits_500r_hits_primers_coverage/"
        
        actual_output_dir = get_output_dir(output_dir, hits, primer_pair = True)
        
        self.assertEqual(actual_output_dir, expected_output_dir)
        
        
    def test_get_primer_pairs(self):
        """ Properly finds all possible forward, reverse primer pairs """
        
        hits_f = ['../200f_hits.txt', '500r_hits.txt', 'primers/1000r_hits.txt']
        
        expected_pairs = [('../200f_hits.txt', '500r_hits.txt'),
                          ('../200f_hits.txt', 'primers/1000r_hits.txt')]
                          
        actual_pairs = get_primer_pairs(hits_f)
        
        self.assertEqual(actual_pairs, expected_pairs)
        
        # Should get empty list with no foward/reverse primer pairs present
        
        hits_f = ['500r_hits.txt', 'primers/1000r_hits.txt']
        
        expected_pairs = []
                          
        actual_pairs = get_primer_pairs(hits_f)
        
        self.assertEqual(actual_pairs, expected_pairs)
        
    def test_graph_taxa_coverage(self):
        """ Main program loop, check text output, creation of graphs """
        
        hits_fps = self.output_dir
        taxa_fp = self.taxa_mapping
        taxa_depth = 2
        primer_pairs = True
        all_files = True
        output_dir = self.output_dir
        score_type = "weighted_score"
        score_threshold = 0.5
        
        graph_taxa_coverage(hits_fps, taxa_fp, taxa_depth, primer_pairs,
         all_files, output_dir, score_type, score_threshold)
         
        # Many files generated, create a list to iterate through to test
        # for creation of all files, see end of this module.
        expected_pdf_fps = self.pdf_filepaths
        
        # Correction for output filepaths due to randomly generated output dir
        for fp_index in range(len(expected_pdf_fps)):
            expected_pdf_fps[fp_index] = output_dir + expected_pdf_fps[fp_index]
            
        
        for curr_fp in expected_pdf_fps:
            if not isfile(curr_fp):
                raise ValueError('Expected pdf %s not found.' % curr_fp)
                
        # Test for proper contents of log, coverage files
        f_log_fp =\
         output_dir + '100f_hits_primer_coverage/100f_hits.log'
        
        f_log_f = open(f_log_fp, "U")
        
        f_log_results = [line.strip() for line in f_log_f]
        
        self.assertEqual(f_log_results, self.expected_100f_log)
        
        f_coverage_fp =\
         output_dir + '100f_hits_primer_coverage/100f_hits_coverage.txt'
         
        f_coverage_f = open(f_coverage_fp, "U")
        
        f_coverage_results = [line.strip() for line in f_coverage_f]
        
        self.assertEqual(f_coverage_results, self.expected_100f_coverage)
        
        # Test reverse primer results
        r_log_fp =\
         output_dir + '500r_hits_primer_coverage/500r_hits.log'
        
        r_log_f = open(r_log_fp, "U")
        
        r_log_results = [line.strip() for line in r_log_f]
        
        self.assertEqual(r_log_results, self.expected_500r_log)
        
        r_coverage_fp =\
         output_dir + '500r_hits_primer_coverage/500r_hits_coverage.txt'
         
        r_coverage_f = open(r_coverage_fp, "U")
        
        r_coverage_results = [line.strip() for line in r_coverage_f]
        
        self.assertEqual(r_coverage_results, self.expected_500r_coverage)
        
        # Test combined results (should be equivalent to reverse primer,
        # which had worse scores
        c_log_fp =\
         output_dir + '100f_hits_500r_hits_primers_coverage/100f_hits_500r_hits.log'
        
        c_log_f = open(c_log_fp, "U")
        
        c_log_results = [line.strip() for line in c_log_f]
        
        self.assertEqual(c_log_results, self.expected_combined_log)
        
        c_coverage_fp =\
         output_dir + '100f_hits_500r_hits_primers_coverage/100f_hits_500r_hits_coverage.txt'
         
        c_coverage_f = open(c_coverage_fp, "U")
        
        c_coverage_results = [line.strip() for line in c_coverage_f]
        
        self.assertEqual(c_coverage_results, self.expected_combined_coverage)

        
        
    def test_write_log_file(self):
        """ Writes log files correctly """

                   
        log_filepath = self.output_dir + "/test.log"
        log_data = self.expected_500r_log
        
        write_log_file(log_filepath, log_data)
        
        actual_log_f = open(log_filepath, "U")
        
        actual_log = [line for line in actual_log_f]
        
        self.assertEqual(actual_log, ["".join(self.expected_500r_log)])
        
    def test_write_primer_coverage(self):
        """ Writes text file of taxonomic coverage correctly """
        
        coverage_fp = self.output_dir + "/test_coverage.txt"
        
        coverage_data = [{('Archaea', 'Archaea'): [1, 1],
         ('Bacteria', 'Bacteria'): [3, 2]},
         {('Bacteria', 'Firmicutes'): [2, 1],
         ('Archaea', 'Euryarchaeota'): [1, 1]},
         {('Bacteria', 'Firmicutes_Uncultured'): [1, 1],
         ('Bacteria', 'Clostridia'): [2, 2],
         ('Archaea', 'Thermoplasmatales'): [1, 1]}]
         
        expected_coverage_file = ['# This file is written in descending levels of taxonomic depth',
         '# Each line beings with the taxonomy level, followed by the first level of', 
         '# taxonomy for a given sequence, generally the domain, followed the taxonomy', 
         '# for the current level.', 
         '# The total sequences for a given taxonomy are listed first, followed by the', 
         '# percentage that are lower to or equal to the threshold score for passing.', 
         '# taxonomy level, first level taxonomy, taxonomy classification for current level, total seqs for given classification, percent seqs passing score', 
         '0,Archaea,Archaea,1,1.0000', 
         '0,Bacteria,Bacteria,3,0.6667', 
         '1,Archaea,Euryarchaeota,1,1.0000', 
         '1,Bacteria,Firmicutes,2,0.5000',
         '2,Archaea,Thermoplasmatales,1,1.0000',
         '2,Bacteria,Firmicutes_Uncultured,1,1.0000',
         '2,Bacteria,Clostridia,2,1.0000']
         
        write_primer_coverage(coverage_fp, coverage_data)
        
        actual_coverage_f = open(coverage_fp, "U")
        
        actual_coverage = [line.strip() for line in actual_coverage_f]
        
        self.assertEqual(actual_coverage, expected_coverage_file) 
        
    def test_write_primer_coverage_graph(self):
        """ Generates output graph files, no data tested """
        
        primer_graph_fp = self.output_dir + '/graph'
        
        graph_title = "unit_tests"
        
        coverage_data = [{('Archaea', 'Archaea'): [1, 1],
         ('Bacteria', 'Bacteria'): [3, 2]},
         {('Bacteria', 'Firmicutes'): [2, 1],
         ('Archaea', 'Euryarchaeota'): [1, 1]},
         {('Bacteria', 'Firmicutes_Uncultured'): [1, 1],
         ('Bacteria', 'Clostridia'): [2, 2],
         ('Archaea', 'Thermoplasmatales'): [1, 1]}]
         
        write_primer_coverage_graph(primer_graph_fp, graph_title, coverage_data)
        
        # Can only test for creation of output graphs at correct filepaths
        
        expected_graph_fps = [
        self.output_dir + 'graph_taxonomy_level_0.pdf',
        self.output_dir + 'graph_Archaea_taxonomy_level_1.pdf',
        self.output_dir + 'graph_Bacteria_taxonomy_level_1.pdf',
        self.output_dir + 'graph_Archaea_taxonomy_level_2.pdf',
        self.output_dir + 'graph_Bacteria_taxonomy_level_2.pdf']
        
        for graph in expected_graph_fps:
            if not isfile(graph):
                raise ValueError('Expected graph %s not created.' % graph)
                
    def test_get_curr_taxa(self):
        """ Returns taxa key correctly, handles 'uncultured' taxa """
        
        # Test without any uncultured groups
        
        taxa = ['Bacteria', 'Firmicutes', 'Clostridria', 'Clostridium']
        taxa_level = 2
        domain_index = 0
        
        expected_domain_key = ('Bacteria', 'Clostridria')
        
        actual_domain_key = get_curr_taxa(taxa, taxa_level, domain_index)
        
        self.assertEqual(actual_domain_key, expected_domain_key)
        
        # Should return None if taxa_level is greater than that provided in
        # taxa
        
        taxa_level = 5
        
        expected_domain_key = None
        
        actual_domain_key = get_curr_taxa(taxa, taxa_level, domain_index)
        
        self.assertEqual(actual_domain_key, expected_domain_key)
        
        # Should give name of higher level taxa followed by _Uncultured
        # if current level given as "uncultured"
        
        taxa = ['Bacteria', 'Firmicutes', 'Uncultured', 'Uncultured']
        taxa_level = 2
        domain_index = 0
        
        expected_domain_key = ('Bacteria', 'Firmicutes_Uncultured')
        
        actual_domain_key = get_curr_taxa(taxa, taxa_level, domain_index)
        
        self.assertEqual(actual_domain_key, expected_domain_key)
        
        # If 'uncultured' already in name of higher level taxa, should just
        # take name of higher level taxa instead
        taxa = ['Bacteria', 'Uncultured firmicutes', 'Uncultured', 'Uncultured']
        taxa_level = 2
        domain_index = 0
        
        expected_domain_key = ('Bacteria', 'Uncultured firmicutes')
        
        actual_domain_key = get_curr_taxa(taxa, taxa_level, domain_index)
        
        self.assertEqual(actual_domain_key, expected_domain_key)
        
# Large data sets at the end of the file for better readability

forward_hits_file = """# Primer: 515f 5'-GTGCCAGCMGCCGCGGTAA-3'
# Input fasta file: combined_v15.fasta
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
AB302407,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,1202,1,0,False,0,0,0.4,False
EF142855,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,1088,1,0,False,0,0,0.4,False
DQ278118,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,350,1,0,False,0,0,0.4,False
DQ278125,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,350,1,0,False,0,0,0.4,False"""

forward_hits_data = """AB302407,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,1202,1,0,False,0,0,0.4,False
EF142855,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,1088,1,0,False,0,0,0.4,False
DQ278118,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,350,1,0,False,0,0,0.4,False
DQ278125,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,350,1,0,False,0,0,0.4,False""".split('\n')

forward_hits_data_taxa_unmapped = """XXXX,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,1202,1,0,False,0,0,0.4,False
EF142855,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,1088,1,0,False,0,0,0.4,False
DQ278118,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,350,1,0,False,0,0,0.4,False
DQ278125,GTGTCAGCCGCCGCGGTAA,GTGCCAGCMGCCGCGGTAA,350,1,0,False,0,0,0.4,False""".split('\n')

reverse_hits_file = """# Primer: 806r 5'-GGACTACVSGGGTATCTAAT-3'
# Input fasta file: combined_v15.fasta
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
AB302407,ATTAGATACCCCTGTAGTCC,ATTAGATACCCSBGTAGTCC,2238,0,0,False,0,0,0.0,False
EF142855,ATTAGATACCCTAGTAGTCC,ATTAGATACCCSBGTAGTCC,1360,2,0,False,0,0,0.8,False
DQ278118,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,622,0,0,False,0,0,0.0,False
DQ278125,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,622,0,0,False,0,0,0.0,False"""

reverse_hits_data = """AB302407,ATTAGATACCCCTGTAGTCC,ATTAGATACCCSBGTAGTCC,2238,0,0,False,0,0,0.0,False
EF142855,ATTAGATACCCTAGTAGTCC,ATTAGATACCCSBGTAGTCC,1360,2,0,False,0,0,0.8,False
DQ278118,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,622,0,0,False,0,0,0.0,False
DQ278125,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,622,0,0,False,0,0,0.0,False""".split('\n')

reverse_hits_data_tp_mm = """AB302407,ATTAGATACCCCTGTAGTCC,ATTAGATACCCSBGTAGTCC,2238,0,0,False,0,0,0.0,False
EF142855,ATTAGATACCCTAGTAGTCC,ATTAGATACCCSBGTAGTCC,1360,2,0,False,0,0,0.8,False
DQ278118,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,622,0,2,False,0,0,0.0,False
DQ278125,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,622,0,3,False,0,0,0.0,False""".split('\n')


reverse_hits_data_not_matched = """XXXXXXX,ATTAGATACCCCTGTAGTCC,ATTAGATACCCSBGTAGTCC,2238,0,0,False,0,0,0.0,False
EF142855,ATTAGATACCCTAGTAGTCC,ATTAGATACCCSBGTAGTCC,1360,2,0,False,0,0,0.8,False
DQ278118,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,622,0,0,False,0,0,0.0,False
DQ278125,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,622,0,0,False,0,0,0.0,False""".split('\n')

reverse_hits_data_taxa_unmapped = """XXXX,ATTAGATACCCCTGTAGTCC,ATTAGATACCCSBGTAGTCC,2238,0,0,False,0,0,0.0,False
EF142855,ATTAGATACCCTAGTAGTCC,ATTAGATACCCSBGTAGTCC,1360,2,0,False,0,0,0.8,False
DQ278118,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,622,0,0,False,0,0,0.0,False
DQ278125,ATTAGATACCCGGGTAGTCC,ATTAGATACCCSBGTAGTCC,622,0,0,False,0,0,0.0,False""".split('\n')

taxa_mapping_file = """AB329763\tArchaea;Crenarchaeota;uncultured;uncultured
AB302407\tBacteria;Firmicutes;uncultured
DQ243733\tArchaea;Crenarchaeota;uncultured;uncultured
EF142855\tBacteria;Firmicutes;Clostridia
DQ278118\tBacteria;Firmicutes;Clostridia
DQ924802\tArchaea;Crenarchaeota;uncultured
EF020872\tArchaea;Crenarchaeota;uncultured;uncultured
EF021345\tArchaea;Crenarchaeota;uncultured;uncultured
EF021657\tArchaea;Crenarchaeota;uncultured;uncultured
EF022107\tArchaea;Crenarchaeota;uncultured;uncultured
DQ278125\tArchaea;Euryarchaeota;Thermoplasmatales;uncultured
EF203630\tArchaea;Crenarchaeota;uncultured;uncultured
FJ404042\tArchaea;Crenarchaeota;uncultured;uncultured
AJ428044\tArchaea;Crenarchaeota;uncultured;uncultured
AY367317\tArchaea;Euryarchaeota;Methanosarcinales;Methanosaeta
AY800210\tArchaea;Euryarchaeota;Halobacteriales;uncultured
DQ186711\tArchaea;Euryarchaeota;Halobacteriales;Halobacteriales;Haloarcula
DQ397381\tArchaea;Euryarchaeota;Halobacteriales;uncultured"""

taxa_mapping_data = {'EF020872': 'Archaea;Crenarchaeota;uncultured;uncultured', 'EF203630': 'Archaea;Crenarchaeota;uncultured;uncultured', 'EF021345': 'Archaea;Crenarchaeota;uncultured;uncultured', 'AB329763': 'Archaea;Crenarchaeota;uncultured;uncultured', 'DQ924802': 'Archaea;Crenarchaeota;uncultured', 'DQ278118': 'Bacteria;Firmicutes;Clostridia', 'AJ428044': 'Archaea;Crenarchaeota;uncultured;uncultured', 'AY800210': 'Archaea;Euryarchaeota;Halobacteriales;uncultured', 'AY367317': 'Archaea;Euryarchaeota;Methanosarcinales;Methanosaeta', 'EF021657': 'Archaea;Crenarchaeota;uncultured;uncultured', 'DQ278125': 'Archaea;Euryarchaeota;Thermoplasmatales;uncultured', 'DQ186711': 'Archaea;Euryarchaeota;Halobacteriales;Halobacteriales;Haloarcula', 'DQ397381': 'Archaea;Euryarchaeota;Halobacteriales;uncultured', 'FJ404042': 'Archaea;Crenarchaeota;uncultured;uncultured', 'DQ243733': 'Archaea;Crenarchaeota;uncultured;uncultured', 'AB302407': 'Bacteria;Firmicutes;uncultured', 'EF142855': 'Bacteria;Firmicutes;Clostridia', 'EF022107': 'Archaea;Crenarchaeota;uncultured;uncultured'}

pdf_filepaths = [
'100f_hits_primer_coverage/100f_hits_coverage_Archaea_taxonomy_level_1.pdf',
'100f_hits_primer_coverage/100f_hits_coverage_Bacteria_taxonomy_level_1.pdf',
'100f_hits_primer_coverage/100f_hits_coverage_taxonomy_level_0.pdf',
'500r_hits_primer_coverage/500r_hits_coverage_Archaea_taxonomy_level_1.pdf',
'500r_hits_primer_coverage/500r_hits_coverage_Bacteria_taxonomy_level_1.pdf',
'500r_hits_primer_coverage/500r_hits_coverage_taxonomy_level_0.pdf',
'100f_hits_500r_hits_primers_coverage/100f_hits_500r_hits_coverage_Archaea_taxonomy_level_1.pdf',
'100f_hits_500r_hits_primers_coverage/100f_hits_500r_hits_coverage_Bacteria_taxonomy_level_1.pdf',                       
'100f_hits_500r_hits_primers_coverage/100f_hits_500r_hits_coverage_taxonomy_level_0.pdf']                       

expected_100f_log = """# Taxonomic coverage report
# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.
# Hits file analyzed: 100f_hits.txt
# Taxonomy mapping file used: taxa_mapping.txt
# Scoring method used: weighted_score
# Score threshold: 0.50
# Total sequences in hits file: 4
# Total sequences with no corresponding taxonomy in taxonomy mapping file: 0
# Sequence IDs lacking taxonomy mapping:""".split('\n')

expected_100f_coverage = """# This file is written in descending levels of taxonomic depth
# Each line beings with the taxonomy level, followed by the first level of
# taxonomy for a given sequence, generally the domain, followed the taxonomy
# for the current level.
# The total sequences for a given taxonomy are listed first, followed by the
# percentage that are lower to or equal to the threshold score for passing.
# taxonomy level, first level taxonomy, taxonomy classification for current level, total seqs for given classification, percent seqs passing score
0,Archaea,Archaea,1,1.0000
0,Bacteria,Bacteria,3,1.0000
1,Archaea,Euryarchaeota,1,1.0000
1,Bacteria,Firmicutes,3,1.0000""".split('\n')

expected_500r_log = """# Taxonomic coverage report
# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.
# Hits file analyzed: 500r_hits.txt
# Taxonomy mapping file used: taxa_mapping.txt
# Scoring method used: weighted_score
# Score threshold: 0.50
# Total sequences in hits file: 4
# Total sequences with no corresponding taxonomy in taxonomy mapping file: 0
# Sequence IDs lacking taxonomy mapping:""".split('\n')

expected_500r_coverage = """# This file is written in descending levels of taxonomic depth
# Each line beings with the taxonomy level, followed by the first level of
# taxonomy for a given sequence, generally the domain, followed the taxonomy
# for the current level.
# The total sequences for a given taxonomy are listed first, followed by the
# percentage that are lower to or equal to the threshold score for passing.
# taxonomy level, first level taxonomy, taxonomy classification for current level, total seqs for given classification, percent seqs passing score
0,Archaea,Archaea,1,1.0000
0,Bacteria,Bacteria,3,0.6667
1,Archaea,Euryarchaeota,1,1.0000
1,Bacteria,Firmicutes,3,0.6667""".split('\n')

expected_combined_log = """# Taxonomic coverage report
# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.
# Hits files analyzed: 100f_hits.txt,500r_hits.txt
# Taxonomy mapping file used: taxa_mapping.txt
# Scoring method used: weighted_score
# Score threshold: 0.50
# Total sequences in hits files: 4
# Total sequences with no corresponding taxonomy in taxonomy mapping file: 0
# Sequence IDs lacking taxonomy mapping:""".split('\n')

expected_combined_coverage = """# This file is written in descending levels of taxonomic depth
# Each line beings with the taxonomy level, followed by the first level of
# taxonomy for a given sequence, generally the domain, followed the taxonomy
# for the current level.
# The total sequences for a given taxonomy are listed first, followed by the
# percentage that are lower to or equal to the threshold score for passing.
# taxonomy level, first level taxonomy, taxonomy classification for current level, total seqs for given classification, percent seqs passing score
0,Archaea,Archaea,1,1.0000
0,Bacteria,Bacteria,3,0.6667
1,Archaea,Euryarchaeota,1,1.0000
1,Bacteria,Firmicutes,3,0.6667""".split('\n')



if __name__ == "__main__":
    main()
