#!/usr/bin/env python

from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

from os.path import basename, isdir, isfile, exists
from os import remove, makedirs
from shutil import rmtree


from numpy import array, nan, isnan
from cogent3.util.unit_test import TestCase, main
from cogent3 import DNA
from cogent3.parse.fasta import MinimalFastaParser
from primerprospector.old_cogent import get_tmp_filename
from primerprospector.cogentutil.misc import remove_files, get_random_directory_name




from primerprospector.analyze_primers import match_scorer_ambigs,\
 get_outf_paths, primer_hit_histogram, write_primer_histogram,\
 primer_to_match_query, pair_hmm_align_unaligned_seqs, local_align_primer_seq,\
 score_primer, hits_seq_end, get_yaxis_max, get_primers, get_hits_data,\
 generate_hits_file_and_histogram, analyze_primers



 
class AnalyzePrimersTests(TestCase):
    
    def setUp(self):
        """ Initialize tmp files for unit tests """
        
        
        self.sample_fasta_file = get_tmp_filename(prefix = "sample_seqs_",
         suffix = ".fasta")
        seq_file = open(self.sample_fasta_file, 'w')
        seq_file.write(sample_bacterial_seqs)
        seq_file.close()
        
        
        self.output_dir = get_random_directory_name(prefix = '/tmp/')
        self.output_dir += '/'
         
        self.expected_bacterial_seqs_hits_file_27f =\
         expected_bacterial_seqs_hits_file_27f 
        self.expected_bacterial_seqs_hits_file2 =\
         expected_bacterial_seqs_hits_file2
         
        self.sample_bacterial_seqs_expected_hits =\
         sample_bacterial_seqs_expected_hits
        self.sample_bacterial_seqs_expected_hist =\
         sample_bacterial_seqs_expected_hist
         
        self.archaeal_16s_seq1 = archaeal_16s_seq1
        
         
        if not exists(self.output_dir):
            makedirs(self.output_dir)
        
        self._files_to_remove =\
         [self.sample_fasta_file]
         
    def tearDown(self):
        remove_files(self._files_to_remove)
        rmtree(self.output_dir)

   
    
    def test_primer_to_match_query(self):
        """primers correctly converted to matchable sequence
        """
        f_primer = DNA.makeSequence('TCGGA',Name='primerf_test')
        expected = 'TCGGA'
        self.assertEqual(primer_to_match_query(f_primer),expected)
        
        r_primer = DNA.makeSequence('ATCGGCC',Name='273r')
        expected = 'GGCCGAT'
        self.assertEqual(primer_to_match_query(r_primer),expected)
        
        invalid_primer = DNA.makeSequence('TCCA',Name='Invalid_primer_name')
        self.assertRaises(ValueError,primer_to_match_query,invalid_primer)
        
        
    def test_match_scorer_ambigs(self):
        """ match_scorer_ambigs identifies matches/mismatches w/ ambig chars
        """
        # random selection of match/mismatch states
        
        sw_scorer=match_scorer_ambigs(1, -1)
        
        # matches
        self.assertEqual(sw_scorer('A', 'A'), 1)
        self.assertEqual(sw_scorer('A', 'R'), 1)
        self.assertEqual(sw_scorer('R', 'A'), 1)
        self.assertEqual(sw_scorer('W', 'A'), 1)
        self.assertEqual(sw_scorer('A', 'W'), 1)
        self.assertEqual(sw_scorer('M', 'A'), 1)
        self.assertEqual(sw_scorer('D', 'A'), 1)
        self.assertEqual(sw_scorer('H', 'A'), 1)
        self.assertEqual(sw_scorer('V', 'A'), 1)
        self.assertEqual(sw_scorer('A', 'V'), 1)
        self.assertEqual(sw_scorer('C', 'C'), 1)
        self.assertEqual(sw_scorer('Y', 'C'), 1)
        self.assertEqual(sw_scorer('C', 'Y'), 1)
        self.assertEqual(sw_scorer('C', 'C'), 1)
        self.assertEqual(sw_scorer('B', 'C'), 1)
        self.assertEqual(sw_scorer('T', 'T'), 1)
        self.assertEqual(sw_scorer('G', 'G'), 1)
        self.assertEqual(sw_scorer('-', '-'), 1)
        
        # mismatches
        self.assertEqual(sw_scorer('A', 'T'), -1)
        self.assertEqual(sw_scorer('T', 'A'), -1)
        self.assertEqual(sw_scorer('T', 'G'), -1)
        self.assertEqual(sw_scorer('T', 'C'), -1)
        self.assertEqual(sw_scorer('W', 'C'), -1)
        self.assertEqual(sw_scorer('D', 'C'), -1)
        self.assertEqual(sw_scorer('Y', 'A'), -1)
        self.assertEqual(sw_scorer('R', 'T'), -1)
        self.assertEqual(sw_scorer('-', 'T'), -1)
        self.assertEqual(sw_scorer('T', '-'), -1)
        
        # unknown chars result in ValueError
        self.assertRaises(ValueError, sw_scorer, 'x', 'x')
        self.assertRaises(ValueError, sw_scorer, 'x', 'T')
        self.assertRaises(ValueError, sw_scorer, 'T', 'x')
        
    def test_local_align_primer_seq_fwd_rev_match(self):
        "local_align function can handle fwd/rev primers with no mismatches"
        # forward primer
        primer = DNA.makeSequence('TAGC', Name='5f')
        seq = 'TAGC'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('TAGC', 'TAGC', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        
        primer = DNA.makeSequence('TAGC', Name='5f')
        seq = 'TAGCCCCC'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('TAGC', 'TAGC', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        
        primer = DNA.makeSequence('TAGC', Name='5f')
        seq = 'CCCTAGCCCCC'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('TAGC', 'TAGC', 3)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

        # different length primer
        primer = DNA.makeSequence('GTTTAGC', Name='5f')
        seq = 'GTTTAGC'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('GTTTAGC', 'GTTTAGC', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        
        # reverse primer
        primer = DNA.makeSequence('TAGC', Name='5r')
        seq = 'TAGC'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('TAGC', 'TAGC', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        
        primer = DNA.makeSequence('TAGC', Name='5r')
        seq = 'TAGCCCCC'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('TAGC', 'TAGC', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        
        primer = DNA.makeSequence('TAGC', Name='5r')
        seq = 'CCCTAGCCCCC'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('TAGC', 'TAGC', 3)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        
    def test_local_align_primer_seq_fwd_rev_match_ambig(self):
        "local_align function can handle fwd/rev primers with ambigs"
        primer = DNA.makeSequence('TASC', Name='5f')
        seq = 'TAGC'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('TASC', 'TAGC', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)

    def test_local_align_primer_seq_mm(self):
        "local_align function can handle fwd/rev primers with mismatches"
        # forward primer
        primer = DNA.makeSequence('AAAAACTTTTT', Name='5f')
        seq = 'AAAAAGTTTTT'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('AAAAACTTTTT', 'AAAAAGTTTTT', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        
        # forward primer
        primer = DNA.makeSequence('AAAACCTTTTT', Name='5f')
        seq = 'AAAAAGTTTTT'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('AAAACCTTTTT','AAAAAGTTTTT', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        
    def test_local_align_primer_seq_indels_middle(self):
        "local_align function can handle fwd/rev primers with indels in middle of seq"

        # Insertion in target sequence
        primer = DNA.makeSequence('CGAATCGCTATCG', Name='5f')
        seq = 'CGAATCTGCTATCG'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('CGAATC-GCTATCG','CGAATCTGCTATCG', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        
        # Deletion in target sequence
        primer = DNA.makeSequence('CGAATCGCTATCG', Name='5f')
        seq = 'CGAATGCTATCG'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('CGAATCGCTATCG','CGAAT-GCTATCG', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        

        
    def test_local_align_primer_seq_multiple_mismatch_indel(self):
        "local_align function can handle fwd/rev primers with indels and mismatches"
        # multiple insertions
        primer = DNA.makeSequence('ATCGGGCGATCATT', Name='5f')
        seq = 'ATCGGGTTCGATCATT'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('ATCGGG--CGATCATT', 'ATCGGGTTCGATCATT', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        # Test that primer returned same length as input primer
        self.assertEqual(len(str(actual[0].replace("-",""))), len(str(primer)))
        
        # two deletions
        primer = DNA.makeSequence('ACGGTACAGTGG', Name='5f')
        seq = 'ACGGCAGTGG'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('ACGGTACAGTGG', 'ACGG--CAGTGG', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
  
        # deletion and mismatch
        primer = DNA.makeSequence('CATCGTCGATCA', Name='5f')
        seq = 'CCTCGTGATCA'
        # primer_hit, target, mismatch_count, hit_start
        expected = ('CATCGTCGATCA', 'CCTCGT-GATCA', 0)
        actual = local_align_primer_seq(primer, seq)
        self.assertEqual(actual, expected)
        
    def test_local_align_real_seq_primer_pairs(self):
        """local_align return correct results with real seq/primer pairs"""
        # these tests are based on manual review of the sequences

        seq = self.archaeal_16s_seq1

        f_primer = 'GYGCASCAGKCGMGAAW'
        expected = ('GYGCASCAGKCGMGAAW','ACGCAGCAGGCGAGAAA',360)
        self.assertEqual(local_align_primer_seq(f_primer,seq),expected)

        r_primer = 'AATTGGAKTCAACGCCGGR'
        expected = ('AATTGGAKTCAACGCCGGR','AATTAGATTCAACGCTGGA',998)
        self.assertEqual(local_align_primer_seq(r_primer,seq),expected)
        

         
    def test_get_primers_errors(self):
        """ Properly raises errors with invalid primer data 
        
        This test is for the primer parameters, that can be a complete text
        file of primers, a single primer name that is reference along with
        the text file, or a single primer name and manually specified sequence
        """
        
        # Should raise a value error with no data passed
        
        self.assertRaises(ValueError, get_primers)
        
        # Should raise an error if only name or sequence passed
        primer_name = "217f"
        primer_sequence = "AACTGARGSN"
        
        self.assertRaises(ValueError, get_primers, primer_name = primer_name)
        self.assertRaises(ValueError, get_primers,
         primer_sequence = primer_sequence)
         
        bad_primer_name = "F217"
        
        self.assertRaises(ValueError, get_primers,
         primer_name = bad_primer_name, primer_sequence = primer_sequence)
         
        # Raise error if primers not named correctly in primers file
        self.assertRaises(ValueError, get_primers,
         primers_data = bad_primers_file)
         
        # Raise error if primer name specified but not in file
        self.assertRaises(ValueError, get_primers,
         primers_data = primers_file, primer_name = primer_name)
        
    def test_get_primers(self):
        """ Properly builds primer data
        
        This test is for the primer parameters, that can be a complete text
        file of primers, a single primer name that is reference along with
        the text file, or a single primer name and manually specified sequence
        """
        
        # Test for manually specified name and sequence
        primer_name = "217f"
        primer_sequence = "AACTGARGSN"
        
        expected_result = [DNA.makeSequence("AACTGARGSN", Name = "217f")]
        
        actual_result = get_primers(primer_name = primer_name,
         primer_sequence = primer_sequence)
         
        self.assertEqual(actual_result, expected_result)
        
        # Test for primer file input results
        expected_result = [DNA.makeSequence('ACTAAGCCATG', Name = "65f"),
                           DNA.makeSequence('GCCTRCTAGSGATTCCGGCT', 
                            Name = "497r")]
                            
        actual_result = get_primers(primers_data = primers_file)
        
        self.assertEqual(actual_result, expected_result)
        
        # Test for correct primer when specifying name and file
        
        expected_result = [DNA.makeSequence('ACTAAGCCATG', Name = "65f")]
        
        actual_result = get_primers(primers_data = primers_file, 
         primer_name = "65f")
         
        self.assertEqual(actual_result, expected_result)
        
        
    def test_get_outf_paths(self):
        """ Tests for proper creation of output filepaths """
        
        # No actual file IO, just string manipulation
        output_dir = "analyzed_primers/"
        primer = DNA.makeSequence("ACCACGT", Name = "657r")
        fasta_filepath = "../eukaryotic_seqs.fasta"
        
        expected_graph_fp = "analyzed_primers//657r_eukaryotic_seqs.ps"
        expected_hits_fp = "analyzed_primers//657r_eukaryotic_seqs_hits.txt"
        
        actual_graph_fp, actual_hits_fp = get_outf_paths(output_dir,
         primer, fasta_filepath)
         
        self.assertEqual(actual_graph_fp, expected_graph_fp)
        self.assertEqual(actual_hits_fp, expected_hits_fp)
        
    def test_analyze_primers(self):
        """ Test for overall program module """
        
                    
        fasta_fps = self.sample_fasta_file
        verbose = False
        output_dir = self.output_dir
        primers_filepath = None
        primer_name = '27f'
        primer_sequence = 'AGAGTTTGATCMTGGCTCAG'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        
        hits_file_output = self.output_dir +\
         "27f_" + basename(self.sample_fasta_file).split('.')[0] + "_hits.txt"

        histogram_file_output = self.output_dir +\
         "27f_" + basename(self.sample_fasta_file).split('.')[0] + ".ps"
         
        analyze_primers(fasta_fps, verbose, output_dir, primers_filepath,
         primer_name, primer_sequence, tp_len, last_base_mm, tp_mm,
         non_tp_mm, tp_gap, non_tp_gap)
         
        
        self._files_to_remove.append(hits_file_output)
        self._files_to_remove.append(histogram_file_output)
        
        expected_hits = self.expected_bacterial_seqs_hits_file_27f
        # Have to manually append tmp file name
        expected_hits[1] += basename(self.sample_fasta_file)
        
        # Test hits result
        hits_result_f = open(hits_file_output, "U")
        actual_hits_result = [line.strip() for line in hits_result_f]
        
        self.assertEqual(actual_hits_result,\
         self.expected_bacterial_seqs_hits_file_27f)
        
        # Can only test existance of .ps output file
        self.assertEqual(isfile(histogram_file_output), True)
        
    def test_score_primer_perfect_match(self):
        """ Test for proper mm, gap, and score calcuations for primers """
        
        
        # Test for perfect score
        primer = DNA.makeSequence('ACGTACGTAA', Name='477f')
        primer_hit = 'ACGTACGTAA'
        target_hit = 'ACGTACGTAA'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        expected_weighted_score = 0
        expected_non_tp_gaps = 0
        expected_tp_gaps = 0
        expected_non_tp_mismatches = 0
        expected_tp_mismatches = 0
        expected_last_base_mismatches = False
        
        weighted_score, non_tp_gaps, tp_gaps, non_tp_mismatches,\
         tp_mismatches, last_base_mismatches = score_primer(primer,
         primer_hit, target_hit, tp_len, last_base_mm, tp_mm, non_tp_mm,
         tp_gap, non_tp_gap)
                 
        self.assertEqual(weighted_score, expected_weighted_score)
        self.assertEqual(non_tp_gaps, expected_non_tp_gaps)
        self.assertEqual(tp_gaps, expected_tp_gaps)
        self.assertEqual(non_tp_mismatches, expected_non_tp_mismatches)
        self.assertEqual(tp_mismatches, expected_tp_mismatches)
        self.assertEqual(last_base_mismatches, expected_last_base_mismatches)
        
    def test_score_primer_non_tp_mm(self):
        """ Test for proper mm, gap, and score calcuations for primers """
        
        
        # Test for perfect score
        primer = DNA.makeSequence('ACGTACGTAA', Name='477f')
        primer_hit = 'AGGTACGTAA'
        target_hit = 'ACGTACGTAA'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        expected_weighted_score = 0.4
        expected_non_tp_gaps = 0
        expected_tp_gaps = 0
        expected_non_tp_mismatches = 1
        expected_tp_mismatches = 0
        expected_last_base_mismatches = False
        
        weighted_score, non_tp_gaps, tp_gaps, non_tp_mismatches,\
         tp_mismatches, last_base_mismatches = score_primer(primer,
         primer_hit, target_hit, tp_len, last_base_mm, tp_mm, non_tp_mm,
         tp_gap, non_tp_gap)
                 
        self.assertEqual(weighted_score, expected_weighted_score)
        self.assertEqual(non_tp_gaps, expected_non_tp_gaps)
        self.assertEqual(tp_gaps, expected_tp_gaps)
        self.assertEqual(non_tp_mismatches, expected_non_tp_mismatches)
        self.assertEqual(tp_mismatches, expected_tp_mismatches)
        self.assertEqual(last_base_mismatches, expected_last_base_mismatches)
        
    def test_score_primer_tp_mm(self):
        """ Test for proper mm, gap, and score calcuations for primers """
        
        
        # Test for perfect score
        primer = DNA.makeSequence('ACGTACGTAA', Name='477f')
        primer_hit = 'ACGTACGTTA'
        target_hit = 'ACGTACGTAA'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        expected_weighted_score = 1.0
        expected_non_tp_gaps = 0
        expected_tp_gaps = 0
        expected_non_tp_mismatches = 0
        expected_tp_mismatches = 1
        expected_last_base_mismatches = False
        
        weighted_score, non_tp_gaps, tp_gaps, non_tp_mismatches,\
         tp_mismatches, last_base_mismatches = score_primer(primer,
         primer_hit, target_hit, tp_len, last_base_mm, tp_mm, non_tp_mm,
         tp_gap, non_tp_gap)
                 
        self.assertEqual(weighted_score, expected_weighted_score)
        self.assertEqual(non_tp_gaps, expected_non_tp_gaps)
        self.assertEqual(tp_gaps, expected_tp_gaps)
        self.assertEqual(non_tp_mismatches, expected_non_tp_mismatches)
        self.assertEqual(tp_mismatches, expected_tp_mismatches)
        self.assertEqual(last_base_mismatches, expected_last_base_mismatches)
        
    def test_score_primer_last_base_mm(self):
        """ Test for proper mm, gap, and score calcuations for primers """
        
        
        # Test for perfect score
        primer = DNA.makeSequence('ACGTACGTAA', Name='477f')
        primer_hit = 'ACGTACGTAA'
        target_hit = 'ACGTACGTAC'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        expected_weighted_score = 3.0
        expected_non_tp_gaps = 0
        expected_tp_gaps = 0
        expected_non_tp_mismatches = 0
        expected_tp_mismatches = 0
        expected_last_base_mismatches = True
        
        weighted_score, non_tp_gaps, tp_gaps, non_tp_mismatches,\
         tp_mismatches, last_base_mismatches = score_primer(primer,
         primer_hit, target_hit, tp_len, last_base_mm, tp_mm, non_tp_mm,
         tp_gap, non_tp_gap)
                 
        self.assertEqual(weighted_score, expected_weighted_score)
        self.assertEqual(non_tp_gaps, expected_non_tp_gaps)
        self.assertEqual(tp_gaps, expected_tp_gaps)
        self.assertEqual(non_tp_mismatches, expected_non_tp_mismatches)
        self.assertEqual(tp_mismatches, expected_tp_mismatches)
        self.assertEqual(last_base_mismatches, expected_last_base_mismatches)
        
        
    def test_score_primer_non_tp_gaps(self):
        """ Test for proper mm, gap, and score calcuations for primers """
        
        
        # Test for perfect score
        primer = DNA.makeSequence('ACGTACGTAA', Name='477f')
        primer_hit = 'A-GTACGTAA'
        target_hit = 'ACG-ACGTAA'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        expected_weighted_score = 2.0
        expected_non_tp_gaps = 2
        expected_tp_gaps = 0
        expected_non_tp_mismatches = 0
        expected_tp_mismatches = 0
        expected_last_base_mismatches = False
        
        weighted_score, non_tp_gaps, tp_gaps, non_tp_mismatches,\
         tp_mismatches, last_base_mismatches = score_primer(primer,
         primer_hit, target_hit, tp_len, last_base_mm, tp_mm, non_tp_mm,
         tp_gap, non_tp_gap)
                 
        self.assertEqual(weighted_score, expected_weighted_score)
        self.assertEqual(non_tp_gaps, expected_non_tp_gaps)
        self.assertEqual(tp_gaps, expected_tp_gaps)
        self.assertEqual(non_tp_mismatches, expected_non_tp_mismatches)
        self.assertEqual(tp_mismatches, expected_tp_mismatches)
        self.assertEqual(last_base_mismatches, expected_last_base_mismatches)
        
    def test_score_primer_tp_gap(self):
        """ Test for proper mm, gap, and score calcuations for primers """
        
        
        # Test for perfect score
        primer = DNA.makeSequence('ACGTACGTAA', Name='477f')
        primer_hit = 'ACGTACGT-A'
        target_hit = 'ACGTACGTAA'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        expected_weighted_score = 3.0
        expected_non_tp_gaps = 0
        expected_tp_gaps = 1
        expected_non_tp_mismatches = 0
        expected_tp_mismatches = 0
        expected_last_base_mismatches = False
        
        weighted_score, non_tp_gaps, tp_gaps, non_tp_mismatches,\
         tp_mismatches, last_base_mismatches = score_primer(primer,
         primer_hit, target_hit, tp_len, last_base_mm, tp_mm, non_tp_mm,
         tp_gap, non_tp_gap)
                 
        self.assertEqual(weighted_score, expected_weighted_score)
        self.assertEqual(non_tp_gaps, expected_non_tp_gaps)
        self.assertEqual(tp_gaps, expected_tp_gaps)
        self.assertEqual(non_tp_mismatches, expected_non_tp_mismatches)
        self.assertEqual(tp_mismatches, expected_tp_mismatches)
        self.assertEqual(last_base_mismatches, expected_last_base_mismatches)
        
    def test_score_primer_perfect_match_degen(self):
        """ Test for proper mm, gap, and score calcuations for primers """
        
        
        # Test for perfect score
        primer = DNA.makeSequence('ACGTACGTAA', Name='477f')
        primer_hit = 'AYGTRCGTWA'
        target_hit = 'ACGTACGTAA'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        expected_weighted_score = 0
        expected_non_tp_gaps = 0
        expected_tp_gaps = 0
        expected_non_tp_mismatches = 0
        expected_tp_mismatches = 0
        expected_last_base_mismatches = False
        
        weighted_score, non_tp_gaps, tp_gaps, non_tp_mismatches,\
         tp_mismatches, last_base_mismatches = score_primer(primer,
         primer_hit, target_hit, tp_len, last_base_mm, tp_mm, non_tp_mm,
         tp_gap, non_tp_gap)
                 
        self.assertEqual(weighted_score, expected_weighted_score)
        self.assertEqual(non_tp_gaps, expected_non_tp_gaps)
        self.assertEqual(tp_gaps, expected_tp_gaps)
        self.assertEqual(non_tp_mismatches, expected_non_tp_mismatches)
        self.assertEqual(tp_mismatches, expected_tp_mismatches)
        self.assertEqual(last_base_mismatches, expected_last_base_mismatches)
        
    def test_score_primer_perfect_match_reverse(self):
        """ Test for proper mm, gap, and score calcuations for primers """
        
        
        # Test for perfect score
        primer = DNA.makeSequence('ACGTACGTAA', Name='825r')
        primer_hit = 'ACGTACGTAA'
        target_hit = 'ACGTACGTAA'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        expected_weighted_score = 0
        expected_non_tp_gaps = 0
        expected_tp_gaps = 0
        expected_non_tp_mismatches = 0
        expected_tp_mismatches = 0
        expected_last_base_mismatches = False
        
        weighted_score, non_tp_gaps, tp_gaps, non_tp_mismatches,\
         tp_mismatches, last_base_mismatches = score_primer(primer,
         primer_hit, target_hit, tp_len, last_base_mm, tp_mm, non_tp_mm,
         tp_gap, non_tp_gap)
                 
        self.assertEqual(weighted_score, expected_weighted_score)
        self.assertEqual(non_tp_gaps, expected_non_tp_gaps)
        self.assertEqual(tp_gaps, expected_tp_gaps)
        self.assertEqual(non_tp_mismatches, expected_non_tp_mismatches)
        self.assertEqual(tp_mismatches, expected_tp_mismatches)
        self.assertEqual(last_base_mismatches, expected_last_base_mismatches)
        
    def test_score_primer_tp_mm_reverse(self):
        """ Test for proper mm, gap, and score calcuations for primers """
        
        
        # Test for perfect score
        primer = DNA.makeSequence('ACGTACGTAA', Name='825r')
        primer_hit = 'ACGTACGTAA'
        target_hit = 'ACTTACGTAA'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        expected_weighted_score = 1.0
        expected_non_tp_gaps = 0
        expected_tp_gaps = 0
        expected_non_tp_mismatches = 0
        expected_tp_mismatches = 1
        expected_last_base_mismatches = False
        
        weighted_score, non_tp_gaps, tp_gaps, non_tp_mismatches,\
         tp_mismatches, last_base_mismatches = score_primer(primer,
         primer_hit, target_hit, tp_len, last_base_mm, tp_mm, non_tp_mm,
         tp_gap, non_tp_gap)
                 
        self.assertEqual(weighted_score, expected_weighted_score)
        self.assertEqual(non_tp_gaps, expected_non_tp_gaps)
        self.assertEqual(tp_gaps, expected_tp_gaps)
        self.assertEqual(non_tp_mismatches, expected_non_tp_mismatches)
        self.assertEqual(tp_mismatches, expected_tp_mismatches)
        self.assertEqual(last_base_mismatches, expected_last_base_mismatches)
        
    def test_score_primer_non_tp_gap_reverse(self):
        """ Test for proper mm, gap, and score calcuations for primers """
        
        
        # Test for perfect score
        primer = DNA.makeSequence('ACGTACGTAA', Name='825r')
        primer_hit = 'ACGTACG-AA'
        target_hit = 'ACGTACGTAA'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        expected_weighted_score = 1.0
        expected_non_tp_gaps = 1
        expected_tp_gaps = 0
        expected_non_tp_mismatches = 0
        expected_tp_mismatches = 0
        expected_last_base_mismatches = False
        
        weighted_score, non_tp_gaps, tp_gaps, non_tp_mismatches,\
         tp_mismatches, last_base_mismatches = score_primer(primer,
         primer_hit, target_hit, tp_len, last_base_mm, tp_mm, non_tp_mm,
         tp_gap, non_tp_gap)
                 
        self.assertEqual(weighted_score, expected_weighted_score)
        self.assertEqual(non_tp_gaps, expected_non_tp_gaps)
        self.assertEqual(tp_gaps, expected_tp_gaps)
        self.assertEqual(non_tp_mismatches, expected_non_tp_mismatches)
        self.assertEqual(tp_mismatches, expected_tp_mismatches)
        self.assertEqual(last_base_mismatches, expected_last_base_mismatches)
        
    def test_score_primer_last_base_mm_reverse(self):
        """ Test for proper mm, gap, and score calcuations for primers """
        
        
        # Test for perfect score
        primer = DNA.makeSequence('ACGTACGTAA', Name='825r')
        primer_hit = 'ACGTACGTAA'
        target_hit = 'GCGTACGTAA'
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        expected_weighted_score = 3.0
        expected_non_tp_gaps = 0
        expected_tp_gaps = 0
        expected_non_tp_mismatches = 0
        expected_tp_mismatches = 0
        expected_last_base_mismatches = True
        
        weighted_score, non_tp_gaps, tp_gaps, non_tp_mismatches,\
         tp_mismatches, last_base_mismatches = score_primer(primer,
         primer_hit, target_hit, tp_len, last_base_mm, tp_mm, non_tp_mm,
         tp_gap, non_tp_gap)
                 
        self.assertEqual(weighted_score, expected_weighted_score)
        self.assertEqual(non_tp_gaps, expected_non_tp_gaps)
        self.assertEqual(tp_gaps, expected_tp_gaps)
        self.assertEqual(non_tp_mismatches, expected_non_tp_mismatches)
        self.assertEqual(tp_mismatches, expected_tp_mismatches)
        self.assertEqual(last_base_mismatches, expected_last_base_mismatches)

        
    def test_hits_seq_end(self):
        """ Test to detection of primer abutting sequence end """

        
        # Test for hitting start of sequence
        test_seq = "AAACTGACCTAGACCTAGGGTATTACG"
        test_hit_start = 0
        primer_len = 15
        
        expected_result = True
        
        # Starts at zero, so should return True
        self.assertEqual(hits_seq_end(test_seq, test_hit_start, primer_len),
         expected_result)
         
        # Test for hitting end of sequence
        test_seq = "AAACTGACCTAGACCTAGGGTATTACG"
        test_hit_start = 20
        primer_len = 15
        
        expected_result = True
        
        # Start position + length of primer hits right end of sequence
        # so should return True
        self.assertEqual(hits_seq_end(test_seq, test_hit_start, primer_len),
         expected_result)
         
        # Test for not hitting end of sequence
        test_seq = "AAACTGACCTAGACCTAGGGTATTACG"
        test_hit_start = 5
        primer_len = 15
        
        expected_result = False
        
        # Should be in middle of sequence, return False
        self.assertEqual(hits_seq_end(test_seq, test_hit_start, primer_len),
         expected_result)
        
        
        
        
    def test_get_yaxis_max(self):
        """ Gets largest bin size from list of lists """
        
        sample_bin_list = [[1, 2, 1], [0, 0, 0], [1, 0, 0], [0, 0, 0], 
         [1.40, 0.80, 0.40], [0, 0, 0]]
         
        # largest bin is 3 (3 instances of 0 in multiple bins)
        expected_max_size = 3
        
        actual_result = get_yaxis_max(sample_bin_list)
        
        self.assertEqual(actual_result, expected_max_size)
    
    def test_primer_hit_histogram(self):
        """ Test building of primer hit scores histogram
        
        Not really possible to test this function-all data that is passed
        to pylab is for building histogram, no returned data.  Will call
        function with normal data to ensure no errors raised when calling
        function normally."""
        
        sample_mm_data = [1, 0, 1, 2, 0, 0, 0]
        figure_title = "test title"
        x_label = "sample mismatches"
        primer_hit_histogram(sample_mm_data, figure_title, x_label)
        
    def test_write_primer_histogram(self):
        """ Test writing of primer histograms 
        
        Not really possible to test this module, except for creation of 
        histogram file.  Data are generated before being passed to this file,
        which sorts data into various components of the histogram."""
        
        histogram_file_output = self.output_dir +\
         "test_histogram_file.ps"
        
        self._files_to_remove.append(histogram_file_output)
        
        
        write_primer_histogram(self.sample_bacterial_seqs_expected_hist,
         histogram_file_output)
        
        # Can only test existance of .ps output file
        self.assertEqual(isfile(histogram_file_output), True)
        
    def test_pair_hmm_align_unaligned_seqs_exact(self):
        """ Properly returns local pairwise alignment of 2 seqs """
        
        # Test exact match of primer to the test archaeal sequence
        primer = "GAAAGACCTCGG"
        
        # since primer is exact match, result should yield exact match to 
        # primer sequence
        expected_primer_hit = "GAAAGACCTCGG"
        expected_target_hit = "GAAAGACCTCGG"
        
        result = pair_hmm_align_unaligned_seqs([primer, self.archaeal_16s_seq1])
        
        primer_hit = str(result.Seqs[0])
        target_hit = str(result.Seqs[1])
        
        self.assertEqual(primer_hit, expected_primer_hit)
        self.assertEqual(target_hit, expected_target_hit)
        
    def test_pair_hmm_align_unaligned_seqs_degen(self):
        """ Properly returns local pairwise alignment of 2 seqs """
        # Test degenerate match of primer
        primer = "GAAARACYTCGG"
        
        # since primer is exact match, result should yield exact match to 
        # primer sequence
        expected_primer_hit = "GAAARACYTCGG"
        expected_target_hit = "GAAAGACCTCGG"
        
        result = pair_hmm_align_unaligned_seqs([primer, self.archaeal_16s_seq1])
        
        primer_hit = str(result.Seqs[0])
        target_hit = str(result.Seqs[1])
        
        self.assertEqual(primer_hit, expected_primer_hit)
        self.assertEqual(target_hit, expected_target_hit)
        
    def test_pair_hmm_align_unaligned_seqs_mismatch(self):
        """ Properly returns local pairwise alignment of 2 seqs """
        # Test mismatch of primer
        primer = "GATAGACCTCGG"
        
        # Should return same position, with one mismatch in primer hit
        expected_primer_hit = "GATAGACCTCGG"
        expected_target_hit = "GAAAGACCTCGG"
        
        result = pair_hmm_align_unaligned_seqs([primer, self.archaeal_16s_seq1])
        
        primer_hit = str(result.Seqs[0])
        target_hit = str(result.Seqs[1])
        
        self.assertEqual(primer_hit, expected_primer_hit)
        self.assertEqual(target_hit, expected_target_hit)
        
    def test_pair_hmm_align_unaligned_seqs_gap_primer(self):
        """ Properly returns local pairwise alignment of 2 seqs """
        
        # Test exact match of primer, with one deletion
        # to the test archaeal sequence.  Because of the fairly high
        # penalty for starting a gap, a longer sequence is needed here
        primer = "TCGGCGAATTGCCCGTAATACGTAG"
        
        # since primer is exact match but with one deletion in the primer
        # should return gap in primer sequence
        expected_primer_hit = "TCGGCGAATTGC-CCGTAATACGTAG"
        expected_target_hit = "TCGGCGAATTGCTCCGTAATACGTAG"
        
        result = pair_hmm_align_unaligned_seqs([primer, self.archaeal_16s_seq1])
        
        primer_hit = str(result.Seqs[0])
        target_hit = str(result.Seqs[1])
        
        self.assertEqual(primer_hit, expected_primer_hit)
        self.assertEqual(target_hit, expected_target_hit)
        
    def test_pair_hmm_align_unaligned_seqs_insertion_primer(self):
        """ Properly returns local pairwise alignment of 2 seqs """
        
        # Test exact match of primer, with one insertion
        # relative to the test archaeal sequence.  Because of the fairly high
        # penalty for starting a gap, a longer sequence is needed here
        primer = "TCGGCGAATTGCTCTCGTAATACGTAG"
        
        # since primer is exact match but with one insertion in the primer
        # should return gap in target sequence
        expected_primer_hit = "TCGGCGAATTGCTCTCGTAATACGTAG"
        expected_target_hit = "TCGGCGAATTGCTC-CGTAATACGTAG"
        
        result = pair_hmm_align_unaligned_seqs([primer, self.archaeal_16s_seq1])
        
        primer_hit = str(result.Seqs[0])
        target_hit = str(result.Seqs[1])
        
        
        self.assertEqual(primer_hit, expected_primer_hit)
        self.assertEqual(target_hit, expected_target_hit)
        
        
        
    def test_get_hits_data(self):
        """ Properly generates hits data for text and histogram files """
        
        # Using real primer, real sequences to test results
        primer = DNA.makeSequence('GGACTACVSGGGTATCTAAT', Name="806r")
        primer_name = "806r"
        fasta_data = open(self.sample_fasta_file, "U")
        # Use default settings for scoring
        tp_len = 5
        last_base_mm = 3.0
        tp_mm = 1.0
        non_tp_mm = 0.4
        tp_gap = 3.0
        non_tp_gap = 1.0
        
        actual_hits_data, actual_hist_data = get_hits_data(primer, primer_name,
         fasta_data, tp_len, last_base_mm, tp_mm, non_tp_mm, tp_gap,
         non_tp_gap)
         
        expected_hits = self.sample_bacterial_seqs_expected_hits
        expected_hist = self.sample_bacterial_seqs_expected_hist
        
        # Add in tmp fasta file name, which is stored in the output file
        expected_hits[1] += basename(self.sample_fasta_file)
        expected_hist[6] += basename(self.sample_fasta_file)
         
        self.assertEqual(actual_hits_data, expected_hits)
        self.assertEqual(actual_hist_data, expected_hist)
        
    def test_generate_hits_file_and_histogram(self):
        """ Generates proper hits file, creates histogram file """
        
        primers = [DNA.makeSequence('GTGCCAGCMGCCGCGGTAA', '515f')]
        primer_ids = ['515f']
        fasta_filepaths = [self.sample_fasta_file]
        max_mismatches = 4
        output_dir = self.output_dir
        verbose = False
        tp_len = 5
        last_base_mm = 3
        tp_mm = 1
        non_tp_mm = 0.4
        tp_gap = 3
        non_tp_gap = 1
        
         
        hits_file_output = self.output_dir +\
         "515f_" + basename(self.sample_fasta_file).split('.')[0] + "_hits.txt"

        histogram_file_output = self.output_dir +\
         "515f_" + basename(self.sample_fasta_file).split('.')[0] + ".ps"
        
        self._files_to_remove.append(hits_file_output)
        self._files_to_remove.append(histogram_file_output)
        
        generate_hits_file_and_histogram(primers, primer_ids, fasta_filepaths,
         max_mismatches, output_dir, verbose, tp_len, last_base_mm, tp_mm,
         non_tp_mm, tp_gap, non_tp_gap)
         
        expected_hits = self.expected_bacterial_seqs_hits_file2
        # Have to manually append tmp file name
        expected_hits[1] += basename(self.sample_fasta_file)
         
        # Test hits result
        hits_result_f = open(hits_file_output, "U")
        actual_hits_result = [line.strip() for line in hits_result_f]
        self.assertEqual(actual_hits_result,\
         self.expected_bacterial_seqs_hits_file2)
        
        # Can only test existance of .ps output file
        self.assertEqual(isfile(histogram_file_output), True)
        
    
# Placing large strings at the end for better readability
archaeal_16s_seq1 = "AGACCATTGCTATCCGAGTTCTACTAAGCCATGCGAGCTGAGAGGTGAAAGACCTCGGCGAATTGCTCCGTAATACGTAGTCAACTTGCCCTCAAGTGGGAGATAATCTCGGGAAACTGAGGCTAATATCCCATAGTTTTACATTACTGGAATGTGTGTAAAGCGAAAGCCACGGCGCTTGAGGATAGGACTGCGCCTGATTAGGCTGTTGTTGGTGTAAAGATACCACAATATGGCTTTATAGAACCAGCAAACCGATATGGAGGTAAGGGTATTCTTTCCTCCACGGCAAAATCAGTACGGGCTATGGAAGTAGGAGCCCGGAGATGGATTCTGAGACATGAATCCAGGTACTACGGTACGCAGCAGGCGAGAAACCTTTGCAATCGGTTAACGCCGGACAGGGGAACTTGGAGTGTTCTAGTTTTACTAGAACTTTTGCTCACTGAAAACGGGTGAGTGAATAAGGGCTGGCCAAGACGGGTGCCAGCCGCCACGGGTAATACCCGCAGCCCAAGTGATGGTCACGATTATTGGGTCTAAAGCGTCCGTAGCCGGTCTAACAAGTTCCTGGTGAAATCTTACAGCAAACTGTAAGGCTTGCTGGGGATACTGTTAGACTTGAGACCGGGAGAGGACAGAGGTACTCGTAGGGTAGGGGTTAAATCCTATAATCCTACGGGGACCACCTGTGGCGAAAGCGTCTGTCTAGAACGGGTCTGACGGTGAGGGACGAAAGCTAGGGGAGCGATCGGGATTAGATATGGAGAAGGAAAAAGTAATTCAATAAACGAAAGTTTAGAGAATCACCTATTCCAGACCCACTAAAAACCCCGGTAGTCCTAGCTGTAAACATTGCCCGCTTGGTGTTGCAGACCTCTTGGGGGTTTGCAGTGCCGGAGCGTAGGTGTTAAGCGGGCCACCTGGGGAGTACAGTCGCAAGGCCGAAACTTAAGGGAATTGGCGGGGGAGCACACAAGGGGTGGGCAGTGCGGTTCAATTAGATTCAACGCTGGAAATCTTACCAGGGGCGACAGCAGGATGAGGGTCAGTCTGAAGGGCTTACCAGACAAGCTGAGAGGTGGTGCATGGCCATCGTCAGCTCGTACCGTGAGGCGTGCCGTTAAGTCGGTTAACGAGCGAGACCCACATTTCATGTTGCAACTATACTTTCCGAAGTATAGGCACTCATGAGAGACCGCTGGTGATAAACCAGAGGAAGGTGTGGGCGACGGTAGGTCTGTATGCCCCGGATCTCCTGGGCTACACGCGCTGCACAATGCGTGCCACAATGGGAAGCAACTCCGAGAGGAGAAGCGAATCCCCTAAAAGCACGCTTAGTTCAGATTGAGGGTTGCAACTCACCCTCATGAAGCCGGAATCCCTAGTAGGCGAATGTCACTAAGTAGTAGACTGCTAATGAATTACATGCACCACGCTTGATTAAAGATGATAGACAGAGGGGAAAGAATCAGTCTGCTATCGTGAAATCGTTCGCCGAATACGTCCCTGCTCC"
        

primers_file = ['#Test primers file',
                '65f\tACTAAGCCATG',
                '497r\tGCCTRCTAGSGATTCCGGCT']
                
bad_primers_file = ['#Test primers file',
                    '65f\tACTAAGCCATG',
                    'R497\tGCCTRCTAGSGATTCCGGCT']
                    

# 3 bacterial sequences from the Silva database.
sample_bacterial_seqs = """>AB240712 1 917 Bacteria/Actinobacteria/Acidimicrobiaceae/Microthrix/uncultured/uncultured
AGAGTTTGATCCTGGCTCAGGGTGAACGCTGGCGGCGGGCTTAATGCATGCAAGTCGAGCGTGCCTGTGACTTCGGTCACAATAGGAAAGCGGCAGACGGCTGAGTAACGCGTGAGTAACTTGCCCTTTGGGGGGGTATAGCCTTGTGAAAACGAGGATAATCCCGCATAAGATCCCCAAGCCCTGGTTTGGGGATGAAAGCCTTCGGGCGCCAGAGGAGAGGCTCGCGTCCTATCAGCTAGTTGGAGAGGTAACGGCTCACCAAGGCATCGACGGGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCATACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGAAAGCCCGATGCAGCGACGCCGCGTGCGGGAAGAAGGCCCTAGGGTTGTAAACCGCTTTCAGTAGGGAAGAAAATGACGGTACCTACAGAAGAAGGTGCGGCCAACTACGTGCCAGCAGCCGCGGTGACACGTAGGCACCAAGCGTTATCCGGATTTATTGGGCGTAAAGAGCTCGTAGGCGGCTCGGTAAGTCGGGTGTGAAAACTCTGGGCTCAACCCAGAGAGGCCACCCGATACTGCTGTGGCTAGAGTACGGTAGGGGAGCGGGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCAGCGGCGAAGGCGCCGCTCTGGGCCGTTACTGACGCTGAGGAGCGAAGCATGGGTAGCAAACAGGATTAGATACCTCTGGTAGTCCATGCCGTAAACGTTGGGCACTAAGTGTGGGTCTCAACACAACGAGAATCTGCGCCGTCGCTAACGCATTAAGATGCCCCGCCTGGGGGAGTACGGTCGCAGACTAAAACTCAAGGGATTTGACGCGGGGCCCGCACAAGCGGGT
>EF446229 1 937 Bacteria/Alphaproteobacteria/Rickettsiales mitochondria/mitochondria
TACGGTTACCCTTGTTACGACTTAAAGCCCTAACTTTAAATAAGCTGTAATACAGCTTTACTATGAGTTTTATACAGAAGAGTAGAATTTTATATGAAAGGGTGAAATCTGGAAATATATAAAGGAATATCAATAGCGAAGGCAACTCTTTAGTAAAAACTGACGTTGAGGAACGAAAGTGTAGGTAGCAAACAGGATTAGATACCCTGGTAGTCTACACTGTAAATTCTGAGTGCTGTAATTTAATGTAAAAGTTAATTCTTTTATTTAGATTTTTAAGTTAACACCATAAGCACTCCGCCTGAGGAGTATGATCGCAAGGTTGGAACTCAAAGGAATTGACGGGGGCTAGGACTAGTGGTGGAGTATGTGGTTTAATGCGATAATCCGCGCAAAACCTTACCAATTTTTGACATAATAGCTAAAATTTTATTTTACTTTGAGGTAAAATAATATATTTTAACTATTACAGGTGTTGCACGACTGTCGTCAGCTCGTGTTGTGAGATGTTTGGTTAATTCCTTTTAAACGAGCGTAACTCTTATCTTTAATTAAATAATTTACTATTTTTAAAAATTTTTTTTTAAAAGTTTTATTTGAAACTGCTAATAAATAAGTTAGAGGAAGGTGAGAACTAAGCCAAGTCAATGTGACCCTTATAAATTGGGCTACACACGTGCTACAATGGATACTACAGAAAGTTACAAAAATTCAGTTTTTGTTAATCTATAAAAGTCTCCCCAGTTCGGATTGTAGGCTGAAACTCGCTTACATGAAGTTGGAATCACTAGTAATCGCGAATCAGAACGTCGTGGTGAATATGTAAACTAGCCTTGTACACACCGCCCGTCACATGCAAAGAGTTAGCTTTTTTTGAAAGTTTAACTTTGTTAAGCTGATTTAAAAGGTAGTAATTTGGATGAAGTCGTAACAAGGT
>EU434431 1 1200 Bacteria/Beta Gammaproteobacteria/Gammaproteobacteria_1/Pseudomonadales_4/Pseudomonas fluorescens et rel./Pseudomonas tolaasii et rel.
AATCTGCACACACCTTCTCAGACTTCCTACGGAGCCAGCATTTGCGGATTATTTGGACAATCGGGCCGAACGCCTTGATTCCAGCCATGCGGGTGTGTGAATGATGGCTCTTCGAATTTCTAACACGCACTTTAAGGTAGCGAGGAAAGCGTTGTAGCATTCATTACTCCTCTATTTGATGTTCACCGACAGAACTAAGCACTCGGTTAATCTCTGTGCTCAGCAGCCGCGGTTAATACAGAGGGGTGCAAGTGTTTAATCGGAAATTACTGGGCATAACAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGGCCTCCACGTGGGAACTGCATTCAAAACTGACTGACTAGAGTATGGTAGAGGGTGGTGGAATCTCCTGTGTCGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCAACTAGCCGTTGGGAGCCTTGAGCTCTTAGTGGCGCAGCTAACGCATTAAGTTGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGCCTTGACATCCAATGAACTTTCTAGAGATAGATTGGTGCCTTCGGGAACATTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGTAACGAGCGCAACCCTTGTCCTTAGTTACCAGCACGTTATGGTGGGCACTCTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGGCCTGGGCTACACACGTGCTACAATGGTCGGTACAGAGGGTTGCCAAGCCGCGAGGTGGAGCTAATCCCAGAAAACCGATTGTAGTCCGGATCGCAGTCTGCAACTCGACTGCGTGAAGTCGGAATCGCTAGTAATCGCGAATCAGAATGTCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCACCAGAAGTAGCTAGTCTAACCTTCGGGAGGACGGTATACATCAGTATAATCATCGCTTTGAG"""

sample_bacterial_seqs_expected_hits =\
 ["# Primer: 806r 5'-GGACTACVSGGGTATCTAAT-3'", 
 '# Input fasta file: ', 
 '# Parameters', "# 3' length: 5", 
 "# non 3' mismatch penalty: 0.40 per mismatch", 
 "# 3' mismatch penalty: 1.00 per mismatch", 
 '# last base mismatch penalty: 3.00', 
 "# non 3' gap penalty: 1.00 per gap", 
 "# 3' gap penalty: 3.00 per gap", 
 "# Note - seq hit and primer hit are the best local pairwise alignment results for a given sequence and primer pair.  A gap in seq hit represents a deletion in the sequence, whereas a gap in the primer hit signifies an insertion in the target sequence.\n#\n# seq ID, seq hit, primer hit, hit start position, non 3' mismatches, 3' mismatches (except last base), last base mismatch, non 3' gaps, 3' gaps, overall weighted score, hits sequence end ", 
 'AB240712,ATTAGATACCTCTGGTAGTCC,ATTAGATACC-CSBGTAGTCC,752,1,0,False,1,0,1.4,False', 
 'EF446229,ATTAGATACCCTGGTAGTCT,ATTAGATACCCSBGTAGTCC,196,2,0,False,0,0,0.8,False', 
 'EU434431,ATTAGATACCCTGGTAGTCC,ATTAGATACCCSBGTAGTCC,495,1,0,False,0,0,0.4,False']
 
sample_bacterial_seqs_expected_hist =\
 [[1, 2, 1], [0, 0, 0], [1, 0, 0], [0, 0, 0], 
 [1.40, 0.80, 0.40], 
 [0, 0, 0], 
 "806r; Degeneracy: 6; GC content 0.45 - 0.50\n5'-GGACTACVSGGGTATCTAAT-3'\nSequences tested: ", "3' length: 5 nucleotides\nWeighted score = non-3' mismatches * 0.40 + 3' mismatches * 1.00 + non 3' gaps * 1.00 + 3' gaps * 3.00\nAn additional 3.00 penalty is assigned if the final 3' base mismatches\nWeighted score is rounded to the nearest whole number in this graphical display"]
 
# expected hits file result when using the 515f primer
expected_bacterial_seqs_hits_file2 =\
"""# Primer: 515f 5'-GTGCCAGCMGCCGCGGTAA-3'
# Input fasta file: 
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
AB240712,GTGCCAGCAGCCGCGGTGA,GTGCCAGCMGCCGCGGTAA,480,0,1,False,0,0,1.0,False
EF446229,ATGCGATAATCCGCGCAAA,GTGCCAGCMGCCGCGGTAA,377,5,2,False,0,0,4.0,False
EU434431,GTGCTCAGCAGCCGCGGTTA,GTGC-CAGCMGCCGCGGTAA,215,0,1,False,1,0,2.0,False""".split('\n')

expected_bacterial_seqs_hits_file_27f =\
"""# Primer: 27f 5'-AGAGTTTGATCMTGGCTCAG-3'
# Input fasta file: 
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
AB240712,AGAGTTTGATCCTGGCTCAG,AGAGTTTGATCMTGGCTCAG,0,0,0,False,0,0,0.0,True
EF446229,AGGAATTGACGGGGGCTAGG,AGAGTTTGATCMTGGCTCAG,333,7,2,False,0,0,4.8,False
EU434431,TCAGTATAATCATCGCTTTG,AGAGTTTGATCMTGGCTCAG,1178,5,2,False,0,0,4.0,False""".split('\n')

if __name__ == "__main__":
    main()
