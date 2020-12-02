#!/usr/bin/env python
#test_generate_primers_denovo.py
#author William Walters
#created 7-27-09

from __future__ import division

__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

from os import remove

from cogent3.core.alphabet import AlphabetError
from cogent3.util.unit_test import TestCase, main
from primerprospector.old_cogent import get_tmp_filename
from cogent3.util.misc import remove_files

from primerprospector.generate_primers_denovo import ProspectivePrimer, \
get_number_seqs_for_primer,build_seq_data,convert_to_numeric,\
get_corrected_index,del_primers,append_primer_hit,\
find_specific_primer_matches,get_deletion_threshold,\
find_sensitive_primer_matches,calculate_percent_match, search_sequences



class GeneratePrimersDeNovoTests(TestCase):
    """ Unit tests for the sequence_searcher.py module """

    def setUp(self):
        self.initial_primers={}
        
        
    def test_append_primer_hit_raises_errors(self):
        """ append_primer_hit raises properly creates primer object """
        test_primer=ProspectivePrimer("ATCGA",15459,700)
        bad_primer="ATCGA"
        label="label"
        hit_index=5
        region_slice=15
        overall_length=25
        unaligned_seq="ACCCGACTADMYRNACCTCGATTCKM"
        primer_len=5
        self.assertRaises(TypeError,append_primer_hit)
        self.assertRaises(AttributeError,append_primer_hit,bad_primer,label,
         hit_index,region_slice,overall_length,unaligned_seq,primer_len)
        
    def test_append_primer_hit_appends_data(self):    
        """ append_primer_hit appends valid data correctly to primer object """
        test_primer=ProspectivePrimer("ATCGA",15459,700)
        label="label"
        hit_index=5
        region_slice=20
        overall_length=25
        unaligned_seq="ACCCGATCGAMYRNACCTCGATTCKG"
        primer_len=5
        expected_primer_labels=["label"]
        expected_primer_match_count=1
        expected_primer_percent_match=0
        expected_primer_aligned_index=15459
        expected_primer_numeric_seq=[1,2,4,8,1]
        expected_primer_seq="ATCGA"
        expected_primer_upstream_regions=["---------------ACCCGATCGA"]
        expected_primer_downstream_regions=["ATCGAMYRNACCTCGATTCKG----"]
        append_primer_hit(test_primer,label,hit_index,region_slice,
         overall_length,unaligned_seq,primer_len)
        self.assertEqual(test_primer.labels,expected_primer_labels)
        self.assertEqual(test_primer.match_count,expected_primer_match_count)
        self.assertEqual(test_primer.percent_match,
         expected_primer_percent_match)
        self.assertEqual(test_primer.aligned_index,
         expected_primer_aligned_index)
        self.assertEqual(test_primer.numeric_seq,expected_primer_numeric_seq)
        self.assertEqual(test_primer.seq,expected_primer_seq)
        self.assertEqual(test_primer.upstream_regions,
         expected_primer_upstream_regions)
        self.assertEqual(test_primer.downstream_regions,
         expected_primer_downstream_regions)
    
    def test_append_primer_hit_no_filler_chars(self):    
        """ Test again with small region, should get no "N" filler characters"""
        test_primer=ProspectivePrimer("ATCGA",15459,700)
        label="label"
        hit_index=5
        region_slice=3
        overall_length=8
        unaligned_seq="ACCCGATCGAMYRNACCTCGATTCKG"
        primer_len=5
        expected_primer_labels=["label"]
        expected_primer_match_count=1
        expected_primer_percent_match=0
        expected_primer_aligned_index=15459
        expected_primer_numeric_seq=[1,2,4,8,1]
        expected_primer_seq="ATCGA"
        expected_primer_upstream_regions=["CCGATCGA"]
        expected_primer_downstream_regions=["ATCGAMYR"]
        append_primer_hit(test_primer,label,hit_index,region_slice,
         overall_length,unaligned_seq,primer_len)
        self.assertEqual(test_primer.labels,expected_primer_labels)
        self.assertEqual(test_primer.match_count,expected_primer_match_count)
        self.assertEqual(test_primer.percent_match,
         expected_primer_percent_match)
        self.assertEqual(test_primer.aligned_index,
         expected_primer_aligned_index)
        self.assertEqual(test_primer.numeric_seq,expected_primer_numeric_seq)
        self.assertEqual(test_primer.seq,expected_primer_seq)
        self.assertEqual(test_primer.upstream_regions,
         expected_primer_upstream_regions)
        self.assertEqual(test_primer.downstream_regions,
         expected_primer_downstream_regions)
        
    def test_del_primers(self):
        """ del_primers properly deletes primer objects from list """
        primers=[ProspectivePrimer("ATTA",500,56),
         ProspectivePrimer("CGGT",400,50), ProspectivePrimer("ATAT",350,45)]
        primers_to_delete=[0,2]
        # Should only retain the ProspectivePrimer("CGGT",400) after deletion
        del_primers(primers,primers_to_delete)
        expected_primers_len=1
        expected_primers_seq="CGGT"
        expected_primers_aligned_index=400
        self.assertEqual(len(primers), expected_primers_len)
        self.assertEqual(primers[0].seq,expected_primers_seq)
        self.assertEqual(primers[0].aligned_index,
         expected_primers_aligned_index)
        
    def test_get_number_seqs_for_primer(self):
        """ get_number_seqs_for_primer gets proper number seqs for sens value"""
        percent_match=0.50
        seq_count=10
        expected_result=5
        self.assertEqual(get_number_seqs_for_primer(percent_match,seq_count),
         expected_result)
        percent_match=0.60
        seq_count=100
        expected_result=40
        self.assertEqual(get_number_seqs_for_primer(percent_match,seq_count),
         expected_result)
        
    def test_build_seq_data(self):
        """ build_seq_data properly initializes primers dictionary """
        primers={}
        bad_seq="AA-T--TCXJ"
        seq_length=5 
        search_range = None
        self.assertRaises(AlphabetError,build_seq_data,bad_seq,seq_length,
         primers, search_range)
        """ build_seq_data builds list of primer objects with valid data """
        single_seq="aA--C--GT--A--C---"
        seq_length=5
        
        # Should return list with 3 primer objects
        
        expected_primers={('CGTAC', 4): 2, ('AACGT', 0): 0, ('ACGTA', 1): 1}
        primers=build_seq_data(single_seq,seq_length,primers, search_range)
        self.assertEqual(primers,expected_primers)
        
    def test_build_seq_data_search_range(self):
        """ Test search range parameter, should narrow range of primers built"""
        
        single_seq="AA--C--GT--A--C---"
        search_range = "1:15"
        primers={}
        seq_length=5
        expected_primers={('CGTAC', 4): 2, ('ACGTA', 1): 1}
        primers=build_seq_data(single_seq,seq_length,primers, search_range)
        self.assertEqual(primers,expected_primers)
        
        
    def test_convert_to_numeric(self):
        """ convert_to_numeric gets correct decimal value for bitwise compare"""
        bad_primer="AATXOB"
        self.assertRaises(KeyError,convert_to_numeric,bad_primer)
        """ convert_to_numeric converts IUPAC characters to numeric list.
        These lists are used for bitwise comparisons """
        sample_seq="ATUCGNRYMKWSBDHV"
        expected_list=[1,2,2,4,8,15,9,6,5,10,3,12,14,11,7,13]
        self.assertEqual(convert_to_numeric(sample_seq),expected_list)

    def test_get_corrected_index(self):
        """ get_corrected_index returns correct unaligned index """
        
        self.assertRaises(TypeError,get_corrected_index)
        self.assertRaises(TypeError,get_corrected_index,10)
        self.assertRaises(TypeError,get_corrected_index,"string","10")
        self.assertRaises(TypeError,get_corrected_index,10,10)
        """ get_corrected_index returns corrected unaligned index """
        seq="----A---CTG--AA----GGC---"
        aligned_index=10
        expected_unaligned_index=3
        self.assertEqual(get_corrected_index(seq,aligned_index),
         expected_unaligned_index)

    def test_find_specific_primer_matches_raises_errors(self):
        """ correctly identifies sequences at specificity level """
        
        unaligned_seq=""
        integer_mapped_seq=[2,1,8,1,1,4,8,2,1,1,2]
        deletion_threshold=0
        seq_count=1
        sequence_length=5
        label="test sequence"
        region_slice=2
        seq="--TA--G-AA--C--GT--A--AT---"
        primers=[]
        self.assertRaises(ValueError,find_specific_primer_matches,primers,
         integer_mapped_seq,deletion_threshold,seq_count,
         sequence_length,label,unaligned_seq,region_slice,seq)
        
    def test_equals(self):
        """ == is overwritten as expected """
        p1 = ProspectivePrimer("AACGT", 8, 3)
        p2 = ProspectivePrimer("AACGT", 8, 3)
        p3 = ProspectivePrimer("TACGT", 8, 3)
        p4 = ProspectivePrimer("TACGT", 9, 4)
        p5 = ProspectivePrimer("TACGT", 9, 3)
        
        self.assertEqual(p1, p2,
         "Equal primers do not evaluate as equal with ==")
        self.assertNotEqual(p1, p3)
        self.assertNotEqual(p3, p4)
        self.assertNotEqual(p4, p5)
        self.assertNotEqual(p3, p5)
        
        
    def test_find_specific_primer_matches_finds_primers(self):
        """find_primer_matches returns list of primer objects """

        primers=[ProspectivePrimer("AACGT",8,3),
         ProspectivePrimer("ACGTA",9,4), ProspectivePrimer("CGTAC",12,5)]
        # Sequence before conversion to numberic:"TAGAACGTAAT"
        # Should have perfect match to "AACGT" and "ACGTA" but not
        # "CGTAC".  Should delete "AACGT" and "ACGTA", but not "CGTAC"
        integer_mapped_seq=[2,1,8,1,1,4,8,2,1,1,2]
        deletion_threshold=0
        seq_count=1
        sequence_length=5
        label="test sequence"
        unaligned_seq="TAGAACGTAAT"
        # length of included up/downstream BPs from primer
        region_slice=2
        seq="--TA--G-AA--C--GT--A--AT---"
        expected_data=[ProspectivePrimer("CGTAC",12,5)]
        expected_data[0].non_specific_hits=0
        primers=find_specific_primer_matches(primers,integer_mapped_seq,
         deletion_threshold,seq_count,sequence_length,
         label,unaligned_seq,region_slice,seq)
        self.assertEqual(primers, expected_data,
         "Equal primers do not evaluate as equal with ==")
        
        
    def test_find_specific_primer_matches_finds_and_deletes(self):
        """ Should find a perfect match to all 3 sequences, and delete all """
        primers=[ProspectivePrimer("AACGT",8,6),
         ProspectivePrimer("ACGTA",9,8), ProspectivePrimer("CGTAC",12,9)]
        # Sequence before conversion to numberic:"TAGAACGTAAT"
        # Should find a perfect match to all 3 sequences, and delete all 3
        integer_mapped_seq=[2,1,8,1,1,4,8,2,1,4,2]
        deletion_threshold=0
        seq_count=1
        sequence_length=5
        label="test sequence"
        unaligned_seq="TAGAACGTACT"
        # length of included up/downstream BPs from primer
        region_slice=2
        seq="--TA--G-AA--C--GT--A--CT---"
        expected_data=[]
        primers=find_specific_primer_matches(primers,integer_mapped_seq,
         deletion_threshold,seq_count,sequence_length,
         label,unaligned_seq,region_slice,seq)
        self.assertEqual(primers,expected_data,
         "Equal primers do not evaluate as equal with ==")

        
    def test_find_specific_primer_matches_different_alignment(self):
        """ Does not delete seqs that match but have different alignment """
        primers=[ProspectivePrimer("AACGT",9,6),
         ProspectivePrimer("ACGTA",10,8),ProspectivePrimer("CGTAC",13,9)]
        # Even though exact matches are in the sequence, should not delete 
        # sequences if the alignment is off from the exact match.
        integer_mapped_seq=[2,1,8,1,1,4,8,2,1,4,2]
        deletion_threshold=0
        seq_count=1
        # should only look at bases at the exact aligned index.
        sequence_length=5
        label="test sequence"
        unaligned_seq="TAGAACGTACT"
        # length of included up/downstream BPs from primer
        region_slice=2
        seq="--TA--G-AA--C--GT--A--CT---"
        expected_data=[ProspectivePrimer("AACGT",9,6),
         ProspectivePrimer("ACGTA",10,8), ProspectivePrimer("CGTAC",13,9)]
        primers=find_specific_primer_matches(primers,integer_mapped_seq,
         deletion_threshold,seq_count,sequence_length,
         label,unaligned_seq,region_slice,seq)
        # Should not delete any primers, should have a total of 3
        self.assertEqual(len(primers),3)
        # all primers should be equal to expected primers
        for n in range(len(primers)):
            self.assertEqual(primers[n], expected_data[n],
             "Equal primers do not evaluate as equal with ==")

            

    def test_find_sensitive_primer_matches_raises_error(self):
        """ correctly identifies sequences at specificity level """
        
        unaligned_seq=""
        integer_mapped_seq=[2,1,8,1,1,4,8,2,1,1,2]
        deletion_threshold=0
        seq_count=1
        sequence_length=5
        label="test sequence"
        region_slice=2
        seq="--TA--G-AA--C--GT--A--AT---"
        primers=[]
        self.assertRaises(ValueError,find_sensitive_primer_matches,primers,
         integer_mapped_seq,deletion_threshold,seq_count,
         sequence_length,label,unaligned_seq,region_slice,seq)
        
    def test_find_sensitive_primer_matches_finds_primers(self):
        """ find_primer_matches returns list of primer objects found """

        primers=[ProspectivePrimer("AACGT",8,6),
         ProspectivePrimer("ACGTA",9,8),ProspectivePrimer("CGTAC",12,9)]
        # Sequence before conversion to numberic:"TAGAACGTAAT"
        # Should have perfect match to "AACGT" and "ACGTA" but not
        # "CGTAC", which should be deleted
        integer_mapped_seq=[2,1,8,1,1,4,8,2,1,1,2]
        deletion_threshold=0
        seq_count=1
        sequence_length=5
        label="test sequence"
        unaligned_seq="TAGAACGTAAT"
        # length of included up/downstream BPs from primer
        region_slice=2
        seq="--TA--G-AA--C--GT--A--AT---"
        expected_data=[ProspectivePrimer("AACGT",8,6),
         ProspectivePrimer("ACGTA",9,8)]
        expected_data[0].labels.append("test")
        expected_data[0].upstream_regions.append('AGAACGT')
        expected_data[0].downstream_regions.append('AACGTAA')
        expected_data[0].match_count=1
        expected_data[1].labels.append("test")
        expected_data[1].upstream_regions.append('GAACGTA')
        expected_data[1].downstream_regions.append('ACGTAAT')
        expected_data[1].match_count=1
        primers=find_sensitive_primer_matches(primers,integer_mapped_seq,
         deletion_threshold,seq_count,sequence_length,
         label,unaligned_seq,region_slice,seq)
        
        for n in range(len(primers)):
            self.assertEqual(primers[n], expected_data[n],
             "Equal primers do not evaluate as equal with ==")
            
           
    def test_find_sensitive_primer_matches_exceeds_end(self):
        """ Test with region slice that exceeds end of sequence """
        # Should append "N" for these 'unknown' base pairs
        primers=[ProspectivePrimer("AACGT",8,7),
         ProspectivePrimer("ACGTA",9,6),ProspectivePrimer("CGTAC",12,10)]
        # Sequence before conversion to numberic:"TAGAACGTAAT"
        # Should have perfect match to "AACGT" and "ACGTA" but not
        # "CGTAC", which should be deleted
        integer_mapped_seq=[2,1,8,1,1,4,8,2,1,1,2]
        deletion_threshold=0
        seq_count=1
        sequence_length=5
        label="test"
        unaligned_seq="TAGAACGTAAT"
        # length of included up/downstream BPs from primer
        region_slice=5
        seq="--TA--G-AA--C--GT--A--AT---"
        expected_data=[ProspectivePrimer("AACGT",8,7),
         ProspectivePrimer("ACGTA",9,6)]
        expected_data[0].labels.append("test")
        expected_data[0].upstream_regions.append('--TAGAACGT')
        expected_data[0].downstream_regions.append('AACGTAAT--')
        expected_data[0].match_count=1
        expected_data[1].labels.append("test")
        expected_data[1].upstream_regions.append('-TAGAACGTA')
        expected_data[1].downstream_regions.append('ACGTAAT---')
        expected_data[1].match_count=1
        primers=find_sensitive_primer_matches(primers,integer_mapped_seq,
         deletion_threshold,seq_count,sequence_length,
         label,unaligned_seq,region_slice,seq)
        for n in range(len(primers)):
            self.assertEqual(primers[n], expected_data[n],
             "Equal primers do not evaluate as equal with ==")

        
    def test_find_sensitive_primer_matches_only_one_hit(self):
        """ Test to ensure that only one "hit" per dictionary key """
        primers=[ProspectivePrimer("AACGT",15,8)]
        integer_mapped_seq=[1,1,4,8,2,2,1,8,1,1,4,8,2,1,1,2]
        deletion_threshold=0
        seq_count=1
        sequence_length=5
        label="test_sequence"
        unaligned_seq="AACGTTAGAACGTAAT"
        region_slice=5
        seq="--AA-C-GTTA--G-AA--C--GT--A--AT---"
        primers=find_sensitive_primer_matches(primers,integer_mapped_seq,
         deletion_threshold,seq_count,sequence_length,
         label,unaligned_seq,region_slice,seq)
        expected_data=[ProspectivePrimer("AACGT",15,8)]
        expected_data[0].labels.append("test_sequence")
        expected_data[0].upstream_regions.append('GTTAGAACGT')
        expected_data[0].downstream_regions.append('AACGTAAT--')
        expected_data[0].match_count=1
        
        for n in range(len(primers)):
            self.assertEqual(primers[n], expected_data[n],
             "Equal primers do not evaluate as equal with ==")

    def test_get_deletion_threshold(self):
        """ get_deletion_threshold gets correct seq number for del threshold """
        
        self.assertRaises(TypeError,get_deletion_threshold)
        self.assertRaises(TypeError,get_deletion_threshold,"A")
        self.assertRaises(TypeError,get_deletion_threshold,10)
        self.assertRaises(TypeError,get_deletion_threshold,10,"A")
        """ get_deletion_treshold returns expected result """
        # In order to have >90% matches, the max number of sequences
        # that can be missed for 100 sequences is 9
        expected_result=9
        total_seqs=100
        percent_match=0.90
        self.assertEqual(get_deletion_threshold(percent_match,
         total_seqs),expected_result)
        expected_result=40
        total_seqs=100
        percent_match=0.60
        self.assertEqual(get_deletion_threshold(percent_match,
         total_seqs),expected_result)

    def test_calculate_percent_match(self):
        """ calculates correct target,exclude percent match """
        primers=[ProspectivePrimer("AACGT",8,6),
         ProspectivePrimer("ACGTA",9,8),ProspectivePrimer("CGTAC",12,9)]
        primers[0].non_specific_hits=10
        primers[0].match_count=5
        primers[1].non_specific_hits=50
        primers[1].match_count=10
        primers[2].non_specific_hits=25
        primers[2].match_count=15
        target_seqs=100
        non_target_seqs=50
        expected_data=[ProspectivePrimer("AACGT",8,6),
         ProspectivePrimer("ACGTA",9,8), ProspectivePrimer("CGTAC",12,9)]
        expected_data[0].non_specific_percent=0.20
        expected_data[0].percent_match=0.05
        expected_data[1].non_specific_percent=1.0
        expected_data[1].percent_match=0.10
        expected_data[2].non_specific_percent=0.50
        expected_data[2].percent_match=0.15
        calculate_percent_match(primers,target_seqs,non_target_seqs)
        for n in range(len(primers)):
            self.assertEqual(primers[n].seq,expected_data[n].seq)
            self.assertEqual(primers[n].non_specific_percent,
             expected_data[n].non_specific_percent)
            self.assertEqual(primers[n].percent_match,
             expected_data[n].percent_match)
             
class OverallInputOutputTests(TestCase):
    """ Tests overall input/output functionality of module """
    
    def setUp(self):
        # create the temporary input files that will be used with the
        # search_sequences function
        self.test_aligned_fasta_data = test_aligned_fasta_data
        
        self.target_seqs_filepath = get_tmp_filename(\
         prefix='target_seqs',
         suffix='.fasta')
        seq_file = open(self.target_seqs_filepath,'w')
        seq_file.write(self.test_aligned_fasta_data)
        seq_file.close()        
        
        self.test_aligned_specificity_data = test_aligned_specificity_data
        
        self.specificity_seqs_filepath = get_tmp_filename(\
         prefix='specificity_data_',
         suffix='.fasta')
        seq_file = open(self.specificity_seqs_filepath,'w')
        seq_file.write(self.test_aligned_specificity_data)
        seq_file.close()
        
        self.test_std_alignment_data = test_std_alignment_data
        
        self.standard_aln_filepath = get_tmp_filename(\
         prefix='std_alignment_seq_',
         suffix='.fasta')
        seq_file = open(self.standard_aln_filepath,'w')
        seq_file.write(self.test_std_alignment_data)
        seq_file.close()
        
        self.expected_output_with_search_range =\
         expected_output_with_search_range
        self.expected_output_with_std_aligned_and_len_10 =\
         expected_output_with_std_aligned_and_len_10
        self.expected_output_xmer_len_6 = expected_output_xmer_len_6
        self.expected_output_primer_len_10 = expected_output_primer_len_10
        self.expected_output_with_spec_file_sixty_spec =\
         expected_output_with_spec_file_sixty_spec
        self.expected_output_with_spec_file_default =\
         expected_output_with_spec_file_default
        self.expected_output_eighty_percent_sens =\
         expected_output_eighty_percent_sens
        self.expected_output_default_settings =\
         expected_output_default_settings
        self.expected_log_file_default_settings =\
         expected_log_file_default_settings
        
        self._files_to_remove =\
         [self.target_seqs_filepath, self.specificity_seqs_filepath, 
         self.standard_aln_filepath]
        
    def tearDown(self):
        remove_files(self._files_to_remove)
    
    def test_search_sequences_default_settings(self):
        """ Tests main function of generate_primers_denovo, default settings """
        
        default_settings_output_filepath = \
         get_tmp_filename(prefix='default_output_', suffix='.txt')
         
        self._files_to_remove.append(default_settings_output_filepath)
         
         
        search_sequences(self.target_seqs_filepath, sequence_length = 5,
         exclude_fasta_filepath = None, verbose=False, percent_match = 0.60,
         full_primer_length = 20, output_f = default_settings_output_filepath,
         specificity_threshold = 0.01, log_filepath = None,
         standard_index_file = None, search_range = None)
         
       
        output_f = open(default_settings_output_filepath, "U")
        
        default_settings_result = [l.strip('\n') for l in output_f]
            
        self.assertEqual(default_settings_result,
         self.expected_output_default_settings)
         
    def test_search_sequences_log_file(self):
        """ Tests log file output """
        
        default_settings_output_filepath = \
         get_tmp_filename(prefix='default_output_', suffix='.txt')
         
        log_output_filepath = \
         get_tmp_filename(prefix='default_output_', suffix='.log')
         
        self._files_to_remove.append(default_settings_output_filepath)
        self._files_to_remove.append(log_output_filepath)
         
         
        search_sequences(self.target_seqs_filepath, sequence_length = 5,
         exclude_fasta_filepath = None, verbose=False, percent_match = 0.60,
         full_primer_length = 20, output_f = default_settings_output_filepath,
         specificity_threshold = 0.01, log_filepath = log_output_filepath,
         standard_index_file = None, search_range = None)
         
       
        log_output_f = open(log_output_filepath, "U")
        
        default_settings_result = [l.strip('\n') for l in log_output_f]
            
        self.assertEqual(default_settings_result,
         self.expected_log_file_default_settings)
         
    def test_search_sequences_high_sens(self):
        """ Should find nothing at 80% sensitivity setting """
        
        eighty_perc_sens_settings_output_filepath = \
         get_tmp_filename(prefix='eighty_perc_sens_output_', suffix='.txt')
         
        self._files_to_remove.append(eighty_perc_sens_settings_output_filepath)
         
         
        search_sequences(self.target_seqs_filepath, sequence_length = 5,
         exclude_fasta_filepath = None, verbose=False, percent_match = 0.80,
         full_primer_length = 20, 
         output_f = eighty_perc_sens_settings_output_filepath,
         specificity_threshold = 0.01, log_filepath = None,
         standard_index_file = None, search_range = None)
         
       
        output_f = open(eighty_perc_sens_settings_output_filepath, "U")
        
        eighty_perc_sens_settings_result = [l.strip('\n') for l in output_f]
            
        self.assertEqual(eighty_perc_sens_settings_result,
         self.expected_output_eighty_percent_sens)
         
    def test_search_sequences_spec_file_default(self):
        """ Tests main function of generate_primers_denovo 
        
        This test includes a specifity file with the default specificity 
        threshold of 1%, which should result in no conserved Xmers being
        retained."""
        
        spec_file_default_output_filepath = \
         get_tmp_filename(prefix='default_output_spec_', suffix='.txt')
         
        self._files_to_remove.append(spec_file_default_output_filepath)
         
         
        search_sequences(self.target_seqs_filepath, sequence_length = 5,
         exclude_fasta_filepath = self.specificity_seqs_filepath, verbose=False,
         percent_match = 0.60, full_primer_length = 20, 
         output_f = spec_file_default_output_filepath,
         specificity_threshold = 0.01, log_filepath = None,
         standard_index_file = None, search_range = None)
         
       
        output_f = open(spec_file_default_output_filepath, "U")
        
        spec_default_settings_result = [l.strip('\n') for l in output_f]
            
        self.assertEqual(spec_default_settings_result,
         self.expected_output_with_spec_file_default)
         
    def test_search_sequences_spec_file_sixty_percent(self):
        """ Tests main function of generate_primers_denovo 
        
        This test reduces the specificity threshold from the default 1% to
        a much more lenient 60% setting, allowing retention of conserved
        Xmers."""
        
        spec_file_sixty_perc_output_filepath = \
         get_tmp_filename(prefix='default_output_spec_', suffix='.txt')
         
        self._files_to_remove.append(spec_file_sixty_perc_output_filepath)
         
         
        search_sequences(self.target_seqs_filepath, sequence_length = 5,
         exclude_fasta_filepath = self.specificity_seqs_filepath, verbose=False,
         percent_match = 0.60, full_primer_length = 20, 
         output_f = spec_file_sixty_perc_output_filepath,
         specificity_threshold = 0.60, log_filepath = None,
         standard_index_file = None, search_range = None)
         
       
        output_f = open(spec_file_sixty_perc_output_filepath, "U")
        
        spec_sixty_settings_result = [l.strip('\n') for l in output_f]
            
        self.assertEqual(spec_sixty_settings_result,
         self.expected_output_with_spec_file_sixty_spec)
         
         
    def test_search_sequences_primer_len_ten(self):
        """ Tests main function of generate_primers_denovo 
        
        This test reduces the overall primer size to 10 from 20, so only 5
        bases are sliced out up and downstream instead of 15 from the default 
        setting."""
        
        primer_len10_output_filepath = \
         get_tmp_filename(prefix='default_output_', suffix='.txt')
         
        self._files_to_remove.append(primer_len10_output_filepath)
         
         
        search_sequences(self.target_seqs_filepath, sequence_length = 5,
         exclude_fasta_filepath = None, verbose=False, percent_match = 0.60,
         full_primer_length = 10, output_f = primer_len10_output_filepath,
         specificity_threshold = 0.01, log_filepath = None,
         standard_index_file = None, search_range = None)
         
       
        output_f = open(primer_len10_output_filepath, "U")
        
        primer_len10_result = [l.strip('\n') for l in output_f]
            
        self.assertEqual(primer_len10_result,
         self.expected_output_primer_len_10)
         
    def test_search_sequences_xmer_len_6(self):
        """ Tests main function of generate_primers_denovo 
        
        This test increases the size of the conserved Xmer from 5 to 6."""
        
        len6_xmer_output_filepath = \
         get_tmp_filename(prefix='default_output_', suffix='.txt')
         
        self._files_to_remove.append(len6_xmer_output_filepath)
         
         
        search_sequences(self.target_seqs_filepath, sequence_length = 6,
         exclude_fasta_filepath = None, verbose=False, percent_match = 0.60,
         full_primer_length = 20, output_f = len6_xmer_output_filepath,
         specificity_threshold = 0.01, log_filepath = None,
         standard_index_file = None, search_range = None)
         
       
        output_f = open(len6_xmer_output_filepath, "U")
        
        len6_xmer_result = [l.strip('\n') for l in output_f]
            
        self.assertEqual(len6_xmer_result,
         self.expected_output_xmer_len_6)
         
    def test_search_sequences_std_alignment(self):
        """ Tests main function of generate_primers_denovo 
        
        This test includes a standard alignment, which should result in 
        forward and reverse indexes being added to the output file."""
        
        std_alignment_output_filepath = \
         get_tmp_filename(prefix='std_alignment_', suffix='.txt')
         
        self._files_to_remove.append(std_alignment_output_filepath)
         
         
        search_sequences(self.target_seqs_filepath, sequence_length = 5,
         exclude_fasta_filepath = None, verbose=False, percent_match = 0.60,
         full_primer_length = 10, output_f = std_alignment_output_filepath,
         specificity_threshold = 0.01, log_filepath = None,
         standard_index_file = self.standard_aln_filepath, search_range = None)
         
       
        output_f = open(std_alignment_output_filepath, "U")
        
        std_alignment_result = [l.strip('\n') for l in output_f]
            
        self.assertEqual(std_alignment_result,
         self.expected_output_with_std_aligned_and_len_10)
         
    def test_search_sequences_search_range(self):
        """ Tests main function of generate_primers_denovo 
        
        This test limits the range for which conserved 5mers will be found,
        only one conserved site should be found."""
        
        search_range_output_filepath = \
         get_tmp_filename(prefix='search_range_', suffix='.txt')
         
        self._files_to_remove.append(search_range_output_filepath)
         
         
        search_sequences(self.target_seqs_filepath, sequence_length = 5,
         exclude_fasta_filepath = None, verbose=False, percent_match = 0.60,
         full_primer_length = 20, output_f = search_range_output_filepath,
         specificity_threshold = 0.01, log_filepath = None,
         standard_index_file = None, search_range = "5:15")
         
       
        output_f = open(search_range_output_filepath, "U")
        
        search_range_result = [l.strip('\n') for l in output_f]
            
        self.assertEqual(search_range_result,
         self.expected_output_with_search_range)
        


# Start input file data
test_aligned_fasta_data = """>Seq1
--------tTCT-CTTT-----TA--T-GC-T-G-G
>Seq2
--------ATCT-GTTT-----GA--T-CC-T-G-G
>Seq3
--------ATCT-GTTT-----GA--T-CC-T-G-G
>Seq4
--------TTCT-CTTT-----GA--T-CC-T-G-G"""

test_aligned_specificity_data = """>Seq1
--------TTCT-CTTT-----TA--T-GC-T-G-G
>Seq2
--------ATCT-GTTT-----GA--T-CC-T-G-G"""

# Need to use primer_len 10 setting due to artificially short sequence
# Has 5 extra bases inserted at beginning of the sequence
test_std_alignment_data = """>std_alignment_seq
--ATCGA-TTCT-CTTT-----TA--T-GC-T-G-G"""

# Start expected output file data
expected_log_file_default_settings = """Building prosective primers
Constructing primer objects
Counting sequences for excluded fasta file(s)
Total number of target sequences: 4
Finding sensitive primer regions.
Module complete""".split('\n')

expected_output_default_settings = """# Matches that are found in at least 60.00% target sequences
# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches
C,ATCCT,9,23,3,75.0%,0.0%
M,23,------ATCTGTTTGATCCT,ATCCTGG-------------,Seq2
M,23,------ATCTGTTTGATCCT,ATCCTGG-------------,Seq3
M,23,------TTCTCTTTGATCCT,ATCCTGG-------------,Seq4
C,TTTGA,5,14,3,75.0%,0.0%
M,14,----------ATCTGTTTGA,TTTGATCCTGG---------,Seq2
M,14,----------ATCTGTTTGA,TTTGATCCTGG---------,Seq3
M,14,----------TTCTCTTTGA,TTTGATCCTGG---------,Seq4
C,GATCC,8,22,3,75.0%,0.0%
M,22,-------ATCTGTTTGATCC,GATCCTGG------------,Seq2
M,22,-------ATCTGTTTGATCC,GATCCTGG------------,Seq3
M,22,-------TTCTCTTTGATCC,GATCCTGG------------,Seq4
C,CCTGG,11,28,3,75.0%,0.0%
M,28,----ATCTGTTTGATCCTGG,CCTGG---------------,Seq2
M,28,----ATCTGTTTGATCCTGG,CCTGG---------------,Seq3
M,28,----TTCTCTTTGATCCTGG,CCTGG---------------,Seq4
C,TCCTG,10,26,3,75.0%,0.0%
M,26,-----ATCTGTTTGATCCTG,TCCTGG--------------,Seq2
M,26,-----ATCTGTTTGATCCTG,TCCTGG--------------,Seq3
M,26,-----TTCTCTTTGATCCTG,TCCTGG--------------,Seq4
C,TTGAT,6,15,3,75.0%,0.0%
M,15,---------ATCTGTTTGAT,TTGATCCTGG----------,Seq2
M,15,---------ATCTGTTTGAT,TTGATCCTGG----------,Seq3
M,15,---------TTCTCTTTGAT,TTGATCCTGG----------,Seq4
C,TGATC,7,16,3,75.0%,0.0%
M,16,--------ATCTGTTTGATC,TGATCCTGG-----------,Seq2
M,16,--------ATCTGTTTGATC,TGATCCTGG-----------,Seq3
M,16,--------TTCTCTTTGATC,TGATCCTGG-----------,Seq4""".split('\n')

expected_output_eighty_percent_sens = """# Matches that are found in at least 80.00% target sequences
# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches""".split('\n')

expected_output_with_spec_file_default = """# Matches that are found in at least 60.00% target sequences and specific at or below the 1.00% threshold 
# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches""".split('\n')

expected_output_with_spec_file_sixty_spec = """# Matches that are found in at least 60.00% target sequences and specific at or below the 60.00% threshold 
# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches
C,ATCCT,9,23,3,75.0%,50.0%
M,23,------ATCTGTTTGATCCT,ATCCTGG-------------,Seq2
M,23,------ATCTGTTTGATCCT,ATCCTGG-------------,Seq3
M,23,------TTCTCTTTGATCCT,ATCCTGG-------------,Seq4
C,TTTGA,5,14,3,75.0%,50.0%
M,14,----------ATCTGTTTGA,TTTGATCCTGG---------,Seq2
M,14,----------ATCTGTTTGA,TTTGATCCTGG---------,Seq3
M,14,----------TTCTCTTTGA,TTTGATCCTGG---------,Seq4
C,GATCC,8,22,3,75.0%,50.0%
M,22,-------ATCTGTTTGATCC,GATCCTGG------------,Seq2
M,22,-------ATCTGTTTGATCC,GATCCTGG------------,Seq3
M,22,-------TTCTCTTTGATCC,GATCCTGG------------,Seq4
C,CCTGG,11,28,3,75.0%,50.0%
M,28,----ATCTGTTTGATCCTGG,CCTGG---------------,Seq2
M,28,----ATCTGTTTGATCCTGG,CCTGG---------------,Seq3
M,28,----TTCTCTTTGATCCTGG,CCTGG---------------,Seq4
C,TCCTG,10,26,3,75.0%,50.0%
M,26,-----ATCTGTTTGATCCTG,TCCTGG--------------,Seq2
M,26,-----ATCTGTTTGATCCTG,TCCTGG--------------,Seq3
M,26,-----TTCTCTTTGATCCTG,TCCTGG--------------,Seq4
C,TTGAT,6,15,3,75.0%,50.0%
M,15,---------ATCTGTTTGAT,TTGATCCTGG----------,Seq2
M,15,---------ATCTGTTTGAT,TTGATCCTGG----------,Seq3
M,15,---------TTCTCTTTGAT,TTGATCCTGG----------,Seq4
C,TGATC,7,16,3,75.0%,50.0%
M,16,--------ATCTGTTTGATC,TGATCCTGG-----------,Seq2
M,16,--------ATCTGTTTGATC,TGATCCTGG-----------,Seq3
M,16,--------TTCTCTTTGATC,TGATCCTGG-----------,Seq4""".split('\n')

expected_output_primer_len_10 = """# Matches that are found in at least 60.00% target sequences
# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches
C,ATCCT,9,23,3,75.0%,0.0%
M,23,GTTTGATCCT,ATCCTGG---,Seq2
M,23,GTTTGATCCT,ATCCTGG---,Seq3
M,23,CTTTGATCCT,ATCCTGG---,Seq4
C,TTTGA,5,14,3,75.0%,0.0%
M,14,ATCTGTTTGA,TTTGATCCTG,Seq2
M,14,ATCTGTTTGA,TTTGATCCTG,Seq3
M,14,TTCTCTTTGA,TTTGATCCTG,Seq4
C,GATCC,8,22,3,75.0%,0.0%
M,22,TGTTTGATCC,GATCCTGG--,Seq2
M,22,TGTTTGATCC,GATCCTGG--,Seq3
M,22,TCTTTGATCC,GATCCTGG--,Seq4
C,CCTGG,11,28,3,75.0%,0.0%
M,28,TTGATCCTGG,CCTGG-----,Seq2
M,28,TTGATCCTGG,CCTGG-----,Seq3
M,28,TTGATCCTGG,CCTGG-----,Seq4
C,TCCTG,10,26,3,75.0%,0.0%
M,26,TTTGATCCTG,TCCTGG----,Seq2
M,26,TTTGATCCTG,TCCTGG----,Seq3
M,26,TTTGATCCTG,TCCTGG----,Seq4
C,TTGAT,6,15,3,75.0%,0.0%
M,15,TCTGTTTGAT,TTGATCCTGG,Seq2
M,15,TCTGTTTGAT,TTGATCCTGG,Seq3
M,15,TCTCTTTGAT,TTGATCCTGG,Seq4
C,TGATC,7,16,3,75.0%,0.0%
M,16,CTGTTTGATC,TGATCCTGG-,Seq2
M,16,CTGTTTGATC,TGATCCTGG-,Seq3
M,16,CTCTTTGATC,TGATCCTGG-,Seq4""".split('\n')

expected_output_xmer_len_6 = """# Matches that are found in at least 60.00% target sequences
# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches
C,ATCCTG,9,23,3,75.0%,0.0%
M,23,-----ATCTGTTTGATCCTG,ATCCTGG-------------,Seq2
M,23,-----ATCTGTTTGATCCTG,ATCCTGG-------------,Seq3
M,23,-----TTCTCTTTGATCCTG,ATCCTGG-------------,Seq4
C,TTGATC,6,15,3,75.0%,0.0%
M,15,--------ATCTGTTTGATC,TTGATCCTGG----------,Seq2
M,15,--------ATCTGTTTGATC,TTGATCCTGG----------,Seq3
M,15,--------TTCTCTTTGATC,TTGATCCTGG----------,Seq4
C,TGATCC,7,16,3,75.0%,0.0%
M,16,-------ATCTGTTTGATCC,TGATCCTGG-----------,Seq2
M,16,-------ATCTGTTTGATCC,TGATCCTGG-----------,Seq3
M,16,-------TTCTCTTTGATCC,TGATCCTGG-----------,Seq4
C,TTTGAT,5,14,3,75.0%,0.0%
M,14,---------ATCTGTTTGAT,TTTGATCCTGG---------,Seq2
M,14,---------ATCTGTTTGAT,TTTGATCCTGG---------,Seq3
M,14,---------TTCTCTTTGAT,TTTGATCCTGG---------,Seq4
C,TCCTGG,10,26,3,75.0%,0.0%
M,26,----ATCTGTTTGATCCTGG,TCCTGG--------------,Seq2
M,26,----ATCTGTTTGATCCTGG,TCCTGG--------------,Seq3
M,26,----TTCTCTTTGATCCTGG,TCCTGG--------------,Seq4
C,GATCCT,8,22,3,75.0%,0.0%
M,22,------ATCTGTTTGATCCT,GATCCTGG------------,Seq2
M,22,------ATCTGTTTGATCCT,GATCCTGG------------,Seq3
M,22,------TTCTCTTTGATCCT,GATCCTGG------------,Seq4""".split('\n')

expected_output_with_std_aligned_and_len_10 = """# Matches that are found in at least 60.00% target sequences
# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches,forward primer standard index,reverse primer standard index
C,ATCCT,9,23,3,75.0%,0.0%,9,24
M,23,GTTTGATCCT,ATCCTGG---,Seq2
M,23,GTTTGATCCT,ATCCTGG---,Seq3
M,23,CTTTGATCCT,ATCCTGG---,Seq4
C,TTTGA,5,14,3,75.0%,0.0%,5,20
M,14,ATCTGTTTGA,TTTGATCCTG,Seq2
M,14,ATCTGTTTGA,TTTGATCCTG,Seq3
M,14,TTCTCTTTGA,TTTGATCCTG,Seq4
C,GATCC,8,22,3,75.0%,0.0%,8,23
M,22,TGTTTGATCC,GATCCTGG--,Seq2
M,22,TGTTTGATCC,GATCCTGG--,Seq3
M,22,TCTTTGATCC,GATCCTGG--,Seq4
C,CCTGG,11,28,3,75.0%,0.0%,11,26
M,28,TTGATCCTGG,CCTGG-----,Seq2
M,28,TTGATCCTGG,CCTGG-----,Seq3
M,28,TTGATCCTGG,CCTGG-----,Seq4
C,TCCTG,10,26,3,75.0%,0.0%,10,25
M,26,TTTGATCCTG,TCCTGG----,Seq2
M,26,TTTGATCCTG,TCCTGG----,Seq3
M,26,TTTGATCCTG,TCCTGG----,Seq4
C,TTGAT,6,15,3,75.0%,0.0%,6,21
M,15,TCTGTTTGAT,TTGATCCTGG,Seq2
M,15,TCTGTTTGAT,TTGATCCTGG,Seq3
M,15,TCTCTTTGAT,TTGATCCTGG,Seq4
C,TGATC,7,16,3,75.0%,0.0%,7,22
M,16,CTGTTTGATC,TGATCCTGG-,Seq2
M,16,CTGTTTGATC,TGATCCTGG-,Seq3
M,16,CTCTTTGATC,TGATCCTGG-,Seq4""".split('\n')

# With search range parameter 5:15
expected_output_with_search_range = """# Matches that are found in at least 60.00% target sequences
# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches
C,TTTGA,5,14,3,75.0%,0.0%
M,14,----------ATCTGTTTGA,TTTGATCCTGG---------,Seq2
M,14,----------ATCTGTTTGA,TTTGATCCTGG---------,Seq3
M,14,----------TTCTCTTTGA,TTTGATCCTGG---------,Seq4""".split('\n')
        
if __name__ == "__main__":
    main()
