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
 

from cogent3.util.unit_test import TestCase, main
from cogent3 import LoadSeqs, DNA

from primerprospector.format import write_formatted_primers, \
 write_overlap_primers, write_primers_summary, write_denovo_output_file, \
 generate_denovo_output_file, organize_overlapping_primers, get_amplicon_pairs,\
 write_amplicon_pairs, format_linker_data, write_assigned_taxa,\
 write_accuracy_report, format_base_freq_data
from primerprospector.sort_denovo_primers import Primer
from primerprospector.generate_primers_denovo import ProspectivePrimer
from primerprospector.util import FakeOutFile


class FormatTests(TestCase):
    """ Unit tests for the format.py module """
    def setUp(self):
        """ Set up variables for use in unit tests """
        
        self.f_primer_data = f_primer_data
        self.r_primer_data = r_primer_data
        
        self.expected_primers_formatted = expected_primers_formatted
        
        self.expected_data_overlap = expected_data_overlap
        self.expected_primers_report = expected_primer_report
        self.expected_denovo_output_data = expected_denovo_output_data
        self.expected_output_header_data = expected_output_header_data
        self.expected_linker_data = expected_linker_data
        self.expected_taxa_accuracy_report = expected_taxa_accuracy_report
        self.expected_taxa_assignments = expected_taxa_assignments
        self.formatted_base_freq = formatted_base_freq
        
    def test_format_base_freq_data(self):
        """ Formats base frequency for file output correctly """
        
        base_freqs = [{'A':0.8, 'C':0.1, 'T':0.1, 'G':0},
                      {'A':0.6, 'C':0.4, 'T':0, 'G':0},
                      {'A':1.0, 'C':0, 'T':0, 'G':0}]
        primer_seq = "AMA"
        hits_fp = "unit_test_hits"
        score_type = "weighted_score"
        score_threshold = 1.0
        
        actual_results = format_base_freq_data(base_freqs, primer_seq,
         hits_fp, score_type, score_threshold)
         
        expected_result = self.formatted_base_freq
        
        self.assertEqual(actual_results, expected_result)
    
    def test_write_formatted_primers(self):
        """ writes output in proper format """
        formatted_primers_f = FakeOutFile()
        
        write_formatted_primers(self.f_primer_data, self.r_primer_data, 
         formatted_primers_f)
        self.assertEqual(formatted_primers_f.data, 
         self.expected_primers_formatted)
        
    def test_write_overlap_primers(self):
        """ Writes details about overlap w/ known primers in proper format"""
        
        overlap_primers_f = FakeOutFile()
        
        small_seqs_f=LoadSeqs(data=['>Seq1', 'ACCCGATC', '>Seq2', 'ACCCGGGG',
         '>Seq3', 'ACCTGGGG', '>Seq4', 'ACCCGGGG', '>Seq5', 'ACCTGGGG'], 
         moltype=DNA)
        small_seqs_r=LoadSeqs(data=['>Seq1', 'GATCTGGT', '>Seq2', 'GGGGTGGT',
         '>Seq3', 'GGGGTGGT', '>Seq4', 'GGGGTGGT', '>Seq5', 'GGGGTGGT'], 
         moltype=DNA)
        
        sample_primer = [ Primer("GATC,51,2041,9,90.0%,10.0%,102,137", 
         small_seqs_f, small_seqs_r)]
        sample_primer[0].f_unique = True
        sample_primer[0].f_partial_overlap = []
        sample_primer[0].f_3prime_match = []
        sample_primer[0].r_unique = False
        sample_primer[0].r_partial_overlap = \
         ['138r\tGGGYAKRTTRHCYACGTGTT\ttest_known_primer_r\tGTTAACTACGTGTTA\tGTTAACTACG\n']
        sample_primer[0].r_3prime_match = \
         ['137r\tGGYAKRTTRHCYACGTGTTA\ttest_known_primer_r\tGTTAACTACGTGTTA\n']
         
        write_overlap_primers(sample_primer, overlap_primers_f)
        
        self.assertEqual(overlap_primers_f.data, self.expected_data_overlap)
        
    def test_organize_overlapping_primers(self):
        """ Organizes overlapping primer data """
        
        small_seqs_f=LoadSeqs(data=['>Seq1', 'ACCCGATC', '>Seq2', 'ACCCGGGG',
         '>Seq3', 'ACCTGGGG', '>Seq4', 'ACCCGGGG', '>Seq5', 'ACCTGGGG'], 
         moltype=DNA)
        small_seqs_r=LoadSeqs(data=['>Seq1', 'GATCTGGT', '>Seq2', 'GGGGTGGT',
         '>Seq3', 'GGGGTGGT', '>Seq4', 'GGGGTGGT', '>Seq5', 'GGGGTGGT'], 
         moltype=DNA)
        
        sample_primer = [ Primer("GATC,51,2041,9,90.0%,10.0%,102,137", 
         small_seqs_f, small_seqs_r)]
        sample_primer[0].f_unique = True
        sample_primer[0].f_partial_overlap = []
        sample_primer[0].f_3prime_match = []
        sample_primer[0].r_unique = False
        sample_primer[0].r_partial_overlap = \
         ['138r\tGGGYAKRTTRHCYACGTGTT\ttest_known_primer_r\tGTTAACTACGTGTTA\tGTTAACTACG\n']
        sample_primer[0].r_3prime_match = \
         ['137r\tGGYAKRTTRHCYACGTGTTA\ttest_known_primer_r\tGTTAACTACGTGTTA\n']
         
        unique_primers, overlapping_primers, three_prime_match =\
         organize_overlapping_primers(sample_primer)
         
        expected_unique_primers = ['102f\tACCYGRKS\n'] 
        expected_overlapping_primers = ['138r\tGGGYAKRTTRHCYACGTGTT\ttest_known_primer_r\tGTTAACTACGTGTTA\tGTTAACTACG\n'] 
        expected_three_prime_match = ['137r\tGGYAKRTTRHCYACGTGTTA\ttest_known_primer_r\tGTTAACTACGTGTTA\n']
        
        self.assertEqual(unique_primers, expected_unique_primers)
        self.assertEqual(overlapping_primers, expected_overlapping_primers)
        self.assertEqual(three_prime_match, expected_three_prime_match)
        
    def test_write_primers_summary(self):
        """ Writes summary of primer data in proper format """
        
        primer_report_f = FakeOutFile()
        
        small_seqs_f=LoadSeqs(data=['>Seq1', 'ACCCGATC', '>Seq2', 'ACCCGGGG',
         '>Seq3', 'ACCTGGGG', '>Seq4', 'ACCCGGGG', '>Seq5', 'ACCTGGGG'], 
         moltype=DNA)
        small_seqs_r=LoadSeqs(data=['>Seq1', 'GATCTGGT', '>Seq2', 'GGGGTGGT',
         '>Seq3', 'GGGGTGGT', '>Seq4', 'GGGGTGGT', '>Seq5', 'GGGGTGGT'], 
         moltype=DNA)
        
        sample_primer = [ Primer("GATC,51,2041,9,90.0%,10.0%,102,137", 
         small_seqs_f, small_seqs_r)]
        sample_primer[0].f_unique = True
        sample_primer[0].f_partial_overlap = []
        sample_primer[0].f_3prime_match = []
        sample_primer[0].r_unique = False
        sample_primer[0].r_partial_overlap = \
         ['138r\tGGGYAKRTTRHCYACGTGTT\ttest_known_primer_r\tGTTAACTACGTGTTA\tGTTAACTACG\n']
        sample_primer[0].r_3prime_match = \
         ['137r\tGGYAKRTTRHCYACGTGTTA\ttest_known_primer_r\tGTTAACTACGTGTTA\n']
         
        test_variable_pos_freq = 0.20
        
        write_primers_summary(sample_primer, primer_report_f, \
         test_variable_pos_freq)
         
        self.assertEqual(primer_report_f.data, self.expected_primers_report)
        
    def test_write_denovo_output_file(self):
        """ Writes details about Xmers and upstream/downstream sequences """
        
        denovo_output_f = FakeOutFile()
        test_region_slice = 25
        
        primer_data=[ProspectivePrimer("AACGT",8,7),
         ProspectivePrimer("ACGTA",9,6)]
        primer_data[0].labels.append("test")
        primer_data[0].upstream_regions.append('--TAGAACGT')
        primer_data[0].downstream_regions.append('AACGTAAT--')
        primer_data[0].match_count=1
        primer_data[1].labels.append("test")
        primer_data[1].upstream_regions.append('-TAGAACGTA')
        primer_data[1].downstream_regions.append('ACGTAAT---')
        primer_data[1].match_count=1
        
        

        
        write_denovo_output_file(primer_data, denovo_output_f, 
         test_region_slice)
         
        self.assertEqual(denovo_output_f.data, self.expected_denovo_output_data)
        
        
    def test_generate_denovo_output_file(self):
        """ Organizes details about de novo primer data correctly """
        
        
        test_primers = []
        output_filepath = FakeOutFile()
        test_specificity_threshold = 0.10
        test_region_slice = 25
        test_standard_index_seq = None
        test_percent_match = 0.75
        test_specificity_test = True
        
        generate_denovo_output_file(test_primers, output_filepath, 
         test_specificity_threshold, test_region_slice, test_standard_index_seq,
         test_percent_match, test_specificity_test)
         
        self.assertEqual(output_filepath.data, self.expected_output_header_data)
        
    def test_get_amplicon_pairs(self):
        """ Organizes output data for denovo primer pairs within amp_len """
        
        amplicon_len = "10:100"
        f_primers_sorted = ["25f\tAACCGATC","48f\tACRRCGTN"]
        r_primers_sorted = ["77r\tACGGACTT","200r\tRVCCTG"]
        
        # Should find 2 potential amplicons in the 10-100 bp range
        expected_amplicon_data =\
         ['# Primer pairs within estimated amplicon size\n', 
         '# Min size: 10\n', '# Max size: 100\n', '\n', '25f\tAACCGATC', 
         '77r\tACGGACTT', 'Estimated amplicon size: 36\n', '\n', 
         '48f\tACRRCGTN', '77r\tACGGACTT', 'Estimated amplicon size: 13\n']
        
        actual_amplicon_data = get_amplicon_pairs(amplicon_len,
         f_primers_sorted, r_primers_sorted)
        self.assertEqual(actual_amplicon_data, expected_amplicon_data)
        
        
        
    def test_write_amplicon_pairs(self):
        """ Properly writes amplicon restricted primer pair data """
        
        output_filepath = FakeOutFile()
        
        amplicon_data = ['# Primer pairs within estimated amplicon size\n', 
         '# Min size: 10\n', '# Max size: 100\n', '\n', '25f\tAACCGATC', 
         '77r\tACGGACTT', 'Estimated amplicon size: 52\n', '\n', 
         '48f\tACRRCGTN', '77r\tACGGACTT', 'Estimated amplicon size: 29\n']
         
        # Should write out the lines to the file object
        
        write_amplicon_pairs(amplicon_data, output_filepath)
        
        self.assertEqual(output_filepath.data, "".join(amplicon_data))
        
    def test_format_linker_data(self):
        """ Properly formats linker data for output """
        
        
        linker_summary = """Base position 1\nA: 0.033 T: 0.000 C: 0.967 G: 0.000\nBase position 2\nA: 0.000 T: 0.767 C: 0.133 G: 0.100\n"""
        best_seq = "TA"
        worst_seq = "CT"
        hits = "/home/user/primers/forward_primer_hits.txt"
        score_type = "weighted_score"
        score_threshold = 4
        
        actual_linker_data = format_linker_data(linker_summary, best_seq,
         worst_seq, hits, score_type, score_threshold)
         
        self.assertEqual(actual_linker_data, self.expected_linker_data)
        
            
    def test_write_accuracy_report(self):
        """ Writes accuracy reports correctly """
        
        
        
        accuracy_values = [(0.90,50),(0.85,45),(0.66,35)]
        seqs_lacking_taxa_mapping = 1
        seqs_assigned_root = 3
        report_f = FakeOutFile()
        fasta_fp = "/fasta_source/255f_forward_reads.fasta"
        assignment_method = "rdp"
        training_data_fp = None
        
        write_accuracy_report(accuracy_values, seqs_lacking_taxa_mapping,
         seqs_assigned_root, report_f, fasta_fp, assignment_method,
         training_data_fp)
         
        
        self.assertEqual(report_f.data, self.expected_taxa_accuracy_report)
        
        
        
    def test_write_assigned_taxa(self):
        """ Writes taxonomic assignments correctly """
        
        assigned_taxa_f = FakeOutFile()
        
                
        taxa_assignments = {"seq1":("Root;Archaea;Euryarchaeota;Halobacteria",0.8800),
                            "seq2":("Root;Bacteria;Firmicutes",0.9900),
                            "seq3":("Root;Eukarya;Metazoa",0.75)}
                            
        write_assigned_taxa(taxa_assignments, assigned_taxa_f)
        
        
        # Because data is stored in a dict, order of writing can be random
        # Sort to results to test for equivalency
        
        actual_results = []
        
        for line in assigned_taxa_f.data.split("\n"):
            # Avoid appending empty list element due to last split on \n
            if len(line):
                actual_results.append(line)

        actual_results.sort()
            
        
        self.assertEqual(actual_results, self.expected_taxa_assignments)
        




        
        

f_primer_data = ['104f\tGCRDACGGCTCAGTAACACG\n', '103f\tGGCRDACGGCTCAGTAACAC\n', '102f\tHGGCRDACGGCTCAGTAACA\n', '101f\tKHGGCRDACGGCTCAGTAAC\n', '105f\tCRWACGGCTCAGTAACACGT\n', '99f\tVCKHGGCRDACGGCTCAGTA\n', '98f\tDVCKHGGCRDACGGCTCAGT\n', '97f\tRDVCKHGGCRDACGGCTCAG\n']
r_primer_data = ['139r\tWGGGYAKRTTRHCYACGTGT\n', '138r\tGGGYAKRTTRHCYACGTGTT\n', '137r\tGGYAKRTTRHCYACGTGTTA\n', '136r\tGYAKRTTRHCYACGTGTTAC\n', '140r\tHWGGGYAKRTTRHCYACGTG\n', '134r\tAKRTTRMCYACGTGTTACTG\n', '133r\tKRTTRMCYACGTGTTACTGA\n', '132r\tRTTRMCYACGTGTTACTGAG\n']

expected_primers_formatted = """# primer_id <tab> primer sequence (5'->3')\n104f\tGCRDACGGCTCAGTAACACG\n103f\tGGCRDACGGCTCAGTAACAC\n102f\tHGGCRDACGGCTCAGTAACA\n101f\tKHGGCRDACGGCTCAGTAAC\n105f\tCRWACGGCTCAGTAACACGT\n99f\tVCKHGGCRDACGGCTCAGTA\n98f\tDVCKHGGCRDACGGCTCAGT\n97f\tRDVCKHGGCRDACGGCTCAG\n139r\tWGGGYAKRTTRHCYACGTGT\n138r\tGGGYAKRTTRHCYACGTGTT\n137r\tGGYAKRTTRHCYACGTGTTA\n136r\tGYAKRTTRHCYACGTGTTAC\n140r\tHWGGGYAKRTTRHCYACGTG\n134r\tAKRTTRMCYACGTGTTACTG\n133r\tKRTTRMCYACGTGTTACTGA\n132r\tRTTRMCYACGTGTTACTGAG\n"""

expected_data_overlap = """#Unique primers that do not overlap or share 3' ends with known primers
#primer name<tab>primer sequence
102f\tACCYGRKS
#Primers overlapping with known primers.
#primer name<tab>primer sequence<tab>known primer name<tab>known primer sequence<tab>overlapping sequence
138r\tGGGYAKRTTRHCYACGTGTT\ttest_known_primer_r\tGTTAACTACGTGTTA\tGTTAACTACG
#Primers sharing 3' ends with known primers.
#primer name<tab>primer sequence<tab>known primer name<tab>known primer sequence
137r\tGGYAKRTTRHCYACGTGTTA\ttest_known_primer_r\tGTTAACTACGTGTTA
"""

expected_primer_report = """# Primer Report
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
GATC,51,2041,9,90.0%,10.0%,102,137,ACCYGRKS,ACCYGRKS,ACCCGGGG,ACCASMYC,ACCASMYC,ACCACCCC
"""

expected_denovo_output_data = """C,AACGT,7,8,1,0%,0%
M,8,--TAGAACGT,AACGTAAT--,test
C,ACGTA,6,9,1,0%,0%
M,9,-TAGAACGTA,ACGTAAT---,test
"""
expected_output_header_data = """# Matches that are found in at least 75.00% target sequences and specific at or below the 10.00% threshold 
# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches
"""

expected_linker_data = """# Summary data for suggested linkers
# Note-degenerate bases are ignored for linker results, as are linkers that
# exceed the length of the input fasta sequence.  Base position starts at the
# 5' position of the linker, so position 0 would be the 5' most base in the
# linker, position 1 would be the next base, and so on.
Hits file used to generate linkers: forward_primer_hits.txt
Score type used: weighted_score
Score threshold: 4

Base position 1
A: 0.033 T: 0.000 C: 0.967 G: 0.000
Base position 2
A: 0.000 T: 0.767 C: 0.133 G: 0.100

Suggested linker sequence:
TA
Worst linker sequence:
CT"""

expected_taxa_assignments = ['seq1\tRoot;Archaea;Euryarchaeota;Halobacteria\t0.880',
'seq2\tRoot;Bacteria;Firmicutes\t0.990',
'seq3\tRoot;Eukarya;Metazoa\t0.750']

expected_taxa_accuracy_report = """# Taxonomy Assignment Report File
# This report starts at the highest level of taxa (i.e., Domain level)
# and lists the accuracy of the assignment and the number of sequences
# that had a taxonomic assignment and taxonomic mapping to that depth
# Fasta file used for taxonomic assignments: 255f_forward_reads.fasta
# Assignment method: rdp
# Start report accuracy data
# Taxa level, percent accurate assignment, number of sequences with taxa defined at this level
0,0.900,50
1,0.850,45
2,0.660,35
Sequences lacking corresponding ID in taxonomy mapping file: 1
Sequences assigned only as 'Root': 3
"""

formatted_base_freq = """# Base frequency report for optimizing primers
# This file is tab separated for easy importation into Excel or other spreadsheets
# Degenerate DNA codes (listed here for convenience): R=AG, Y=CT, M=AC, K=GT, W=AT, S=CG, B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT
# The primer is listed in 5' to 3' orientation.
# Hits file used to generate base frequency data: unit_test_hits
# Score type used: weighted_score
# Score threshold: 1
# Primer sequence: AMA
#
Primer\tA\tM\tA
Base
A\t0.80\t0.60\t1.00
T\t0.10\t0.00\t0.00
C\t0.10\t0.40\t0.00
G\t0.00\t0.00\t0.00
"""


if __name__ == "__main__":
    main()
