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
from cogent3.util.misc import remove_files, create_dir, get_random_directory_name

from primerprospector.taxa_assignment_report import generate_taxa_report,\
 get_accuracy_report, get_report_output_fp, get_taxa_assignment_output_fp,\
 write_accuracy_report, write_assigned_taxa, parse_lineage, build_tree,\
 generate_training_seqs, generate_training_files




class TaxaCoverageTests(TestCase):
    """unit tests for taxa_coverage """
    
    def setUp(self):
        """ Generate temporary files for use in unit testing"""
        
        self.output_d = get_random_directory_name(prefix = '/tmp/')
        self.output_d += '/'
        
        create_dir(self.output_d)
        
        self.expected_report_header = expected_report_header
        self.expected_report_header_retraining = \
         expected_report_header_retraining
        self.taxonomy_mapping_file = taxonomy_mapping_file
        self.training_file = training_file
        self.fasta_file = fasta_file
        self.RDP_taxonomy_mapping_file = RDP_taxonomy_mapping_file
        self.expected_RDP_formatted_taxa = \
         expected_RDP_formatted_taxa
        self.rdp_expected_training_seqs = rdp_expected_training_seqs
        self.expected_training_seqs = expected_training_seqs
        
        self.mapping_fp = self.output_d + "taxa_mapping.txt"
        taxa_map = open(self.mapping_fp, "w")
        taxa_map.write(self.taxonomy_mapping_file)
        taxa_map.close()
        
        self.RDP_taxonomy_mapping_fp = self.output_d + "_RDP_taxa_mapping.txt"
        taxa_map = open(self.RDP_taxonomy_mapping_fp, "w")
        taxa_map.write(self.RDP_taxonomy_mapping_file)
        taxa_map.close()
        
        self.fasta_fp = self.output_d + "seqs.fasta"
        fasta_f = open(self.fasta_fp, "w")
        fasta_f.write(self.fasta_file)
        fasta_f.close()
        
        self.training_fp = self.output_d + "_training_seqs.fasta"
        fasta_f = open(self.training_fp, "w")
        fasta_f.write(self.training_file)
        fasta_f.close()
        
        

        self._files_to_remove = [self.mapping_fp, self.fasta_fp]
        
    def tearDown(self):
        remove_files(self._files_to_remove)
        if exists(self.output_d ):
            rmtree(self.output_d )
        
            
    def test_generate_taxa_report(self):
        """ Overall test for taxa assignment report """
        
        taxa_mapping_fp = self.mapping_fp
        fasta_fp = self.fasta_fp
        taxa_depth = 4
        output_dir = self.output_d
        assignment_method = 'rdp'
        min_confidence = 0.80
        training_data_fp = None
        

        
        generate_taxa_report(taxa_mapping_fp, fasta_fp, taxa_depth, output_dir,
         assignment_method, min_confidence, training_data_fp)
         
        actual_report_file =\
         open(self.output_d + 'seqs_accuracy_report.txt', 'U')
         
        # Can only test for proper header since RDP gives somewhat random result
        actual_report = [line.strip() for line in actual_report_file][0:8]
        
        self.assertEqual(actual_report, self.expected_report_header)
        
        actual_taxa_assignment_file =\
         open(self.output_d + 'seqs_assignments.txt', 'U')
         
        # RDP gives random results, so only testing for identify of input
        # sequence labels in output file.
        
        expected_seqs_labels = ['seq1', 'seq2', 'seq3', 'seq4']
        
        actual_seq_labels =\
         [line.split('\t')[0] for line in actual_taxa_assignment_file]
         
        actual_seq_labels.sort()
        
        self.assertEqual(actual_seq_labels, expected_seqs_labels)
        
    def test_generate_taxa_report_retraining(self):
        """ Overall test for taxa assignment report with retraining """
        
        taxa_mapping_fp = self.RDP_taxonomy_mapping_fp
        fasta_fp = self.fasta_fp
        taxa_depth = 4
        output_dir = self.output_d
        assignment_method = 'rdp'
        min_confidence = 0.80
        training_data_fp = self.training_fp
        

        
        generate_taxa_report(taxa_mapping_fp, fasta_fp, taxa_depth, output_dir,
         assignment_method, min_confidence, training_data_fp)
         
        actual_report_file =\
         open(self.output_d + 'seqs_accuracy_report.txt', 'U')
         
        # Can only test for proper header since RDP gives somewhat random result
        actual_report = [line.strip() for line in actual_report_file][0:8]
        
        self.assertEqual(actual_report, self.expected_report_header_retraining)
        
        actual_taxa_assignment_file =\
         open(self.output_d + 'seqs_assignments.txt', 'U')
         
        # RDP gives random results, so only testing for identify of input
        # sequence labels in output file.
        
        expected_seqs_labels = ['seq1', 'seq2', 'seq3', 'seq4']
        
        actual_seq_labels =\
         [line.split('\t')[0] for line in actual_taxa_assignment_file]
         
        actual_seq_labels.sort()
        
        self.assertEqual(actual_seq_labels, expected_seqs_labels)
        
    def test_get_accuracy_report(self):
        """ Properly assigns accuracy percentages for each level of taxa """
        
        
        
        taxa_map = {'seq1':'Archaea;Crenarchaeota;Pyrobaculum',
                    'seq2':'Bacteria;Firmicutes;Clostridia',
                    'seq3':'Bacteria;Proteobacteria;Enterococcus',
                    'seq4':'Eukarya;Metazoa;Hemichordata'}
                    
        assigned_taxa =\
         {'seq1':('Root;Archaea;Euryarchaeota;Halobacteriales', 0.8800),
          'seq2':('Root;Bacteria;Firmicutes;Clostridia;Clostridiales', 0.950),
          'seq3':('Root', 0.50),
          'seq5':('Root;Bacteria;Spirocheates;Treponema', 0.755)}
        taxa_depth = 4
        
        accuracy_values, seqs_lacking_taxa_mapping, seqs_assigned_root =\
         get_accuracy_report(taxa_map, assigned_taxa, taxa_depth)
         
        expected_accuracy_values = [(100.0, 2), (50.0, 2), (50.0, 2)]
        expected_lacking_mapping = 1
        expected_assigned_root = 1
        
        self.assertEqual(accuracy_values, expected_accuracy_values)
        self.assertEqual(seqs_lacking_taxa_mapping, expected_lacking_mapping)
        self.assertEqual(seqs_assigned_root, expected_assigned_root)
        
        
        
    def test_get_report_output_fp(self):
        """ Gets reports output filepath correctly """
        
        output_dir = '../test_output'
        fasta_fp = '/source_fasta/227f_339r_r_reads_150.fasta'
        
        # Not written to, just checks for proper name of output file based on
        # input file analyzed.
        expected_report_output_fp =\
         '../test_output/227f_339r_r_reads_150_accuracy_report.txt'
         
        actual_report_output_fp = get_report_output_fp(output_dir, fasta_fp)
        
        self.assertEqual(actual_report_output_fp, expected_report_output_fp)
        
    def test_get_taxa_assignment_output_fp(self):
        """ Gets taxonomic assignments filepath correctly """
        
        output_dir = '/test_output/'
        fasta_fp = '../227f_339r_r_reads_150.fasta'
        
        expected_assignments_output_fp =\
         '/test_output/227f_339r_r_reads_150_assignments.txt'
         
        actual_assignments_output_fp =\
         get_taxa_assignment_output_fp(output_dir, fasta_fp)
        
        self.assertEqual(actual_assignments_output_fp,
         expected_assignments_output_fp)
         
    def test_parse_lineage(self):
        """Lineage in csv format is correctly parsed to a list
        
        Modified from QIIME assign_taxonomy code by Kyle Bittinger.
        """
        str =\
         'Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.;a;b'
        actual = parse_lineage(str)
        expected = ['Archaea', 'Euryarchaeota', 'Methanomicrobiales', 
                    'Methanomicrobium et rel.', 'a', 'b']
        self.assertEqual(actual, expected)
        
    def test_build_tree(self):
        """build_tree() gets tree with Rdp-format taxonomy
        
        Modified from QIIME assign_taxonomy code by Kyle Bittinger.
        """
        taxa_mapping_fp = open(self.RDP_taxonomy_mapping_fp, "U")
        
        tree = build_tree(taxa_mapping_fp)
        actual = tree.rdp_taxonomy().strip()
        # The order of the lines in this file depends on python's
        # dict() implementation, so we should ideally build two sets
        # of lines and check that their contents match.
        
        expected = self.expected_RDP_formatted_taxa.strip()
        self.assertEqual(actual, expected)
        
    def test_generate_training_seqs(self):
        """ Generates training sequences for RDP retraining properly """
        mapping_file = open(self.RDP_taxonomy_mapping_fp, "U")
        training_data = open(self.training_fp, "U")  
        seqs = generate_training_seqs(training_data, mapping_file)
        
        actual_result = [seq for seq in seqs]
         
        self.assertEqual(actual_result, self.rdp_expected_training_seqs)
        
    def test_generate_training_files(self):
        """ Generates training files for RDP retraining properly """
        
        rdp_taxonomy_file, rdp_training_seqs_file =\
         generate_training_files(self.training_fp,
         self.RDP_taxonomy_mapping_fp)
         
        actual_taxonomy_result =\
         "\n".join([line.strip() for line in rdp_taxonomy_file])
        actual_rdp_training_seqs = [line for line in rdp_training_seqs_file]
        
        self.assertEqual(actual_taxonomy_result,
         self.expected_RDP_formatted_taxa)
        self.assertEqual(actual_rdp_training_seqs, self.expected_training_seqs)


        
        
fasta_file = """>seq1
TACGGTTACCCTTGTTACGACTTAAAGCCCTAACTTTAAATAAGCTGTAATACAGCTTTACTATGAGTTTTATACAGAAGAGTAGAATTTTATATGAAAGGGTGAAATCTGGAAATATATAAAGGAATATCAATAGCGAAGGCAACTCTTTAGTAAAAACTGACGTTGAGGAACGAAAGTGTAGGTAGCAAACAGGATTAGATACCCTGGTAGTCTACACTGTAAATTCTGAGTGCTGTAATTTAATGTAAAAGTTAATTCTTTTATTTAGATTTTTAAGTTAACACCATAAGCACTCCGCCTGAGGAGTATGATCGCAAGGTTGGAACTCAAAGGAATTGACGGGGGCTAGGACTAGTGGTGGAGTATGTGGTTTAATGCGATAATCCGCGCAAAACCTTACCAATTTTTGACATAATAGCTAAAATTTTATTTTACTTTGAGGTAAAATAATATATTTTAACTATTACAGGTGTTGCACGACTGTCGTCAGCTCGTGTTGTGAGATGTTTGGTTAATTCCTTTTAAACGAGCGTAACTCTTATCTTTAATTAAATAATTTACTATTTTTAAAAATTTTTTTTTAAAAGTTTTATTTGAAACTGCTAATAAATAAGTTAGAGGAAGGTGAGAACTAAGCCAAGTCAATGTGACCCTTATAAATTGGGCTACACACGTGCTACAATGGATACTACAGAAAGTTACAAAAATTCAGTTTTTGTTAATCTATAAAAGTCTCCCCAGTTCGGATTGTAGGCTGAAACTCGCTTACATGAAGTTGGAATCACTAGTAATCGCGAATCAGAACGTCGTGGTGAATATGTAAACTAGCCTTGTACACACCGCCCGTCACATGCAAAGAGTTAGCTTTTTTTGAAAGTTTAACTTTGTTAAGCTGATTTAAAAGGTAGTAATTTGGATGAAGTCGTAACAAGGT
>seq2
AATCTGCACACACCTTCTCAGACTTCCTACGGAGCCAGCATTTGCGGATTATTTGGACAATCGGGCCGAACGCCTTGATTCCAGCCATGCGGGTGTGTGAATGATGGCTCTTCGAATTTCTAACACGCACTTTAAGGTAGCGAGGAAAGCGTTGTAGCATTCATTACTCCTCTATTTGATGTTCACCGACAGAACTAAGCACTCGGTTAATCTCTGTGCCCAGCAGCCGCGGTTAATACAGAGGGGTGCAAGTGTTTAATCGGAAATTACTGGGCATAACAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGGCCTCCACGTGGGAACTGCATTCAAAACTGACTGACTAGAGTATGGTAGAGGGTGGTGGAATCTCCTGTGTCGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCAACTAGCCGTTGGGAGCCTTGAGCTCTTAGTGGCGCAGCTAACGCATTAAGTTGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGAAGCAACGCGAAGAACCTTACCAGGCCTTGACATCCAATGAACTTTCTAGAGATAGATTGGTGCCTTCGGGAACATTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGTAACGAGCGCAACCCTTGTCCTTAGTTACCAGCACGTTATGGTGGGCACTCTAAGGAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGGCCTGGGCTACACACGTGCTACAATGGTCGGTACAGAGGGTTGCCAAGCCGCGAGGTGGAGCTAATCCCAGAAAACCGATTGTAGTCCGGATCGCAGTCTGCAACTCGACTGCGTGAAGTCGGAATCGCTAGTAATCGCGAATCAGAATGTCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCACCAGAAGTAGCTAGTCTAACCTTCGGGAGGACGGTATACATCAGTATAATCATCGCTTTGAG
>seq3
GGCAGCTCAAGTGATGGCCAATCTTATTGGGCCTAAAGCGTCCGTAGCTGGCCGTGAAAGTTCGTCGGGAAATCCATCCGCTCAACGGATGGGCGTCCGGCGAAAACTTCACGGCTTGGGACCGGAAGGCTCGAGGGGTACGTCTGGGGTAGGAGTGAAATCCCGTAATCCTGGACGGACCGCCGATGGCGAAAGCACCTCGAGAAGACGGATCCGACAGTGAGGGACGAAAGCTAGGGTCTCGAACCGG
>seq4
CTGCAGAAGATACCGCGGGACTAGGAGGCGGGAGAGGTGGACGGTACTCCACGTGTAGGGGTAAAATCCTTTGATCCGTGGAAGACCACCAGTGGCGAAGGCGGTCCACCAGAACGCG"""

training_file = """>seq2
ACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACTGTGCCATGCAGCACTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCTGTATGAGGCAGGTT
>seq3
TCGCCTGAGATGCATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACCGTACCATGCGGTACTGTGCGCTTATGCGGTATTAGCAGCCGTTTCCAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCG"""

taxonomy_mapping_file = """# Test taxonomy mapping file
seq2\tBacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Flavimonas
seq3\tArchaea;Crenarchaeota;Sulfolobales;Sulfolobacaea
seq4\tArchaea"""

RDP_taxonomy_mapping_file = """seq1\tRoot;Bacteria;Firmicutes;Clostridiales;Clostridia;Clostridium
seq2\tRoot;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae
seq3\tRoot;Archaea;Crenarchaeota;Sulfolobales;Sulfolobacaea;Sulfolobacaea
seq4\tRoot;Bacteria;unclassified;unclassified;unclassified;unclassified"""

expected_RDP_formatted_taxa = """1*Root*0*0*domain
2*Archaea*1*1*phylum
3*Crenarchaeota*2*2*class
4*Sulfolobales*3*3*order
5*Sulfolobacaea*4*4*family
6*Sulfolobacaea*5*5*genus
7*Bacteria*1*1*phylum
8*Firmicutes*7*2*class
9*Clostridiales*8*3*order
10*Clostridia*9*4*family
11*Clostridium*10*5*genus
12*Proteobacteria*7*2*class
13*Gammaproteobacteria*12*3*order
14*Pseudomonadales*13*4*family
15*Pseudomonadaceae*14*5*genus
16*unclassified*7*2*class
17*unclassified*16*3*order
18*unclassified*17*4*family
19*unclassified*18*5*genus"""

expected_report_header = """# Taxonomy Assignment Report File
# This report starts at the highest level of taxa (i.e., Domain level)
# and lists the accuracy of the assignment and the number of sequences
# that had a taxonomic assignment and taxonomic mapping to that depth
# Fasta file used for taxonomic assignments: seqs.fasta
# Assignment method: rdp
# Start report accuracy data
# Taxa level, percent accurate assignment, number of sequences with taxa defined at this level""".split('\n')

expected_report_header_retraining = ['# Taxonomy Assignment Report File', '# This report starts at the highest level of taxa (i.e., Domain level)', '# and lists the accuracy of the assignment and the number of sequences', '# that had a taxonomic assignment and taxonomic mapping to that depth', '# Fasta file used for taxonomic assignments: seqs.fasta', '# Assignment method: rdp', '# Training data filepath for RDP classifier: _training_seqs.fasta', '# Start report accuracy data']

rdp_expected_training_seqs = [('seq2 Root;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae', 'ACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACTGTGCCATGCAGCACTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCTGTATGAGGCAGGTT'), ('seq3 Root;Archaea;Crenarchaeota;Sulfolobales;Sulfolobacaea;Sulfolobacaea', 'TCGCCTGAGATGCATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACCGTACCATGCGGTACTGTGCGCTTATGCGGTATTAGCAGCCGTTTCCAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCG')]

expected_training_seqs = ['>seq2 Root;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae\n', 'ACTGATCGTCGCCTTGGTGAGCCGTTACCTCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACTGTGCCATGCAGCACTGTGCGCTTATGCGGTATTAGCAGTCATTTCTAACTGTTATCCCCTGTATGAGGCAGGTT\n', '>seq3 Root;Archaea;Crenarchaeota;Sulfolobales;Sulfolobacaea;Sulfolobacaea\n', 'TCGCCTGAGATGCATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGGTCACCCTCTCAGGTCGGCTACTGATCGTCGCCTTGGTAGGCCGTTACCCCACCAACTAGCTAATCAGACGCGGGTCCATCTCATACCACCGGAGTTTTTCACACCGTACCATGCGGTACTGTGCGCTTATGCGGTATTAGCAGCCGTTTCCAACTGTTATCCCCCTGTATGAGGCAGGTTACCCACGCGTTACTCACCCGTCCGCCG\n']

if __name__ == "__main__":
    main()

