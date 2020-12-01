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
from os import remove, makedirs
from shutil import rmtree
from glob import glob

from cogent.util.unit_test import TestCase, main
from cogent.app.util import get_tmp_filename
from cogent.util.misc import remove_files, get_random_directory_name

from primerprospector.check_primer_barcode_dimers import expand_degeneracies,\
 get_energy_res, check_barcodes_and_primers
 
class CheckPrimerBarcodeTests(TestCase):
    
    def setUp(self):
        """ Initialize tmp files for unit tests """
        
        
        self.sample_bc_file = get_tmp_filename(prefix = "sample_bc_",
         suffix = ".txt")
        bc_file = open(self.sample_bc_file, 'w')
        bc_file.write(sample_low_energy_barcodes)
        bc_file.close()
        
        self.dna_energies_file = get_tmp_filename(prefix = "dna_energies",
         suffix = ".txt")
        e_file = open(self.dna_energies_file, 'w')
        e_file.write(dna_parameters)
        e_file.close()
        
        
        self.output_dir = get_random_directory_name(prefix = '/tmp/')
        self.output_dir += '/'
        
        self.expected_results_no_flags = expected_results_no_flags
        self.expected_filtered_bcs_no_flags = expected_filtered_bcs_no_flags
        
        self.expected_results_low_temp = expected_results_low_temp
        self.expected_filtered_bcs_low_temp = expected_filtered_bcs_low_temp
        
        self.expected_results_degen = expected_results_degen
        self.expected_filtered_bcs_degen = expected_filtered_bcs_degen
        
        self.expected_results_paired_end = expected_results_paired_end
        self.expected_filtered_bcs_paired_end = expected_filtered_bcs_paired_end
        
         
        if not exists(self.output_dir):
            makedirs(self.output_dir)
        
        self._files_to_remove =\
         [self.sample_bc_file, self.dna_energies_file]
         
        
         
    def tearDown(self):
        remove_files(self._files_to_remove)
        rmtree(self.output_dir)
        
        
    def test_expand_degeneracies(self):
        """ Properly gets all degenerate possibilities for a DNA sequence """
        
        # Basic example of purine/pyramidine degeneracies
        
        sample_seq = "ATACGRATC"
        expected_seqs = ['ATACGAATC', 'ATACGGATC']
        expected_seqs.sort()
        actual_results = expand_degeneracies(sample_seq)
        actual_results.sort()
        self.assertEqual(actual_results, expected_seqs)
        
        # Test other degenerate codes
        
        sample_seq = "N"
        expected_seqs = ['A', 'T', 'C', 'G']
        expected_seqs.sort()
        actual_results = expand_degeneracies(sample_seq)
        actual_results.sort()
        self.assertEqual(actual_results, expected_seqs)
        
        sample_seq = "YW"
        expected_seqs = ['CA', 'CT', 'TA', 'TT']
        expected_seqs.sort()
        actual_results = expand_degeneracies(sample_seq)
        actual_results.sort()
        self.assertEqual(actual_results, expected_seqs)
        
        sample_seq = "BH"
        expected_seqs = ['CA', 'CC', 'CT', 'GA', 'GC', 'GT', 'TA', 'TC', 'TT']
        expected_seqs.sort()
        actual_results = expand_degeneracies(sample_seq)
        actual_results.sort()
        self.assertEqual(actual_results, expected_seqs)
        
    def test_get_energy_res(self):
        """ Pulls energy value from string correctly """
        
        sample_rna_fold_output = "TACGACCATACGAAAAACCCCCGGGGG----------CCCCCGGGGGTTTTT,..((......)).((((((((((((((..........))))))))))))))., (-22.39)"
        
        expected_energy = -22.39
        
        actual_energy = get_energy_res(sample_rna_fold_output)
        
        self.assertEqual(actual_energy, expected_energy)
        
    def test_check_barcodes_and_primers_default(self):
        """ Flags barcodes/seqs correctly for low energy structures """
        
        # Cannot test data in output .ps files, but can test for these files
        # being generated properly
        
        primer1 = "ACAGACAATAGGAATACG"
        primer2 = "CGAACGAGAATAGGACAGAGAT"
        
        check_barcodes_and_primers(self.sample_bc_file, primer1, primer2,
         self.dna_energies_file, output_dir = self.output_dir)
         
        # Should be no .ps files created in output directory
        
        glob_results = glob(self.output_dir + "/*.ps")
        
        self.assertEqual(len(glob_results), 0)
        
        actual_results_f = open(self.output_dir + "/barcode_results.txt", "U")
        
        actual_results = "\n".join([line.strip() for line in actual_results_f])
        
        self.assertEqual(actual_results, self.expected_results_no_flags)
        
        actual_filtered_f = open(self.output_dir + "/filtered_barcodes.txt",
         "U")
         
        actual_results = "\n".join([line.strip() for line in actual_filtered_f])
        
        self.assertEqual(actual_results, self.expected_filtered_bcs_no_flags)
        
    def test_check_barcodes_and_primers_low_temp(self):
        """ Finds more secondary structure at lower temperatures """
        
        # Cannot test data in output .ps files, but can test for these files
        # being generated properly
        
        primer1 = "ACAGACAATAGGAATACG"
        primer2 = "CGAACGAGAATAGGACAGAGAT"
        
        temp = 15
        
        check_barcodes_and_primers(self.sample_bc_file, primer1, primer2,
         self.dna_energies_file, output_dir = self.output_dir,
         annealing_temp = temp)
         
        # Should create 3 .ps files in output directory
        
        glob_results = glob(self.output_dir + "/*.ps")
        
        self.assertEqual(len(glob_results), 3)
        
        actual_results_f = open(self.output_dir + "/barcode_results.txt", "U")
        
        actual_results = "\n".join([line.strip() for line in actual_results_f])
        
        self.assertEqual(actual_results, self.expected_results_low_temp)
        
        actual_filtered_f = open(self.output_dir + "/filtered_barcodes.txt",
         "U")
         
        actual_results = "\n".join([line.strip() for line in actual_filtered_f])
        
        self.assertEqual(actual_results, self.expected_filtered_bcs_low_temp)
        
    def test_check_barcodes_and_primers_handles_degeneracies(self):
        """ Handles/flags degeneracy in primers correctly """
        
        # Cannot test data in output .ps files, but can test for these files
        # being generated properly
        
        primer1 = "ACAGACAARYYGAATACG"
        primer2 = "CGAACGAGAATAGGACAGAGAT"
        
        check_barcodes_and_primers(self.sample_bc_file, primer1, primer2,
         self.dna_energies_file, output_dir = self.output_dir)
         
        # Should create 1 .ps files created in output directory
        
        glob_results = glob(self.output_dir + "/*.ps")
        
        self.assertEqual(len(glob_results), 3)
        
        actual_results_f = open(self.output_dir + "/barcode_results.txt", "U")
        
        actual_results = "\n".join([line.strip() for line in actual_results_f])
        
        self.assertEqual(actual_results, self.expected_results_degen)
        
        actual_filtered_f = open(self.output_dir + "/filtered_barcodes.txt",
         "U")
         
        actual_results = "\n".join([line.strip() for line in actual_filtered_f])
        
        self.assertEqual(actual_results, self.expected_filtered_bcs_degen)
        
        
    def test_check_barcodes_and_primers_paired_end(self):
        """ Correctly tests for paired end barcodes secondary structure """
        
        # Cannot test data in output .ps files, but can test for these files
        # being generated properly
        
        primer1 = "ACAGACAATAGGAATACG"
        # Primer 2 altered to it will pair with barcode 2
        primer2 = "CGAACGAGAACGGCGAGAGAT"
        
        check_barcodes_and_primers(self.sample_bc_file, primer1, primer2,
         self.dna_energies_file, output_dir = self.output_dir,
         paired_end_barcodes = True)
         
        # Should create 2 .ps files created in output directory, with primer1 
        # and primer 2 binding, and primer2 binding itself
        
        glob_results = glob(self.output_dir + "/*.ps")
        
        self.assertEqual(len(glob_results), 2)
        
        actual_results_f = open(self.output_dir + "/barcode_results.txt", "U")
        
        actual_results = "\n".join([line.strip() for line in actual_results_f])
        
        self.assertEqual(actual_results, self.expected_results_paired_end)
        
        actual_filtered_f = open(self.output_dir + "/filtered_barcodes.txt",
         "U")
         
        actual_results = "\n".join([line.strip() for line in actual_filtered_f])
        
        self.assertEqual(actual_results, self.expected_filtered_bcs_paired_end)
        
        
         
        
        
        
# Large strings placed at end for better readability
        
sample_low_energy_barcodes = """# unit test barcode data
CAGCTTGAGATT
ATATCGCCGTTA
"""

expected_results_no_flags = """# Barcode/Primer combinations that fall below the energy threshold.
# primer1 and primer2 will be listed in all non-degenerate forms that result in energy values below threshold.
# If primer2 is not barcoded, it will be tested against itself once, and listed first with the barcode index and sequence fields left empty.
# line number of barcode from input barcode file, barcode sequence, primer1 identity, primer2 identity, combined sequence, secondary structure, Gibbs energy in kcal/mol
No barcodes/primer combination was found below the specified energy threshold."""

expected_filtered_bcs_no_flags = """# The following barcodes were not flagged for any potential secondary structure with the given primers.
CAGCTTGAGATT
ATATCGCCGTTA"""

expected_results_low_temp = """# Barcode/Primer combinations that fall below the energy threshold.
# primer1 and primer2 will be listed in all non-degenerate forms that result in energy values below threshold.
# If primer2 is not barcoded, it will be tested against itself once, and listed first with the barcode index and sequence fields left empty.
# line number of barcode from input barcode file, barcode sequence, primer1 identity, primer2 identity, combined sequence, secondary structure, Gibbs energy in kcal/mol
0,CAGCTTGAGATT,Primer1,Primer2,CAGCTTGAGATTACAGACAATAGGAATACG----------CGAACGAGAATAGGACAGAGAT,....(((..........)))........((..............))................,-2.66
1,ATATCGCCGTTA,Primer1,Primer1,ATATCGCCGTTAACAGACAATAGGAATACG----------ATATCGCCGTTAACAGACAATAGGAATACG,........((((((..................................))))))................,-4.61
1,ATATCGCCGTTA,Primer1,Primer2,ATATCGCCGTTAACAGACAATAGGAATACG----------CGAACGAGAATAGGACAGAGAT,........(((....)))..........((..............))................,-2.72"""

expected_filtered_bcs_low_temp = """# All barcodes were flagged for secondary structure with the given primers."""

expected_results_degen = """# Barcode/Primer combinations that fall below the energy threshold.
# primer1 and primer2 will be listed in all non-degenerate forms that result in energy values below threshold.
# If primer2 is not barcoded, it will be tested against itself once, and listed first with the barcode index and sequence fields left empty.
# line number of barcode from input barcode file, barcode sequence, primer1 identity, primer2 identity, combined sequence, secondary structure, Gibbs energy in kcal/mol
0,CAGCTTGAGATT,Primer1,Primer1,CAGCTTGAGATTACAGACAAGCTGAATACG----------CAGCTTGAGATTACAGACAAGCTGAATACG,(((((((..........)))))))................(((((((..........)))))))......,-7.79
0,CAGCTTGAGATT,Primer1,Primer2,CAGCTTGAGATTACAGACAAGCTGAATACG----------CGAACGAGAATAGGACAGAGAT,(((((((..........)))))))......................................,-3.66
0,CAGCTTGAGATT,Primer1,Primer1,CAGCTTGAGATTACAGACAAGCCGAATACG----------CAGCTTGAGATTACAGACAAGCCGAATACG,..(((((..........)))))..................(.(((((..........))))).)......,-2.55"""

expected_filtered_bcs_degen = """# The following barcodes were not flagged for any potential secondary structure with the given primers.
ATATCGCCGTTA"""

expected_results_paired_end = """# Barcode/Primer combinations that fall below the energy threshold.
# primer1 and primer2 will be listed in all non-degenerate forms that result in energy values below threshold.
# If primer2 is not barcoded, it will be tested against itself once, and listed first with the barcode index and sequence fields left empty.
# line number of barcode from input barcode file, barcode sequence, primer1 identity, primer2 identity, combined sequence, secondary structure, Gibbs energy in kcal/mol
1,ATATCGCCGTTA,Primer1,Primer2,ATATCGCCGTTAACAGACAATAGGAATACG----------ATATCGCCGTTACGAACGAGAACGGCGAGAGAT,............................................(((((((.........)))))))......,-6.36
1,ATATCGCCGTTA,Primer2,Primer2,ATATCGCCGTTACGAACGAGAACGGCGAGAGAT----------ATATCGCCGTTACGAACGAGAACGGCGAGAGAT,....(((((((.........)))))))....................(((((((.........)))))))......,-12.72"""

expected_filtered_bcs_paired_end = """# The following barcodes were not flagged for any potential secondary structure with the given primers.
CAGCTTGAGATT"""

dna_parameters = """## RNAfold parameter file
/* This file contains energy parameters for folding DNA  */
/* The parameters were made available by David Mathews */
/* and are identical to the ones distributed with RNAstructure */

/* WARNING!! This file is still incomplete and largely untested */
/* It does not contain: */
/* terminal mismatch energies for 1xn and 2x3 loops */
/* we currently still use dangling ends instead of terminal mismatches */
/* in multi-loops and the exterior loop. */
/* There's no special parameters for co-axial stacking */

# stack_energies
/*  CG    GC    GU    UG    AU    UA    @  */
  -220  -190   -20   -60  -130  -160    0
  -190  -220   -50     0  -160  -130    0
   -20   -50   130    80   -10    40    0
   -60     0    80    20    70    40    0
  -130  -160   -10    70   -80  -100    0
  -160  -130    40    40  -100   -60    0
     0     0     0     0     0     0    0 


# stack_enthalpies
/*  CG    GC    GU    UG    AU    UA    @  */
  -980  -750  -240  -430  -580  -990    0
  -750  -790  -240   430  -780  -850    0
  -240  -240   500   800  -170   -90    0
  -430   430   800  -670   370   -90    0
  -580  -780  -170   370  -530  -720    0
  -990  -850   -90   -90  -720  -670    0
     0     0     0     0     0     0    0 


# dangle5
/*            N     A     C     G     T */
/* NP */     INF   INF   INF   INF   INF
/* CG */     -50  -100   -50   -60   -80
/* GC */     -40   -60   -40   -60   -60
/* GT */     -40   -60   -40   -60   -60
/* TG */     -50  -100   -50   -60   -80
/* AT */     -40   -40   -40   -60   -80
/* TA */      70   -30    10    70    10
/* NN */      70   -30    10    70    10

# dangle3
/*            N     A     C     G     T */
/* NP */     INF   INF   INF   INF   INF
/* CG */     -30   -90   -30   -50   -40
/* GC */       0  -100   -30     0   -50
/* GT */       0  -100   -30     0   -50
/* TG */     -30   -90   -30   -50   -40
/* AT */     -20   -40   -20   -50   -30
/* TA */      50    10    50    10    20
/* NN */      50    10    50    10    20

# dangle5_enthalpies
/*            N     A     C     G     T */
/* NP */     INF   INF   INF   INF   INF
/* CG */     -90  -620  -330  -510   -90
/* GC */    -420  -460  -420  -470  -500
/* GT */    -420  -460  -420  -470  -500
/* TG */     -90  -620  -330  -510   -90
/* AT */     -60   -90   -60  -250  -900
/* TA */     200    80   200   150   110
/* NN */     200    80   200   150   110

# dangle3_enthalpies
/*            N     A     C     G     T */
/* NP */     INF   INF   INF   INF   INF
/* CG */     -20  -190   -20  -400  -550
/* GC */    -250  -800  -250  -310  -390
/* GT */    -250  -800  -250  -310  -390
/* TG */     -20  -190   -20  -400  -550
/* AT */     200  -530   200  -480    90
/* TA */     630   430   630   250   260
/* NN */     630   430   630   250   260

# int11_energies
/* CG..CG */
/* N */    160    90   160    90    90
/* A */     90    70    90   -20    90
/* C */    160    90   160    90    60
/* G */     90   -20    90    20    90
/* T */     90    90    60    90     0
/* CG..GC */
/* N */    130   100   130   100    80
/* A */     50    50    10  -120    10
/* C */    130   100   130   100    80
/* G */     10   -40    10   -60    10
/* T */    100   100    70   100    70
/* CG..GT */
/* N */    130   100   130   100    80
/* A */     50    50    10  -120    10
/* C */    130   100   130   100    80
/* G */     10   -40    10   -60    10
/* T */    100   100    70   100    70
/* CG..TG */
/* N */    160    90   160    90    90
/* A */     90    70    90   -20    90
/* C */    160    90   160    90    60
/* G */     90   -20    90    20    90
/* T */     90    90    60    90     0
/* CG..AT */
/* N */    180   120   180   120   110
/* A */    120   120   110   -80   110
/* C */    180   120   180   120   -60
/* G */    110   -10   110   -70   110
/* T */    120   120     0   120   -40
/* CG..TA */
/* N */    220   140   220   140   170
/* A */    170   140   170    70   170
/* C */    220   140   220   140   120
/* G */    170    40   170    30   170
/* T */    140   140   120   140    80
/* CG..NN */
/* N */    220   140   220   140   170
/* A */    170   140   170    70   170
/* C */    220   140   220   140   120
/* G */    170    40   170    30   170
/* T */    140   140   120   140    80
/* GC..CG */
/* N */    130    50   130    10   100
/* A */    100    50   100   -40   100
/* C */    130    10   130    10    70
/* G */    100  -120   100   -60   100
/* T */     80    10    80    10    70
/* GC..GC */
/* N */    100    40   100    40    90
/* A */     40    20    40  -150    40
/* C */    100    40   100    40    70
/* G */     40  -150    40  -180    40
/* T */     90    40    70    40    90
/* GC..GT */
/* N */    100    40   100    40    90
/* A */     40    20    40  -150    40
/* C */    100    40   100    40    70
/* G */     40  -150    40  -180    40
/* T */     90    40    70    40    90
/* GC..TG */
/* N */    130    50   130    10   100
/* A */    100    50   100   -40   100
/* C */    130    10   130    10    70
/* G */    100  -120   100   -60   100
/* T */     80    10    80    10    70
/* GC..AT */
/* N */    260   140   260   130   160
/* A */    140   140   120   -20   120
/* C */    260   130   260   130   160
/* G */    120    10   120  -110   120
/* T */    130   130    90   130    90
/* GC..TA */
/* N */    200   120   200   120   170
/* A */    170    70   170    70   170
/* C */    200   120   200   120   130
/* G */    170    40   170    50   170
/* T */    200   120   200   120    40
/* GC..NN */
/* N */    260   140   260   130   170
/* A */    170   140   170    70   170
/* C */    260   130   260   130   160
/* G */    170    40   170    50   170
/* T */    200   130   200   130    90
/* GT..CG */
/* N */    130    50   130    10   100
/* A */    100    50   100   -40   100
/* C */    130    10   130    10    70
/* G */    100  -120   100   -60   100
/* T */     80    10    80    10    70
/* GT..GC */
/* N */    100    40   100    40    90
/* A */     40    20    40  -150    40
/* C */    100    40   100    40    70
/* G */     40  -150    40  -180    40
/* T */     90    40    70    40    90
/* GT..GT */
/* N */    100    40   100    40    90
/* A */     40    20    40  -150    40
/* C */    100    40   100    40    70
/* G */     40  -150    40  -180    40
/* T */     90    40    70    40    90
/* GT..TG */
/* N */    130    50   130    10   100
/* A */    100    50   100   -40   100
/* C */    130    10   130    10    70
/* G */    100  -120   100   -60   100
/* T */     80    10    80    10    70
/* GT..AT */
/* N */    260   140   260   130   160
/* A */    140   140   120   -20   120
/* C */    260   130   260   130   160
/* G */    120    10   120  -110   120
/* T */    130   130    90   130    90
/* GT..TA */
/* N */    200   120   200   120   170
/* A */    170    70   170    70   170
/* C */    200   120   200   120   130
/* G */    170    40   170    50   170
/* T */    200   120   200   120    40
/* GT..NN */
/* N */    260   140   260   130   170
/* A */    170   140   170    70   170
/* C */    260   130   260   130   160
/* G */    170    40   170    50   170
/* T */    200   130   200   130    90
/* TG..CG */
/* N */    160    90   160    90    90
/* A */     90    70    90   -20    90
/* C */    160    90   160    90    60
/* G */     90   -20    90    20    90
/* T */     90    90    60    90     0
/* TG..GC */
/* N */    130   100   130   100    80
/* A */     50    50    10  -120    10
/* C */    130   100   130   100    80
/* G */     10   -40    10   -60    10
/* T */    100   100    70   100    70
/* TG..GT */
/* N */    130   100   130   100    80
/* A */     50    50    10  -120    10
/* C */    130   100   130   100    80
/* G */     10   -40    10   -60    10
/* T */    100   100    70   100    70
/* TG..TG */
/* N */    160    90   160    90    90
/* A */     90    70    90   -20    90
/* C */    160    90   160    90    60
/* G */     90   -20    90    20    90
/* T */     90    90    60    90     0
/* TG..AT */
/* N */    180   120   180   120   110
/* A */    120   120   110   -80   110
/* C */    180   120   180   120   -60
/* G */    110   -10   110   -70   110
/* T */    120   120     0   120   -40
/* TG..TA */
/* N */    220   140   220   140   170
/* A */    170   140   170    70   170
/* C */    220   140   220   140   120
/* G */    170    40   170    30   170
/* T */    140   140   120   140    80
/* TG..NN */
/* N */    220   140   220   140   170
/* A */    170   140   170    70   170
/* C */    220   140   220   140   120
/* G */    170    40   170    30   170
/* T */    140   140   120   140    80
/* AT..CG */
/* N */    180   120   180   110   120
/* A */    120   120   120   -10   120
/* C */    180   110   180   110     0
/* G */    120   -80   120   -70   120
/* T */    110   110   -60   110   -40
/* AT..GC */
/* N */    260   140   260   120   130
/* A */    140   140   130    10   130
/* C */    260   120   260   120    90
/* G */    130   -20   130  -110   130
/* T */    160   120   160   120    90
/* AT..GT */
/* N */    260   140   260   120   130
/* A */    140   140   130    10   130
/* C */    260   120   260   120    90
/* G */    130   -20   130  -110   130
/* T */    160   120   160   120    90
/* AT..TG */
/* N */    180   120   180   110   120
/* A */    120   120   120   -10   120
/* C */    180   110   180   110     0
/* G */    120   -80   120   -70   120
/* T */    110   110   -60   110   -40
/* AT..AT */
/* N */    290   120   290   100   140
/* A */    120   120   100   -40   120
/* C */    290   100   290   100    60
/* G */    100   -40   100   -10   100
/* T */    140   120    60   100   140
/* AT..TA */
/* N */    230   130   230   120   150
/* A */    150   130   150    50   150
/* C */    230   120   230   120   120
/* G */    150     0   150    40   150
/* T */    130   120   130   120    30
/* AT..NN */
/* N */    290   140   290   120   150
/* A */    150   140   150    50   150
/* C */    290   120   290   120   120
/* G */    150     0   150    40   150
/* T */    160   120   160   120   140
/* TA..CG */
/* N */    220   170   220   170   140
/* A */    140   140   140    40   140
/* C */    220   170   220   170   120
/* G */    140    70   140    30   140
/* T */    170   170   120   170    80
/* TA..GC */
/* N */    200   170   200   170   200
/* A */    120    70   120    40   120
/* C */    200   170   200   170   200
/* G */    120    70   120    50   120
/* T */    170   170   130   170    40
/* TA..GT */
/* N */    200   170   200   170   200
/* A */    120    70   120    40   120
/* C */    200   170   200   170   200
/* G */    120    70   120    50   120
/* T */    170   170   130   170    40
/* TA..TG */
/* N */    220   170   220   170   140
/* A */    140   140   140    40   140
/* C */    220   170   220   170   120
/* G */    140    70   140    30   140
/* T */    170   170   120   170    80
/* TA..AT */
/* N */    230   150   230   150   130
/* A */    130   130   120     0   120
/* C */    230   150   230   150   130
/* G */    120    50   120    40   120
/* T */    150   150   120   150    30
/* TA..TA */
/* N */    220   180   220   180   180
/* A */    180   160   180    90   180
/* C */    220   180   220   180   110
/* G */    180    90   180   100   180
/* T */    180   180   110   180   120
/* TA..NN */
/* N */    230   180   230   180   200
/* A */    180   160   180    90   180
/* C */    230   180   230   180   200
/* G */    180    90   180   100   180
/* T */    180   180   130   180   120
/* NN..CG */
/* N */    220   170   220   170   140
/* A */    140   140   140    40   140
/* C */    220   170   220   170   120
/* G */    140    70   140    30   140
/* T */    170   170   120   170    80
/* NN..GC */
/* N */    260   170   260   170   200
/* A */    140   140   130    40   130
/* C */    260   170   260   170   200
/* G */    130    70   130    50   130
/* T */    170   170   160   170    90
/* NN..GT */
/* N */    260   170   260   170   200
/* A */    140   140   130    40   130
/* C */    260   170   260   170   200
/* G */    130    70   130    50   130
/* T */    170   170   160   170    90
/* NN..TG */
/* N */    220   170   220   170   140
/* A */    140   140   140    40   140
/* C */    220   170   220   170   120
/* G */    140    70   140    30   140
/* T */    170   170   120   170    80
/* NN..AT */
/* N */    290   150   290   150   160
/* A */    140   140   120     0   120
/* C */    290   150   290   150   160
/* G */    120    50   120    40   120
/* T */    150   150   120   150   140
/* NN..TA */
/* N */    230   180   230   180   180
/* A */    180   160   180    90   180
/* C */    230   180   230   180   130
/* G */    180    90   180   100   180
/* T */    200   180   200   180   120
/* NN..NN */
/* N */    290   180   290   180   200
/* A */    180   160   180    90   180
/* C */    290   180   290   180   200
/* G */    180    90   180   100   180
/* T */    200   180   200   180   140

# int11_enthalpies
/* CG..CG */
/* N */    610   610   510   610   510
/* A */    610   400   510   610   510
/* C */    510   510   410   510   160
/* G */    610   610   510    10   510
/* T */    510   510   160   510  -270
/* CG..GC */
/* N */    930   930   750   930   460
/* A */    460   -60   460   130   460
/* C */    930   930   750   930   190
/* G */    460   390   460   -90   460
/* T */    930   930    50   930   210
/* CG..GT */
/* N */    930   930   750   930   460
/* A */    460   -60   460   130   460
/* C */    930   930   750   930   190
/* G */    460   390   460   -90   460
/* T */    930   930    50   930   210
/* CG..TG */
/* N */    610   610   510   610   510
/* A */    610   400   510   610   510
/* C */    510   510   410   510   160
/* G */    610   610   510    10   510
/* T */    510   510   160   510  -270
/* CG..AT */
/* N */    530   500   530   390   530
/* A */    530   500   530  -220   530
/* C */    390   390   260   390 -1230
/* G */    530  -630   530  -990   530
/* T */    390   390  -490   390  -750
/* CG..TA */
/* N */   1320  1000  1320  1130  1320
/* A */   1320  1000  1320  1130  1320
/* C */    960   960   590   960    10
/* G */   1320   690  1320  -110  1320
/* T */    960   960   300   960    70
/* CG..NN */
/* N */   1320  1000  1320  1130  1320
/* A */   1320  1000  1320  1130  1320
/* C */    960   960   750   960   190
/* G */   1320   690  1320    10  1320
/* T */    960   960   300   960   210
/* GC..CG */
/* N */    930   460   930   460   930
/* A */    930   -60   930   390   930
/* C */    750   460   750   460    50
/* G */    930   130   930   -90   930
/* T */    460   460   190   460   210
/* GC..GC */
/* N */    600   390   600   390   390
/* A */    390   260   390   240   390
/* C */    600   390   600   390   -60
/* G */    390   240   390  -280   390
/* T */    390   390   -60   390    30
/* GC..GT */
/* N */    600   390   600   390   390
/* A */    390   260   390   240   390
/* C */    600   390   600   390   -60
/* G */    390   240   390  -280   390
/* T */    390   390   -60   390    30
/* GC..TG */
/* N */    930   460   930   460   930
/* A */    930   -60   930   390   930
/* C */    750   460   750   460    50
/* G */    930   130   930   -90   930
/* T */    460   460   190   460   210
/* GC..AT */
/* N */   1280  1090  1280  1090  1160
/* A */   1000  1000   410   450   410
/* C */   1280  1090  1280  1090  1160
/* G */    890   890   410  -500   410
/* T */   1090  1090   130  1090    50
/* GC..TA */
/* N */   1320   690  1320  1130  1320
/* A */   1320  -140  1320  1130  1320
/* C */   1290   420  1290   420   -40
/* G */   1320   690  1320   340  1320
/* T */    420   420   370   420 -1000
/* GC..NN */
/* N */   1320  1090  1320  1130  1320
/* A */   1320  1000  1320  1130  1320
/* C */   1290  1090  1290  1090  1160
/* G */   1320   890  1320   340  1320
/* T */   1090  1090   370  1090   210
/* GT..CG */
/* N */    930   460   930   460   930
/* A */    930   -60   930   390   930
/* C */    750   460   750   460    50
/* G */    930   130   930   -90   930
/* T */    460   460   190   460   210
/* GT..GC */
/* N */    600   390   600   390   390
/* A */    390   260   390   240   390
/* C */    600   390   600   390   -60
/* G */    390   240   390  -280   390
/* T */    390   390   -60   390    30
/* GT..GT */
/* N */    600   390   600   390   390
/* A */    390   260   390   240   390
/* C */    600   390   600   390   -60
/* G */    390   240   390  -280   390
/* T */    390   390   -60   390    30
/* GT..TG */
/* N */    930   460   930   460   930
/* A */    930   -60   930   390   930
/* C */    750   460   750   460    50
/* G */    930   130   930   -90   930
/* T */    460   460   190   460   210
/* GT..AT */
/* N */   1280  1090  1280  1090  1160
/* A */   1000  1000   410   450   410
/* C */   1280  1090  1280  1090  1160
/* G */    890   890   410  -500   410
/* T */   1090  1090   130  1090    50
/* GT..TA */
/* N */   1320   690  1320  1130  1320
/* A */   1320  -140  1320  1130  1320
/* C */   1290   420  1290   420   -40
/* G */   1320   690  1320   340  1320
/* T */    420   420   370   420 -1000
/* GT..NN */
/* N */   1320  1090  1320  1130  1320
/* A */   1320  1000  1320  1130  1320
/* C */   1290  1090  1290  1090  1160
/* G */   1320   890  1320   340  1320
/* T */   1090  1090   370  1090   210
/* TG..CG */
/* N */    610   610   510   610   510
/* A */    610   400   510   610   510
/* C */    510   510   410   510   160
/* G */    610   610   510    10   510
/* T */    510   510   160   510  -270
/* TG..GC */
/* N */    930   930   750   930   460
/* A */    460   -60   460   130   460
/* C */    930   930   750   930   190
/* G */    460   390   460   -90   460
/* T */    930   930    50   930   210
/* TG..GT */
/* N */    930   930   750   930   460
/* A */    460   -60   460   130   460
/* C */    930   930   750   930   190
/* G */    460   390   460   -90   460
/* T */    930   930    50   930   210
/* TG..TG */
/* N */    610   610   510   610   510
/* A */    610   400   510   610   510
/* C */    510   510   410   510   160
/* G */    610   610   510    10   510
/* T */    510   510   160   510  -270
/* TG..AT */
/* N */    530   500   530   390   530
/* A */    530   500   530  -220   530
/* C */    390   390   260   390 -1230
/* G */    530  -630   530  -990   530
/* T */    390   390  -490   390  -750
/* TG..TA */
/* N */   1320  1000  1320  1130  1320
/* A */   1320  1000  1320  1130  1320
/* C */    960   960   590   960    10
/* G */   1320   690  1320  -110  1320
/* T */    960   960   300   960    70
/* TG..NN */
/* N */   1320  1000  1320  1130  1320
/* A */   1320  1000  1320  1130  1320
/* C */    960   960   750   960   190
/* G */   1320   690  1320    10  1320
/* T */    960   960   300   960   210
/* AT..CG */
/* N */    530   530   390   530   390
/* A */    500   500   390  -630   390
/* C */    530   530   260   530  -490
/* G */    390  -220   390  -990   390
/* T */    530   530 -1230   530  -750
/* AT..GC */
/* N */   1280  1000  1280   890  1090
/* A */   1090  1000  1090   890  1090
/* C */   1280   410  1280   410   130
/* G */   1090   450  1090  -500  1090
/* T */   1160   410  1160   410    50
/* AT..GT */
/* N */   1280  1000  1280   890  1090
/* A */   1090  1000  1090   890  1090
/* C */   1280   410  1280   410   130
/* G */   1090   450  1090  -500  1090
/* T */   1160   410  1160   410    50
/* AT..TG */
/* N */    530   530   390   530   390
/* A */    500   500   390  -630   390
/* C */    530   530   260   530  -490
/* G */    390  -220   390  -990   390
/* T */    530   530 -1230   530  -750
/* AT..AT */
/* N */    980   980   980   980   980
/* A */    980   490   980   710   980
/* C */    980   980   750   980   430
/* G */    980   710   980   410   980
/* T */    980   980   430   980   570
/* AT..TA */
/* N */   1470  1290  1470  1070   150
/* A */   1290  1290   150  1070   150
/* C */   1470  1040  1470  1040  -160
/* G */    180   180   150   160   150
/* T */   1040  1040   350  1040  -480
/* AT..NN */
/* N */   1470  1290  1470  1070  1090
/* A */   1290  1290  1090  1070  1090
/* C */   1470  1040  1470  1040   430
/* G */   1090   710  1090   410  1090
/* T */   1160  1040  1160  1040   570
/* TA..CG */
/* N */   1320  1320   960  1320   960
/* A */   1000  1000   960   690   960
/* C */   1320  1320   590  1320   300
/* G */   1130  1130   960  -110   960
/* T */   1320  1320    10  1320    70
/* TA..GC */
/* N */   1320  1320  1290  1320   420
/* A */    690  -140   420   690   420
/* C */   1320  1320  1290  1320   370
/* G */   1130  1130   420   340   420
/* T */   1320  1320   -40  1320 -1000
/* TA..GT */
/* N */   1320  1320  1290  1320   420
/* A */    690  -140   420   690   420
/* C */   1320  1320  1290  1320   370
/* G */   1130  1130   420   340   420
/* T */   1320  1320   -40  1320 -1000
/* TA..TG */
/* N */   1320  1320   960  1320   960
/* A */   1000  1000   960   690   960
/* C */   1320  1320   590  1320   300
/* G */   1130  1130   960  -110   960
/* T */   1320  1320    10  1320    70
/* TA..AT */
/* N */   1470  1290  1470   180  1040
/* A */   1290  1290  1040   180  1040
/* C */   1470   150  1470   150   350
/* G */   1070  1070  1040   160  1040
/* T */    150   150  -160   150  -480
/* TA..TA */
/* N */   1740  1430  1740  1430  1430
/* A */   1430  1210  1430  1190  1430
/* C */   1740  1430  1740  1430   250
/* G */   1430  1190  1430   840  1430
/* T */   1430  1430   250  1430   620
/* TA..NN */
/* N */   1740  1430  1740  1430  1430
/* A */   1430  1290  1430  1190  1430
/* C */   1740  1430  1740  1430   370
/* G */   1430  1190  1430   840  1430
/* T */   1430  1430   250  1430   620
/* NN..CG */
/* N */   1320  1320   960  1320   960
/* A */   1000  1000   960   690   960
/* C */   1320  1320   750  1320   300
/* G */   1130  1130   960    10   960
/* T */   1320  1320   190  1320   210
/* NN..GC */
/* N */   1320  1320  1290  1320  1090
/* A */   1090  1000  1090   890  1090
/* C */   1320  1320  1290  1320   370
/* G */   1130  1130  1090   340  1090
/* T */   1320  1320  1160  1320   210
/* NN..GT */
/* N */   1320  1320  1290  1320  1090
/* A */   1090  1000  1090   890  1090
/* C */   1320  1320  1290  1320   370
/* G */   1130  1130  1090   340  1090
/* T */   1320  1320  1160  1320   210
/* NN..TG */
/* N */   1320  1320   960  1320   960
/* A */   1000  1000   960   690   960
/* C */   1320  1320   750  1320   300
/* G */   1130  1130   960    10   960
/* T */   1320  1320   190  1320   210
/* NN..AT */
/* N */   1470  1290  1470  1090  1160
/* A */   1290  1290  1040   710  1040
/* C */   1470  1090  1470  1090  1160
/* G */   1070  1070  1040   410  1040
/* T */   1090  1090   430  1090   570
/* NN..TA */
/* N */   1740  1430  1740  1430  1430
/* A */   1430  1290  1430  1190  1430
/* C */   1740  1430  1740  1430   250
/* G */   1430  1190  1430   840  1430
/* T */   1430  1430   370  1430   620
/* NN..NN */
/* N */   1740  1430  1740  1430  1430
/* A */   1430  1290  1430  1190  1430
/* C */   1740  1430  1740  1430  1160
/* G */   1430  1190  1430   840  1430
/* T */   1430  1430  1160  1430   620

# int21_energies
/* CG N ..CG */
   540   470   540   470   470 /* N */
   470   470   470   470   450 /* A */
   540   470   540   470   470 /* C */
   470   470   470   470   440 /* G */
   470   450   470   440   470 /* T */
/* CG A ..CG */
   470   450   470   360   470 /* N */
   450   450   450   360   450 /* A */
   470   450   470   360   470 /* C */
   360   360   360   360   360 /* G */
   470   450   470   360   470 /* T */
/* CG C ..CG */
   540   470   540   470   440 /* N */
   470   470   470   470   440 /* A */
   540   470   540   470   440 /* C */
   470   470   470   470   440 /* G */
   440   440   440   440   440 /* T */
/* CG G ..CG */
   470   360   470   400   470 /* N */
   360   360   360   360   360 /* A */
   470   360   470   400   470 /* C */
   400   360   400   400   400 /* G */
   470   360   470   400   470 /* T */
/* CG T ..CG */
   470   470   440   470   380 /* N */
   470   470   440   470   380 /* A */
   440   440   440   440   380 /* C */
   470   470   440   470   380 /* G */
   380   380   380   380   380 /* T */
/* CG N ..GC */
   510   480   510   480   460 /* N */
   480   480   480   480   460 /* A */
   510   480   510   480   460 /* C */
   480   480   480   480   460 /* G */
   460   460   460   460   460 /* T */
/* CG A ..GC */
   430   430   390   260   390 /* N */
   430   430   390   260   390 /* A */
   390   390   390   260   390 /* C */
   260   260   260   260   260 /* G */
   390   390   390   260   390 /* T */
/* CG C ..GC */
   510   480   510   480   460 /* N */
   480   480   480   480   460 /* A */
   510   480   510   480   460 /* C */
   480   480   480   480   460 /* G */
   460   460   460   460   460 /* T */
/* CG G ..GC */
   390   340   390   320   390 /* N */
   340   340   340   320   340 /* A */
   390   340   390   320   390 /* C */
   320   320   320   320   320 /* G */
   390   340   390   320   390 /* T */
/* CG T ..GC */
   480   480   450   480   450 /* N */
   480   480   450   480   450 /* A */
   450   450   450   450   450 /* C */
   480   480   450   480   450 /* G */
   450   450   450   450   450 /* T */
/* CG N ..GT */
   510   480   510   480   460 /* N */
   480   480   480   480   460 /* A */
   510   480   510   480   460 /* C */
   480   480   480   480   460 /* G */
   460   460   460   460   460 /* T */
/* CG A ..GT */
   430   430   390   260   390 /* N */
   430   430   390   260   390 /* A */
   390   390   390   260   390 /* C */
   260   260   260   260   260 /* G */
   390   390   390   260   390 /* T */
/* CG C ..GT */
   510   480   510   480   460 /* N */
   480   480   480   480   460 /* A */
   510   480   510   480   460 /* C */
   480   480   480   480   460 /* G */
   460   460   460   460   460 /* T */
/* CG G ..GT */
   390   340   390   320   390 /* N */
   340   340   340   320   340 /* A */
   390   340   390   320   390 /* C */
   320   320   320   320   320 /* G */
   390   340   390   320   390 /* T */
/* CG T ..GT */
   480   480   450   480   450 /* N */
   480   480   450   480   450 /* A */
   450   450   450   450   450 /* C */
   480   480   450   480   450 /* G */
   450   450   450   450   450 /* T */
/* CG N ..TG */
   540   470   540   470   470 /* N */
   470   470   470   470   450 /* A */
   540   470   540   470   470 /* C */
   470   470   470   470   440 /* G */
   470   450   470   440   470 /* T */
/* CG A ..TG */
   470   450   470   360   470 /* N */
   450   450   450   360   450 /* A */
   470   450   470   360   470 /* C */
   360   360   360   360   360 /* G */
   470   450   470   360   470 /* T */
/* CG C ..TG */
   540   470   540   470   440 /* N */
   470   470   470   470   440 /* A */
   540   470   540   470   440 /* C */
   470   470   470   470   440 /* G */
   440   440   440   440   440 /* T */
/* CG G ..TG */
   470   360   470   400   470 /* N */
   360   360   360   360   360 /* A */
   470   360   470   400   470 /* C */
   400   360   400   400   400 /* G */
   470   360   470   400   470 /* T */
/* CG T ..TG */
   470   470   440   470   380 /* N */
   470   470   440   470   380 /* A */
   440   440   440   440   380 /* C */
   470   470   440   470   380 /* G */
   380   380   380   380   380 /* T */
/* CG N ..AT */
   560   500   560   500   490 /* N */
   500   500   500   500   490 /* A */
   560   500   560   500   490 /* C */
   500   500   500   500   340 /* G */
   490   490   490   340   490 /* T */
/* CG A ..AT */
   500   500   490   300   490 /* N */
   500   500   490   300   490 /* A */
   490   490   490   300   490 /* C */
   300   300   300   300   300 /* G */
   490   490   490   300   490 /* T */
/* CG C ..AT */
   560   500   560   500   320 /* N */
   500   500   500   500   320 /* A */
   560   500   560   500   320 /* C */
   500   500   500   500   320 /* G */
   320   320   320   320   320 /* T */
/* CG G ..AT */
   490   370   490   310   490 /* N */
   370   370   370   310   370 /* A */
   490   370   490   310   490 /* C */
   310   310   310   310   310 /* G */
   490   370   490   310   490 /* T */
/* CG T ..AT */
   500   500   380   500   340 /* N */
   500   500   380   500   340 /* A */
   380   380   380   380   340 /* C */
   500   500   380   500   340 /* G */
   340   340   340   340   340 /* T */
/* CG N ..TA */
   600   520   600   520   550 /* N */
   520   520   520   520   520 /* A */
   600   520   600   520   550 /* C */
   520   520   520   520   500 /* G */
   550   520   550   500   550 /* T */
/* CG A ..TA */
   550   520   550   450   550 /* N */
   520   520   520   450   520 /* A */
   550   520   550   450   550 /* C */
   450   450   450   450   450 /* G */
   550   520   550   450   550 /* T */
/* CG C ..TA */
   600   520   600   520   500 /* N */
   520   520   520   520   500 /* A */
   600   520   600   520   500 /* C */
   520   520   520   520   500 /* G */
   500   500   500   500   500 /* T */
/* CG G ..TA */
   550   420   550   410   550 /* N */
   420   420   420   410   420 /* A */
   550   420   550   410   550 /* C */
   410   410   410   410   410 /* G */
   550   420   550   410   550 /* T */
/* CG T ..TA */
   520   520   500   520   460 /* N */
   520   520   500   520   460 /* A */
   500   500   500   500   460 /* C */
   520   520   500   520   460 /* G */
   460   460   460   460   460 /* T */
/* CG N ..NN */
   600   520   600   520   550 /* N */
   520   520   520   520   520 /* A */
   600   520   600   520   550 /* C */
   520   520   520   520   500 /* G */
   550   520   550   500   550 /* T */
/* CG A ..NN */
   550   520   550   450   550 /* N */
   520   520   520   450   520 /* A */
   550   520   550   450   550 /* C */
   450   450   450   450   450 /* G */
   550   520   550   450   550 /* T */
/* CG C ..NN */
   600   520   600   520   500 /* N */
   520   520   520   520   500 /* A */
   600   520   600   520   500 /* C */
   520   520   520   520   500 /* G */
   500   500   500   500   500 /* T */
/* CG G ..NN */
   550   420   550   410   550 /* N */
   420   420   420   410   420 /* A */
   550   420   550   410   550 /* C */
   410   410   410   410   410 /* G */
   550   420   550   410   550 /* T */
/* CG T ..NN */
   520   520   500   520   460 /* N */
   520   520   500   520   460 /* A */
   500   500   500   500   460 /* C */
   520   520   500   520   460 /* G */
   460   460   460   460   460 /* T */
/* GC N ..CG */
   510   430   510   390   480 /* N */
   430   430   430   390   430 /* A */
   510   430   510   390   480 /* C */
   390   390   390   390   390 /* G */
   480   430   480   390   480 /* T */
/* GC A ..CG */
   480   430   480   340   480 /* N */
   430   430   430   340   430 /* A */
   480   430   480   340   480 /* C */
   340   340   340   340   340 /* G */
   480   430   480   340   480 /* T */
/* GC C ..CG */
   510   390   510   390   450 /* N */
   390   390   390   390   390 /* A */
   510   390   510   390   450 /* C */
   390   390   390   390   390 /* G */
   450   390   450   390   450 /* T */
/* GC G ..CG */
   480   260   480   320   480 /* N */
   260   260   260   260   260 /* A */
   480   260   480   320   480 /* C */
   320   260   320   320   320 /* G */
   480   260   480   320   480 /* T */
/* GC T ..CG */
   460   390   460   390   450 /* N */
   390   390   390   390   390 /* A */
   460   390   460   390   450 /* C */
   390   390   390   390   390 /* G */
   450   390   450   390   450 /* T */
/* GC N ..GC */
   480   420   480   420   470 /* N */
   420   420   420   420   420 /* A */
   480   420   480   420   450 /* C */
   420   420   420   420   420 /* G */
   470   420   450   420   470 /* T */
/* GC A ..GC */
   420   400   420   230   420 /* N */
   400   400   400   230   400 /* A */
   420   400   420   230   420 /* C */
   230   230   230   230   230 /* G */
   420   400   420   230   420 /* T */
/* GC C ..GC */
   480   420   480   420   450 /* N */
   420   420   420   420   420 /* A */
   480   420   480   420   450 /* C */
   420   420   420   420   420 /* G */
   450   420   450   420   450 /* T */
/* GC G ..GC */
   420   230   420   200   420 /* N */
   230   230   230   200   230 /* A */
   420   230   420   200   420 /* C */
   200   200   200   200   200 /* G */
   420   230   420   200   420 /* T */
/* GC T ..GC */
   470   420   450   420   470 /* N */
   420   420   420   420   420 /* A */
   450   420   450   420   450 /* C */
   420   420   420   420   420 /* G */
   470   420   450   420   470 /* T */
/* GC N ..GT */
   480   420   480   420   470 /* N */
   420   420   420   420   420 /* A */
   480   420   480   420   450 /* C */
   420   420   420   420   420 /* G */
   470   420   450   420   470 /* T */
/* GC A ..GT */
   420   400   420   230   420 /* N */
   400   400   400   230   400 /* A */
   420   400   420   230   420 /* C */
   230   230   230   230   230 /* G */
   420   400   420   230   420 /* T */
/* GC C ..GT */
   480   420   480   420   450 /* N */
   420   420   420   420   420 /* A */
   480   420   480   420   450 /* C */
   420   420   420   420   420 /* G */
   450   420   450   420   450 /* T */
/* GC G ..GT */
   420   230   420   200   420 /* N */
   230   230   230   200   230 /* A */
   420   230   420   200   420 /* C */
   200   200   200   200   200 /* G */
   420   230   420   200   420 /* T */
/* GC T ..GT */
   470   420   450   420   470 /* N */
   420   420   420   420   420 /* A */
   450   420   450   420   450 /* C */
   420   420   420   420   420 /* G */
   470   420   450   420   470 /* T */
/* GC N ..TG */
   510   430   510   390   480 /* N */
   430   430   430   390   430 /* A */
   510   430   510   390   480 /* C */
   390   390   390   390   390 /* G */
   480   430   480   390   480 /* T */
/* GC A ..TG */
   480   430   480   340   480 /* N */
   430   430   430   340   430 /* A */
   480   430   480   340   480 /* C */
   340   340   340   340   340 /* G */
   480   430   480   340   480 /* T */
/* GC C ..TG */
   510   390   510   390   450 /* N */
   390   390   390   390   390 /* A */
   510   390   510   390   450 /* C */
   390   390   390   390   390 /* G */
   450   390   450   390   450 /* T */
/* GC G ..TG */
   480   260   480   320   480 /* N */
   260   260   260   260   260 /* A */
   480   260   480   320   480 /* C */
   320   260   320   320   320 /* G */
   480   260   480   320   480 /* T */
/* GC T ..TG */
   460   390   460   390   450 /* N */
   390   390   390   390   390 /* A */
   460   390   460   390   450 /* C */
   390   390   390   390   390 /* G */
   450   390   450   390   450 /* T */
/* GC N ..AT */
   640   520   640   510   540 /* N */
   520   520   510   510   510 /* A */
   640   510   640   510   540 /* C */
   510   510   510   510   510 /* G */
   540   510   540   510   540 /* T */
/* GC A ..AT */
   520   520   500   360   500 /* N */
   520   520   500   360   500 /* A */
   500   500   500   360   500 /* C */
   360   360   360   360   360 /* G */
   500   500   500   360   500 /* T */
/* GC C ..AT */
   640   510   640   510   540 /* N */
   510   510   510   510   510 /* A */
   640   510   640   510   540 /* C */
   510   510   510   510   510 /* G */
   540   510   540   510   540 /* T */
/* GC G ..AT */
   500   390   500   270   500 /* N */
   390   390   390   270   390 /* A */
   500   390   500   270   500 /* C */
   270   270   270   270   270 /* G */
   500   390   500   270   500 /* T */
/* GC T ..AT */
   510   510   470   510   470 /* N */
   510   510   470   510   470 /* A */
   470   470   470   470   470 /* C */
   510   510   470   510   470 /* G */
   470   470   470   470   470 /* T */
/* GC N ..TA */
   580   500   580   500   550 /* N */
   500   500   500   500   500 /* A */
   580   500   580   500   550 /* C */
   500   500   500   500   500 /* G */
   550   500   550   500   550 /* T */
/* GC A ..TA */
   550   450   550   450   550 /* N */
   450   450   450   450   450 /* A */
   550   450   550   450   550 /* C */
   450   450   450   450   450 /* G */
   550   450   550   450   550 /* T */
/* GC C ..TA */
   580   500   580   500   510 /* N */
   500   500   500   500   500 /* A */
   580   500   580   500   510 /* C */
   500   500   500   500   500 /* G */
   510   500   510   500   510 /* T */
/* GC G ..TA */
   550   420   550   430   550 /* N */
   420   420   420   420   420 /* A */
   550   420   550   430   550 /* C */
   430   420   430   430   430 /* G */
   550   420   550   430   550 /* T */
/* GC T ..TA */
   580   500   580   500   420 /* N */
   500   500   500   500   420 /* A */
   580   500   580   500   420 /* C */
   500   500   500   500   420 /* G */
   420   420   420   420   420 /* T */
/* GC N ..NN */
   640   520   640   510   550 /* N */
   520   520   510   510   510 /* A */
   640   510   640   510   550 /* C */
   510   510   510   510   510 /* G */
   550   510   550   510   550 /* T */
/* GC A ..NN */
   550   520   550   450   550 /* N */
   520   520   500   450   500 /* A */
   550   500   550   450   550 /* C */
   450   450   450   450   450 /* G */
   550   500   550   450   550 /* T */
/* GC C ..NN */
   640   510   640   510   540 /* N */
   510   510   510   510   510 /* A */
   640   510   640   510   540 /* C */
   510   510   510   510   510 /* G */
   540   510   540   510   540 /* T */
/* GC G ..NN */
   550   420   550   430   550 /* N */
   420   420   420   420   420 /* A */
   550   420   550   430   550 /* C */
   430   420   430   430   430 /* G */
   550   420   550   430   550 /* T */
/* GC T ..NN */
   580   510   580   510   470 /* N */
   510   510   500   510   470 /* A */
   580   500   580   500   470 /* C */
   510   510   500   510   470 /* G */
   470   470   470   470   470 /* T */
/* GT N ..CG */
   510   430   510   390   480 /* N */
   430   430   430   390   430 /* A */
   510   430   510   390   480 /* C */
   390   390   390   390   390 /* G */
   480   430   480   390   480 /* T */
/* GT A ..CG */
   480   430   480   340   480 /* N */
   430   430   430   340   430 /* A */
   480   430   480   340   480 /* C */
   340   340   340   340   340 /* G */
   480   430   480   340   480 /* T */
/* GT C ..CG */
   510   390   510   390   450 /* N */
   390   390   390   390   390 /* A */
   510   390   510   390   450 /* C */
   390   390   390   390   390 /* G */
   450   390   450   390   450 /* T */
/* GT G ..CG */
   480   260   480   320   480 /* N */
   260   260   260   260   260 /* A */
   480   260   480   320   480 /* C */
   320   260   320   320   320 /* G */
   480   260   480   320   480 /* T */
/* GT T ..CG */
   460   390   460   390   450 /* N */
   390   390   390   390   390 /* A */
   460   390   460   390   450 /* C */
   390   390   390   390   390 /* G */
   450   390   450   390   450 /* T */
/* GT N ..GC */
   480   420   480   420   470 /* N */
   420   420   420   420   420 /* A */
   480   420   480   420   450 /* C */
   420   420   420   420   420 /* G */
   470   420   450   420   470 /* T */
/* GT A ..GC */
   420   400   420   230   420 /* N */
   400   400   400   230   400 /* A */
   420   400   420   230   420 /* C */
   230   230   230   230   230 /* G */
   420   400   420   230   420 /* T */
/* GT C ..GC */
   480   420   480   420   450 /* N */
   420   420   420   420   420 /* A */
   480   420   480   420   450 /* C */
   420   420   420   420   420 /* G */
   450   420   450   420   450 /* T */
/* GT G ..GC */
   420   230   420   200   420 /* N */
   230   230   230   200   230 /* A */
   420   230   420   200   420 /* C */
   200   200   200   200   200 /* G */
   420   230   420   200   420 /* T */
/* GT T ..GC */
   470   420   450   420   470 /* N */
   420   420   420   420   420 /* A */
   450   420   450   420   450 /* C */
   420   420   420   420   420 /* G */
   470   420   450   420   470 /* T */
/* GT N ..GT */
   480   420   480   420   470 /* N */
   420   420   420   420   420 /* A */
   480   420   480   420   450 /* C */
   420   420   420   420   420 /* G */
   470   420   450   420   470 /* T */
/* GT A ..GT */
   420   400   420   230   420 /* N */
   400   400   400   230   400 /* A */
   420   400   420   230   420 /* C */
   230   230   230   230   230 /* G */
   420   400   420   230   420 /* T */
/* GT C ..GT */
   480   420   480   420   450 /* N */
   420   420   420   420   420 /* A */
   480   420   480   420   450 /* C */
   420   420   420   420   420 /* G */
   450   420   450   420   450 /* T */
/* GT G ..GT */
   420   230   420   200   420 /* N */
   230   230   230   200   230 /* A */
   420   230   420   200   420 /* C */
   200   200   200   200   200 /* G */
   420   230   420   200   420 /* T */
/* GT T ..GT */
   470   420   450   420   470 /* N */
   420   420   420   420   420 /* A */
   450   420   450   420   450 /* C */
   420   420   420   420   420 /* G */
   470   420   450   420   470 /* T */
/* GT N ..TG */
   510   430   510   390   480 /* N */
   430   430   430   390   430 /* A */
   510   430   510   390   480 /* C */
   390   390   390   390   390 /* G */
   480   430   480   390   480 /* T */
/* GT A ..TG */
   480   430   480   340   480 /* N */
   430   430   430   340   430 /* A */
   480   430   480   340   480 /* C */
   340   340   340   340   340 /* G */
   480   430   480   340   480 /* T */
/* GT C ..TG */
   510   390   510   390   450 /* N */
   390   390   390   390   390 /* A */
   510   390   510   390   450 /* C */
   390   390   390   390   390 /* G */
   450   390   450   390   450 /* T */
/* GT G ..TG */
   480   260   480   320   480 /* N */
   260   260   260   260   260 /* A */
   480   260   480   320   480 /* C */
   320   260   320   320   320 /* G */
   480   260   480   320   480 /* T */
/* GT T ..TG */
   460   390   460   390   450 /* N */
   390   390   390   390   390 /* A */
   460   390   460   390   450 /* C */
   390   390   390   390   390 /* G */
   450   390   450   390   450 /* T */
/* GT N ..AT */
   640   520   640   510   540 /* N */
   520   520   510   510   510 /* A */
   640   510   640   510   540 /* C */
   510   510   510   510   510 /* G */
   540   510   540   510   540 /* T */
/* GT A ..AT */
   520   520   500   360   500 /* N */
   520   520   500   360   500 /* A */
   500   500   500   360   500 /* C */
   360   360   360   360   360 /* G */
   500   500   500   360   500 /* T */
/* GT C ..AT */
   640   510   640   510   540 /* N */
   510   510   510   510   510 /* A */
   640   510   640   510   540 /* C */
   510   510   510   510   510 /* G */
   540   510   540   510   540 /* T */
/* GT G ..AT */
   500   390   500   270   500 /* N */
   390   390   390   270   390 /* A */
   500   390   500   270   500 /* C */
   270   270   270   270   270 /* G */
   500   390   500   270   500 /* T */
/* GT T ..AT */
   510   510   470   510   470 /* N */
   510   510   470   510   470 /* A */
   470   470   470   470   470 /* C */
   510   510   470   510   470 /* G */
   470   470   470   470   470 /* T */
/* GT N ..TA */
   580   500   580   500   550 /* N */
   500   500   500   500   500 /* A */
   580   500   580   500   550 /* C */
   500   500   500   500   500 /* G */
   550   500   550   500   550 /* T */
/* GT A ..TA */
   550   450   550   450   550 /* N */
   450   450   450   450   450 /* A */
   550   450   550   450   550 /* C */
   450   450   450   450   450 /* G */
   550   450   550   450   550 /* T */
/* GT C ..TA */
   580   500   580   500   510 /* N */
   500   500   500   500   500 /* A */
   580   500   580   500   510 /* C */
   500   500   500   500   500 /* G */
   510   500   510   500   510 /* T */
/* GT G ..TA */
   550   420   550   430   550 /* N */
   420   420   420   420   420 /* A */
   550   420   550   430   550 /* C */
   430   420   430   430   430 /* G */
   550   420   550   430   550 /* T */
/* GT T ..TA */
   580   500   580   500   420 /* N */
   500   500   500   500   420 /* A */
   580   500   580   500   420 /* C */
   500   500   500   500   420 /* G */
   420   420   420   420   420 /* T */
/* GT N ..NN */
   640   520   640   510   550 /* N */
   520   520   510   510   510 /* A */
   640   510   640   510   550 /* C */
   510   510   510   510   510 /* G */
   550   510   550   510   550 /* T */
/* GT A ..NN */
   550   520   550   450   550 /* N */
   520   520   500   450   500 /* A */
   550   500   550   450   550 /* C */
   450   450   450   450   450 /* G */
   550   500   550   450   550 /* T */
/* GT C ..NN */
   640   510   640   510   540 /* N */
   510   510   510   510   510 /* A */
   640   510   640   510   540 /* C */
   510   510   510   510   510 /* G */
   540   510   540   510   540 /* T */
/* GT G ..NN */
   550   420   550   430   550 /* N */
   420   420   420   420   420 /* A */
   550   420   550   430   550 /* C */
   430   420   430   430   430 /* G */
   550   420   550   430   550 /* T */
/* GT T ..NN */
   580   510   580   510   470 /* N */
   510   510   500   510   470 /* A */
   580   500   580   500   470 /* C */
   510   510   500   510   470 /* G */
   470   470   470   470   470 /* T */
/* TG N ..CG */
   540   470   540   470   470 /* N */
   470   470   470   470   450 /* A */
   540   470   540   470   470 /* C */
   470   470   470   470   440 /* G */
   470   450   470   440   470 /* T */
/* TG A ..CG */
   470   450   470   360   470 /* N */
   450   450   450   360   450 /* A */
   470   450   470   360   470 /* C */
   360   360   360   360   360 /* G */
   470   450   470   360   470 /* T */
/* TG C ..CG */
   540   470   540   470   440 /* N */
   470   470   470   470   440 /* A */
   540   470   540   470   440 /* C */
   470   470   470   470   440 /* G */
   440   440   440   440   440 /* T */
/* TG G ..CG */
   470   360   470   400   470 /* N */
   360   360   360   360   360 /* A */
   470   360   470   400   470 /* C */
   400   360   400   400   400 /* G */
   470   360   470   400   470 /* T */
/* TG T ..CG */
   470   470   440   470   380 /* N */
   470   470   440   470   380 /* A */
   440   440   440   440   380 /* C */
   470   470   440   470   380 /* G */
   380   380   380   380   380 /* T */
/* TG N ..GC */
   510   480   510   480   460 /* N */
   480   480   480   480   460 /* A */
   510   480   510   480   460 /* C */
   480   480   480   480   460 /* G */
   460   460   460   460   460 /* T */
/* TG A ..GC */
   430   430   390   260   390 /* N */
   430   430   390   260   390 /* A */
   390   390   390   260   390 /* C */
   260   260   260   260   260 /* G */
   390   390   390   260   390 /* T */
/* TG C ..GC */
   510   480   510   480   460 /* N */
   480   480   480   480   460 /* A */
   510   480   510   480   460 /* C */
   480   480   480   480   460 /* G */
   460   460   460   460   460 /* T */
/* TG G ..GC */
   390   340   390   320   390 /* N */
   340   340   340   320   340 /* A */
   390   340   390   320   390 /* C */
   320   320   320   320   320 /* G */
   390   340   390   320   390 /* T */
/* TG T ..GC */
   480   480   450   480   450 /* N */
   480   480   450   480   450 /* A */
   450   450   450   450   450 /* C */
   480   480   450   480   450 /* G */
   450   450   450   450   450 /* T */
/* TG N ..GT */
   510   480   510   480   460 /* N */
   480   480   480   480   460 /* A */
   510   480   510   480   460 /* C */
   480   480   480   480   460 /* G */
   460   460   460   460   460 /* T */
/* TG A ..GT */
   430   430   390   260   390 /* N */
   430   430   390   260   390 /* A */
   390   390   390   260   390 /* C */
   260   260   260   260   260 /* G */
   390   390   390   260   390 /* T */
/* TG C ..GT */
   510   480   510   480   460 /* N */
   480   480   480   480   460 /* A */
   510   480   510   480   460 /* C */
   480   480   480   480   460 /* G */
   460   460   460   460   460 /* T */
/* TG G ..GT */
   390   340   390   320   390 /* N */
   340   340   340   320   340 /* A */
   390   340   390   320   390 /* C */
   320   320   320   320   320 /* G */
   390   340   390   320   390 /* T */
/* TG T ..GT */
   480   480   450   480   450 /* N */
   480   480   450   480   450 /* A */
   450   450   450   450   450 /* C */
   480   480   450   480   450 /* G */
   450   450   450   450   450 /* T */
/* TG N ..TG */
   540   470   540   470   470 /* N */
   470   470   470   470   450 /* A */
   540   470   540   470   470 /* C */
   470   470   470   470   440 /* G */
   470   450   470   440   470 /* T */
/* TG A ..TG */
   470   450   470   360   470 /* N */
   450   450   450   360   450 /* A */
   470   450   470   360   470 /* C */
   360   360   360   360   360 /* G */
   470   450   470   360   470 /* T */
/* TG C ..TG */
   540   470   540   470   440 /* N */
   470   470   470   470   440 /* A */
   540   470   540   470   440 /* C */
   470   470   470   470   440 /* G */
   440   440   440   440   440 /* T */
/* TG G ..TG */
   470   360   470   400   470 /* N */
   360   360   360   360   360 /* A */
   470   360   470   400   470 /* C */
   400   360   400   400   400 /* G */
   470   360   470   400   470 /* T */
/* TG T ..TG */
   470   470   440   470   380 /* N */
   470   470   440   470   380 /* A */
   440   440   440   440   380 /* C */
   470   470   440   470   380 /* G */
   380   380   380   380   380 /* T */
/* TG N ..AT */
   560   500   560   500   490 /* N */
   500   500   500   500   490 /* A */
   560   500   560   500   490 /* C */
   500   500   500   500   340 /* G */
   490   490   490   340   490 /* T */
/* TG A ..AT */
   500   500   490   300   490 /* N */
   500   500   490   300   490 /* A */
   490   490   490   300   490 /* C */
   300   300   300   300   300 /* G */
   490   490   490   300   490 /* T */
/* TG C ..AT */
   560   500   560   500   320 /* N */
   500   500   500   500   320 /* A */
   560   500   560   500   320 /* C */
   500   500   500   500   320 /* G */
   320   320   320   320   320 /* T */
/* TG G ..AT */
   490   370   490   310   490 /* N */
   370   370   370   310   370 /* A */
   490   370   490   310   490 /* C */
   310   310   310   310   310 /* G */
   490   370   490   310   490 /* T */
/* TG T ..AT */
   500   500   380   500   340 /* N */
   500   500   380   500   340 /* A */
   380   380   380   380   340 /* C */
   500   500   380   500   340 /* G */
   340   340   340   340   340 /* T */
/* TG N ..TA */
   600   520   600   520   550 /* N */
   520   520   520   520   520 /* A */
   600   520   600   520   550 /* C */
   520   520   520   520   500 /* G */
   550   520   550   500   550 /* T */
/* TG A ..TA */
   550   520   550   450   550 /* N */
   520   520   520   450   520 /* A */
   550   520   550   450   550 /* C */
   450   450   450   450   450 /* G */
   550   520   550   450   550 /* T */
/* TG C ..TA */
   600   520   600   520   500 /* N */
   520   520   520   520   500 /* A */
   600   520   600   520   500 /* C */
   520   520   520   520   500 /* G */
   500   500   500   500   500 /* T */
/* TG G ..TA */
   550   420   550   410   550 /* N */
   420   420   420   410   420 /* A */
   550   420   550   410   550 /* C */
   410   410   410   410   410 /* G */
   550   420   550   410   550 /* T */
/* TG T ..TA */
   520   520   500   520   460 /* N */
   520   520   500   520   460 /* A */
   500   500   500   500   460 /* C */
   520   520   500   520   460 /* G */
   460   460   460   460   460 /* T */
/* TG N ..NN */
   600   520   600   520   550 /* N */
   520   520   520   520   520 /* A */
   600   520   600   520   550 /* C */
   520   520   520   520   500 /* G */
   550   520   550   500   550 /* T */
/* TG A ..NN */
   550   520   550   450   550 /* N */
   520   520   520   450   520 /* A */
   550   520   550   450   550 /* C */
   450   450   450   450   450 /* G */
   550   520   550   450   550 /* T */
/* TG C ..NN */
   600   520   600   520   500 /* N */
   520   520   520   520   500 /* A */
   600   520   600   520   500 /* C */
   520   520   520   520   500 /* G */
   500   500   500   500   500 /* T */
/* TG G ..NN */
   550   420   550   410   550 /* N */
   420   420   420   410   420 /* A */
   550   420   550   410   550 /* C */
   410   410   410   410   410 /* G */
   550   420   550   410   550 /* T */
/* TG T ..NN */
   520   520   500   520   460 /* N */
   520   520   500   520   460 /* A */
   500   500   500   500   460 /* C */
   520   520   500   520   460 /* G */
   460   460   460   460   460 /* T */
/* AT N ..CG */
   560   500   560   490   500 /* N */
   500   500   500   490   500 /* A */
   560   500   560   490   500 /* C */
   490   490   490   490   380 /* G */
   500   500   500   380   500 /* T */
/* AT A ..CG */
   500   500   500   370   500 /* N */
   500   500   500   370   500 /* A */
   500   500   500   370   500 /* C */
   370   370   370   370   370 /* G */
   500   500   500   370   500 /* T */
/* AT C ..CG */
   560   490   560   490   380 /* N */
   490   490   490   490   380 /* A */
   560   490   560   490   380 /* C */
   490   490   490   490   380 /* G */
   380   380   380   380   380 /* T */
/* AT G ..CG */
   500   300   500   310   500 /* N */
   300   300   300   300   300 /* A */
   500   300   500   310   500 /* C */
   310   300   310   310   310 /* G */
   500   300   500   310   500 /* T */
/* AT T ..CG */
   490   490   320   490   340 /* N */
   490   490   320   490   340 /* A */
   320   320   320   320   320 /* C */
   490   490   320   490   340 /* G */
   340   340   320   340   340 /* T */
/* AT N ..GC */
   640   520   640   500   510 /* N */
   520   520   510   500   510 /* A */
   640   510   640   500   510 /* C */
   500   500   500   500   470 /* G */
   510   510   510   470   510 /* T */
/* AT A ..GC */
   520   520   510   390   510 /* N */
   520   520   510   390   510 /* A */
   510   510   510   390   510 /* C */
   390   390   390   390   390 /* G */
   510   510   510   390   510 /* T */
/* AT C ..GC */
   640   500   640   500   470 /* N */
   500   500   500   500   470 /* A */
   640   500   640   500   470 /* C */
   500   500   500   500   470 /* G */
   470   470   470   470   470 /* T */
/* AT G ..GC */
   510   360   510   270   510 /* N */
   360   360   360   270   360 /* A */
   510   360   510   270   510 /* C */
   270   270   270   270   270 /* G */
   510   360   510   270   510 /* T */
/* AT T ..GC */
   540   500   540   500   470 /* N */
   500   500   500   500   470 /* A */
   540   500   540   500   470 /* C */
   500   500   500   500   470 /* G */
   470   470   470   470   470 /* T */
/* AT N ..GT */
   640   520   640   500   510 /* N */
   520   520   510   500   510 /* A */
   640   510   640   500   510 /* C */
   500   500   500   500   470 /* G */
   510   510   510   470   510 /* T */
/* AT A ..GT */
   520   520   510   390   510 /* N */
   520   520   510   390   510 /* A */
   510   510   510   390   510 /* C */
   390   390   390   390   390 /* G */
   510   510   510   390   510 /* T */
/* AT C ..GT */
   640   500   640   500   470 /* N */
   500   500   500   500   470 /* A */
   640   500   640   500   470 /* C */
   500   500   500   500   470 /* G */
   470   470   470   470   470 /* T */
/* AT G ..GT */
   510   360   510   270   510 /* N */
   360   360   360   270   360 /* A */
   510   360   510   270   510 /* C */
   270   270   270   270   270 /* G */
   510   360   510   270   510 /* T */
/* AT T ..GT */
   540   500   540   500   470 /* N */
   500   500   500   500   470 /* A */
   540   500   540   500   470 /* C */
   500   500   500   500   470 /* G */
   470   470   470   470   470 /* T */
/* AT N ..TG */
   560   500   560   490   500 /* N */
   500   500   500   490   500 /* A */
   560   500   560   490   500 /* C */
   490   490   490   490   380 /* G */
   500   500   500   380   500 /* T */
/* AT A ..TG */
   500   500   500   370   500 /* N */
   500   500   500   370   500 /* A */
   500   500   500   370   500 /* C */
   370   370   370   370   370 /* G */
   500   500   500   370   500 /* T */
/* AT C ..TG */
   560   490   560   490   380 /* N */
   490   490   490   490   380 /* A */
   560   490   560   490   380 /* C */
   490   490   490   490   380 /* G */
   380   380   380   380   380 /* T */
/* AT G ..TG */
   500   300   500   310   500 /* N */
   300   300   300   300   300 /* A */
   500   300   500   310   500 /* C */
   310   300   310   310   310 /* G */
   500   300   500   310   500 /* T */
/* AT T ..TG */
   490   490   320   490   340 /* N */
   490   490   320   490   340 /* A */
   320   320   320   320   320 /* C */
   490   490   320   490   340 /* G */
   340   340   320   340   340 /* T */
/* AT N ..AT */
   670   500   670   480   520 /* N */
   500   500   480   480   500 /* A */
   670   480   670   480   480 /* C */
   480   480   480   480   480 /* G */
   520   500   480   480   520 /* T */
/* AT A ..AT */
   500   500   480   340   500 /* N */
   500   500   480   340   500 /* A */
   480   480   480   340   480 /* C */
   340   340   340   340   340 /* G */
   500   500   480   340   500 /* T */
/* AT C ..AT */
   670   480   670   480   440 /* N */
   480   480   480   480   440 /* A */
   670   480   670   480   440 /* C */
   480   480   480   480   440 /* G */
   440   440   440   440   440 /* T */
/* AT G ..AT */
   480   340   480   370   480 /* N */
   340   340   340   340   340 /* A */
   480   340   480   370   480 /* C */
   370   340   370   370   370 /* G */
   480   340   480   370   480 /* T */
/* AT T ..AT */
   520   480   440   480   520 /* N */
   480   480   440   480   480 /* A */
   440   440   440   440   440 /* C */
   480   480   440   480   480 /* G */
   520   480   440   480   520 /* T */
/* AT N ..TA */
   610   510   610   500   530 /* N */
   510   510   510   500   510 /* A */
   610   510   610   500   530 /* C */
   500   500   500   500   500 /* G */
   530   510   530   500   530 /* T */
/* AT A ..TA */
   530   510   530   430   530 /* N */
   510   510   510   430   510 /* A */
   530   510   530   430   530 /* C */
   430   430   430   430   430 /* G */
   530   510   530   430   530 /* T */
/* AT C ..TA */
   610   500   610   500   500 /* N */
   500   500   500   500   500 /* A */
   610   500   610   500   500 /* C */
   500   500   500   500   500 /* G */
   500   500   500   500   500 /* T */
/* AT G ..TA */
   530   380   530   420   530 /* N */
   380   380   380   380   380 /* A */
   530   380   530   420   530 /* C */
   420   380   420   420   420 /* G */
   530   380   530   420   530 /* T */
/* AT T ..TA */
   510   500   510   500   410 /* N */
   500   500   500   500   410 /* A */
   510   500   510   500   410 /* C */
   500   500   500   500   410 /* G */
   410   410   410   410   410 /* T */
/* AT N ..NN */
   670   520   670   500   530 /* N */
   520   520   510   500   510 /* A */
   670   510   670   500   530 /* C */
   500   500   500   500   500 /* G */
   530   510   530   500   530 /* T */
/* AT A ..NN */
   530   520   530   430   530 /* N */
   520   520   510   430   510 /* A */
   530   510   530   430   530 /* C */
   430   430   430   430   430 /* G */
   530   510   530   430   530 /* T */
/* AT C ..NN */
   670   500   670   500   500 /* N */
   500   500   500   500   500 /* A */
   670   500   670   500   500 /* C */
   500   500   500   500   500 /* G */
   500   500   500   500   500 /* T */
/* AT G ..NN */
   530   380   530   420   530 /* N */
   380   380   380   380   380 /* A */
   530   380   530   420   530 /* C */
   420   380   420   420   420 /* G */
   530   380   530   420   530 /* T */
/* AT T ..NN */
   540   500   540   500   520 /* N */
   500   500   500   500   480 /* A */
   540   500   540   500   470 /* C */
   500   500   500   500   480 /* G */
   520   480   470   480   520 /* T */
/* TA N ..CG */
   600   550   600   550   520 /* N */
   550   550   550   550   520 /* A */
   600   550   600   550   520 /* C */
   550   550   550   550   500 /* G */
   520   520   520   500   520 /* T */
/* TA A ..CG */
   520   520   520   420   520 /* N */
   520   520   520   420   520 /* A */
   520   520   520   420   520 /* C */
   420   420   420   420   420 /* G */
   520   520   520   420   520 /* T */
/* TA C ..CG */
   600   550   600   550   500 /* N */
   550   550   550   550   500 /* A */
   600   550   600   550   500 /* C */
   550   550   550   550   500 /* G */
   500   500   500   500   500 /* T */
/* TA G ..CG */
   520   450   520   410   520 /* N */
   450   450   450   410   450 /* A */
   520   450   520   410   520 /* C */
   410   410   410   410   410 /* G */
   520   450   520   410   520 /* T */
/* TA T ..CG */
   550   550   500   550   460 /* N */
   550   550   500   550   460 /* A */
   500   500   500   500   460 /* C */
   550   550   500   550   460 /* G */
   460   460   460   460   460 /* T */
/* TA N ..GC */
   580   550   580   550   580 /* N */
   550   550   550   550   550 /* A */
   580   550   580   550   580 /* C */
   550   550   550   550   550 /* G */
   580   550   580   550   580 /* T */
/* TA A ..GC */
   500   450   500   420   500 /* N */
   450   450   450   420   450 /* A */
   500   450   500   420   500 /* C */
   420   420   420   420   420 /* G */
   500   450   500   420   500 /* T */
/* TA C ..GC */
   580   550   580   550   580 /* N */
   550   550   550   550   550 /* A */
   580   550   580   550   580 /* C */
   550   550   550   550   550 /* G */
   580   550   580   550   580 /* T */
/* TA G ..GC */
   500   450   500   430   500 /* N */
   450   450   450   430   450 /* A */
   500   450   500   430   500 /* C */
   430   430   430   430   430 /* G */
   500   450   500   430   500 /* T */
/* TA T ..GC */
   550   550   510   550   420 /* N */
   550   550   510   550   420 /* A */
   510   510   510   510   420 /* C */
   550   550   510   550   420 /* G */
   420   420   420   420   420 /* T */
/* TA N ..GT */
   580   550   580   550   580 /* N */
   550   550   550   550   550 /* A */
   580   550   580   550   580 /* C */
   550   550   550   550   550 /* G */
   580   550   580   550   580 /* T */
/* TA A ..GT */
   500   450   500   420   500 /* N */
   450   450   450   420   450 /* A */
   500   450   500   420   500 /* C */
   420   420   420   420   420 /* G */
   500   450   500   420   500 /* T */
/* TA C ..GT */
   580   550   580   550   580 /* N */
   550   550   550   550   550 /* A */
   580   550   580   550   580 /* C */
   550   550   550   550   550 /* G */
   580   550   580   550   580 /* T */
/* TA G ..GT */
   500   450   500   430   500 /* N */
   450   450   450   430   450 /* A */
   500   450   500   430   500 /* C */
   430   430   430   430   430 /* G */
   500   450   500   430   500 /* T */
/* TA T ..GT */
   550   550   510   550   420 /* N */
   550   550   510   550   420 /* A */
   510   510   510   510   420 /* C */
   550   550   510   550   420 /* G */
   420   420   420   420   420 /* T */
/* TA N ..TG */
   600   550   600   550   520 /* N */
   550   550   550   550   520 /* A */
   600   550   600   550   520 /* C */
   550   550   550   550   500 /* G */
   520   520   520   500   520 /* T */
/* TA A ..TG */
   520   520   520   420   520 /* N */
   520   520   520   420   520 /* A */
   520   520   520   420   520 /* C */
   420   420   420   420   420 /* G */
   520   520   520   420   520 /* T */
/* TA C ..TG */
   600   550   600   550   500 /* N */
   550   550   550   550   500 /* A */
   600   550   600   550   500 /* C */
   550   550   550   550   500 /* G */
   500   500   500   500   500 /* T */
/* TA G ..TG */
   520   450   520   410   520 /* N */
   450   450   450   410   450 /* A */
   520   450   520   410   520 /* C */
   410   410   410   410   410 /* G */
   520   450   520   410   520 /* T */
/* TA T ..TG */
   550   550   500   550   460 /* N */
   550   550   500   550   460 /* A */
   500   500   500   500   460 /* C */
   550   550   500   550   460 /* G */
   460   460   460   460   460 /* T */
/* TA N ..AT */
   610   530   610   530   510 /* N */
   530   530   530   530   510 /* A */
   610   530   610   530   510 /* C */
   530   530   530   530   510 /* G */
   510   510   510   510   510 /* T */
/* TA A ..AT */
   510   510   500   380   500 /* N */
   510   510   500   380   500 /* A */
   500   500   500   380   500 /* C */
   380   380   380   380   380 /* G */
   500   500   500   380   500 /* T */
/* TA C ..AT */
   610   530   610   530   510 /* N */
   530   530   530   530   510 /* A */
   610   530   610   530   510 /* C */
   530   530   530   530   510 /* G */
   510   510   510   510   510 /* T */
/* TA G ..AT */
   500   430   500   420   500 /* N */
   430   430   430   420   430 /* A */
   500   430   500   420   500 /* C */
   420   420   420   420   420 /* G */
   500   430   500   420   500 /* T */
/* TA T ..AT */
   530   530   500   530   410 /* N */
   530   530   500   530   410 /* A */
   500   500   500   500   410 /* C */
   530   530   500   530   410 /* G */
   410   410   410   410   410 /* T */
/* TA N ..TA */
   600   560   600   560   560 /* N */
   560   560   560   560   540 /* A */
   600   560   600   560   560 /* C */
   560   560   560   560   500 /* G */
   560   540   560   500   560 /* T */
/* TA A ..TA */
   560   540   560   470   560 /* N */
   540   540   540   470   540 /* A */
   560   540   560   470   560 /* C */
   470   470   470   470   470 /* G */
   560   540   560   470   560 /* T */
/* TA C ..TA */
   600   560   600   560   490 /* N */
   560   560   560   560   490 /* A */
   600   560   600   560   490 /* C */
   560   560   560   560   490 /* G */
   490   490   490   490   490 /* T */
/* TA G ..TA */
   560   470   560   480   560 /* N */
   470   470   470   470   470 /* A */
   560   470   560   480   560 /* C */
   480   470   480   480   480 /* G */
   560   470   560   480   560 /* T */
/* TA T ..TA */
   560   560   490   560   500 /* N */
   560   560   490   560   500 /* A */
   490   490   490   490   490 /* C */
   560   560   490   560   500 /* G */
   500   500   490   500   500 /* T */
/* TA N ..NN */
   610   560   610   560   580 /* N */
   560   560   560   560   550 /* A */
   610   560   610   560   580 /* C */
   560   560   560   560   550 /* G */
   580   550   580   550   580 /* T */
/* TA A ..NN */
   560   540   560   470   560 /* N */
   540   540   540   470   540 /* A */
   560   540   560   470   560 /* C */
   470   470   470   470   470 /* G */
   560   540   560   470   560 /* T */
/* TA C ..NN */
   610   560   610   560   580 /* N */
   560   560   560   560   550 /* A */
   610   560   610   560   580 /* C */
   560   560   560   560   550 /* G */
   580   550   580   550   580 /* T */
/* TA G ..NN */
   560   470   560   480   560 /* N */
   470   470   470   470   470 /* A */
   560   470   560   480   560 /* C */
   480   470   480   480   480 /* G */
   560   470   560   480   560 /* T */
/* TA T ..NN */
   560   560   510   560   500 /* N */
   560   560   510   560   500 /* A */
   510   510   510   510   490 /* C */
   560   560   510   560   500 /* G */
   500   500   490   500   500 /* T */
/* NN N ..CG */
   600   550   600   550   520 /* N */
   550   550   550   550   520 /* A */
   600   550   600   550   520 /* C */
   550   550   550   550   500 /* G */
   520   520   520   500   520 /* T */
/* NN A ..CG */
   520   520   520   420   520 /* N */
   520   520   520   420   520 /* A */
   520   520   520   420   520 /* C */
   420   420   420   420   420 /* G */
   520   520   520   420   520 /* T */
/* NN C ..CG */
   600   550   600   550   500 /* N */
   550   550   550   550   500 /* A */
   600   550   600   550   500 /* C */
   550   550   550   550   500 /* G */
   500   500   500   500   500 /* T */
/* NN G ..CG */
   520   450   520   410   520 /* N */
   450   450   450   410   450 /* A */
   520   450   520   410   520 /* C */
   410   410   410   410   410 /* G */
   520   450   520   410   520 /* T */
/* NN T ..CG */
   550   550   500   550   460 /* N */
   550   550   500   550   460 /* A */
   500   500   500   500   460 /* C */
   550   550   500   550   460 /* G */
   460   460   460   460   460 /* T */
/* NN N ..GC */
   640   550   640   550   580 /* N */
   550   550   550   550   550 /* A */
   640   550   640   550   580 /* C */
   550   550   550   550   550 /* G */
   580   550   580   550   580 /* T */
/* NN A ..GC */
   520   520   510   420   510 /* N */
   520   520   510   420   510 /* A */
   510   510   510   420   510 /* C */
   420   420   420   420   420 /* G */
   510   510   510   420   510 /* T */
/* NN C ..GC */
   640   550   640   550   580 /* N */
   550   550   550   550   550 /* A */
   640   550   640   550   580 /* C */
   550   550   550   550   550 /* G */
   580   550   580   550   580 /* T */
/* NN G ..GC */
   510   450   510   430   510 /* N */
   450   450   450   430   450 /* A */
   510   450   510   430   510 /* C */
   430   430   430   430   430 /* G */
   510   450   510   430   510 /* T */
/* NN T ..GC */
   550   550   540   550   470 /* N */
   550   550   510   550   470 /* A */
   540   510   540   510   470 /* C */
   550   550   510   550   470 /* G */
   470   470   470   470   470 /* T */
/* NN N ..GT */
   640   550   640   550   580 /* N */
   550   550   550   550   550 /* A */
   640   550   640   550   580 /* C */
   550   550   550   550   550 /* G */
   580   550   580   550   580 /* T */
/* NN A ..GT */
   520   520   510   420   510 /* N */
   520   520   510   420   510 /* A */
   510   510   510   420   510 /* C */
   420   420   420   420   420 /* G */
   510   510   510   420   510 /* T */
/* NN C ..GT */
   640   550   640   550   580 /* N */
   550   550   550   550   550 /* A */
   640   550   640   550   580 /* C */
   550   550   550   550   550 /* G */
   580   550   580   550   580 /* T */
/* NN G ..GT */
   510   450   510   430   510 /* N */
   450   450   450   430   450 /* A */
   510   450   510   430   510 /* C */
   430   430   430   430   430 /* G */
   510   450   510   430   510 /* T */
/* NN T ..GT */
   550   550   540   550   470 /* N */
   550   550   510   550   470 /* A */
   540   510   540   510   470 /* C */
   550   550   510   550   470 /* G */
   470   470   470   470   470 /* T */
/* NN N ..TG */
   600   550   600   550   520 /* N */
   550   550   550   550   520 /* A */
   600   550   600   550   520 /* C */
   550   550   550   550   500 /* G */
   520   520   520   500   520 /* T */
/* NN A ..TG */
   520   520   520   420   520 /* N */
   520   520   520   420   520 /* A */
   520   520   520   420   520 /* C */
   420   420   420   420   420 /* G */
   520   520   520   420   520 /* T */
/* NN C ..TG */
   600   550   600   550   500 /* N */
   550   550   550   550   500 /* A */
   600   550   600   550   500 /* C */
   550   550   550   550   500 /* G */
   500   500   500   500   500 /* T */
/* NN G ..TG */
   520   450   520   410   520 /* N */
   450   450   450   410   450 /* A */
   520   450   520   410   520 /* C */
   410   410   410   410   410 /* G */
   520   450   520   410   520 /* T */
/* NN T ..TG */
   550   550   500   550   460 /* N */
   550   550   500   550   460 /* A */
   500   500   500   500   460 /* C */
   550   550   500   550   460 /* G */
   460   460   460   460   460 /* T */
/* NN N ..AT */
   670   530   670   530   540 /* N */
   530   530   530   530   510 /* A */
   670   530   670   530   540 /* C */
   530   530   530   530   510 /* G */
   540   510   540   510   540 /* T */
/* NN A ..AT */
   520   520   500   380   500 /* N */
   520   520   500   380   500 /* A */
   500   500   500   380   500 /* C */
   380   380   380   380   380 /* G */
   500   500   500   380   500 /* T */
/* NN C ..AT */
   670   530   670   530   540 /* N */
   530   530   530   530   510 /* A */
   670   530   670   530   540 /* C */
   530   530   530   530   510 /* G */
   540   510   540   510   540 /* T */
/* NN G ..AT */
   500   430   500   420   500 /* N */
   430   430   430   420   430 /* A */
   500   430   500   420   500 /* C */
   420   420   420   420   420 /* G */
   500   430   500   420   500 /* T */
/* NN T ..AT */
   530   530   500   530   520 /* N */
   530   530   500   530   480 /* A */
   500   500   500   500   470 /* C */
   530   530   500   530   480 /* G */
   520   480   470   480   520 /* T */
/* NN N ..TA */
   610   560   610   560   560 /* N */
   560   560   560   560   540 /* A */
   610   560   610   560   560 /* C */
   560   560   560   560   500 /* G */
   560   540   560   500   560 /* T */
/* NN A ..TA */
   560   540   560   470   560 /* N */
   540   540   540   470   540 /* A */
   560   540   560   470   560 /* C */
   470   470   470   470   470 /* G */
   560   540   560   470   560 /* T */
/* NN C ..TA */
   610   560   610   560   510 /* N */
   560   560   560   560   500 /* A */
   610   560   610   560   510 /* C */
   560   560   560   560   500 /* G */
   510   500   510   500   510 /* T */
/* NN G ..TA */
   560   470   560   480   560 /* N */
   470   470   470   470   470 /* A */
   560   470   560   480   560 /* C */
   480   470   480   480   480 /* G */
   560   470   560   480   560 /* T */
/* NN T ..TA */
   580   560   580   560   500 /* N */
   560   560   500   560   500 /* A */
   580   500   580   500   490 /* C */
   560   560   500   560   500 /* G */
   500   500   490   500   500 /* T */
/* NN N ..NN */
   670   560   670   560   580 /* N */
   560   560   560   560   550 /* A */
   670   560   670   560   580 /* C */
   560   560   560   560   550 /* G */
   580   550   580   550   580 /* T */
/* NN A ..NN */
   560   540   560   470   560 /* N */
   540   540   540   470   540 /* A */
   560   540   560   470   560 /* C */
   470   470   470   470   470 /* G */
   560   540   560   470   560 /* T */
/* NN C ..NN */
   670   560   670   560   580 /* N */
   560   560   560   560   550 /* A */
   670   560   670   560   580 /* C */
   560   560   560   560   550 /* G */
   580   550   580   550   580 /* T */
/* NN G ..NN */
   560   470   560   480   560 /* N */
   470   470   470   470   470 /* A */
   560   470   560   480   560 /* C */
   480   470   480   480   480 /* G */
   560   470   560   480   560 /* T */
/* NN T ..NN */
   580   560   580   560   520 /* N */
   560   560   510   560   500 /* A */
   580   510   580   510   490 /* C */
   560   560   510   560   500 /* G */
   520   500   490   500   520 /* T */

# int21_enthalpies
/* CG N ..CG */
  2500  2500  2500  2500  2500 /* N */
  2500  2500  2500  2500  2500 /* A */
  2500  2500  2400  2500  2400 /* C */
  2500  2500  2500  2500  2500 /* G */
  2500  2500  2400  2500  2400 /* T */
/* CG A ..CG */
  2500  2500  2500  2500  2500 /* N */
  2500  2290  2290  2500  2290 /* A */
  2500  2290  2400  2500  2400 /* C */
  2500  2500  2500  2500  2500 /* G */
  2500  2290  2400  2500  2400 /* T */
/* CG C ..CG */
  2400  2400  2400  2400  2050 /* N */
  2400  2400  2400  2400  2050 /* A */
  2400  2400  2300  2400  2050 /* C */
  2400  2400  2400  2400  2050 /* G */
  2050  2050  2050  2050  2050 /* T */
/* CG G ..CG */
  2500  2500  2500  2500  2500 /* N */
  2500  2500  2500  2500  2500 /* A */
  2500  2500  2400  1900  2400 /* C */
  2500  2500  1900  1900  1900 /* G */
  2500  2500  2400  1900  2400 /* T */
/* CG T ..CG */
  2400  2400  2050  2400  1620 /* N */
  2400  2400  2050  2400  1620 /* A */
  2050  2050  2050  2050  1620 /* C */
  2400  2400  2050  2400  1620 /* G */
  1620  1620  1620  1620  1620 /* T */
/* CG N ..GC */
  2820  2820  2820  2820  2350 /* N */
  2820  2820  2820  2820  2350 /* A */
  2820  2820  2640  2820  2350 /* C */
  2820  2820  2820  2820  2100 /* G */
  2350  2350  2350  2100  2350 /* T */
/* CG A ..GC */
  2350  2350  2350  2020  2350 /* N */
  2350  1830  2350  2020  2350 /* A */
  2350  2350  2350  2020  2350 /* C */
  2020  2020  2020  2020  2020 /* G */
  2350  2350  2350  2020  2350 /* T */
/* CG C ..GC */
  2820  2820  2820  2820  2080 /* N */
  2820  2820  2820  2820  2080 /* A */
  2820  2820  2640  2820  2080 /* C */
  2820  2820  2820  2820  2080 /* G */
  2080  2080  2080  2080  2080 /* T */
/* CG G ..GC */
  2350  2280  2350  1800  2350 /* N */
  2280  2280  2280  1800  2280 /* A */
  2350  2280  2350  1800  2350 /* C */
  1800  1800  1800  1800  1800 /* G */
  2350  2280  2350  1800  2350 /* T */
/* CG T ..GC */
  2820  2820  2100  2820  2100 /* N */
  2820  2820  1940  2820  2100 /* A */
  1940  1940  1940  1940  1940 /* C */
  2820  2820  1940  2820  2100 /* G */
  2100  2100  2100  2100  2100 /* T */
/* CG N ..GT */
  2820  2820  2820  2820  2350 /* N */
  2820  2820  2820  2820  2350 /* A */
  2820  2820  2640  2820  2350 /* C */
  2820  2820  2820  2820  2100 /* G */
  2350  2350  2350  2100  2350 /* T */
/* CG A ..GT */
  2350  2350  2350  2020  2350 /* N */
  2350  1830  2350  2020  2350 /* A */
  2350  2350  2350  2020  2350 /* C */
  2020  2020  2020  2020  2020 /* G */
  2350  2350  2350  2020  2350 /* T */
/* CG C ..GT */
  2820  2820  2820  2820  2080 /* N */
  2820  2820  2820  2820  2080 /* A */
  2820  2820  2640  2820  2080 /* C */
  2820  2820  2820  2820  2080 /* G */
  2080  2080  2080  2080  2080 /* T */
/* CG G ..GT */
  2350  2280  2350  1800  2350 /* N */
  2280  2280  2280  1800  2280 /* A */
  2350  2280  2350  1800  2350 /* C */
  1800  1800  1800  1800  1800 /* G */
  2350  2280  2350  1800  2350 /* T */
/* CG T ..GT */
  2820  2820  2100  2820  2100 /* N */
  2820  2820  1940  2820  2100 /* A */
  1940  1940  1940  1940  1940 /* C */
  2820  2820  1940  2820  2100 /* G */
  2100  2100  2100  2100  2100 /* T */
/* CG N ..TG */
  2500  2500  2500  2500  2500 /* N */
  2500  2500  2500  2500  2500 /* A */
  2500  2500  2400  2500  2400 /* C */
  2500  2500  2500  2500  2500 /* G */
  2500  2500  2400  2500  2400 /* T */
/* CG A ..TG */
  2500  2500  2500  2500  2500 /* N */
  2500  2290  2290  2500  2290 /* A */
  2500  2290  2400  2500  2400 /* C */
  2500  2500  2500  2500  2500 /* G */
  2500  2290  2400  2500  2400 /* T */
/* CG C ..TG */
  2400  2400  2400  2400  2050 /* N */
  2400  2400  2400  2400  2050 /* A */
  2400  2400  2300  2400  2050 /* C */
  2400  2400  2400  2400  2050 /* G */
  2050  2050  2050  2050  2050 /* T */
/* CG G ..TG */
  2500  2500  2500  2500  2500 /* N */
  2500  2500  2500  2500  2500 /* A */
  2500  2500  2400  1900  2400 /* C */
  2500  2500  1900  1900  1900 /* G */
  2500  2500  2400  1900  2400 /* T */
/* CG T ..TG */
  2400  2400  2050  2400  1620 /* N */
  2400  2400  2050  2400  1620 /* A */
  2050  2050  2050  2050  1620 /* C */
  2400  2400  2050  2400  1620 /* G */
  1620  1620  1620  1620  1620 /* T */
/* CG N ..AT */
  2420  2420  2420  2280  2420 /* N */
  2420  2390  2420  2280  2420 /* A */
  2420  2420  2420  2280  2420 /* C */
  2280  2280  2280  2280  1670 /* G */
  2420  2420  2420  1670  2420 /* T */
/* CG A ..AT */
  2420  2420  2420  1670  2420 /* N */
  2420  2390  2420  1670  2420 /* A */
  2420  2420  2420  1670  2420 /* C */
  1670  1670  1670  1670  1670 /* G */
  2420  2420  2420  1670  2420 /* T */
/* CG C ..AT */
  2280  2280  2280  2280   660 /* N */
  2280  2280  2280  2280   660 /* A */
  2280  2280  2150  2280   660 /* C */
  2280  2280  2280  2280   660 /* G */
   660   660   660   660   660 /* T */
/* CG G ..AT */
  2420  1260  2420   900  2420 /* N */
  1260  1260  1260   900  1260 /* A */
  2420  1260  2420   900  2420 /* C */
   900   900   900   900   900 /* G */
  2420  1260  2420   900  2420 /* T */
/* CG T ..AT */
  2280  2280  1400  2280  1140 /* N */
  2280  2280  1400  2280  1140 /* A */
  1400  1400  1400  1400  1140 /* C */
  2280  2280  1400  2280  1140 /* G */
  1140  1140  1140  1140  1140 /* T */
/* CG N ..TA */
  3210  3020  3210  3020  3210 /* N */
  3020  2890  2890  3020  2890 /* A */
  3210  2890  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2890  3210  3020  3210 /* T */
/* CG A ..TA */
  3210  3020  3210  3020  3210 /* N */
  3020  2890  2890  3020  2890 /* A */
  3210  2890  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2890  3210  3020  3210 /* T */
/* CG C ..TA */
  2850  2850  2850  2850  1900 /* N */
  2850  2850  2850  2850  1900 /* A */
  2850  2850  2480  2850  1900 /* C */
  2850  2850  2850  2850  1900 /* G */
  1900  1900  1900  1900  1900 /* T */
/* CG G ..TA */
  3210  2580  3210  1780  3210 /* N */
  2580  2580  2580  1780  2580 /* A */
  3210  2580  3210  1780  3210 /* C */
  1780  1780  1780  1780  1780 /* G */
  3210  2580  3210  1780  3210 /* T */
/* CG T ..TA */
  2850  2850  2190  2850  1960 /* N */
  2850  2850  2190  2850  1960 /* A */
  2190  2190  2190  2190  1960 /* C */
  2850  2850  2190  2850  1960 /* G */
  1960  1960  1960  1960  1960 /* T */
/* CG N ..NN */
  3210  3020  3210  3020  3210 /* N */
  3020  2890  2890  3020  2890 /* A */
  3210  2890  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2890  3210  3020  3210 /* T */
/* CG A ..NN */
  3210  3020  3210  3020  3210 /* N */
  3020  2890  2890  3020  2890 /* A */
  3210  2890  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2890  3210  3020  3210 /* T */
/* CG C ..NN */
  2850  2850  2850  2850  2080 /* N */
  2850  2850  2850  2850  2080 /* A */
  2850  2850  2640  2850  2080 /* C */
  2850  2850  2850  2850  2080 /* G */
  2080  2080  2080  2080  2080 /* T */
/* CG G ..NN */
  3210  2580  3210  2500  3210 /* N */
  2580  2580  2580  2500  2580 /* A */
  3210  2580  3210  1900  3210 /* C */
  2500  2500  1900  1900  1900 /* G */
  3210  2580  3210  1900  3210 /* T */
/* CG T ..NN */
  2850  2850  2190  2850  2100 /* N */
  2850  2850  2190  2850  2100 /* A */
  2190  2190  2190  2190  1960 /* C */
  2850  2850  2190  2850  2100 /* G */
  2100  2100  2100  2100  2100 /* T */
/* GC N ..CG */
  2820  2350  2820  2350  2820 /* N */
  2350  2350  2350  2350  2350 /* A */
  2820  2350  2820  2350  2820 /* C */
  2350  2350  2350  2350  2350 /* G */
  2820  2350  2820  2350  2820 /* T */
/* GC A ..CG */
  2820  2280  2820  2280  2820 /* N */
  2280  1830  1830  2280  1830 /* A */
  2820  1830  2820  2280  2820 /* C */
  2280  2280  2280  2280  2280 /* G */
  2820  1830  2820  2280  2820 /* T */
/* GC C ..CG */
  2640  2350  2640  2350  2350 /* N */
  2350  2350  2350  2350  2350 /* A */
  2640  2350  2640  2350  1940 /* C */
  2350  2350  2350  2350  2350 /* G */
  2350  2350  1940  2350  1940 /* T */
/* GC G ..CG */
  2820  2020  2820  2020  2820 /* N */
  2020  2020  2020  2020  2020 /* A */
  2820  2020  2820  1800  2820 /* C */
  2020  2020  1800  1800  1800 /* G */
  2820  2020  2820  1800  2820 /* T */
/* GC T ..CG */
  2350  2350  2350  2350  2350 /* N */
  2350  2350  2350  2350  2350 /* A */
  2350  2350  2080  2350  2100 /* C */
  2350  2350  2350  2350  2350 /* G */
  2350  2350  2100  2350  2100 /* T */
/* GC N ..GC */
  2490  2280  2490  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2490  2280  2490  2280  2280 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  2280  2280  2280 /* T */
/* GC A ..GC */
  2280  2150  2280  2130  2280 /* N */
  2150  2150  2150  2130  2150 /* A */
  2280  2150  2280  2130  2280 /* C */
  2130  2130  2130  2130  2130 /* G */
  2280  2150  2280  2130  2280 /* T */
/* GC C ..GC */
  2490  2280  2490  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2490  2280  2490  2280  1830 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  1830  2280  1830 /* T */
/* GC G ..GC */
  2280  2130  2280  1610  2280 /* N */
  2130  2130  2130  1610  2130 /* A */
  2280  2130  2280  1610  2280 /* C */
  1610  1610  1610  1610  1610 /* G */
  2280  2130  2280  1610  2280 /* T */
/* GC T ..GC */
  2280  2280  2280  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2280  2280  1830  2280  1830 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  1830  2280  1920 /* T */
/* GC N ..GT */
  2490  2280  2490  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2490  2280  2490  2280  2280 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  2280  2280  2280 /* T */
/* GC A ..GT */
  2280  2150  2280  2130  2280 /* N */
  2150  2150  2150  2130  2150 /* A */
  2280  2150  2280  2130  2280 /* C */
  2130  2130  2130  2130  2130 /* G */
  2280  2150  2280  2130  2280 /* T */
/* GC C ..GT */
  2490  2280  2490  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2490  2280  2490  2280  1830 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  1830  2280  1830 /* T */
/* GC G ..GT */
  2280  2130  2280  1610  2280 /* N */
  2130  2130  2130  1610  2130 /* A */
  2280  2130  2280  1610  2280 /* C */
  1610  1610  1610  1610  1610 /* G */
  2280  2130  2280  1610  2280 /* T */
/* GC T ..GT */
  2280  2280  2280  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2280  2280  1830  2280  1830 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  1830  2280  1920 /* T */
/* GC N ..TG */
  2820  2350  2820  2350  2820 /* N */
  2350  2350  2350  2350  2350 /* A */
  2820  2350  2820  2350  2820 /* C */
  2350  2350  2350  2350  2350 /* G */
  2820  2350  2820  2350  2820 /* T */
/* GC A ..TG */
  2820  2280  2820  2280  2820 /* N */
  2280  1830  1830  2280  1830 /* A */
  2820  1830  2820  2280  2820 /* C */
  2280  2280  2280  2280  2280 /* G */
  2820  1830  2820  2280  2820 /* T */
/* GC C ..TG */
  2640  2350  2640  2350  2350 /* N */
  2350  2350  2350  2350  2350 /* A */
  2640  2350  2640  2350  1940 /* C */
  2350  2350  2350  2350  2350 /* G */
  2350  2350  1940  2350  1940 /* T */
/* GC G ..TG */
  2820  2020  2820  2020  2820 /* N */
  2020  2020  2020  2020  2020 /* A */
  2820  2020  2820  1800  2820 /* C */
  2020  2020  1800  1800  1800 /* G */
  2820  2020  2820  1800  2820 /* T */
/* GC T ..TG */
  2350  2350  2350  2350  2350 /* N */
  2350  2350  2350  2350  2350 /* A */
  2350  2350  2080  2350  2100 /* C */
  2350  2350  2350  2350  2350 /* G */
  2350  2350  2100  2350  2100 /* T */
/* GC N ..AT */
  3170  2980  3170  2980  3050 /* N */
  2980  2980  2980  2980  2980 /* A */
  3170  2980  3170  2980  3050 /* C */
  2980  2980  2980  2980  2980 /* G */
  3050  2980  3050  2980  3050 /* T */
/* GC A ..AT */
  2890  2890  2340  2340  2340 /* N */
  2890  2890  2300  2340  2300 /* A */
  2340  2300  2300  2340  2300 /* C */
  2340  2340  2340  2340  2340 /* G */
  2340  2300  2300  2340  2300 /* T */
/* GC C ..AT */
  3170  2980  3170  2980  3050 /* N */
  2980  2980  2980  2980  2980 /* A */
  3170  2980  3170  2980  3050 /* C */
  2980  2980  2980  2980  2980 /* G */
  3050  2980  3050  2980  3050 /* T */
/* GC G ..AT */
  2780  2780  2780  1390  2780 /* N */
  2780  2780  2780  1390  2780 /* A */
  2780  2780  2300  1390  2300 /* C */
  1390  1390  1390  1390  1390 /* G */
  2780  2780  2300  1390  2300 /* T */
/* GC T ..AT */
  2980  2980  2020  2980  2020 /* N */
  2980  2980  2020  2980  1940 /* A */
  2020  2020  2020  2020  2020 /* C */
  2980  2980  2020  2980  1940 /* G */
  1940  1940  1940  1940  1940 /* T */
/* GC N ..TA */
  3210  3020  3210  3020  3210 /* N */
  2580  2580  2580  2580  2580 /* A */
  3210  2580  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2580  3210  3020  3210 /* T */
/* GC A ..TA */
  3210  3020  3210  3020  3210 /* N */
  1750  1750  1750  1750  1750 /* A */
  3210  1750  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  1750  3210  3020  3210 /* T */
/* GC C ..TA */
  3180  2310  3180  2310  2310 /* N */
  2310  2310  2310  2310  2310 /* A */
  3180  2310  3180  2310  1850 /* C */
  2310  2310  2310  2310  2310 /* G */
  2310  2310  1850  2310  1850 /* T */
/* GC G ..TA */
  3210  2580  3210  2580  3210 /* N */
  2580  2580  2580  2580  2580 /* A */
  3210  2580  3210  2230  3210 /* C */
  2580  2580  2230  2230  2230 /* G */
  3210  2580  3210  2230  3210 /* T */
/* GC T ..TA */
  2310  2310  2310  2310   890 /* N */
  2310  2310  2310  2310   890 /* A */
  2310  2310  2260  2310   890 /* C */
  2310  2310  2310  2310   890 /* G */
   890   890   890   890   890 /* T */
/* GC N ..NN */
  3210  3020  3210  3020  3210 /* N */
  2980  2980  2980  2980  2980 /* A */
  3210  2980  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2980  3210  3020  3210 /* T */
/* GC A ..NN */
  3210  3020  3210  3020  3210 /* N */
  2890  2890  2300  2340  2300 /* A */
  3210  2300  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2300  3210  3020  3210 /* T */
/* GC C ..NN */
  3180  2980  3180  2980  3050 /* N */
  2980  2980  2980  2980  2980 /* A */
  3180  2980  3180  2980  3050 /* C */
  2980  2980  2980  2980  2980 /* G */
  3050  2980  3050  2980  3050 /* T */
/* GC G ..NN */
  3210  2780  3210  2580  3210 /* N */
  2780  2780  2780  2580  2780 /* A */
  3210  2780  3210  2230  3210 /* C */
  2580  2580  2230  2230  2230 /* G */
  3210  2780  3210  2230  3210 /* T */
/* GC T ..NN */
  2980  2980  2350  2980  2350 /* N */
  2980  2980  2350  2980  2350 /* A */
  2350  2350  2260  2350  2100 /* C */
  2980  2980  2350  2980  2350 /* G */
  2350  2350  2100  2350  2100 /* T */
/* GT N ..CG */
  2820  2350  2820  2350  2820 /* N */
  2350  2350  2350  2350  2350 /* A */
  2820  2350  2820  2350  2820 /* C */
  2350  2350  2350  2350  2350 /* G */
  2820  2350  2820  2350  2820 /* T */
/* GT A ..CG */
  2820  2280  2820  2280  2820 /* N */
  2280  1830  1830  2280  1830 /* A */
  2820  1830  2820  2280  2820 /* C */
  2280  2280  2280  2280  2280 /* G */
  2820  1830  2820  2280  2820 /* T */
/* GT C ..CG */
  2640  2350  2640  2350  2350 /* N */
  2350  2350  2350  2350  2350 /* A */
  2640  2350  2640  2350  1940 /* C */
  2350  2350  2350  2350  2350 /* G */
  2350  2350  1940  2350  1940 /* T */
/* GT G ..CG */
  2820  2020  2820  2020  2820 /* N */
  2020  2020  2020  2020  2020 /* A */
  2820  2020  2820  1800  2820 /* C */
  2020  2020  1800  1800  1800 /* G */
  2820  2020  2820  1800  2820 /* T */
/* GT T ..CG */
  2350  2350  2350  2350  2350 /* N */
  2350  2350  2350  2350  2350 /* A */
  2350  2350  2080  2350  2100 /* C */
  2350  2350  2350  2350  2350 /* G */
  2350  2350  2100  2350  2100 /* T */
/* GT N ..GC */
  2490  2280  2490  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2490  2280  2490  2280  2280 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  2280  2280  2280 /* T */
/* GT A ..GC */
  2280  2150  2280  2130  2280 /* N */
  2150  2150  2150  2130  2150 /* A */
  2280  2150  2280  2130  2280 /* C */
  2130  2130  2130  2130  2130 /* G */
  2280  2150  2280  2130  2280 /* T */
/* GT C ..GC */
  2490  2280  2490  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2490  2280  2490  2280  1830 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  1830  2280  1830 /* T */
/* GT G ..GC */
  2280  2130  2280  1610  2280 /* N */
  2130  2130  2130  1610  2130 /* A */
  2280  2130  2280  1610  2280 /* C */
  1610  1610  1610  1610  1610 /* G */
  2280  2130  2280  1610  2280 /* T */
/* GT T ..GC */
  2280  2280  2280  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2280  2280  1830  2280  1830 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  1830  2280  1920 /* T */
/* GT N ..GT */
  2490  2280  2490  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2490  2280  2490  2280  2280 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  2280  2280  2280 /* T */
/* GT A ..GT */
  2280  2150  2280  2130  2280 /* N */
  2150  2150  2150  2130  2150 /* A */
  2280  2150  2280  2130  2280 /* C */
  2130  2130  2130  2130  2130 /* G */
  2280  2150  2280  2130  2280 /* T */
/* GT C ..GT */
  2490  2280  2490  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2490  2280  2490  2280  1830 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  1830  2280  1830 /* T */
/* GT G ..GT */
  2280  2130  2280  1610  2280 /* N */
  2130  2130  2130  1610  2130 /* A */
  2280  2130  2280  1610  2280 /* C */
  1610  1610  1610  1610  1610 /* G */
  2280  2130  2280  1610  2280 /* T */
/* GT T ..GT */
  2280  2280  2280  2280  2280 /* N */
  2280  2280  2280  2280  2280 /* A */
  2280  2280  1830  2280  1830 /* C */
  2280  2280  2280  2280  2280 /* G */
  2280  2280  1830  2280  1920 /* T */
/* GT N ..TG */
  2820  2350  2820  2350  2820 /* N */
  2350  2350  2350  2350  2350 /* A */
  2820  2350  2820  2350  2820 /* C */
  2350  2350  2350  2350  2350 /* G */
  2820  2350  2820  2350  2820 /* T */
/* GT A ..TG */
  2820  2280  2820  2280  2820 /* N */
  2280  1830  1830  2280  1830 /* A */
  2820  1830  2820  2280  2820 /* C */
  2280  2280  2280  2280  2280 /* G */
  2820  1830  2820  2280  2820 /* T */
/* GT C ..TG */
  2640  2350  2640  2350  2350 /* N */
  2350  2350  2350  2350  2350 /* A */
  2640  2350  2640  2350  1940 /* C */
  2350  2350  2350  2350  2350 /* G */
  2350  2350  1940  2350  1940 /* T */
/* GT G ..TG */
  2820  2020  2820  2020  2820 /* N */
  2020  2020  2020  2020  2020 /* A */
  2820  2020  2820  1800  2820 /* C */
  2020  2020  1800  1800  1800 /* G */
  2820  2020  2820  1800  2820 /* T */
/* GT T ..TG */
  2350  2350  2350  2350  2350 /* N */
  2350  2350  2350  2350  2350 /* A */
  2350  2350  2080  2350  2100 /* C */
  2350  2350  2350  2350  2350 /* G */
  2350  2350  2100  2350  2100 /* T */
/* GT N ..AT */
  3170  2980  3170  2980  3050 /* N */
  2980  2980  2980  2980  2980 /* A */
  3170  2980  3170  2980  3050 /* C */
  2980  2980  2980  2980  2980 /* G */
  3050  2980  3050  2980  3050 /* T */
/* GT A ..AT */
  2890  2890  2340  2340  2340 /* N */
  2890  2890  2300  2340  2300 /* A */
  2340  2300  2300  2340  2300 /* C */
  2340  2340  2340  2340  2340 /* G */
  2340  2300  2300  2340  2300 /* T */
/* GT C ..AT */
  3170  2980  3170  2980  3050 /* N */
  2980  2980  2980  2980  2980 /* A */
  3170  2980  3170  2980  3050 /* C */
  2980  2980  2980  2980  2980 /* G */
  3050  2980  3050  2980  3050 /* T */
/* GT G ..AT */
  2780  2780  2780  1390  2780 /* N */
  2780  2780  2780  1390  2780 /* A */
  2780  2780  2300  1390  2300 /* C */
  1390  1390  1390  1390  1390 /* G */
  2780  2780  2300  1390  2300 /* T */
/* GT T ..AT */
  2980  2980  2020  2980  2020 /* N */
  2980  2980  2020  2980  1940 /* A */
  2020  2020  2020  2020  2020 /* C */
  2980  2980  2020  2980  1940 /* G */
  1940  1940  1940  1940  1940 /* T */
/* GT N ..TA */
  3210  3020  3210  3020  3210 /* N */
  2580  2580  2580  2580  2580 /* A */
  3210  2580  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2580  3210  3020  3210 /* T */
/* GT A ..TA */
  3210  3020  3210  3020  3210 /* N */
  1750  1750  1750  1750  1750 /* A */
  3210  1750  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  1750  3210  3020  3210 /* T */
/* GT C ..TA */
  3180  2310  3180  2310  2310 /* N */
  2310  2310  2310  2310  2310 /* A */
  3180  2310  3180  2310  1850 /* C */
  2310  2310  2310  2310  2310 /* G */
  2310  2310  1850  2310  1850 /* T */
/* GT G ..TA */
  3210  2580  3210  2580  3210 /* N */
  2580  2580  2580  2580  2580 /* A */
  3210  2580  3210  2230  3210 /* C */
  2580  2580  2230  2230  2230 /* G */
  3210  2580  3210  2230  3210 /* T */
/* GT T ..TA */
  2310  2310  2310  2310   890 /* N */
  2310  2310  2310  2310   890 /* A */
  2310  2310  2260  2310   890 /* C */
  2310  2310  2310  2310   890 /* G */
   890   890   890   890   890 /* T */
/* GT N ..NN */
  3210  3020  3210  3020  3210 /* N */
  2980  2980  2980  2980  2980 /* A */
  3210  2980  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2980  3210  3020  3210 /* T */
/* GT A ..NN */
  3210  3020  3210  3020  3210 /* N */
  2890  2890  2300  2340  2300 /* A */
  3210  2300  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2300  3210  3020  3210 /* T */
/* GT C ..NN */
  3180  2980  3180  2980  3050 /* N */
  2980  2980  2980  2980  2980 /* A */
  3180  2980  3180  2980  3050 /* C */
  2980  2980  2980  2980  2980 /* G */
  3050  2980  3050  2980  3050 /* T */
/* GT G ..NN */
  3210  2780  3210  2580  3210 /* N */
  2780  2780  2780  2580  2780 /* A */
  3210  2780  3210  2230  3210 /* C */
  2580  2580  2230  2230  2230 /* G */
  3210  2780  3210  2230  3210 /* T */
/* GT T ..NN */
  2980  2980  2350  2980  2350 /* N */
  2980  2980  2350  2980  2350 /* A */
  2350  2350  2260  2350  2100 /* C */
  2980  2980  2350  2980  2350 /* G */
  2350  2350  2100  2350  2100 /* T */
/* TG N ..CG */
  2500  2500  2500  2500  2500 /* N */
  2500  2500  2500  2500  2500 /* A */
  2500  2500  2400  2500  2400 /* C */
  2500  2500  2500  2500  2500 /* G */
  2500  2500  2400  2500  2400 /* T */
/* TG A ..CG */
  2500  2500  2500  2500  2500 /* N */
  2500  2290  2290  2500  2290 /* A */
  2500  2290  2400  2500  2400 /* C */
  2500  2500  2500  2500  2500 /* G */
  2500  2290  2400  2500  2400 /* T */
/* TG C ..CG */
  2400  2400  2400  2400  2050 /* N */
  2400  2400  2400  2400  2050 /* A */
  2400  2400  2300  2400  2050 /* C */
  2400  2400  2400  2400  2050 /* G */
  2050  2050  2050  2050  2050 /* T */
/* TG G ..CG */
  2500  2500  2500  2500  2500 /* N */
  2500  2500  2500  2500  2500 /* A */
  2500  2500  2400  1900  2400 /* C */
  2500  2500  1900  1900  1900 /* G */
  2500  2500  2400  1900  2400 /* T */
/* TG T ..CG */
  2400  2400  2050  2400  1620 /* N */
  2400  2400  2050  2400  1620 /* A */
  2050  2050  2050  2050  1620 /* C */
  2400  2400  2050  2400  1620 /* G */
  1620  1620  1620  1620  1620 /* T */
/* TG N ..GC */
  2820  2820  2820  2820  2350 /* N */
  2820  2820  2820  2820  2350 /* A */
  2820  2820  2640  2820  2350 /* C */
  2820  2820  2820  2820  2100 /* G */
  2350  2350  2350  2100  2350 /* T */
/* TG A ..GC */
  2350  2350  2350  2020  2350 /* N */
  2350  1830  2350  2020  2350 /* A */
  2350  2350  2350  2020  2350 /* C */
  2020  2020  2020  2020  2020 /* G */
  2350  2350  2350  2020  2350 /* T */
/* TG C ..GC */
  2820  2820  2820  2820  2080 /* N */
  2820  2820  2820  2820  2080 /* A */
  2820  2820  2640  2820  2080 /* C */
  2820  2820  2820  2820  2080 /* G */
  2080  2080  2080  2080  2080 /* T */
/* TG G ..GC */
  2350  2280  2350  1800  2350 /* N */
  2280  2280  2280  1800  2280 /* A */
  2350  2280  2350  1800  2350 /* C */
  1800  1800  1800  1800  1800 /* G */
  2350  2280  2350  1800  2350 /* T */
/* TG T ..GC */
  2820  2820  2100  2820  2100 /* N */
  2820  2820  1940  2820  2100 /* A */
  1940  1940  1940  1940  1940 /* C */
  2820  2820  1940  2820  2100 /* G */
  2100  2100  2100  2100  2100 /* T */
/* TG N ..GT */
  2820  2820  2820  2820  2350 /* N */
  2820  2820  2820  2820  2350 /* A */
  2820  2820  2640  2820  2350 /* C */
  2820  2820  2820  2820  2100 /* G */
  2350  2350  2350  2100  2350 /* T */
/* TG A ..GT */
  2350  2350  2350  2020  2350 /* N */
  2350  1830  2350  2020  2350 /* A */
  2350  2350  2350  2020  2350 /* C */
  2020  2020  2020  2020  2020 /* G */
  2350  2350  2350  2020  2350 /* T */
/* TG C ..GT */
  2820  2820  2820  2820  2080 /* N */
  2820  2820  2820  2820  2080 /* A */
  2820  2820  2640  2820  2080 /* C */
  2820  2820  2820  2820  2080 /* G */
  2080  2080  2080  2080  2080 /* T */
/* TG G ..GT */
  2350  2280  2350  1800  2350 /* N */
  2280  2280  2280  1800  2280 /* A */
  2350  2280  2350  1800  2350 /* C */
  1800  1800  1800  1800  1800 /* G */
  2350  2280  2350  1800  2350 /* T */
/* TG T ..GT */
  2820  2820  2100  2820  2100 /* N */
  2820  2820  1940  2820  2100 /* A */
  1940  1940  1940  1940  1940 /* C */
  2820  2820  1940  2820  2100 /* G */
  2100  2100  2100  2100  2100 /* T */
/* TG N ..TG */
  2500  2500  2500  2500  2500 /* N */
  2500  2500  2500  2500  2500 /* A */
  2500  2500  2400  2500  2400 /* C */
  2500  2500  2500  2500  2500 /* G */
  2500  2500  2400  2500  2400 /* T */
/* TG A ..TG */
  2500  2500  2500  2500  2500 /* N */
  2500  2290  2290  2500  2290 /* A */
  2500  2290  2400  2500  2400 /* C */
  2500  2500  2500  2500  2500 /* G */
  2500  2290  2400  2500  2400 /* T */
/* TG C ..TG */
  2400  2400  2400  2400  2050 /* N */
  2400  2400  2400  2400  2050 /* A */
  2400  2400  2300  2400  2050 /* C */
  2400  2400  2400  2400  2050 /* G */
  2050  2050  2050  2050  2050 /* T */
/* TG G ..TG */
  2500  2500  2500  2500  2500 /* N */
  2500  2500  2500  2500  2500 /* A */
  2500  2500  2400  1900  2400 /* C */
  2500  2500  1900  1900  1900 /* G */
  2500  2500  2400  1900  2400 /* T */
/* TG T ..TG */
  2400  2400  2050  2400  1620 /* N */
  2400  2400  2050  2400  1620 /* A */
  2050  2050  2050  2050  1620 /* C */
  2400  2400  2050  2400  1620 /* G */
  1620  1620  1620  1620  1620 /* T */
/* TG N ..AT */
  2420  2420  2420  2280  2420 /* N */
  2420  2390  2420  2280  2420 /* A */
  2420  2420  2420  2280  2420 /* C */
  2280  2280  2280  2280  1670 /* G */
  2420  2420  2420  1670  2420 /* T */
/* TG A ..AT */
  2420  2420  2420  1670  2420 /* N */
  2420  2390  2420  1670  2420 /* A */
  2420  2420  2420  1670  2420 /* C */
  1670  1670  1670  1670  1670 /* G */
  2420  2420  2420  1670  2420 /* T */
/* TG C ..AT */
  2280  2280  2280  2280   660 /* N */
  2280  2280  2280  2280   660 /* A */
  2280  2280  2150  2280   660 /* C */
  2280  2280  2280  2280   660 /* G */
   660   660   660   660   660 /* T */
/* TG G ..AT */
  2420  1260  2420   900  2420 /* N */
  1260  1260  1260   900  1260 /* A */
  2420  1260  2420   900  2420 /* C */
   900   900   900   900   900 /* G */
  2420  1260  2420   900  2420 /* T */
/* TG T ..AT */
  2280  2280  1400  2280  1140 /* N */
  2280  2280  1400  2280  1140 /* A */
  1400  1400  1400  1400  1140 /* C */
  2280  2280  1400  2280  1140 /* G */
  1140  1140  1140  1140  1140 /* T */
/* TG N ..TA */
  3210  3020  3210  3020  3210 /* N */
  3020  2890  2890  3020  2890 /* A */
  3210  2890  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2890  3210  3020  3210 /* T */
/* TG A ..TA */
  3210  3020  3210  3020  3210 /* N */
  3020  2890  2890  3020  2890 /* A */
  3210  2890  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2890  3210  3020  3210 /* T */
/* TG C ..TA */
  2850  2850  2850  2850  1900 /* N */
  2850  2850  2850  2850  1900 /* A */
  2850  2850  2480  2850  1900 /* C */
  2850  2850  2850  2850  1900 /* G */
  1900  1900  1900  1900  1900 /* T */
/* TG G ..TA */
  3210  2580  3210  1780  3210 /* N */
  2580  2580  2580  1780  2580 /* A */
  3210  2580  3210  1780  3210 /* C */
  1780  1780  1780  1780  1780 /* G */
  3210  2580  3210  1780  3210 /* T */
/* TG T ..TA */
  2850  2850  2190  2850  1960 /* N */
  2850  2850  2190  2850  1960 /* A */
  2190  2190  2190  2190  1960 /* C */
  2850  2850  2190  2850  1960 /* G */
  1960  1960  1960  1960  1960 /* T */
/* TG N ..NN */
  3210  3020  3210  3020  3210 /* N */
  3020  2890  2890  3020  2890 /* A */
  3210  2890  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2890  3210  3020  3210 /* T */
/* TG A ..NN */
  3210  3020  3210  3020  3210 /* N */
  3020  2890  2890  3020  2890 /* A */
  3210  2890  3210  3020  3210 /* C */
  3020  3020  3020  3020  3020 /* G */
  3210  2890  3210  3020  3210 /* T */
/* TG C ..NN */
  2850  2850  2850  2850  2080 /* N */
  2850  2850  2850  2850  2080 /* A */
  2850  2850  2640  2850  2080 /* C */
  2850  2850  2850  2850  2080 /* G */
  2080  2080  2080  2080  2080 /* T */
/* TG G ..NN */
  3210  2580  3210  2500  3210 /* N */
  2580  2580  2580  2500  2580 /* A */
  3210  2580  3210  1900  3210 /* C */
  2500  2500  1900  1900  1900 /* G */
  3210  2580  3210  1900  3210 /* T */
/* TG T ..NN */
  2850  2850  2190  2850  2100 /* N */
  2850  2850  2190  2850  2100 /* A */
  2190  2190  2190  2190  1960 /* C */
  2850  2850  2190  2850  2100 /* G */
  2100  2100  2100  2100  2100 /* T */
/* AT N ..CG */
  2420  2420  2420  2420  2390 /* N */
  2420  2420  2420  2420  2390 /* A */
  2420  2420  2280  2420  2280 /* C */
  2420  2420  2420  2420  1400 /* G */
  2280  2280  2280  1400  2280 /* T */
/* AT A ..CG */
  2390  2390  2390  1260  2390 /* N */
  2390  2390  2390  1260  2390 /* A */
  2280  2280  2280  1260  2280 /* C */
  1260  1260  1260  1260  1260 /* G */
  2280  2280  2280  1260  2280 /* T */
/* AT C ..CG */
  2420  2420  2420  2420  1400 /* N */
  2420  2420  2420  2420  1400 /* A */
  2420  2420  2150  2420  1400 /* C */
  2420  2420  2420  2420  1400 /* G */
  1400  1400  1400  1400  1400 /* T */
/* AT G ..CG */
  2280  1670  2280  1670  2280 /* N */
  1670  1670  1670  1670  1670 /* A */
  2280  1670  2280   900  2280 /* C */
  1670  1670   900   900   900 /* G */
  2280  1670  2280   900  2280 /* T */
/* AT T ..CG */
  2420  2420   660  2420  1140 /* N */
  2420  2420   660  2420  1140 /* A */
   660   660   660   660   660 /* C */
  2420  2420   660  2420  1140 /* G */
  1140  1140   660  1140  1140 /* T */
/* AT N ..GC */
  3170  2980  3170  2780  2980 /* N */
  2980  2890  2980  2780  2980 /* A */
  3170  2980  3170  2780  2980 /* C */
  2780  2780  2780  2780  2780 /* G */
  2980  2980  2980  2780  2980 /* T */
/* AT A ..GC */
  2980  2980  2980  2780  2980 /* N */
  2980  2890  2980  2780  2980 /* A */
  2980  2980  2980  2780  2980 /* C */
  2780  2780  2780  2780  2780 /* G */
  2980  2980  2980  2780  2980 /* T */
/* AT C ..GC */
  3170  2300  3170  2300  2020 /* N */
  2300  2300  2300  2300  2020 /* A */
  3170  2300  3170  2300  2020 /* C */
  2300  2300  2300  2300  2020 /* G */
  2020  2020  2020  2020  2020 /* T */
/* AT G ..GC */
  2980  2340  2980  1390  2980 /* N */
  2340  2340  2340  1390  2340 /* A */
  2980  2340  2980  1390  2980 /* C */
  1390  1390  1390  1390  1390 /* G */
  2980  2340  2980  1390  2980 /* T */
/* AT T ..GC */
  3050  2300  3050  2300  1940 /* N */
  2300  2300  2300  2300  1940 /* A */
  3050  2300  3050  2300  1940 /* C */
  2300  2300  2300  2300  1940 /* G */
  1940  1940  1940  1940  1940 /* T */
/* AT N ..GT */
  3170  2980  3170  2780  2980 /* N */
  2980  2890  2980  2780  2980 /* A */
  3170  2980  3170  2780  2980 /* C */
  2780  2780  2780  2780  2780 /* G */
  2980  2980  2980  2780  2980 /* T */
/* AT A ..GT */
  2980  2980  2980  2780  2980 /* N */
  2980  2890  2980  2780  2980 /* A */
  2980  2980  2980  2780  2980 /* C */
  2780  2780  2780  2780  2780 /* G */
  2980  2980  2980  2780  2980 /* T */
/* AT C ..GT */
  3170  2300  3170  2300  2020 /* N */
  2300  2300  2300  2300  2020 /* A */
  3170  2300  3170  2300  2020 /* C */
  2300  2300  2300  2300  2020 /* G */
  2020  2020  2020  2020  2020 /* T */
/* AT G ..GT */
  2980  2340  2980  1390  2980 /* N */
  2340  2340  2340  1390  2340 /* A */
  2980  2340  2980  1390  2980 /* C */
  1390  1390  1390  1390  1390 /* G */
  2980  2340  2980  1390  2980 /* T */
/* AT T ..GT */
  3050  2300  3050  2300  1940 /* N */
  2300  2300  2300  2300  1940 /* A */
  3050  2300  3050  2300  1940 /* C */
  2300  2300  2300  2300  1940 /* G */
  1940  1940  1940  1940  1940 /* T */
/* AT N ..TG */
  2420  2420  2420  2420  2390 /* N */
  2420  2420  2420  2420  2390 /* A */
  2420  2420  2280  2420  2280 /* C */
  2420  2420  2420  2420  1400 /* G */
  2280  2280  2280  1400  2280 /* T */
/* AT A ..TG */
  2390  2390  2390  1260  2390 /* N */
  2390  2390  2390  1260  2390 /* A */
  2280  2280  2280  1260  2280 /* C */
  1260  1260  1260  1260  1260 /* G */
  2280  2280  2280  1260  2280 /* T */
/* AT C ..TG */
  2420  2420  2420  2420  1400 /* N */
  2420  2420  2420  2420  1400 /* A */
  2420  2420  2150  2420  1400 /* C */
  2420  2420  2420  2420  1400 /* G */
  1400  1400  1400  1400  1400 /* T */
/* AT G ..TG */
  2280  1670  2280  1670  2280 /* N */
  1670  1670  1670  1670  1670 /* A */
  2280  1670  2280   900  2280 /* C */
  1670  1670   900   900   900 /* G */
  2280  1670  2280   900  2280 /* T */
/* AT T ..TG */
  2420  2420   660  2420  1140 /* N */
  2420  2420   660  2420  1140 /* A */
   660   660   660   660   660 /* C */
  2420  2420   660  2420  1140 /* G */
  1140  1140   660  1140  1140 /* T */
/* AT N ..AT */
  2870  2870  2870  2870  2870 /* N */
  2870  2870  2870  2870  2870 /* A */
  2870  2870  2870  2870  2870 /* C */
  2870  2870  2870  2870  2870 /* G */
  2870  2870  2870  2870  2870 /* T */
/* AT A ..AT */
  2870  2870  2870  2600  2870 /* N */
  2870  2380  2870  2600  2380 /* A */
  2870  2870  2870  2600  2870 /* C */
  2600  2600  2600  2600  2600 /* G */
  2870  2870  2870  2600  2870 /* T */
/* AT C ..AT */
  2870  2870  2870  2870  2320 /* N */
  2870  2870  2870  2870  2320 /* A */
  2870  2870  2640  2870  2320 /* C */
  2870  2870  2870  2870  2320 /* G */
  2320  2320  2320  2320  2320 /* T */
/* AT G ..AT */
  2870  2600  2870  2600  2870 /* N */
  2600  2600  2600  2600  2600 /* A */
  2870  2600  2870  2300  2870 /* C */
  2600  2600  2300  2300  2300 /* G */
  2870  2600  2870  2300  2870 /* T */
/* AT T ..AT */
  2870  2870  2320  2870  2870 /* N */
  2870  2870  2320  2870  2870 /* A */
  2320  2320  2320  2320  2320 /* C */
  2870  2870  2320  2870  2870 /* G */
  2870  2870  2320  2870  2460 /* T */
/* AT N ..TA */
  3360  3180  3360  2960  3180 /* N */
  3180  3180  3180  2960  3180 /* A */
  3360  3180  3360  2960  2040 /* C */
  2960  2960  2960  2960  2960 /* G */
  3180  3180  2040  2960  2040 /* T */
/* AT A ..TA */
  3180  3180  3180  2960  3180 /* N */
  3180  3180  3180  2960  3180 /* A */
  3180  3180  2040  2960  2040 /* C */
  2960  2960  2960  2960  2960 /* G */
  3180  3180  2040  2960  2040 /* T */
/* AT C ..TA */
  3360  2930  3360  2930  2930 /* N */
  2930  2930  2930  2930  2930 /* A */
  3360  2930  3360  2930  1730 /* C */
  2930  2930  2930  2930  2930 /* G */
  1730  1730  1730  1730  1730 /* T */
/* AT G ..TA */
  2070  2070  2070  2070  2070 /* N */
  2070  2070  2070  2070  2070 /* A */
  2070  2070  2040  2050  2040 /* C */
  2070  2070  2050  2050  2050 /* G */
  2070  2070  2040  2050  2040 /* T */
/* AT T ..TA */
  2930  2930  2930  2930  1410 /* N */
  2930  2930  2930  2930  1410 /* A */
  2930  2930  2240  2930  1410 /* C */
  2930  2930  2930  2930  1410 /* G */
  1410  1410  1410  1410  1410 /* T */
/* AT N ..NN */
  3360  3180  3360  2960  3180 /* N */
  3180  3180  3180  2960  3180 /* A */
  3360  3180  3360  2960  2980 /* C */
  2960  2960  2960  2960  2960 /* G */
  3180  3180  2980  2960  2980 /* T */
/* AT A ..NN */
  3180  3180  3180  2960  3180 /* N */
  3180  3180  3180  2960  3180 /* A */
  3180  3180  2980  2960  2980 /* C */
  2960  2960  2960  2960  2960 /* G */
  3180  3180  2980  2960  2980 /* T */
/* AT C ..NN */
  3360  2930  3360  2930  2930 /* N */
  2930  2930  2930  2930  2930 /* A */
  3360  2930  3360  2930  2320 /* C */
  2930  2930  2930  2930  2930 /* G */
  2320  2320  2320  2320  2320 /* T */
/* AT G ..NN */
  2980  2600  2980  2600  2980 /* N */
  2600  2600  2600  2600  2600 /* A */
  2980  2600  2980  2300  2980 /* C */
  2600  2600  2300  2300  2300 /* G */
  2980  2600  2980  2300  2980 /* T */
/* AT T ..NN */
  3050  2930  3050  2930  2870 /* N */
  2930  2930  2930  2930  2870 /* A */
  3050  2930  3050  2930  2320 /* C */
  2930  2930  2930  2930  2870 /* G */
  2870  2870  2320  2870  2460 /* T */
/* TA N ..CG */
  3210  3210  3210  3210  3020 /* N */
  3210  3210  3210  3210  3020 /* A */
  3210  3210  2850  3210  2850 /* C */
  3210  3210  3210  3210  2580 /* G */
  3020  3020  2850  2580  2850 /* T */
/* TA A ..CG */
  2890  2890  2890  2580  2890 /* N */
  2890  2890  2890  2580  2890 /* A */
  2850  2850  2850  2580  2850 /* C */
  2580  2580  2580  2580  2580 /* G */
  2850  2850  2850  2580  2850 /* T */
/* TA C ..CG */
  3210  3210  3210  3210  2190 /* N */
  3210  3210  3210  3210  2190 /* A */
  3210  3210  2480  3210  2190 /* C */
  3210  3210  3210  3210  2190 /* G */
  2190  2190  2190  2190  2190 /* T */
/* TA G ..CG */
  3020  3020  3020  1780  3020 /* N */
  3020  3020  3020  1780  3020 /* A */
  3020  3020  2850  1780  2850 /* C */
  1780  1780  1780  1780  1780 /* G */
  3020  3020  2850  1780  2850 /* T */
/* TA T ..CG */
  3210  3210  1960  3210  1960 /* N */
  3210  3210  1900  3210  1960 /* A */
  1960  1900  1900  1900  1960 /* C */
  3210  3210  1900  3210  1960 /* G */
  1960  1960  1960  1960  1960 /* T */
/* TA N ..GC */
  3210  3210  3210  3210  3210 /* N */
  3210  3210  3210  3210  3210 /* A */
  3210  3210  3180  3210  3180 /* C */
  3210  3210  3210  3210  3210 /* G */
  3210  3210  2310  3210  2310 /* T */
/* TA A ..GC */
  2580  2580  2580  2580  2580 /* N */
  2580  1750  1750  2580  1750 /* A */
  2580  1750  2310  2580  2310 /* C */
  2580  2580  2580  2580  2580 /* G */
  2580  1750  2310  2580  2310 /* T */
/* TA C ..GC */
  3210  3210  3210  3210  3210 /* N */
  3210  3210  3210  3210  3210 /* A */
  3210  3210  3180  3210  3180 /* C */
  3210  3210  3210  3210  3210 /* G */
  3210  3210  2260  3210  2260 /* T */
/* TA G ..GC */
  3020  3020  3020  2230  3020 /* N */
  3020  3020  3020  2230  3020 /* A */
  3020  3020  2310  2230  2310 /* C */
  2230  2230  2230  2230  2230 /* G */
  3020  3020  2310  2230  2310 /* T */
/* TA T ..GC */
  3210  3210  1850  3210   890 /* N */
  3210  3210  1850  3210   890 /* A */
  1850  1850  1850  1850   890 /* C */
  3210  3210  1850  3210   890 /* G */
   890   890   890   890   890 /* T */
/* TA N ..GT */
  3210  3210  3210  3210  3210 /* N */
  3210  3210  3210  3210  3210 /* A */
  3210  3210  3180  3210  3180 /* C */
  3210  3210  3210  3210  3210 /* G */
  3210  3210  2310  3210  2310 /* T */
/* TA A ..GT */
  2580  2580  2580  2580  2580 /* N */
  2580  1750  1750  2580  1750 /* A */
  2580  1750  2310  2580  2310 /* C */
  2580  2580  2580  2580  2580 /* G */
  2580  1750  2310  2580  2310 /* T */
/* TA C ..GT */
  3210  3210  3210  3210  3210 /* N */
  3210  3210  3210  3210  3210 /* A */
  3210  3210  3180  3210  3180 /* C */
  3210  3210  3210  3210  3210 /* G */
  3210  3210  2260  3210  2260 /* T */
/* TA G ..GT */
  3020  3020  3020  2230  3020 /* N */
  3020  3020  3020  2230  3020 /* A */
  3020  3020  2310  2230  2310 /* C */
  2230  2230  2230  2230  2230 /* G */
  3020  3020  2310  2230  2310 /* T */
/* TA T ..GT */
  3210  3210  1850  3210   890 /* N */
  3210  3210  1850  3210   890 /* A */
  1850  1850  1850  1850   890 /* C */
  3210  3210  1850  3210   890 /* G */
   890   890   890   890   890 /* T */
/* TA N ..TG */
  3210  3210  3210  3210  3020 /* N */
  3210  3210  3210  3210  3020 /* A */
  3210  3210  2850  3210  2850 /* C */
  3210  3210  3210  3210  2580 /* G */
  3020  3020  2850  2580  2850 /* T */
/* TA A ..TG */
  2890  2890  2890  2580  2890 /* N */
  2890  2890  2890  2580  2890 /* A */
  2850  2850  2850  2580  2850 /* C */
  2580  2580  2580  2580  2580 /* G */
  2850  2850  2850  2580  2850 /* T */
/* TA C ..TG */
  3210  3210  3210  3210  2190 /* N */
  3210  3210  3210  3210  2190 /* A */
  3210  3210  2480  3210  2190 /* C */
  3210  3210  3210  3210  2190 /* G */
  2190  2190  2190  2190  2190 /* T */
/* TA G ..TG */
  3020  3020  3020  1780  3020 /* N */
  3020  3020  3020  1780  3020 /* A */
  3020  3020  2850  1780  2850 /* C */
  1780  1780  1780  1780  1780 /* G */
  3020  3020  2850  1780  2850 /* T */
/* TA T ..TG */
  3210  3210  1960  3210  1960 /* N */
  3210  3210  1900  3210  1960 /* A */
  1960  1900  1900  1900  1960 /* C */
  3210  3210  1900  3210  1960 /* G */
  1960  1960  1960  1960  1960 /* T */
/* TA N ..AT */
  3360  3180  3360  2240  2960 /* N */
  3180  3180  2960  2070  2960 /* A */
  3360  2960  3360  2070  2930 /* C */
  2240  2070  2070  2070  2240 /* G */
  2960  2960  2930  2240  2930 /* T */
/* TA A ..AT */
  3180  3180  2930  2070  2930 /* N */
  3180  3180  2930  2070  2930 /* A */
  2930  2930  2930  2070  2930 /* C */
  2070  2070  2070  2070  2070 /* G */
  2930  2930  2930  2070  2930 /* T */
/* TA C ..AT */
  3360  2240  3360  2240  2240 /* N */
  2240  2040  2040  2040  2240 /* A */
  3360  2040  3360  2040  2240 /* C */
  2240  2040  2040  2040  2240 /* G */
  2240  2240  2240  2240  2240 /* T */
/* TA G ..AT */
  2960  2960  2960  2050  2960 /* N */
  2960  2960  2960  2050  2960 /* A */
  2960  2960  2930  2050  2930 /* C */
  2050  2050  2050  2050  2050 /* G */
  2960  2960  2930  2050  2930 /* T */
/* TA T ..AT */
  2040  2040  1730  2040  1410 /* N */
  2040  2040  1730  2040  1410 /* A */
  1730  1730  1730  1730  1410 /* C */
  2040  2040  1730  2040  1410 /* G */
  1410  1410  1410  1410  1410 /* T */
/* TA N ..TA */
  3630  3320  3630  3320  3320 /* N */
  3320  3320  3320  3320  3100 /* A */
  3630  3320  3630  3320  3320 /* C */
  3320  3320  3320  3320  3080 /* G */
  3320  3100  3320  3080  3320 /* T */
/* TA A ..TA */
  3320  3100  3320  3080  3320 /* N */
  3100  3100  3100  3080  3100 /* A */
  3320  3100  3320  3080  3320 /* C */
  3080  3080  3080  3080  3080 /* G */
  3320  3100  3320  3080  3320 /* T */
/* TA C ..TA */
  3630  3320  3630  3320  2140 /* N */
  3320  3320  3320  3320  2140 /* A */
  3630  3320  3630  3320  2140 /* C */
  3320  3320  3320  3320  2140 /* G */
  2140  2140  2140  2140  2140 /* T */
/* TA G ..TA */
  3320  3080  3320  3080  3320 /* N */
  3080  3080  3080  3080  3080 /* A */
  3320  3080  3320  2730  3320 /* C */
  3080  3080  2730  2730  2730 /* G */
  3320  3080  3320  2730  3320 /* T */
/* TA T ..TA */
  3320  3320  2140  3320  2510 /* N */
  3320  3320  2140  3320  2510 /* A */
  2140  2140  2140  2140  2140 /* C */
  3320  3320  2140  3320  2510 /* G */
  2510  2510  2140  2510  2510 /* T */
/* TA N ..NN */
  3630  3320  3630  3320  3320 /* N */
  3320  3320  3320  3320  3210 /* A */
  3630  3320  3630  3320  3320 /* C */
  3320  3320  3320  3320  3210 /* G */
  3320  3210  3320  3210  3320 /* T */
/* TA A ..NN */
  3320  3180  3320  3080  3320 /* N */
  3180  3180  3100  3080  3100 /* A */
  3320  3100  3320  3080  3320 /* C */
  3080  3080  3080  3080  3080 /* G */
  3320  3100  3320  3080  3320 /* T */
/* TA C ..NN */
  3630  3320  3630  3320  3210 /* N */
  3320  3320  3320  3320  3210 /* A */
  3630  3320  3630  3320  3180 /* C */
  3320  3320  3320  3320  3210 /* G */
  3210  3210  2260  3210  2260 /* T */
/* TA G ..NN */
  3320  3080  3320  3080  3320 /* N */
  3080  3080  3080  3080  3080 /* A */
  3320  3080  3320  2730  3320 /* C */
  3080  3080  2730  2730  2730 /* G */
  3320  3080  3320  2730  3320 /* T */
/* TA T ..NN */
  3320  3320  2140  3320  2510 /* N */
  3320  3320  2140  3320  2510 /* A */
  2140  2140  2140  2140  2140 /* C */
  3320  3320  2140  3320  2510 /* G */
  2510  2510  2140  2510  2510 /* T */
/* NN N ..CG */
  3210  3210  3210  3210  3020 /* N */
  3210  3210  3210  3210  3020 /* A */
  3210  3210  2850  3210  2850 /* C */
  3210  3210  3210  3210  2580 /* G */
  3020  3020  2850  2580  2850 /* T */
/* NN A ..CG */
  2890  2890  2890  2580  2890 /* N */
  2890  2890  2890  2580  2890 /* A */
  2850  2850  2850  2580  2850 /* C */
  2580  2580  2580  2580  2580 /* G */
  2850  2850  2850  2580  2850 /* T */
/* NN C ..CG */
  3210  3210  3210  3210  2350 /* N */
  3210  3210  3210  3210  2350 /* A */
  3210  3210  2640  3210  2190 /* C */
  3210  3210  3210  3210  2350 /* G */
  2350  2350  2190  2350  2190 /* T */
/* NN G ..CG */
  3020  3020  3020  2500  3020 /* N */
  3020  3020  3020  2500  3020 /* A */
  3020  3020  2850  1900  2850 /* C */
  2500  2500  1900  1900  1900 /* G */
  3020  3020  2850  1900  2850 /* T */
/* NN T ..CG */
  3210  3210  2350  3210  2350 /* N */
  3210  3210  2350  3210  2350 /* A */
  2350  2350  2080  2350  2100 /* C */
  3210  3210  2350  3210  2350 /* G */
  2350  2350  2100  2350  2100 /* T */
/* NN N ..GC */
  3210  3210  3210  3210  3210 /* N */
  3210  3210  3210  3210  3210 /* A */
  3210  3210  3180  3210  3180 /* C */
  3210  3210  3210  3210  3210 /* G */
  3210  3210  2980  3210  2980 /* T */
/* NN A ..GC */
  2980  2980  2980  2780  2980 /* N */
  2980  2890  2980  2780  2980 /* A */
  2980  2980  2980  2780  2980 /* C */
  2780  2780  2780  2780  2780 /* G */
  2980  2980  2980  2780  2980 /* T */
/* NN C ..GC */
  3210  3210  3210  3210  3210 /* N */
  3210  3210  3210  3210  3210 /* A */
  3210  3210  3180  3210  3180 /* C */
  3210  3210  3210  3210  3210 /* G */
  3210  3210  2260  3210  2260 /* T */
/* NN G ..GC */
  3020  3020  3020  2230  3020 /* N */
  3020  3020  3020  2230  3020 /* A */
  3020  3020  2980  2230  2980 /* C */
  2230  2230  2230  2230  2230 /* G */
  3020  3020  2980  2230  2980 /* T */
/* NN T ..GC */
  3210  3210  3050  3210  2280 /* N */
  3210  3210  2300  3210  2280 /* A */
  3050  2300  3050  2300  1940 /* C */
  3210  3210  2300  3210  2280 /* G */
  2280  2280  2100  2280  2100 /* T */
/* NN N ..GT */
  3210  3210  3210  3210  3210 /* N */
  3210  3210  3210  3210  3210 /* A */
  3210  3210  3180  3210  3180 /* C */
  3210  3210  3210  3210  3210 /* G */
  3210  3210  2980  3210  2980 /* T */
/* NN A ..GT */
  2980  2980  2980  2780  2980 /* N */
  2980  2890  2980  2780  2980 /* A */
  2980  2980  2980  2780  2980 /* C */
  2780  2780  2780  2780  2780 /* G */
  2980  2980  2980  2780  2980 /* T */
/* NN C ..GT */
  3210  3210  3210  3210  3210 /* N */
  3210  3210  3210  3210  3210 /* A */
  3210  3210  3180  3210  3180 /* C */
  3210  3210  3210  3210  3210 /* G */
  3210  3210  2260  3210  2260 /* T */
/* NN G ..GT */
  3020  3020  3020  2230  3020 /* N */
  3020  3020  3020  2230  3020 /* A */
  3020  3020  2980  2230  2980 /* C */
  2230  2230  2230  2230  2230 /* G */
  3020  3020  2980  2230  2980 /* T */
/* NN T ..GT */
  3210  3210  3050  3210  2280 /* N */
  3210  3210  2300  3210  2280 /* A */
  3050  2300  3050  2300  1940 /* C */
  3210  3210  2300  3210  2280 /* G */
  2280  2280  2100  2280  2100 /* T */
/* NN N ..TG */
  3210  3210  3210  3210  3020 /* N */
  3210  3210  3210  3210  3020 /* A */
  3210  3210  2850  3210  2850 /* C */
  3210  3210  3210  3210  2580 /* G */
  3020  3020  2850  2580  2850 /* T */
/* NN A ..TG */
  2890  2890  2890  2580  2890 /* N */
  2890  2890  2890  2580  2890 /* A */
  2850  2850  2850  2580  2850 /* C */
  2580  2580  2580  2580  2580 /* G */
  2850  2850  2850  2580  2850 /* T */
/* NN C ..TG */
  3210  3210  3210  3210  2350 /* N */
  3210  3210  3210  3210  2350 /* A */
  3210  3210  2640  3210  2190 /* C */
  3210  3210  3210  3210  2350 /* G */
  2350  2350  2190  2350  2190 /* T */
/* NN G ..TG */
  3020  3020  3020  2500  3020 /* N */
  3020  3020  3020  2500  3020 /* A */
  3020  3020  2850  1900  2850 /* C */
  2500  2500  1900  1900  1900 /* G */
  3020  3020  2850  1900  2850 /* T */
/* NN T ..TG */
  3210  3210  2350  3210  2350 /* N */
  3210  3210  2350  3210  2350 /* A */
  2350  2350  2080  2350  2100 /* C */
  3210  3210  2350  3210  2350 /* G */
  2350  2350  2100  2350  2100 /* T */
/* NN N ..AT */
  3360  3180  3360  2980  3050 /* N */
  3180  3180  2980  2980  2980 /* A */
  3360  2980  3360  2980  3050 /* C */
  2980  2980  2980  2980  2980 /* G */
  3050  2980  3050  2980  3050 /* T */
/* NN A ..AT */
  3180  3180  2930  2600  2930 /* N */
  3180  3180  2930  2600  2930 /* A */
  2930  2930  2930  2600  2930 /* C */
  2600  2600  2600  2600  2600 /* G */
  2930  2930  2930  2600  2930 /* T */
/* NN C ..AT */
  3360  2980  3360  2980  3050 /* N */
  2980  2980  2980  2980  2980 /* A */
  3360  2980  3360  2980  3050 /* C */
  2980  2980  2980  2980  2980 /* G */
  3050  2980  3050  2980  3050 /* T */
/* NN G ..AT */
  2960  2960  2960  2600  2960 /* N */
  2960  2960  2960  2600  2960 /* A */
  2960  2960  2930  2300  2930 /* C */
  2600  2600  2300  2300  2300 /* G */
  2960  2960  2930  2300  2930 /* T */
/* NN T ..AT */
  2980  2980  2320  2980  2870 /* N */
  2980  2980  2320  2980  2870 /* A */
  2320  2320  2320  2320  2320 /* C */
  2980  2980  2320  2980  2870 /* G */
  2870  2870  2320  2870  2460 /* T */
/* NN N ..TA */
  3630  3320  3630  3320  3320 /* N */
  3320  3320  3320  3320  3180 /* A */
  3630  3320  3630  3320  3320 /* C */
  3320  3320  3320  3320  3080 /* G */
  3320  3180  3320  3080  3320 /* T */
/* NN A ..TA */
  3320  3180  3320  3080  3320 /* N */
  3180  3180  3180  3080  3180 /* A */
  3320  3180  3320  3080  3320 /* C */
  3080  3080  3080  3080  3080 /* G */
  3320  3180  3320  3080  3320 /* T */
/* NN C ..TA */
  3630  3320  3630  3320  2930 /* N */
  3320  3320  3320  3320  2930 /* A */
  3630  3320  3630  3320  2140 /* C */
  3320  3320  3320  3320  2930 /* G */
  2310  2310  2140  2310  2140 /* T */
/* NN G ..TA */
  3320  3080  3320  3080  3320 /* N */
  3080  3080  3080  3080  3080 /* A */
  3320  3080  3320  2730  3320 /* C */
  3080  3080  2730  2730  2730 /* G */
  3320  3080  3320  2730  3320 /* T */
/* NN T ..TA */
  3320  3320  2930  3320  2510 /* N */
  3320  3320  2930  3320  2510 /* A */
  2930  2930  2260  2930  2140 /* C */
  3320  3320  2930  3320  2510 /* G */
  2510  2510  2140  2510  2510 /* T */
/* NN N ..NN */
  3630  3320  3630  3320  3320 /* N */
  3320  3320  3320  3320  3210 /* A */
  3630  3320  3630  3320  3320 /* C */
  3320  3320  3320  3320  3210 /* G */
  3320  3210  3320  3210  3320 /* T */
/* NN A ..NN */
  3320  3180  3320  3080  3320 /* N */
  3180  3180  3180  3080  3180 /* A */
  3320  3180  3320  3080  3320 /* C */
  3080  3080  3080  3080  3080 /* G */
  3320  3180  3320  3080  3320 /* T */
/* NN C ..NN */
  3630  3320  3630  3320  3210 /* N */
  3320  3320  3320  3320  3210 /* A */
  3630  3320  3630  3320  3180 /* C */
  3320  3320  3320  3320  3210 /* G */
  3210  3210  3050  3210  3050 /* T */
/* NN G ..NN */
  3320  3080  3320  3080  3320 /* N */
  3080  3080  3080  3080  3080 /* A */
  3320  3080  3320  2730  3320 /* C */
  3080  3080  2730  2730  2730 /* G */
  3320  3080  3320  2730  3320 /* T */
/* NN T ..NN */
  3320  3320  3050  3320  2870 /* N */
  3320  3320  2930  3320  2870 /* A */
  3050  2930  3050  2930  2320 /* C */
  3320  3320  2930  3320  2870 /* G */
  2870  2870  2320  2870  2510 /* T */

# int22_energies
/* CG A A ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG A C ..CG */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG A G ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG A T ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG C A ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG C C ..CG */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG C G ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG C T ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG G A ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG G C ..CG */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG G G ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG G T ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG T A ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG T C ..CG */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG T G ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG T T ..CG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG A A ..GC */
    90   110    90    90 /* A */
    90   110    90    90 /* C */
    90   110    90    90 /* G */
    90   110    90    90 /* T */
/* CG A C ..GC */
   140   160   140   140 /* A */
   140   160   140   140 /* C */
   140   160   140   140 /* G */
   140   160   140   140 /* T */
/* CG A G ..GC */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG A T ..GC */
   110   130   110   110 /* A */
   110   130   110   110 /* C */
   110   130   110   110 /* G */
   110   130   110   110 /* T */
/* CG C A ..GC */
    90   110    90    90 /* A */
    90   110    90    90 /* C */
    90   110    90    90 /* G */
    90   110    90    90 /* T */
/* CG C C ..GC */
   140   160   140   140 /* A */
   140   160   140   140 /* C */
   140   160   140   140 /* G */
   140   160   140   140 /* T */
/* CG C G ..GC */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG C T ..GC */
   110   130   110   110 /* A */
   110   130   110   110 /* C */
   110   130   110   110 /* G */
   110   130   110   110 /* T */
/* CG G A ..GC */
    90   110    90    90 /* A */
    90   110    90    90 /* C */
    90   110    90    90 /* G */
    90   110    90    90 /* T */
/* CG G C ..GC */
   140   160   140   140 /* A */
   140   160   140   140 /* C */
   140   160   140   140 /* G */
   140   160   140   140 /* T */
/* CG G G ..GC */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG G T ..GC */
   110   130   110   110 /* A */
   110   130   110   110 /* C */
   110   130   110   110 /* G */
   110   130   110   110 /* T */
/* CG T A ..GC */
    90   110    90    90 /* A */
    90   110    90    90 /* C */
    90   110    90    90 /* G */
    90   110    90    90 /* T */
/* CG T C ..GC */
   140   160   140   140 /* A */
   140   160   140   140 /* C */
   140   160   140   140 /* G */
   140   160   140   140 /* T */
/* CG T G ..GC */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG T T ..GC */
   110   130   110   110 /* A */
   110   130   110   110 /* C */
   110   130   110   110 /* G */
   110   130   110   110 /* T */
/* CG A A ..GT */
   110   130   110   110 /* A */
   110   130   110   110 /* C */
   110   130   110   110 /* G */
   110   130   110   110 /* T */
/* CG A C ..GT */
   160   180   160   160 /* A */
   160   180   160   160 /* C */
   160   180   160   160 /* G */
   160   180   160   160 /* T */
/* CG A G ..GT */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG A T ..GT */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG C A ..GT */
   110   130   110   110 /* A */
   110   130   110   110 /* C */
   110   130   110   110 /* G */
   110   130   110   110 /* T */
/* CG C C ..GT */
   160   180   160   160 /* A */
   160   180   160   160 /* C */
   160   180   160   160 /* G */
   160   180   160   160 /* T */
/* CG C G ..GT */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG C T ..GT */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG G A ..GT */
   110   130   110   110 /* A */
   110   130   110   110 /* C */
   110   130   110   110 /* G */
   110   130   110   110 /* T */
/* CG G C ..GT */
   160   180   160   160 /* A */
   160   180   160   160 /* C */
   160   180   160   160 /* G */
   160   180   160   160 /* T */
/* CG G G ..GT */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG G T ..GT */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG T A ..GT */
   110   130   110   110 /* A */
   110   130   110   110 /* C */
   110   130   110   110 /* G */
   110   130   110   110 /* T */
/* CG T C ..GT */
   160   180   160   160 /* A */
   160   180   160   160 /* C */
   160   180   160   160 /* G */
   160   180   160   160 /* T */
/* CG T G ..GT */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG T T ..GT */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG A A ..TG */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG A C ..TG */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG A G ..TG */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG A T ..TG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG C A ..TG */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG C C ..TG */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG C G ..TG */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG C T ..TG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG G A ..TG */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG G C ..TG */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG G G ..TG */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG G T ..TG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG T A ..TG */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG T C ..TG */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG T G ..TG */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG T T ..TG */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG A A ..AT */
   180   200   180   180 /* A */
   180   200   180   180 /* C */
   180   200   180   180 /* G */
   180   200   180   180 /* T */
/* CG A C ..AT */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG A G ..AT */
   280   300   280   280 /* A */
   280   300   280   280 /* C */
   280   300   280   280 /* G */
   280   300   280   280 /* T */
/* CG A T ..AT */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG C A ..AT */
   180   200   180   180 /* A */
   180   200   180   180 /* C */
   180   200   180   180 /* G */
   180   200   180   180 /* T */
/* CG C C ..AT */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG C G ..AT */
   280   300   280   280 /* A */
   280   300   280   280 /* C */
   280   300   280   280 /* G */
   280   300   280   280 /* T */
/* CG C T ..AT */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG G A ..AT */
   180   200   180   180 /* A */
   180   200   180   180 /* C */
   180   200   180   180 /* G */
   180   200   180   180 /* T */
/* CG G C ..AT */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG G G ..AT */
   280   300   280   280 /* A */
   280   300   280   280 /* C */
   280   300   280   280 /* G */
   280   300   280   280 /* T */
/* CG G T ..AT */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG T A ..AT */
   180   200   180   180 /* A */
   180   200   180   180 /* C */
   180   200   180   180 /* G */
   180   200   180   180 /* T */
/* CG T C ..AT */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG T G ..AT */
   280   300   280   280 /* A */
   280   300   280   280 /* C */
   280   300   280   280 /* G */
   280   300   280   280 /* T */
/* CG T T ..AT */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG A A ..TA */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG A C ..TA */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG A G ..TA */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG A T ..TA */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG C A ..TA */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG C C ..TA */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG C G ..TA */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG C T ..TA */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG G A ..TA */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG G C ..TA */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG G G ..TA */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG G T ..TA */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG T A ..TA */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG T C ..TA */
   170   190   170   170 /* A */
   170   190   170   170 /* C */
   170   190   170   170 /* G */
   170   190   170   170 /* T */
/* CG T G ..TA */
   150   170   150   150 /* A */
   150   170   150   150 /* C */
   150   170   150   150 /* G */
   150   170   150   150 /* T */
/* CG T T ..TA */
   130   150   130   130 /* A */
   130   150   130   130 /* C */
   130   150   130   130 /* G */
   130   150   130   130 /* T */
/* CG A A ..NN */
   180   200   180   180 /* A */
   180   200   180   180 /* C */
   180   200   180   180 /* G */
   180   200   180   180 /* T */
/* CG A C ..NN */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG A G ..NN */
   280   300   280   280 /* A */
   280   300   280   280 /* C */
   280   300   280   280 /* G */
   280   300   280   280 /* T */
/* CG A T ..NN */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG C A ..NN */
   180   200   180   180 /* A */
   180   200   180   180 /* C */
   180   200   180   180 /* G */
   180   200   180   180 /* T */
/* CG C C ..NN */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG C G ..NN */
   280   300   280   280 /* A */
   280   300   280   280 /* C */
   280   300   280   280 /* G */
   280   300   280   280 /* T */
/* CG C T ..NN */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG G A ..NN */
   180   200   180   180 /* A */
   180   200   180   180 /* C */
   180   200   180   180 /* G */
   180   200   180   180 /* T */
/* CG G C ..NN */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG G G ..NN */
   280   300   280   280 /* A */
   280   300   280   280 /* C */
   280   300   280   280 /* G */
   280   300   280   280 /* T */
/* CG G T ..NN */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG T A ..NN */
   180   200   180   180 /* A */
   180   200   180   180 /* C */
   180   200   180   180 /* G */
   180   200   180   180 /* T */
/* CG T C ..NN */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* CG T G ..NN */
   280   300   280   280 /* A */
   280   300   280   280 /* C */
   280   300   280   280 /* G */
   280   300   280   280 /* T */
/* CG T T ..NN */
   220   240   220   220 /* A */
   220   240   220   220 /* C */
   220   240   220   220 /* G */
   220   240   220   220 /* T */
/* GC A A ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC A C ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC A G ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC A T ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC C A ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC C C ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC C G ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC C T ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC G A ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC G C ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC G G ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC G T ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC T A ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC T C ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC T G ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC T T ..CG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC A A ..GC */
    50   100    90    70 /* A */
    50   100    90    70 /* C */
    50   100    90    70 /* G */
    50   100    90    70 /* T */
/* GC A C ..GC */
   100   150   140   120 /* A */
   100   150   140   120 /* C */
   100   150   140   120 /* G */
   100   150   140   120 /* T */
/* GC A G ..GC */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC A T ..GC */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GC C A ..GC */
    50   100    90    70 /* A */
    50   100    90    70 /* C */
    50   100    90    70 /* G */
    50   100    90    70 /* T */
/* GC C C ..GC */
   100   150   140   120 /* A */
   100   150   140   120 /* C */
   100   150   140   120 /* G */
   100   150   140   120 /* T */
/* GC C G ..GC */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC C T ..GC */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GC G A ..GC */
    50   100    90    70 /* A */
    50   100    90    70 /* C */
    50   100    90    70 /* G */
    50   100    90    70 /* T */
/* GC G C ..GC */
   100   150   140   120 /* A */
   100   150   140   120 /* C */
   100   150   140   120 /* G */
   100   150   140   120 /* T */
/* GC G G ..GC */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC G T ..GC */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GC T A ..GC */
    50   100    90    70 /* A */
    50   100    90    70 /* C */
    50   100    90    70 /* G */
    50   100    90    70 /* T */
/* GC T C ..GC */
   100   150   140   120 /* A */
   100   150   140   120 /* C */
   100   150   140   120 /* G */
   100   150   140   120 /* T */
/* GC T G ..GC */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC T T ..GC */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GC A A ..GT */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GC A C ..GT */
   120   170   160   140 /* A */
   120   170   160   140 /* C */
   120   170   160   140 /* G */
   120   170   160   140 /* T */
/* GC A G ..GT */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC A T ..GT */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC C A ..GT */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GC C C ..GT */
   120   170   160   140 /* A */
   120   170   160   140 /* C */
   120   170   160   140 /* G */
   120   170   160   140 /* T */
/* GC C G ..GT */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC C T ..GT */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC G A ..GT */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GC G C ..GT */
   120   170   160   140 /* A */
   120   170   160   140 /* C */
   120   170   160   140 /* G */
   120   170   160   140 /* T */
/* GC G G ..GT */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC G T ..GT */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC T A ..GT */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GC T C ..GT */
   120   170   160   140 /* A */
   120   170   160   140 /* C */
   120   170   160   140 /* G */
   120   170   160   140 /* T */
/* GC T G ..GT */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC T T ..GT */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC A A ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC A C ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC A G ..TG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC A T ..TG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC C A ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC C C ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC C G ..TG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC C T ..TG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC G A ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC G C ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC G G ..TG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC G T ..TG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC T A ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC T C ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC T G ..TG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC T T ..TG */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC A A ..AT */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GC A C ..AT */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC A G ..AT */
   240   290   280   260 /* A */
   240   290   280   260 /* C */
   240   290   280   260 /* G */
   240   290   280   260 /* T */
/* GC A T ..AT */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC C A ..AT */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GC C C ..AT */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC C G ..AT */
   240   290   280   260 /* A */
   240   290   280   260 /* C */
   240   290   280   260 /* G */
   240   290   280   260 /* T */
/* GC C T ..AT */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC G A ..AT */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GC G C ..AT */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC G G ..AT */
   240   290   280   260 /* A */
   240   290   280   260 /* C */
   240   290   280   260 /* G */
   240   290   280   260 /* T */
/* GC G T ..AT */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC T A ..AT */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GC T C ..AT */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC T G ..AT */
   240   290   280   260 /* A */
   240   290   280   260 /* C */
   240   290   280   260 /* G */
   240   290   280   260 /* T */
/* GC T T ..AT */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC A A ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC A C ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC A G ..TA */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC A T ..TA */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC C A ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC C C ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC C G ..TA */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC C T ..TA */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC G A ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC G C ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC G G ..TA */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC G T ..TA */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC T A ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC T C ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GC T G ..TA */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GC T T ..TA */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GC A A ..NN */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GC A C ..NN */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC A G ..NN */
   240   290   280   260 /* A */
   240   290   280   260 /* C */
   240   290   280   260 /* G */
   240   290   280   260 /* T */
/* GC A T ..NN */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC C A ..NN */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GC C C ..NN */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC C G ..NN */
   240   290   280   260 /* A */
   240   290   280   260 /* C */
   240   290   280   260 /* G */
   240   290   280   260 /* T */
/* GC C T ..NN */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC G A ..NN */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GC G C ..NN */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC G G ..NN */
   240   290   280   260 /* A */
   240   290   280   260 /* C */
   240   290   280   260 /* G */
   240   290   280   260 /* T */
/* GC G T ..NN */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC T A ..NN */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GC T C ..NN */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GC T G ..NN */
   240   290   280   260 /* A */
   240   290   280   260 /* C */
   240   290   280   260 /* G */
   240   290   280   260 /* T */
/* GC T T ..NN */
   180   230   220   200 /* A */
   180   230   220   200 /* C */
   180   230   220   200 /* G */
   180   230   220   200 /* T */
/* GT A A ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT A C ..CG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT A G ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT A T ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT C A ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT C C ..CG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT C G ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT C T ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT G A ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT G C ..CG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT G G ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT G T ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT T A ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT T C ..CG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT T G ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT T T ..CG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT A A ..GC */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GT A C ..GC */
   120   170   160   140 /* A */
   120   170   160   140 /* C */
   120   170   160   140 /* G */
   120   170   160   140 /* T */
/* GT A G ..GC */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT A T ..GC */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GT C A ..GC */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GT C C ..GC */
   120   170   160   140 /* A */
   120   170   160   140 /* C */
   120   170   160   140 /* G */
   120   170   160   140 /* T */
/* GT C G ..GC */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT C T ..GC */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GT G A ..GC */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GT G C ..GC */
   120   170   160   140 /* A */
   120   170   160   140 /* C */
   120   170   160   140 /* G */
   120   170   160   140 /* T */
/* GT G G ..GC */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT G T ..GC */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GT T A ..GC */
    70   120   110    90 /* A */
    70   120   110    90 /* C */
    70   120   110    90 /* G */
    70   120   110    90 /* T */
/* GT T C ..GC */
   120   170   160   140 /* A */
   120   170   160   140 /* C */
   120   170   160   140 /* G */
   120   170   160   140 /* T */
/* GT T G ..GC */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT T T ..GC */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GT A A ..GT */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GT A C ..GT */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GT A G ..GT */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT A T ..GT */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT C A ..GT */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GT C C ..GT */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GT C G ..GT */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT C T ..GT */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT G A ..GT */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GT G C ..GT */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GT G G ..GT */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT G T ..GT */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT T A ..GT */
    90   140   130   110 /* A */
    90   140   130   110 /* C */
    90   140   130   110 /* G */
    90   140   130   110 /* T */
/* GT T C ..GT */
   140   190   180   160 /* A */
   140   190   180   160 /* C */
   140   190   180   160 /* G */
   140   190   180   160 /* T */
/* GT T G ..GT */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT T T ..GT */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT A A ..TG */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT A C ..TG */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT A G ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT A T ..TG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT C A ..TG */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT C C ..TG */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT C G ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT C T ..TG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT G A ..TG */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT G C ..TG */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT G G ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT G T ..TG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT T A ..TG */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT T C ..TG */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT T G ..TG */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT T T ..TG */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT A A ..AT */
   160   210   200   180 /* A */
   160   210   200   180 /* C */
   160   210   200   180 /* G */
   160   210   200   180 /* T */
/* GT A C ..AT */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT A G ..AT */
   260   310   300   280 /* A */
   260   310   300   280 /* C */
   260   310   300   280 /* G */
   260   310   300   280 /* T */
/* GT A T ..AT */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT C A ..AT */
   160   210   200   180 /* A */
   160   210   200   180 /* C */
   160   210   200   180 /* G */
   160   210   200   180 /* T */
/* GT C C ..AT */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT C G ..AT */
   260   310   300   280 /* A */
   260   310   300   280 /* C */
   260   310   300   280 /* G */
   260   310   300   280 /* T */
/* GT C T ..AT */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT G A ..AT */
   160   210   200   180 /* A */
   160   210   200   180 /* C */
   160   210   200   180 /* G */
   160   210   200   180 /* T */
/* GT G C ..AT */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT G G ..AT */
   260   310   300   280 /* A */
   260   310   300   280 /* C */
   260   310   300   280 /* G */
   260   310   300   280 /* T */
/* GT G T ..AT */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT T A ..AT */
   160   210   200   180 /* A */
   160   210   200   180 /* C */
   160   210   200   180 /* G */
   160   210   200   180 /* T */
/* GT T C ..AT */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT T G ..AT */
   260   310   300   280 /* A */
   260   310   300   280 /* C */
   260   310   300   280 /* G */
   260   310   300   280 /* T */
/* GT T T ..AT */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT A A ..TA */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT A C ..TA */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT A G ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT A T ..TA */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT C A ..TA */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT C C ..TA */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT C G ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT C T ..TA */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT G A ..TA */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT G C ..TA */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT G G ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT G T ..TA */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT T A ..TA */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT T C ..TA */
   150   200   190   170 /* A */
   150   200   190   170 /* C */
   150   200   190   170 /* G */
   150   200   190   170 /* T */
/* GT T G ..TA */
   130   180   170   150 /* A */
   130   180   170   150 /* C */
   130   180   170   150 /* G */
   130   180   170   150 /* T */
/* GT T T ..TA */
   110   160   150   130 /* A */
   110   160   150   130 /* C */
   110   160   150   130 /* G */
   110   160   150   130 /* T */
/* GT A A ..NN */
   160   210   200   180 /* A */
   160   210   200   180 /* C */
   160   210   200   180 /* G */
   160   210   200   180 /* T */
/* GT A C ..NN */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT A G ..NN */
   260   310   300   280 /* A */
   260   310   300   280 /* C */
   260   310   300   280 /* G */
   260   310   300   280 /* T */
/* GT A T ..NN */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT C A ..NN */
   160   210   200   180 /* A */
   160   210   200   180 /* C */
   160   210   200   180 /* G */
   160   210   200   180 /* T */
/* GT C C ..NN */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT C G ..NN */
   260   310   300   280 /* A */
   260   310   300   280 /* C */
   260   310   300   280 /* G */
   260   310   300   280 /* T */
/* GT C T ..NN */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT G A ..NN */
   160   210   200   180 /* A */
   160   210   200   180 /* C */
   160   210   200   180 /* G */
   160   210   200   180 /* T */
/* GT G C ..NN */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT G G ..NN */
   260   310   300   280 /* A */
   260   310   300   280 /* C */
   260   310   300   280 /* G */
   260   310   300   280 /* T */
/* GT G T ..NN */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT T A ..NN */
   160   210   200   180 /* A */
   160   210   200   180 /* C */
   160   210   200   180 /* G */
   160   210   200   180 /* T */
/* GT T C ..NN */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* GT T G ..NN */
   260   310   300   280 /* A */
   260   310   300   280 /* C */
   260   310   300   280 /* G */
   260   310   300   280 /* T */
/* GT T T ..NN */
   200   250   240   220 /* A */
   200   250   240   220 /* C */
   200   250   240   220 /* G */
   200   250   240   220 /* T */
/* TG A A ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG A C ..CG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG A G ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG A T ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG C A ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG C C ..CG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG C G ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG C T ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG G A ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG G C ..CG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG G G ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG G T ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG T A ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG T C ..CG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG T G ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG T T ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG A A ..GC */
   130   130   110    90 /* A */
   130   130   110    90 /* C */
   130   130   110    90 /* G */
   130   130   110    90 /* T */
/* TG A C ..GC */
   180   180   160   140 /* A */
   180   180   160   140 /* C */
   180   180   160   140 /* G */
   180   180   160   140 /* T */
/* TG A G ..GC */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG A T ..GC */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TG C A ..GC */
   130   130   110    90 /* A */
   130   130   110    90 /* C */
   130   130   110    90 /* G */
   130   130   110    90 /* T */
/* TG C C ..GC */
   180   180   160   140 /* A */
   180   180   160   140 /* C */
   180   180   160   140 /* G */
   180   180   160   140 /* T */
/* TG C G ..GC */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG C T ..GC */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TG G A ..GC */
   130   130   110    90 /* A */
   130   130   110    90 /* C */
   130   130   110    90 /* G */
   130   130   110    90 /* T */
/* TG G C ..GC */
   180   180   160   140 /* A */
   180   180   160   140 /* C */
   180   180   160   140 /* G */
   180   180   160   140 /* T */
/* TG G G ..GC */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG G T ..GC */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TG T A ..GC */
   130   130   110    90 /* A */
   130   130   110    90 /* C */
   130   130   110    90 /* G */
   130   130   110    90 /* T */
/* TG T C ..GC */
   180   180   160   140 /* A */
   180   180   160   140 /* C */
   180   180   160   140 /* G */
   180   180   160   140 /* T */
/* TG T G ..GC */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG T T ..GC */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TG A A ..GT */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TG A C ..GT */
   200   200   180   160 /* A */
   200   200   180   160 /* C */
   200   200   180   160 /* G */
   200   200   180   160 /* T */
/* TG A G ..GT */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG A T ..GT */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG C A ..GT */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TG C C ..GT */
   200   200   180   160 /* A */
   200   200   180   160 /* C */
   200   200   180   160 /* G */
   200   200   180   160 /* T */
/* TG C G ..GT */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG C T ..GT */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG G A ..GT */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TG G C ..GT */
   200   200   180   160 /* A */
   200   200   180   160 /* C */
   200   200   180   160 /* G */
   200   200   180   160 /* T */
/* TG G G ..GT */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG G T ..GT */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG T A ..GT */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TG T C ..GT */
   200   200   180   160 /* A */
   200   200   180   160 /* C */
   200   200   180   160 /* G */
   200   200   180   160 /* T */
/* TG T G ..GT */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG T T ..GT */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG A A ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG A C ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG A G ..TG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG A T ..TG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG C A ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG C C ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG C G ..TG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG C T ..TG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG G A ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG G C ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG G G ..TG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG G T ..TG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG T A ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG T C ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG T G ..TG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG T T ..TG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG A A ..AT */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TG A C ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG A G ..AT */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TG A T ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG C A ..AT */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TG C C ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG C G ..AT */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TG C T ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG G A ..AT */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TG G C ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG G G ..AT */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TG G T ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG T A ..AT */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TG T C ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG T G ..AT */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TG T T ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG A A ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG A C ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG A G ..TA */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG A T ..TA */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG C A ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG C C ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG C G ..TA */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG C T ..TA */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG G A ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG G C ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG G G ..TA */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG G T ..TA */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG T A ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG T C ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TG T G ..TA */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TG T T ..TA */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TG A A ..NN */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TG A C ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG A G ..NN */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TG A T ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG C A ..NN */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TG C C ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG C G ..NN */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TG C T ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG G A ..NN */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TG G C ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG G G ..NN */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TG G T ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG T A ..NN */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TG T C ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TG T G ..NN */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TG T T ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* AT A A ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT A C ..CG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT A G ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT A T ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT C A ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT C C ..CG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT C G ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT C T ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT G A ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT G C ..CG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT G G ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT G T ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT T A ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT T C ..CG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT T G ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT T T ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT A A ..GC */
   140   180   240   180 /* A */
   140   180   240   180 /* C */
   140   180   240   180 /* G */
   140   180   240   180 /* T */
/* AT A C ..GC */
   190   230   290   230 /* A */
   190   230   290   230 /* C */
   190   230   290   230 /* G */
   190   230   290   230 /* T */
/* AT A G ..GC */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT A T ..GC */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* AT C A ..GC */
   140   180   240   180 /* A */
   140   180   240   180 /* C */
   140   180   240   180 /* G */
   140   180   240   180 /* T */
/* AT C C ..GC */
   190   230   290   230 /* A */
   190   230   290   230 /* C */
   190   230   290   230 /* G */
   190   230   290   230 /* T */
/* AT C G ..GC */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT C T ..GC */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* AT G A ..GC */
   140   180   240   180 /* A */
   140   180   240   180 /* C */
   140   180   240   180 /* G */
   140   180   240   180 /* T */
/* AT G C ..GC */
   190   230   290   230 /* A */
   190   230   290   230 /* C */
   190   230   290   230 /* G */
   190   230   290   230 /* T */
/* AT G G ..GC */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT G T ..GC */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* AT T A ..GC */
   140   180   240   180 /* A */
   140   180   240   180 /* C */
   140   180   240   180 /* G */
   140   180   240   180 /* T */
/* AT T C ..GC */
   190   230   290   230 /* A */
   190   230   290   230 /* C */
   190   230   290   230 /* G */
   190   230   290   230 /* T */
/* AT T G ..GC */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT T T ..GC */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* AT A A ..GT */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* AT A C ..GT */
   210   250   310   250 /* A */
   210   250   310   250 /* C */
   210   250   310   250 /* G */
   210   250   310   250 /* T */
/* AT A G ..GT */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT A T ..GT */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT C A ..GT */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* AT C C ..GT */
   210   250   310   250 /* A */
   210   250   310   250 /* C */
   210   250   310   250 /* G */
   210   250   310   250 /* T */
/* AT C G ..GT */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT C T ..GT */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT G A ..GT */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* AT G C ..GT */
   210   250   310   250 /* A */
   210   250   310   250 /* C */
   210   250   310   250 /* G */
   210   250   310   250 /* T */
/* AT G G ..GT */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT G T ..GT */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT T A ..GT */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* AT T C ..GT */
   210   250   310   250 /* A */
   210   250   310   250 /* C */
   210   250   310   250 /* G */
   210   250   310   250 /* T */
/* AT T G ..GT */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT T T ..GT */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT A A ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT A C ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT A G ..TG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT A T ..TG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT C A ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT C C ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT C G ..TG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT C T ..TG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT G A ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT G C ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT G G ..TG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT G T ..TG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT T A ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT T C ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT T G ..TG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT T T ..TG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT A A ..AT */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* AT A C ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT A G ..AT */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* AT A T ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT C A ..AT */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* AT C C ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT C G ..AT */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* AT C T ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT G A ..AT */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* AT G C ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT G G ..AT */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* AT G T ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT T A ..AT */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* AT T C ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT T G ..AT */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* AT T T ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT A A ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT A C ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT A G ..TA */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT A T ..TA */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT C A ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT C C ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT C G ..TA */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT C T ..TA */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT G A ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT G C ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT G G ..TA */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT G T ..TA */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT T A ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT T C ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* AT T G ..TA */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* AT T T ..TA */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* AT A A ..NN */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* AT A C ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT A G ..NN */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* AT A T ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT C A ..NN */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* AT C C ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT C G ..NN */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* AT C T ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT G A ..NN */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* AT G C ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT G G ..NN */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* AT G T ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT T A ..NN */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* AT T C ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* AT T G ..NN */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* AT T T ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* TA A A ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA A C ..CG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA A G ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA A T ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA C A ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA C C ..CG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA C G ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA C T ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA G A ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA G C ..CG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA G G ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA G T ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA T A ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA T C ..CG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA T G ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA T T ..CG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA A A ..GC */
   130   130   110    90 /* A */
   130   130   110    90 /* C */
   130   130   110    90 /* G */
   130   130   110    90 /* T */
/* TA A C ..GC */
   180   180   160   140 /* A */
   180   180   160   140 /* C */
   180   180   160   140 /* G */
   180   180   160   140 /* T */
/* TA A G ..GC */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA A T ..GC */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TA C A ..GC */
   130   130   110    90 /* A */
   130   130   110    90 /* C */
   130   130   110    90 /* G */
   130   130   110    90 /* T */
/* TA C C ..GC */
   180   180   160   140 /* A */
   180   180   160   140 /* C */
   180   180   160   140 /* G */
   180   180   160   140 /* T */
/* TA C G ..GC */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA C T ..GC */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TA G A ..GC */
   130   130   110    90 /* A */
   130   130   110    90 /* C */
   130   130   110    90 /* G */
   130   130   110    90 /* T */
/* TA G C ..GC */
   180   180   160   140 /* A */
   180   180   160   140 /* C */
   180   180   160   140 /* G */
   180   180   160   140 /* T */
/* TA G G ..GC */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA G T ..GC */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TA T A ..GC */
   130   130   110    90 /* A */
   130   130   110    90 /* C */
   130   130   110    90 /* G */
   130   130   110    90 /* T */
/* TA T C ..GC */
   180   180   160   140 /* A */
   180   180   160   140 /* C */
   180   180   160   140 /* G */
   180   180   160   140 /* T */
/* TA T G ..GC */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA T T ..GC */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TA A A ..GT */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TA A C ..GT */
   200   200   180   160 /* A */
   200   200   180   160 /* C */
   200   200   180   160 /* G */
   200   200   180   160 /* T */
/* TA A G ..GT */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA A T ..GT */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA C A ..GT */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TA C C ..GT */
   200   200   180   160 /* A */
   200   200   180   160 /* C */
   200   200   180   160 /* G */
   200   200   180   160 /* T */
/* TA C G ..GT */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA C T ..GT */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA G A ..GT */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TA G C ..GT */
   200   200   180   160 /* A */
   200   200   180   160 /* C */
   200   200   180   160 /* G */
   200   200   180   160 /* T */
/* TA G G ..GT */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA G T ..GT */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA T A ..GT */
   150   150   130   110 /* A */
   150   150   130   110 /* C */
   150   150   130   110 /* G */
   150   150   130   110 /* T */
/* TA T C ..GT */
   200   200   180   160 /* A */
   200   200   180   160 /* C */
   200   200   180   160 /* G */
   200   200   180   160 /* T */
/* TA T G ..GT */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA T T ..GT */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA A A ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA A C ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA A G ..TG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA A T ..TG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA C A ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA C C ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA C G ..TG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA C T ..TG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA G A ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA G C ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA G G ..TG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA G T ..TG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA T A ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA T C ..TG */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA T G ..TG */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA T T ..TG */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA A A ..AT */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TA A C ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA A G ..AT */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TA A T ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA C A ..AT */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TA C C ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA C G ..AT */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TA C T ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA G A ..AT */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TA G C ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA G G ..AT */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TA G T ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA T A ..AT */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TA T C ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA T G ..AT */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TA T T ..AT */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA A A ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA A C ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA A G ..TA */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA A T ..TA */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA C A ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA C C ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA C G ..TA */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA C T ..TA */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA G A ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA G C ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA G G ..TA */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA G T ..TA */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA T A ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA T C ..TA */
   210   210   190   170 /* A */
   210   210   190   170 /* C */
   210   210   190   170 /* G */
   210   210   190   170 /* T */
/* TA T G ..TA */
   190   190   170   150 /* A */
   190   190   170   150 /* C */
   190   190   170   150 /* G */
   190   190   170   150 /* T */
/* TA T T ..TA */
   170   170   150   130 /* A */
   170   170   150   130 /* C */
   170   170   150   130 /* G */
   170   170   150   130 /* T */
/* TA A A ..NN */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TA A C ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA A G ..NN */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TA A T ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA C A ..NN */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TA C C ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA C G ..NN */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TA C T ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA G A ..NN */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TA G C ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA G G ..NN */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TA G T ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA T A ..NN */
   220   220   200   180 /* A */
   220   220   200   180 /* C */
   220   220   200   180 /* G */
   220   220   200   180 /* T */
/* TA T C ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* TA T G ..NN */
   320   320   300   280 /* A */
   320   320   300   280 /* C */
   320   320   300   280 /* G */
   320   320   300   280 /* T */
/* TA T T ..NN */
   260   260   240   220 /* A */
   260   260   240   220 /* C */
   260   260   240   220 /* G */
   260   260   240   220 /* T */
/* NN A A ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN A C ..CG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN A G ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN A T ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN C A ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN C C ..CG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN C G ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN C T ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN G A ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN G C ..CG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN G G ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN G T ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN T A ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN T C ..CG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN T G ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN T T ..CG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN A A ..GC */
   140   180   240   180 /* A */
   140   180   240   180 /* C */
   140   180   240   180 /* G */
   140   180   240   180 /* T */
/* NN A C ..GC */
   190   230   290   230 /* A */
   190   230   290   230 /* C */
   190   230   290   230 /* G */
   190   230   290   230 /* T */
/* NN A G ..GC */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN A T ..GC */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* NN C A ..GC */
   140   180   240   180 /* A */
   140   180   240   180 /* C */
   140   180   240   180 /* G */
   140   180   240   180 /* T */
/* NN C C ..GC */
   190   230   290   230 /* A */
   190   230   290   230 /* C */
   190   230   290   230 /* G */
   190   230   290   230 /* T */
/* NN C G ..GC */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN C T ..GC */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* NN G A ..GC */
   140   180   240   180 /* A */
   140   180   240   180 /* C */
   140   180   240   180 /* G */
   140   180   240   180 /* T */
/* NN G C ..GC */
   190   230   290   230 /* A */
   190   230   290   230 /* C */
   190   230   290   230 /* G */
   190   230   290   230 /* T */
/* NN G G ..GC */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN G T ..GC */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* NN T A ..GC */
   140   180   240   180 /* A */
   140   180   240   180 /* C */
   140   180   240   180 /* G */
   140   180   240   180 /* T */
/* NN T C ..GC */
   190   230   290   230 /* A */
   190   230   290   230 /* C */
   190   230   290   230 /* G */
   190   230   290   230 /* T */
/* NN T G ..GC */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN T T ..GC */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* NN A A ..GT */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* NN A C ..GT */
   210   250   310   250 /* A */
   210   250   310   250 /* C */
   210   250   310   250 /* G */
   210   250   310   250 /* T */
/* NN A G ..GT */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN A T ..GT */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN C A ..GT */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* NN C C ..GT */
   210   250   310   250 /* A */
   210   250   310   250 /* C */
   210   250   310   250 /* G */
   210   250   310   250 /* T */
/* NN C G ..GT */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN C T ..GT */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN G A ..GT */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* NN G C ..GT */
   210   250   310   250 /* A */
   210   250   310   250 /* C */
   210   250   310   250 /* G */
   210   250   310   250 /* T */
/* NN G G ..GT */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN G T ..GT */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN T A ..GT */
   160   200   260   200 /* A */
   160   200   260   200 /* C */
   160   200   260   200 /* G */
   160   200   260   200 /* T */
/* NN T C ..GT */
   210   250   310   250 /* A */
   210   250   310   250 /* C */
   210   250   310   250 /* G */
   210   250   310   250 /* T */
/* NN T G ..GT */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN T T ..GT */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN A A ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN A C ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN A G ..TG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN A T ..TG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN C A ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN C C ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN C G ..TG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN C T ..TG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN G A ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN G C ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN G G ..TG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN G T ..TG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN T A ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN T C ..TG */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN T G ..TG */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN T T ..TG */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN A A ..AT */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* NN A C ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN A G ..AT */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* NN A T ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN C A ..AT */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* NN C C ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN C G ..AT */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* NN C T ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN G A ..AT */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* NN G C ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN G G ..AT */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* NN G T ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN T A ..AT */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* NN T C ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN T G ..AT */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* NN T T ..AT */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN A A ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN A C ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN A G ..TA */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN A T ..TA */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN C A ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN C C ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN C G ..TA */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN C T ..TA */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN G A ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN G C ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN G G ..TA */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN G T ..TA */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN T A ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN T C ..TA */
   220   260   320   260 /* A */
   220   260   320   260 /* C */
   220   260   320   260 /* G */
   220   260   320   260 /* T */
/* NN T G ..TA */
   200   240   300   240 /* A */
   200   240   300   240 /* C */
   200   240   300   240 /* G */
   200   240   300   240 /* T */
/* NN T T ..TA */
   180   220   280   220 /* A */
   180   220   280   220 /* C */
   180   220   280   220 /* G */
   180   220   280   220 /* T */
/* NN A A ..NN */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* NN A C ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN A G ..NN */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* NN A T ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN C A ..NN */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* NN C C ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN C G ..NN */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* NN C T ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN G A ..NN */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* NN G C ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN G G ..NN */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* NN G T ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN T A ..NN */
   230   270   330   270 /* A */
   230   270   330   270 /* C */
   230   270   330   270 /* G */
   230   270   330   270 /* T */
/* NN T C ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */
/* NN T G ..NN */
   330   370   430   370 /* A */
   330   370   430   370 /* C */
   330   370   430   370 /* G */
   330   370   430   370 /* T */
/* NN T T ..NN */
   270   310   370   310 /* A */
   270   310   370   310 /* C */
   270   310   370   310 /* G */
   270   310   370   310 /* T */

# int22_enthalpies
/* CG A A ..CG */
  -920  -880  -930  -960 /* A */
  -920  -880  -930  -960 /* C */
  -920  -880  -930  -960 /* G */
  -920  -880  -930  -960 /* T */
/* CG A C ..CG */
  -880  -840  -890  -920 /* A */
  -880  -840  -890  -920 /* C */
  -880  -840  -890  -920 /* G */
  -880  -840  -890  -920 /* T */
/* CG A G ..CG */
  -930  -890  -940  -970 /* A */
  -930  -890  -940  -970 /* C */
  -930  -890  -940  -970 /* G */
  -930  -890  -940  -970 /* T */
/* CG A T ..CG */
  -960  -920  -970 -1000 /* A */
  -960  -920  -970 -1000 /* C */
  -960  -920  -970 -1000 /* G */
  -960  -920  -970 -1000 /* T */
/* CG C A ..CG */
  -920  -880  -930  -960 /* A */
  -920  -880  -930  -960 /* C */
  -920  -880  -930  -960 /* G */
  -920  -880  -930  -960 /* T */
/* CG C C ..CG */
  -880  -840  -890  -920 /* A */
  -880  -840  -890  -920 /* C */
  -880  -840  -890  -920 /* G */
  -880  -840  -890  -920 /* T */
/* CG C G ..CG */
  -930  -890  -940  -970 /* A */
  -930  -890  -940  -970 /* C */
  -930  -890  -940  -970 /* G */
  -930  -890  -940  -970 /* T */
/* CG C T ..CG */
  -960  -920  -970 -1000 /* A */
  -960  -920  -970 -1000 /* C */
  -960  -920  -970 -1000 /* G */
  -960  -920  -970 -1000 /* T */
/* CG G A ..CG */
  -920  -880  -930  -960 /* A */
  -920  -880  -930  -960 /* C */
  -920  -880  -930  -960 /* G */
  -920  -880  -930  -960 /* T */
/* CG G C ..CG */
  -880  -840  -890  -920 /* A */
  -880  -840  -890  -920 /* C */
  -880  -840  -890  -920 /* G */
  -880  -840  -890  -920 /* T */
/* CG G G ..CG */
  -930  -890  -940  -970 /* A */
  -930  -890  -940  -970 /* C */
  -930  -890  -940  -970 /* G */
  -930  -890  -940  -970 /* T */
/* CG G T ..CG */
  -960  -920  -970 -1000 /* A */
  -960  -920  -970 -1000 /* C */
  -960  -920  -970 -1000 /* G */
  -960  -920  -970 -1000 /* T */
/* CG T A ..CG */
  -920  -880  -930  -960 /* A */
  -920  -880  -930  -960 /* C */
  -920  -880  -930  -960 /* G */
  -920  -880  -930  -960 /* T */
/* CG T C ..CG */
  -880  -840  -890  -920 /* A */
  -880  -840  -890  -920 /* C */
  -880  -840  -890  -920 /* G */
  -880  -840  -890  -920 /* T */
/* CG T G ..CG */
  -930  -890  -940  -970 /* A */
  -930  -890  -940  -970 /* C */
  -930  -890  -940  -970 /* G */
  -930  -890  -940  -970 /* T */
/* CG T T ..CG */
  -960  -920  -970 -1000 /* A */
  -960  -920  -970 -1000 /* C */
  -960  -920  -970 -1000 /* G */
  -960  -920  -970 -1000 /* T */
/* CG A A ..GC */
 -1080 -1040 -1090 -1120 /* A */
 -1080 -1040 -1090 -1120 /* C */
 -1080 -1040 -1090 -1120 /* G */
 -1080 -1040 -1090 -1120 /* T */
/* CG A C ..GC */
  -790  -750  -800  -830 /* A */
  -790  -750  -800  -830 /* C */
  -790  -750  -800  -830 /* G */
  -790  -750  -800  -830 /* T */
/* CG A G ..GC */
  -970  -930  -980 -1010 /* A */
  -970  -930  -980 -1010 /* C */
  -970  -930  -980 -1010 /* G */
  -970  -930  -980 -1010 /* T */
/* CG A T ..GC */
  -550  -510  -560  -590 /* A */
  -550  -510  -560  -590 /* C */
  -550  -510  -560  -590 /* G */
  -550  -510  -560  -590 /* T */
/* CG C A ..GC */
 -1080 -1040 -1090 -1120 /* A */
 -1080 -1040 -1090 -1120 /* C */
 -1080 -1040 -1090 -1120 /* G */
 -1080 -1040 -1090 -1120 /* T */
/* CG C C ..GC */
  -790  -750  -800  -830 /* A */
  -790  -750  -800  -830 /* C */
  -790  -750  -800  -830 /* G */
  -790  -750  -800  -830 /* T */
/* CG C G ..GC */
  -970  -930  -980 -1010 /* A */
  -970  -930  -980 -1010 /* C */
  -970  -930  -980 -1010 /* G */
  -970  -930  -980 -1010 /* T */
/* CG C T ..GC */
  -550  -510  -560  -590 /* A */
  -550  -510  -560  -590 /* C */
  -550  -510  -560  -590 /* G */
  -550  -510  -560  -590 /* T */
/* CG G A ..GC */
 -1080 -1040 -1090 -1120 /* A */
 -1080 -1040 -1090 -1120 /* C */
 -1080 -1040 -1090 -1120 /* G */
 -1080 -1040 -1090 -1120 /* T */
/* CG G C ..GC */
  -790  -750  -800  -830 /* A */
  -790  -750  -800  -830 /* C */
  -790  -750  -800  -830 /* G */
  -790  -750  -800  -830 /* T */
/* CG G G ..GC */
  -970  -930  -980 -1010 /* A */
  -970  -930  -980 -1010 /* C */
  -970  -930  -980 -1010 /* G */
  -970  -930  -980 -1010 /* T */
/* CG G T ..GC */
  -550  -510  -560  -590 /* A */
  -550  -510  -560  -590 /* C */
  -550  -510  -560  -590 /* G */
  -550  -510  -560  -590 /* T */
/* CG T A ..GC */
 -1080 -1040 -1090 -1120 /* A */
 -1080 -1040 -1090 -1120 /* C */
 -1080 -1040 -1090 -1120 /* G */
 -1080 -1040 -1090 -1120 /* T */
/* CG T C ..GC */
  -790  -750  -800  -830 /* A */
  -790  -750  -800  -830 /* C */
  -790  -750  -800  -830 /* G */
  -790  -750  -800  -830 /* T */
/* CG T G ..GC */
  -970  -930  -980 -1010 /* A */
  -970  -930  -980 -1010 /* C */
  -970  -930  -980 -1010 /* G */
  -970  -930  -980 -1010 /* T */
/* CG T T ..GC */
  -550  -510  -560  -590 /* A */
  -550  -510  -560  -590 /* C */
  -550  -510  -560  -590 /* G */
  -550  -510  -560  -590 /* T */
/* CG A A ..GT */
  -760  -720  -770  -800 /* A */
  -760  -720  -770  -800 /* C */
  -760  -720  -770  -800 /* G */
  -760  -720  -770  -800 /* T */
/* CG A C ..GT */
  -470  -430  -480  -510 /* A */
  -470  -430  -480  -510 /* C */
  -470  -430  -480  -510 /* G */
  -470  -430  -480  -510 /* T */
/* CG A G ..GT */
  -650  -610  -660  -690 /* A */
  -650  -610  -660  -690 /* C */
  -650  -610  -660  -690 /* G */
  -650  -610  -660  -690 /* T */
/* CG A T ..GT */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG C A ..GT */
  -760  -720  -770  -800 /* A */
  -760  -720  -770  -800 /* C */
  -760  -720  -770  -800 /* G */
  -760  -720  -770  -800 /* T */
/* CG C C ..GT */
  -470  -430  -480  -510 /* A */
  -470  -430  -480  -510 /* C */
  -470  -430  -480  -510 /* G */
  -470  -430  -480  -510 /* T */
/* CG C G ..GT */
  -650  -610  -660  -690 /* A */
  -650  -610  -660  -690 /* C */
  -650  -610  -660  -690 /* G */
  -650  -610  -660  -690 /* T */
/* CG C T ..GT */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG G A ..GT */
  -760  -720  -770  -800 /* A */
  -760  -720  -770  -800 /* C */
  -760  -720  -770  -800 /* G */
  -760  -720  -770  -800 /* T */
/* CG G C ..GT */
  -470  -430  -480  -510 /* A */
  -470  -430  -480  -510 /* C */
  -470  -430  -480  -510 /* G */
  -470  -430  -480  -510 /* T */
/* CG G G ..GT */
  -650  -610  -660  -690 /* A */
  -650  -610  -660  -690 /* C */
  -650  -610  -660  -690 /* G */
  -650  -610  -660  -690 /* T */
/* CG G T ..GT */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG T A ..GT */
  -760  -720  -770  -800 /* A */
  -760  -720  -770  -800 /* C */
  -760  -720  -770  -800 /* G */
  -760  -720  -770  -800 /* T */
/* CG T C ..GT */
  -470  -430  -480  -510 /* A */
  -470  -430  -480  -510 /* C */
  -470  -430  -480  -510 /* G */
  -470  -430  -480  -510 /* T */
/* CG T G ..GT */
  -650  -610  -660  -690 /* A */
  -650  -610  -660  -690 /* C */
  -650  -610  -660  -690 /* G */
  -650  -610  -660  -690 /* T */
/* CG T T ..GT */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG A A ..TG */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG A C ..TG */
  -200  -160  -210  -240 /* A */
  -200  -160  -210  -240 /* C */
  -200  -160  -210  -240 /* G */
  -200  -160  -210  -240 /* T */
/* CG A G ..TG */
  -390  -350  -400  -430 /* A */
  -390  -350  -400  -430 /* C */
  -390  -350  -400  -430 /* G */
  -390  -350  -400  -430 /* T */
/* CG A T ..TG */
 -1040 -1000 -1050 -1080 /* A */
 -1040 -1000 -1050 -1080 /* C */
 -1040 -1000 -1050 -1080 /* G */
 -1040 -1000 -1050 -1080 /* T */
/* CG C A ..TG */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG C C ..TG */
  -200  -160  -210  -240 /* A */
  -200  -160  -210  -240 /* C */
  -200  -160  -210  -240 /* G */
  -200  -160  -210  -240 /* T */
/* CG C G ..TG */
  -390  -350  -400  -430 /* A */
  -390  -350  -400  -430 /* C */
  -390  -350  -400  -430 /* G */
  -390  -350  -400  -430 /* T */
/* CG C T ..TG */
 -1040 -1000 -1050 -1080 /* A */
 -1040 -1000 -1050 -1080 /* C */
 -1040 -1000 -1050 -1080 /* G */
 -1040 -1000 -1050 -1080 /* T */
/* CG G A ..TG */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG G C ..TG */
  -200  -160  -210  -240 /* A */
  -200  -160  -210  -240 /* C */
  -200  -160  -210  -240 /* G */
  -200  -160  -210  -240 /* T */
/* CG G G ..TG */
  -390  -350  -400  -430 /* A */
  -390  -350  -400  -430 /* C */
  -390  -350  -400  -430 /* G */
  -390  -350  -400  -430 /* T */
/* CG G T ..TG */
 -1040 -1000 -1050 -1080 /* A */
 -1040 -1000 -1050 -1080 /* C */
 -1040 -1000 -1050 -1080 /* G */
 -1040 -1000 -1050 -1080 /* T */
/* CG T A ..TG */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG T C ..TG */
  -200  -160  -210  -240 /* A */
  -200  -160  -210  -240 /* C */
  -200  -160  -210  -240 /* G */
  -200  -160  -210  -240 /* T */
/* CG T G ..TG */
  -390  -350  -400  -430 /* A */
  -390  -350  -400  -430 /* C */
  -390  -350  -400  -430 /* G */
  -390  -350  -400  -430 /* T */
/* CG T T ..TG */
 -1040 -1000 -1050 -1080 /* A */
 -1040 -1000 -1050 -1080 /* C */
 -1040 -1000 -1050 -1080 /* G */
 -1040 -1000 -1050 -1080 /* T */
/* CG A A ..AT */
   -60   -20   -70  -100 /* A */
   -60   -20   -70  -100 /* C */
   -60   -20   -70  -100 /* G */
   -60   -20   -70  -100 /* T */
/* CG A C ..AT */
    60   100    50    20 /* A */
    60   100    50    20 /* C */
    60   100    50    20 /* G */
    60   100    50    20 /* T */
/* CG A G ..AT */
    10    50     0   -30 /* A */
    10    50     0   -30 /* C */
    10    50     0   -30 /* G */
    10    50     0   -30 /* T */
/* CG A T ..AT */
   -30    10   -40   -70 /* A */
   -30    10   -40   -70 /* C */
   -30    10   -40   -70 /* G */
   -30    10   -40   -70 /* T */
/* CG C A ..AT */
   -60   -20   -70  -100 /* A */
   -60   -20   -70  -100 /* C */
   -60   -20   -70  -100 /* G */
   -60   -20   -70  -100 /* T */
/* CG C C ..AT */
    60   100    50    20 /* A */
    60   100    50    20 /* C */
    60   100    50    20 /* G */
    60   100    50    20 /* T */
/* CG C G ..AT */
    10    50     0   -30 /* A */
    10    50     0   -30 /* C */
    10    50     0   -30 /* G */
    10    50     0   -30 /* T */
/* CG C T ..AT */
   -30    10   -40   -70 /* A */
   -30    10   -40   -70 /* C */
   -30    10   -40   -70 /* G */
   -30    10   -40   -70 /* T */
/* CG G A ..AT */
   -60   -20   -70  -100 /* A */
   -60   -20   -70  -100 /* C */
   -60   -20   -70  -100 /* G */
   -60   -20   -70  -100 /* T */
/* CG G C ..AT */
    60   100    50    20 /* A */
    60   100    50    20 /* C */
    60   100    50    20 /* G */
    60   100    50    20 /* T */
/* CG G G ..AT */
    10    50     0   -30 /* A */
    10    50     0   -30 /* C */
    10    50     0   -30 /* G */
    10    50     0   -30 /* T */
/* CG G T ..AT */
   -30    10   -40   -70 /* A */
   -30    10   -40   -70 /* C */
   -30    10   -40   -70 /* G */
   -30    10   -40   -70 /* T */
/* CG T A ..AT */
   -60   -20   -70  -100 /* A */
   -60   -20   -70  -100 /* C */
   -60   -20   -70  -100 /* G */
   -60   -20   -70  -100 /* T */
/* CG T C ..AT */
    60   100    50    20 /* A */
    60   100    50    20 /* C */
    60   100    50    20 /* G */
    60   100    50    20 /* T */
/* CG T G ..AT */
    10    50     0   -30 /* A */
    10    50     0   -30 /* C */
    10    50     0   -30 /* G */
    10    50     0   -30 /* T */
/* CG T T ..AT */
   -30    10   -40   -70 /* A */
   -30    10   -40   -70 /* C */
   -30    10   -40   -70 /* G */
   -30    10   -40   -70 /* T */
/* CG A A ..TA */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG A C ..TA */
  -200  -160  -210  -240 /* A */
  -200  -160  -210  -240 /* C */
  -200  -160  -210  -240 /* G */
  -200  -160  -210  -240 /* T */
/* CG A G ..TA */
  -390  -350  -400  -430 /* A */
  -390  -350  -400  -430 /* C */
  -390  -350  -400  -430 /* G */
  -390  -350  -400  -430 /* T */
/* CG A T ..TA */
 -1040 -1000 -1050 -1080 /* A */
 -1040 -1000 -1050 -1080 /* C */
 -1040 -1000 -1050 -1080 /* G */
 -1040 -1000 -1050 -1080 /* T */
/* CG C A ..TA */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG C C ..TA */
  -200  -160  -210  -240 /* A */
  -200  -160  -210  -240 /* C */
  -200  -160  -210  -240 /* G */
  -200  -160  -210  -240 /* T */
/* CG C G ..TA */
  -390  -350  -400  -430 /* A */
  -390  -350  -400  -430 /* C */
  -390  -350  -400  -430 /* G */
  -390  -350  -400  -430 /* T */
/* CG C T ..TA */
 -1040 -1000 -1050 -1080 /* A */
 -1040 -1000 -1050 -1080 /* C */
 -1040 -1000 -1050 -1080 /* G */
 -1040 -1000 -1050 -1080 /* T */
/* CG G A ..TA */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG G C ..TA */
  -200  -160  -210  -240 /* A */
  -200  -160  -210  -240 /* C */
  -200  -160  -210  -240 /* G */
  -200  -160  -210  -240 /* T */
/* CG G G ..TA */
  -390  -350  -400  -430 /* A */
  -390  -350  -400  -430 /* C */
  -390  -350  -400  -430 /* G */
  -390  -350  -400  -430 /* T */
/* CG G T ..TA */
 -1040 -1000 -1050 -1080 /* A */
 -1040 -1000 -1050 -1080 /* C */
 -1040 -1000 -1050 -1080 /* G */
 -1040 -1000 -1050 -1080 /* T */
/* CG T A ..TA */
  -230  -190  -240  -270 /* A */
  -230  -190  -240  -270 /* C */
  -230  -190  -240  -270 /* G */
  -230  -190  -240  -270 /* T */
/* CG T C ..TA */
  -200  -160  -210  -240 /* A */
  -200  -160  -210  -240 /* C */
  -200  -160  -210  -240 /* G */
  -200  -160  -210  -240 /* T */
/* CG T G ..TA */
  -390  -350  -400  -430 /* A */
  -390  -350  -400  -430 /* C */
  -390  -350  -400  -430 /* G */
  -390  -350  -400  -430 /* T */
/* CG T T ..TA */
 -1040 -1000 -1050 -1080 /* A */
 -1040 -1000 -1050 -1080 /* C */
 -1040 -1000 -1050 -1080 /* G */
 -1040 -1000 -1050 -1080 /* T */
/* CG A A ..NN */
   -60   -20   -70  -100 /* A */
   -60   -20   -70  -100 /* C */
   -60   -20   -70  -100 /* G */
   -60   -20   -70  -100 /* T */
/* CG A C ..NN */
    60   100    50    20 /* A */
    60   100    50    20 /* C */
    60   100    50    20 /* G */
    60   100    50    20 /* T */
/* CG A G ..NN */
    10    50     0   -30 /* A */
    10    50     0   -30 /* C */
    10    50     0   -30 /* G */
    10    50     0   -30 /* T */
/* CG A T ..NN */
   -30    10   -40   -70 /* A */
   -30    10   -40   -70 /* C */
   -30    10   -40   -70 /* G */
   -30    10   -40   -70 /* T */
/* CG C A ..NN */
   -60   -20   -70  -100 /* A */
   -60   -20   -70  -100 /* C */
   -60   -20   -70  -100 /* G */
   -60   -20   -70  -100 /* T */
/* CG C C ..NN */
    60   100    50    20 /* A */
    60   100    50    20 /* C */
    60   100    50    20 /* G */
    60   100    50    20 /* T */
/* CG C G ..NN */
    10    50     0   -30 /* A */
    10    50     0   -30 /* C */
    10    50     0   -30 /* G */
    10    50     0   -30 /* T */
/* CG C T ..NN */
   -30    10   -40   -70 /* A */
   -30    10   -40   -70 /* C */
   -30    10   -40   -70 /* G */
   -30    10   -40   -70 /* T */
/* CG G A ..NN */
   -60   -20   -70  -100 /* A */
   -60   -20   -70  -100 /* C */
   -60   -20   -70  -100 /* G */
   -60   -20   -70  -100 /* T */
/* CG G C ..NN */
    60   100    50    20 /* A */
    60   100    50    20 /* C */
    60   100    50    20 /* G */
    60   100    50    20 /* T */
/* CG G G ..NN */
    10    50     0   -30 /* A */
    10    50     0   -30 /* C */
    10    50     0   -30 /* G */
    10    50     0   -30 /* T */
/* CG G T ..NN */
   -30    10   -40   -70 /* A */
   -30    10   -40   -70 /* C */
   -30    10   -40   -70 /* G */
   -30    10   -40   -70 /* T */
/* CG T A ..NN */
   -60   -20   -70  -100 /* A */
   -60   -20   -70  -100 /* C */
   -60   -20   -70  -100 /* G */
   -60   -20   -70  -100 /* T */
/* CG T C ..NN */
    60   100    50    20 /* A */
    60   100    50    20 /* C */
    60   100    50    20 /* G */
    60   100    50    20 /* T */
/* CG T G ..NN */
    10    50     0   -30 /* A */
    10    50     0   -30 /* C */
    10    50     0   -30 /* G */
    10    50     0   -30 /* T */
/* CG T T ..NN */
   -30    10   -40   -70 /* A */
   -30    10   -40   -70 /* C */
   -30    10   -40   -70 /* G */
   -30    10   -40   -70 /* T */
/* GC A A ..CG */
 -1080  -790  -970  -550 /* A */
 -1080  -790  -970  -550 /* C */
 -1080  -790  -970  -550 /* G */
 -1080  -790  -970  -550 /* T */
/* GC A C ..CG */
 -1040  -750  -930  -510 /* A */
 -1040  -750  -930  -510 /* C */
 -1040  -750  -930  -510 /* G */
 -1040  -750  -930  -510 /* T */
/* GC A G ..CG */
 -1090  -800  -980  -560 /* A */
 -1090  -800  -980  -560 /* C */
 -1090  -800  -980  -560 /* G */
 -1090  -800  -980  -560 /* T */
/* GC A T ..CG */
 -1120  -830 -1010  -590 /* A */
 -1120  -830 -1010  -590 /* C */
 -1120  -830 -1010  -590 /* G */
 -1120  -830 -1010  -590 /* T */
/* GC C A ..CG */
 -1080  -790  -970  -550 /* A */
 -1080  -790  -970  -550 /* C */
 -1080  -790  -970  -550 /* G */
 -1080  -790  -970  -550 /* T */
/* GC C C ..CG */
 -1040  -750  -930  -510 /* A */
 -1040  -750  -930  -510 /* C */
 -1040  -750  -930  -510 /* G */
 -1040  -750  -930  -510 /* T */
/* GC C G ..CG */
 -1090  -800  -980  -560 /* A */
 -1090  -800  -980  -560 /* C */
 -1090  -800  -980  -560 /* G */
 -1090  -800  -980  -560 /* T */
/* GC C T ..CG */
 -1120  -830 -1010  -590 /* A */
 -1120  -830 -1010  -590 /* C */
 -1120  -830 -1010  -590 /* G */
 -1120  -830 -1010  -590 /* T */
/* GC G A ..CG */
 -1080  -790  -970  -550 /* A */
 -1080  -790  -970  -550 /* C */
 -1080  -790  -970  -550 /* G */
 -1080  -790  -970  -550 /* T */
/* GC G C ..CG */
 -1040  -750  -930  -510 /* A */
 -1040  -750  -930  -510 /* C */
 -1040  -750  -930  -510 /* G */
 -1040  -750  -930  -510 /* T */
/* GC G G ..CG */
 -1090  -800  -980  -560 /* A */
 -1090  -800  -980  -560 /* C */
 -1090  -800  -980  -560 /* G */
 -1090  -800  -980  -560 /* T */
/* GC G T ..CG */
 -1120  -830 -1010  -590 /* A */
 -1120  -830 -1010  -590 /* C */
 -1120  -830 -1010  -590 /* G */
 -1120  -830 -1010  -590 /* T */
/* GC T A ..CG */
 -1080  -790  -970  -550 /* A */
 -1080  -790  -970  -550 /* C */
 -1080  -790  -970  -550 /* G */
 -1080  -790  -970  -550 /* T */
/* GC T C ..CG */
 -1040  -750  -930  -510 /* A */
 -1040  -750  -930  -510 /* C */
 -1040  -750  -930  -510 /* G */
 -1040  -750  -930  -510 /* T */
/* GC T G ..CG */
 -1090  -800  -980  -560 /* A */
 -1090  -800  -980  -560 /* C */
 -1090  -800  -980  -560 /* G */
 -1090  -800  -980  -560 /* T */
/* GC T T ..CG */
 -1120  -830 -1010  -590 /* A */
 -1120  -830 -1010  -590 /* C */
 -1120  -830 -1010  -590 /* G */
 -1120  -830 -1010  -590 /* T */
/* GC A A ..GC */
 -1240  -950 -1130  -710 /* A */
 -1240  -950 -1130  -710 /* C */
 -1240  -950 -1130  -710 /* G */
 -1240  -950 -1130  -710 /* T */
/* GC A C ..GC */
  -950  -660  -840  -420 /* A */
  -950  -660  -840  -420 /* C */
  -950  -660  -840  -420 /* G */
  -950  -660  -840  -420 /* T */
/* GC A G ..GC */
 -1130  -840 -1020  -600 /* A */
 -1130  -840 -1020  -600 /* C */
 -1130  -840 -1020  -600 /* G */
 -1130  -840 -1020  -600 /* T */
/* GC A T ..GC */
  -710  -420  -600  -180 /* A */
  -710  -420  -600  -180 /* C */
  -710  -420  -600  -180 /* G */
  -710  -420  -600  -180 /* T */
/* GC C A ..GC */
 -1240  -950 -1130  -710 /* A */
 -1240  -950 -1130  -710 /* C */
 -1240  -950 -1130  -710 /* G */
 -1240  -950 -1130  -710 /* T */
/* GC C C ..GC */
  -950  -660  -840  -420 /* A */
  -950  -660  -840  -420 /* C */
  -950  -660  -840  -420 /* G */
  -950  -660  -840  -420 /* T */
/* GC C G ..GC */
 -1130  -840 -1020  -600 /* A */
 -1130  -840 -1020  -600 /* C */
 -1130  -840 -1020  -600 /* G */
 -1130  -840 -1020  -600 /* T */
/* GC C T ..GC */
  -710  -420  -600  -180 /* A */
  -710  -420  -600  -180 /* C */
  -710  -420  -600  -180 /* G */
  -710  -420  -600  -180 /* T */
/* GC G A ..GC */
 -1240  -950 -1130  -710 /* A */
 -1240  -950 -1130  -710 /* C */
 -1240  -950 -1130  -710 /* G */
 -1240  -950 -1130  -710 /* T */
/* GC G C ..GC */
  -950  -660  -840  -420 /* A */
  -950  -660  -840  -420 /* C */
  -950  -660  -840  -420 /* G */
  -950  -660  -840  -420 /* T */
/* GC G G ..GC */
 -1130  -840 -1020  -600 /* A */
 -1130  -840 -1020  -600 /* C */
 -1130  -840 -1020  -600 /* G */
 -1130  -840 -1020  -600 /* T */
/* GC G T ..GC */
  -710  -420  -600  -180 /* A */
  -710  -420  -600  -180 /* C */
  -710  -420  -600  -180 /* G */
  -710  -420  -600  -180 /* T */
/* GC T A ..GC */
 -1240  -950 -1130  -710 /* A */
 -1240  -950 -1130  -710 /* C */
 -1240  -950 -1130  -710 /* G */
 -1240  -950 -1130  -710 /* T */
/* GC T C ..GC */
  -950  -660  -840  -420 /* A */
  -950  -660  -840  -420 /* C */
  -950  -660  -840  -420 /* G */
  -950  -660  -840  -420 /* T */
/* GC T G ..GC */
 -1130  -840 -1020  -600 /* A */
 -1130  -840 -1020  -600 /* C */
 -1130  -840 -1020  -600 /* G */
 -1130  -840 -1020  -600 /* T */
/* GC T T ..GC */
  -710  -420  -600  -180 /* A */
  -710  -420  -600  -180 /* C */
  -710  -420  -600  -180 /* G */
  -710  -420  -600  -180 /* T */
/* GC A A ..GT */
  -920  -630  -810  -390 /* A */
  -920  -630  -810  -390 /* C */
  -920  -630  -810  -390 /* G */
  -920  -630  -810  -390 /* T */
/* GC A C ..GT */
  -630  -340  -520  -100 /* A */
  -630  -340  -520  -100 /* C */
  -630  -340  -520  -100 /* G */
  -630  -340  -520  -100 /* T */
/* GC A G ..GT */
  -810  -520  -700  -280 /* A */
  -810  -520  -700  -280 /* C */
  -810  -520  -700  -280 /* G */
  -810  -520  -700  -280 /* T */
/* GC A T ..GT */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC C A ..GT */
  -920  -630  -810  -390 /* A */
  -920  -630  -810  -390 /* C */
  -920  -630  -810  -390 /* G */
  -920  -630  -810  -390 /* T */
/* GC C C ..GT */
  -630  -340  -520  -100 /* A */
  -630  -340  -520  -100 /* C */
  -630  -340  -520  -100 /* G */
  -630  -340  -520  -100 /* T */
/* GC C G ..GT */
  -810  -520  -700  -280 /* A */
  -810  -520  -700  -280 /* C */
  -810  -520  -700  -280 /* G */
  -810  -520  -700  -280 /* T */
/* GC C T ..GT */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC G A ..GT */
  -920  -630  -810  -390 /* A */
  -920  -630  -810  -390 /* C */
  -920  -630  -810  -390 /* G */
  -920  -630  -810  -390 /* T */
/* GC G C ..GT */
  -630  -340  -520  -100 /* A */
  -630  -340  -520  -100 /* C */
  -630  -340  -520  -100 /* G */
  -630  -340  -520  -100 /* T */
/* GC G G ..GT */
  -810  -520  -700  -280 /* A */
  -810  -520  -700  -280 /* C */
  -810  -520  -700  -280 /* G */
  -810  -520  -700  -280 /* T */
/* GC G T ..GT */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC T A ..GT */
  -920  -630  -810  -390 /* A */
  -920  -630  -810  -390 /* C */
  -920  -630  -810  -390 /* G */
  -920  -630  -810  -390 /* T */
/* GC T C ..GT */
  -630  -340  -520  -100 /* A */
  -630  -340  -520  -100 /* C */
  -630  -340  -520  -100 /* G */
  -630  -340  -520  -100 /* T */
/* GC T G ..GT */
  -810  -520  -700  -280 /* A */
  -810  -520  -700  -280 /* C */
  -810  -520  -700  -280 /* G */
  -810  -520  -700  -280 /* T */
/* GC T T ..GT */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC A A ..TG */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC A C ..TG */
  -360   -70  -250   170 /* A */
  -360   -70  -250   170 /* C */
  -360   -70  -250   170 /* G */
  -360   -70  -250   170 /* T */
/* GC A G ..TG */
  -550  -260  -440   -20 /* A */
  -550  -260  -440   -20 /* C */
  -550  -260  -440   -20 /* G */
  -550  -260  -440   -20 /* T */
/* GC A T ..TG */
 -1200  -910 -1090  -670 /* A */
 -1200  -910 -1090  -670 /* C */
 -1200  -910 -1090  -670 /* G */
 -1200  -910 -1090  -670 /* T */
/* GC C A ..TG */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC C C ..TG */
  -360   -70  -250   170 /* A */
  -360   -70  -250   170 /* C */
  -360   -70  -250   170 /* G */
  -360   -70  -250   170 /* T */
/* GC C G ..TG */
  -550  -260  -440   -20 /* A */
  -550  -260  -440   -20 /* C */
  -550  -260  -440   -20 /* G */
  -550  -260  -440   -20 /* T */
/* GC C T ..TG */
 -1200  -910 -1090  -670 /* A */
 -1200  -910 -1090  -670 /* C */
 -1200  -910 -1090  -670 /* G */
 -1200  -910 -1090  -670 /* T */
/* GC G A ..TG */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC G C ..TG */
  -360   -70  -250   170 /* A */
  -360   -70  -250   170 /* C */
  -360   -70  -250   170 /* G */
  -360   -70  -250   170 /* T */
/* GC G G ..TG */
  -550  -260  -440   -20 /* A */
  -550  -260  -440   -20 /* C */
  -550  -260  -440   -20 /* G */
  -550  -260  -440   -20 /* T */
/* GC G T ..TG */
 -1200  -910 -1090  -670 /* A */
 -1200  -910 -1090  -670 /* C */
 -1200  -910 -1090  -670 /* G */
 -1200  -910 -1090  -670 /* T */
/* GC T A ..TG */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC T C ..TG */
  -360   -70  -250   170 /* A */
  -360   -70  -250   170 /* C */
  -360   -70  -250   170 /* G */
  -360   -70  -250   170 /* T */
/* GC T G ..TG */
  -550  -260  -440   -20 /* A */
  -550  -260  -440   -20 /* C */
  -550  -260  -440   -20 /* G */
  -550  -260  -440   -20 /* T */
/* GC T T ..TG */
 -1200  -910 -1090  -670 /* A */
 -1200  -910 -1090  -670 /* C */
 -1200  -910 -1090  -670 /* G */
 -1200  -910 -1090  -670 /* T */
/* GC A A ..AT */
  -220    70  -110   310 /* A */
  -220    70  -110   310 /* C */
  -220    70  -110   310 /* G */
  -220    70  -110   310 /* T */
/* GC A C ..AT */
  -100   190    10   430 /* A */
  -100   190    10   430 /* C */
  -100   190    10   430 /* G */
  -100   190    10   430 /* T */
/* GC A G ..AT */
  -150   140   -40   380 /* A */
  -150   140   -40   380 /* C */
  -150   140   -40   380 /* G */
  -150   140   -40   380 /* T */
/* GC A T ..AT */
  -190   100   -80   340 /* A */
  -190   100   -80   340 /* C */
  -190   100   -80   340 /* G */
  -190   100   -80   340 /* T */
/* GC C A ..AT */
  -220    70  -110   310 /* A */
  -220    70  -110   310 /* C */
  -220    70  -110   310 /* G */
  -220    70  -110   310 /* T */
/* GC C C ..AT */
  -100   190    10   430 /* A */
  -100   190    10   430 /* C */
  -100   190    10   430 /* G */
  -100   190    10   430 /* T */
/* GC C G ..AT */
  -150   140   -40   380 /* A */
  -150   140   -40   380 /* C */
  -150   140   -40   380 /* G */
  -150   140   -40   380 /* T */
/* GC C T ..AT */
  -190   100   -80   340 /* A */
  -190   100   -80   340 /* C */
  -190   100   -80   340 /* G */
  -190   100   -80   340 /* T */
/* GC G A ..AT */
  -220    70  -110   310 /* A */
  -220    70  -110   310 /* C */
  -220    70  -110   310 /* G */
  -220    70  -110   310 /* T */
/* GC G C ..AT */
  -100   190    10   430 /* A */
  -100   190    10   430 /* C */
  -100   190    10   430 /* G */
  -100   190    10   430 /* T */
/* GC G G ..AT */
  -150   140   -40   380 /* A */
  -150   140   -40   380 /* C */
  -150   140   -40   380 /* G */
  -150   140   -40   380 /* T */
/* GC G T ..AT */
  -190   100   -80   340 /* A */
  -190   100   -80   340 /* C */
  -190   100   -80   340 /* G */
  -190   100   -80   340 /* T */
/* GC T A ..AT */
  -220    70  -110   310 /* A */
  -220    70  -110   310 /* C */
  -220    70  -110   310 /* G */
  -220    70  -110   310 /* T */
/* GC T C ..AT */
  -100   190    10   430 /* A */
  -100   190    10   430 /* C */
  -100   190    10   430 /* G */
  -100   190    10   430 /* T */
/* GC T G ..AT */
  -150   140   -40   380 /* A */
  -150   140   -40   380 /* C */
  -150   140   -40   380 /* G */
  -150   140   -40   380 /* T */
/* GC T T ..AT */
  -190   100   -80   340 /* A */
  -190   100   -80   340 /* C */
  -190   100   -80   340 /* G */
  -190   100   -80   340 /* T */
/* GC A A ..TA */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC A C ..TA */
  -360   -70  -250   170 /* A */
  -360   -70  -250   170 /* C */
  -360   -70  -250   170 /* G */
  -360   -70  -250   170 /* T */
/* GC A G ..TA */
  -550  -260  -440   -20 /* A */
  -550  -260  -440   -20 /* C */
  -550  -260  -440   -20 /* G */
  -550  -260  -440   -20 /* T */
/* GC A T ..TA */
 -1200  -910 -1090  -670 /* A */
 -1200  -910 -1090  -670 /* C */
 -1200  -910 -1090  -670 /* G */
 -1200  -910 -1090  -670 /* T */
/* GC C A ..TA */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC C C ..TA */
  -360   -70  -250   170 /* A */
  -360   -70  -250   170 /* C */
  -360   -70  -250   170 /* G */
  -360   -70  -250   170 /* T */
/* GC C G ..TA */
  -550  -260  -440   -20 /* A */
  -550  -260  -440   -20 /* C */
  -550  -260  -440   -20 /* G */
  -550  -260  -440   -20 /* T */
/* GC C T ..TA */
 -1200  -910 -1090  -670 /* A */
 -1200  -910 -1090  -670 /* C */
 -1200  -910 -1090  -670 /* G */
 -1200  -910 -1090  -670 /* T */
/* GC G A ..TA */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC G C ..TA */
  -360   -70  -250   170 /* A */
  -360   -70  -250   170 /* C */
  -360   -70  -250   170 /* G */
  -360   -70  -250   170 /* T */
/* GC G G ..TA */
  -550  -260  -440   -20 /* A */
  -550  -260  -440   -20 /* C */
  -550  -260  -440   -20 /* G */
  -550  -260  -440   -20 /* T */
/* GC G T ..TA */
 -1200  -910 -1090  -670 /* A */
 -1200  -910 -1090  -670 /* C */
 -1200  -910 -1090  -670 /* G */
 -1200  -910 -1090  -670 /* T */
/* GC T A ..TA */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GC T C ..TA */
  -360   -70  -250   170 /* A */
  -360   -70  -250   170 /* C */
  -360   -70  -250   170 /* G */
  -360   -70  -250   170 /* T */
/* GC T G ..TA */
  -550  -260  -440   -20 /* A */
  -550  -260  -440   -20 /* C */
  -550  -260  -440   -20 /* G */
  -550  -260  -440   -20 /* T */
/* GC T T ..TA */
 -1200  -910 -1090  -670 /* A */
 -1200  -910 -1090  -670 /* C */
 -1200  -910 -1090  -670 /* G */
 -1200  -910 -1090  -670 /* T */
/* GC A A ..NN */
  -220    70  -110   310 /* A */
  -220    70  -110   310 /* C */
  -220    70  -110   310 /* G */
  -220    70  -110   310 /* T */
/* GC A C ..NN */
  -100   190    10   430 /* A */
  -100   190    10   430 /* C */
  -100   190    10   430 /* G */
  -100   190    10   430 /* T */
/* GC A G ..NN */
  -150   140   -40   380 /* A */
  -150   140   -40   380 /* C */
  -150   140   -40   380 /* G */
  -150   140   -40   380 /* T */
/* GC A T ..NN */
  -190   100   -80   340 /* A */
  -190   100   -80   340 /* C */
  -190   100   -80   340 /* G */
  -190   100   -80   340 /* T */
/* GC C A ..NN */
  -220    70  -110   310 /* A */
  -220    70  -110   310 /* C */
  -220    70  -110   310 /* G */
  -220    70  -110   310 /* T */
/* GC C C ..NN */
  -100   190    10   430 /* A */
  -100   190    10   430 /* C */
  -100   190    10   430 /* G */
  -100   190    10   430 /* T */
/* GC C G ..NN */
  -150   140   -40   380 /* A */
  -150   140   -40   380 /* C */
  -150   140   -40   380 /* G */
  -150   140   -40   380 /* T */
/* GC C T ..NN */
  -190   100   -80   340 /* A */
  -190   100   -80   340 /* C */
  -190   100   -80   340 /* G */
  -190   100   -80   340 /* T */
/* GC G A ..NN */
  -220    70  -110   310 /* A */
  -220    70  -110   310 /* C */
  -220    70  -110   310 /* G */
  -220    70  -110   310 /* T */
/* GC G C ..NN */
  -100   190    10   430 /* A */
  -100   190    10   430 /* C */
  -100   190    10   430 /* G */
  -100   190    10   430 /* T */
/* GC G G ..NN */
  -150   140   -40   380 /* A */
  -150   140   -40   380 /* C */
  -150   140   -40   380 /* G */
  -150   140   -40   380 /* T */
/* GC G T ..NN */
  -190   100   -80   340 /* A */
  -190   100   -80   340 /* C */
  -190   100   -80   340 /* G */
  -190   100   -80   340 /* T */
/* GC T A ..NN */
  -220    70  -110   310 /* A */
  -220    70  -110   310 /* C */
  -220    70  -110   310 /* G */
  -220    70  -110   310 /* T */
/* GC T C ..NN */
  -100   190    10   430 /* A */
  -100   190    10   430 /* C */
  -100   190    10   430 /* G */
  -100   190    10   430 /* T */
/* GC T G ..NN */
  -150   140   -40   380 /* A */
  -150   140   -40   380 /* C */
  -150   140   -40   380 /* G */
  -150   140   -40   380 /* T */
/* GC T T ..NN */
  -190   100   -80   340 /* A */
  -190   100   -80   340 /* C */
  -190   100   -80   340 /* G */
  -190   100   -80   340 /* T */
/* GT A A ..CG */
  -760  -470  -650  -230 /* A */
  -760  -470  -650  -230 /* C */
  -760  -470  -650  -230 /* G */
  -760  -470  -650  -230 /* T */
/* GT A C ..CG */
  -720  -430  -610  -190 /* A */
  -720  -430  -610  -190 /* C */
  -720  -430  -610  -190 /* G */
  -720  -430  -610  -190 /* T */
/* GT A G ..CG */
  -770  -480  -660  -240 /* A */
  -770  -480  -660  -240 /* C */
  -770  -480  -660  -240 /* G */
  -770  -480  -660  -240 /* T */
/* GT A T ..CG */
  -800  -510  -690  -270 /* A */
  -800  -510  -690  -270 /* C */
  -800  -510  -690  -270 /* G */
  -800  -510  -690  -270 /* T */
/* GT C A ..CG */
  -760  -470  -650  -230 /* A */
  -760  -470  -650  -230 /* C */
  -760  -470  -650  -230 /* G */
  -760  -470  -650  -230 /* T */
/* GT C C ..CG */
  -720  -430  -610  -190 /* A */
  -720  -430  -610  -190 /* C */
  -720  -430  -610  -190 /* G */
  -720  -430  -610  -190 /* T */
/* GT C G ..CG */
  -770  -480  -660  -240 /* A */
  -770  -480  -660  -240 /* C */
  -770  -480  -660  -240 /* G */
  -770  -480  -660  -240 /* T */
/* GT C T ..CG */
  -800  -510  -690  -270 /* A */
  -800  -510  -690  -270 /* C */
  -800  -510  -690  -270 /* G */
  -800  -510  -690  -270 /* T */
/* GT G A ..CG */
  -760  -470  -650  -230 /* A */
  -760  -470  -650  -230 /* C */
  -760  -470  -650  -230 /* G */
  -760  -470  -650  -230 /* T */
/* GT G C ..CG */
  -720  -430  -610  -190 /* A */
  -720  -430  -610  -190 /* C */
  -720  -430  -610  -190 /* G */
  -720  -430  -610  -190 /* T */
/* GT G G ..CG */
  -770  -480  -660  -240 /* A */
  -770  -480  -660  -240 /* C */
  -770  -480  -660  -240 /* G */
  -770  -480  -660  -240 /* T */
/* GT G T ..CG */
  -800  -510  -690  -270 /* A */
  -800  -510  -690  -270 /* C */
  -800  -510  -690  -270 /* G */
  -800  -510  -690  -270 /* T */
/* GT T A ..CG */
  -760  -470  -650  -230 /* A */
  -760  -470  -650  -230 /* C */
  -760  -470  -650  -230 /* G */
  -760  -470  -650  -230 /* T */
/* GT T C ..CG */
  -720  -430  -610  -190 /* A */
  -720  -430  -610  -190 /* C */
  -720  -430  -610  -190 /* G */
  -720  -430  -610  -190 /* T */
/* GT T G ..CG */
  -770  -480  -660  -240 /* A */
  -770  -480  -660  -240 /* C */
  -770  -480  -660  -240 /* G */
  -770  -480  -660  -240 /* T */
/* GT T T ..CG */
  -800  -510  -690  -270 /* A */
  -800  -510  -690  -270 /* C */
  -800  -510  -690  -270 /* G */
  -800  -510  -690  -270 /* T */
/* GT A A ..GC */
  -920  -630  -810  -390 /* A */
  -920  -630  -810  -390 /* C */
  -920  -630  -810  -390 /* G */
  -920  -630  -810  -390 /* T */
/* GT A C ..GC */
  -630  -340  -520  -100 /* A */
  -630  -340  -520  -100 /* C */
  -630  -340  -520  -100 /* G */
  -630  -340  -520  -100 /* T */
/* GT A G ..GC */
  -810  -520  -700  -280 /* A */
  -810  -520  -700  -280 /* C */
  -810  -520  -700  -280 /* G */
  -810  -520  -700  -280 /* T */
/* GT A T ..GC */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GT C A ..GC */
  -920  -630  -810  -390 /* A */
  -920  -630  -810  -390 /* C */
  -920  -630  -810  -390 /* G */
  -920  -630  -810  -390 /* T */
/* GT C C ..GC */
  -630  -340  -520  -100 /* A */
  -630  -340  -520  -100 /* C */
  -630  -340  -520  -100 /* G */
  -630  -340  -520  -100 /* T */
/* GT C G ..GC */
  -810  -520  -700  -280 /* A */
  -810  -520  -700  -280 /* C */
  -810  -520  -700  -280 /* G */
  -810  -520  -700  -280 /* T */
/* GT C T ..GC */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GT G A ..GC */
  -920  -630  -810  -390 /* A */
  -920  -630  -810  -390 /* C */
  -920  -630  -810  -390 /* G */
  -920  -630  -810  -390 /* T */
/* GT G C ..GC */
  -630  -340  -520  -100 /* A */
  -630  -340  -520  -100 /* C */
  -630  -340  -520  -100 /* G */
  -630  -340  -520  -100 /* T */
/* GT G G ..GC */
  -810  -520  -700  -280 /* A */
  -810  -520  -700  -280 /* C */
  -810  -520  -700  -280 /* G */
  -810  -520  -700  -280 /* T */
/* GT G T ..GC */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GT T A ..GC */
  -920  -630  -810  -390 /* A */
  -920  -630  -810  -390 /* C */
  -920  -630  -810  -390 /* G */
  -920  -630  -810  -390 /* T */
/* GT T C ..GC */
  -630  -340  -520  -100 /* A */
  -630  -340  -520  -100 /* C */
  -630  -340  -520  -100 /* G */
  -630  -340  -520  -100 /* T */
/* GT T G ..GC */
  -810  -520  -700  -280 /* A */
  -810  -520  -700  -280 /* C */
  -810  -520  -700  -280 /* G */
  -810  -520  -700  -280 /* T */
/* GT T T ..GC */
  -390  -100  -280   140 /* A */
  -390  -100  -280   140 /* C */
  -390  -100  -280   140 /* G */
  -390  -100  -280   140 /* T */
/* GT A A ..GT */
  -600  -310  -490   -70 /* A */
  -600  -310  -490   -70 /* C */
  -600  -310  -490   -70 /* G */
  -600  -310  -490   -70 /* T */
/* GT A C ..GT */
  -310   -20  -200   220 /* A */
  -310   -20  -200   220 /* C */
  -310   -20  -200   220 /* G */
  -310   -20  -200   220 /* T */
/* GT A G ..GT */
  -490  -200  -380    40 /* A */
  -490  -200  -380    40 /* C */
  -490  -200  -380    40 /* G */
  -490  -200  -380    40 /* T */
/* GT A T ..GT */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT C A ..GT */
  -600  -310  -490   -70 /* A */
  -600  -310  -490   -70 /* C */
  -600  -310  -490   -70 /* G */
  -600  -310  -490   -70 /* T */
/* GT C C ..GT */
  -310   -20  -200   220 /* A */
  -310   -20  -200   220 /* C */
  -310   -20  -200   220 /* G */
  -310   -20  -200   220 /* T */
/* GT C G ..GT */
  -490  -200  -380    40 /* A */
  -490  -200  -380    40 /* C */
  -490  -200  -380    40 /* G */
  -490  -200  -380    40 /* T */
/* GT C T ..GT */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT G A ..GT */
  -600  -310  -490   -70 /* A */
  -600  -310  -490   -70 /* C */
  -600  -310  -490   -70 /* G */
  -600  -310  -490   -70 /* T */
/* GT G C ..GT */
  -310   -20  -200   220 /* A */
  -310   -20  -200   220 /* C */
  -310   -20  -200   220 /* G */
  -310   -20  -200   220 /* T */
/* GT G G ..GT */
  -490  -200  -380    40 /* A */
  -490  -200  -380    40 /* C */
  -490  -200  -380    40 /* G */
  -490  -200  -380    40 /* T */
/* GT G T ..GT */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT T A ..GT */
  -600  -310  -490   -70 /* A */
  -600  -310  -490   -70 /* C */
  -600  -310  -490   -70 /* G */
  -600  -310  -490   -70 /* T */
/* GT T C ..GT */
  -310   -20  -200   220 /* A */
  -310   -20  -200   220 /* C */
  -310   -20  -200   220 /* G */
  -310   -20  -200   220 /* T */
/* GT T G ..GT */
  -490  -200  -380    40 /* A */
  -490  -200  -380    40 /* C */
  -490  -200  -380    40 /* G */
  -490  -200  -380    40 /* T */
/* GT T T ..GT */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT A A ..TG */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT A C ..TG */
   -40   250    70   490 /* A */
   -40   250    70   490 /* C */
   -40   250    70   490 /* G */
   -40   250    70   490 /* T */
/* GT A G ..TG */
  -230    60  -120   300 /* A */
  -230    60  -120   300 /* C */
  -230    60  -120   300 /* G */
  -230    60  -120   300 /* T */
/* GT A T ..TG */
  -880  -590  -770  -350 /* A */
  -880  -590  -770  -350 /* C */
  -880  -590  -770  -350 /* G */
  -880  -590  -770  -350 /* T */
/* GT C A ..TG */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT C C ..TG */
   -40   250    70   490 /* A */
   -40   250    70   490 /* C */
   -40   250    70   490 /* G */
   -40   250    70   490 /* T */
/* GT C G ..TG */
  -230    60  -120   300 /* A */
  -230    60  -120   300 /* C */
  -230    60  -120   300 /* G */
  -230    60  -120   300 /* T */
/* GT C T ..TG */
  -880  -590  -770  -350 /* A */
  -880  -590  -770  -350 /* C */
  -880  -590  -770  -350 /* G */
  -880  -590  -770  -350 /* T */
/* GT G A ..TG */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT G C ..TG */
   -40   250    70   490 /* A */
   -40   250    70   490 /* C */
   -40   250    70   490 /* G */
   -40   250    70   490 /* T */
/* GT G G ..TG */
  -230    60  -120   300 /* A */
  -230    60  -120   300 /* C */
  -230    60  -120   300 /* G */
  -230    60  -120   300 /* T */
/* GT G T ..TG */
  -880  -590  -770  -350 /* A */
  -880  -590  -770  -350 /* C */
  -880  -590  -770  -350 /* G */
  -880  -590  -770  -350 /* T */
/* GT T A ..TG */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT T C ..TG */
   -40   250    70   490 /* A */
   -40   250    70   490 /* C */
   -40   250    70   490 /* G */
   -40   250    70   490 /* T */
/* GT T G ..TG */
  -230    60  -120   300 /* A */
  -230    60  -120   300 /* C */
  -230    60  -120   300 /* G */
  -230    60  -120   300 /* T */
/* GT T T ..TG */
  -880  -590  -770  -350 /* A */
  -880  -590  -770  -350 /* C */
  -880  -590  -770  -350 /* G */
  -880  -590  -770  -350 /* T */
/* GT A A ..AT */
   100   390   210   630 /* A */
   100   390   210   630 /* C */
   100   390   210   630 /* G */
   100   390   210   630 /* T */
/* GT A C ..AT */
   220   510   330   750 /* A */
   220   510   330   750 /* C */
   220   510   330   750 /* G */
   220   510   330   750 /* T */
/* GT A G ..AT */
   170   460   280   700 /* A */
   170   460   280   700 /* C */
   170   460   280   700 /* G */
   170   460   280   700 /* T */
/* GT A T ..AT */
   130   420   240   660 /* A */
   130   420   240   660 /* C */
   130   420   240   660 /* G */
   130   420   240   660 /* T */
/* GT C A ..AT */
   100   390   210   630 /* A */
   100   390   210   630 /* C */
   100   390   210   630 /* G */
   100   390   210   630 /* T */
/* GT C C ..AT */
   220   510   330   750 /* A */
   220   510   330   750 /* C */
   220   510   330   750 /* G */
   220   510   330   750 /* T */
/* GT C G ..AT */
   170   460   280   700 /* A */
   170   460   280   700 /* C */
   170   460   280   700 /* G */
   170   460   280   700 /* T */
/* GT C T ..AT */
   130   420   240   660 /* A */
   130   420   240   660 /* C */
   130   420   240   660 /* G */
   130   420   240   660 /* T */
/* GT G A ..AT */
   100   390   210   630 /* A */
   100   390   210   630 /* C */
   100   390   210   630 /* G */
   100   390   210   630 /* T */
/* GT G C ..AT */
   220   510   330   750 /* A */
   220   510   330   750 /* C */
   220   510   330   750 /* G */
   220   510   330   750 /* T */
/* GT G G ..AT */
   170   460   280   700 /* A */
   170   460   280   700 /* C */
   170   460   280   700 /* G */
   170   460   280   700 /* T */
/* GT G T ..AT */
   130   420   240   660 /* A */
   130   420   240   660 /* C */
   130   420   240   660 /* G */
   130   420   240   660 /* T */
/* GT T A ..AT */
   100   390   210   630 /* A */
   100   390   210   630 /* C */
   100   390   210   630 /* G */
   100   390   210   630 /* T */
/* GT T C ..AT */
   220   510   330   750 /* A */
   220   510   330   750 /* C */
   220   510   330   750 /* G */
   220   510   330   750 /* T */
/* GT T G ..AT */
   170   460   280   700 /* A */
   170   460   280   700 /* C */
   170   460   280   700 /* G */
   170   460   280   700 /* T */
/* GT T T ..AT */
   130   420   240   660 /* A */
   130   420   240   660 /* C */
   130   420   240   660 /* G */
   130   420   240   660 /* T */
/* GT A A ..TA */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT A C ..TA */
   -40   250    70   490 /* A */
   -40   250    70   490 /* C */
   -40   250    70   490 /* G */
   -40   250    70   490 /* T */
/* GT A G ..TA */
  -230    60  -120   300 /* A */
  -230    60  -120   300 /* C */
  -230    60  -120   300 /* G */
  -230    60  -120   300 /* T */
/* GT A T ..TA */
  -880  -590  -770  -350 /* A */
  -880  -590  -770  -350 /* C */
  -880  -590  -770  -350 /* G */
  -880  -590  -770  -350 /* T */
/* GT C A ..TA */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT C C ..TA */
   -40   250    70   490 /* A */
   -40   250    70   490 /* C */
   -40   250    70   490 /* G */
   -40   250    70   490 /* T */
/* GT C G ..TA */
  -230    60  -120   300 /* A */
  -230    60  -120   300 /* C */
  -230    60  -120   300 /* G */
  -230    60  -120   300 /* T */
/* GT C T ..TA */
  -880  -590  -770  -350 /* A */
  -880  -590  -770  -350 /* C */
  -880  -590  -770  -350 /* G */
  -880  -590  -770  -350 /* T */
/* GT G A ..TA */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT G C ..TA */
   -40   250    70   490 /* A */
   -40   250    70   490 /* C */
   -40   250    70   490 /* G */
   -40   250    70   490 /* T */
/* GT G G ..TA */
  -230    60  -120   300 /* A */
  -230    60  -120   300 /* C */
  -230    60  -120   300 /* G */
  -230    60  -120   300 /* T */
/* GT G T ..TA */
  -880  -590  -770  -350 /* A */
  -880  -590  -770  -350 /* C */
  -880  -590  -770  -350 /* G */
  -880  -590  -770  -350 /* T */
/* GT T A ..TA */
   -70   220    40   460 /* A */
   -70   220    40   460 /* C */
   -70   220    40   460 /* G */
   -70   220    40   460 /* T */
/* GT T C ..TA */
   -40   250    70   490 /* A */
   -40   250    70   490 /* C */
   -40   250    70   490 /* G */
   -40   250    70   490 /* T */
/* GT T G ..TA */
  -230    60  -120   300 /* A */
  -230    60  -120   300 /* C */
  -230    60  -120   300 /* G */
  -230    60  -120   300 /* T */
/* GT T T ..TA */
  -880  -590  -770  -350 /* A */
  -880  -590  -770  -350 /* C */
  -880  -590  -770  -350 /* G */
  -880  -590  -770  -350 /* T */
/* GT A A ..NN */
   100   390   210   630 /* A */
   100   390   210   630 /* C */
   100   390   210   630 /* G */
   100   390   210   630 /* T */
/* GT A C ..NN */
   220   510   330   750 /* A */
   220   510   330   750 /* C */
   220   510   330   750 /* G */
   220   510   330   750 /* T */
/* GT A G ..NN */
   170   460   280   700 /* A */
   170   460   280   700 /* C */
   170   460   280   700 /* G */
   170   460   280   700 /* T */
/* GT A T ..NN */
   130   420   240   660 /* A */
   130   420   240   660 /* C */
   130   420   240   660 /* G */
   130   420   240   660 /* T */
/* GT C A ..NN */
   100   390   210   630 /* A */
   100   390   210   630 /* C */
   100   390   210   630 /* G */
   100   390   210   630 /* T */
/* GT C C ..NN */
   220   510   330   750 /* A */
   220   510   330   750 /* C */
   220   510   330   750 /* G */
   220   510   330   750 /* T */
/* GT C G ..NN */
   170   460   280   700 /* A */
   170   460   280   700 /* C */
   170   460   280   700 /* G */
   170   460   280   700 /* T */
/* GT C T ..NN */
   130   420   240   660 /* A */
   130   420   240   660 /* C */
   130   420   240   660 /* G */
   130   420   240   660 /* T */
/* GT G A ..NN */
   100   390   210   630 /* A */
   100   390   210   630 /* C */
   100   390   210   630 /* G */
   100   390   210   630 /* T */
/* GT G C ..NN */
   220   510   330   750 /* A */
   220   510   330   750 /* C */
   220   510   330   750 /* G */
   220   510   330   750 /* T */
/* GT G G ..NN */
   170   460   280   700 /* A */
   170   460   280   700 /* C */
   170   460   280   700 /* G */
   170   460   280   700 /* T */
/* GT G T ..NN */
   130   420   240   660 /* A */
   130   420   240   660 /* C */
   130   420   240   660 /* G */
   130   420   240   660 /* T */
/* GT T A ..NN */
   100   390   210   630 /* A */
   100   390   210   630 /* C */
   100   390   210   630 /* G */
   100   390   210   630 /* T */
/* GT T C ..NN */
   220   510   330   750 /* A */
   220   510   330   750 /* C */
   220   510   330   750 /* G */
   220   510   330   750 /* T */
/* GT T G ..NN */
   170   460   280   700 /* A */
   170   460   280   700 /* C */
   170   460   280   700 /* G */
   170   460   280   700 /* T */
/* GT T T ..NN */
   130   420   240   660 /* A */
   130   420   240   660 /* C */
   130   420   240   660 /* G */
   130   420   240   660 /* T */
/* TG A A ..CG */
  -230  -200  -390 -1040 /* A */
  -230  -200  -390 -1040 /* C */
  -230  -200  -390 -1040 /* G */
  -230  -200  -390 -1040 /* T */
/* TG A C ..CG */
  -190  -160  -350 -1000 /* A */
  -190  -160  -350 -1000 /* C */
  -190  -160  -350 -1000 /* G */
  -190  -160  -350 -1000 /* T */
/* TG A G ..CG */
  -240  -210  -400 -1050 /* A */
  -240  -210  -400 -1050 /* C */
  -240  -210  -400 -1050 /* G */
  -240  -210  -400 -1050 /* T */
/* TG A T ..CG */
  -270  -240  -430 -1080 /* A */
  -270  -240  -430 -1080 /* C */
  -270  -240  -430 -1080 /* G */
  -270  -240  -430 -1080 /* T */
/* TG C A ..CG */
  -230  -200  -390 -1040 /* A */
  -230  -200  -390 -1040 /* C */
  -230  -200  -390 -1040 /* G */
  -230  -200  -390 -1040 /* T */
/* TG C C ..CG */
  -190  -160  -350 -1000 /* A */
  -190  -160  -350 -1000 /* C */
  -190  -160  -350 -1000 /* G */
  -190  -160  -350 -1000 /* T */
/* TG C G ..CG */
  -240  -210  -400 -1050 /* A */
  -240  -210  -400 -1050 /* C */
  -240  -210  -400 -1050 /* G */
  -240  -210  -400 -1050 /* T */
/* TG C T ..CG */
  -270  -240  -430 -1080 /* A */
  -270  -240  -430 -1080 /* C */
  -270  -240  -430 -1080 /* G */
  -270  -240  -430 -1080 /* T */
/* TG G A ..CG */
  -230  -200  -390 -1040 /* A */
  -230  -200  -390 -1040 /* C */
  -230  -200  -390 -1040 /* G */
  -230  -200  -390 -1040 /* T */
/* TG G C ..CG */
  -190  -160  -350 -1000 /* A */
  -190  -160  -350 -1000 /* C */
  -190  -160  -350 -1000 /* G */
  -190  -160  -350 -1000 /* T */
/* TG G G ..CG */
  -240  -210  -400 -1050 /* A */
  -240  -210  -400 -1050 /* C */
  -240  -210  -400 -1050 /* G */
  -240  -210  -400 -1050 /* T */
/* TG G T ..CG */
  -270  -240  -430 -1080 /* A */
  -270  -240  -430 -1080 /* C */
  -270  -240  -430 -1080 /* G */
  -270  -240  -430 -1080 /* T */
/* TG T A ..CG */
  -230  -200  -390 -1040 /* A */
  -230  -200  -390 -1040 /* C */
  -230  -200  -390 -1040 /* G */
  -230  -200  -390 -1040 /* T */
/* TG T C ..CG */
  -190  -160  -350 -1000 /* A */
  -190  -160  -350 -1000 /* C */
  -190  -160  -350 -1000 /* G */
  -190  -160  -350 -1000 /* T */
/* TG T G ..CG */
  -240  -210  -400 -1050 /* A */
  -240  -210  -400 -1050 /* C */
  -240  -210  -400 -1050 /* G */
  -240  -210  -400 -1050 /* T */
/* TG T T ..CG */
  -270  -240  -430 -1080 /* A */
  -270  -240  -430 -1080 /* C */
  -270  -240  -430 -1080 /* G */
  -270  -240  -430 -1080 /* T */
/* TG A A ..GC */
  -390  -360  -550 -1200 /* A */
  -390  -360  -550 -1200 /* C */
  -390  -360  -550 -1200 /* G */
  -390  -360  -550 -1200 /* T */
/* TG A C ..GC */
  -100   -70  -260  -910 /* A */
  -100   -70  -260  -910 /* C */
  -100   -70  -260  -910 /* G */
  -100   -70  -260  -910 /* T */
/* TG A G ..GC */
  -280  -250  -440 -1090 /* A */
  -280  -250  -440 -1090 /* C */
  -280  -250  -440 -1090 /* G */
  -280  -250  -440 -1090 /* T */
/* TG A T ..GC */
   140   170   -20  -670 /* A */
   140   170   -20  -670 /* C */
   140   170   -20  -670 /* G */
   140   170   -20  -670 /* T */
/* TG C A ..GC */
  -390  -360  -550 -1200 /* A */
  -390  -360  -550 -1200 /* C */
  -390  -360  -550 -1200 /* G */
  -390  -360  -550 -1200 /* T */
/* TG C C ..GC */
  -100   -70  -260  -910 /* A */
  -100   -70  -260  -910 /* C */
  -100   -70  -260  -910 /* G */
  -100   -70  -260  -910 /* T */
/* TG C G ..GC */
  -280  -250  -440 -1090 /* A */
  -280  -250  -440 -1090 /* C */
  -280  -250  -440 -1090 /* G */
  -280  -250  -440 -1090 /* T */
/* TG C T ..GC */
   140   170   -20  -670 /* A */
   140   170   -20  -670 /* C */
   140   170   -20  -670 /* G */
   140   170   -20  -670 /* T */
/* TG G A ..GC */
  -390  -360  -550 -1200 /* A */
  -390  -360  -550 -1200 /* C */
  -390  -360  -550 -1200 /* G */
  -390  -360  -550 -1200 /* T */
/* TG G C ..GC */
  -100   -70  -260  -910 /* A */
  -100   -70  -260  -910 /* C */
  -100   -70  -260  -910 /* G */
  -100   -70  -260  -910 /* T */
/* TG G G ..GC */
  -280  -250  -440 -1090 /* A */
  -280  -250  -440 -1090 /* C */
  -280  -250  -440 -1090 /* G */
  -280  -250  -440 -1090 /* T */
/* TG G T ..GC */
   140   170   -20  -670 /* A */
   140   170   -20  -670 /* C */
   140   170   -20  -670 /* G */
   140   170   -20  -670 /* T */
/* TG T A ..GC */
  -390  -360  -550 -1200 /* A */
  -390  -360  -550 -1200 /* C */
  -390  -360  -550 -1200 /* G */
  -390  -360  -550 -1200 /* T */
/* TG T C ..GC */
  -100   -70  -260  -910 /* A */
  -100   -70  -260  -910 /* C */
  -100   -70  -260  -910 /* G */
  -100   -70  -260  -910 /* T */
/* TG T G ..GC */
  -280  -250  -440 -1090 /* A */
  -280  -250  -440 -1090 /* C */
  -280  -250  -440 -1090 /* G */
  -280  -250  -440 -1090 /* T */
/* TG T T ..GC */
   140   170   -20  -670 /* A */
   140   170   -20  -670 /* C */
   140   170   -20  -670 /* G */
   140   170   -20  -670 /* T */
/* TG A A ..GT */
   -70   -40  -230  -880 /* A */
   -70   -40  -230  -880 /* C */
   -70   -40  -230  -880 /* G */
   -70   -40  -230  -880 /* T */
/* TG A C ..GT */
   220   250    60  -590 /* A */
   220   250    60  -590 /* C */
   220   250    60  -590 /* G */
   220   250    60  -590 /* T */
/* TG A G ..GT */
    40    70  -120  -770 /* A */
    40    70  -120  -770 /* C */
    40    70  -120  -770 /* G */
    40    70  -120  -770 /* T */
/* TG A T ..GT */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG C A ..GT */
   -70   -40  -230  -880 /* A */
   -70   -40  -230  -880 /* C */
   -70   -40  -230  -880 /* G */
   -70   -40  -230  -880 /* T */
/* TG C C ..GT */
   220   250    60  -590 /* A */
   220   250    60  -590 /* C */
   220   250    60  -590 /* G */
   220   250    60  -590 /* T */
/* TG C G ..GT */
    40    70  -120  -770 /* A */
    40    70  -120  -770 /* C */
    40    70  -120  -770 /* G */
    40    70  -120  -770 /* T */
/* TG C T ..GT */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG G A ..GT */
   -70   -40  -230  -880 /* A */
   -70   -40  -230  -880 /* C */
   -70   -40  -230  -880 /* G */
   -70   -40  -230  -880 /* T */
/* TG G C ..GT */
   220   250    60  -590 /* A */
   220   250    60  -590 /* C */
   220   250    60  -590 /* G */
   220   250    60  -590 /* T */
/* TG G G ..GT */
    40    70  -120  -770 /* A */
    40    70  -120  -770 /* C */
    40    70  -120  -770 /* G */
    40    70  -120  -770 /* T */
/* TG G T ..GT */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG T A ..GT */
   -70   -40  -230  -880 /* A */
   -70   -40  -230  -880 /* C */
   -70   -40  -230  -880 /* G */
   -70   -40  -230  -880 /* T */
/* TG T C ..GT */
   220   250    60  -590 /* A */
   220   250    60  -590 /* C */
   220   250    60  -590 /* G */
   220   250    60  -590 /* T */
/* TG T G ..GT */
    40    70  -120  -770 /* A */
    40    70  -120  -770 /* C */
    40    70  -120  -770 /* G */
    40    70  -120  -770 /* T */
/* TG T T ..GT */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG A A ..TG */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG A C ..TG */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TG A G ..TG */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TG A T ..TG */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TG C A ..TG */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG C C ..TG */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TG C G ..TG */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TG C T ..TG */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TG G A ..TG */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG G C ..TG */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TG G G ..TG */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TG G T ..TG */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TG T A ..TG */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG T C ..TG */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TG T G ..TG */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TG T T ..TG */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TG A A ..AT */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TG A C ..AT */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TG A G ..AT */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TG A T ..AT */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TG C A ..AT */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TG C C ..AT */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TG C G ..AT */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TG C T ..AT */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TG G A ..AT */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TG G C ..AT */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TG G G ..AT */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TG G T ..AT */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TG T A ..AT */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TG T C ..AT */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TG T G ..AT */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TG T T ..AT */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TG A A ..TA */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG A C ..TA */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TG A G ..TA */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TG A T ..TA */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TG C A ..TA */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG C C ..TA */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TG C G ..TA */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TG C T ..TA */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TG G A ..TA */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG G C ..TA */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TG G G ..TA */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TG G T ..TA */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TG T A ..TA */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TG T C ..TA */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TG T G ..TA */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TG T T ..TA */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TG A A ..NN */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TG A C ..NN */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TG A G ..NN */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TG A T ..NN */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TG C A ..NN */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TG C C ..NN */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TG C G ..NN */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TG C T ..NN */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TG G A ..NN */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TG G C ..NN */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TG G G ..NN */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TG G T ..NN */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TG T A ..NN */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TG T C ..NN */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TG T G ..NN */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TG T T ..NN */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* AT A A ..CG */
   -60    60    10   -30 /* A */
   -60    60    10   -30 /* C */
   -60    60    10   -30 /* G */
   -60    60    10   -30 /* T */
/* AT A C ..CG */
   -20   100    50    10 /* A */
   -20   100    50    10 /* C */
   -20   100    50    10 /* G */
   -20   100    50    10 /* T */
/* AT A G ..CG */
   -70    50     0   -40 /* A */
   -70    50     0   -40 /* C */
   -70    50     0   -40 /* G */
   -70    50     0   -40 /* T */
/* AT A T ..CG */
  -100    20   -30   -70 /* A */
  -100    20   -30   -70 /* C */
  -100    20   -30   -70 /* G */
  -100    20   -30   -70 /* T */
/* AT C A ..CG */
   -60    60    10   -30 /* A */
   -60    60    10   -30 /* C */
   -60    60    10   -30 /* G */
   -60    60    10   -30 /* T */
/* AT C C ..CG */
   -20   100    50    10 /* A */
   -20   100    50    10 /* C */
   -20   100    50    10 /* G */
   -20   100    50    10 /* T */
/* AT C G ..CG */
   -70    50     0   -40 /* A */
   -70    50     0   -40 /* C */
   -70    50     0   -40 /* G */
   -70    50     0   -40 /* T */
/* AT C T ..CG */
  -100    20   -30   -70 /* A */
  -100    20   -30   -70 /* C */
  -100    20   -30   -70 /* G */
  -100    20   -30   -70 /* T */
/* AT G A ..CG */
   -60    60    10   -30 /* A */
   -60    60    10   -30 /* C */
   -60    60    10   -30 /* G */
   -60    60    10   -30 /* T */
/* AT G C ..CG */
   -20   100    50    10 /* A */
   -20   100    50    10 /* C */
   -20   100    50    10 /* G */
   -20   100    50    10 /* T */
/* AT G G ..CG */
   -70    50     0   -40 /* A */
   -70    50     0   -40 /* C */
   -70    50     0   -40 /* G */
   -70    50     0   -40 /* T */
/* AT G T ..CG */
  -100    20   -30   -70 /* A */
  -100    20   -30   -70 /* C */
  -100    20   -30   -70 /* G */
  -100    20   -30   -70 /* T */
/* AT T A ..CG */
   -60    60    10   -30 /* A */
   -60    60    10   -30 /* C */
   -60    60    10   -30 /* G */
   -60    60    10   -30 /* T */
/* AT T C ..CG */
   -20   100    50    10 /* A */
   -20   100    50    10 /* C */
   -20   100    50    10 /* G */
   -20   100    50    10 /* T */
/* AT T G ..CG */
   -70    50     0   -40 /* A */
   -70    50     0   -40 /* C */
   -70    50     0   -40 /* G */
   -70    50     0   -40 /* T */
/* AT T T ..CG */
  -100    20   -30   -70 /* A */
  -100    20   -30   -70 /* C */
  -100    20   -30   -70 /* G */
  -100    20   -30   -70 /* T */
/* AT A A ..GC */
  -220  -100  -150  -190 /* A */
  -220  -100  -150  -190 /* C */
  -220  -100  -150  -190 /* G */
  -220  -100  -150  -190 /* T */
/* AT A C ..GC */
    70   190   140   100 /* A */
    70   190   140   100 /* C */
    70   190   140   100 /* G */
    70   190   140   100 /* T */
/* AT A G ..GC */
  -110    10   -40   -80 /* A */
  -110    10   -40   -80 /* C */
  -110    10   -40   -80 /* G */
  -110    10   -40   -80 /* T */
/* AT A T ..GC */
   310   430   380   340 /* A */
   310   430   380   340 /* C */
   310   430   380   340 /* G */
   310   430   380   340 /* T */
/* AT C A ..GC */
  -220  -100  -150  -190 /* A */
  -220  -100  -150  -190 /* C */
  -220  -100  -150  -190 /* G */
  -220  -100  -150  -190 /* T */
/* AT C C ..GC */
    70   190   140   100 /* A */
    70   190   140   100 /* C */
    70   190   140   100 /* G */
    70   190   140   100 /* T */
/* AT C G ..GC */
  -110    10   -40   -80 /* A */
  -110    10   -40   -80 /* C */
  -110    10   -40   -80 /* G */
  -110    10   -40   -80 /* T */
/* AT C T ..GC */
   310   430   380   340 /* A */
   310   430   380   340 /* C */
   310   430   380   340 /* G */
   310   430   380   340 /* T */
/* AT G A ..GC */
  -220  -100  -150  -190 /* A */
  -220  -100  -150  -190 /* C */
  -220  -100  -150  -190 /* G */
  -220  -100  -150  -190 /* T */
/* AT G C ..GC */
    70   190   140   100 /* A */
    70   190   140   100 /* C */
    70   190   140   100 /* G */
    70   190   140   100 /* T */
/* AT G G ..GC */
  -110    10   -40   -80 /* A */
  -110    10   -40   -80 /* C */
  -110    10   -40   -80 /* G */
  -110    10   -40   -80 /* T */
/* AT G T ..GC */
   310   430   380   340 /* A */
   310   430   380   340 /* C */
   310   430   380   340 /* G */
   310   430   380   340 /* T */
/* AT T A ..GC */
  -220  -100  -150  -190 /* A */
  -220  -100  -150  -190 /* C */
  -220  -100  -150  -190 /* G */
  -220  -100  -150  -190 /* T */
/* AT T C ..GC */
    70   190   140   100 /* A */
    70   190   140   100 /* C */
    70   190   140   100 /* G */
    70   190   140   100 /* T */
/* AT T G ..GC */
  -110    10   -40   -80 /* A */
  -110    10   -40   -80 /* C */
  -110    10   -40   -80 /* G */
  -110    10   -40   -80 /* T */
/* AT T T ..GC */
   310   430   380   340 /* A */
   310   430   380   340 /* C */
   310   430   380   340 /* G */
   310   430   380   340 /* T */
/* AT A A ..GT */
   100   220   170   130 /* A */
   100   220   170   130 /* C */
   100   220   170   130 /* G */
   100   220   170   130 /* T */
/* AT A C ..GT */
   390   510   460   420 /* A */
   390   510   460   420 /* C */
   390   510   460   420 /* G */
   390   510   460   420 /* T */
/* AT A G ..GT */
   210   330   280   240 /* A */
   210   330   280   240 /* C */
   210   330   280   240 /* G */
   210   330   280   240 /* T */
/* AT A T ..GT */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT C A ..GT */
   100   220   170   130 /* A */
   100   220   170   130 /* C */
   100   220   170   130 /* G */
   100   220   170   130 /* T */
/* AT C C ..GT */
   390   510   460   420 /* A */
   390   510   460   420 /* C */
   390   510   460   420 /* G */
   390   510   460   420 /* T */
/* AT C G ..GT */
   210   330   280   240 /* A */
   210   330   280   240 /* C */
   210   330   280   240 /* G */
   210   330   280   240 /* T */
/* AT C T ..GT */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT G A ..GT */
   100   220   170   130 /* A */
   100   220   170   130 /* C */
   100   220   170   130 /* G */
   100   220   170   130 /* T */
/* AT G C ..GT */
   390   510   460   420 /* A */
   390   510   460   420 /* C */
   390   510   460   420 /* G */
   390   510   460   420 /* T */
/* AT G G ..GT */
   210   330   280   240 /* A */
   210   330   280   240 /* C */
   210   330   280   240 /* G */
   210   330   280   240 /* T */
/* AT G T ..GT */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT T A ..GT */
   100   220   170   130 /* A */
   100   220   170   130 /* C */
   100   220   170   130 /* G */
   100   220   170   130 /* T */
/* AT T C ..GT */
   390   510   460   420 /* A */
   390   510   460   420 /* C */
   390   510   460   420 /* G */
   390   510   460   420 /* T */
/* AT T G ..GT */
   210   330   280   240 /* A */
   210   330   280   240 /* C */
   210   330   280   240 /* G */
   210   330   280   240 /* T */
/* AT T T ..GT */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT A A ..TG */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT A C ..TG */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* AT A G ..TG */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* AT A T ..TG */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* AT C A ..TG */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT C C ..TG */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* AT C G ..TG */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* AT C T ..TG */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* AT G A ..TG */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT G C ..TG */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* AT G G ..TG */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* AT G T ..TG */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* AT T A ..TG */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT T C ..TG */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* AT T G ..TG */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* AT T T ..TG */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* AT A A ..AT */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* AT A C ..AT */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* AT A G ..AT */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* AT A T ..AT */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* AT C A ..AT */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* AT C C ..AT */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* AT C G ..AT */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* AT C T ..AT */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* AT G A ..AT */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* AT G C ..AT */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* AT G G ..AT */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* AT G T ..AT */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* AT T A ..AT */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* AT T C ..AT */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* AT T G ..AT */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* AT T T ..AT */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* AT A A ..TA */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT A C ..TA */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* AT A G ..TA */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* AT A T ..TA */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* AT C A ..TA */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT C C ..TA */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* AT C G ..TA */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* AT C T ..TA */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* AT G A ..TA */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT G C ..TA */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* AT G G ..TA */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* AT G T ..TA */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* AT T A ..TA */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* AT T C ..TA */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* AT T G ..TA */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* AT T T ..TA */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* AT A A ..NN */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* AT A C ..NN */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* AT A G ..NN */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* AT A T ..NN */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* AT C A ..NN */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* AT C C ..NN */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* AT C G ..NN */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* AT C T ..NN */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* AT G A ..NN */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* AT G C ..NN */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* AT G G ..NN */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* AT G T ..NN */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* AT T A ..NN */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* AT T C ..NN */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* AT T G ..NN */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* AT T T ..NN */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* TA A A ..CG */
  -230  -200  -390 -1040 /* A */
  -230  -200  -390 -1040 /* C */
  -230  -200  -390 -1040 /* G */
  -230  -200  -390 -1040 /* T */
/* TA A C ..CG */
  -190  -160  -350 -1000 /* A */
  -190  -160  -350 -1000 /* C */
  -190  -160  -350 -1000 /* G */
  -190  -160  -350 -1000 /* T */
/* TA A G ..CG */
  -240  -210  -400 -1050 /* A */
  -240  -210  -400 -1050 /* C */
  -240  -210  -400 -1050 /* G */
  -240  -210  -400 -1050 /* T */
/* TA A T ..CG */
  -270  -240  -430 -1080 /* A */
  -270  -240  -430 -1080 /* C */
  -270  -240  -430 -1080 /* G */
  -270  -240  -430 -1080 /* T */
/* TA C A ..CG */
  -230  -200  -390 -1040 /* A */
  -230  -200  -390 -1040 /* C */
  -230  -200  -390 -1040 /* G */
  -230  -200  -390 -1040 /* T */
/* TA C C ..CG */
  -190  -160  -350 -1000 /* A */
  -190  -160  -350 -1000 /* C */
  -190  -160  -350 -1000 /* G */
  -190  -160  -350 -1000 /* T */
/* TA C G ..CG */
  -240  -210  -400 -1050 /* A */
  -240  -210  -400 -1050 /* C */
  -240  -210  -400 -1050 /* G */
  -240  -210  -400 -1050 /* T */
/* TA C T ..CG */
  -270  -240  -430 -1080 /* A */
  -270  -240  -430 -1080 /* C */
  -270  -240  -430 -1080 /* G */
  -270  -240  -430 -1080 /* T */
/* TA G A ..CG */
  -230  -200  -390 -1040 /* A */
  -230  -200  -390 -1040 /* C */
  -230  -200  -390 -1040 /* G */
  -230  -200  -390 -1040 /* T */
/* TA G C ..CG */
  -190  -160  -350 -1000 /* A */
  -190  -160  -350 -1000 /* C */
  -190  -160  -350 -1000 /* G */
  -190  -160  -350 -1000 /* T */
/* TA G G ..CG */
  -240  -210  -400 -1050 /* A */
  -240  -210  -400 -1050 /* C */
  -240  -210  -400 -1050 /* G */
  -240  -210  -400 -1050 /* T */
/* TA G T ..CG */
  -270  -240  -430 -1080 /* A */
  -270  -240  -430 -1080 /* C */
  -270  -240  -430 -1080 /* G */
  -270  -240  -430 -1080 /* T */
/* TA T A ..CG */
  -230  -200  -390 -1040 /* A */
  -230  -200  -390 -1040 /* C */
  -230  -200  -390 -1040 /* G */
  -230  -200  -390 -1040 /* T */
/* TA T C ..CG */
  -190  -160  -350 -1000 /* A */
  -190  -160  -350 -1000 /* C */
  -190  -160  -350 -1000 /* G */
  -190  -160  -350 -1000 /* T */
/* TA T G ..CG */
  -240  -210  -400 -1050 /* A */
  -240  -210  -400 -1050 /* C */
  -240  -210  -400 -1050 /* G */
  -240  -210  -400 -1050 /* T */
/* TA T T ..CG */
  -270  -240  -430 -1080 /* A */
  -270  -240  -430 -1080 /* C */
  -270  -240  -430 -1080 /* G */
  -270  -240  -430 -1080 /* T */
/* TA A A ..GC */
  -390  -360  -550 -1200 /* A */
  -390  -360  -550 -1200 /* C */
  -390  -360  -550 -1200 /* G */
  -390  -360  -550 -1200 /* T */
/* TA A C ..GC */
  -100   -70  -260  -910 /* A */
  -100   -70  -260  -910 /* C */
  -100   -70  -260  -910 /* G */
  -100   -70  -260  -910 /* T */
/* TA A G ..GC */
  -280  -250  -440 -1090 /* A */
  -280  -250  -440 -1090 /* C */
  -280  -250  -440 -1090 /* G */
  -280  -250  -440 -1090 /* T */
/* TA A T ..GC */
   140   170   -20  -670 /* A */
   140   170   -20  -670 /* C */
   140   170   -20  -670 /* G */
   140   170   -20  -670 /* T */
/* TA C A ..GC */
  -390  -360  -550 -1200 /* A */
  -390  -360  -550 -1200 /* C */
  -390  -360  -550 -1200 /* G */
  -390  -360  -550 -1200 /* T */
/* TA C C ..GC */
  -100   -70  -260  -910 /* A */
  -100   -70  -260  -910 /* C */
  -100   -70  -260  -910 /* G */
  -100   -70  -260  -910 /* T */
/* TA C G ..GC */
  -280  -250  -440 -1090 /* A */
  -280  -250  -440 -1090 /* C */
  -280  -250  -440 -1090 /* G */
  -280  -250  -440 -1090 /* T */
/* TA C T ..GC */
   140   170   -20  -670 /* A */
   140   170   -20  -670 /* C */
   140   170   -20  -670 /* G */
   140   170   -20  -670 /* T */
/* TA G A ..GC */
  -390  -360  -550 -1200 /* A */
  -390  -360  -550 -1200 /* C */
  -390  -360  -550 -1200 /* G */
  -390  -360  -550 -1200 /* T */
/* TA G C ..GC */
  -100   -70  -260  -910 /* A */
  -100   -70  -260  -910 /* C */
  -100   -70  -260  -910 /* G */
  -100   -70  -260  -910 /* T */
/* TA G G ..GC */
  -280  -250  -440 -1090 /* A */
  -280  -250  -440 -1090 /* C */
  -280  -250  -440 -1090 /* G */
  -280  -250  -440 -1090 /* T */
/* TA G T ..GC */
   140   170   -20  -670 /* A */
   140   170   -20  -670 /* C */
   140   170   -20  -670 /* G */
   140   170   -20  -670 /* T */
/* TA T A ..GC */
  -390  -360  -550 -1200 /* A */
  -390  -360  -550 -1200 /* C */
  -390  -360  -550 -1200 /* G */
  -390  -360  -550 -1200 /* T */
/* TA T C ..GC */
  -100   -70  -260  -910 /* A */
  -100   -70  -260  -910 /* C */
  -100   -70  -260  -910 /* G */
  -100   -70  -260  -910 /* T */
/* TA T G ..GC */
  -280  -250  -440 -1090 /* A */
  -280  -250  -440 -1090 /* C */
  -280  -250  -440 -1090 /* G */
  -280  -250  -440 -1090 /* T */
/* TA T T ..GC */
   140   170   -20  -670 /* A */
   140   170   -20  -670 /* C */
   140   170   -20  -670 /* G */
   140   170   -20  -670 /* T */
/* TA A A ..GT */
   -70   -40  -230  -880 /* A */
   -70   -40  -230  -880 /* C */
   -70   -40  -230  -880 /* G */
   -70   -40  -230  -880 /* T */
/* TA A C ..GT */
   220   250    60  -590 /* A */
   220   250    60  -590 /* C */
   220   250    60  -590 /* G */
   220   250    60  -590 /* T */
/* TA A G ..GT */
    40    70  -120  -770 /* A */
    40    70  -120  -770 /* C */
    40    70  -120  -770 /* G */
    40    70  -120  -770 /* T */
/* TA A T ..GT */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA C A ..GT */
   -70   -40  -230  -880 /* A */
   -70   -40  -230  -880 /* C */
   -70   -40  -230  -880 /* G */
   -70   -40  -230  -880 /* T */
/* TA C C ..GT */
   220   250    60  -590 /* A */
   220   250    60  -590 /* C */
   220   250    60  -590 /* G */
   220   250    60  -590 /* T */
/* TA C G ..GT */
    40    70  -120  -770 /* A */
    40    70  -120  -770 /* C */
    40    70  -120  -770 /* G */
    40    70  -120  -770 /* T */
/* TA C T ..GT */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA G A ..GT */
   -70   -40  -230  -880 /* A */
   -70   -40  -230  -880 /* C */
   -70   -40  -230  -880 /* G */
   -70   -40  -230  -880 /* T */
/* TA G C ..GT */
   220   250    60  -590 /* A */
   220   250    60  -590 /* C */
   220   250    60  -590 /* G */
   220   250    60  -590 /* T */
/* TA G G ..GT */
    40    70  -120  -770 /* A */
    40    70  -120  -770 /* C */
    40    70  -120  -770 /* G */
    40    70  -120  -770 /* T */
/* TA G T ..GT */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA T A ..GT */
   -70   -40  -230  -880 /* A */
   -70   -40  -230  -880 /* C */
   -70   -40  -230  -880 /* G */
   -70   -40  -230  -880 /* T */
/* TA T C ..GT */
   220   250    60  -590 /* A */
   220   250    60  -590 /* C */
   220   250    60  -590 /* G */
   220   250    60  -590 /* T */
/* TA T G ..GT */
    40    70  -120  -770 /* A */
    40    70  -120  -770 /* C */
    40    70  -120  -770 /* G */
    40    70  -120  -770 /* T */
/* TA T T ..GT */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA A A ..TG */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA A C ..TG */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TA A G ..TG */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TA A T ..TG */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TA C A ..TG */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA C C ..TG */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TA C G ..TG */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TA C T ..TG */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TA G A ..TG */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA G C ..TG */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TA G G ..TG */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TA G T ..TG */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TA T A ..TG */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA T C ..TG */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TA T G ..TG */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TA T T ..TG */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TA A A ..AT */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TA A C ..AT */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TA A G ..AT */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TA A T ..AT */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TA C A ..AT */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TA C C ..AT */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TA C G ..AT */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TA C T ..AT */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TA G A ..AT */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TA G C ..AT */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TA G G ..AT */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TA G T ..AT */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TA T A ..AT */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TA T C ..AT */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TA T G ..AT */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TA T T ..AT */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TA A A ..TA */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA A C ..TA */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TA A G ..TA */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TA A T ..TA */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TA C A ..TA */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA C C ..TA */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TA C G ..TA */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TA C T ..TA */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TA G A ..TA */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA G C ..TA */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TA G G ..TA */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TA G T ..TA */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TA T A ..TA */
   460   490   300  -350 /* A */
   460   490   300  -350 /* C */
   460   490   300  -350 /* G */
   460   490   300  -350 /* T */
/* TA T C ..TA */
   490   520   330  -320 /* A */
   490   520   330  -320 /* C */
   490   520   330  -320 /* G */
   490   520   330  -320 /* T */
/* TA T G ..TA */
   300   330   140  -510 /* A */
   300   330   140  -510 /* C */
   300   330   140  -510 /* G */
   300   330   140  -510 /* T */
/* TA T T ..TA */
  -350  -320  -510 -1160 /* A */
  -350  -320  -510 -1160 /* C */
  -350  -320  -510 -1160 /* G */
  -350  -320  -510 -1160 /* T */
/* TA A A ..NN */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TA A C ..NN */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TA A G ..NN */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TA A T ..NN */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TA C A ..NN */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TA C C ..NN */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TA C G ..NN */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TA C T ..NN */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TA G A ..NN */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TA G C ..NN */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TA G G ..NN */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TA G T ..NN */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* TA T A ..NN */
   630   660   470  -180 /* A */
   630   660   470  -180 /* C */
   630   660   470  -180 /* G */
   630   660   470  -180 /* T */
/* TA T C ..NN */
   750   780   590   -60 /* A */
   750   780   590   -60 /* C */
   750   780   590   -60 /* G */
   750   780   590   -60 /* T */
/* TA T G ..NN */
   700   730   540  -110 /* A */
   700   730   540  -110 /* C */
   700   730   540  -110 /* G */
   700   730   540  -110 /* T */
/* TA T T ..NN */
   660   690   500  -150 /* A */
   660   690   500  -150 /* C */
   660   690   500  -150 /* G */
   660   690   500  -150 /* T */
/* NN A A ..CG */
   -60    60    10   -30 /* A */
   -60    60    10   -30 /* C */
   -60    60    10   -30 /* G */
   -60    60    10   -30 /* T */
/* NN A C ..CG */
   -20   100    50    10 /* A */
   -20   100    50    10 /* C */
   -20   100    50    10 /* G */
   -20   100    50    10 /* T */
/* NN A G ..CG */
   -70    50     0   -40 /* A */
   -70    50     0   -40 /* C */
   -70    50     0   -40 /* G */
   -70    50     0   -40 /* T */
/* NN A T ..CG */
  -100    20   -30   -70 /* A */
  -100    20   -30   -70 /* C */
  -100    20   -30   -70 /* G */
  -100    20   -30   -70 /* T */
/* NN C A ..CG */
   -60    60    10   -30 /* A */
   -60    60    10   -30 /* C */
   -60    60    10   -30 /* G */
   -60    60    10   -30 /* T */
/* NN C C ..CG */
   -20   100    50    10 /* A */
   -20   100    50    10 /* C */
   -20   100    50    10 /* G */
   -20   100    50    10 /* T */
/* NN C G ..CG */
   -70    50     0   -40 /* A */
   -70    50     0   -40 /* C */
   -70    50     0   -40 /* G */
   -70    50     0   -40 /* T */
/* NN C T ..CG */
  -100    20   -30   -70 /* A */
  -100    20   -30   -70 /* C */
  -100    20   -30   -70 /* G */
  -100    20   -30   -70 /* T */
/* NN G A ..CG */
   -60    60    10   -30 /* A */
   -60    60    10   -30 /* C */
   -60    60    10   -30 /* G */
   -60    60    10   -30 /* T */
/* NN G C ..CG */
   -20   100    50    10 /* A */
   -20   100    50    10 /* C */
   -20   100    50    10 /* G */
   -20   100    50    10 /* T */
/* NN G G ..CG */
   -70    50     0   -40 /* A */
   -70    50     0   -40 /* C */
   -70    50     0   -40 /* G */
   -70    50     0   -40 /* T */
/* NN G T ..CG */
  -100    20   -30   -70 /* A */
  -100    20   -30   -70 /* C */
  -100    20   -30   -70 /* G */
  -100    20   -30   -70 /* T */
/* NN T A ..CG */
   -60    60    10   -30 /* A */
   -60    60    10   -30 /* C */
   -60    60    10   -30 /* G */
   -60    60    10   -30 /* T */
/* NN T C ..CG */
   -20   100    50    10 /* A */
   -20   100    50    10 /* C */
   -20   100    50    10 /* G */
   -20   100    50    10 /* T */
/* NN T G ..CG */
   -70    50     0   -40 /* A */
   -70    50     0   -40 /* C */
   -70    50     0   -40 /* G */
   -70    50     0   -40 /* T */
/* NN T T ..CG */
  -100    20   -30   -70 /* A */
  -100    20   -30   -70 /* C */
  -100    20   -30   -70 /* G */
  -100    20   -30   -70 /* T */
/* NN A A ..GC */
  -220  -100  -150  -190 /* A */
  -220  -100  -150  -190 /* C */
  -220  -100  -150  -190 /* G */
  -220  -100  -150  -190 /* T */
/* NN A C ..GC */
    70   190   140   100 /* A */
    70   190   140   100 /* C */
    70   190   140   100 /* G */
    70   190   140   100 /* T */
/* NN A G ..GC */
  -110    10   -40   -80 /* A */
  -110    10   -40   -80 /* C */
  -110    10   -40   -80 /* G */
  -110    10   -40   -80 /* T */
/* NN A T ..GC */
   310   430   380   340 /* A */
   310   430   380   340 /* C */
   310   430   380   340 /* G */
   310   430   380   340 /* T */
/* NN C A ..GC */
  -220  -100  -150  -190 /* A */
  -220  -100  -150  -190 /* C */
  -220  -100  -150  -190 /* G */
  -220  -100  -150  -190 /* T */
/* NN C C ..GC */
    70   190   140   100 /* A */
    70   190   140   100 /* C */
    70   190   140   100 /* G */
    70   190   140   100 /* T */
/* NN C G ..GC */
  -110    10   -40   -80 /* A */
  -110    10   -40   -80 /* C */
  -110    10   -40   -80 /* G */
  -110    10   -40   -80 /* T */
/* NN C T ..GC */
   310   430   380   340 /* A */
   310   430   380   340 /* C */
   310   430   380   340 /* G */
   310   430   380   340 /* T */
/* NN G A ..GC */
  -220  -100  -150  -190 /* A */
  -220  -100  -150  -190 /* C */
  -220  -100  -150  -190 /* G */
  -220  -100  -150  -190 /* T */
/* NN G C ..GC */
    70   190   140   100 /* A */
    70   190   140   100 /* C */
    70   190   140   100 /* G */
    70   190   140   100 /* T */
/* NN G G ..GC */
  -110    10   -40   -80 /* A */
  -110    10   -40   -80 /* C */
  -110    10   -40   -80 /* G */
  -110    10   -40   -80 /* T */
/* NN G T ..GC */
   310   430   380   340 /* A */
   310   430   380   340 /* C */
   310   430   380   340 /* G */
   310   430   380   340 /* T */
/* NN T A ..GC */
  -220  -100  -150  -190 /* A */
  -220  -100  -150  -190 /* C */
  -220  -100  -150  -190 /* G */
  -220  -100  -150  -190 /* T */
/* NN T C ..GC */
    70   190   140   100 /* A */
    70   190   140   100 /* C */
    70   190   140   100 /* G */
    70   190   140   100 /* T */
/* NN T G ..GC */
  -110    10   -40   -80 /* A */
  -110    10   -40   -80 /* C */
  -110    10   -40   -80 /* G */
  -110    10   -40   -80 /* T */
/* NN T T ..GC */
   310   430   380   340 /* A */
   310   430   380   340 /* C */
   310   430   380   340 /* G */
   310   430   380   340 /* T */
/* NN A A ..GT */
   100   220   170   130 /* A */
   100   220   170   130 /* C */
   100   220   170   130 /* G */
   100   220   170   130 /* T */
/* NN A C ..GT */
   390   510   460   420 /* A */
   390   510   460   420 /* C */
   390   510   460   420 /* G */
   390   510   460   420 /* T */
/* NN A G ..GT */
   210   330   280   240 /* A */
   210   330   280   240 /* C */
   210   330   280   240 /* G */
   210   330   280   240 /* T */
/* NN A T ..GT */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN C A ..GT */
   100   220   170   130 /* A */
   100   220   170   130 /* C */
   100   220   170   130 /* G */
   100   220   170   130 /* T */
/* NN C C ..GT */
   390   510   460   420 /* A */
   390   510   460   420 /* C */
   390   510   460   420 /* G */
   390   510   460   420 /* T */
/* NN C G ..GT */
   210   330   280   240 /* A */
   210   330   280   240 /* C */
   210   330   280   240 /* G */
   210   330   280   240 /* T */
/* NN C T ..GT */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN G A ..GT */
   100   220   170   130 /* A */
   100   220   170   130 /* C */
   100   220   170   130 /* G */
   100   220   170   130 /* T */
/* NN G C ..GT */
   390   510   460   420 /* A */
   390   510   460   420 /* C */
   390   510   460   420 /* G */
   390   510   460   420 /* T */
/* NN G G ..GT */
   210   330   280   240 /* A */
   210   330   280   240 /* C */
   210   330   280   240 /* G */
   210   330   280   240 /* T */
/* NN G T ..GT */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN T A ..GT */
   100   220   170   130 /* A */
   100   220   170   130 /* C */
   100   220   170   130 /* G */
   100   220   170   130 /* T */
/* NN T C ..GT */
   390   510   460   420 /* A */
   390   510   460   420 /* C */
   390   510   460   420 /* G */
   390   510   460   420 /* T */
/* NN T G ..GT */
   210   330   280   240 /* A */
   210   330   280   240 /* C */
   210   330   280   240 /* G */
   210   330   280   240 /* T */
/* NN T T ..GT */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN A A ..TG */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN A C ..TG */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* NN A G ..TG */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* NN A T ..TG */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* NN C A ..TG */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN C C ..TG */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* NN C G ..TG */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* NN C T ..TG */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* NN G A ..TG */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN G C ..TG */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* NN G G ..TG */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* NN G T ..TG */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* NN T A ..TG */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN T C ..TG */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* NN T G ..TG */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* NN T T ..TG */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* NN A A ..AT */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* NN A C ..AT */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* NN A G ..AT */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* NN A T ..AT */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* NN C A ..AT */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* NN C C ..AT */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* NN C G ..AT */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* NN C T ..AT */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* NN G A ..AT */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* NN G C ..AT */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* NN G G ..AT */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* NN G T ..AT */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* NN T A ..AT */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* NN T C ..AT */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* NN T G ..AT */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* NN T T ..AT */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* NN A A ..TA */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN A C ..TA */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* NN A G ..TA */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* NN A T ..TA */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* NN C A ..TA */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN C C ..TA */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* NN C G ..TA */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* NN C T ..TA */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* NN G A ..TA */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN G C ..TA */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* NN G G ..TA */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* NN G T ..TA */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* NN T A ..TA */
   630   750   700   660 /* A */
   630   750   700   660 /* C */
   630   750   700   660 /* G */
   630   750   700   660 /* T */
/* NN T C ..TA */
   660   780   730   690 /* A */
   660   780   730   690 /* C */
   660   780   730   690 /* G */
   660   780   730   690 /* T */
/* NN T G ..TA */
   470   590   540   500 /* A */
   470   590   540   500 /* C */
   470   590   540   500 /* G */
   470   590   540   500 /* T */
/* NN T T ..TA */
  -180   -60  -110  -150 /* A */
  -180   -60  -110  -150 /* C */
  -180   -60  -110  -150 /* G */
  -180   -60  -110  -150 /* T */
/* NN A A ..NN */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* NN A C ..NN */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* NN A G ..NN */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* NN A T ..NN */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* NN C A ..NN */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* NN C C ..NN */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* NN C G ..NN */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* NN C T ..NN */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* NN G A ..NN */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* NN G C ..NN */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* NN G G ..NN */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* NN G T ..NN */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */
/* NN T A ..NN */
   800   920   870   830 /* A */
   800   920   870   830 /* C */
   800   920   870   830 /* G */
   800   920   870   830 /* T */
/* NN T C ..NN */
   920  1040   990   950 /* A */
   920  1040   990   950 /* C */
   920  1040   990   950 /* G */
   920  1040   990   950 /* T */
/* NN T G ..NN */
   870   990   940   900 /* A */
   870   990   940   900 /* C */
   870   990   940   900 /* G */
   870   990   940   900 /* T */
/* NN T T ..NN */
   830   950   900   860 /* A */
   830   950   900   860 /* C */
   830   950   900   860 /* G */
   830   950   900   860 /* T */

# hairpin
   INF   INF   INF   340   340   350   420   420   420   430
   440   450   460   470   480   500   510   520   530   540
   550   560   570   580   590   600   610   620   630   640
   650

# bulge
   INF   380   280   320   350   380   400   420   430   440
   450   460   470   480   490   500   500   510   520   520
   530   530   540   540   550   550   560   560   560   570
   570

# internal_loop
   INF   INF   INF   INF   250   270   290   310   320   340
   350   360   370   380   390   390   400   410   410   420
   420   430   430   440   440   450   450   460   460   460
   470

# ML_params
/* the numbers below are really meant to be used with co-axial stacking */
/* without the penalty is certainly too large */
/* F = cu*n_unpaired + cc + ci*loop_degree (+TermAU) */
/*      cu      cc      ci   TerminalAU DuplexInit*/
         0    1280     -60      20        210

# NINIO
/* Ninio = MIN(max, m*|n1-n2| */
/*       m   max              */
     60  300

# Tetraloops
    GGUCAC   -80  /* 160 */
    GGCUAC   -70  /* 170 */
    GGAAAC     0  /* 240 */
    CUUUUG   -60  /* 220 */
    CAAAAG     0  /* 280 */
    AUUUUU  -380  /*  10 */
    AAAAAU   -70  /* 280 */
    UUAAUA  -120  /* 180 */
    GAAAAC   -60  /* 180 */
    GUAAUC   -80  /* 180 */
    CUUUUG  -200  /*  80 */
        GGTCAC   -80  /* 160 */
        GGCTAC   -70  /* 170 */
        GGAAAC     0  /* 240 */
        CTTTTG   -60  /* 220 */
        CAAAAG     0  /* 280 */
        ATTTTT  -380  /*  10 */
        AAAAAT   -70  /* 280 */
        TTAATA  -120  /* 180 */
        GAAAAC   -60  /* 180 */
        GTAATC   -80  /* 180 */
        CTTTTG  -200  /*  80 */
        CAAAAG   -90  /* 190 */
    CAAAAG   -90  /* 190 */

# Triloops
    CGCAG   -390  /* -50 */
    CUUUG   -290  /*  50 */
    AUAAU     80  /* 460 */
        CGCAG   -390  /* -50 */
        CTTTG   -290  /*  50 */
        ATAAT     80  /* 460 */

# mismatch_hairpin
/* CG.. */
   -40   -60   -40   -60   -60
   -40   -60   -40   -60   -60
   -40   -60   -40   -60   -60
   -40   -60   -40   -60   -60
   -40   -60   -40   -60   -60
/* GC.. */
   -50  -100   -50   -60   -80
   -50  -100   -50   -60   -80
   -50  -100   -50   -60   -80
   -50  -100   -50   -60   -80
   -50  -100   -50   -60   -80
/* GT.. */
   -30   -80   -30   -40   -60
   -30   -80   -30   -40   -60
   -30   -80   -30   -40   -60
   -30   -80   -30   -40   -60
   -30   -80   -30   -40   -60
/* TG.. */
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
/* AT.. */
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30
/* TA.. */
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
/* NN.. */
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30

# mismatch_interior
/* CG.. */
   -40   -60   -40   -60   -60
   -40   -60   -40   -60   -60
   -40   -60   -40   -60   -60
   -40   -60   -40   -60   -60
   -40   -60   -40   -60   -60
/* GC.. */
   -50  -100   -50   -60   -80
   -50  -100   -50   -60   -80
   -50  -100   -50   -60   -80
   -50  -100   -50   -60   -80
   -50  -100   -50   -60   -80
/* GT.. */
   -30   -80   -30   -40   -60
   -30   -80   -30   -40   -60
   -30   -80   -30   -40   -60
   -30   -80   -30   -40   -60
   -30   -80   -30   -40   -60
/* TG.. */
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
/* AT.. */
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30
/* TA.. */
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
   -20   -20   -20   -40   -60
/* NN.. */
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30
    90   -10    30    90    30

# mismatch_enthalpies
/* CG.. */
  -420  -460  -420  -470  -500
  -420  -460  -420  -470  -500
  -420  -460  -420  -470  -500
  -420  -460  -420  -470  -500
  -420  -460  -420  -470  -500
/* GC.. */
   -90  -620  -330  -510   -90
   -90  -620  -330  -510   -90
   -90  -620  -330  -510   -90
   -90  -620  -330  -510   -90
   -90  -620  -330  -510   -90
/* GT.. */
   230  -300   -10  -190   230
   230  -300   -10  -190   230
   230  -300   -10  -190   230
   230  -300   -10  -190   230
   230  -300   -10  -190   230
/* TG.. */
   230   230   170    70  -580
   230   230   170    70  -580
   230   230   170    70  -580
   230   230   170    70  -580
   230   230   170    70  -580
/* AT.. */
   520   400   520   470   430
   520   400   520   470   430
   520   400   520   470   430
   520   400   520   470   430
   520   400   520   470   430
/* TA.. */
   230   230   170    70  -580
   230   230   170    70  -580
   230   230   170    70  -580
   230   230   170    70  -580
   230   230   170    70  -580
/* NN.. */
   520   400   520   470   430
   520   400   520   470   430
   520   400   520   470   430
   520   400   520   470   430
   520   400   520   470   430"""
        
if __name__ == "__main__":
    main()
