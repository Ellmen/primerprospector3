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
 

from cogent.util.unit_test import TestCase, main

from primerprospector.util import get_DNA_to_numeric, get_numeric_to_DNA,\
 correct_primer_name, get_primer_direction

class UtilTests(TestCase):
    """ Unit tests for the util.py module of primer prospector"""
    
    def setUp(self):
        pass
        
    def test_sequence_sanity(self):
        """ NT should be converted back to same value from numeric """
        
        nucleotides = "ATCGYRNSDHVBW-"
        
        for n in nucleotides:
            curr_nt = n
            test_nt = get_numeric_to_DNA()[get_DNA_to_numeric()[n]]
            self.assertEqual(test_nt, curr_nt)
            
    def test_correct_primer_name(self):
        """ Should autocorrect caps errors for f and r on primer names """
        
        # Should not alter primer that has lowercase suffix
        expected_primer = '238f'
        
        actual_primer = correct_primer_name('238f')
        
        self.assertEqual(actual_primer, expected_primer)
        
        # Should correct if caps at end of primer
        
        actual_primer = correct_primer_name('238F')
        
        self.assertEqual(actual_primer, expected_primer)
        
    def test_get_primer_direction(self):
        """ Returns primer direction from primer hits filename correctly """
        
        # Test for forward hits file in current directory
        hits_filepath = "23f_hits.txt"
        expected_direction = "f"
        
        actual_direction = get_primer_direction(hits_filepath)
        
        self.assertEqual(actual_direction, expected_direction)
        
        # Test for relative directories
        hits_filepath = "../323r_hits.txt"
        expected_direction = "r"
        
        actual_direction = get_primer_direction(hits_filepath)
        
        self.assertEqual(actual_direction, expected_direction)
        
        # Should raise error if incorrect hits file naming format
        
        bad_hits_filepath = "F422_hits.txt"
        
        self.assertRaises(ValueError, get_primer_direction, bad_hits_filepath)
            

    

if __name__ == "__main__":
    main()
