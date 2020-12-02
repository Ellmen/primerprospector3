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
 

lower = str.lower
from os.path import basename

""" Contains dictionaries and data for usage in primer design and analysis,
 as well as other miscellaneous functions. """

class FakeOutFile(object):
    """ Class to test file output """
    
    def __init__(self):
        self.data = ""
    
    def write(self,s):
        self.data += s

def get_DNA_to_numeric():
    return {'A':1, 'T':2, 'U':2, 'C':4, 'G':8, 'N':15, 'R':9,\
    'Y':6, 'M':5, 'K':10, 'W':3, 'S':12, 'B':14, 'D':11, 'H':7,\
    'V':13, '-':0}
        
def get_numeric_to_DNA():
    return {1:'A', 2:'T', 4:'C', 8:'G', 15:'N', 9:'R', 6:'Y',
    5:'M', 10:'K', 3:'W', 12:'S', 14:'B', 11:'D', 7:'H', 13:'V', 0:'-'}
    
def correct_primer_name(primer_name):
    "Tries to autocorrect primer name to have lowercase f or r at end """
    
    # Fix primer name if not in correct format (followed by lower case f or
    # r.  Leave other components unchanged.
    primer_prefix = primer_name[:-1]
    primer_suffix = lower(primer_name[-1])
    primer_name = primer_prefix + primer_suffix
    
    return primer_name
    
def get_primer_direction(hits_file):
    """ Returns F or R direction for primer based on hits file name 
    
    hits_file: Current hits file filepath"""
    
    # Split this way in case multiple directories in filepath.
    hits_name=basename(hits_file).split("_")[0]
    if hits_name.endswith("f"):
        return "f"
    elif hits_name.endswith("r"):
        return "r"
    else:
        raise ValueError('%s not named correctly, all hits files '%hits_name+\
         'must start with the numeric value indicating position of the '+\
         'primer followed by "f" or "r".  Example: 219f_bacteria_hits.txt')
    

    

