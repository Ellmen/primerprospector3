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

""" get_linker is designed to suggest primer linkers based upon the least
frequently occuring base pairs (default 2) immediately upstream from the 
5' region of the primer.

A single _hits.txt file, or a directory of _hits.txt files can be passed as
input (see the archaeal_primers_util_modified.py file for information about
the hits.txt file format).  The fasta files used to generate the _hits.txt
files are also required input for this function.

The output file will contine the following information for each _hits.txt file
passed:  The percentage occurance of each base at all positions in the linker,
and a suggested linker based on complementary bases to the least frequently
occuring bases.
"""

from os.path import basename

from cogent3 import DNA

from primerprospector.parse import get_primer_seq, get_hits_field,\
 get_primer_hits_data, parse_hits_data
from primerprospector.util import get_primer_direction
from primerprospector.format import format_base_freq_data

    
def generate_primers_pos_data(hits_fp,
                              output_dir = ".",
                              score_type = "weighted_score",
                              score_threshold = 2.0):
    """ Main program function for listing primer base position data
    hits_fp: filepath for primer hits to be tested.
    output_dir: output directory
    score_type: The type of scoring used to determine if a sequence should
     be considered for base position frequency.  Can be weighted_score, 
     tp_mismatches, or overall_mismatches
    score_threshold: Value at which or below a sequence score is considered
     acceptable to be considered (i.e., that sequence is
     likely to amplify during PCR).
    
    """
    
    # Get primer sequence from hits file
    primer_seq = get_primer_seq(hits_fp)
    
    # Get primer direction, based on hits file name
    primer_direction = get_primer_direction(hits_fp)
    
    # Get hits lines, parse out lines that fail to meet score threshold
    raw_hits_data = get_primer_hits_data(hits_fp)
    hits_data = parse_hits_data(raw_hits_data, score_threshold, score_type)
    
    # Pull out sequence hits fields from the hits_data
    sequence_hit_index = 1
    sequence_hits = get_hits_field(hits_data, sequence_hit_index)
    
    # Get reverse complement of sequence hits if direction is reverse
    if primer_direction == "r":
        sequence_hits = reverse_sequence_hits(sequence_hits)
        
    # Filter out sequence hits that are not the same length as the primer
    # Should not happen with the current local alignment module, but just to 
    # be sure, testing here
    sequence_hits = filter_uneven_sequence_hits(sequence_hits, primer_seq)
    
    base_freqs = get_base_freqs(sequence_hits, primer_seq)
    
    freq_file_data = format_base_freq_data(base_freqs, primer_seq, hits_fp,
     score_type, score_threshold)
     
    base_freqs_fp = get_primer_base_freqs_fp(output_dir, hits_fp)
    
    base_freqs_f = open(base_freqs_fp, "w")
    
    base_freqs_f.write(freq_file_data)
     
    
    
        
    
def filter_uneven_sequence_hits(sequence_hits,
                                primer_seq):
    """ Filters out sequence hits that do not match primer length
    
    sequence_hits: list of sequence hits.
    primer_seq: sequence of primer from hits file.
    """
    
    filtered_seqs = []
    
    primer_len = len(primer_seq)
    
    for seq in sequence_hits:
        if len(seq) == primer_len:
            filtered_seqs.append(seq)
            
    return filtered_seqs
    
def reverse_sequence_hits(sequence_hits):
    """ Returns reverse complement of sequences in given list
    
    sequence_hits: list of sequences in IUPAC format
    """
    
    reverse_seqs = []
    
    for seq in sequence_hits:
        reverse_seqs.append(DNA.rc(seq))
        
    return reverse_seqs
    
def get_base_freqs(sequence_hits,
                   primer):
    """ Generates list ofdict of percent base occurace for each base in linker
    
    sequence_hits:  List of sequence hits from each input sequence
     that passed scoring threshold.
    primer: primer sequence from hits file
    """
    
    # Get total number of sequences considered for percent calcuations
    total_seqs = len(sequence_hits)
    
    primer_len = len(primer)

    # Will ignore degenerate bases, gaps for purposes of getting base freqs
    nucleotides = ['A', 'T', 'C', 'G']
    
    base_counts = []
    base_freqs = []

    for base_pos in range(primer_len):
        base_counts.append({'A':0, 'T':0, 'C':0, 'G':0})
        base_freqs.append({'A':0, 'T':0, 'C':0, 'G':0})
        for seq in sequence_hits:
            
            # Add counts of each base for given base index, skip degeneracies
            if seq[base_pos] not in nucleotides:
                continue
            
            base_counts[base_pos][seq[base_pos]] += 1
                

        # Get percentages for each base
        for base in base_counts[base_pos].keys():
            base_freqs[base_pos][base] =\
             base_counts[base_pos][base] / total_seqs

    return base_freqs
    
    
    
    
def get_primer_base_freqs_fp(output_dir,
                             hits):
    """ Generates base frequency output filepath
    
    hits: filepath of input hits file, used to generate output filename
    output_dir: output directory """
    
    if not output_dir.endswith('/'):
        output_dir += '/'
        
    base_freqs_fp = output_dir + basename(hits).split('.')[0] + \
     "_base_frequencies.txt"
     
    return base_freqs_fp
    
    
