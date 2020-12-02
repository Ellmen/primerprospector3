#!/usr/bin/env python
# Author: William Walters
# sort_denovo_primers.py

from __future__ import division

from future.utils import raise_
__author__ = "William Walters"
__copyright__ = "Copyright 2010, The Primer Prospector Project"
__credits__ =  ["William Walters"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

""" Description
File created 22 July 2009.

The purpose of this module is to analyze the conserved sequence hits file
(see generate_primers_denovo.py) to determine if the upstream or downstream
sequences from the perfectly (or near so) matching 3' primer end are 
reasonably conserved and suitable for use in designing primers.

This module uses the sequences given in the hits file to calculate
the shannon entropies for each position in the primers.  The output
from this module will contain the overall consensus sequence for each
potential primer (going upstream or downstream), an degenerate 
IUPAC sequence that considers all bases found in each position, a 
filtered degenerate IUPAC sequence that does not consider bases whose
percentage is under a specified value, and the shannon entropy scores for
the overall sequences.

Note-the non-filtered IUPAC sequence output will contain a "." for positions
that contain a "-" character (this results from filling in gap characters at
unknown bases that exceed the beginning or end of a sequence in the 
generate_primers_denovo.py program).

Following this initial analysis, the primers are then sorted into a summary
file containing information about each prospective primer and a primer file 
that is formatted for use with the analyze_primers.py module.  
The primers can be sorted according to sensitivity (greatest to least),
specificity (most to least), or shannon entropy of the overall primer.

Finally, a list of preexisting primers can be passed via the -k option to 
compare to the de novo primers generated by generate_primers_denovo.py.
The de novo primers are compared to the primers passed, and flags any
primers that overlap (considers matching degenerate characters as well).
The 'primers_overlap.txt' file contains information about the 
overlapping primer for the entire primer set.  This file contains a section
showing primers that have a 'match' to the supplied primers, meaning that the
primers overlap and match at the 3' end.  An 'overlap' section shows details 
about primers that overlap with the given primers but do not match at the 
3' end.  Finally, a 'unique' primers section shows details about primers that
do not overlap with the supplied primer set.

These formatted primer files are in the following format:
primer_id <tab> primer sequence (5'->3')  <tab> primer citation
Any comments are preceeded by the pound (#) symbol. 

If a standard alignment was used to record indices in the 
generate_primers_denovo.py module (-a option), this module will detect the 
presence of the standard aligned indices.  If absent, the primers' numeric 
names will be based on the initial unaligned index of the sequence they were 
found in.
"""

from operator import itemgetter
from os import remove
from copy import deepcopy

from cogent import LoadSeqs, DNA
from cogent.core.alphabet import AlphabetError
from cogent.core.moltype import IUPAC_DNA_ambiguities, IUPAC_DNA_chars,\
 IUPAC_gap, make_pairs, DnaStandardPairs
from numpy import array, bitwise_or, bitwise_and

from primerprospector.util import get_DNA_to_numeric, get_numeric_to_DNA
from primerprospector.parse import parse_primer_data, get_known_primers,\
 get_primers_data
from primerprospector.format import write_formatted_primers, \
 write_overlap_primers, write_primers_summary, organize_overlapping_primers,\
 get_amplicon_pairs, write_amplicon_pairs


class Primer(object):
    """ Stores, determines information about primer sequences """
    def __init__(self, 
                 header, 
                 f_primer_seqs, 
                 r_primer_seqs,
                 variable_pos_freq=0.20, 
                 primer_name = "", 
                 cmp_trunc_len = 5):
        """ initialize information with PyCogent alignment functions 
        
        header: header line from primer hits file, contains information
         regarding primer sensitivity, specificity, index, and standard
         aligned index if these data were generated.
        f_primer_seqs: sequence collection object of the conserved
         Xmer and upstream sequences.
        r_primer_seqs: sequence collection object of the conserved 
         Xmer and downstream sequences.
        variable_pos_freq: frequency at which a particular base has to appear
         to be considered when making the filtered degenerate version of new
         primers.  Lower to make more degenerate, raise to limit degeneracy.
        primer_name: optional string that will be inserted into primer name
         if specified.
        cmp_trunc_len: length of 3' end of known primers (see -k option) to 
         compare to the new primers to determine if primer is duplicate of
         already known primer."""
        
        # f_ and r_ refer to forward and reverse primers
        self.header=header
        self.f_primer_seqs=f_primer_seqs
        self.r_primer_seqs=r_primer_seqs
        self.f_consensus="".join(f_primer_seqs.majorityConsensus())
        # Generate Reverse complement for the reverse primer
        self.r_consensus=DNA.rc("".join(r_primer_seqs.majorityConsensus()))
        self.f_varied=f_primer_seqs.variablePositions()
        self.r_varied=r_primer_seqs.variablePositions()
        self.f_base_freq=f_primer_seqs.columnProbs()
        self.r_base_freq=r_primer_seqs.columnProbs()
        self.f_uncertainties=f_primer_seqs.uncertainties()
        self.r_uncertainties=r_primer_seqs.uncertainties()
        # Generate filtered results
        self.filtered_f_base_freq=\
         self.filter_base_freq(self.f_base_freq,variable_pos_freq)
        self.filtered_r_base_freq=\
         self.filter_base_freq(self.r_base_freq,variable_pos_freq)
        self.f_degenerate_seq=\
         self.get_filtered_degenerate_seq(self.filtered_f_base_freq)
        # Get reverse complement
        self.r_degenerate_seq=DNA.rc(\
         self.get_filtered_degenerate_seq(self.filtered_r_base_freq))
        self.f_overall_score=sum(self.f_uncertainties)
        self.r_overall_score=sum(self.r_uncertainties)
        
        self.sensitivity = float(self.header.split(',')[4].replace("%",""))
        self.specificity = float(self.header.split(',')[5].replace("%",""))
        # Detect if standard alignment file was supplied
        try:
            self.primer_f_name = primer_name + self.header.split(',')[6] + "f"
            self.primer_r_name = primer_name + self.header.split(',')[7] + "r"
        except IndexError:
            self.primer_f_name = primer_name + self.header.split(',')[1] + "f"
            self.primer_r_name = primer_name + self.header.split(',')[1] + "r"

        
        # IUPAC consensus results, generally very degenerate
        self.f_IUPAC=f_primer_seqs.IUPACConsensus().replace("?",".")
        self.r_IUPAC=DNA.rc(r_primer_seqs.IUPACConsensus().replace("?","."))
        
        # normal, complement, reverse, and reverse-complement forms of the 
        # forward and reverse primer for testing with overlap with known primers
        # In integer form for bitwise comparison
        
        self.f_seq = convert_to_numeric(self.f_degenerate_seq)
        self.f_complement =\
         convert_to_numeric(DNA.complement(self.f_degenerate_seq))
        self.f_reverse = convert_to_numeric(self.f_degenerate_seq[::-1])
        self.f_rc = convert_to_numeric(DNA.rc(self.f_degenerate_seq))
        
        self.r_seq = convert_to_numeric(self.r_degenerate_seq)
        self.r_complement =\
         convert_to_numeric(DNA.complement(self.r_degenerate_seq))
        self.r_reverse = convert_to_numeric(self.r_degenerate_seq[::-1])
        self.r_rc = convert_to_numeric(DNA.rc(self.r_degenerate_seq))
        
        # Store sequence of 3' end for comparison with known primers
        self.f_trunc_seq =\
         convert_to_numeric(self.f_degenerate_seq[-cmp_trunc_len:])
        self.r_trunc_seq =\
         convert_to_numeric(self.r_degenerate_seq[-cmp_trunc_len:])        
        # Store information about overlapping primers
        
        self.f_unique = True
        self.f_partial_overlap = []
        self.f_3prime_match = []
        self.r_unique = True
        self.r_partial_overlap = []
        self.r_3prime_match = []

        
    def filter_base_freq(self, base_freq, variable_pos_freq=0.20):
        """ Creates filtered base freq with rare occurances removed 

        Will retain the base(s) with the highest frequency, regardless
        of variable_pos_freq value
        
        base_freq: contains base frequencies, generated from columnProbs()
         method of sequence collection object
        variable_pos_freq: frequency at which a particular base has to appear
         to be considered when making the filtered degenerate version of new
         primers.  Lower to make more degenerate, raise to limit degeneracy."""

        pos_results=deepcopy(base_freq)

        for base_pos in pos_results:
            highest_freq=0
            # Determine highest frequency for retention
            for base in base_pos.keys():
                if base_pos[base]>highest_freq:
                    highest_freq=base_pos[base]
            for base in base_pos.keys():
                if base_pos[base]<variable_pos_freq and \
                 base_pos[base]<highest_freq:
                    del base_pos[base]

        return pos_results
        

    def get_filtered_degenerate_seq(self, filtered_base_freq):
        """ Create degenerate sequence based on filtered base freq 
        
        Creates a bit pattern by using bitwise_or of all bases in each
        position, then uses this bitmasked value in a dictionary of 
        degenerate IUPAC codes to determine the appropriate IUPAC code
        for the base position.
        
        filtered_base_freq: contains filtered version of columnProbs()
         generated base frequencies to remove low frequency occuring bases."""

        
        sequence=[]
        DNA_to_numeric = get_DNA_to_numeric()
        numeric_to_DNA = get_numeric_to_DNA()
        
        for base in filtered_base_freq:
            base_list=[]
            degen_pos=array(0)
            for degen_base in base.keys():
                base_list.append(DNA_to_numeric[degen_base])
            for n in base_list:
                degen_pos=bitwise_or(degen_pos,n)
            sequence.append(numeric_to_DNA[degen_pos])
                            
        return "".join(sequence)

class KnownPrimer(object):
    """ contains data for supplied known primers
    
    These known primers will be used to compare with prospective primers
    to determine if the prospective primers are unique, have overlapping
    regions, or are duplicates (share 3' ends) with the supplied known
    primers"""
    def __init__(self, seq, name, trunc_len):
        """ Data for known primers 
        
        seq: sequence of known primer, 5'->3'
        name: known primer name
        trunc_len: length of 3' end of the primer, used for determine
         if a new primer is a complete duplication of a known primer """
    
        self.name = name
        self.seq = DNA.makeSequence(seq,Name=name)
        self.trunc_seq = convert_to_numeric(self.seq[-trunc_len:])
        if len(self.trunc_seq) < trunc_len:
            raise ValueError("Primer in known primers file %s " % name +\
             "has a length less than the specified 3' length %d " % trunc_len)
        self.numeric_seq = convert_to_numeric(seq)



def analyze_primers(hits_file,
                    output_dir,
                    verbose,
                    variable_pos_freq,
                    sort_method, 
                    known_primers_filepath,
                    primer_name,
                    match_len, 
                    cmp_truncate_len,
                    amplicon_len):
    """ Analyzes the sequences given in the primers hits file 
    
    hits_file: output filepath from generate_primers_denovo.py
    output_dir: output directory
    verbose: verbose output
    variable_pos_freq: frequency at which a particular base has to appear
     to be considered when making the filtered degenerate version of new
     primers.  Lower to make more degenerate, raise to limit degeneracy.
    sort_method: sorts primers in the formatted_primers.txt file to sort by
     sensitivity, specificity, or overall shannon entropy.
    known_primers_filepath: optional primers file of known primers to test for
     overlapping or matching primers with new primers.
    primer_name: optional name to append to the primers generated.  For 
     instance, a primer with the name "Bacterial" would be given a name 
     starting with the index, followed by the primer_name, then 'f' or 'r'.
    match_len: the length of overlap that a known primer must have with a new
     primer to be considered overlapping.
    cmp_truncate_len: length of the 3' end of the primers, used for determining
     if a known and new primer are duplicates, perhaps with variable overall
     lengths.
    amplicon_len: If enabled, will generate an additional file showing forward
     and reverse primer pairs that amplicons with an estimated range in X:Y
     given by this parameter as int values separated by a colon.  Requires
     a standard alignment file was used to generate the primer indices with
     generate_primers_denovo."""

    if not output_dir.endswith("/"):
        output_dir += "/"
        
        
    # Read, generate primer information from given primer hits file
    primers_f = open(hits_file,'U')
    primers_data = get_primers_data(primers_f)
    primers = build_primers(primers_data, variable_pos_freq, primer_name,
     cmp_truncate_len)
     
    if amplicon_len:
        check_for_std_alignment(primers_data)
    
    # Write Primer summary file
    analyzed_primers_f = output_dir+"primers_details.txt"
    primer_report = open(analyzed_primers_f, 'w')
    write_primers_summary(primers, primer_report, variable_pos_freq)
    
    # Write out sorted primers in a valid format for the primer analysis
    # pipeline
    formatted_primers_f = output_dir + "formatted_primers.txt"
    formatted_primers = open(formatted_primers_f, 'w')
    f_primers_sorted, r_primers_sorted = sort_primers_data(primers,
     sort_method)

    write_formatted_primers(f_primers_sorted, r_primers_sorted,
     formatted_primers)
     
    # If a list of known, or universal primers are specified, find out if and
    # how much the de novo primers overlap, output to primers_overlap.txt file
    if known_primers_filepath:
        known_primers = open(known_primers_filepath, "U")
        primers_overlap = open(output_dir + "primers_overlap.txt","w")
        primers = find_matches(primers, known_primers, match_len,
         cmp_truncate_len)
        write_overlap_primers(primers, primers_overlap)
        
    # Generate file with primer pairs that have amplicons in the estimated range
    if amplicon_len:
        amplicon_len_pairs = output_dir + "amplicon_len_pairs.txt"
        amplicon_pairs_f = open(amplicon_len_pairs, "w")
        amplicon_data = get_amplicon_pairs(amplicon_len, f_primers_sorted,
         r_primers_sorted)
        write_amplicon_pairs(amplicon_data, amplicon_pairs_f)


def build_primers(primers_data, 
                  variable_pos_freq, 
                  primer_name,
                  cmp_truncate_len):
    """ Constructs primer objects from primer dictionaries 
    
    primers_data: list of lines, from primers data generated with 
     generate_primers_denovo.py
    variable_pos_freq: frequency at which a particular base has to appear
     to be considered when making the filtered degenerate version of new
     primers.  Lower to make more degenerate, raise to limit degeneracy.
    primer_name: optional primer name to add to numeric index and 'f' or 'r'.
    cmp_truncate_len: length of the 3' end of the primers, used for determining
     if a known and new primer are duplicates, perhaps with variable overall
     lengths."""
    
    
    forward_primers, reverse_primers = parse_primer_data(primers_data)
    
    primers = []
    
    for n in forward_primers:
        header = reverse_primers[n][0]
        forward_seqs = forward_primers[n][1:]
        reverse_seqs = reverse_primers[n][1:]
        primers.append(Primer(header,
         LoadSeqs(data=forward_seqs,moltype=DNA,label_to_name=lambda x: x),
         LoadSeqs(data=reverse_seqs,moltype=DNA,label_to_name=lambda x: x),
         variable_pos_freq, primer_name, cmp_truncate_len))
    
    return primers

    

def check_for_std_alignment(primers_data):
    """ Checks for standard alignment data from primers data """
    
    for line in primers_data:
        if line.startswith('C'):
            try:
                std_index = line.split(",")[7]
            except IndexError:
                raise IndexError('Usage of the -a option requires that a '+\
                 'standard alignment option (see generate_primers_denovo.py) '+\
                 'was used to generate accurate indices of de novo primers.')
            
    
    
def sort_primers_data(primers, sort_method):
    """ sort primers according to analysis type
    
    Returns a list of strings with the primer name <tab> primer sequence
    for the purpose of writing a formatted primer file for the primer
    analysis pipeline 
    
    primers: list of Primer objects
    sort_method: method by which primers will be sorted-sensitivity, 
     specificity, or degeneracy (shannon entropy)."""
    

    
    f_primers_sorted = []
    r_primers_sorted = []
    # Sensitivity, specificity will be equal for prospective forward and reverse
    # primers.  Overall shannon entropy scores, will likely be differenty
    if sort_method == 'S':
        f_primers = []
        r_primers = []
        for p in primers:
            # Generate tuples for sorting via itemgetter
            f_primers.append((p.primer_f_name, p.f_degenerate_seq,
             p.sensitivity))
            r_primers.append((p.primer_r_name, p.r_degenerate_seq,
             p.sensitivity))
        # Sort primers according to sensitivity, from greatest to least
        f_primers = sorted(f_primers, key = itemgetter(2), reverse = True)
        r_primers = sorted(r_primers, key = itemgetter(2), reverse = True)
    elif sort_method == 'O':
        f_primers = []
        r_primers = []
        for p in primers:
            # Generate tuples for sorting via itemgetter
            f_primers.append((p.primer_f_name, p.f_degenerate_seq,
             p.f_overall_score))
            r_primers.append((p.primer_r_name, p.r_degenerate_seq,
             p.r_overall_score))
        # Sort primers according to overall shannon entropy, least to greatest
        f_primers = sorted(f_primers, key = itemgetter(2))
        r_primers = sorted(r_primers, key = itemgetter(2))
    elif sort_method == 'P':
        f_primers = []
        r_primers = []
        for p in primers:
            # Generate tuples for sorting via itemgetter
            f_primers.append((p.primer_f_name, p.f_degenerate_seq,
             p.specificity))
            r_primers.append((p.primer_r_name, p.r_degenerate_seq,
             p.specificity))
        # Sort primers according to specificity, from smallest to largest
        f_primers = sorted(f_primers, key = itemgetter(2))
        r_primers = sorted(r_primers, key = itemgetter(2))
    else:
        raise_(ValueError,('Sort method must be "S", "P", or "O"'))
        
         
    # Generate list in format of a list of strings for writing to file
    for n in f_primers:
        f_primers_sorted.append(str(n[0]) + '\t' + str(n[1]) + '\n')
    for n in r_primers:
        r_primers_sorted.append(str(n[0]) + '\t' + str(n[1]) + '\n')
    
    
    return f_primers_sorted, r_primers_sorted
    
    




def convert_to_numeric(sequence):
    """ convert BP codes to numeric values for bitwise comparisons
    
    returns a numeric list corresponding to the nucleotide sequence
    
    sequence: IUPAC DNA sequence"""
    
    int_mapped_seq=[]
    DNA_to_numeric = get_DNA_to_numeric()
    for n in sequence:
        int_mapped_seq.append(DNA_to_numeric[n])
    return int_mapped_seq
    
def convert_to_DNA(numeric_sequence):
    """ Converts numeric sequence to DNA 
    
    numeric_sequence: DNA sequence that has been converted to a 4 bit version
     for the purpose of doing bitwise operations."""
    
    DNA_seq=""
    numeric_to_DNA = get_numeric_to_DNA()
    
    for n in numeric_sequence:
        DNA_seq+=numeric_to_DNA[n]
        
    return DNA_seq

def find_overlap(match_len,
                 seq1,
                 seq2):
    """ Searching for match_len of perfect matches between seq1 and seq2 
    
    match_len: size, in base pairs, of overlap between known and new primer
     to be considered overlapping.
    seq1: new primer sequence
    seq2: known primer sequence"""
    
    for n in range(len(seq1)-match_len+1):
        seq_slice=seq1[n:n+match_len]
        for m in range(len(seq2)-match_len+1):
            compare_slice=seq2[m:m+match_len]
            seq_bitwise=bitwise_and(seq_slice,compare_slice)
            if len(seq_bitwise.nonzero()[0])==match_len:
                return convert_to_DNA(seq_bitwise.tolist())
                
    return False
            


    
def find_3prime_match(test_primer,
                      known_primer,
                      cmp_truncate_len):
    """ Test for perfect 3' match, if so, append match to primer object 
    
    test_primer: new primer Primer object
    known_primer: known primer CmpPrimer object
    cmp_truncate_len: length of the 3' end of the primers, used for determining
     if a known and new primer are duplicates, perhaps with variable overall
     lengths.
    """
        
    # only test forward to forward, reverse to reverse primers
    # If perfect match at the 3' end, append string to primer object
    # for overlapping primer output later
    if known_primer.name.endswith("f"):
        seq_bitwise=bitwise_and(test_primer.f_trunc_seq, known_primer.trunc_seq)
        if len(seq_bitwise.nonzero()[0])==cmp_truncate_len:
            test_primer.f_unique = False
            test_primer.f_3prime_match.append("\t".join(
                                        (test_primer.primer_f_name,
                                         test_primer.f_degenerate_seq,
                                         known_primer.name,
                                         str(known_primer.seq) + "\n")))
    elif known_primer.name.endswith("r"):
        seq_bitwise=bitwise_and(test_primer.r_trunc_seq, known_primer.trunc_seq)
        if len(seq_bitwise.nonzero()[0])==cmp_truncate_len:
            test_primer.r_unique = False
            test_primer.r_3prime_match.append("\t".join(
                                        (test_primer.primer_r_name,
                                         test_primer.r_degenerate_seq, 
                                         known_primer.name, 
                                         str(known_primer.seq) + "\n")))  
    else:
        raise ValueError('Known primer %s named incorrectly, ' %\
         known_primer.name + 'should end with "f" or "r"')
    

def build_known_primers(k_primers,
                        cmp_truncate_len):
    """ Builds list of KnownPrimers objects from list of known primers
    
    k_primers: list of known primers with tuples of seq, name
    cmp_truncate_len: length of the 3' end of the primers, used for determining
     if a known and new primer are duplicates, perhaps with variable overall
     lengths."""
    
    known_primers = []
    
    # k_primers is a list of tuples, where element 0 is the seq,
    # element 1 is the name of each known primer
    for p in k_primers:
        known_primers.append(KnownPrimer(p[0], p[1], cmp_truncate_len))
        
    return known_primers

def find_matches(primers, 
                 known_primers_data, 
                 match_len,
                 cmp_truncate_len):
    """ Finds and returns info about overlapping primers and known seqs 
    
    primers: list of Primer objects for new primers
    known_primers_data: list of lines of data for known primers
    match_len: length of overlap between known and new primers to be 
     considered overlapping.
    cmp_truncate_len: length of the 3' end of the primers, used for determining
     if a known and new primer are duplicates, perhaps with variable overall
     lengths."""
    
    # Build list of known KnownPrimer objects from the supplied known primers
    k_primers= get_known_primers(known_primers_data)
    known_primers = build_known_primers(k_primers, cmp_truncate_len)
        
    
    # For each prospective primer, compare all four possible orientations of
    # the primer to each known primer for overlap and for a 3' match.
    for p in primers:
        f_overlap = False
        r_overlap = False
        for k in known_primers:
            find_3prime_match(p, k, cmp_truncate_len)
            k_seq = k.numeric_seq
            f_seq = p.f_seq
            f_comp = p.f_complement
            f_rev = p.f_reverse
            f_rc = p.f_rc
            r_seq = p.r_seq
            r_comp = p.r_complement
            r_rev = p.r_reverse
            r_rc = p.r_rc
            # Should only be one overlap for a given primer
            f_overlap = find_overlap(match_len, f_seq, k_seq) or\
             find_overlap(match_len, f_comp, k_seq) or\
             find_overlap(match_len, f_rev, k_seq) or\
             find_overlap(match_len, f_rc, k_seq)
            r_overlap = find_overlap(match_len, r_seq, k_seq) or\
             find_overlap(match_len, r_comp, k_seq) or\
             find_overlap(match_len, r_rev, k_seq) or\
             find_overlap(match_len, r_rc, k_seq)
            if f_overlap:
                p.f_unique = False
                p.f_partial_overlap.append("\t".join((p.primer_f_name,
                                                      p.f_degenerate_seq, 
                                                      k.name, 
                                                      str(k.seq),
                                                      f_overlap +"\n")))
            if r_overlap:
                p.r_unique = False
                p.r_partial_overlap.append("\t".join((p.primer_r_name,
                                                      p.r_degenerate_seq, 
                                                      k.name, 
                                                      str(k.seq), 
                                                      r_overlap +"\n")))
   
    return primers
    


