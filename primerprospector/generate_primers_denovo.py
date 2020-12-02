#!/usr/bin/env python
#generate_primers_denovo.py
#author William Walters
#created 7-2-09

from __future__ import division
from __future__ import print_function

from future.utils import raise_
__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

""" generate_primers_denovo is designed to take input aligned
fasta file(s) and search for all Xmers of a certain sequence length 
(default 5) that are common to the sequences in the target
fasta files.
If excluded fasta files are supplied (-e parameter), specificity is tested 
against the excluded fasta file(s).  Any Xmers that are found to exist above
a certain threshold (1% is the default threshold) in the 
excluded fasta sequences are not included in the output for the 
specifity output files.
Next sensitivity is tested.  If the percent of perfectly matching Xmers falls 
below a certain threshold, this particular Xmer is deleted from the list of 
primer objects and is not considered for later analysis.  
"""

from os.path import isdir
from copy import deepcopy

from cogent3.core.alphabet import AlphabetError
from cogent3 import DNA
from optparse import OptionParser
from numpy import array, bitwise_and
from cogent3.parse.fasta import MinimalFastaParser

from primerprospector.util import get_DNA_to_numeric
from primerprospector.format import generate_denovo_output_file, \
 write_denovo_output_file

class ProspectivePrimer(object):
    """ contains data for prospective primers
    
    This object includes the IUPAC sequence of the primer, the aligned
    index where the primer was first encountered, the 
    numeric equivalent of the IUPAC codes (for using numpy's
    bitwise_and), the number of sequences that have an exact match for
    this primer, the percentage of matches overall, and lists of the
    upstream, downstream, and taxonomic information for the primer."""
    def __init__(self,
                 seq,
                 aligned_index,
                 unaligned_index):
        """ Initializes data for ProspectivePrimer object
        
        seq: Xmer sequence being analyzed
        aligned_index: index of first base pair in aligned sequence
        unaligned_index: index of first base pair after degapping """
    
        self.seq=seq
        self.aligned_index=aligned_index
        self.unaligned_index=unaligned_index
        self.numeric_seq=convert_to_numeric(self.seq)
        self.upstream_regions=[]
        self.downstream_regions=[]
        self.labels=[]
        self.match_count=0
        self.percent_match=0
        self.non_specific_hits=0
        self.non_specific_percent=0
        
        self.std_index = False
        self.f_std_index = None
        self.r_std_index = None

        
    def __eq__(self, 
               other):
        """ Override equality operator for testing primer objects """
        return (self.seq == other.seq and
                self.aligned_index == other.aligned_index and
                self.unaligned_index == other.unaligned_index and
                self.numeric_seq == other.numeric_seq and
                self.upstream_regions == other.upstream_regions and
                self.downstream_regions == other.downstream_regions and
                self.labels == other.labels and
                self.match_count == other.match_count and
                self.percent_match == other.percent_match and
                self.non_specific_hits == other.non_specific_hits and
                self.non_specific_percent == other.non_specific_percent)
    

def get_sequence_count(input_fasta_files):
    """ returns count of sequences in given fasta file(s) 
    
    The input_fasta_files is a list of fasta filepaths """
    
    # Correction for the case that only one file passed
    if type(input_fasta_files)==str:
        input_fasta_files=[input_fasta_files]
        
    count=0
    for n in input_fasta_files:
        fasta_f=open(n,'U')
        for label,seq in MinimalFastaParser(fasta_f):
            count+=1
        fasta_f.close()
    return count
    
def get_number_seqs_for_primer(percent_match,
                               seq_count):
    """ Calculates number of sequences needed to include in primer objects
    
    percent_match: sensitivity threshold
    seq_count: total number of seqs in target set"""
    
    total_seq_use=int((1-percent_match)*seq_count)
    
    return total_seq_use

def iterate_target_sequences(input_fasta_files,
                             sequence_length,
                             percent_match,
                             search_range):
    """ Iterates over given fasta files to build target dictionary 
    
    input_fasta_files: List of fasta file paths to search for conserved seqs
    sequence_length: size of the Xmer being analyzed
    percent_match: sensitivity threshold for Xmer
    search_range: optional restrictive search range in the aligned seq for
      conserved sites."""
    
    initial_primers={}

    for n in input_fasta_files:
        # seq_count and total_seq_use based on percent_match parameter to
        # limit the number of sequences searched and optimize performance.
        analyzed_count=0
        seq_count=get_sequence_count(n)
        total_seq_use=get_number_seqs_for_primer(percent_match, seq_count)
        fasta_f=open(n,'U')
        for label,seq in MinimalFastaParser(fasta_f):
            if analyzed_count>total_seq_use:
                break
            analyzed_count+=1
            seq = seq.upper()
            initial_primers=build_seq_data(seq.replace("U","T"),
             sequence_length,initial_primers, search_range)
        fasta_f.close()
    if len(initial_primers)==0:
        raise ValueError('Cannot find any primers from the given fasta '+\
         'files, please check file format, sensitivity/specificity, '+\
         'and search_range parameters.')
    return initial_primers
    
def construct_primers(initial_primers):
    """ Builds list of primer objects from initial_primers
    
    initial_primers: dictionaries containing primers data, used first
     for faster processing. """
    

    primers=[]
    for n in initial_primers:
        primers.append(ProspectivePrimer(n[0],n[1],initial_primers[n]))
    
    return primers
        

def build_seq_data(seq,
                   sequence_length,
                   initial_primers,
                   search_range):
    """ builds dictionary of Xmers from given fasta files
    
    Uses sequences to build a list of objects of Xmers with related
    information: aligned index, numeric version of
    the Xmer for bitwise comparison, and several values initialized to
    zero for later data collection.
    
        
    Because of the inevitable occurance of duplicate Xmers in a 
    sequence, the key had to be made with additional information to 
    ensure its uniqueness.  The total length of the Xmer and up or down
    stream sequences are by default 20 base pairs, and can be altered
    with the -s parameter
    
    seq: current fasta sequence being searched for new Xmers
    sequence_length: size of the Xmer
    initial_primers: dictionary of nascent Xmers
    search_range: optional restrictive search range in the aligned seq for
      conserved sites."""
    
    aligned_seq=DNA.makeSequence(seq)
    # remove gap characters
    unaligned_seq=str(DNA.makeSequence(seq).degap())
    gaps=aligned_seq.gapMaps()
    
    if search_range:
        primer_start = get_corrected_index(seq,int(search_range.split(":")[0]))
        primer_end = get_corrected_index(seq,int(search_range.split(":")[1]))
        # Correct in case end index is close to the end of the sequence
        if primer_end + sequence_length > len(unaligned_seq):
            primer_end = len(unaligned_seq)-sequence_length+1

    else:
        primer_start = 0
        primer_end = len(unaligned_seq)-sequence_length+1
    
    for n in range(primer_start, primer_end):
        seq_slice=unaligned_seq[n:n+sequence_length]
        aligned_index=gaps[0][n]
        unaligned_index=n
        init_key=(seq_slice,aligned_index)
        initial_primers[init_key]=unaligned_index
    
    return initial_primers
    
    

def convert_to_numeric(sequence):
    """ convert DNA codes to numeric values for bitwise comparisons
    
    returns a numeric list corresponding to the nucleotide sequence
    
    sequence: IUPAC DNA sequence"""
    
    int_mapped_seq=[]
    DNA_to_numeric = get_DNA_to_numeric()
    
    for n in sequence:
        int_mapped_seq.append(DNA_to_numeric[n])
    return int_mapped_seq
    
def get_corrected_index(seq,
                        aligned_index):
    """ returns a corrected unaligned index based on aligned index 
    
    seq: aligned IUPAC DNA sequence
    aligned_index: index of Xmer in aligned sequence
    """
    
    # Counts the number of nucleotides in aligned sequence, returns
    # count of nucleotides occuring before aligned index reached
    slice_seq=seq[0:aligned_index]
    # If different gap characters used, may need to modify this
    # In current form, it is optimized for speed
    corrected_index=\
     aligned_index - (slice_seq.count("-") + slice_seq.count("."))
     

     
    return corrected_index
    
    
    
def del_primers(primers,
                deletions):
    """ Deletes primer objects given a list of indices to delete 
    
    primers: list of ProspectivePrimer objects
    deletions: list of indices in the primers list for deletion"""
    
    # Sort primers in reverse order so indices remain correct during deletion
    deletions.sort(reverse=True)
    for n in deletions:
        del primers[n]
    return primers
    
def append_primer_hit(primer, 
                      label,
                      hit_index,
                      region_slice,
                      overall_length,
                      unaligned_seq,
                      primer_len):
    """ Appends upstream and downstream sequence information for primer hit 
    
    Because some sequences may be hit near the 5' or 3' end of sequence read,
    it is necessary to append N's to the upstream or downstream region.  This
    makes both visual inspection of the primers easier and allows for
    alignment objects to be loaded given a list of primers.
    
    primer: conserved Xmer
    label: sequence label from fasta file
    hit_index: unaligned index where Xmer was first found
    region_slice: size of the upstream/downstream slice from conserved Xmer
    overall_length: size of slice + Xmer size
    unaligned_seq: degapped version of the sequence
    primer_len: length of Xmer"""
    
    
    primer.match_count+=1
    primer.labels.append(label.split()[0])
    # Fill in 'N' for incomplete sequences
    # Set primer_index to 0 in case slicing left end of sequence
    primer_index=hit_index-region_slice
    if primer_index<0:
        primer_index=0
    unknown_bases=overall_length-len(unaligned_seq[primer_index:hit_index+
     primer_len])
    if unknown_bases>0:
        filler="-"*unknown_bases
    else:
        filler=""
    upstream_region=filler+unaligned_seq[primer_index:hit_index+primer_len]
    primer.upstream_regions.append(upstream_region)
    unknown_bases=overall_length-len(unaligned_seq[hit_index:hit_index+
     primer_len+region_slice])
    if unknown_bases>0:
        filler="-"*unknown_bases
    else:
        filler=""
    downstream_region=unaligned_seq[hit_index:hit_index +
     primer_len+region_slice]+filler
    primer.downstream_regions.append(downstream_region)
    return
    
def find_specific_primer_matches(primers,
                                 integer_mapped_seq,
                                 deletion_threshold,
                                 seq_count,
                                 sequence_length,
                                 label,
                                 unaligned_seq,
                                 region_slice,
                                 seq):
    """ searches through integer mapped sequence to find specific matches
    
    This function does not append data from sequences, rather its purpose is
    to eliminate non-specific primers before the sensitive primers (along with
    the associated sequence data) are built.
    
    primers: list of ProspectivePrimer objects
    integer_mapped_seq: seq converted to numeric values for bitwise compare
    deletion_threshold: threshold to delete based on number of sequences and
     sensitivity value
    seq_count: number of seqs already analyzed
    label: fasta label for seq
    unaligned_seq: degapped seq
    region_slice: length of seq to slice out upstream/downstream of Xmer
    seq: fasta sequence
    """
   
    primer_len=sequence_length
    overall_length=region_slice+primer_len
    bad_primers=[]
    seq_length=len(integer_mapped_seq)
    
    if len(unaligned_seq)==0:
        raise_(ValueError,('unaligned sequence contains no data.'))
    
    for p in range(len(primers)):
        corrected_index = get_corrected_index(seq,primers[p].aligned_index)
        start_index = corrected_index
        end_index = corrected_index + primer_len
        
        
        # skip test if testing beyond the end of the sequence
        if end_index > seq_length:
            continue
        # Will return all non-zeros with perfect base pair matching
        seq_bitwise = bitwise_and(primers[p].numeric_seq,
         integer_mapped_seq[start_index:end_index])
        if len(seq_bitwise.nonzero()[0])==primer_len:
            primers[p].non_specific_hits +=1
        if primers[p].non_specific_hits>deletion_threshold:
            bad_primers.append(p)

            
    del_primers(primers,bad_primers)
    return primers


def find_sensitive_primer_matches(primers,
                                  integer_mapped_seq,
                                  deletion_threshold,
                                  seq_count,
                                  sequence_length,
                                  label,
                                  unaligned_seq,
                                  region_slice,
                                  seq):
    """ searches through integer mapped sequence to find matches
    
    primers: list of ProspectivePrimer objects
    integer_mapped_seq: seq converted to numeric values for bitwise compare
    deletion_threshold: threshold to delete based on number of sequences and
     specificity value
    seq_count: number of seqs already analyzed
    label: fasta label for seq
    unaligned_seq: degapped seq
    region_slice: length of seq to slice out upstream/downstream of Xmer
    seq: fasta sequence"""
   
    quality_threshold=seq_count-deletion_threshold
    primer_len=sequence_length
    overall_length=region_slice+primer_len
    
    
    bad_primers=[]
    seq_length=len(integer_mapped_seq)
    if len(unaligned_seq)==0:
        raise_(ValueError,('unaligned_seq contains no data.'))
    
    for p in range(len(primers)):
        corrected_index = get_corrected_index(seq,primers[p].aligned_index)
        start_index = corrected_index
        end_index = corrected_index + primer_len
        
        # skip test if testing beyond the end of the sequence
        if end_index > seq_length:
            # This counts as a miss, so do miss check
            if primers[p].match_count<quality_threshold:
                bad_primers.append(p)
            continue
        
        seq_bitwise = bitwise_and(primers[p].numeric_seq,
         integer_mapped_seq[start_index:end_index])
        if len(seq_bitwise.nonzero()[0])==primer_len:
            append_primer_hit(primers[p],label,start_index,region_slice,
             overall_length,unaligned_seq,primer_len)
        if primers[p].match_count<quality_threshold:
            bad_primers.append(p)

    del_primers(primers,bad_primers)
    
    return primers
            
    
def get_deletion_threshold(percent_match,
                           seq_total):
    """ returns max number of missed seqs to be above threshold 
    
    percent_match: Threshold of required percent matching sequences
    seq_total: total number of sequences being tested."""
    

    return int(float(1-percent_match)*seq_total)

def get_specific_hits(primers,
                      exclude_fasta_files,
                      specificity_max,
                      sequence_length,
                      region_slice,
                      seq_total_exclude):
    """ Iterates through excluded sequences, deletes non-specific primers 
    
    primers: list of ProspectivePrimer objects
    exclude_fasta_files: list of fasta filepaths of sequences where exact 
     match with Xmer is not desired
    specificity_max: specificity threshold for unwanted sequences
    sequence_length: size of Xmer
    region_slice: length of seq to slice out upstream/downstream of Xmer
    seq_total_exclude: total number of seqs in exclusion files"""
    
    seq_count=0
    # Once sequence is found deletion_threshold number of times in excluded
    # fasta sequences, delete the primer as being nonspecific
    deletion_threshold=int(round(specificity_max*seq_total_exclude))
    for n in exclude_fasta_files:
        fasta_f=open(n,'U')
        for label,seq in MinimalFastaParser(fasta_f):
            seq_count+=1
            unaligned_seq = seq.replace("-","")
            unaligned_seq = unaligned_seq.replace(".","")
            unaligned_seq = unaligned_seq.replace("U","T")
            unaligned_seq = unaligned_seq.upper()
            integer_mapped_seq = convert_to_numeric(unaligned_seq)
            primers=find_specific_primer_matches(primers, integer_mapped_seq,
             deletion_threshold, seq_count, sequence_length,
             label, unaligned_seq, region_slice,seq)
        fasta_f.close()
    
    return primers
    
def get_sensitive_hits(primers,
                       input_fasta_files,
                       percent_match,
                       sequence_length,
                       region_slice):
    """ Retains and appends information to primer objects w/ sensitive hits 
    
    primers: list of ProspectivePrimer objects
    input_fasta_files: list of fasta filepaths of sequences where exact 
     match with Xmer is desired
    percent_match: sensitivity threshold
    sequence_length: size of Xmer
    region_slice: length of seq to slice out upstream/downstream of Xmer"""

    seq_count=0
    for n in input_fasta_files:
        seq_total_target=get_sequence_count(n)
        deletion_threshold=get_deletion_threshold(percent_match,
         seq_total_target)
        fasta_f=open(n,'U')
        for label,seq in MinimalFastaParser(fasta_f):
            seq_count+=1
            unaligned_seq = seq.replace("-","")
            unaligned_seq = unaligned_seq.replace(".","")
            unaligned_seq = unaligned_seq.upper()
            unaligned_seq = unaligned_seq.replace("U","T")
            integer_mapped_seq = convert_to_numeric(unaligned_seq)
            primers=find_sensitive_primer_matches(primers, integer_mapped_seq,
             deletion_threshold, seq_count, sequence_length,
             label,unaligned_seq, region_slice, seq)
        fasta_f.close()
    
    return primers
    
    
    
def calculate_percent_match(primers,
                            seq_count,
                            exclude_seq_count=1):
    """ Iterates list of primer objects, calculates percent matches """
    # Calculate percent of sequences that are 'hit' by each primer
    for n in range(len(primers)):
        # Calculate percent perfect match
        primers[n].percent_match=float(primers[n].match_count/seq_count)
        primers[n].non_specific_percent=\
        float(primers[n].non_specific_hits/exclude_seq_count)
    
    return primers

def append_std_aligned_index(primers,
                             standard_index_seq,
                             region_slice):
    """ Appends standard aligned index value to ProspectivePrimer objects
    
    primers: list of ProspectivePrimer objects
    standard_index_seq: aligned standard index sequence for giving meaningfule
     index value for forward and reverse primers
    region_slice: number of base pairs up and downstream to include in output"""
    
    for n in primers:
        n.std_index = True
        standard_unaligned_index = get_corrected_index(standard_index_seq,
         n.aligned_index)
        # 5' for forward primer would be upstream of the Xmer by the
        # number of bases in the region slice
        n.f_std_index = standard_unaligned_index - region_slice
        # 5' for reverse primer is the length of the Xmer plus the number
        # of bases in the region slice.
        n.r_std_index = standard_unaligned_index + len(n.seq) + region_slice
        
    return primers

def search_sequences(input_fasta_filepath, 
                     sequence_length,
                     exclude_fasta_filepath,
                     verbose,
                     percent_match,
                     full_primer_length,
                     output_f,
                     specificity_threshold,
                     log_filepath, 
                     standard_index_file, 
                     search_range):
    """ Search sequences for occurences of Xmers of sequence_length 
    
    input_fasta_filepath: fasta file(s) containing target sequences for 
     finding conserved sites
    sequence_length: size of the conserved Xmer to find
    exclude_fasta_filepath: fasta file(s) containing sequences where a perfect
     match to the conserved Xmer is not desired
    verbose: enables printing of status to stdout
    percent_match: sensitivity threshold
    full_primer_length: size of the overall primers, including conserved Xmer
     and upstream/downstream slice
    output_f: filepath for output prospective primers data
    specificity_threshold: threshold in unwanted sequences, sequences found
     above this value will result in primer being discarded
    log_filepath: optional log output, potentially used on a cluster
     environment where stdout does not yield useful results
    standard_index_file: Aligned fasta file that will be used to give 
     numeric value to forward and reverse primers.  For example, if the 
     aligned index of a particular primer matches the unaligne index of 500 in
     the standard index file, and the overall primer size is 20, the forward
     primer will be named 480f and the reverse primer will be named 520r
    search_range: optional limitation to the range in the aligned input 
     sequences where the primers will be searched for.  If one desired to 
     only search for primers at the beginning of a group of 16S sequences, this
     parameter could be used to restrict the search to this range.
    """
    
    # Check input and output files before generating data

    if isdir(output_f):
        raise IOError('%s is a directory, please specify a file path.' \
         % output_f)
    
    try:
        output_filepath=open(output_f, 'w')
    except IOError:
        raise IOError('Unabled to open output filepath %s' %\
         output_f)
        
    if standard_index_file:
        try:
            test_alignment_file = open(standard_index_file, "U")
            test_alignment_file.close()
        except IOError:
            raise IOError('Unable to open standard index file %s'%\
             standard_index_file)
             
    if log_filepath:
        if isdir(log_filepath):
            raise IOError('log_filepath %s is a directory, please specify '+\
             'a filepath.' % log_filepath)
        try:
            test_log_f = open(log_filepath, 'w')
        except IOError:
            raise IOError('Unable to open log file %s' %\
             log_filepath)
    
    region_slice=full_primer_length-sequence_length
       
    
    if log_filepath:
        log_f = open(log_filepath, 'w')
    if verbose:
        print("Building prospective primers")
    if log_filepath:
        log_f.write("Building prosective primers\n")
        
    input_fasta_files=input_fasta_filepath.split(":")
    initial_primers=iterate_target_sequences(input_fasta_files,sequence_length,\
    percent_match, search_range)
    
    if verbose:
        print("Constructing primer objects")
    if log_filepath:
        log_f.write("Constructing primer objects\n")

    primers=construct_primers(initial_primers)

    if exclude_fasta_filepath:
        exclude_fasta_files=exclude_fasta_filepath.split(":")
    else:
        if not exclude_fasta_filepath:
            # Setting variable to 1 in case no exclusion files
            # Limits need for redundant functions
            seq_total_exclude=1
  
    if verbose and exclude_fasta_filepath:
        print("Counting sequences for excluded fasta file(s)")
    if log_filepath:
        log_f.write("Counting sequences for excluded fasta file(s)\n")

    if exclude_fasta_filepath:
        seq_total_exclude=get_sequence_count(exclude_fasta_files)
    if verbose and exclude_fasta_filepath:
        print("Total sequences: %d" % seq_total_exclude)
    if log_filepath and exclude_fasta_filepath:
        log_f.write("Total sequences: %d\n" % seq_total_exclude)
        
    if verbose and exclude_fasta_filepath:
        print("Finding specific hits")
    if log_filepath and exclude_fasta_filepath:
        log_f.write("Finding specific hits\n")
    
    if exclude_fasta_filepath:
        primers=get_specific_hits(primers,exclude_fasta_files,\
        specificity_threshold,sequence_length,region_slice,\
        seq_total_exclude)
    
    seq_total_target=get_sequence_count(input_fasta_files)
    if verbose:
        print("Total number of target sequences: %d" % seq_total_target)
    if log_filepath:
        log_f.write("Total number of target sequences: %d\n" \
        % seq_total_target)

    if verbose:
        print("Finding sensitive primer regions.")
    if log_filepath:
        log_f.write("Finding sensitive primer regions.\n")
       
    primers=get_sensitive_hits(primers,input_fasta_files,\
    percent_match,sequence_length,region_slice)
    primers=calculate_percent_match(primers,seq_total_target,seq_total_exclude)
    
    if standard_index_file:
        standard_index_fasta = open(standard_index_file, "U")
        # Only read first file
        for label, seq in MinimalFastaParser(standard_index_fasta):
            standard_index_seq = seq
            break
        primers = append_std_aligned_index(primers, standard_index_seq,
         region_slice)
        
    else:
        standard_index_seq = None
    
    
    generate_denovo_output_file(primers,output_filepath,\
    specificity_threshold, region_slice, standard_index_seq, percent_match,
    bool(exclude_fasta_filepath))
    
    if verbose:
        print("Module complete")
    if log_filepath:
        log_f.write("Module complete\n")
    



    
        
        
    
