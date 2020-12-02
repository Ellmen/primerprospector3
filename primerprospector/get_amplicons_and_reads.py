#!/usr/bin/env python
# Author: William Walters (William.A.Walters@colorado.edu)
# get_amplicons_and_reads.py

from __future__ import division

__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

""" The get_amplicons_and_reads module is designed to generate amplicons from
specified _hits.txt file or pair of files (generated by the analyze_primers.py 
module).  By default, this module also generates
a sequence "read" of 250 base pairs from the 3' position of the reverse
primer.

If a single hits file is specified, the "amplicon" is considered to be
the sequence of the fasta files starting at the 3' end of the primer and going
to the end of the sequence (direction depends on whether primer is forward
or reverse)

The amplicon files will be generated in a fasta compatable format, with the
sequence ID as the label.  Only hits that have a score equal or lower than the
score threshold (default 1.0) will have the sequence ID and sequence appended.
If a single hits file is specified, the amplicon starts at the 3' end of the
primer and goes to the end of the fasta sequence.
For primer pairs, this module will append the entire length of the predicted 
amplicon between the pair, again with the restrictions in place considering 
the primer hit score for both primers, and the minimal_seq_len.

The read files will contain the sequence ID and sequence.  By default, a 
reverse read of 250 base pairs will be generated (if the amplicons are 
shorted than this, the read length will match the amplicon sequence).

For the read files generated, the files will be named according to the 
primer hits file(s) used to generate the amplicon, with a _f_reads.fasta or
_r_reads.fasta to indicate direction.  If paired end reads are specified, 
both files will be generated.
"""

from os.path import basename

from cogent.parse.fasta import MinimalFastaParser
            
from primerprospector.parse import get_fasta_filepaths, build_fasta_data,\
 get_primer_hits_data_pair, get_primer_hits_data
from primerprospector.util import get_primer_direction
    
def generate_paired_amplicons(primer_hits_data,
                             score_threshold,
                             min_seq_len,
                             fasta_data,
                             score_type):
    """ Returns predicted amplicons from given primer pair 
    
    primer_hits_data: 2 dimensional list containing hits data for a given
     primer pair
    score_threshold: Value at which or below a primer is considering to
     function during PCR.
    min_seq_len: The minimum sequence length required for inclusion of the
     sequence in the output amplicons and reads.
    fasta_data: dictionary of fasta label: sequence for all fasta files
     passed.
    score_type: Defines which values from the primer hits file will be used
     to assess primer extension (e.g., weighted_score, overall_mismatches...)
    """
    
    
    # Current listing of indices used for scoring primer.  Will need to add
    # other indices when Gibbs energy scores are added.
    weighted_score_index = 9
    label_index = 0
    seq_hits_index = 3
    primer_seq_index = 2
    non_tp_mm_index = 4
    tp_mm_index = 5
    
    amplicon_data=[]
    
    valid_score_types = ['overall_mismatches', 'tp_mismatches', 
     'weighted_score']
     
    if score_type not in valid_score_types:
        raise ValueError('score_type %s not a valid score type.' % score_type)
        
    if score_type == 'weighted_score':
        target_score_index = [weighted_score_index]
    elif score_type == 'tp_mismatches':
        target_score_index = [tp_mm_index]
    elif score_type == 'overall_mismatches':
        target_score_index = [non_tp_mm_index, tp_mm_index]
    
    for n in range(len(primer_hits_data[0])):
        
        f_hits = primer_hits_data[0][n]
        r_hits = primer_hits_data[1][n]
        line_f = f_hits.split(",")
        line_r = r_hits.split(",")
        
        # Ensure that fasta labels match between hits files
        if line_f[label_index] != line_r[label_index]:
            raise ValueError('Hits files must be generated from the same '+\
             'fasta file soures.  The following fasta labels are different: '+\
             '%s, %s' % (line_f[label_index], line_r[label_index]))
             
        label = line_f[label_index]
        
        score_f = 0
        score_r = 0
        
        for score_index in target_score_index:
            score_f += float(line_f[score_index])
            score_r += float(line_r[score_index])
            
        # Don't include fasta sequence if either primer over given score
        if score_f > score_threshold or score_r > score_threshold:
            continue
            
        # Need correction for starting the amplicon at the end of the primer
        # Only affects forward primer, reverse primer 3' matches primer hit
        # index
        f_primer_correction=len(line_f[primer_seq_index])
        
        hit_index_f = int(line_f[seq_hits_index]) + f_primer_correction
        hit_index_r = int(line_r[seq_hits_index])


        # Check to see if amplicon meets minimal length requirements
        if (hit_index_r - hit_index_f) < min_seq_len:
            continue
            
        # Slice out amplicon
        try:
            amplicon=fasta_data[label][hit_index_f:hit_index_r]
        except IndexError:
            raise IndexError('Found fasta label %s ' % label + 'in hits '+\
             'file that is not found in the input fasta file(s), please '+\
             'check fasta files to ensure they match those used to generate '+\
             'the hits files.')
             
        # If no amplicon generated, skip
        if not len(amplicon):
            continue
            
        # If all criteria are matched will append fasta label + seq
        amplicon_data.append(">" + label)
        amplicon_data.append(amplicon)
        
    return amplicon_data
    
    
    
def generate_unidirectional_amplicons(primer_hits_data,
                                      score_threshold,
                                      min_seq_len,
                                      fasta_data,
                                      score_type,
                                      primer_direction):
    """ Returns amplicons (to sequence end) from given primer 
    
    primer_hits_data: list containing hits data for a given
     primer
    score_threshold: Value at which or below a primer is considering to
     function during PCR.
    min_seq_len: The minimum sequence length required for inclusion of the
     sequence in the output amplicons and reads.
    fasta_data: dictionary of fasta label: sequence for all fasta files
     passed.
    score_type: Defines which values from the primer hits file will be used
     to assess primer extension (e.g., weighted_score, overall_mismatches...)
    primer_direction: Either 'f' or 'r' to indicate forward or reverse primer
    """
    
    # Current listing of indices used for scoring primer.  Will need to add
    # other indices when Gibbs energy scores are added.
    weighted_score_index = 9
    label_index = 0
    seq_hits_index = 3
    primer_seq_index = 2
    non_tp_mm_index = 4
    tp_mm_index = 5
    
    amplicon_data=[]
    
    valid_score_types = ['overall_mismatches', 'tp_mismatches', 
     'weighted_score']
     
    if score_type not in valid_score_types:
        raise ValueError('score_type %s not a valid score type.' % score_type)
        
    if score_type == 'weighted_score':
        target_score_index = [weighted_score_index]
    elif score_type == 'tp_mismatches':
        target_score_index = [tp_mm_index]
    elif score_type == 'overall_mismatches':
        target_score_index = [non_tp_mm_index, tp_mm_index]
    
    for hits in primer_hits_data:
        
        line = hits.split(',')
        
             
        label = line[label_index]
        
        score = 0
        
        for score_index in target_score_index:
            score += float(line[score_index])

            
        # Don't include fasta sequence if primer over given score
        if score > score_threshold:
            continue
            
        # Need correction for starting the amplicon at the end of the primer
        # Only affects forward primer, reverse primer 3' matches primer hit
        # index
        
        if primer_direction == 'f':
            f_primer_correction=len(line[primer_seq_index])
            hit_index = int(line[seq_hits_index]) + f_primer_correction
        else:
            hit_index = int(line[seq_hits_index])
            
        # Slice out amplicon
        if primer_direction == 'f':
            try:
                amplicon=fasta_data[label][hit_index:]
            except IndexError:
                raise IndexError('Found fasta label %s ' % label + 'in hits '+\
                 'file that is not found in the input fasta file(s), please '+\
                 'check fasta files to ensure they match those used to  '+\
                 'generate the hits files.')
        else:
            try:
                amplicon=fasta_data[label][0:hit_index]
            except IndexError:
                raise IndexError('Found fasta label %s ' % label + 'in hits '+\
                 'file that is not found in the input fasta file(s), please '+\
                 'check fasta files to ensure they match those used to '+\
                 'generate the hits files.')
             
        # If no amplicon generated, skip
        if not len(amplicon):
            continue
            
        # If amplicon below minimum length, skip
        if len(amplicon) < min_seq_len:
            continue
            
        # If all criteria are matched will append fasta label + seq
        amplicon_data.append(">" + label)
        amplicon_data.append(amplicon)
        
    return amplicon_data
    


def get_output_name_single_primer(primer_hit,
                                  output_dir):
    """ Returns output_filepath name for unidirectional primer
    
    primer_hit: primer hits filepath
    output_dir: output directory"""
    
    if not output_dir.endswith('/'):
        output_dir += '/'
    
    amplicons_fp = output_dir + basename(primer_hit).split('_')[0] + "_" +\
     "amplicons.fasta"
    
    
    return amplicons_fp
    
def get_output_name_primer_pair(primer_pair,
                                output_dir):
    """ Returns output_filepath name for primer pair amplicons 
    
    primer_pair: tuple of forward primer hits filepath, 
     reverse primer hits filepath
    output_dir: output directory"""
    
    if not output_dir.endswith('/'):
        output_dir += '/'
    
    forward_name = basename(primer_pair[0]).split('_')[0] + "_"
    reverse_name = basename(primer_pair[1]).split('_')[0] + "_"
    
    amplicons_fp = output_dir + forward_name + reverse_name + "amplicons.fasta"
    
    return amplicons_fp
    


    

def generate_amplicons_from_hits_files(hits_files, 
                                       fasta_data,
                                       output_dir,
                                       score_type,
                                       score_threshold,
                                       min_seq_len):
    """ Matches forward and reverse hits files, generates amplicon fasta file
    
    hits_files: list of hits file filepaths (generated by analyze_primers)
    fasta_data: dictionary of fasta label: sequence
    output_dir: output directory
    score_type: Defines which values from the primer hits file will be used
     to assess primer extension (e.g., weighted_score, overall_mismatches...)
    score_threshold: Value at which or below a primer is considering to
     function during PCR.
    min_seq_len: The minimum sequence length required for inclusion of the
     sequence in the output amplicons and reads.
    """
    
    
    forward_hits = None
    reverse_hits = None
    
    for hit_file in hits_files:
        # Will only get 'f' or 'r' or raised error from get_primer_direction
        primer_direction = get_primer_direction(hit_file)
        if primer_direction == 'f':
            forward_hits = hit_file
        else:
            reverse_hits = hit_file
            
            
    # Flag for single hits file (i.e., 'amplicon' is sequence starting at
    # 3' end of primer and goes to the end of the sequence)
    # Also check to make sure a forward and reverse
    # pair was found if 2 hits files were specified
    if len(hits_files) == 1:
        single_hits = True
    else:
        single_hits = False
        if not forward_hits and reverse_hits:
            raise ValueError('Must specify a forward and reverse primer '+\
             'hits file when passing a pair of primer hits files.')
        
    if single_hits:
        # Pass whichever hits file was found
        primer_hits_data = get_primer_hits_data(forward_hits or reverse_hits)
        
        # Generate list of amplicons, from primer hit to end if above score
        amplicon_data = generate_unidirectional_amplicons(primer_hits_data,
         score_threshold, min_seq_len, fasta_data, score_type, primer_direction)
         
        # Get output amplicon filepath based on input name
        amplicons_fp =\
         get_output_name_single_primer(forward_hits or reverse_hits, output_dir)
        
    else:
        # Read hits data, check for equivalent amounts of data
        primer_hits_data = get_primer_hits_data_pair((forward_hits,
         reverse_hits))
        # Generate list of amplicons to write in fasta format, according
        # to the score threshold for a given primer pair.
        amplicon_data = generate_paired_amplicons(primer_hits_data,\
         score_threshold, min_seq_len, fasta_data, score_type)
        # Get output amplicon filepath based on input names
        amplicons_fp =\
         get_output_name_primer_pair((forward_hits, reverse_hits), output_dir)
         
        
    # Write amplicons_data
    out_f = open(amplicons_fp, "w")
    out_f.write('\n'.join(line for line in amplicon_data))
    out_f.close()
    
    return amplicons_fp
    

def get_output_name_reads(amplicons_fp, 
                          read_direction, 
                          read_len,
                          output_dir):
    """ Generates list of output names based on filename, direction, and len
    
    amplicons_fp: filepath of amplicons file that was written by this module
    read_direction: forward, reverse, or pair end reads
    read_len: length of read to generate in base pairs
    output_dir: output directory
    """
    
    read_output_names = []
    
    if not output_dir.endswith('/'):
        output_dir += '/'
    
    # Cut off 'amplicons' from output name, retain the primer names
    root_name = "_".join(basename(amplicons_fp).split("_")[:-1]) + "_"
    
    if read_direction == 'f':
        direction_name = ['f_']
    elif read_direction == 'r':
        direction_name = ['r_']
    else:
        direction_name = ['f_', 'r_']
        
    len_name = str(read_len) + "_"
    
    suffix = "reads.fasta"
    
    
    for n in direction_name:
        read_output_names.append(output_dir + root_name + n + len_name + suffix)
        
    return read_output_names
        
        
    


def generate_reads(amplicons_fp, 
                   read_direction,
                   read_len,
                   output_dir):
    """ Generate reads from fasta seqs 
    
    amplicons_fp: filepath of amplicons file generated by this module
    read_direction: forward, reverse, or pair end reads (f, r, or p)
    read_len: length of reads in base pairs
    output_dir: output directory
    """
    
    out_f = []
    
    amplicon_f = open(amplicons_fp, "U")
    
    reads_outf = get_output_name_reads(amplicons_fp, read_direction, read_len,
     output_dir)
     
    for reads_out in reads_outf:
        out_f.append(open(reads_out, "w"))
         
    for label, seq in MinimalFastaParser(amplicon_f):
        
        fasta_label = ">" + label + "\n"
        if read_direction == 'f':
            out_f[0].write(fasta_label)
            out_f[0].write(seq[0:read_len] + '\n')
        elif read_direction == 'r':
            out_f[0].write(fasta_label)
            out_f[0].write(seq[-read_len:] + '\n')
        else:
            out_f[0].write(fasta_label)
            out_f[0].write(seq[0:read_len] + '\n')
            out_f[1].write(fasta_label)
            out_f[1].write(seq[-read_len:] + '\n')
            
    
    
def get_hits_files(primer_hits):
    """ Generates list of hits files for analysis 
    
    primer_hits: primer hits filepaths (if multiple, separated by a colon)
    """
    
    # Generate list of hits file paths
    if primer_hits.count(":"):
        hits_files=primer_hits.split(":")
        if len(hits_files) != 2:
            raise ValueError('This module can only accept a single or pair '+\
             'of primer hits files.')
    else:
        hits_files=[primer_hits]
        
    # Test to make sure all filepaths are legitimate
    for f in hits_files:
        try:
            test_f = open(f, "U")
            test_f.close()
        except IOError:
            raise IOError('Incorrect filepath for %s hits file.' % f)
            
    return hits_files

    
def get_amplicons_and_reads(primer_hits, 
                            fasta_fps,
                            output_dir,
                            score_type,
                            score_threshold,
                            min_seq_len,
                            read_direction,
                            read_len):
    """ Main function to generating amplicons and reads 
    
    primer_hits: primer hits filepaths (if two, separated by a colon)
    fasta_fps: Fasta filepaths.  Need all fasta files used to generate the
     primer hits file.  Seperate multiple fasta files with a colon.
    output_dir: Directory where amplicons and reads will be written.
    score_type: Defines which values from the primer hits file will be used
     to assess primer extension (e.g., weighted_score, overall_mismatches...)
    score_threshold: Value at which or below a primer is considering to
     function during PCR.
    min_seq_len: The minimum sequence length required for inclusion of the
     sequence in the output amplicons.
    read_direction: Direction (from 3' end of primer) to generate reads.
    read_len: Length of read to generate."""
    
    
    # Test fasta files for correct format, return list of fasta filepaths
    fasta_fps = get_fasta_filepaths(fasta_fps)
    
    # Build dictionary of fasta label: fasta seq
    # Loading all data this way so there is not a dependency on the order
    # of the sequences used to generate the hits files.
    fasta_data = build_fasta_data(fasta_fps)
    
    # Get list of hits files, test for validity
    hits_files = get_hits_files(primer_hits)
    
    # Generate the amplicons based on primer hit index, score
    amplicons_fp = generate_amplicons_from_hits_files(hits_files, fasta_data,
     output_dir, score_type, score_threshold, min_seq_len)
     
    generate_reads(amplicons_fp, read_direction, read_len, output_dir)
    
  
    


    
