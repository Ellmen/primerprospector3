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

from cogent import DNA

from primerprospector.parse import get_fasta_filepaths, build_fasta_data,\
 get_primer_hits_data, get_hits_files, parse_hits_data
from primerprospector.util import get_primer_direction
from primerprospector.format import format_linker_data

    
def sort_linker_freqs(base_freqs):
    """ Returns parsed data, best/worst linkers from base_freqs data 
    
    base_freqs:  A list of dictionaries, each containing the 4 nucleotides
     and their frequency of occurance for the given base position.  The first
     dictionary in this list represents the most 5' position in the linker, and
     each consecutive position represents the next base in the linker."""
    
    # Ignores degenerate codes, puts data in summarized and readable format
    linker_summary = ""

    bases = ['A', 'T', 'C', 'G']
    best_seq = ""
    worst_seq = ""
    
    for base_pos in range(len(base_freqs)):

        linker_summary += "Base position %d\n" % base_pos
        A_content = "A: %1.3f " % base_freqs[base_pos]['A']
        T_content = "T: %1.3f " % base_freqs[base_pos]['T']
        C_content = "C: %1.3f " % base_freqs[base_pos]['C']
        G_content = "G: %1.3f " % base_freqs[base_pos]['G']
        
        linker_summary += A_content + T_content + C_content + G_content + "\n"
        
        # find least frequently occuring non-degenerate base
        smallest_freq = min(base_freqs[base_pos]['A'], 
         base_freqs[base_pos]['T'], base_freqs[base_pos]['C'], 
         base_freqs[base_pos]['G'])
        
        largest_freq = max(base_freqs[base_pos]['A'],
         base_freqs[base_pos]['T'], base_freqs[base_pos]['C'],
         base_freqs[base_pos]['G'])
        
        for base in base_freqs[base_pos]:
            if base_freqs[base_pos][base] == largest_freq:
                worst_seq += base
                break
        
        for base in base_freqs[base_pos]:
            if base_freqs[base_pos][base] == smallest_freq:
                best_seq += base
                break
                
    return linker_summary, best_seq, worst_seq
    

        
def get_linker_site_seqs(hits_data,
                         fasta_data,
                         primer_direction,
                         linker_len):
    """ Returns list of seqs from the prospective linker site
    
    hits_data: list of lines of data from hits file.  Should already be parsed
     to remove hits that do not meet score threshold settings
    fasta_data: dictionary of fasta label: sequence
    primer_direction: either 'f' or 'r' to indicate direction of primer.  Used
     to find linker site (upstream or downstream of primer hit site from hits
     file)
    linker_len: size of the linker in base pairs.
    """
    
    if primer_direction not in ['f', 'r']:
        raise ValueError,('Primer direction must be "f" or "r".  Got '+\
         '%s instead.' % primer_direction)
    
    seq_hit_index = 1
    fasta_site_index = 3
    primer_hit_index = 2
    fasta_label_index = 0
    
    linker_site_seqs = []
    
    for hit in hits_data:
        line = hit.split(',')
        
        primer_hit = line[primer_hit_index]
        seq_hit = line[seq_hit_index]
        
        # Skip any hits that have gaps
        if primer_hit.count('-') or seq_hit.count('-'):
            continue
        
        fasta_label = line[fasta_label_index]
        fasta_site = int(line[fasta_site_index])
        
        if primer_direction == 'f':
            linker_start = fasta_site - linker_len
            linker_end = fasta_site
            
            # Test for linker going past the end of the sequence, skip if so
            if linker_start < 0:
                continue
                
            # Slice out linker sequence and append to list
            linker = fasta_data[fasta_label][linker_start:linker_end]
            linker_site_seqs.append(linker)
        else:
            # Have to count length of the primer to find 5' position for linker
            linker_start = fasta_site + len(line[seq_hit_index])
            linker_end = linker_start + linker_len
            
            # Make sure we're not slicing off past the end of the sequence
            if linker_end > len(fasta_data[fasta_label]):
                continue
                
            # Slice out linker sequence and append to list
            # In this case, we want the reverse complement, so we can write
            # suggested linkers in a 5'->3' direction and complementing the
            # appropriate DNA strand
            linker = fasta_data[fasta_label][linker_start:linker_end]
            linker_site_seqs.append(DNA.rc(linker))
        
    
    
    return linker_site_seqs
    
    
def calc_base_freqs(linker_site_seqs,
                    linker_len):
    """ Generates list ofdict of percent base occurace for each base in linker
    
    linker_site_seqs:  List of linker sequence sites from each input sequence
     that passed scoring threshold.
    linker_len: size of linker in base pairs
    """
    
    # Get total number of sequences considered for percent calcuations
    total_seqs = len(linker_site_seqs)

    nucleotides = ['A', 'T', 'C', 'G']
    
    base_counts = []
    
    base_freqs = []
    
    
    for base_pos in range(linker_len):
        base_counts.append({'A':0, 'T':0, 'C':0, 'G':0})
        base_freqs.append({'A':0, 'T':0, 'C':0, 'G':0})
        for seq in linker_site_seqs:
            
            # Add counts of each base for given base index, skip degeneracies
            if seq[base_pos] not in nucleotides:
                continue
            
            base_counts[base_pos][seq[base_pos]] += 1
                
                
        
        # Get percentages for each base
        for base in base_counts[base_pos].keys():
            base_freqs[base_pos][base] =\
             base_counts[base_pos][base] / total_seqs


    
    return base_freqs
    
    
    
    
def get_linkers_fp(output_dir,
                   hits):
    """ Generates suggested linkers output filepath
    
    hits: filepath of input hits file, used to generate output filename
    output_dir: output directory """
    
    if not output_dir.endswith('/'):
        output_dir += '/'
        
    linkers_fp = output_dir + basename(hits).split('.')[0] + \
     "_suggested_linkers.txt"
     
    return linkers_fp
    
    
def generate_linkers(hits,
                     fasta_data,
                     linker_len,
                     output_dir,
                     score_type,
                     score_threshold):
    """ Finds linkers for given hits data, writes suggested linkers file
    
    hits: filepath of input hits file, used to generate output filename
    fasta_data: dict of fasta label: sequence
    linker_len: size of linker in base pairs
    output_dir: output directory
    score_type: weighted_score, tp_mismatches, or overall_mismatches.  Used
     for recording score type in output suggested linkers file.
    score_threshold: Score at which or below a sequence is considered for
     generating linkers.  Logged in suggested linkers file.
    """
    
    # Get hits lines, parse out lines that fail to meet score threshold
    raw_hits_data = get_primer_hits_data(hits)
    hits_data = parse_hits_data(raw_hits_data, score_threshold, score_type)
    
    primer_direction = get_primer_direction(hits)
    
    # Get list of linkers
    linker_site_seqs = get_linker_site_seqs(hits_data, fasta_data, 
     primer_direction, linker_len)
     
    # Get percentage of base frequences in each position of linker
    base_freqs = calc_base_freqs(linker_site_seqs, linker_len)
    
    # Sort results, get best and worst linker sequence
    linker_summary, best_seq, worst_seq = sort_linker_freqs(base_freqs)
    
    # Get output filename
    linkers_fp = get_linkers_fp(output_dir, hits)
    
    
    formatted_linker_data = format_linker_data(linker_summary, best_seq, 
     worst_seq, hits, score_type, score_threshold)
     
    
    linkers_f = open(linkers_fp, "w")
    linkers_f.write(formatted_linker_data)
    linkers_f.close()
    
    

     

    

    


def get_linkers(hits_fps,
                fasta_fps,
                linker_len,
                all_files,
                output_dir,
                score_type,
                score_threshold):
    """ Main program function for generating linkers 
    
    hits_fps: filepath(s) for primer hits to be tested.  Multiple files can
     be separated with a colon.  If the all_files flag is True, this points to
     the directory where all files ending with _hits.txt will be tested.
    fasta_fps: fasta filepaths.  Must include all fasta files used to generate
     the hits files.
    linker_len: length of linker in nucleotides.
    all_files: Set to True to test all _hits file in directory specified with
     the hits_fps parameter.
    output_dir: output directory
    score_type: The type of scoring used to determine if a sequence should
     be considered for making linkers.  Can be weighted_score, tp_mismatches,
     or overall_mismatches
    score_threshold: Value at which or below a sequence score is considered
     acceptable to be considered for making a linker (i.e., that sequence is
     likely to amplify during PCR).
    """
    
    # Check fasta files, build dict of seq ID: sequence
    fasta_fps = get_fasta_filepaths(fasta_fps)
    
    fasta_data = build_fasta_data(fasta_fps)
    
    hits_fps = get_hits_files(hits_fps, all_files)
    
    
    for hits in hits_fps:
        generate_linkers(hits, fasta_data, linker_len, output_dir,
         score_type, score_threshold)

    

    

