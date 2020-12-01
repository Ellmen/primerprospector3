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

from glob import glob
from string import lower, upper

from cogent.parse.fasta import MinimalFastaParser

from primerprospector.util import correct_primer_name

""" Contains parse functions for the various modules of PrimerProspector """
 
def parse_primer_data(primers_data):
    """ Reads primer sequences from list of lines, assembles primer dictionaries
    
    The primer header is read from the list of lines from the input file. 
    These sequences and headers are used to generate dictionaries"""
    
    forward_primers = {}
    reverse_primers = {}
    
    
    for line in primers_data:
        fields = line.split(',')
        # Get primer index, sensitivity, specificity
        if fields[0] == 'C':
            forward_primers[fields[3]] = [",".join(fields[1:])]
            reverse_primers[fields[3]] = [",".join(fields[1:])]
        # Append upstream/downstream sequences for testing shannon entropies
        elif fields[0] == 'M':
            forward_primers[fields[1]].append(fields[2])
            reverse_primers[fields[1]].append(fields[3])
        else:
            raise ValueError,('Invalid primer hits file, see line %s' % line) 
            
    return forward_primers, reverse_primers
    
def get_known_primers(known_primers_data):
    """ Read data from known primers file, builds list of tuples
    
    tuple contains (seq, name) for each primer """
    
    k_primers = []
    
    for n in known_primers_data:
        if n.startswith("#") or not(n.strip()):
            continue
        primer_name = (n.split('\t')[0].strip())
        # Fix primer name if not in correct format (followed by lower case f
        # or r.  Leave other components unchanged.
        primer_name = correct_primer_name(primer_name)
        if not (primer_name.endswith('f') or primer_name.endswith('r')):
            raise ValueError,('primer %s from known primers ' % primer_name +
             'file should end with "f" or "r" for forward or reverse primer')
        # Caps for DNA sequences
        seq = upper(n.split('\t')[1].strip())
        k_primers.append((seq, primer_name))
        
    return k_primers
    
def get_primers_data(primers_infile):
    """ Returns list of lines from specificity_infile """
    
    primers_data=[]
    # Parse infile, skip comments or empty lines and remove white space
    for n in primers_infile:
        if n.startswith("#") or len(n.strip())==0:
            continue
        primers_data.append(n.strip())
        
    return primers_data
    
def parse_formatted_primers_data(primers_data):
    """ Parses primer data, returns list of tuples (primer_name, primer_seq)
    
     primers_data: open file object of formatted primers file"""

    raw_primers = []
    
    # tuples were used instead of a dictionary, as people may be interested
    # in having a particular order of primers analyzed.
    for n in primers_data:
        # Skip comments, empty lines
        if n.startswith("#") or not n.strip():
            continue
        try:
            primer_name=(n.split('\t')[0].strip())
            # Fix primer name if not in correct format (followed by lower case f
            # or r.  Leave other components unchanged.
            primer_name = correct_primer_name(primer_name)
            # Stripped in case leading/trailing spaces are present
            primer_seq=upper(n.split('\t')[1].strip())
        except IndexError:
            raise IndexError,('Incorrect file format for input primers file. '+\
             'Comments should be preceeded by "#" and the primers should be '+\
             'given as primer_name<tab>primer_sequence.')
        raw_primers.append((primer_name, primer_seq))

    return raw_primers
    
def get_fasta_filepaths(fasta_fps):
    """ Tests all filepaths, returns list of fasta filepaths 
    
    fasta_fps: input fasta filepath(s), if multiple paths, each should be
     separated by a colon."""
    
    fasta_filepaths = []

    
    # If multiple filepaths were passed, they are separated by a colon
    fps = fasta_fps.split(":")
    # append all filepaths to the lists
    for f in fps:
        fasta_filepaths.append(f)

        
    # Test each file to make sure it can be opened
    for fp in fasta_filepaths:
        try:
            test_file = open(fp, "U")
            test_file.close()
        except IOError:
            raise IOError, ("Unable to open %s, please check filepath." % fp)
            
    for fp in fasta_filepaths:
        f = open(fp, "U")
        for label, seq in MinimalFastaParser(f):
            # Halt and raise error if gaps or "U" nucleotides in sequences
            if "U" in seq:
                raise ValueError,("Fasta file %s " % fp +\
                "contains 'U' characters, please replace with 'T' "+\
                "characters.  Try running the module clean_fasta.py.")
            if "-" in seq or "." in seq:
                raise ValueError,("Fasta file %s " % fp +\
                "contains gap characters ('-' or '.').  These characters can "+\
                "be removed with the module clean_fasta.py")
        f.close()

            
    return fasta_filepaths
    
    
def build_fasta_data(fasta_fps):
    """ Reads fasta files into fasta_data dictionary """
    
    fasta_data={}
    for fasta_file in fasta_fps:
        f=open(fasta_file,'r')
        for label,seq in MinimalFastaParser(f):
            # Stripping white space from fasta label so will match label
            # in hits files
            fasta_data[label.split()[0]]=seq
            
    return fasta_data
    
def get_primer_hits_data_pair(primer_pair):
    """ Reads primer hits data into 2D list 
    
    primer_pair: 2 element list containing filepath for forward and reverse
     hits file."""
    
    f_primer_hits = open(primer_pair[0], 'U')
    r_primer_hits = open(primer_pair[1], 'U')
    
    hits_data = [[],[]]
    
    for n in f_primer_hits:
        # Skip comments
        if n.startswith("#") or not n.strip():
            continue
        else:
            hits_data[0].append(n)
    for n in r_primer_hits:
        # Skip comments
        if n.startswith("#") or not n.strip():
            continue
        else:
            hits_data[1].append(n)
            
    # Test hits data to make sure forward and reverse hits files match
    if len(hits_data[0]) != len(hits_data[1]):
        raise ValueError, ('Hits file pair does not have equal number of '+\
         'sequences.  Please check that %s and %s have the same input ' %\
         (primer_pair[0], primer_pair[1]) + 'fasta sequences.')
         
    # Make sure all labels match for the primer pair
    for n in range(len(hits_data[0])):
        f_primer_seq_label = hits_data[0][n].split(',')[0]
        r_primer_seq_label = hits_data[1][n].split(',')[0]
        if f_primer_seq_label != r_primer_seq_label:
            raise ValueError,('Mismatched sequence labels for files '+\
             '%s and %s ' % (primer_pair[0], primer_pair[1]) + ' at labels '+\
             '%s and %s ' % (f_primer_seq_label, r_primer_seq_label))
         
        
    return hits_data
    
def get_primer_hits_data(primer_hits):
    """ Reads primer hits data into list 
    
    primer_hits: primer hits filepath"""
    
    try:
        hits_f = open(primer_hits, "U")
    except IOError:
        raise IOError,("Can not open specified hits filepath %s.  " %\
         primer_hits)
    
    hits_data = []
    
    for line in hits_f:
        if line.startswith("#") or not line.strip():
            continue
        else:
            hits_data.append(line)
            
    hits_f.close()
            
    return hits_data
    
def get_primer_seq(hits_fp):
    """ Pulls primer sequence from hits file
    
    hits_fp: primer hits filepath
    """
    
    try:
        hits_f = open(hits_fp, "U")
    except IOError:
        raise IOError,("Can not open specified hits filepath %s.  " %\
         hits_fp)
    
    primer_line_start = "# Primer:"
    
    for line in hits_f:
        if line.startswith(primer_line_start):
            raw_primer_seq = line.split(' ')[-1]
            primer_seq = raw_primer_seq.split('-')[1]
            
    hits_f.close()
            
    return primer_seq
    
def get_hits_field(hits_data,
                   hits_index):
    """ Generates list of lines from hits data, specified by hits_index
    
    hits_data:  Lines of data from hits file.
    hits_index:  Index of data to be returned (data will be split on commas)
    """
    
    # Do initial sanity test in case hits_index illogical
    try:
        first_line = hits_data[0].split(',')[hits_index]
    except IndexError:
        raise IndexError,('Unable to parse hits data, please check the '+\
         'format of the hits file and the hits_index utilized')
         
    hits_field_data = []
    
    for line in hits_data:
        hits_field_data.append(line.split(',')[hits_index])
        
    return hits_field_data
        
    
def parse_taxa_mapping_file(taxa_lines):
    """ Parses list of lines, returns dictionary of seq ID: taxa
    
    taxa_lines: list of lines, e.g., open file object of taxa mapping file.
     taxa must be in the format of seq ID<tab>taxa(separated by semicolons
     in order of decending taxa, starting at domain.
     Example (taxonomy taken from the Silva97 database):
     AB175604<tab>Archaea;Euryarchaeota;Halobacteriales;uncultured
     """
     
    taxa_map = {}
     
    seq_id_index = 0
    taxa_index = 1
     
    for line in taxa_lines:
        
        # skip empty or comment lines
        if line.startswith("#") or not line.strip():
            continue
        
        curr_line = line.split('\t')
        # Check for correct tab separation-only 1 tab per line
        if len(curr_line) != 2:
            raise ValueError,('taxa mapping file not formatted correctly, '+\
             'line %s not properly separated with a single tab.' % line)
             
        # Remove "Root" from taxa if present
        taxa = curr_line[taxa_index].split(';')
        if lower(taxa[0].strip()) == 'root':
            taxa = ";".join(taxa[1:])
        else:
            taxa = ";".join(taxa)
       
        # Strip off leading or trailing whitespace that may be in taxa file
        taxa_map[curr_line[seq_id_index].strip()] = taxa.strip()
    
    return taxa_map
    
def get_amplicons_filepaths(amplicons_filepath,
                            all_files=False):
    """ Returns list of amplicons filepaths
    
    amplicons_filepath: amplicon filepaths-if multiple passed, separate with 
     a colon
    all_files: True/False flag, if True, use amplicons_filepath as a directory
     to search for all files ending with _amplicons.fasta to generate a list
     of filepaths.
    """
    
    
    if not all_files:
        amp_filepaths = amplicons_filepath.split(":")
    else:
        if not amplicons_filepath.endswith("/"):
            amplicons_filepath += "/"
        amp_filepaths = glob(amplicons_filepath + "*_amplicons.fasta")
    
    for amp_file in amp_filepaths:
        try:
            current_f = open(amp_file, "U")
        except IOError:
            raise IOError,('Unable to open %s, please check filepaths.' %\
             amp_file)
        current_f.close()
        
    return amp_filepaths
        
        
def get_hits_files(hits_fps,
                   all_files = False):
    """ Generates list of hits files for analysis 
    
    hits_fps: hits filepaths.  If multiple hits filepaths, separated with a
     colon.
    all_files: If True, test directory specified by the hits_fps parameter 
     for all files ending with _hits.txt"""
    
    
    if all_files:
        if not hits_fps.endswith('/'):
            hits_fps += '/'
        hits_files = glob(hits_fps + "*_hits.txt")
    elif hits_fps.count(":"):
        hits_files = hits_fps.split(":")
    else:
        hits_files = [hits_fps]
    return hits_files
    
def parse_hits_data(primer_hits_data,
                    score_threshold,
                    score_type):
    """ Returns list of lines of hits data whose scores are <= score_threshold
    
    Returning a list instead of a dict, so the order of hits can be 
     maintained.
    
    primer_hits_data: list containing hits data for a given
     primer
    score_threshold: Value at which or below a primer is considering to
     function during PCR.
    score_type: Defines which values from the primer hits file will be used
     to assess primer extension (e.g., weighted_score, overall_mismatches...)
    """
    
    # Current listing of indices used for scoring primer.  Will need to add
    # other indices when Gibbs energy scores are added.
    weighted_score_index = 9
    non_tp_mm_index = 4
    tp_mm_index = 5
    
    hits_results = []
    
    valid_score_types = ['overall_mismatches', 'tp_mismatches', 
     'weighted_score']
     
    if score_type not in valid_score_types:
        raise ValueError,('score_type %s not a valid score type.' % score_type)
        
    if score_type == 'weighted_score':
        target_score_index = [weighted_score_index]
    elif score_type == 'tp_mismatches':
        target_score_index = [tp_mm_index]
    elif score_type == 'overall_mismatches':
        target_score_index = [non_tp_mm_index, tp_mm_index]
    
    for hits in primer_hits_data:
        
        line = hits.split(',')
        
        score = 0
        
        for score_index in target_score_index:
            score += float(line[score_index])

            
        # Don't include hits line if over given score
        if score > score_threshold:
            continue
        else:
            hits_results.append(hits)

        
    return hits_results
    
def get_barcodes(barcodes_f):
    """ Returns list of barcodes from open file object.
    
    barcodes_f: open file of barcodes.  Comments (preceeded by "#") and empty
     lines are ignored.
    """
    
    barcodes = []
    
    for line in barcodes_f:
        curr_line = line.strip()
        
        if len(curr_line) == 0 or curr_line.startswith("#"):
            continue
        barcodes.append(curr_line)
        
    return barcodes
