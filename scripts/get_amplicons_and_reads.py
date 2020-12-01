#!/usr/bin/env python
# File created on 02 Jun 2010
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2010, The Primer Prospector Project"
__credits__ =  ["William Walters"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

from os.path import isdir

from cogent.util.misc import parse_command_line_parameters, create_dir
from optparse import make_option

from primerprospector.get_amplicons_and_reads import get_amplicons_and_reads



script_info={}
script_info['brief_description']="""Generates predicted amplicons and reads from primers"""
script_info['script_description']="""Using hits files (generated with analyze_primers.py), this module will look at an individual primer or a specified primer pair to generate amplicons and reads.  By default, the weighted score is used to determine if a particular primer will amplify, but other results (3' mismatches, overall mismatches) can be used to determine if a primer will amplify.  Every fasta file used to generate the hits file must be passed to this module (multiple fastas should be separated by a colon).  If a single hits file is specified, the amplicon will be the entire sequence following the 3' end of the primer."""

script_info['script_usage']=[]
script_info['script_usage'].append(("""Standard Example usage (examine pair of primer hits files):""","""""","""%prog [options] {-f fasta_filepath -i primer_hits_filepath1:primer_hits_filepath2 [required] }"""))
script_info['script_usage'].append(("""Use total mismatches instead of weighted score to determine if a given primer will amplify:""","""""","""get_amplicons_and_reads.py -f fasta_filepath -i primer_hits_filepath1:primer_hits_filepath2 -s total_mismatches\n\n"""))

script_info['output_description']="""An amplicon file (in fasta format) will be written for the forward and reverse primer pair, with a name based upon the name of the primers (taken from hits file name, e.g. 27f_338r_amplicons.fasta).  The reads file will be generated in the same manner, with an "_reads" instead of "_amplicons"."""


script_info['required_options']=[\
    make_option('-i', '--primer_hits',
        help='Target primer hits files.  Separate multiple files '+\
        'with a colon.'),
        
    make_option('-f', '--fasta_fps',
        help='Fasta filepaths.  Must match the fasta files used in '+\
        'the analyze_primers module.  Multiple fasta files can be '+\
        'passed, separated with a colon.  Order not important.' )]


script_info['optional_options']=[\
    make_option('-o', '--output_dir',
        help='Specify output directory for amplicons and reads. '+\
        '[default: %default]', default="."),
        
    make_option('-s', '--score_type',
    help='Value to use from primer hits file to determine a given '+\
        'primer\'s amplification success.  Valid choices are weighted_score, '+\
        'overall_mismatches, tp_mismatches.  Gibbs energy scores '+\
        'not currently implemented [default: %default]', 
        default='weighted_score'),

    make_option('-t','--score_threshold',\
        help='If primer has score at or below this parameter, the primer '+\
        'amplification is considered to be successful [default: %default]',
        default=1.0),
        
    make_option('-m', '--min_seq_len',\
        help='Sets minimum sequence length of amplicon to be included in the '+\
        'output amplicon file [default: %default]', default=100),
        
    make_option('-d', '--read_direction',\
        help='Direction of reads generated. Can be forward (f), reverse (r),'+\
        ' or paired end (p).  [default: %default]', default='r'),
        
    make_option('-R', '--read_len',\
        help='Length of reads to generate.  Should be set according to '+\
        'sequencing technology/reagents used.  [default: %default]',
        default = 250)]

script_info['version'] = __version__

def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    primer_hits = opts.primer_hits
    fasta_fps = opts.fasta_fps
    output_dir = opts.output_dir
    score_type = opts.score_type
    score_threshold = float(opts.score_threshold)
    min_seq_len = int(opts.min_seq_len)
    read_direction = opts.read_direction
    read_len = int(opts.read_len)
    
    # Check primer score type
    valid_score_types = ['overall_mismatches', 'tp_mismatches', 
     'weighted_score']
    
    if score_type not in valid_score_types:
        raise ValueError, ('Score type must be one of the following: '+\
         '%s' % ", ".join(valid_score_types))
         
    # Check read direction parameter
    valid_read_directions = ['f', 'r', 'p']
    
    if read_direction not in valid_read_directions:
        raise ValueError, ('Read direction must be one of the following: '+\
         '%s' % ", ".join(valid_read_directions))
    
    
    # Create output directory if it does not exist
    create_dir(output_dir)
    
    get_amplicons_and_reads(primer_hits, fasta_fps, output_dir, 
     score_type, score_threshold, min_seq_len, read_direction, read_len)
    
    
    
    
if __name__ == "__main__":
    main()

