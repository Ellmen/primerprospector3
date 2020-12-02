#!/usr/bin/env python

__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

from primerprospector.cogentutil.misc import parse_command_line_parameters, create_dir
from optparse import make_option


from primerprospector.taxa_coverage import graph_taxa_coverage


script_info = {}
script_info['brief_description'] = """ Generates taxa coverage graphs, reports about given primers or primer pairs """
script_info['script_description'] = """
Taxa coverage graphs and text file summaries can be generated by this module 
for a primer, multiple primers, or primer pairs.  These graphs are generated by 
using primer data generated by the analyze_primers.py module, to determine if a 
primer will amplify a given sequence.

A single _hits.txt file, multiple hits files, or a directory of hits.txt 
files can be passed as input (see the analyze_primers.py file for 
information about the hits.txt file format).

A taxonomy mapping file is required.  Sequence IDs in the hits files that do
not have corresponding taxonomy mapping data will be binned in a category
of "Unknown" sequences under the domain level output of the graph and 
text summary files.  Taxonomy mapping files have the format of the sequence
ID<tab>taxonomy, delineated by semicolons, starting with domain.

Example (taxonomy comes from the Silva 97 reference set):\n
EU884952<tab>Bacteria;Bacteroidetes-Chlorobi;Bacteroidetes;Rikenella

Forward and reverse primer pairs can be tested (-p option).  In these cases,
amplification for the primers are decided by the poorest score of two primers
tested.  *Warning* This module will not check for logical combinations of
primers, so all forward and reverse primers, including those that are not 
positioned to generate amplicons (i.e., 910f and 495r) will be tested if the
-p option is enabled.
"""

script_info['script_usage']=[]
script_info['script_usage'].append(("""Standard Example usage:""","""""","""%prog [options] {-i hits_filepath [required] -t taxa_mapping_filepath [required]}"""))
script_info['script_usage'].append(("""Change taxa depth reported to 5, manually pass 2 hits files to test, use overall mismatches for scoring:""","""""","""taxa_coverage.py -i primer1f_hits.txt:primer2r_hits.txt -t taxa_map.txt -d 5 -s overall_mismatches """))
script_info['script_usage'].append(("""Test all hits file in a target directory, get primer pair results, put output in taxa_coverage directory:""","""""","""taxa_coverage.py -i hits_dir/ -t taxa_map.txt -p -o taxa_coverage/"""))

script_info['output_description']="""A log file will be generated for each hits file tested in the output directory.  For each primer or primer pair tested, a subdirectory in the main output directory will be generated with taxonomy graphs for each level of taxonomy tested as well as a text file summary of the taxonomy coverage."""


script_info['required_options']=[\
    # hits filepaths
    make_option('-i', '--hits_fps',
        help='Target primer hits files to generate linkers against.  '+\
        'Separate multiple files with a colon.'),
                    
    make_option('-T', '--taxa_fp',
        help='Taxonomy mapping file.')]
                 
script_info['optional_options']=[\
    # Linker size
    make_option('-d', '--taxa_depth', 
        help='Depth of taxa to generate graphs and summaries for, starting '+\
        'with domain. [default: %default]', default = 3),
        
    # Analyze all _hits.txt files in directory specified with -i
    make_option('-r', '--all_files', 
    help='Test all _hits.txt files in directory specified with -i.  '+\
        ' [default: %default]', default = False, action = 'store_true'),
        
    # Look for primer pairs to test
    make_option('-p', '--primer_pairs', 
    help='Test primer pairs.  Will test all input hits files that are forward'+\
        ' and reverse primers.  Hits files must have matching sequences.  '+\
        'The worse scoring primer of the pair dictates amplification success. '+\
        '[default: %default]', default = False, action = 'store_true'),
        
    # Specify output directory
    make_option('-o', '--output_dir',
        help='Specify base output directory for taxa summary.  A log file '+\
        'be output to this directory.  Taxonomy graphs and text summaries '+\
        'will be generated in separated subdirectories from the main output '+\
        'directory. [default: %default]', default = "."),
        
    make_option('-s', '--score_type',
    help='Value to use from primer hits file to determine a given'+\
        'primer\'s amplification success.  Valid choices are weighted_score, '+\
        'overall_mismatches, tp_mismatches.  Gibbs energy scores '+\
        'not currently implemented [default: %default]', 
        default='weighted_score'),

    make_option('-t','--score_threshold',\
        help='If primer has score at or below this parameter, the primer '+\
        'amplification is considered to be successful [default: %default]',
        default=1.0)
]




script_info['version'] = __version__

def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    hits_fps = opts.hits_fps
    taxa_fp = opts.taxa_fp
    taxa_depth = int(opts.taxa_depth)
    primer_pairs = opts.primer_pairs
    all_files = opts.all_files
    output_dir = opts.output_dir
    score_type = opts.score_type
    score_threshold = float(opts.score_threshold)
    
    # Check primer score type
    valid_score_types = ['overall_mismatches', 'tp_mismatches', 
     'weighted_score']
    
    if score_type not in valid_score_types:
        raise ValueError('Score type must be one of the following: '+\
         '%s' % ", ".join(valid_score_types))
         
    if taxa_depth not in range(1, 10):
        raise ValueError('Taxa depth should be a positive integer in the '+\
         'range of 1 to 10.')
    
    # create output directory if not present
    create_dir(output_dir)
    
    graph_taxa_coverage(hits_fps, taxa_fp, taxa_depth, primer_pairs, all_files,
     output_dir, score_type, score_threshold)
    
    
if __name__ == "__main__":
    main()

