#!/usr/bin/env python

__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

from cogent.util.misc import parse_command_line_parameters, create_dir
from optparse import make_option


from primerprospector.optimize_primers import generate_primers_pos_data

script_info = {}
script_info['brief_description'] = """ Records base frequency for each position of a primer for optimization by modulating degeneracy. """
script_info['script_description'] = """

Optimize primers takes an input primer hits file (created by analyze_primers.py)
and generates an output file containing the base frequency for each position 
in the primer.

Only primer hits that exceed a given score threshold will be considered for
the base frequency output.
"""

script_info['script_usage']=[]
script_info['script_usage'].append(("""Standard Example usage:""","""""","""%prog [options] {-i hits_filepath [required]}"""))
script_info['script_usage'].append(("""Use overall mismatches for scoring, set score threshold to only consider primer hits that are equal to or less than two mismatches:""","""""","""optimize_primers.py -i primer2r_hits.txt -s overall_mismatches -t 2"""))

script_info['output_description']="""A tab-seperated text file will be generated
containing the primer sequence, in 5' to 3' orientation, and the base 
frequencies for each nucleotide. """


script_info['required_options']=[\
    # hits filepath
    make_option('-i', '--hits_fp',
        help='Target primer hits file to generate base frequencies with.')]
                 
script_info['optional_options']=[\
        
    # Specify output directory
    make_option('-o', '--output_dir',
        help='Specify output directory for linkers summary. '+\
        '[default: %default]', default = "."),
        
    make_option('-s', '--score_type',
    help='Value to use from primer hits file to determine a given '+\
        'primer\'s amplification success.  Valid choices are weighted_score, '+\
        'overall_mismatches, tp_mismatches.  Gibbs energy scores '+\
        'not currently implemented [default: %default]', 
        default='weighted_score'),

    make_option('-t','--score_threshold',\
        help='If primer has score at or below this parameter, the primer '+\
        'amplification is considered to be successful. [default: %default]',
        default=2.0)
]




script_info['version'] = __version__

def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    hits_fp = opts.hits_fp
    output_dir = opts.output_dir
    score_type = opts.score_type
    score_threshold = float(opts.score_threshold)
    
    # Check primer score type
    valid_score_types = ['overall_mismatches', 'tp_mismatches', 
     'weighted_score']
    
    if score_type not in valid_score_types:
        raise ValueError, ('Score type must be one of the following: '+\
         '%s' % ", ".join(valid_score_types))
    
    # create output directory if not present
    create_dir(output_dir)
    
    generate_primers_pos_data(hits_fp, output_dir, score_type, score_threshold)
    
    
if __name__ == "__main__":
    main()

    
