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


from primerprospector.generate_linkers import get_linkers

script_info = {}
script_info['brief_description'] = """ Find the best linkers to use for the given primers """
script_info['script_description'] = """

generate_linkers is designed to suggest primer linkers based upon the least frequently occuring base pairs (default 2) immediately upstream from the 5' region of the primer.

A single _hits.txt file, multiple hits files, or a directory of hits.txt 
files can be passed as input (see the analyze_primers.py file for 
information about the hits.txt file format).  The fasta files used to 
generate the _hits.txt files are also required input for this module.

The output file will contain the percentage occurance of each base at 
all positions in the linker and a suggested linker based on complementarity 
to the least frequently occurring bases.
"""

script_info['script_usage']=[]
script_info['script_usage'].append(("""Standard Example usage:""","""""","""%prog [options] {-i hits_filepath [required] -f input_fasta_filepath [required]}"""))
script_info['script_usage'].append(("""Change linker size to 3, manually pass 2 hits files to test, use overall mismatches for scoring:""","""""","""generate_linkers.py -i primer1f_hits.txt:primer2r_hits.txt -f seqs.fasta -l 3 -s overall_mismatches """))
script_info['script_usage'].append(("""Test all hits file in a target directory, put output in linkers directory:""","""""","""generate_linkers.py -i hits_dir/ -f seqs.fasta -o linkers/"""))

script_info['output_description']="""A text file will be generated for each hits file containing details about the base frequency for each position of the linker and the suggested linker sequence."""


script_info['required_options']=[\
    # hits filepaths
    make_option('-i', '--hits_fps',
        help='Target primer hits files to generate linkers with.  '+\
        'Separate multiple files with a colon.'),
                    
    make_option('-f', '--fasta_fps',
        help='Fasta filepath(s).  Must include all fasta sequences used to '+\
        'generate the hits files.  Separate multiple files with a '+\
        'colon.')]
                 
script_info['optional_options']=[\
    # Linker size
    make_option('-l', '--linker_len', 
        help='Size of linker in base pairs. [default: %default]', default = 2),
        
    # Analyze all _hits.txt files in directory specified with -i
    make_option('-r', '--all_files', 
    help='Test all _hits.txt files in directory specified with -i.  '+\
        ' [default: %default]', default = False, action = 'store_true'),
        
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
        default=1.0)
]




script_info['version'] = __version__

def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    hits_fps = opts.hits_fps
    fasta_fps = opts.fasta_fps
    linker_len = int(opts.linker_len)
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
    
    # create output directory if not present
    create_dir(output_dir)
    
    get_linkers(hits_fps, fasta_fps, linker_len, all_files, output_dir,
     score_type, score_threshold)
    
    
if __name__ == "__main__":
    main()

    
