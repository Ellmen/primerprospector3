#!/usr/bin/env python
# File created 18 Mar 2009
from __future__ import division

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

from primerprospector.analyze_primers import analyze_primers


script_info={}
script_info['brief_description']=""" Score primers for binding to a set of target sequences. """
script_info['script_description']="""

This module performs a local alignment for the input primer(s) against target 
sequences, and then determines a weighted score based upon the number of 
mismatches and gaps.  A summary graph showing mismatches and weighted scores is 
generated for all input fasta files, as is a hits file containing details about 
the primer mismatches, index in the sequence, and other details about primer 
binding.

This module takes an input primer file and one or more fasta files.  Each
primer is tested against every sequence to find the best local alignment.
Mismatches and gaps are calculated for the primer, along with a weighting score 
which gives larger penalties to gaps and mismatches in the 3' end of the primer.

An output hits file is generated for each primer, recording information about
the primer hit site, mismatches, and overall weighted score (a perfect score
starts at zero and increases as penalties are added).  A graph is
also generated, showing mismatches/gaps and overall score information for the 
primer and the target sequences.

The primers input file should be generated with the following format:
Comments are preceeded by a pound "#" symbol.
The primer data are tab delineated with the primer name first, such as
"349_v2r", the actual nucleotide sequence next, listed in a 5' to 3' sequence, 
(example: "AATCGRACGNTYA"), and finally a comment 
or citation, if any, can be listed.  Forward primers should be followed by a 
"f" while reverse primers are followed by "r".
A complete example line is listed below.

815_v34f	GTGGCCNATRRCYAGAACGC	Darrow,Scopes,Bryan et al. 1926

The input sequences should be in fasta format.  If more than one file is 
supplied, they should be separated by a colon.  Each fasta file passed will
have its sequence coverage displayed in separate output graphics and hits
files.
"""
script_info['script_usage']=[]


script_info['script_usage'].append(("""Example""","""Standard Example:""","""%prog [options] {-P input_primers_filepath [required] -f input_fasta_filepath [required]}"""))
script_info['script_usage'].append(("""Manually specify a primer name and sequence:""","""Note - primer name must end with 'f' or 'r':""","""analyze_primers.py -p "primer_name_f" -s "ACCTGACRGGTAATC" -f input_fasta_filepath"""))
script_info['script_usage'].append(("""Use multiple target files, change scoring parameters:""","""Pass a primers file, two target fasta files, change the size of the 3' region from the default 5 bases to 7 bases, and lower the 3' mismatch penalty from the default 1.0 to 0.6:""","""analyze_primers.py -P primers.txt -f bacterial_seqs.fasta:eukaryotic_seqs.fasta -e 7 -t 0.6\n\n"""))

script_info['output_description']=("""For each primer tested, an output _hits.txt file containing information about the primer hit to each sequence is generated, as well a .ps file showing overviews for the mismatches and weighted score for the primer and target sequences.  Both output files are named by the primer and fasta file tested""")

script_info['required_options']=[\
    # input fasta filepath(s)
    make_option('-f', '--fasta_seqs',
                help='Target fasta file(s) to score input primer(s) against.'+\
                    ' Separate multiple files with a colon.')]
                    
script_info['optional_options']=[\
    # Primers file.  Must have primers file, or specify a single primer and seq
    make_option('-P', '--primers_filepath', 
        help='Path to input primers file.  This file is tab delineated, with'+\
        " the first column being the primer name, which must end with an 'f'"+\
        " or a 'r'.  The second column contains the primer sequence in 5' to "+\
        " 3' format. One must supply either a primer file or a primer name "+\
        " (-p parameter) and primer sequence (-s parameter). "+\
        "[default: %default]", default=None),
        
    # Option to specify a single primer to analyze
    make_option('-p', '--primer_name', 
    help='Specify a single primer to analyze.  One can either specify a'+\
        ' single primer that is listed in a primers file (-P parameter) or '+\
        ' specify a sequence with the -s parameter.  Primer name must '+\
        'end with a "f" or "r" to designate orientation. [default: %default]', 
        default=None),
        
    # Option to specify a single primer sequence
    make_option('-s', '--primer_sequence',
        help='Primer sequence if using the -p option.  Must be specified'+\
        ' if not passing a primers file with the -P option.  If both -P and '+\
        '-p parameters are passed, the sequence given with this option will '+\
        'be taken rather than sequences in the -P primers file. '+\
        '[default: %default]', default=None),
        
    # Specify output directory
    make_option('-o', '--output_dir',
        help='Specify output directory for hits files and primer summary '+\
        'graphs. [default: %default]', default="."),
   
    # specify length of 3' primer region
    make_option('-e','--three_prime_len',\
        help="Length of primer considered to be part of the 3' region "+\
        'for the purpose of giving a weighted score for mismatches and/or '+\
        'gaps. [default: %default]', default=5),
        
    # penalty for last base mismatch
    make_option('-l', '--last_base_mismatch',\
        help="Sets penalty for mismatch in final base of 3' end of the "+\
        'primer. [default: %default]', default=3),
        
    # Specify penalty for all 3' mismatches except very last base
    make_option('-t', '--three_prime_mismatch',\
        help="Penalty for all 3' mismatches except final base."+\
        '[default: %default]', default=1),
        
    # Specify penalty for all non-3' mismatches
    make_option('-T', '--non_three_prime_mismatch',\
        help= "Penalty for all non-3' mismatches. [default: %default]",
        default=0.4),
        
    # Set penalty for gaps in 3' region of primer
    make_option('-g', '--three_prime_gap',\
        help="Penalty for gaps in the 3' region of the primer. "+\
        " [default: %default]",default=3),
        
    # Set penalty for gaps in non 3' region of primer
    make_option('-G', '--non_three_prime_gap',\
        help="Penalty for non 3' gaps. [default: %default]", default=1)]



script_info['version'] = __version__

def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    # For sake of reducing variable name length, tp = "three_prime"
    # mm = mismatch
    fasta_seqs = opts.fasta_seqs
    verbose = opts.verbose
    output_dir = opts.output_dir
    primers_filepath = opts.primers_filepath
    primer_name = opts.primer_name
    primer_sequence = opts.primer_sequence
    tp_len = int(opts.three_prime_len)
    last_base_mm = float(opts.last_base_mismatch)
    tp_mm = float(opts.three_prime_mismatch)
    non_tp_mm = float(opts.non_three_prime_mismatch)
    tp_gap = float(opts.three_prime_gap)
    non_tp_gap = float(opts.non_three_prime_gap)
    
    # Create output directory if it does not exist
    create_dir(output_dir)
            
    analyze_primers(fasta_seqs, verbose, output_dir, primers_filepath, 
     primer_name, primer_sequence, tp_len, last_base_mm, tp_mm, non_tp_mm, 
     tp_gap, non_tp_gap)
    


if __name__ == "__main__":
    main()
    




