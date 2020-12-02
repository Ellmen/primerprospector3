#!/usr/bin/env python
# created 11-24-09

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

from primerprospector.amplicons_histograms import generate_histograms


script_info={}
script_info['brief_description']="""Create histograms of predicted amplicons for primer pairs"""

script_info['script_description']="""Using amplicons files (generated with get_amplicons_and_reads.py), this module will generate histogram(s) showing the predicted amplicon sizes.  If a taxonomy mapping file is passed, the histograms will be separated according to domain (archaea, bacteria, and eukaryotic)."""

script_info['script_usage']=[]
script_info['script_usage'].append(("""""","""Standard Example usage (create a histogram from an amplicon file, output to the current directory):""","""%prog [options] {-f amplicons_filepath }"""))
script_info['script_usage'].append(("""""","""Test for all _amplicons.fasta files in current directory, pass a taxonomy mapping file so histograms are plotted according to domain, and output to amplicons_graph directory:""","""amplicons_histograms.py -f . -r -t taxonomy_mapping.txt -o amplicons_graph"""))

script_info['output_description']="""There will be an output graph in postscript format with the same root name as the input _amplicons.fasta file for each _amplicons.fasta file passed."""


script_info['required_options']=[\
    make_option('-f', '--amplicons_filepath',
        help='Target amplicons files.  Separate multiple files '+\
        'with a colon.')]
        
script_info['optional_options']=[\
    make_option('-o', '--output_dir',
        help='Specify output directory for amplicons and reads. '+\
        '[default: %default]', default="."),
        
    make_option('-r', '--all_files',
    help='Generate histograms for all files ending with _amplicons.fasta '+\
        'in the directory specified with the -f parameter [default: %default]', 
        default=False, action = "store_true"),

    make_option('-t','--taxa_map',\
        help='Filepath to taxonomy mapping file, used to separate graphs '+\
        'according to domain. [default: %default]',
        default=None)]
        
        
script_info['version'] = __version__


def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    amplicons_filepath = opts.amplicons_filepath
    output_dir = opts.output_dir
    all_files = opts.all_files
    taxa_map_filepath = opts.taxa_map
    
    # Test to make sure directory is passed for amplicons_filepath if all_files
    # is True
    if all_files:
        if not isdir(amplicons_filepath):
            raise ValueError('-f needs to specify a directory ' +\
             'if -r option enabled.')
    
    # Create output directory if it does not exist
    create_dir(output_dir)
    
    generate_histograms(amplicons_filepath, output_dir, all_files, 
     taxa_map_filepath)
    
    
    
    
if __name__ == "__main__":
    main()
    
