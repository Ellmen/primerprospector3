#!/usr/bin/env python
# File created on 02 Jun 2010
from __future__ import division

from future.utils import raise_
__author__ = "William Walters"
__copyright__ = "Copyright 2010, The Primer Prospector Project"
__credits__ =  ["William Walters"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"
 
from os.path import isdir

from cogent.util.misc import parse_command_line_parameters
from optparse import make_option

from primerprospector.generate_primers_denovo import search_sequences


script_info={}
script_info['brief_description']="""Find conserved and/or specific DNA sequences for use as the 3' end of a primer"""
script_info['script_description']="""The generate_primers_denovo module is designed to take an input sequence length (default 5) and aligned fasta file(s) to search for all Xmers of the given sequence length that are conserved in the target fasta files.  Optionally, any Xmers that are found to exist above a certain threshold (1% is the default) in the excluded fasta sequences are discarded.  The remaining Xmers, along with their upstream and downstream sequences are written to an output file."""

script_info['script_usage']=[]
script_info['script_usage'].append(("""Standard Example usage:""","""""","""%prog [options] {-i include_fasta_filepath(s) -o output_primers_filepath}"""))
script_info['script_usage'].append(("""Look for common 5mers in a given aligned fasta file, bact_sample.fasta, that are conserved at least 60% (default setting) of the time.  Output to primers.txt:""","""""","""generate_primers_denovo.py -i bact_sample.fasta -o primers.txt"""))
script_info['script_usage'].append(("""Look for common 6mers in two input aligned fasta files (bact_sample.fasta, arch_sample.fasta), that are not found in the aligned fasta file euk_sample.fasta:""","""""","""generate_primers_denovo.py -i bact_sample:arch_sample.fasta -x 6 -e euk_sample.fasta -o primers.txt"""))

script_info['output_description']="""The output file from generate_primers_denovo is a text file containing information about the target Xmers that met the sensitivity (and specificity if given) threshold(s).  For each Xmer, the calculated sensitivity and specificity values are given, as well as the upstream and downstream sequences from the Xmer for all of the sequences with perfect matches."""

script_info['required_options']=[\
    make_option('-i', '--target_seqs', dest='target_seqs',    #input_fasta_filepath
                help='Target aligned fasta sequence files to find conserved '+\
                    'sites for primer design.  Separate multiple files with a '+\
                    'colon.'),
    make_option('-o', '--output_filepath', dest='output_filepath',  #output_root_name
                help='name of output filepath to write details about '+\
                    'conserved sequence sites.' )]
                
script_info['optional_options']=[\
    make_option('-e', '--exclude_fasta', dest='exclude_fasta_filepath', 
        help='Excluded aligned fasta file(s).  To pass '+\
        "multiple files, separate each file with a colon.  Example: -e "+\
        "test1.fasta:test2.fasta.  If not specified, will skip exclusion "+\
        "step [default: %default]", default=None),
        
    # Percentage of samples that must have perfect matches
    make_option('-p', '--percent_match', dest='percent_match',
    help='Percentage of sequence matches to '+\
        'primer that must match in order to retain prospective '+\
        'sequence in dictionary. [default: %default]', default=0.60),
        
    # Number of base pairs up and downstream to include with sequence
    make_option('-s', '--full_primer_length', dest="full_primer_length",
        help='Overall primer length to '+\
        'retrieve from sequences. [default: %default]', default=20),

    # Number of conserved bases to search for
    make_option('-x','--xmer_length',\
        dest="sequence_length",help='Xmer length to search for in target '+\
        'fasta sequence(s). [default: %default]', default=5),
        
    make_option('-S', '--specificity_threshold',\
        dest="specificity_threshold",\
        help='Sets specificity threshold for excluded fasta sequences. '+\
        '[default: %default]', default=0.01),
        
    make_option('-l', '--log_file',\
        dest='log_file',help='log filepath. If not specified, no log file '+\
        'will be written.  [default: %default]',
        default=None),
        
    make_option('-a', '--standard_index_file',\
        dest='standard_index_file', help='Aligned sequence file with which '+\
        'to assign prospective primer indices to.  The alignment where a '+\
        'conserved sequence is found will be used to determine the unaligned '+\
        'index in the supplied file (for instance an E. coli sequence) '+\
        'and will be recorded in the output file for the purpose of giving '+\
        'a meaningful name to prospective primers.  Only the first sequence '+\
        'in the file will be used for determining an index [default: %default]',
        default=None),
        
    make_option('-r', '--search_range',\
        dest='search_range', help='Range of nucleotides in the supplied '+\
        'aligned target sequences to search for primers.  Supply the '+\
        'starting index and end index separated by a colon.  Example -r 1500:'+\
        '2700  Enable this option to generate primers that target certain '+\
        'regions. [default: %default]', default=None)]


    
script_info['version'] = __version__


def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if float(opts.specificity_threshold<0.001) or \
     float(opts.specificity_threshold)>0.9:
        raise ValueError('Specificity Thresholds expected to be in '+\
        'range 0.001 to 0.90')
        
    if int(opts.sequence_length) not in range(1,15):
        raise_(ValueError,\
    ('Sequence length argument expected to be in range 1-15'))
    if float(opts.percent_match) < 0.01 or float(opts.percent_match)>1.0:
        raise_(ValueError,\
    ('Percent match argument expected to be in range 0.1-1.0'))

    if int(opts.full_primer_length) not in range (10,40):
        raise_(ValueError,('primer length expected to be 10-40 base pairs.'))
        
    if opts.search_range:
        try:
            start_index = int(opts.search_range.split(":")[0])
            end_index = int(opts.search_range.split(":")[1])
        except IndexError:
            raise IndexError('search_range option needs to be two integers '+\
             'separated by a colon.  Example:  -r 1000:2500')
        
    search_sequences(input_fasta_filepath = opts.target_seqs,
     sequence_length = int(opts.sequence_length),
     exclude_fasta_filepath = opts.exclude_fasta_filepath,
     verbose = opts.verbose,
     percent_match = float(opts.percent_match),
     full_primer_length = int(opts.full_primer_length),
     output_f = opts.output_filepath,
     specificity_threshold = float(opts.specificity_threshold),
     log_filepath = opts.log_file,
     standard_index_file = opts.standard_index_file,
     search_range = opts.search_range)


               
if __name__ == "__main__":
    main()


