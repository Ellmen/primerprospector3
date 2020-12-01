#!/usr/bin/env python


__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

from os.path import isfile

from cogent.util.misc import parse_command_line_parameters, create_dir
from optparse import make_option

from primerprospector.check_primer_barcode_dimers import \
 check_barcodes_and_primers

script_info={}
script_info['brief_description']=""" Tests input combinations of barcodes and primers for potential secondary structure."""
script_info['script_description']="""

This module takes an input text file containing barcodes, line separated, and
listed in 5' to 3' direction.  Two primers are specified at the command line.
The first, primer 1, is the primer sequence (and linker) that is associated
with a barcode.  Primer 2 by default has no barcode associated with it, 
although barcode association with both primers can be enabled.

Each barcode/primer combination that potentially will form secondary 
structures will be listed with line number, barcode sequence, secondary 
structure, and energy in an output text file.  Additionally, a postscript file
showing graphical output for the secondary structure will be generated in a 
subdirectory of the output directory.  Degenerate primers will be listed in 
the form(s) that are likely to have secondary structure.
"""
script_info['script_usage']=[]


script_info['script_usage'].append(("""Example""","""Standard Example:""","""%prog [options] {-b barcode_filepath [required] -p primer_sequence_1 [required] -P primer_sequence_2 [required] -e DNA_energy_parameters_filepath [required]}"""))
script_info['script_usage'].append(("""Set the annealing temperature to 65, set score threshold to -10, and point to the energy parameters file in the directory /home/BigErn/pprospector/DNA_parameters:""","""""","""check_primer_barcode_dimers.py -b barcodes.txt -p "ACCTGACRGGTAATC" -P "CGACCAGTACG" -t 65 -s -10 -e /home/BigErn/pprospector/DNA_parameters/dna_DM.par"""))

script_info['output_description']=("""Barcode/primer combinations that are likely to form secondary structures will be listed in the output 'flagged_barcodes.txt' file.  Any such barcode/primer combination will also be present in the graphs/ directory, and named for the line number that the barcode was found on and the barcode sequence.""")

script_info['required_options']=[\
    # Barcodes file
    make_option('-b', '--barcodes',
                help='Filepath of barcodes to score input primer(s) against.'),
    
    make_option('-p', '--primer1',
                help="Primer, written in 5' to 3', that is linked to barcodes"+\
                " tested.  If linker sequence is present between primer and "+\
                "barcode, include it with this sequences."),
                
    make_option('-P', '--primer2',
                help="Second primer, written in 5' to 3' orientation.  This "+\
                "primer by default is not associated with any barcodes."),
                
    make_option('-e', '--energy_parameters',
        help='Specify energy parameters file for predicting secondary '+\
        'structures.  A DNA parameters file, dna_DM.par, is found in the '+\
        'DNA_parameters folder of Primer Prospector, and should be pointed '+\
        'to with this parameter.  If an incorrect file is used, the '+\
        'Vienna software will use default parameters, which are for RNA '+\
        'folding, and could give misleading results.  The provided DNA '+\
        "parameters file is a modified form of the DNA parameters from "+\
        " David Mathews' RNAstructure program.")]
                    
script_info['optional_options']=[\
        
    make_option('-t', '--annealing_temp', 
                help='Specify an annealing temperature in degrees Celsius.'+\
                ' [default: %default]', default=50),
        
    make_option('-s', '--score_threshold',
        help='Specify a score threshold for the Gibbs energy calculation, '+\
        'below which a barcode/primer combination is flagged for potential '+\
        'secondary structure.  [default: %default]', default=-10.0),
        
    # Specify output directory
    make_option('-o', '--output_dir',
        help='Specify output directory for barcode/primer secondary structure'+\
        ' summary and graphs. [default: %default]', default="."),
        
    
    make_option('-B', '--paired_end_barcodes',
        help='If enabled, barcodes will be appended to both primer 1 and '+\
        'primer 2.  [default: %default]', default = False, 
        action = 'store_true'),
        
    make_option('-g', '--suppress_graphs',
        help='Suppress retention of output postscript graphs. '+\
        '[default: %default]', default = False, action = 'store_true')
        
    
]



script_info['version'] = __version__

def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)
    

    barcodes_fp  = opts.barcodes
    primer1 = opts.primer1
    primer2 = opts.primer2
    output_dir = opts.output_dir
    annealing_temp = float(opts.annealing_temp)
    score_threshold = float(opts.score_threshold)
    paired_end_barcodes = opts.paired_end_barcodes
    energy_parameters = opts.energy_parameters
    suppress_graphs = opts.suppress_graphs
    
    if not isfile(energy_parameters):
        raise ValueError,('Specified path %s ' % energy_parameters +\
         'for DNA energy parameters not a file.')

    
    # Create output directory if it does not exist
    create_dir(output_dir)
    
    check_barcodes_and_primers(barcodes_fp, primer1, primer2, energy_parameters,
     output_dir, annealing_temp, score_threshold, paired_end_barcodes,
     suppress_graphs)

    


if __name__ == "__main__":
    main() 
