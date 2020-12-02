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


from primerprospector.cogentutil.misc import parse_command_line_parameters, create_dir
from optparse import make_option

from primerprospector.taxa_assignment_report import generate_taxa_report



script_info = {}
script_info['brief_description'] = """ Assign taxonomies for input fasta sequences, generate report on accuracy of assignment by comparing to known taxonomies """
script_info['script_description'] = """

This module will attempt to assign taxonomies to the input fasta files, then 
compare these results to the known taxonomies supplied for each sequence ID 
in the input taxa mapping file.  Currently only the RDP classifier is 
implemented for assigning taxonomies.

The primary purpose of this module is to test amplicons and/or reads for a
prospective primer pair (see get_amplicons_and_reads.py) for their 
phylogenetic usefulness.  Very short reads, or reads in a region of high
conservation will be less accurate in terms of taxonomic assignment.

This module will save the taxonomic assignments, as well as a report detailing
the percentage that were accurately assigned for each level of taxonomic depth.
An assignment of Archaea;Euryarchaeota;Thermoplasmatales for a sequence that
was actually Archaea;Euryarchaeota;Halobacteriales is accurate to the second
level of taxonomy, but not the third, and will be recorded as such.  The 
default depth of taxa to test is 3.
"""

script_info['script_usage']=[]
script_info['script_usage'].append(("""Standard Example usage:""","""""","""%prog [options] {-t taxa_mapping_filepath [required] -f input_fasta_filepath [required]}"""))
script_info['script_usage'].append(("""Change taxa depth to 4, change output directory to taxa_report/""","""""","""taxa_assignment_report.py -t taxa_mapping_filepath -f input_fasta_filepath  -d 4 -o taxa_report/ """))

script_info['output_description']=""" A log file from the taxonomy classifier will be generated, along with a taxonomic assignment file.  Additionally, a taxa_assignment_report.txt file will be generated detailing the accuracy of assignment for each level of taxonomic depth tested."""


script_info['version'] = __version__

script_info['required_options']=[\

    make_option('-t', '--taxa_mapping_fp',
        help='Taxonomy mapping filepath'),
                    
    make_option('-f', '--fasta_fp',
        help='Fasta sequence file.')]
                 
script_info['optional_options']=[\

    make_option('-d', '--taxa_depth', 
        help='Depth of taxonomy to test for accuracy.  Depth that exceeds '+\
        'specifications in the taxonomy mapping file or report will be '+\
        'ignored [default: %default]', default = 3),
        
    make_option('-o', '--output_dir',
        help='Specify output directory for reports, log. '+\
        '[default: %default]', default = "."),
        
    make_option('-m', '--assignment_method',
        help=' Taxonomic assignment method.  Currently only RDP classifier '+\
        'implemented. [default: %default]', default='rdp'),
        
    make_option('-c', '--min_confidence',
        help='Minimum confidence for taxonomic assignment.  '+\
        '[default: %default]', default =0.80),
        
    make_option('-T', '--training_data_fp',
        help='Training sequence data filepath for rdp classifier. '+\
        '[default: %default]', default = None)
]

def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    taxa_mapping_fp = opts.taxa_mapping_fp
    fasta_fp = opts.fasta_fp
    taxa_depth = int(opts.taxa_depth)
    output_dir = opts.output_dir
    assignment_method = opts.assignment_method
    min_confidence = float(opts.min_confidence)
    training_data_fp = opts.training_data_fp

    # Check assignment method
    valid_assignment_methods = ['rdp']
    
    if assignment_method not in valid_assignment_methods:
        raise ValueError('Assignment method must be one of the following: '+\
         '%s' % ",".join(valid_assignment_methods))
         
    
    # create output directory if not present
    create_dir(output_dir)
    
    generate_taxa_report(taxa_mapping_fp, fasta_fp, taxa_depth,
     output_dir, assignment_method, min_confidence, training_data_fp)
    
    
if __name__ == "__main__":
    main()


