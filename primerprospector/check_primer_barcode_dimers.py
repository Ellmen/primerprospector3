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

""" This module contains functions for finding the folding energy of 
barcode-primer/primer combinations, creating a text file of these combinations
whose energy falls below a particular threshold, and generation of graphical
files displaying the secondary structure of these flagged sequences"""

from os import rename, remove
from os.path import isfile
upper = str.upper

from cogent3.app.vienna_package import RNAfold
from cogent3.seqsim.sequence_generators import SequenceGenerator, IUPAC_DNA
from cogent3 import DNA

from primerprospector.parse import get_barcodes

def expand_degeneracies(raw_primer):
    """ Returns all non-degenerate versions of a given primer sequence """
    
    primers=SequenceGenerator(template=raw_primer,alphabet=IUPAC_DNA)
    expanded_primers=[]
    for primer in primers:
        expanded_primers.append(primer)
        
    return expanded_primers
    
    
def get_energy_res(text_out):
    """ Pulls Gibbs energy result from text output of RNAfold
    
    text_out: String of text output from RNAfold, combined into single line
    """
    
    text_out = text_out.replace(')', '')
    text_out = text_out.replace('(', '')
    text_out = text_out.strip()
    
    energy = float(text_out.split(' ')[-1])
    
    return energy
    

def check_barcodes_and_primers(barcodes_fp,
                               primer1,
                               primer2,
                               energy_parameters,
                               output_dir = '.',
                               annealing_temp = 55.0,
                               score_threshold = -2.0,
                               paired_end_barcodes = False,
                               suppress_graphs = False):
    """ Main function for checking barcodes/primers for secondary structure.
    
    barcodes_fp: filepath for barcodes.  This file has all barcodes to be used,
     line separated, and listed in 5' to 3' direction.
    primer1: primer associated with barcodes (can be forward or reverse primer).
     This primer is listed in 5' to 3' direction.
    primer2: second primer.  By default, is not associated with barcodes.
    energy_parameters: filepath containing DNA parameters file.
    output_dir: output directory.
    annealing_temp: Annealing temperature-this should be the lowest 
     temperature during the PCR process.  In degrees Celcius.
    score_threshold: Value at or below which barcode/primer combinations will
     be flagged for potential secondary structure.  Gibbs energy in kcal/mol.
    paired_end_barcodes: If True, primer2 will also have barcode appended to 
     its 5' end for dimer/secondary structure calculations.
    """
    
    primer1 = upper(primer1)
    primer2 = upper(primer2)
    
    f_primers = expand_degeneracies(primer1)
    r_primers = expand_degeneracies(primer2)
    
    connector_seq = "----------"
    
    barcodes_f = open(barcodes_fp, "U")
    barcodes = get_barcodes(barcodes_f)
    
    # Save a list of good barcodes to create filtered output file
    good_barcodes = []
    
    text_output = [text_output_header]
    barcodes_flagged = False
    
    if not output_dir.endswith("/"):
        output_dir += "/"

        
    # Single test for reverse primer dimers if not using paired end barcodes
    # Second test-Reverse primer to reverse primer pairing
                    
    if not paired_end_barcodes:
        
        for r_primer in r_primers:
            
            curr_seq = [r_primer + connector_seq + r_primer]

            r = RNAfold(WorkingDir = output_dir)
            r.Parameters['-T'].on(annealing_temp)
            r.Parameters['-P'].on(energy_parameters)
           
            result = r(curr_seq)
            
            # Create text output line containing all sequences and index
            # of the barcode (line number) it was found in the input 
            # barcode file.  Also, replace "U" with "T" for output so there
            # is no confusion about sequences.
            
            raw_text = ",".join([l.strip().replace("U", "T") \
             for l in result['StdOut'].readlines()])
            
            energy_res = get_energy_res(raw_text)
            
            text_out = "%s," % ",,Primer2,Primer2," +\
             raw_text.split(' ')[0] + ",%2.2f\n" % energy_res
            
            if energy_res <= score_threshold:
                
                barcodes_flagged = True
                good_bc = False
                text_output.append(text_out)
                
                ps_outf_name =\
                 output_dir + "primer2V%s_primer2V%s.ps" %\
                 (r_primers.index(r_primer), r_primers.index(r_primer))
                 
                if suppress_graphs:
                    if isfile(output_dir + "rna.ps"):
                        remove(output_dir + "rna.ps")
                else:
                    rename(output_dir + "rna.ps", ps_outf_name)
    
    # Test the forward and reverse primers against each other, as well as the
    # forward against forward, reverse against reverse.
    for curr_bc in barcodes:
        
        good_bc = True

        for f_primer in f_primers:

            # Test forward primer and itself
            curr_seq = [curr_bc + f_primer + connector_seq + curr_bc + f_primer]
            
            r = RNAfold(WorkingDir = output_dir)
            r.Parameters['-T'].on(annealing_temp)
            r.Parameters['-P'].on(energy_parameters)
           
            result = r(curr_seq)
            
            # Create text output line containing all sequences and index
            # of the barcode (line number) it was found in the input 
            # barcode file.  Also, replace "U" with "T" for output so there
            # is no confusion about sequences.
            
            raw_text = ",".join([l.strip().replace("U", "T") \
             for l in result['StdOut'].readlines()])
            
            energy_res = get_energy_res(raw_text)
            
            text_out = "%s," % barcodes.index(curr_bc) + curr_bc +\
             ",Primer1,Primer1," +\
             raw_text.split(' ')[0] + ",%2.2f\n" % energy_res
            
            if energy_res <= score_threshold:
                
                barcodes_flagged = True
                good_bc = False
                text_output.append(text_out)
                
                ps_outf_name =\
                 output_dir + "Line%s_%s_primer1V%s_primer1V%s.ps" %\
                 (barcodes.index(curr_bc), curr_bc,  
                 f_primers.index(f_primer), f_primers.index(f_primer))
                 
                if suppress_graphs:
                    if isfile(output_dir + "rna.ps"):
                        remove(output_dir + "rna.ps")
                else:
                    rename(output_dir + "rna.ps", ps_outf_name)
                
            # Remove default postscript file if present
            if isfile(output_dir + "rna.ps"):
                remove(output_dir + "rna.ps")
            
            for r_primer in r_primers:
                
                if paired_end_barcodes:
                    curr_seq = [curr_bc + f_primer + connector_seq +\
                     curr_bc + r_primer]
                else:
                    curr_seq = [curr_bc + f_primer + connector_seq + r_primer]
                
                r = RNAfold(WorkingDir = output_dir)
                r.Parameters['-T'].on(annealing_temp)
                r.Parameters['-P'].on(energy_parameters)
               
                result = r(curr_seq)
                
                # Create text output line containing all sequences and index
                # of the barcode (line number) it was found in the input 
                # barcode file.  Also, replace "U" with "T" for output so there
                # is no confusion about sequences.
                
                raw_text = ",".join([l.strip().replace("U", "T") \
                 for l in result['StdOut'].readlines()])
                
                energy_res = get_energy_res(raw_text)
                
                text_out = "%s," % barcodes.index(curr_bc) + curr_bc +\
                 ",Primer1,Primer2," +\
                 raw_text.split(' ')[0] + ",%2.2f\n" % energy_res
                
                if energy_res <= score_threshold:
                    
                    barcodes_flagged = True
                    good_bc = False
                    text_output.append(text_out)
                    
                    ps_outf_name =\
                     output_dir + "Line%s_%s_primer1V%s_primer2V%s.ps" %\
                     (barcodes.index(curr_bc), curr_bc,  
                     f_primers.index(f_primer), r_primers.index(r_primer))
                     
                    if suppress_graphs:
                        if isfile(output_dir + "rna.ps"):
                            remove(output_dir + "rna.ps")
                    else:
                        rename(output_dir + "rna.ps", ps_outf_name)
                        
                    # Second test-Reverse primer to reverse primer pairing
                    
                    if paired_end_barcodes:
                        curr_seq = [curr_bc + r_primer + connector_seq +\
                         curr_bc + r_primer]
                    else:
                        # Do single test for reverse primers before to save
                        # computation time
                        continue
                    
                    r = RNAfold(WorkingDir = output_dir)
                    r.Parameters['-T'].on(annealing_temp)
                    r.Parameters['-P'].on(energy_parameters)
                   
                    result = r(curr_seq)
                    
                    # Create text output line containing all sequences and index
                    # of the barcode (line number) it was found in the input 
                    # barcode file.  Also, replace "U" with "T" for output so there
                    # is no confusion about sequences.
                    
                    raw_text = ",".join([l.strip().replace("U", "T") \
                     for l in result['StdOut'].readlines()])
                    
                    energy_res = get_energy_res(raw_text)
                    
                    text_out = "%s," % barcodes.index(curr_bc) + curr_bc +\
                     ",Primer2,Primer2," +\
                     raw_text.split(' ')[0] + ",%2.2f\n" % energy_res
                    
                    if energy_res <= score_threshold:
                        
                        barcodes_flagged = True
                        good_bc = False
                        text_output.append(text_out)
                        
                        ps_outf_name =\
                         output_dir + "Line%s_%s_primer2V%s_primer2V%s.ps" %\
                         (barcodes.index(curr_bc), curr_bc,  
                         f_primers.index(f_primer), r_primers.index(r_primer))
                         
                        if suppress_graphs:
                            if isfile(output_dir + "rna.ps"):
                                remove(output_dir + "rna.ps")
                        else:
                            rename(output_dir + "rna.ps", ps_outf_name)
                    
                # Remove default postscript file if present
                if isfile(output_dir + "rna.ps"):
                    remove(output_dir + "rna.ps")
        
        if good_bc:
            good_barcodes.append(curr_bc)
    
    if not barcodes_flagged:
        text_output.append('No barcodes/primer combination was found below '+\
         'the specified energy threshold.')
    
    text_f = open(output_dir + "barcode_results.txt", "w")
    
    for line in text_output:
        text_f.write(line)
        
    filtered_bc_f = open(output_dir + "filtered_barcodes.txt", "w")
    
    if len(good_barcodes):
        filtered_bc_f.write("# The following barcodes were not flagged for "+\
         "any potential secondary structure with the given primers.\n")
        for line in good_barcodes:
            filtered_bc_f.write(line + "\n")
    else:
        filtered_bc_f.write("# All barcodes were flagged for secondary "+\
         "structure with the given primers.")
        
    
    

text_output_header = """# Barcode/Primer combinations that fall below the energy threshold.
# primer1 and primer2 will be listed in all non-degenerate forms that result in energy values below threshold.
# If primer2 is not barcoded, it will be tested against itself once, and listed first with the barcode index and sequence fields left empty.
# line number of barcode from input barcode file, barcode sequence, primer1 identity, primer2 identity, combined sequence, secondary structure, Gibbs energy in kcal/mol
"""




