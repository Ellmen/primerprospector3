#!/usr/bin/env python
# File created on 26 Apr 2010
from __future__ import division

__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters", "Greg Caporaso", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"
 
""" Contains file formatting, IO for PrimerProspector modules """

from os.path import basename

def write_formatted_primers(f_primers_sorted, 
                            r_primers_sorted,
                            formatted_primers):
    """ Write sorted primers in proper format for data analysis pipeline 
    
    f_primers_sorted: list of forward primers, sorted according to specified
     method
    r_primers_sorted: list of reverse primers, sorted.
    formatted_primers: open output file object."""
    
    
    formatted_primers.write("# primer_id <tab> primer sequence (5'->3')\n")
    for p in f_primers_sorted:
        formatted_primers.write(p)
    for p in r_primers_sorted:
        formatted_primers.write(p)
    return
    
def organize_overlapping_primers(primers):
    """ Sort primers according to overlap state 
    
    primers: list of Primer objects"""
    
    unique_primers = []
    overlapping_primers = []
    three_prime_match = []
    # two passes, to build up f primer information first, then r primer
    for p in primers:
        if p.f_unique:
            unique_primers.append((p.primer_f_name + '\t' +\
             p.f_degenerate_seq + '\n'))
        else:
            if p.f_partial_overlap:
                for n in p.f_partial_overlap:
                    overlapping_primers.append(n)
            if p.f_3prime_match:
                for n in p.f_3prime_match:
                    three_prime_match.append(n)
    for p in primers:
        if p.r_unique:
            unique_primers.append((p.primer_r_name + '\t' +\
             p.r_degenerate_seq + '\n'))
        else:
            if p.r_partial_overlap:
                for n in p.r_partial_overlap:
                    overlapping_primers.append(n)
            if p.r_3prime_match:
                for n in p.r_3prime_match:
                    three_prime_match.append(n)
                    
    return unique_primers, overlapping_primers, three_prime_match
        

        
def write_overlap_primers(primers,
                          primers_overlap):
    """ Write report about primers overlapping with known primers 
    
    primers: list of Primers objects
    primers_overlap: open output file object for overlapping primers data."""
    

    unique_primers, overlapping_primers, three_prime_match =\
     organize_overlapping_primers(primers)
    
                    
    primers_overlap.write("#Unique primers that do not overlap or share 3' "+\
     "ends with known primers\n")
    primers_overlap.write("#primer name<tab>primer sequence\n")
    for n in unique_primers:
        primers_overlap.write(n)
    primers_overlap.write("#Primers overlapping with known primers.\n")
    primers_overlap.write("#primer name<tab>primer sequence<tab>known primer "+\
     "name<tab>known primer sequence<tab>overlapping sequence\n")
    for n in overlapping_primers:
        primers_overlap.write(n)
    primers_overlap.write("#Primers sharing 3' ends with known primers.\n")
    primers_overlap.write("#primer name<tab>primer sequence<tab>known primer "+\
     "name<tab>known primer sequence\n")
    for n in three_prime_match:
        primers_overlap.write(n)
        
    return


def write_primers_summary(primers,
                          primer_report,
                          variable_pos_freq):
    """ Write summary primers data in comma delineated form 
    
    primers: list of Primer objects
    primer_report: open file object to write report with
    variable_pos_freq: Allowed degeneracy for primer."""
    
    # Write details about data
    primer_report.write(primers_summary_info)
    
    for p in primers:
        header = p.header
        # Check for standard alignment, add empty fields if missing
        if len(header.split(",")) == 6:
            header += ",,"
        primer_data = ','.join([header, 
                                p.f_IUPAC, 
                                p.f_degenerate_seq,
                                p.f_consensus, 
                                p.r_IUPAC, 
                                p.r_degenerate_seq, 
                                p.r_consensus])
        primer_report.write(primer_data)
        primer_report.write("\n")
        
    return


def write_denovo_output_file(output_data,
                             output_filepath,
                             region_slice):
    """ Writes data to output files 
    
    output_data: list of strings of denovo primer data
    output_filepath: open output filepath
    region_slice: slice of sequence up and downstream from conserved site to
     retain
    """
    

    for n in output_data:
        # Build string, data depends on if standard alignment option enabled.
        alignment_line=",".join(("C", n.seq, str(n.unaligned_index),
         str(n.aligned_index), str(n.match_count), str(n.percent_match*100)+"%",
         str(n.non_specific_percent*100)+"%"))
        # If present, find the standard index for the prospective primer
        if n.std_index:
            # format additional alignment data, preceeded by comma
            alignment_line += "," + ",".join((str(n.f_std_index), 
             str(n.r_std_index))) + "\n"
        else:
            alignment_line += "\n"
        
        output_filepath.write(alignment_line)

        label_list=n.labels
        upstream_primer=n.upstream_regions
        downstream_primer=n.downstream_regions
        for m in range(len(label_list)):
            output_filepath.write(",".join(("M", str(n.aligned_index),
             upstream_primer[m], downstream_primer[m], label_list[m] )))
            output_filepath.write("\n")
             


def generate_denovo_output_file(primers, 
                                output_filepath,
                                specificity_threshold,
                                region_slice, 
                                standard_index_seq,
                                percent_match, 
                                specificity_test):
    """ Generates text output files 
    
    primers: list of Primer objects
    output_filepath: open output file object
    specificity_threshold: Max allowed non-specific sequences
    region_slice: size, in base pairs, of up and downstream sequences to 
     retain along with the short conserved site acting as the 3' primer end.
    standard_index_seq: Optional standard alignment (e.g. E. coli) to give
     de novo primers indices that approximate indices of known primers.
    percent_match: sensitivity threshold
    specificity_test: Set to True if specificity requirement when generating
     denovo primers.
    """

    title = "# Matches that are found in at least %1.2f%s target sequences" %\
     (percent_match*100, "%")
    if specificity_test:
        title +=" and specific at or below the %1.2f%s threshold \n" %\
         (specificity_threshold*100, "%")
    else:
        title +="\n"
    output_filepath.write(title)
    primer_info_line = "# record type,sequence,unaligned index,aligned index"+\
        ",number of perfect matches,percent perfect matches,percent "+\
        "non-specific matches"
        
    if standard_index_seq:
        primer_info_line += ",forward primer standard index,reverse primer "+\
         "standard index\n"
    else:
        primer_info_line += "\n"
            
    output_filepath.write(primer_info_line)
    
    write_denovo_output_file(primers,output_filepath,\
    region_slice)
    
def get_amplicon_pairs(amplicon_len, 
                       f_primers_sorted,
                       r_primers_sorted):
    """ Generates formatted data for primer pairs in given amp size range 
    
    amplicon_len:  min, max size of amplicon in format X:Y
    f_primers_sorted: list of forward primer names<tab>sequences
    r_primers_sorted: list of reverse primer names<tab>sequences"""
    
    # Get min and max amplicon sizes
    min_amp_size = int(amplicon_len.split(":")[0])
    max_amp_size = int(amplicon_len.split(":")[1])
    
    amplicon_data = []
    
    amplicon_data.append("# Primer pairs within estimated amplicon size\n")
    amplicon_data.append("# Min size: %d\n" % min_amp_size)
    amplicon_data.append("# Max size: %d\n" % max_amp_size)
    
    found_pair = False
    
    for forward_primer in f_primers_sorted:
        
        # Get 3' index of primer by stripping all alphabet characters
        f_primer_index = forward_primer.split("\t")[0]
        f_primer_index = f_primer_index.strip(uppercase)
        f_primer_index = int(f_primer_index.strip(lowercase))
        
        # Also need to add the length of the primer to get the 3' end site
        f_primer_index += len(forward_primer.split("\t")[1].strip())
        
        for reverse_primer in r_primers_sorted:
            
            r_primer_index = reverse_primer.split("\t")[0]
            r_primer_index = r_primer_index.strip(uppercase)
            r_primer_index = int(r_primer_index.strip(lowercase))
            
            r_primer_index -= len(reverse_primer.split("\t")[1].strip())
            
            est_amp_size = r_primer_index - f_primer_index
            
            if est_amp_size in range(min_amp_size, max_amp_size):
                
                found_pair = True
                
                amplicon_data.append("\n")
                amplicon_data.append(forward_primer)
                amplicon_data.append(reverse_primer)
                amplicon_data.append("Estimated amplicon size: %d\n" %\
                 est_amp_size)
                 
    if not found_pair:
        amplicon_data.append("No primer pairs were found that match amplicon "+\
         "size criteria.")
         
    return amplicon_data
    
    
def write_amplicon_pairs(amplicon_data,
                         amplicon_pairs_f):
    """ Writes primer pair data that fits amplicon size criteria 
    
    amplicon_data:  list of lines containing amplicon pairs to be written
    amplicon_pairs_f:  open file object to write amplicon_data to"""
    
    for line in amplicon_data:
        amplicon_pairs_f.write(line)
        
        
def format_base_freq_data(base_freqs,
                          primer_seq,
                          hits_fp,
                          score_type,
                          score_threshold):
    """ Structures base frequency data into string for writing to file
    
    base_freqs: list of dicts containing base frequency data for each position
     in the sequence hits.
    primer_seq:  primer sequence.
    hits_fp: filepath for primer hits file tested.  Used to log the file name
     of the hits file used to generate base frequency data.
    score_type: The type of scoring used to determine if a sequence should
     be considered for base frequency data.  Can be weighted_score, 
     tp_mismatches, or overall_mismatches
    score_threshold: Value at which or below a sequence score is considered
     acceptable (i.e., that sequence is likely to amplify during PCR).
    """
    
    base_freq_data = base_freq_header
    
    base_freq_data += "# Hits file used to generate base frequency data: %s\n"%\
     basename(hits_fp)
    base_freq_data += "# Score type used: %s\n" % score_type
    base_freq_data += "# Score threshold: %d\n" % score_threshold
    base_freq_data += "# Primer sequence: %s\n" % primer_seq
    base_freq_data += "#\n"
    base_freq_data += "Primer"
    
    for base in primer_seq:
        base_freq_data += "\t%s" % base
    
    base_freq_data += "\n"
    
    base_freq_data += "Base\n"
    
    bases = ['A', 'T', 'C', 'G']
    
    for curr_base in bases:
        
        base_freq_data += curr_base
        
        for base_freq in base_freqs:
            
            base_freq_data += "\t%1.2f" % base_freq[curr_base]
        
        base_freq_data += "\n"
        
    return base_freq_data
    
        
        
    
    
def format_linker_data(linker_summary,
                       best_seq,
                       worst_seq,
                       hits,
                       score_type,
                       score_threshold):
    """ Structures suggested linker data into string for output
    
    linker_summary: Percentages of nucleotides appearing at each base position.
    best_seq: sequence of nucleotides that least frequently appear
    worse_seq: sequence of nucleotides that most often appear
    hits: filepath for primer hits file tested.  Used to log the file name
     of the hits file used to generate linker data.
    score_type: The type of scoring used to determine if a sequence should
     be considered for making linkers.  Can be weighted_score, tp_mismatches,
     or overall_mismatches
    score_threshold: Value at which or below a sequence score is considered
     acceptable to be considered for making a linker (i.e., that sequence is
     likely to amplify during PCR).
    """
    
    summary_data = linker_summary_notes
    
    summary_data += "Hits file used to generate linkers: %s\n" % basename(hits)
    summary_data += "Score type used: %s\n" % score_type
    summary_data += "Score threshold: %d\n" % score_threshold
    summary_data += "\n"
    
    summary_data += linker_summary
    
    summary_data += "\nSuggested linker sequence:\n%s\n" % best_seq
    summary_data += "Worst linker sequence:\n%s" % worst_seq
    
    return summary_data
    
def write_assigned_taxa(assigned_taxa,
                        taxa_assignment_outf):
    """ Writes taxonomy assignments to the filepath taxa_assignment_outf
    
    assigned_taxa: dict of sequence ID: (assigned taxonomy, confidence)
    taxa_assignment_outf: taxa assignment output file object
    """

    
    for assignment in assigned_taxa:
        line = assignment + "\t" + assigned_taxa[assignment][0] + "\t%1.3f\n" %\
         assigned_taxa[assignment][1]
        taxa_assignment_outf.write(line)
        
def write_accuracy_report(accuracy_values,
                          seqs_lacking_taxa_mapping,
                          seqs_assigned_root,
                          report_outf,
                          fasta_fp,
                          assignment_method,
                          training_data_fp):
    """ Writes report about accuracy of assignments, number of unmapped seqs
    
    accuracy_values:  list of tuples of (percent accurate assignment, number
     of seqs for taxa level)
    seqs_lacking_taxa_mapping: Total number of sequences where the sequence
     ID could not be found in the mapping file.
    seqs_assigned_root: Sequences assigned only as 'Root'.  These may or may
     not be real assignments, quantified and logged.
    report_outf: Accuracy report file object
    fasta_fp: filepath of fasta sequences which were assigned to taxonomies.
     Used to log filename in report file.
    assignment_method: method used for assignment.  Used to log in report file.
    training_data_fp: Training data for RDP classifier, if used, logged in
     report file.
    """
    
    
    report = report_header
    
    report += "# Fasta file used for taxonomic assignments: %s\n" %\
     basename(fasta_fp)
    report += "# Assignment method: %s\n" % assignment_method
    if training_data_fp:
        report += "# Training data filepath for RDP classifier: %s\n" %\
         basename(training_data_fp)
         
    report += "# Start report accuracy data\n"
    report += "# Taxa level, percent accurate assignment, number of sequences"+\
     " with taxa defined at this level\n"
    
    
    count = 0
    for taxa_level in accuracy_values:
        report += "%d,%2.3f,%d\n" % (count, taxa_level[0], taxa_level[1])
        count += 1
         
    report +=\
     "Sequences lacking corresponding ID in taxonomy mapping file: %d\n" %\
     seqs_lacking_taxa_mapping
     
    report +=\
     "Sequences assigned only as 'Root': %d\n" % seqs_assigned_root
         
    report_outf.write(report)
    
    
    
linker_summary_notes="""# Summary data for suggested linkers
# Note-degenerate bases are ignored for linker results, as are linkers that
# exceed the length of the input fasta sequence.  Base position starts at the
# 5' position of the linker, so position 0 would be the 5' most base in the
# linker, position 1 would be the next base, and so on.
"""
           

primers_summary_info = """# Primer Report
# This file contains details about the prospective primers given in the formatted_primers.txt file.
# Conserved_Xmer is the conserved site found with the sequence_searcher module that serves as the primer 3' region.
# Unaligned_Index is the first position the X-mer was found when running the sequence_searcher module after degapping.
# Aligned_Index is the first position the X-mer was found when running the sequence_searcher module.
# Number_Hits is the number of sequences that had a perfect match to the X-mer.
# Percent_Match is the percentage of total sequences with a perfect match.
# Nonspecific_Match is the percentage of sequences in exclusion files that had perfect matches.  0.0% if no files were specified.
# 5'_Standard_Index gives a corrected index based on an alignment given with the -a parameter with the sequence_searcher module.  Will be empty if no standard alignment file was specified.
# 3'_Standard_Index gives a corrected index based on an alignment given with the -a parameter with the sequence_searcher module.
# F_Primer_IUPAC_Sequence shows the fully degenerate version of the prospective forward primer.  If any gaps are present, this base will be given as a '.' character.
# F_Primer_Filtered_Sequence shows a filtered version of the degenerate sequence where bases have to be present more than the specified frequency (see -p parameter) to contribute to degeneracy.
# F_Primer_Majority_Consensus gives the most frequently occurring bases.
# R_Primer data are the same as the F_Primer, but for the reverse primer, and have been reverse complemented.
# Start primer data:
# Conserved_Xmer, Unaligned_Index, Aligned_Index, Number_Hits, Percent_Match, Nonspecific_Match, 5'_Standard_Index, 3'_Standard_Index, F_Primer_IUPAC_Sequence, F_Primer_Filtered_Sequence, F_Primer_Majority_Consensus, R_Primer_IUPAC_Sequence, R_Primer_Filtered_Sequence, R_Primer_Majority_Consensus
"""


report_header = """# Taxonomy Assignment Report File
# This report starts at the highest level of taxa (i.e., Domain level)
# and lists the accuracy of the assignment and the number of sequences
# that had a taxonomic assignment and taxonomic mapping to that depth
"""

base_freq_header = """# Base frequency report for optimizing primers
# This file is tab separated for easy importation into Excel or other spreadsheets
# Degenerate DNA codes (listed here for convenience): R=AG, Y=CT, M=AC, K=GT, W=AT, S=CG, B=CGT, D=AGT, H=ACT, V=ACG, N=ACGT
# The primer is listed in 5' to 3' orientation.
"""

