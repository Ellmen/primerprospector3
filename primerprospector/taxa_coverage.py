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


""" taxa_coverage.py is intended to return phylogenetic information about
    the taxa expected to be amplified by a given primer or primer pair.
    (see analyze_primers.py for information regarding primer scoring).  
    This module will generate output text files giving a summary of the 
    percent expected hits in each level of taxa specified (default 3).  The
    depth of taxa that can be covered in this module depends on the level
    specified in the taxa mapping file that is supplied.
    
    The scoring methods can be weighted score (which weights mismatches with
    a harsher penalty that are near the 3' end of the primer), the overall
    mismatches, or three prime mismatches only.
    
    This module can either analyze a single hits file, multiple hits files
    separated with a colon, or all files ending in "_hits.txt" for a given 
    input directory.  Output files will be based upon the given _hits.txt 
    file name(s), and will be given the "_taxa_coverage.txt" or 
    "_taxa_coverage_level_X.pdf" extension.
    
    """


from os.path import basename
capitalize = str.capitalize
lower = str.lower
from math import ceil

from numpy import arange
from matplotlib import use
use('Agg')
from matplotlib.pyplot import ylim
from pylab import figure, savefig, figtext, yticks, xticks
from matplotlib import colors
from matplotlib.pyplot import clf, close, subplots_adjust
from cogent.util.misc import create_dir

from primerprospector.parse import get_primer_hits_data, get_hits_files, \
 parse_hits_data, parse_taxa_mapping_file, get_primer_hits_data_pair
from primerprospector.util import get_primer_direction

def get_coverage_data_primer_pair(pair_hits_data,
                                  taxa_map,
                                  taxa_depth,
                                  score_type,
                                  score_threshold,
                                  hits_f_fp,
                                  hits_r_fp,
                                  taxa_fp):
    """ Sorts coverage data into categories, total seqs, seqs passing score
    
    For primer pairs, the worst scoring primer of the two decides whether a 
    particular sequences is considered passing.
    
    
    pair_hits_data: 2D list of lines of hit data from forward and reverse
     hits file
    taxa_map: dict of sequence ID: taxonomy (semicolon separated, highest level 
     of taxa to lowest level)
    taxa_depth: level of taxa to attempt to categories.  This function will
     stop trying to add new levels to the returned list when the current list
     completes and is empty.
    score_type: Type of score from hits file to utilize.  Valid choices are
     weighted_score, overall_mismatches, and tp_mismatches.
    score_threshold: Value at which or below a particular sequence will be 
     considered to be successfully amplified during PCR.
    hits_f_fp: filepath of forward hits file tested.  Used for log file.
    hits_r_fp: filepath of reverse hits file tested.  Used for log file.
    taxa_fp: filepath of the taxonomy mapping file.  Used for log file.

    """
    
    coverage_data = []
    # log data will contain score type, threshold, total seqs analyzed, total
    # seqs lacking corresponding IDs in the mapping file, and the IDs of these
    # seqs not found in the mapping file.
    log_data = []
    
    f_hits_index = 0
    r_hits_index = 1

    
    # Add data to beginning of log file
    log_data.append(log_file_header)
    log_data.append("# Hits files analyzed: %s,%s\n" %\
     (basename(hits_f_fp), basename(hits_r_fp)))
    log_data.append("# Taxonomy mapping file used: %s\n" % basename(taxa_fp))
    log_data.append("# Scoring method used: %s\n" % score_type)
    log_data.append("# Score threshold: %2.2f\n" % score_threshold)
    
    total_seqs = len(pair_hits_data[f_hits_index])
    
    log_data.append("# Total sequences in hits files: %d\n" % total_seqs)
    
    total_seqs_no_taxa = 0
    seqs_no_tax = ""
    
    domain_index = 0
    
    totals_index = 0
    passing_index = 1
    
    no_taxa_key = ("Unknown","Unknown")
    
    seq_id_index = 0
    weighted_score_index = 9
    non_tp_mm_index = 4
    tp_mm_index = 5
    
    valid_score_types = ['overall_mismatches', 'tp_mismatches', 
     'weighted_score']
     
    if score_type not in valid_score_types:
        raise ValueError('score_type %s not a valid score type.' % score_type)
        
    if score_type == 'weighted_score':
        target_score_index = [weighted_score_index]
    elif score_type == 'tp_mismatches':
        target_score_index = [tp_mm_index]
    elif score_type == 'overall_mismatches':
        target_score_index = [non_tp_mm_index, tp_mm_index]
    
    for taxa_level in range(taxa_depth):
        
        # Add dict to contain categories for current taxa depth
        coverage_data.append({})
        
        
        for hits_index in range(len(pair_hits_data[f_hits_index])):
            
            hit_line_f = pair_hits_data[f_hits_index][hits_index].split(',')
            hit_line_r = pair_hits_data[r_hits_index][hits_index].split(',')

            # f and r data have already been tested for equality when data 
            # were loaded, but test anyway
            seq_id = hit_line_f[seq_id_index]
            seq_id_r = hit_line_r[seq_id_index]
            
            if seq_id != seq_id_r:
                raise ValueError('Forward and reverse sequence IDs do not '+\
                 'match.  See forward ID %s, reverse ID %s.' %\
                 (seq_id, seq_id_r))
            
            score_f = 0
            score_r = 0
            
            for score_index in target_score_index:
                score_f += float(hit_line_f[score_index])
                score_r += float(hit_line_r[score_index])
            
            if score_f <= score_threshold and score_r <= score_threshold:
                passes_score = 1
            else:
                passes_score = 0
            
            try:
                taxa = taxa_map[seq_id].split(';')
            except KeyError:
                # Add to 'Unknown' bin, only count when at domain level
                if taxa_level == domain_index:
                    total_seqs_no_taxa += 1
                    seqs_no_tax += seq_id + "\n"
                    try:
                        coverage_data[taxa_level][no_taxa_key][totals_index]+= 1
                        coverage_data[taxa_level][no_taxa_key][passing_index]+=\
                         passes_score
                    except KeyError:
                        coverage_data[taxa_level][no_taxa_key] = [1, 0]
                        coverage_data[taxa_level][no_taxa_key][passing_index]+=\
                         passes_score
                    continue
                else:
                    # Skip further processing if beyond domain level
                    continue
            
            # Determine current taxa-backstep if uncultured to find the
            # best name
            domain_taxa_key = get_curr_taxa(taxa, taxa_level, domain_index)
                
            if domain_taxa_key == None:
                continue
                
            
            # Add key if doesn't exist, otherwise increment counts
            try:
                coverage_data[taxa_level][domain_taxa_key][totals_index] += 1
                coverage_data[taxa_level][domain_taxa_key][passing_index] +=\
                 passes_score
            except KeyError:
                coverage_data[taxa_level][domain_taxa_key] = [1, 0]
                coverage_data[taxa_level][domain_taxa_key][passing_index] +=\
                 passes_score
                
                
            
        
    log_data.append("# Total sequences with no corresponding taxonomy "+\
     "in taxonomy mapping file: %d\n" % total_seqs_no_taxa)
    log_data.append("# Sequence IDs lacking taxonomy mapping:\n")
    log_data.append(seqs_no_tax)
    
    
    return coverage_data, log_data                       

def get_coverage_data_single_primer(hits_data,
                                    taxa_map,
                                    taxa_depth,
                                    score_type,
                                    score_threshold,
                                    hits,
                                    taxa_fp):
    """ Sorts coverage data into categories, total seqs, seqs passing score
    
    hits_data: lines of hit data from hits file
    taxa_map: dict of sequence ID: taxonomy (comma separated, highest level of
     taxa to lowest level)
    taxa_depth: level of taxa to attempt to categories.  This function will
     stop trying to add new levels to the returned list when the current list
     completes and is empty.
    score_type: Type of score from hits file to utilize.  Valid choices are
     weighted_score, overall_mismatches, and tp_mismatches.
    score_threshold: Value at which or below a particular sequence will be 
     considered to be successfully amplified during PCR.
    hits: filepath of current hits file tested.  Used for log file.
    taxa_fp: filepath of the taxonomy mapping file.  Used for log file.

    """
    
    coverage_data = []
    # log data will contain score type, threshold, total seqs analyzed, total
    # seqs lacking corresponding IDs in the mapping file, and the IDs of these
    # seqs not found in the mapping file.
    log_data = []
    
    # Add data to beginning of log file
    log_data.append(log_file_header)
    log_data.append("# Hits file analyzed: %s\n" % basename(hits))
    log_data.append("# Taxonomy mapping file used: %s\n" % basename(taxa_fp))
    log_data.append("# Scoring method used: %s\n" % score_type)
    log_data.append("# Score threshold: %2.2f\n" % score_threshold)
    
    total_seqs = len(hits_data)
    
    log_data.append("# Total sequences in hits file: %d\n" % total_seqs)
    
    total_seqs_no_taxa = 0
    seqs_no_tax = ""
    
    domain_index = 0
    
    totals_index = 0
    passing_index = 1
    
    no_taxa_key = ("Unknown","Unknown")
    
    seq_id_index = 0
    weighted_score_index = 9
    non_tp_mm_index = 4
    tp_mm_index = 5
    
    valid_score_types = ['overall_mismatches', 'tp_mismatches', 
     'weighted_score']
     
    if score_type not in valid_score_types:
        raise ValueError('score_type %s not a valid score type.' % score_type)
        
    if score_type == 'weighted_score':
        target_score_index = [weighted_score_index]
    elif score_type == 'tp_mismatches':
        target_score_index = [tp_mm_index]
    elif score_type == 'overall_mismatches':
        target_score_index = [non_tp_mm_index, tp_mm_index]
    
    for taxa_level in range(taxa_depth):
        
        # Add dict to contain categories for current taxa depth
        coverage_data.append({})
        
        for hit in hits_data:
            
            hit_line = hit.split(',')

            seq_id = hit_line[seq_id_index]
            
            score = 0
            for score_index in target_score_index:
                score += float(hit_line[score_index])
            
            if score <= score_threshold:
                passes_score = 1
            else:
                passes_score = 0
            
            try:
                taxa = taxa_map[seq_id].split(';')
            except KeyError:
                # Add to 'Unknown' bin, only count when at domain level
                # Used in case some sequences are not in mapping file.
                if taxa_level == domain_index:
                    total_seqs_no_taxa += 1
                    seqs_no_tax += seq_id + "\n"
                    try:
                        coverage_data[taxa_level][no_taxa_key][totals_index]+= 1
                        coverage_data[taxa_level][no_taxa_key][passing_index]+=\
                         passes_score
                    except KeyError:
                        coverage_data[taxa_level][no_taxa_key] = [1, 0]
                        coverage_data[taxa_level][no_taxa_key][passing_index]+=\
                         passes_score
                    continue
                else:
                    # Skip further checking if beyond domain level, only
                    # quanitify seqs lacking mapping once.
                    continue

            # Determine current taxa-backstep if uncultured to find the
            # best name
            domain_taxa_key = get_curr_taxa(taxa, taxa_level, domain_index)
                
            if domain_taxa_key == None:
                continue
                
            
            # Add key if doesn't exist, otherwise increment counts
            try:
                coverage_data[taxa_level][domain_taxa_key][totals_index] += 1
                coverage_data[taxa_level][domain_taxa_key][passing_index] +=\
                 passes_score
            except KeyError:
                coverage_data[taxa_level][domain_taxa_key] = [1, 0]
                coverage_data[taxa_level][domain_taxa_key][passing_index] +=\
                 passes_score
                
                
            
        
    log_data.append("# Total sequences with no corresponding taxonomy "+\
     "in taxonomy mapping file: %d\n" % total_seqs_no_taxa)
    log_data.append("# Sequence IDs lacking taxonomy mapping:\n")
    log_data.append(seqs_no_tax)
    
    
    return coverage_data, log_data
    
def get_curr_taxa(taxa,
                  taxa_level,
                  domain_index):
    """ Returns taxa string for current level, or None if mapping missing
    
    taxa: taxonomy for the current sequence ID, split into a list according to
     taxa depth (0 is domain level, 1 phylum, and so on).
    taxa_level: int with current taxonomy level being examined.
    domain_index: Should be zero, but left as variable in case differing
     taxonomy mapping format used.
    """
    
    try:
        taxa_current_level = taxa[taxa_level].capitalize()
        # Try to distinguish 'Uncultured' microbes according to higher
        # level taxonomy if available.  If uncultured part (but not
        # the whole name) of the higher level taxa, give this level
        # the name that already is partially named uncultured.  Or go
        # to a higher level to find name.  This is done to prevent 
        # redundant "uncultured" appearing in the final output, but 
        # also to separate 'uncultured' items in the current taxa level
        # according to known higher-level taxonomy.
        taxa_prefix = ""
        if taxa_current_level == "Uncultured":
            for back_taxa in range(taxa_level, domain_index, -1):
                curr_back_taxa = taxa[back_taxa].capitalize()
                if curr_back_taxa != "Uncultured":
                    if "uncultured" in curr_back_taxa.lower():
                        taxa_current_level = curr_back_taxa
                        break
                    else:
                        taxa_prefix = curr_back_taxa + "_"
                        break
        # Capitalize so taxonomic names are not case sensitive, still
        # need to be equivalent letters; 'Archaeal' is not equivalent
        # to 'Archaea' and will be binned separately.
        domain_taxa_key = (taxa[domain_index].capitalize(),
         taxa_prefix + taxa_current_level)
    except IndexError:
        # Skip current hit if taxa_level is deeper than taxa depth
        # given in the taxa mapping file
        domain_taxa_key = None
        
    return domain_taxa_key
    
    
def get_output_dir(output_dir,
                   hits,
                   primer_pair = False):
    """ Returns output subdirectory path for single primer results
    
    output_dir: base output directory
    hits: filepath of hits file for current primer analyzed
    primer_pair: If True, will split hits filenames to generate output dir
     based on names of both hits files.
    """
    
    if not output_dir.endswith('/'):
        output_dir += '/'
        
    if primer_pair:
        f_hits_index = 0
        r_hits_index = 1
        curr_primer_output_dir =\
         output_dir + basename(hits[f_hits_index]).split('.')[0] + "_" +\
         basename(hits[r_hits_index]).split('.')[0] + "_primers_coverage" + '/'
    else:
        curr_primer_output_dir = output_dir + basename(hits).split('.')[0] +\
         "_primer_coverage" + '/'
    
    return curr_primer_output_dir
    
def get_log_filepath(curr_primer_output_dir,
                     hits,
                     primer_pair = False):
    """ Returns output filepath for log file
    
    curr_primer_output_dir: subdirectory to output log file to
    hits: filepath of primer/primers, used to generate log file name.
    primer_pair: if True, then hits is a 2 element list containing the forward
     and reverse hits filepaths.
    """
    
    if primer_pair:
        f_hits_index = 0
        r_hits_index = 1
        log_filepath = curr_primer_output_dir +\
         basename(hits[f_hits_index]).split('.')[0] + "_" +\
         basename(hits[r_hits_index]).split('.')[0] + ".log"
        
    else:
        log_filepath = curr_primer_output_dir + basename(hits).split('.')[0] +\
         ".log"
        
    return log_filepath
    
    
    
def write_log_file(log_filepath,
                   log_data):
    """ Writes log file for current primer/primer pair
    
    log_filepath = filepath of output log file
    log_data: list of lines of log data
    """
    
    log_f = open(log_filepath, "w")
    
    for line in log_data:
        log_f.write(line)
        
    log_f.close()
    
def get_coverage_filepath(curr_primer_output_dir,
                          hits,
                          primer_pair = False):
    """ Generates name for text file output of percentage coverage for each taxa
    
    curr_primer_output_dir: subdirectory to output log file to
    hits: filepath of primer/primers, used to generate coverage file name.
    primer_pair: if True, then two hits filepaths have to be passed as a string,
     separated by a colon to generate name based on combination of the two.
    """
    
    if primer_pair:
        
        f_hits_index = 0
        r_hits_index = 1
        primer_coverage_filepath = curr_primer_output_dir +\
         basename(hits[f_hits_index]).split('.')[0] + "_" +\
         basename(hits[r_hits_index]).split('.')[0] + "_coverage.txt"
         
    else:
        primer_coverage_filepath = curr_primer_output_dir +\
         basename(hits).split('.')[0] + "_coverage.txt"
        
    return primer_coverage_filepath
     
    
    
    
    
def write_primer_coverage(primer_coverage_filepath,
                          coverage_data):
    """ Writes primer coverage percentages, in descending order of taxa levels
    
    primer_coverage_filepath: output filepath for text file of coverage.
    coverage_data: list of dictionaries containing taxonomic coverage for
     each level of taxa depth specified/available.
    """
    
    primer_coverage = open(primer_coverage_filepath, "w")
    
    primer_coverage.write(coverage_header)

    for taxa_level in coverage_data:
        curr_taxa_level = coverage_data.index(taxa_level)
        # Put all data in a list for sorting according to domain name
        curr_level_data = []
        for taxa in taxa_level:
            curr_level_data.append(("%d" % curr_taxa_level, taxa[0], taxa[1], 
             "%d" % taxa_level[taxa][0],
             "%1.4f\n" % (taxa_level[taxa][1]/taxa_level[taxa][0])))
             
        curr_level_data = sorted(curr_level_data,
         key = lambda domain: domain[1])
        for line in curr_level_data:
            primer_coverage.write(",".join(line))

        
    
def write_primer_coverage_graph(primer_graph_filepath,
                                graph_title,
                                coverage_data):
    """ Generates graphs for each level of taxonomic output
    
    Additionally, these graphs are separated according to domain (or the first
    level of taxa assigned to each sample).
    
    primer_graph_filepath: base name of graph filepath to be used.  Will have
     the domain and taxa level appended to name in final output.
    graph_title: Title for each graph created, has name of primer hits file
     used to generate data.
    coverage_data: List of dictionaries containing descending levels of 
     taxonomic coverage, along with total sequences and number of sequences
     that are at or below the score threshold for the given primer.
    """
    
    # The top level of taxa will be graphed together.  All subsequent levels
    # will be graphed separately according to the top level taxonomy (usually
    # the bacteria, archaeal, and eukaryotic domains).
    unknown_taxa = ("Unknown", "Unknown")
    
    top_level_taxa = {}
    
    for taxa in coverage_data[0]:
            if not taxa == unknown_taxa:
                top_level_taxa[taxa[0]] = 0
                
    top_level_index = 0
    colors = 'rbgycmk'
    
    figure_subtext = "Numeric values above bins represent\ntotal sequence "+\
     "counts for each set"
     
    y_label = "Percent Coverage"
    ytick_labels = range(0, 100, 20)
    
    for taxa_level in coverage_data:
        
        curr_taxa_level = coverage_data.index(taxa_level)

        # Combine the top level graph, all other graphs will be separated 
        # according to top level taxa
        
        if curr_taxa_level == top_level_index:
            
            curr_graph_filepath = get_graph_fp(primer_graph_filepath,
             curr_taxa_level)
             
            curr_graph_title =\
             graph_title + "Taxonomy Level %d\n" % curr_taxa_level +\
             figure_subtext
            
            fig = figure(figsize=(10, 30))
            
            figtext(0.5, 0.93, curr_graph_title, ha='center', color='black',
             size='medium')

             
            data_len = len(taxa_level) + 1
            ind = arange(2, data_len + 3)
            width = 0.35
            
            ind_count = 0
            ax = fig.add_subplot(30, 1, 1)
            
            graph_correction = 1
            
            x_labels = [""]
            
            for taxa_data in taxa_level:
                
                perc_coverage =\
                 100* taxa_level[taxa_data][1] / taxa_level[taxa_data][0]
                seq_count = taxa_level[taxa_data][0]
                y = ax.bar(ind[ind_count + graph_correction], perc_coverage, 
                 width, color = colors[ind_count%7])
                ax.text(y[0].get_x() + y[0].get_width() / 2,\
                 1.05 * y[0].get_height(), '%d' % seq_count, ha='center', \
                  va='bottom', fontsize=7)
                ind_count += 1
                x_labels.append(taxa_data[1])
                
            ax.set_ylabel(y_label, fontsize = 7)
            
            ax.set_xticks(ind + width / 2)

            
            ax.set_xticklabels(x_labels, rotation=45,\
             fontsize = 7, horizontalalignment='right' )
             
            ax.set_yticklabels(ytick_labels, fontsize = 7)
            
            ymin = 0
            ymax = 100
            ylim(ymin, ymax)
            
            
            savefig(curr_graph_filepath)
            clf()
            close()
            
        else:
            # Split according to domain.  Unknown only displayed at top domain
            
            domains = {}
            for domain in taxa_level:
                # Separate into dictionaries based on domain name, so should
                # have an archaea, bacteria, eukaryotic dictionary (if
                # there are sequences present under those domains.
                try:
                    domains[domain[0]][domain[1]] = taxa_level[domain]  
                except KeyError:
                    domains[domain[0]] = {}
                    domains[domain[0]][domain[1]] = taxa_level[domain]
                    
            for domain in domains:
                
                total_categories = len(domains[domain])
                
                corrected_figure_size = ceil(total_categories / 8) + 11
                if corrected_figure_size < 30:
                    corrected_figure_size = 30
                    
                
                subplot_size = ceil(total_categories / 6) + 11
               
                # Correct to avoid distended graphs for low number categories
                if subplot_size < 15:
                    subplot_size = 15
                
                
                curr_graph_filepath = get_graph_fp(primer_graph_filepath,
                 curr_taxa_level, domain = domain)
                 
                curr_graph_title =\
                 graph_title + "Sequences In Category %s\n" % domain +\
                 "Taxonomy Level %d\n" % curr_taxa_level + figure_subtext
                
                fig = figure(figsize=(10, corrected_figure_size))
                
                figtext(0.5, 0.95, curr_graph_title, ha='center', color='black',
                 size='medium')

                max_categories = 28
                width = 0.35
                
                graph_correction = 1
            
                x_labels = [""]
                
                ind_count = 0
                subplot_index = 1
                
                if corrected_figure_size > 59:
                    subplot_increment = 4
                else:
                    subplot_increment = 3
                
                category_count = 0
                
                ind = arange(max_categories + graph_correction)
                
                
                total_graphed = 0
                
                for taxa_data in domains[domain]:
                    
                    if category_count == 0:
                        ax = fig.add_subplot(subplot_size, 1,\
                         subplot_index)
                        subplot_index += subplot_increment
                    
                    category_count += 1
                    
                    
                    perc_coverage =\
                     100* domains[domain][taxa_data][1] /\
                     domains[domain][taxa_data][0]
                    
                    seq_count = domains[domain][taxa_data][0]
                    
                    y = ax.bar(ind[ind_count + graph_correction], \
                     perc_coverage, width, color = colors[ind_count%7])
                    
                    ax.text(y[0].get_x() + y[0].get_width() / 2,\
                     1.05 * y[0].get_height(), '%d' % seq_count, ha='center', \
                      va='bottom', fontsize=7)
                    
                    #ind_count = (ind_count + 1) % max_categories
                    ind_count += 1
                    
                    x_labels.append(taxa_data)
                    
                    total_graphed += 1
                    
                    if category_count == max_categories or\
                     total_graphed == total_categories:
                        
                        category_count = 0
                        
                        ax.set_ylabel(y_label, fontsize = 7)
                        ax.set_xticks(ind + width / 2)
                        ax.set_xticklabels(x_labels, rotation=45,\
                         fontsize = 7, horizontalalignment='right' )
                 
                        ax.set_yticklabels(ytick_labels, fontsize = 7)
                
                        ymin = 0
                        ymax = 100
                        ylim(ymin, ymax)
                        x_labels = [""]
                        ind_count = 0
                
                subplots_adjust(hspace=1.0)
                
                
                savefig(curr_graph_filepath, dpi = 300)
                clf()
                close()

            
    
    
def get_graph_fp(primer_graph_filepath,
                 curr_taxa_level,
                 domain = None):
    """ Gets graph name for plot at given taxonomy level
    
    primer_graph_filepath: root filepath to append taxonomy name to.
    curr_taxa_level: Depth of taxa, starting at 0 for top level taxa
    """
    
    if curr_taxa_level == 0:
        curr_graph_filepath = primer_graph_filepath +\
         "_taxonomy_level_%d.pdf" % curr_taxa_level
    else:
        curr_graph_filepath = primer_graph_filepath +\
         "_%s_taxonomy_level_%d.pdf" % (domain, curr_taxa_level)
         
    return curr_graph_filepath
    
    
    
    
    
def get_coverage_graph_root_filepath(curr_primer_output_dir,
                                     hits,
                                     primer_pair = False):
    """ Returns output graph filename for primer/primer pair tested
    
    curr_primer_output_dir: subdirectory to output log file to
    hits: filepath of primer/primers, used to generate graph file name.
    primer_pair: if True, then two hits filepaths have to be passed as a string,
     separated by a colon to generate name based on combination of the two.
    """
    
    if primer_pair:
        
        f_hits_index = 0
        r_hits_index = 1
        primer_graph_filepath = curr_primer_output_dir +\
         basename(hits[f_hits_index]).split('.')[0] + "_" +\
         basename(hits[r_hits_index]).split('.')[0] + "_coverage"
        graph_title = "Predicted Taxonomic Coverage\n%s\n" %\
         (basename(hits[f_hits_index]).split('.')[0] + "_" +\
         basename(hits[r_hits_index]).split('.')[0])
    else:
        primer_graph_filepath = curr_primer_output_dir +\
         basename(hits).split('.')[0] + "_coverage"
        graph_title = "Predicted Taxonomic Coverage\n%s\n" %\
         basename(hits).split('.')[0]
        
    return primer_graph_filepath, graph_title
    

def get_primer_pairs(hits_f):
    """ Finds all forward and reverse primer hits paths
    
    This function returns a list of tuples containing every possible forward
    and reverse primer hits file.
    
    hits_f: list of hits filepaths
    """
    
    forward_hits = []
    reverse_hits = []
    
    primer_pairs = []
    
    for hits in hits_f:
        if get_primer_direction(hits) == 'f':
            forward_hits.append(hits)
        else:
            reverse_hits.append(hits)
            
    for f_hits in forward_hits:
        for r_hits in reverse_hits:
            primer_pairs.append((f_hits, r_hits))
            
    return primer_pairs


def graph_taxa_coverage(hits_fps,
                        taxa_fp,
                        taxa_depth,
                        primer_pairs,
                        all_files,
                        output_dir,
                        score_type,
                        score_threshold):
    """ Main program function for generating taxa coverage graphs, summaries
    
    hits_fps: hits file paths.  Multiple hits filepaths should be separated 
     with a colon.  If all_files is True, this should point to a directory of
     hits files.
    taxa_fp: taxonomy mapping filepath.
    taxa_depth: Level of taxa to generate graphs and text summaries.  Depth 
     depends on how specific the taxa are defined in the taxonomy mapping file.
    primer_pairs: If True, will find all forward and reverse primer pairs from
     the input hits files, and generate taxonomy coverage for these primer 
     pairs based on the worst scoring of the two.
    all_files: Set to True to look for all hits files in the directory specified
     with hits_fps.
    output_dir: output directory
    score_type: Changes values taken from the hits file for scoring the primer.
     valid choices are weighted_score, overall_mismatches, and tp_mismatches.
    score_threshold: Value at which or below the primer amplification is 
     considered a success.
    """
    
    # Get dictionary of sequence ID: taxonomy
    taxa_f = open(taxa_fp, "U")
    taxa_map = parse_taxa_mapping_file(taxa_f)
    
    hits_f = get_hits_files(hits_fps, all_files)
    
    for hits in hits_f:
        # Get list of lines from current hits file
        hits_data = get_primer_hits_data(hits)
        
        # Get list of lists.  Top level of lists indicates the depth of taxa to  
        # be plumbed.  The first element of the list is the domain level, the
        # second the phylum, and so on.  Each of these elements is a list of
        # dictionaries, containing tuples of (domain,categories) for said level,
        # the value being a list, with element zero being the total number of
        # sequences for that category, and element one being the total number
        # of sequences that are equal to or under the score threshold for 
        # whichever score_type is being tested.  
        coverage_data, log_data = get_coverage_data_single_primer(hits_data, 
         taxa_map, taxa_depth, score_type, score_threshold, hits, taxa_fp)
         
        curr_primer_output_dir = get_output_dir(output_dir, hits)
         
        create_dir(curr_primer_output_dir)
        
        log_filepath = get_log_filepath(curr_primer_output_dir, hits,
         primer_pair = False)
        
        write_log_file(log_filepath, log_data)
        
        primer_coverage_filepath = get_coverage_filepath(curr_primer_output_dir,
         hits, primer_pair = False)
        
        write_primer_coverage(primer_coverage_filepath, coverage_data)
        
        primer_graph_filepath, graph_title =\
         get_coverage_graph_root_filepath(curr_primer_output_dir, 
         hits, primer_pair = False)
         
        write_primer_coverage_graph(primer_graph_filepath, graph_title,
         coverage_data)
                 
    if primer_pairs:
        
        primer_f_index = 0
        primer_r_index = 1
        
        primer_pairs = get_primer_pairs(hits_f)
        
        if not primer_pairs:
            raise ValueError('No valid forward and reverse primer pairs '+\
             'found, please check input filepaths/directory if using the -p '+\
             'parameter.')
             
        for primer_pair in primer_pairs:
            
            hits_f_fp = primer_pair[primer_f_index]
            hits_r_fp = primer_pair[primer_r_index]
            
            # Returns 2D list of hits data, checks for congruence between hits
            # data for forward and reverse primer.
            pair_hits_data = get_primer_hits_data_pair([hits_f_fp, hits_r_fp])

            
            coverage_data, log_data =\
             get_coverage_data_primer_pair(pair_hits_data,
             taxa_map, taxa_depth, score_type, score_threshold,
             hits_f_fp, hits_r_fp, taxa_fp)
             
            curr_primer_output_dir = get_output_dir(output_dir,
             [hits_f_fp, hits_r_fp], primer_pair = True)
             
            create_dir(curr_primer_output_dir)
            
            log_filepath = get_log_filepath(curr_primer_output_dir, 
             [hits_f_fp, hits_r_fp], primer_pair = True)
            
            write_log_file(log_filepath, log_data)
            
            primer_coverage_filepath =\
             get_coverage_filepath(curr_primer_output_dir, 
             [hits_f_fp, hits_r_fp], primer_pair = True)
            
            write_primer_coverage(primer_coverage_filepath, coverage_data)
            
            primer_graph_filepath, graph_title =\
             get_coverage_graph_root_filepath(curr_primer_output_dir, 
             [hits_f_fp, hits_r_fp], primer_pair = True)
             
            write_primer_coverage_graph(primer_graph_filepath, graph_title,
             coverage_data)
                 

             
        
        
    
    
    

log_file_header = """# Taxonomic coverage report
# Sequence IDs without corresponding IDs in the input taxonomy mapping file are listed at the end of this file.
"""

coverage_header = """# This file is written in descending levels of taxonomic depth
# Each line beings with the taxonomy level, followed by the first level of 
# taxonomy for a given sequence, generally the domain, followed the taxonomy 
# for the current level.
# The total sequences for a given taxonomy are listed first, followed by the
# percentage that are lower to or equal to the threshold score for passing.
# taxonomy level, first level taxonomy, taxonomy classification for current level, total seqs for given classification, percent seqs passing score
"""
