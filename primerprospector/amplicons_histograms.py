#!/usr/bin/env python

from __future__ import division

__author__ = "William A. Walters"
__copyright__ = "Copyright 2010, The PrimerProspector project"
__credits__ = ["William A. Walters"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William A. Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"

""" amplicons_histogram's purpose is to provide a graphical output of
    amplicon sizes.  See the get_amplicons_and_reads.py module for details about
    the generation of the fasta files containing predicted amplicons.  In short,
    the fasta files contain the sequence ID followed by the predicted amplicon 
    between two given primers that are expected to amplify said sequence.
    """
    
from os.path import basename
upper = str.upper

from numpy import arange
from matplotlib import use
use('Agg')
from matplotlib.pyplot import ylim
from pylab import figure, savefig, figtext, yticks, xticks
from cogent.parse.fasta import MinimalFastaParser

from primerprospector.parse import parse_taxa_mapping_file, \
 get_amplicons_filepaths
 

        
def get_amplicons_lens(amplicons):
    """ Reads amplicons file, returns dict of seq ID: seq length 
    
    amplicons: open file object (or list) of amplicons in fasta format"""
    
    
    # Create dictionary with sequence label as key, len of amplicon as value
    raw_amplicons_data={}
    
    for label,seq in MinimalFastaParser(amplicons):
        raw_amplicons_data[label]=len(seq)
        
    return raw_amplicons_data
    
def get_histogram_name(amplicons_filepath, 
                       output_dir):
    """ Generates output filename for histograms
    
    amplicons_filepath: filepath for current amplicons
    output_dir: output directory
    """
    
    histogram_name = output_dir + "/" +\
     basename(amplicons_filepath).split('.')[0] + ".ps"
     
    return histogram_name
    
    
def sort_combined_amplicons(raw_amplicons_data):
    """ Sort amplicons into dictionary of value: counts for that value
    
    raw_amplicons_data: dictionary of seq ID: seq(len)
    """
    
    amplicons_data = {}
    
    for amplicon in raw_amplicons_data:
        try:
            amplicons_data[raw_amplicons_data[amplicon]] += 1
        except KeyError:
            amplicons_data[raw_amplicons_data[amplicon]] = 1
            
    return amplicons_data
    
    
    
def sort_amplicons(raw_amplicons_data,
                   taxa_mapping):
    """ Sorts amplicons into separate domain dictionaries
    
    raw_amplicons_data: dictionary of seq ID: seq(len)
    taxa_mapping: dict of seq ID: taxonomy
    """
    
    arch_amplicons={}
    bact_amplicons={}
    euk_amplicons={}
    other_amplicons={}
    
    # Added the "K__" sections in case using RDP compatable mapping file.

    for amplicon in raw_amplicons_data:
        try:
            domain = taxa_mapping[amplicon].split(';')[0]
        except KeyError:
            raise KeyError('Sequence ID %s not found in ' %\
             amplicon + 'taxonomy mapping file.')
        if upper(domain).count("ARCHAEA"):
            try:
                arch_amplicons[raw_amplicons_data[amplicon]] += 1
            except KeyError:
                arch_amplicons[raw_amplicons_data[amplicon]] = 1
        elif upper(domain).count("EUKARYA") or upper(domain).count("EUKARYOTA"):
            try:
                euk_amplicons[raw_amplicons_data[amplicon]] += 1
            except KeyError:
                euk_amplicons[raw_amplicons_data[amplicon]] = 1
        elif upper(domain).count("BACTERIA"):
            try:
                bact_amplicons[raw_amplicons_data[amplicon]] += 1
            except KeyError:
                bact_amplicons[raw_amplicons_data[amplicon]] = 1
        else:
            try:
                other_amplicons[raw_amplicons_data[amplicon]] += 1
            except KeyError:
                other_amplicons[raw_amplicons_data[amplicon]] = 1

    
    return [arch_amplicons, bact_amplicons, euk_amplicons, other_amplicons]

def generate_histogram_combined(amplicons_data,
                                amp_f,
                                output_dir):
    """ Generates histogram data with all seq lengths in single histogram """
    
    histogram_outf_name = get_histogram_name(amp_f, output_dir)
    
    title = basename(amp_f).split('.')[0]
    fig_num = 0
    
    y_label="Amplicon Count"
    fig = figure(fig_num, figsize = (8, 10))
    fig.text(0.5, 0.93, title, ha="center", size='large')

    width=0.45

    xtitle = title
    plot_number = 1
    
    ax = fig.add_subplot(4, 1, plot_number)
    ind = arange((max(amplicons_data.keys()) + 1))
    for n in amplicons_data:
        y = ax.bar(ind[n], amplicons_data[n], width, color='y')
    ax.set_ylabel(y_label)

    ax.set_xlabel("Amplicon Size in Base Pairs")
    
    max_data_points = max(amplicons_data.values()) + 1
    y_tick_step = int(0.2 * max_data_points) or 2
    yticks(range(0, max_data_points+2, y_tick_step), size = 7)
    
    bin_range = max(amplicons_data.keys()) + 1
    x_tick_step = int(0.2 * bin_range) or 1
    xticks(range(0, bin_range + 100, x_tick_step), size = 7)

    savefig(histogram_outf_name)
    
    fig.clf()
    
    return

def max_value_from_bins(data_sets):
    """ Returns largest single value in list of lists
    
    The purpose of this function is to find the largest value in the lists of
    lists.
    
    data_sets: list of lists containing amplicon size data
    """
    
    counts_all_max = []
    
    for data_set in data_sets:
        
        # skip empty bins
        if not(len(data_set)):
            continue
        
        max_value = max(data_set)
        
        counts_all_max.append(max_value)
        
    return max(counts_all_max)
    
    

def generate_histograms_domains(amplicons_data,
                                amp_f,
                                output_dir):
    """ Plots histogram based on counts in domain dictionaries 
    
    amplicons_data: 4 element list with dictionaries that contains bins of 
     amplicon size: counts.  Dictionaries are ordered in the following:
     [archaea, bacteria, eukarya, unclassified]
    amp_f: current amplicons filepath being analyzed, used to generate graph
     title.
    output_dir: output directory
    """
    
    arch_index = 0
    bact_index = 1
    euk_index = 2
    unclassified_index = 3
    
    # Get separate dictionaries from the input list of dicts, makes later code
    # easier to read
    arch_amplicons = amplicons_data[arch_index]
    bact_amplicons = amplicons_data[bact_index]
    euk_amplicons = amplicons_data[euk_index]
    other_amplicons = amplicons_data[unclassified_index]
    
    max_x_bins = max_value_from_bins([arch_amplicons, bact_amplicons, 
     euk_amplicons, other_amplicons])
    
    
    
    plot_titles=['Archaeal Amplicons' , 'Bacterial Amplicons',
     'Eukaryotic Amplicons' , 'Unclassified Amplicons']
     
    histogram_outf_name = get_histogram_name(amp_f, output_dir)
    
    title = basename(amp_f).split('.')[0]
    

    y_label = "Amplicon Count"
    x_label = "Amplicon Size in Base Pairs"
    
    fig = figure(0, figsize=(8, 10))
    fig.text(0.5, 0.95, title, ha='center', size='large')

    plot_number = 1
    ax = [1, 2, 3, 4, 5]
    width = 0.45
    
    # Plot archaeal amplicon histogram
    ax[plot_number] = fig.add_subplot(5, 1, plot_number)
    
    # Handle case of empty dictionary
    if len(arch_amplicons):
        max_data_points = max(arch_amplicons.values()) + 1
        ind = arange((max(arch_amplicons.keys()) + 1))
    else:
        max_data_points = 10
        ind = 0
    
    for n in arch_amplicons:
        ax[plot_number].bar(ind[n], arch_amplicons[n], width)

    y_tick_step = int(0.2*max_data_points) or 2
    yticks(range(0, max_data_points+2, y_tick_step), size=7)
    
    x_tick_step = int(0.2 * max_x_bins) or 1
    xticks(range(0, max_x_bins + 100, x_tick_step), size = 7)
    
    ax[plot_number].set_ylabel(y_label, size='small')
    xtitle = plot_titles[plot_number - 1]
    ax[plot_number].set_title(xtitle, size='small')
    ax[plot_number].set_xlabel(x_label, size='small')
    
    # Plot bacterial amplicon histogram
    plot_number += 1
    ax[plot_number] = fig.add_subplot(5, 1, plot_number)
    
    if len(bact_amplicons):
        max_data_points = max(bact_amplicons.values()) + 1
        ind = arange((max(bact_amplicons.keys()) + 1))
    else:
        max_data_points = 10
        ind = 0
    
    for n in bact_amplicons:
        ax[plot_number].bar(ind[n], bact_amplicons[n], width)
    
    y_tick_step = int(0.2*max_data_points) or 2
    yticks(range(0, max_data_points+2, y_tick_step), size=7)
    
    x_tick_step = int(0.2 * max_x_bins) or 1
    xticks(range(0, max_x_bins + 100, x_tick_step), size = 7)
    
    ax[plot_number].set_ylabel(y_label, size='small')
    xtitle = plot_titles[plot_number - 1]
    ax[plot_number].set_title(xtitle, size='small')
    ax[plot_number].set_xlabel(x_label, size='small')
    
    # Plot eukaryotic amplicon histogram
    plot_number += 1
    ax[plot_number] = fig.add_subplot(5, 1, plot_number)
    # handle case of empty dictionary
    if len(euk_amplicons):
        max_data_points = max(euk_amplicons.values()) + 1
        ind = arange((max(euk_amplicons.keys()) + 1))
    else:
        max_data_points = 10
        ind = 0
    
    for n in euk_amplicons:
        ax[plot_number].bar(ind[n], euk_amplicons[n], width)
    
    y_tick_step = int(0.2*max_data_points) or 2
    yticks(range(0, max_data_points+2, y_tick_step), size=7)
    
    x_tick_step = int(0.2 * max_x_bins) or 1
    xticks(range(0, max_x_bins + 100, x_tick_step), size = 7)
    
    ax[plot_number].set_ylabel(y_label, size='small')
    xtitle = plot_titles[plot_number - 1]
    ax[plot_number].set_title(xtitle, size='small')
    ax[plot_number].set_xlabel(x_label, size='small')

    # Plot unclassified amplicon histogram
    plot_number += 1
    ax[plot_number] = fig.add_subplot(5, 1, plot_number)

    # Handle case of empty dictionary
    if len(other_amplicons):
        max_data_points = max(other_amplicons.values()) + 1
        ind = arange((max(other_amplicons.keys()) + 1))
    else:
        max_data_points = 10
        ind = 0

    for n in other_amplicons:
        ax[plot_number].bar(ind[n], other_amplicons[n], width)
    
    y_tick_step = int(0.2 * max_data_points) or 2
    yticks(range(0, max_data_points+2, y_tick_step), size=7)
    
    x_tick_step = int(0.2 * max_x_bins) or 1
    xticks(range(0, max_x_bins + 100, x_tick_step), size = 7)
     
    
    ax[plot_number].set_ylabel(y_label, size='small')
    xtitle = plot_titles[plot_number - 1]
    ax[plot_number].set_title(xtitle, size='small')
    ax[plot_number].set_xlabel(x_label, size='small')


    fig.subplots_adjust(hspace = 1.5)
    savefig(histogram_outf_name)
    
    fig.clf()

    
    return

def generate_histograms(amplicons_filepath,
                        output_dir,
                        all_files, 
                        taxa_map_filepath):
    """ Main module for handling generation of histograms from amplicons 
    
    amplicons_filepath: filepath(s) to amplicons files generated with 
     get_amplicons_and_reads module
    output_dir: output directory
    all_files: If True, get all files ending with _amplicons.fasta from 
     directory specified with amplicons_filepath
    taxa_map_filepath: filepath for taxonomy mapping file, used to separate
     amplicons into domains for histograms"""
    
    amp_filepaths = get_amplicons_filepaths(amplicons_filepath, all_files)
    
    if taxa_map_filepath:
        try:
            taxa_map = open(taxa_map_filepath, "U")
        except IOError:
            raise IOError('Unable to open taxa mapping file, please check '+\
             'filepath.')
        
        taxa_mapping = parse_taxa_mapping_file(taxa_map)
        
    
    for amp_f in amp_filepaths:
        amplicons = open(amp_f, "U")
        raw_amplicons_data = get_amplicons_lens(amplicons)
        
        if taxa_map_filepath:
            # Will get a list of 4 dicts containing amplicon size bins for
            # each domain and unclassified
            amplicons_data = sort_amplicons(raw_amplicons_data, taxa_mapping)
            generate_histograms_domains(amplicons_data, amp_f, output_dir)
        else:
            amplicons_data = sort_combined_amplicons(raw_amplicons_data)
            generate_histogram_combined(amplicons_data, amp_f, output_dir)
    
        
    
    
