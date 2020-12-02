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

import re
lower = str.lower
from os.path import basename
from tempfile import NamedTemporaryFile
strip = str.strip
from itertools import count

from cogent3.parse.fasta import MinimalFastaParser
from cogent3.app.rdp_classifier import assign_taxonomy,\
 train_rdp_classifier_and_assign_taxonomy

from primerprospector.parse import parse_taxa_mapping_file
from primerprospector.format import write_assigned_taxa, write_accuracy_report

""" Assigns taxa for given reads/amplicons, generates report on accuracy
 to a given depth of taxa.
 
 Only RDP classifier implemented currently"""


def generate_training_files(training_data_fp,
                            taxa_mapping_fp):
    """Returns a tuple of file objects suitable for passing to the
    RdpTrainer application controller.
    
    Function modified from assign_taxonomy.py QIIME code, written by
    Kyle Bittinger.
    
    training_data_fp: filepath to fasta sequences used for retraining RPD
     classifier.
    taxa_mapping_fp: filepath to taxonomy mapping file.  Must have matching
     IDs to training_data_fp sequences.  In format of sequence ID<tab> six 
     levels of taxonomy, starting with "Root", in descending taxonomic depth,
     separated by semicolons.
    """
    id_to_taxonomy_file = open(taxa_mapping_fp, 'U')
    reference_seqs_file = open(training_data_fp, 'U')

    # Generate taxonomic tree and write to file 
    tree = build_tree(id_to_taxonomy_file)
    id_to_taxonomy_file.seek(0)

    rdp_taxonomy_file = NamedTemporaryFile(
        prefix='RdpTaxonAssigner_taxonomy_', suffix='.txt')
    rdp_taxonomy_file.write(tree.rdp_taxonomy())
    rdp_taxonomy_file.seek(0)

    # Generate a set of training seqs and write to file
    training_seqs = generate_training_seqs(reference_seqs_file,
     id_to_taxonomy_file)

    rdp_training_seqs_file = NamedTemporaryFile(
        prefix='RdpTaxonAssigner_training_seqs_', suffix='.fasta')
    for rdp_id, seq in training_seqs:
        rdp_training_seqs_file.write('>%s\n%s\n' % (rdp_id, seq))
    rdp_training_seqs_file.seek(0)
        

    return rdp_taxonomy_file, rdp_training_seqs_file
    
def generate_training_seqs(reference_seqs_file, id_to_taxonomy_file):
    """Returns an iterator of valid training sequences in
    RDP-compatible format
    
    Function modified from assign_taxonomy.py QIIME code, written by
    Kyle Bittinger.

    Each training sequence is represented by a tuple (rdp_id,
    seq).  The rdp_id consists of two items: the original sequence
    ID with whitespace replaced by underscores, and the lineage
    with taxa separated by semicolons.
    """
    # Rdp requires unique sequence IDs without whitespace.  Can't
    # trust user IDs to not have whitespace, so we replace all
    # whitespace with an underscore.  Classification may fail if
    # the replacement method generates a name collision.

    id_to_taxonomy_map = parse_id_to_taxonomy_file(id_to_taxonomy_file)

    for id, seq in MinimalFastaParser(reference_seqs_file):
        taxonomy = id_to_taxonomy_map[id]
        lineage = parse_lineage(taxonomy)
        rdp_id = '%s %s' % (re.sub('\s', '_', id), ';'.join(lineage))
        yield rdp_id, seq
    
def parse_id_to_taxonomy_file(f):
    """ parse the id_to_taxonomy file into a dict mapping id -> taxonomy
    
    Function modified from assign_taxonomy.py QIIME code, written by
    Kyle Bittinger.
    """
    result = {}
    for line in f:
        line = line.strip()
        if line:
            identifier, taxonomy = map(strip, line.split('\t'))
            result[identifier] = taxonomy
    return result 

def build_tree(id_to_taxonomy_file):
    """Returns an RdpTree object representing the taxonomic tree
    derived from the id_to_taxonomy open file object.
    
    Function modified from assign_taxonomy.py QIIME code, written by
    Kyle Bittinger.
    
    id_to_taxonomy_file: Open file object of taxonomy mapping file.
    """

    
    tree = RdpTree()
    for line in id_to_taxonomy_file:

        id, taxonomy = map(strip, line.split('\t'))
        lineage = parse_lineage(taxonomy)
        tree.insert_lineage(lineage)
    return tree
    
    
def parse_lineage(lineage_str):
    """Returns a list of taxa from the semi-colon-separated
    lineage string of an id_to_taxonomy file.
    
    Function modified from assign_taxonomy.py QIIME code, written by
    Kyle Bittinger.
    """
    lineage = lineage_str.strip().split(';')
    # The RDP Classifier can only deal with a lineage that is 6
    # levels deep.  We detect this problem now to avoid an
    # ApplicationError later on.
    if len(lineage) != 6:
        raise ValueError(
            'Each reference assignment must contain 6 items, specifying '
            'domain, phylum, class, order, family, and genus.  '
            'Detected %s items in "%s": %s.' % \
            (len(lineage), lineage_str, lineage))
    return lineage
        
# The RdpTree class is defined as a nested class to prevent the
# implementation details of creating RDP-compatible training files
# from being exposed at the module level.
class RdpTree(object):
    """Simple, specialized tree class used to generate a taxonomy
    file for the Rdp Classifier.
    
    Function modified from assign_taxonomy.py QIIME code, written by
    Kyle Bittinger.
    """
    taxonomic_ranks = ['domain', 'phylum', 'class', 
                       'order', 'family', 'genus']

    def __init__(self, name='Root', parent=None):
        self.name = name
        self.parent = parent
        if parent is None:
            self.depth = 0
        else:
            self.depth = parent.depth + 1
        self.children = dict()  # name => subtree

    def insert_lineage(self, lineage):
        """Inserts a lineage into the taxonomic tree.
        
        Lineage must support the iterator interface, or provide an
        __iter__() method that returns an iterator.
        """
        lineage = lineage.__iter__()
        try:
            taxon = next(lineage)
            if taxon not in self.children:
                self.children[taxon] = RdpTree(name=taxon, parent=self)
            self.children[taxon].insert_lineage(lineage)
        except StopIteration:
            pass
        return self

    def rdp_taxonomy(self, counter=None):
        """Returns a string, in Rdp-compatible format.
        """
        if counter is None:
            # initialize a new counter to assign IDs
            counter = count(0)

        # Assign ID to current node; used by child nodes
        self.id = next(counter)

        if self.parent is None:
            # do not print line for Root node
            retval = ''
        else:
            # Rdp taxonomy file does not count the root node
            # when considering the depth of the taxon.  We
            # therefore specify a taxon depth which is one
            # less than the tree depth.
            taxon_depth = self.depth - 1

            # In this simplest-possible implementation, we
            # also use the taxon depth to retrieve the string
            # value of the taxonomic rank.  Taxa beyond the
            # depth of our master list are given a rank of ''.
            try:
                taxon_rank = self.taxonomic_ranks[taxon_depth]
            except IndexError:
                taxon_rank = ''

            fields = [self.id, self.name, self.parent.id, taxon_depth,
                      taxon_rank]
            retval = '*'.join(map(str, fields)) + "\n"

        # Recursively append lines from sorted list of subtrees
        child_names = list(self.children.keys())
        child_names.sort()
        subtrees = [self.children[name] for name in child_names]
        for subtree in subtrees:
            retval += subtree.rdp_taxonomy(counter)
        return retval
        
        

def get_taxa_assignment_output_fp(output_dir,
                                  fasta_fp):
    """ Gets output filepath for the taxonomy assignment file
    
    output_dir: output directory
    fasta_fp: input fasta filepath
    """
    
    if not output_dir.endswith("/"):
        output_dir += '/'
        
    taxa_assignment_outf = output_dir + basename(fasta_fp).split('.')[0] +\
     "_assignments.txt"
    
    return taxa_assignment_outf
    
def get_report_output_fp(output_dir,
                         fasta_fp):
    """ Gets output filepath for accuracy report
    
    output_dir: output directory
    fasta_fp: input fasta filepath
    """
    if not output_dir.endswith("/"):
        output_dir += '/'
        
    report_outf = output_dir + basename(fasta_fp).split('.')[0] +\
     "_accuracy_report.txt"
    
    return report_outf
    
    

        

def get_accuracy_report(taxa_map, 
                        assigned_taxa,
                        taxa_depth):
    """ Gets accurate assignemtn percent for depth of taxa specified
    
    taxa_map:  dict of sequence ID : taxonomy (known taxa from taxonomy 
     mapping file)
    assigned_taxa: dict of sequence ID : (assigned taxonomy, confidence)
    taxa_depth: level of assigned taxa to test for accuracy
    """
    
    accuracy_values = []
    
    taxa_index = 0
    
    seqs_lacking_taxa_mapping = 0
    seqs_assigned_root = 0
    
    # taxa assigners add "Root" to the beginning of the taxonomy, go to 
    # next index to correct for this
    root_correction = 1
    
    # Correction for "Root" assignments
    assignment_correction = 0
    
    for curr_level in range(taxa_depth):
        total_seqs = 0
        correct_assignment = 0
        
        
        
        for assignment in assigned_taxa:
            
            try:
                # Taxonomy split on semicolons
                actual_taxonomy = \
                 taxa_map[assignment].split(';')[curr_level +\
                 assignment_correction]
                if actual_taxonomy == "Root":
                    assignment_correction = 1
                    actual_taxonomy = \
                    taxa_map[assignment].split(';')[curr_level +\
                    assignment_correction]
                actual_taxonomy = lower(actual_taxonomy)
            except IndexError:
                # Skip
                continue
            except KeyError:
                # Count up total at domain level
                if curr_level == 0:
                    seqs_lacking_taxa_mapping += 1
                continue
            
            # Classifier splits taxa levels on semicolons
            # skip if assigned taxa not as deep as taxa_depth
            try:
                if assigned_taxa[assignment][taxa_index] == "Root" and\
                 curr_level == 0:
                    seqs_assigned_root += 1
                    continue
                curr_assignment = \
                 assigned_taxa[assignment][taxa_index].split(';')[curr_level +\
                 root_correction]
                # Strip off any quotes from assignment
                curr_assignment = curr_assignment.replace('"', '')
                curr_assignment = lower(curr_assignment)
            except IndexError:
                # Skip
                continue
                
            
                
            total_seqs += 1

            
            if curr_assignment == actual_taxonomy:
                correct_assignment += 1
                
        if total_seqs == 0:
            # end tests if unable to test any assignments at this level
            break
        accuracy_values.append((100 * correct_assignment / total_seqs, 
         total_seqs))
                
    
    return accuracy_values, seqs_lacking_taxa_mapping, seqs_assigned_root
    

def generate_taxa_report(taxa_mapping_fp,
                         fasta_fp,
                         taxa_depth,
                         output_dir,
                         assignment_method = 'rdp',
                         min_confidence = 0.80,
                         training_data_fp = None):
    """ Main program function for assignment and accuracy report of taxa 
    
    taxa_mapping_fp: filepath of taxonomy mapping file
    fasta_fp: filepath of fasta file to test
    taxa_depth: level of taxa in results and taxonomy mapping to compare in
     accuracy report.  Depends on both resulting depth and input depth.
    output_dir: output directory
    assignment_method: taxonomy assignment method
    min_confidence: minimum confidence for taxonomic assignment
    training_data_fp: training data filepath for RDP classifier
    """
    
    taxa_map_f = open(taxa_mapping_fp, "U")
    
    
    taxa_map = parse_taxa_mapping_file(taxa_map_f)

    fasta_seqs = open(fasta_fp, "U")
    
    taxa_assignment_outf = get_taxa_assignment_output_fp(output_dir,
     fasta_fp)
     
    report_outf = get_report_output_fp(output_dir, fasta_fp)
    
    if assignment_method == 'rdp':
        if training_data_fp:
            # Convert files to format compatable with RDP classifier
            taxonomy_file, training_seqs_file =\
             generate_training_files(training_data_fp, taxa_mapping_fp)
            
            assigned_taxa =\
             train_rdp_classifier_and_assign_taxonomy(training_seqs_file,
             taxonomy_file, fasta_seqs, min_confidence = min_confidence)
        else:
            assigned_taxa = assign_taxonomy(fasta_seqs, 
             min_confidence = min_confidence)
        
    taxa_assignment_outf = open(taxa_assignment_outf, "w")
     
    write_assigned_taxa(assigned_taxa, taxa_assignment_outf)
     
    # Get a list of (floats, total seqs at level) for the accuracy of 
    # assignment at each depth of taxa assignment tested, 
    accuracy_values, seqs_lacking_taxa_mapping, seqs_assigned_root =\
     get_accuracy_report(taxa_map, assigned_taxa, taxa_depth)
     
    report_outf = open(report_outf, "w")
     
    write_accuracy_report(accuracy_values, seqs_lacking_taxa_mapping, 
     seqs_assigned_root, report_outf, fasta_fp, assignment_method, 
     training_data_fp)
    

