.. _essential_files:

===================
File Format Details
===================


This section covers the details of file formats for input and output files.  Most files are described in the `tutorial <../tutorial/tutorial.html>`_, however a few files are described here more fully for advanced users who may wish to parse data directly from Primer Prospector output.  Sample files can be found in the ``pprospector/doc/tutorial_files/`` directory.


^^^^^^^^^^^^^^^^^^
Input File Formats
^^^^^^^^^^^^^^^^^^

+++++++++++
Fasta Files
+++++++++++

Fasta files must comply with standard format (label line preceded by the '>' character, followed by sequences).  Modules in the *de novo* pipeline require aligned fasta files, while the analysis pipeline utilizes unaligned (degapped) fasta files.  Aligned fasta files can be degapped, or otherwise cleaned up for use with the `clean_fasta.py <../scripts/clean_fasta.html>`_ module.

+++++++++++++++++++
Primers File Format
+++++++++++++++++++

Primers supplied in text files to Primer Prospector should be tab separated with sequences listed in 5' to 3' direction.  Comments are preceded by a pound (#) character.  See the example ``primers.txt`` file in the ``pprospector/doc/tutorial_files/analysis_files/`` directory.


++++++++++++++++++++++
Taxonomy Mapping Files
++++++++++++++++++++++

Taxonomy mapping files are used by several modules in the primer analysis pipeline.  A few example lines of a taxonomy mapping file are shown below:

::

	103547	Bacteria;Proteobacteria;Gammaproteobacteria;Legionellales;Coxiellaceae;Coxiella
	214073	Bacteria;Cyanobacteria;Cyanobacteria;Chloroplast;Streptophyta
	254376	Bacteria;Firmicutes;"Clostridia";Clostridiales;"Lachnospiraceae";Roseburia

The proper format is the sequence ID, followed by a tab, followed by semicolon (;) separated taxonomies, starting with the highest level (domain) and descending to more specific levels.  See the ``taxonomy_mapping.txt`` file in the ``pprospector/doc/tutorial_files/analysis_files/`` as an example.

A Greengenes mapping file is available in a `zip <http://greengenes.lbl.gov/Download/OTUs/gg_otus_6oct2010.zip>`_ file.

^^^^^^^^^^^^^^^^^^^
Output File Formats
^^^^^^^^^^^^^^^^^^^

+++++++++++++++++
*De Novo* Primers
+++++++++++++++++

The `generate_primers_denovo.py <../scripts/generate_primers_denovo.html>`_ finds short conserved sites and records these and the surrounding sequences.  An example output file is below (see ``pprospector/doc/tutorial_files/de_novo_files/final_denovo_primers.txt`` for the complete file):

::

	# Matches that are found in at least 95.00% target sequences and specific at or below the 1.00% threshold 
	# record type,sequence,unaligned index,aligned index,number of perfect matches,percent perfect matches,percent non-specific matches,forward primer standard index,reverse primer standard index
	C,GCGTG,371,1974,18,100.0%,0.0%,388,423
	M,1974,TGACGCAGCGACGCCGCGTG,GCGTGAGGGATGACGGCCTT,11544
	M,1974,TGATGCAGCAACGCCGCGTG,GCGTGGAGGATGACGCATTT,57484
	M,1974,TGATCCAGCAATGCCGCGTG,GCGTGTGTGAAGAAGGCCCT,103547
	<snip>

The header lines show the sensitivity and specificity thresholds.  In the following lines, the 'C' line represents the conserved short sequence that was found, along with other details about where this sequence was found originally.  In the above example, the unaligned (degapped) index where the sequence GCGTG was original found was 371, while the aligned (gapped) index was 1974.  18 of the input sequences (100% of them) had an exact match for this sequence, and the sequence was found in 0.0% of the unwanted (specificity test) file.  The standard forward and reverse indices was the position of this sequence in an *E. coli* sequence, corrected for the length of the prospective primer (20 base pairs).

Lines that start with "M" indicate the sequences with a perfect match to the short sequence GCGTG.  Because sequences can be repeated, the aligned index, 1974, is also recorded to guarantee uniqueness for each conserved site recorded.  The following two sequences are the upstream and downstream sequences of the conserved site (default setting records 15 base pairs up and downstream).  The final number is the sequence identifier, pulled from the fasta label.

++++++++++++++++++++++++
Sorted *De Novo* Primers
++++++++++++++++++++++++

Most of the output from `sort_denovo_primers.py <../scripts/sort_denovo_primers.html>`_  is covered in the `tutorial <../tutorial/tutorial.html>`_.  The ``primer_details.txt`` file is more fully covered here (an example file is present in the ``pprospector/doc/tutorial_files/de_novo_files/`` directory):

::

	# Primer Report
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
	GGAAT,330,1908,18,100.0%,0.0%,347,382,GGGAGGCAGCAGTVRGGAAT,GGGAGGCAGCAGTGRGGAAT,GGGAGGCAGCAGTGGGGAAT,YSHSCATTGNSSAADATTCC,CSCCCATTGYSCAATATTCC,CGCCCATTGTGCAATATTCC
	CGCGT,370,1972,18,100.0%,0.0%,387,422,BYKRNSSRRCVAYGCCGCGT,CTGAYSCAGCRAYGCCGCGT,CTGATGCAGCAACGCCGCGT,RRDVCBTHNTCVYHCACGCG,ARKMCTTCWTCVYTCACGCG,AGGCCTTCTTCACTCACGCG
	CGGGA,314,1877,18,100.0%,0.0%,331,366,CGGNCCMRACTCCTACGGGA,CGGCCCAGACTCCTACGGGA,CGGCCCAGACTCCTACGGGA,TTCCYBACTGCTGCCTCCCG,TTCCYCACTGCTGCCTCCCG,TTCCCCACTGCTGCCTCCCG
	<snip>

Much of the data are the same that are stored in the initial *de novo* primer output file.  However, it is worth explaining the IUPAC, filtered, and consensus sequences.  As mentioned for the *de novo* output, the up and downstream sequences from the conserved sites are recorded.  These sequences are combined in this sorted output.  The IUPAC sequence contains degeneracy for any nucleotide that occurs.  For the degeneracy to be present in the filtered line, the nucleotide must be present at least 20% of the time (default setting).  The consensus sequence shows the nucleotides that occur most frequently, with no degeneracy.


++++++++++++++++++++++++++
Primer Analysis Hits Files
++++++++++++++++++++++++++

The primer hits files generated by `analyze_primers.py <../scripts/analyze_primers.html>`_ are used by many of the modules in the primer analysis pipeline.  An example hits file is shown below (more examples can be found in the ``pprospector/doc/tutorial_files/analysis_files/`` directory):


::

	# Primer: 27f 5'-AGAGTTTGATCMTGGCTCAG-3'
	# Input fasta file: sample_bacteria_seqs_unaligned.fasta
	# Parameters
	# 3' length: 5
	# non 3' mismatch penalty: 0.40 per mismatch
	# 3' mismatch penalty: 1.00 per mismatch
	# last base mismatch penalty: 3.00
	# non 3' gap penalty: 1.00 per gap
	# 3' gap penalty: 3.00 per gap
	# Note - seq hit and primer hit are the best local pairwise alignment results for a given sequence and primer pair.  A gap in seq hit represents a deletion in the sequence, whereas a gap in the primer hit signifies an insertion in the target sequence.
	#
	# seq ID, seq hit, primer hit, hit start position, non 3' mismatches, 3' mismatches (except last base), last base mismatch, non 3' gaps, 3' gaps, overall weighted score, hits sequence end 
	11544,TTGATCCTGGCTCAG,TTGATCMTGGCTCAG,0,0,0,False,0,0,0.0,True
	57484,AGGGTGCAAGCGTTACTCGG,AGAGTTTGATCMTGGCTCAG,478,8,1,False,0,0,4.2,False
	103547,AGAGTTTGATTCTGGCTCAG,AGAGTTTGATCMTGGCTCAG,5,1,0,False,0,0,0.4,False
	120123,AGTGAAAGCCCGGGGCTCAA,AGAGTTTGATCMTGGCTCAG,551,8,0,True,0,0,6.2,False
	127471,AGAGTTTGATCCTGGCTCAG,AGAGTTTGATCMTGGCTCAG,0,0,0,False,0,0,0.0,True
	<snip>


Header lines contain details about the scoring parameters as well as the primer and sequences tested.  There are a few details in the data lines that may not be immediately apparent.  Both the sequence hit and primer hit sequence are recorded, as an insertion in the target sequence will create a gap in the primer hit, while a deletion in the target sequence will create a gap in the sequence hit.  The hit start position records where the primer binds the target sequence.  Scoring is fully described in the `tutorial <../tutorial/tutorial.html>`_.  Finally, the hits sequence end value is stored as True or False.  This value does not affect primer scoring.  Primers can abut sequence ends if, for example, sequences were generated with this particular primer and submitted without removal of said primer.





