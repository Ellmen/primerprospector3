.. _taxa_assignment_report:

.. index:: taxa_assignment_report.py

*taxa_assignment_report.py* --  Assign taxonomies for input fasta sequences, generate report on accuracy of assignment by comparing to known taxonomies 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**



This module will attempt to assign taxonomies to the input fasta files, then 
compare these results to the known taxonomies supplied for each sequence ID 
in the input taxa mapping file.  Currently only the RDP classifier is 
implemented for assigning taxonomies.

The primary purpose of this module is to test amplicons and/or reads for a
prospective primer pair (see `get_amplicons_and_reads.py <./get_amplicons_and_reads.html>`_) for their 
phylogenetic usefulness.  Very short reads, or reads in a region of high
conservation will be less accurate in terms of taxonomic assignment.

This module will save the taxonomic assignments, as well as a report detailing
the percentage that were accurately assigned for each level of taxonomic depth.
An assignment of Archaea;Euryarchaeota;Thermoplasmatales for a sequence that
was actually Archaea;Euryarchaeota;Halobacteriales is accurate to the second
level of taxonomy, but not the third, and will be recorded as such.  The 
default depth of taxa to test is 3.



**Usage:** :file:`taxa_assignment_report.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-t, `-`-taxa_mapping_fp
		Taxonomy mapping filepath
	-f, `-`-fasta_fp
		Fasta sequence file.
	
	**[OPTIONAL]**
		
	-d, `-`-taxa_depth
		Depth of taxonomy to test for accuracy.  Depth that exceeds specifications in the taxonomy mapping file or report will be ignored [default: 3]
	-o, `-`-output_dir
		Specify output directory for reports, log. [default: .]
	-m, `-`-assignment_method
		 Taxonomic assignment method.  Currently only RDP classifier implemented. [default: rdp]
	-c, `-`-min_confidence
		Minimum confidence for taxonomic assignment.  [default: 0.8]
	-T, `-`-training_data_fp
		Training sequence data filepath for rdp classifier. [default: None]


**Output:**

 A log file from the taxonomy classifier will be generated, along with a taxonomic assignment file.  Additionally, a taxa_assignment_report.txt file will be generated detailing the accuracy of assignment for each level of taxonomic depth tested.


**Standard Example usage:**

::

	taxa_assignment_report.py [options] {-t taxa_mapping_filepath [required] -f input_fasta_filepath [required]}

**Change taxa depth to 4, change output directory to taxa_report/:**

::

	taxa_assignment_report.py -t taxa_mapping_filepath -f input_fasta_filepath  -d 4 -o taxa_report/ 


