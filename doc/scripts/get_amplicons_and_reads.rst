.. _get_amplicons_and_reads:

.. index:: get_amplicons_and_reads.py

*get_amplicons_and_reads.py* -- Generates predicted amplicons and reads from primers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Using hits files (generated with `analyze_primers.py <./analyze_primers.html>`_), this module will look at an individual primer or a specified primer pair to generate amplicons and reads.  By default, the weighted score is used to determine if a particular primer will amplify, but other results (3' mismatches, overall mismatches) can be used to determine if a primer will amplify.  Every fasta file used to generate the hits file must be passed to this module (multiple fastas should be separated by a colon).  If a single hits file is specified, the amplicon will be the entire sequence following the 3' end of the primer.


**Usage:** :file:`get_amplicons_and_reads.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-primer_hits
		Target primer hits files.  Separate multiple files with a colon.
	-f, `-`-fasta_fps
		Fasta filepaths.  Must match the fasta files used in the analyze_primers module.  Multiple fasta files can be passed, separated with a colon.  Order not important.
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Specify output directory for amplicons and reads. [default: .]
	-s, `-`-score_type
		Value to use from primer hits file to determine a given primer's amplification success.  Valid choices are weighted_score, overall_mismatches, tp_mismatches.  Gibbs energy scores not currently implemented [default: weighted_score]
	-t, `-`-score_threshold
		If primer has score at or below this parameter, the primer amplification is considered to be successful [default: 1.0]
	-m, `-`-min_seq_len
		Sets minimum sequence length of amplicon to be included in the output amplicon file [default: 100]
	-d, `-`-read_direction
		Direction of reads generated. Can be forward (f), reverse (r), or paired end (p).  [default: r]
	-R, `-`-read_len
		Length of reads to generate.  Should be set according to sequencing technology/reagents used.  [default: 250]


**Output:**

An amplicon file (in fasta format) will be written for the forward and reverse primer pair, with a name based upon the name of the primers (taken from hits file name, e.g. 27f_338r_amplicons.fasta).  The reads file will be generated in the same manner, with an "_reads" instead of "_amplicons".


**Standard Example usage (examine pair of primer hits files):**

::

	get_amplicons_and_reads.py [options] {-f fasta_filepath -i primer_hits_filepath1:primer_hits_filepath2 [required] }

**Use total mismatches instead of weighted score to determine if a given primer will amplify:**

::

	get_amplicons_and_reads.py -f fasta_filepath -i primer_hits_filepath1:primer_hits_filepath2 -s total_mismatches




