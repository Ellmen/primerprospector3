.. _generate_linkers:

.. index:: generate_linkers.py

*generate_linkers.py* --  Find the best linkers to use for the given primers 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**



generate_linkers is designed to suggest primer linkers based upon the least frequently occuring base pairs (default 2) immediately upstream from the 5' region of the primer.

A single _hits.txt file, multiple hits files, or a directory of hits.txt 
files can be passed as input (see the `analyze_primers.py <./analyze_primers.html>`_ file for 
information about the hits.txt file format).  The fasta files used to 
generate the _hits.txt files are also required input for this module.

The output file will contain the percentage occurance of each base at 
all positions in the linker and a suggested linker based on complementarity 
to the least frequently occurring bases.



**Usage:** :file:`generate_linkers.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-hits_fps
		Target primer hits files to generate linkers with.  Separate multiple files with a colon.
	-f, `-`-fasta_fps
		Fasta filepath(s).  Must include all fasta sequences used to generate the hits files.  Separate multiple files with a colon.
	
	**[OPTIONAL]**
		
	-l, `-`-linker_len
		Size of linker in base pairs. [default: 2]
	-r, `-`-all_files
		Test all _hits.txt files in directory specified with -i.   [default: False]
	-o, `-`-output_dir
		Specify output directory for linkers summary. [default: .]
	-s, `-`-score_type
		Value to use from primer hits file to determine a given primer's amplification success.  Valid choices are weighted_score, overall_mismatches, tp_mismatches.  Gibbs energy scores not currently implemented [default: weighted_score]
	-t, `-`-score_threshold
		If primer has score at or below this parameter, the primer amplification is considered to be successful. [default: 1.0]


**Output:**

A text file will be generated for each hits file containing details about the base frequency for each position of the linker and the suggested linker sequence.


**Standard Example usage:**

::

	generate_linkers.py [options] {-i hits_filepath [required] -f input_fasta_filepath [required]}

**Change linker size to 3, manually pass 2 hits files to test, use overall mismatches for scoring:**

::

	generate_linkers.py -i primer1f_hits.txt:primer2r_hits.txt -f seqs.fasta -l 3 -s overall_mismatches 

**Test all hits file in a target directory, put output in linkers directory:**

::

	generate_linkers.py -i hits_dir/ -f seqs.fasta -o linkers/


