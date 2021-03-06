.. _clean_fasta:

.. index:: clean_fasta.py

*clean_fasta.py* -- Filters input fasta file to remove gaps and 'U' characters.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This module will filter out gap characters ('.' and '-'), spaces, and/or uracil characters ('U') from an input fasta file(s).  The module can also capitalize the characters in the filtered file(s).  This should be done for fasta files to be utilized with `analyze_primers.py <./analyze_primers.html>`_ or any downstream module requring input fasta files. 


**Usage:** :file:`clean_fasta.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-fasta_seqs
		Target fasta file(s) to filter. Separate multiple files with a colon.
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Output directory for filtered fasta file(s). [default: .]
	-g, `-`-gap_chars
		Filter gap characters "." and "-" [default: True]
	-s, `-`-space_chars
		Filter space characters. [default: True]
	-u, `-`-convert_uracil
		Convert Uracil "U" characters to "T". [default: True]
	-c, `-`-cap_seqs
		Capitalize filtered sequences [default: True]


**Output:**

Output fasta file(s) generated by this module will have a '_filtered.fasta' appended to the original fasta file name.


**Standard Example usage:**

::

	clean_fasta.py [options] {-f input_fasta_filepath [required] }

**Keep gap characters in two aligned fasta files, but allow all other filtering:**

::

	clean_fasta.py -f input_aligned_seqs1.fasta:input_aligned_seqs2.fasta -g


