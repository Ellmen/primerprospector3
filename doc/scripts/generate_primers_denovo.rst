.. _generate_primers_denovo:

.. index:: generate_primers_denovo.py

*generate_primers_denovo.py* -- Find conserved and/or specific DNA sequences for use as the 3' end of a primer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The generate_primers_denovo module is designed to take an input sequence length (default 5) and aligned fasta file(s) to search for all Xmers of the given sequence length that are conserved in the target fasta files.  Optionally, any Xmers that are found to exist above a certain threshold (1% is the default) in the excluded fasta sequences are discarded.  The remaining Xmers, along with their upstream and downstream sequences are written to an output file.


**Usage:** :file:`generate_primers_denovo.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-target_seqs
		Target aligned fasta sequence files to find conserved sites for primer design.  Separate multiple files with a colon.
	-o, `-`-output_filepath
		Name of output filepath to write details about conserved sequence sites.
	
	**[OPTIONAL]**
		
	-e, `-`-exclude_fasta
		Excluded aligned fasta file(s).  To pass multiple files, separate each file with a colon.  Example: -e test1.fasta:test2.fasta.  If not specified, will skip exclusion step [default: None]
	-p, `-`-percent_match
		Percentage of sequence matches to primer that must match in order to retain prospective sequence in dictionary. [default: 0.6]
	-s, `-`-full_primer_length
		Overall primer length to retrieve from sequences. [default: 20]
	-x, `-`-xmer_length
		Xmer length to search for in target fasta sequence(s). [default: 5]
	-S, `-`-specificity_threshold
		Sets specificity threshold for excluded fasta sequences. [default: 0.01]
	-l, `-`-log_file
		Log filepath. If not specified, no log file will be written.  [default: None]
	-a, `-`-standard_index_file
		Aligned sequence file with which to assign prospective primer indices to.  The alignment where a conserved sequence is found will be used to determine the unaligned index in the supplied file (for instance an E. coli sequence) and will be recorded in the output file for the purpose of giving a meaningful name to prospective primers.  Only the first sequence in the file will be used for determining an index [default: None]
	-r, `-`-search_range
		Range of nucleotides in the supplied aligned target sequences to search for primers.  Supply the starting index and end index separated by a colon.  Example -r 1500:2700  Enable this option to generate primers that target certain regions. [default: None]


**Output:**

The output file from generate_primers_denovo is a text file containing information about the target Xmers that met the sensitivity (and specificity if given) threshold(s).  For each Xmer, the calculated sensitivity and specificity values are given, as well as the upstream and downstream sequences from the Xmer for all of the sequences with perfect matches.


**Standard Example usage:**

::

	generate_primers_denovo.py [options] {-i include_fasta_filepath(s) -o output_primers_filepath}

**Look for common 5mers in a given aligned fasta file, bact_sample.fasta, that are conserved at least 60% (default setting) of the time.  Output to primers.txt:**

::

	generate_primers_denovo.py -i bact_sample.fasta -o primers.txt

**Look for common 6mers in two input aligned fasta files (bact_sample.fasta, arch_sample.fasta), that are not found in the aligned fasta file euk_sample.fasta:**

::

	generate_primers_denovo.py -i bact_sample:arch_sample.fasta -x 6 -e euk_sample.fasta -o primers.txt


