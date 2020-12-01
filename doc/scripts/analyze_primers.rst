.. _analyze_primers:

.. index:: analyze_primers.py

*analyze_primers.py* --  Score primers for binding to a set of target sequences. 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**



This module performs a local alignment for the input primer(s) against target 
sequences, and then determines a weighted score based upon the number of 
mismatches and gaps.  A summary graph showing mismatches and weighted scores is 
generated for all input fasta files, as is a hits file containing details about 
the primer mismatches, index in the sequence, and other details about primer 
binding.

This module takes an input primer file and one or more fasta files.  Each
primer is tested against every sequence to find the best local alignment.
Mismatches and gaps are calculated for the primer, along with a weighting score 
which gives larger penalties to gaps and mismatches in the 3' end of the primer.

An output hits file is generated for each primer, recording information about
the primer hit site, mismatches, and overall weighted score (a perfect score
starts at zero and increases as penalties are added).  A graph is
also generated, showing mismatches/gaps and overall score information for the 
primer and the target sequences.

The primers input file should be generated with the following format:
Comments are preceeded by a pound "#" symbol.
The primer data are tab delineated with the primer name first, such as
"349_v2r", the actual nucleotide sequence next, listed in a 5' to 3' sequence, 
(example: "AATCGRACGNTYA"), and finally a comment 
or citation, if any, can be listed.  Forward primers should be followed by a 
"f" while reverse primers are followed by "r".
A complete example line is listed below.

815_v34f	GTGGCCNATRRCYAGAACGC	Darrow,Scopes,Bryan et al. 1926

The input sequences should be in fasta format.  If more than one file is 
supplied, they should be separated by a colon.  Each fasta file passed will
have its sequence coverage displayed in separate output graphics and hits
files.



**Usage:** :file:`analyze_primers.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-fasta_seqs
		Target fasta file(s) to score input primer(s) against. Separate multiple files with a colon.
	
	**[OPTIONAL]**
		
	-P, `-`-primers_filepath
		Path to input primers file.  This file is tab delineated, with the first column being the primer name, which must end with an 'f' or a 'r'.  The second column contains the primer sequence in 5' to  3' format. One must supply either a primer file or a primer name  (-p parameter) and primer sequence (-s parameter). [default: None]
	-p, `-`-primer_name
		Specify a single primer to analyze.  One can either specify a single primer that is listed in a primers file (-P parameter) or  specify a sequence with the -s parameter.  Primer name must end with a "f" or "r" to designate orientation. [default: None]
	-s, `-`-primer_sequence
		Primer sequence if using the -p option.  Must be specified if not passing a primers file with the -P option.  If both -P and -p parameters are passed, the sequence given with this option will be taken rather than sequences in the -P primers file. [default: None]
	-o, `-`-output_dir
		Specify output directory for hits files and primer summary graphs. [default: .]
	-e, `-`-three_prime_len
		Length of primer considered to be part of the 3' region for the purpose of giving a weighted score for mismatches and/or gaps. [default: 5]
	-l, `-`-last_base_mismatch
		Sets penalty for mismatch in final base of 3' end of the primer. [default: 3]
	-t, `-`-three_prime_mismatch
		Penalty for all 3' mismatches except final base.[default: 1]
	-T, `-`-non_three_prime_mismatch
		Penalty for all non-3' mismatches. [default: 0.4]
	-g, `-`-three_prime_gap
		Penalty for gaps in the 3' region of the primer.  [default: 3]
	-G, `-`-non_three_prime_gap
		Penalty for non 3' gaps. [default: 1]


**Output:**

For each primer tested, an output _hits.txt file containing information about the primer hit to each sequence is generated, as well a .ps file showing overviews for the mismatches and weighted score for the primer and target sequences.  Both output files are named by the primer and fasta file tested


**Example:**

Standard Example:

::

	analyze_primers.py [options] {-P input_primers_filepath [required] -f input_fasta_filepath [required]}

**Manually specify a primer name and sequence:**

Note - primer name must end with 'f' or 'r':

::

	analyze_primers.py -p "primer_name_f" -s "ACCTGACRGGTAATC" -f input_fasta_filepath

**Use multiple target files, change scoring parameters:**

Pass a primers file, two target fasta files, change the size of the 3' region from the default 5 bases to 7 bases, and lower the 3' mismatch penalty from the default 1.0 to 0.6:

::

	analyze_primers.py -P primers.txt -f bacterial_seqs.fasta:eukaryotic_seqs.fasta -e 7 -t 0.6




