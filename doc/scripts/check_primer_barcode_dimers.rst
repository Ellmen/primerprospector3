.. _check_primer_barcode_dimers:

.. index:: check_primer_barcode_dimers.py

*check_primer_barcode_dimers.py* --  Tests input combinations of barcodes and primers for potential secondary structure.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**



This module takes an input text file containing barcodes, line separated, and
listed in 5' to 3' direction.  Two primers are specified at the command line.
The first, primer 1, is the primer sequence (and linker) that is associated
with a barcode.  Primer 2 by default has no barcode associated with it, 
although barcode association with both primers can be enabled.

Each barcode/primer combination that potentially will form secondary 
structures will be listed with line number, barcode sequence, secondary 
structure, and energy in an output text file.  Additionally, a postscript file
showing graphical output for the secondary structure will be generated in a 
subdirectory of the output directory.  Degenerate primers will be listed in 
the form(s) that are likely to have secondary structure.



**Usage:** :file:`check_primer_barcode_dimers.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-b, `-`-barcodes
		Filepath of barcodes to score input primer(s) against.
	-p, `-`-primer1
		Primer, written in 5' to 3', that is linked to barcodes tested.  If linker sequence is present between primer and barcode, include it with this sequences.
	-P, `-`-primer2
		Second primer, written in 5' to 3' orientation.  This primer by default is not associated with any barcodes.
	-e, `-`-energy_parameters
		Specify energy parameters file for predicting secondary structures.  A DNA parameters file, dna_DM.par, is found in the DNA_parameters folder of Primer Prospector, and should be pointed to with this parameter.  If an incorrect file is used, the Vienna software will use default parameters, which are for RNA folding, and could give misleading results.  The provided DNA parameters file is a modified form of the DNA parameters from  David Mathews' RNAstructure program.
	
	**[OPTIONAL]**
		
	-t, `-`-annealing_temp
		Specify an annealing temperature in degrees Celsius. [default: 50]
	-s, `-`-score_threshold
		Specify a score threshold for the Gibbs energy calculation, below which a barcode/primer combination is flagged for potential secondary structure.  [default: -10.0]
	-o, `-`-output_dir
		Specify output directory for barcode/primer secondary structure summary and graphs. [default: .]
	-B, `-`-paired_end_barcodes
		If enabled, barcodes will be appended to both primer 1 and primer 2.  [default: False]
	-g, `-`-suppress_graphs
		Suppress retention of output postscript graphs. [default: False]


**Output:**

Barcode/primer combinations that are likely to form secondary structures will be listed in the output 'flagged_barcodes.txt' file.  Any such barcode/primer combination will also be present in the graphs/ directory, and named for the line number that the barcode was found on and the barcode sequence.


**Example:**

Standard Example:

::

	check_primer_barcode_dimers.py [options] {-b barcode_filepath [required] -p primer_sequence_1 [required] -P primer_sequence_2 [required] -e DNA_energy_parameters_filepath [required]}

**Set the annealing temperature to 65, set score threshold to -10, and point to the energy parameters file in the directory /home/BigErn/pprospector/DNA_parameters:**

::

	check_primer_barcode_dimers.py -b barcodes.txt -p "ACCTGACRGGTAATC" -P "CGACCAGTACG" -t 65 -s -10 -e /home/BigErn/pprospector/DNA_parameters/dna_DM.par


