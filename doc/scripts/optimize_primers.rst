.. _optimize_primers:

.. index:: optimize_primers.py

*optimize_primers.py* --  Records base frequency for each position of a primer for optimization by modulating degeneracy. 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**



Optimize primers takes an input primer hits file (created by `analyze_primers.py <./analyze_primers.html>`_)
and generates an output file containing the base frequency for each position 
in the primer.

Only primer hits that exceed a given score threshold will be considered for
the base frequency output.



**Usage:** :file:`optimize_primers.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-hits_fp
		Target primer hits file to generate base frequencies with.
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Specify output directory for linkers summary. [default: .]
	-s, `-`-score_type
		Value to use from primer hits file to determine a given primer's amplification success.  Valid choices are weighted_score, overall_mismatches, tp_mismatches.  Gibbs energy scores not currently implemented [default: weighted_score]
	-t, `-`-score_threshold
		If primer has score at or below this parameter, the primer amplification is considered to be successful. [default: 2.0]


**Output:**

A tab-seperated text file will be generated
containing the primer sequence, in 5' to 3' orientation, and the base 
frequencies for each nucleotide. 


**Standard Example usage:**

::

	optimize_primers.py [options] {-i hits_filepath [required]}

**Use overall mismatches for scoring, set score threshold to only consider primer hits that are equal to or less than two mismatches:**

::

	optimize_primers.py -i primer2r_hits.txt -s overall_mismatches -t 2


