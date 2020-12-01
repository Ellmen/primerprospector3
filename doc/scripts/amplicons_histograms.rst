.. _amplicons_histograms:

.. index:: amplicons_histograms.py

*amplicons_histograms.py* -- Create histograms of predicted amplicons for primer pairs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Using amplicons files (generated with `get_amplicons_and_reads.py <./get_amplicons_and_reads.html>`_), this module will generate histogram(s) showing the predicted amplicon sizes.  If a taxonomy mapping file is passed, the histograms will be separated according to domain (archaea, bacteria, and eukaryotic).


**Usage:** :file:`amplicons_histograms.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-amplicons_filepath
		Target amplicons files.  Separate multiple files with a colon.
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Specify output directory for amplicons and reads. [default: .]
	-r, `-`-all_files
		Generate histograms for all files ending with _amplicons.fasta in the directory specified with the -f parameter [default: False]
	-t, `-`-taxa_map
		Filepath to taxonomy mapping file, used to separate graphs according to domain. [default: None]


**Output:**

There will be an output graph in postscript format with the same root name as the input _amplicons.fasta file for each _amplicons.fasta file passed.


Standard Example usage (create a histogram from an amplicon file, output to the current directory):

::

	amplicons_histograms.py [options] {-f amplicons_filepath }

Test for all _amplicons.fasta files in current directory, pass a taxonomy mapping file so histograms are plotted according to domain, and output to amplicons_graph directory:

::

	amplicons_histograms.py -f . -r -t taxonomy_mapping.txt -o amplicons_graph


