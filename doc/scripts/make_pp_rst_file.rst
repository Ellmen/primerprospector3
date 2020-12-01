.. _make_pp_rst_file:

.. index:: make_pp_rst_file.py

*make_pp_rst_file.py* -- Make Sphinx RST file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script will take a script file and convert the usage strings and options to generate a documentation .rst file.


**Usage:** :file:`make_pp_rst_file.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_script
		This is the input script for which to  make a .rst file
	-o, `-`-output_dir
		Output directory for .rst files


**Output:**

This will output a Sphinx rst-formatted file.


**Example:**

::

	make_qiime_rst_file.py -i make_2d_plots.py -o doc/


