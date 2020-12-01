.. _doc_install:
.. Primer Prospector documentation master file, created by William Walters (modified from QIIME install file, by Jesse Stombaugh)
   sphinx-quickstart on Mon Jan 25 12:57:02 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. index:: Installing Primer Prospector

============================
Installing Primer Prospector
============================


Primer Prospector consists of native code and additionally wraps many external applications. 

As a consequence of this 'pipeline' architecture, depending on the features of Primer Prospector that you plan to use, you may or may not need all of the Primer Prospector dependencies.

The following programs are needed to use all of the features of Primer Prospector.

You should follow instructions provided by the package developers to install the dependencies.

Dependencies required for Primer Prospector
-------------------------------------------

To install most of following dependencies you need to have a build environment on your machine. On OS X, this involves installing the `developer tools <http://developer.apple.com/technologies/xcode.html>`_. On Debian-based Linux (e.g., Ubuntu), this involves installing the ``build-essential`` package::

	sudo apt-get install build-essential

The following are required by Primer Prospector:

* Python 2.6 (`src <http://www.python.org/ftp/python/2.6.4/Python-2.6.4.tgz>`_) (license: PSF)
* PyCogent 1.5 (`src <http://sourceforge.net/projects/pycogent/files/PyCogent/1.5/PyCogent-1.5.tgz/download>`_) (license: GPL)
* Numpy 1.3.0 (`src <http://sourceforge.net/projects/numpy/files/NumPy/1.3.0/numpy-1.3.0.tar.gz/download>`_) (license: BSD)
* Matplotlib 0.98.5.3 (`src <http://iweb.dl.sourceforge.net/project/matplotlib/OldFiles/matplotlib-0.98.5.3.tar.gz>`_) (license: BSD)

To use the ``taxa_assignment_report.py`` module, it is necessary to install a taxonomy assignment program.  Currently, only the RDP classifier is implemented.

* rdp_classifier-2.0.1 (`src <http://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.0.1/rdp_classifier_2.0.1.tar.gz>`_) See :ref:`RDP install notes <rdp-install>`. (license: GPL)

To test for secondary structures between barcodes and primers with the ``check_primer_barcode_dimers.py``, we use the `Vienna <http://www.tbi.univie.ac.at/RNA/>`_ RNA secondary structure prediction software  with nearest-neighbor DNA energy values.

* Vienna 1.8.4 (`right click and select save as or download <http://www.tbi.univie.ac.at/~ivo/RNA/ViennaRNA-1.8.4.tar.gz>`_) Vienna package for predicting DNA/RNA secondary structure. (Free software for non-commercial use)

The DNA parameter file, ``dna_DM.par``, can be found in the ``DNA_parameters`` folder of Primer Prospector.  This file a modified form of the DNA parameters from David Mathews’ `RNAstructure program <http://rna.urmc.rochester.edu/RNAstructure.html>`_.

Useful source of reference sequences and taxonomy mapping files(currently only Bacteria and Archaea sequences):

* Greengenes 97% OTUs, taxonomies, and tree (`zip <http://greengenes.lbl.gov/Download/OTUs/gg_otus_6oct2010.zip>`_)

Silva (`src <http://www.arb-silva.de/>`_) also has sources of small subunit sequences, including eukaryotic sequences.  At the current time, however, there are not publicly hosted taxonomy mapping files in a format that is compatable with Primer Prospector.  The latest version of the aligned Silva SSU fasta files can be found here: (`src <http://www.arb-silva.de/typo3conf/ext/myth_repository/secure.php?u=0&file=fileadmin/silva_databases/release_104/Exports/SSURef_104_full_align_term.fasta.tgz&t=1290014689&hash=ac628445cf2b17f01e7c26414946a594>`_).  The `clean_fasta.py` module in Primer Prospector can be used to degap, remove spaces, and convert the uracil "U" to thymadine "T" characters, which is necessary for files such as the Silva fasta file listed above.

If you plan to build the Primer Prospector documentation locally:

* Sphinx 1.0.4 (`src <http://pypi.python.org/pypi/Sphinx>`_) (license: BSD)



License information for external dependencies
---------------------------------------------
We have attempted to provide accurate licensing information for the above dependencies for the convenience of our users. This information is by no means definitive and may contain errors. Any questions about licenses or the legality of specific uses of these software packages should be directed to the authors of the software. Do not rely solely on the license information presented above!

Shortcuts in this document
--------------------------
For simplicity throughout this document, we assume that you have downloaded Primer Prospector in ``/home/pprospector/``. You should consider all occurrences of ``/home/pprospector/`` in the remainder of this document as references to the directory which contains the Primer Prospector directory which you'll have after downloading and unpacking Primer Prospector.


Stable Release
^^^^^^^^^^^^^^^^^^
Currently the most stable version of Primer Prospector is our 1.0.1 release, which you can download from `here <http://sourceforge.net/projects/pprospector/files/pprospector-1.0.1.tar.gz/download>`_.

Latest Development Version
^^^^^^^^^^^^^^^^^^^^^^^^^^
To get the latest development version of Primer Prospector, you access our Sourceforge repository. While this code is subject to minor changes in interface, it will provide access to the latest and greatest features. The official web documentation is likely to be out-of-date with respect to the development software. You should instead refer to the svn documentation in pprospector/doc. Check out the latest version of Primer Prospector using svn with the command::

	svn co https://pprospector.svn.sourceforge.net/svnroot/pprospector/trunk pprospector

svn users should periodically update Primer Prospector by using the following command::

	svn update /home/pprospector/

Unpacking Primer Prospector (release only)
------------------------------------------
After downloading the Primer Prospector release tar file you'll need to unpack the code. For simplicity in this document, we will assume that you have downloaded Primer Prospector to the directory ``/home/pprospector/``. 

Unpack the release Primer Prospector tar file with the commands::

	cd /home/pprospector
	tar -xvzf pprospector-1.0.1.tar.gz

	
If you have downloaded from svn, Primer Prospector is already unpacked.
	
Installing Primer Prospector
----------------------------
Primer Prospector consists of library code (in ``pprospector/primerprospector``), test code (in ``pprospector/tests``), documentation (in ``pprospector/doc``), and scripts (in ``pprospector/scripts``). Installing Primer Prospector consists of running the tests (optional, but highly recommend), installing the library code in a place where python knows where to find it, and installing the scripts in a place where the shell looks for executable files.



Installing the library code and scripts with setup.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Using ``pprospector/setup.py`` (and thereby python's ``distutils`` package) is the recommended way of installing the Primer Prospector library code and scripts. You can optionally specify where the library code and scripts should be installed -- depending on your setup, you may want to do this. By default, the Primer Prospector library code will be placed under python's ``site-packages``, and the Primer Prospector scripts will be place in ``/usr/local/bin/``. You may need to run ``setup.py`` using ``sudo`` if you do not have permission to place files in the default locations. 

First, ensure that you are in the top-level QIIME directory::
	
	cd /home/pprospector/

By default the Primer Prospector scripts will be installed in ``/usr/local/bin``. As there are a lot of Primer Prospector scripts, we recommend customizing the script directory to keep your system organized. This can be customized with the ``--install_scripts`` option::
	
	python setup.py install --install-scripts=/home/pprospector/bin/
	
You can similarly install the library code in an alternate location using the ``--install-purelib`` option::
	
	python setup.py install --install-purelib=/home/pprospector/lib/


Combine these options as follows::
	
	python setup.py install --install-scripts=/home/pprospector/bin/ --install-purelib=/home/pprospector/lib/

For a complete discussion of customizations related to the setup.py script, `see this page <http://docs.python.org/install/index.html#alternate-installation-the-home-scheme>`_.

If you used default values for ``--install-scripts`` and ``--install-purelib`` (by not specifying them), your installation should be complete. If you specified an alternate value for ``--install-scripts``, you'll need to ensure that the shell knows where to look for the scripts. If you are using the bash shell and the locations specified in the examples above, you can do this with the following command::
	
	echo "export PATH=/home/pprospector/bin/:$PATH" >> /home/pprospector/.bashrc

If you specified an alternate value for ``--install-purelib``, you'll need to be sure that python knows where to look for Primer Prospector. If you are using the bash shell and the locations specified in the examples above, you can do this with the following command::
	
	echo "export PYTHONPATH=/home/pprospector/lib/:$PYTHONPATH" >> /home/pprospector/.bashrc
	
The source your ``.bashrc``::

	source /home/prospector/.bashrc


Running the test suite
----------------------
Next you should run the test suite. Execute the following commands::
	
	cd /home/pprospector/tests/
	python all_tests.py

You will see test output on the terminal indicating test successes and failures. Some failures are OK. The ``all_tests.py`` command will complete with a summary of test failures. Some tests may fail due to missing external applications -- these will be noted separately from other test failures. If these are related to features of Primer Prospector that you are not using, this is acceptable. Otherwise, you'll need to ensure that you have the external applications installed correctly (and the correct versions), and re-run the tests. 

Testing your Primer Prospector installation
-------------------------------------------
If Primer Prospector is installed correctly, you should be able to run the Primer Prospector scripts. Try the following::
	
	cd
	analyze_primers.py -h
	
This should give you help text describing the interface to the analyze_primers.py script. (Note that if you do not have a /home/pprospector/.bashrc you may get an error at the ``source`` step. If you did not specify alternate values for ``--install-purelib`` or ``--install-scripts`` this shouldn't be a problem.)

External application install notes
----------------------------------

PATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^

External applications used by Primer Prospector need to be visible to the shell by existing in executable search path (i.e., listed in the ``$PATH`` environment variable). For example, if you plan to use the RDP classifier, and have the RDP classifier executables installed in ``/home/pprospector/bin`` you can add this directory to your system path with the commands::
	
	echo "export PATH=/home/pprospector/bin/:$PATH" >> /home/pprospector/.bashrc
	source /home/pprospector/.bashrc

PYTHONPATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Primer Prospector, PyCogent, Matplotlib, and NumPy must be visible to python for all features of Primer Prospector.  If you have used these, you should not need to modify your PYTHONPATH to make the library code visible. If you haven't used the respective setup.py scripts, or if you specified an alternate value for ``--install-purelib``, you may need to add the locations of these libraries to your PYTHONPATH environment variable. 

For example, if you've installed RDP in ``/home/pprospector/RDP`` you can add this to your PYTHONPATH with the commands::
	
	echo "export PYTHONPATH=/home/pprospector/RDP/:$PYTHONPATH" >> /home/pprospector/.bashrc
	source /home/pprospector/.bashrc


RDP_JAR_PATH Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _rdp-install:

If you plan to use the RDP classifier for taxonomy assessment reports, you must also define an RDP_JAR_PATH variable. If you have the RDP classifier jar file (``rdp_classifier-2.0.1.jar``) in ``/home/rdp/app`` you can do this with the following command::

	echo "export RDP_JAR_PATH=/home/rdp/app/rdp_classifier-2.0.1.jar" >> /home/rdp/.bashrc
	source /home/rdp/.bashrc
	
Building The Primer Prospector Documentation
--------------------------------------------

.. _build-primer-prospectorr-docs:

If you are using the svn version of Primer Prospector, you may want to build the documentation locally for access to the latest version. You can change to the ``pprospector/doc`` directory and run (you may need to preceed the command with `sudo`)::

	make html
	
We try to update the documentation as we update the code, but svn users may notice some discrepancies. After building the documentation, you can view it in a web browser by opening the file ``pprospector/doc/_build/html/index.html``. You may want to bookmark that page for easy access. 
