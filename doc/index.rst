
.. Primer Prospector documentation master file, created by
   sphinx-quickstart on Mon Jan 25 12:57:02 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#################
Primer Prospector
#################

Primer Prospector is a pipeline of programs to design and analyze PCR primers.

`De novo` primers can be generated from a set of aligned fasta sequences and these primers, or known primers from literature, can be scored against target sequences.  Predicted primer performance can be visualized in summary graphs showing overall primer performance or across specific taxonomies.

Primer Prospector is built in Python using the open-source PyCogent_ toolkit.


Download Primer Prospector
==========================

 * Stable Release: The current release is 1.0.1, which you can `download here <http://sourceforge.net/projects/pprospector/files/pprospector-1.0.1.tar.gz/download>`_.

 * Development Version: Primer Prospector is under development.  To get the latest development version of Primer Prospector, you access our Sourceforge repository. While this code is subject to minor changes in interface, it will provide access to the latest and greatest features. The official web documentation is likely to be out-of-date with respect to the development software. You should instead refer to the svn documentation in pprospector/doc. Check out the latest version of Primer Prospector using svn with the command::

	svn co https://pprospector.svn.sourceforge.net/svnroot/pprospector/trunk pprospector

Installing and using Primer Prospector
======================================
New users should begin with the `Primer Prospector installation guide <./install/install.html>`_. After installing Primer Prospector, you should move on to the `Primer Prospector overview tutorial <./tutorial/tutorial.html>`_ to generate and analyze several primers. For general information about Primer Prospector, you'll want to refer to the `Primer Prospector documentation <./documentation/index.html>`_ 

Contact Us
===========

The primary contact for Primer Prospector is William Walters (use @ for AT and . for DOT:  williamDOTaDOTwaltersATcoloradoDOTedu).

Users can submit `bug reports <http://sourceforge.net/tracker/?group_id=318992&atid=1341293>`_ and `feature requests <http://sourceforge.net/tracker/?group_id=318992&atid=1341296>`_ via Sourceforge


Citing Primer Prospector
========================
If you use Primer Prospector for any published research, please include the following citation:

	* `PrimerProspector: de novo design and taxonomic analysis of PCR primers.` William A. Walters (1,6), J. Gregory Caporaso (2,6), Christian L. Lauber (3), Donna Berg-Lyons (3), Noah Fierer (4), and Rob Knight (2,5).  Submitted (as of December 2010)
	1. Department of Molecular, Cellular, and Developmental Biology.  University of Colorado at Boulder, Boulder, CO, USA.
	2. Department of Chemistry and Biochemistry.  University of Colorado at Boulder, Boulder, CO, USA.
	3. Cooperative Institute for Research in Environmental Sciences.  University of Colorado at Boulder, Boulder, CO, USA.
	4. Department of Ecology and Evolutionary Biology.  University of Colorado at Boulder, Boulder, CO, USA.
	5. Howard Hughes Medial Institute, Boulder, CO, USA.
	6. These authors contributed equally to this work.
	

.. _PyCogent: http://pycogent.sourceforge.net


