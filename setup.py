#!/usr/bin/env python
# File created on 17 Feb 2010, modified 11-16-2010 by W. Walters
from __future__ import division
from __future__ import print_function
from distutils.core import setup
from os import chdir, getcwd, listdir
from os.path import join, abspath
from subprocess import call
from glob import glob
import re

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The Primer Prospector project"
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "William Walters"]
__license__ = "GPL"
__version__ = "1.0.1-release"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Release"
 
 
long_description = """Primer Prospector: de novo design and taxonomic analysis of PCR primers.
"""

doc_imports_failed = False
try:
    import sphinx
except ImportError:
    doc_imports_failed = True

def build_html():
    """ Build the sphinx documentation 
    
    The code for building sphinx documentation is based on 
    PyCogent's setup.py.
    
    """
    cwd = getcwd()
    doc_dir = join(cwd,'doc')
    chdir(doc_dir)
    call(["make", "html"])
    chdir(cwd)
    index_html_path = join(abspath(doc_dir),'_build','html','index.html')
    print("Local documentation built with Sphinx. "+\
          "Open to following path with a web browser:\n%s" %\
            index_html_path)

# try:
#     import cogent3
# except ImportError:
#     print("cogent3 not installed but required. (Is it installed? Is it in the current user's $PYTHONPATH or site-packages?) See http://pycogent.sourceforge.net.")
#     exit(1)

# pycogent_version = tuple([int(v) \
#         for v in re.split("[^\d]", cogent3.__version__) if v.isdigit()])
        
# if pycogent_version < (1,4):
#     print("PyCogent >= 1.4.0 required, but %s is installed." % cogent.__version__)
#     exit(1)
    
setup(name='PrimerProspector',
      version=__version__,
      description='Primer Prospector',
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='http://pprospector.sourceforge.net/',
      packages=['primerprospector', 'primerprospector.cogentutil'],
      scripts=glob('scripts/*py'),
      #package_data = {'primerprospector': ['']},
      long_description=long_description
)

if doc_imports_failed:
    print("Sphinx not installed, so cannot build local html documentation.")
else:
    build_html()
    
'''package_data={'qiime':\
                   ['support_files/qiime_config',\
                    'support_files/css/*css',\
                    'support_files/html_templates/*html',\
                    'support_files/images/*png',\
                    'support_files/jar/*jar',\
                    'support_files/js/*js',\
                    'support_files/sra_xml_templates/*xml']},'''
