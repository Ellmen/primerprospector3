# Primer Prospector (3)

This is a clone of the excellent Primer Prospector developed at the University of Colorado, modified as minimally as possible to support Python 3. You can read more about primer prospector here: http://pprospector.sourceforge.net/

Because pycogent is heavily dependent on Python 2, I had to replace the dependency to cogent3 which has a slightly different API.

## Current status

So far I have all functionality in the tutorial (http://pprospector.sourceforge.net/tutorial/tutorial.html) working up to "Assessing Taxonomic Usefulness of Reads". The tests are failing because they have the heaviest dependency on old pycogent utils which I am still working to replace.

## Installation

Clone the repository, install cogent3 (`pip install cogent3`), and run `pip install .`
