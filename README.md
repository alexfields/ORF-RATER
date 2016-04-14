# ORF-RATER
ORF-RATER (Open Reading Frame - Regression Algorithm for Translational Evaluation of Ribosome-protected footprints) comprises a series of scripts for coding sequence annotation based on ribosome profiling data.

The software was created at [Jonathan Weissman's lab at UCSF](http://weissmanlab.ucsf.edu/) and is described in [Fields, Rodriguez, et al., "A regression-based analysis of ribosome-profiling data reveals a conserved complexity to mammalian translation", *Molecular Cell* **60**, 816-827 (2015).](http://dx.doi.org/10.1016/j.molcel.2015.11.013)

Usage information can be found in the Detailed Protocol included in the paper's supplemental materials, or by running each script with the --help/-h flag.

Required packages include [numpy](http://www.numpy.org), [scipy](http://www.scipy.org), [pysam](https://github.com/pysam-developers/pysam), [biopython](http://www.biopython.org), [pandas](http://pandas.pydata.org/), [tables](http://http://www.pytables.org/), [scikit-learn](http://scikit-learn.org/), [pybedtools](https://pythonhosted.org/pybedtools/), and [plastid](https://pypi.python.org/pypi/plastid), all of which are available through [PyPI](https://pypi.python.org/pypi).

Some features require the [multiisotonic](https://github.com/alexfields/multiisotonic) package, which must be downloaded manually. Multiisotonic additionally requires [python-igraph](https://github.com/igraph/python-igraph).

Contact Alex Fields for further information or assistance.
