# ORF-RATER
ORF-RATER (Open Reading Frame - Regression Algorithm for Translational Evaluation of Ribosome-protected footprints) comprises a series of scripts for coding sequence annotation based on ribosome profiling data.

The software was created at [Jonathan Weissman's lab at UCSF](http://weissmanlab.ucsf.edu/) and is described in [Fields, Rodriguez, et al., "A regression-based analysis of ribosome-profiling data reveals a conserved complexity to mammalian translation", *Molecular Cell* **60**, 816-827 (2015).](http://dx.doi.org/10.1016/j.molcel.2015.11.013)

Usage information can be found in the Detailed Protocol included in the paper's supplemental materials, or by running each script with the --help/-h flag.

Required packages include [numpy](http://www.numpy.org), [scipy](http://www.scipy.org), [pysam](https://github.com/pysam-developers/pysam), [biopython](http://www.biopython.org), [pandas](http://pandas.pydata.org/), [tables](http://www.pytables.org/), [scikit-learn](http://scikit-learn.org/), [pybedtools](https://pythonhosted.org/pybedtools/), and [plastid](https://pypi.python.org/pypi/plastid), all of which are available through [PyPI](https://pypi.python.org/pypi).

Some features require the [multiisotonic](https://github.com/alexfields/multiisotonic) package, which must be downloaded manually. Multiisotonic additionally requires [python-igraph](https://github.com/igraph/python-igraph).

Transcripts must be presented in UCSC's [BED12 format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). The most reliable method I've found to convert from GTF to BED12 involves first converting to [genePred format](https://genome.ucsc.edu/FAQ/FAQformat.html#format9), making use of UCSC's "gtfToGenePred" and "genePredToBed" scripts, which are available [here](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/). The full command is `gtfToGenePred INPUT_GTFFILE.gtf stdout | genePredToBed stdin OUTPUT_BEDFILE.bed`. Similarly, a BED file can be converted to a GTF using the command `bedToGenePred INPUT_BEDFILE.bed stdout | genePredToGtf file stdin OUTPUT_GTFFILE.gtf`.

Contact Alex Fields for further information or assistance.
