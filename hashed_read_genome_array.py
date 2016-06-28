from plastid.genomics.genome_array import BAMGenomeArray
import re
import numpy
import itertools

# Regular expressions for identifying the number of mismatches at the 5' end of a read
NMIS_RE_PLUS = re.compile(r'^(0[ACGTN])*')
NMIS_RE_MINUS = re.compile(r'([ACGTN]0)*$')


def read_length_nmis(read):
    """Get the 5'-mismatch-trimmed read length and number of mismatches trimmed
    Can be used as a read_key_fun

    Parameters
    ----------
    read : :py:class:`pysam.AlignedSegment`
        Read to assign a key

    Returns
    -------
    rdlen : int
        Length of trimmed read

    nmis : int
        Number of 5' mismatches trimmed
    """
    nmis = len((NMIS_RE_MINUS if read.is_reverse else NMIS_RE_PLUS).search(read.opt('MD')).group())/2
    return len(read.positions)-nmis, nmis


def ReadKeyMapFactory(offset_dict, read_key_fun=lambda read: len(read.positions)):
    """Returns a mapping function for HashedReadBAMGenomeNDArray. Reads are mapped
    at a specified offset from the 5' end of the alignment, which can vary with
    the value returned by read_key_fun (e.g. with position or sequence). This
    function is NOT compatible with BAMGenomeArray.

    Parameters
    ----------
    offset_dict : dict
        Dictionary mapping read keys (from read_key_fun) to offsets from the 5' end
        of reads. Reads whose keys are not in offset_dict are silently dropped

    read_key_fun : function
        Function to assign appropriate key for each read.
        Assigned keys should match those in offset_dict

    Returns
    -------
    function
        Mapping function
    """

    # docstring of function we will return.
    docstring = """Returns reads covering a region, and a dict of count vectors mapping
        reads to specific positions in the region, mapping reads at possibly varying
        offsets from the 5' end of each read.

        Parameters
        ----------
        seg : |GenomicSegment|
            Region of interest

        reads : list of :py:class:`pysam.AlignedSegment`
            Reads to map

        Returns
        -------
        dict<:py:class:`numpy.ndarray`>
            Dict of count vector at each position in `seg`, keyed according to
            self.read_key_fun
        """

    def map_func(reads,seg):
        # reads_out = []
        count_array = {k: numpy.zeros(len(seg)) for k in offset_dict}
        for read in reads:
            read_key = read_key_fun(read)
            if read_key in offset_dict:
                offset = offset_dict[read_key]
                if not read.is_reverse:
                    p_site = read.positions[offset]
                else:
                    p_site = read.positions[-offset - 1]

                if p_site >= seg.start and p_site < seg.end:
                    # reads_out.append(read)
                    count_array[read_key][p_site - seg.start] += 1
        return count_array
        # return reads_out, count_array
    map_func.read_keys = offset_dict.keys()
    map_func.read_key_fun = read_key_fun
    map_func.__doc__ = docstring
    return map_func


class HashedReadBAMGenomeArray(BAMGenomeArray):
    """Identical to a BAMGenomeArray, but keeps reads classified by a hash
    function (e.g. read length)
    """

    def __init__(self,bamfiles,mapping):
        """Create HashedReadBAMGenomeArray

        Parameters
        ----------
        bamfile : list
            An list of open :py:class:`pysam.AlignmentFile` s. Note: the
            corresponding `BAM`_ files must be sorted and indexed by `samtools`_.

        mapping : func
            Mapping function that determines how each read alignment is mapped to a
            count at a genomic position. Returns a list of pysam.AlignedReads and a
            dict of count vectors corresponding to the number of mapped read counts
            at each position, keyed according to a supplied function. Must have a
            list of valid keys stored as mapping.read_keys.
            Typically generated using ReadKeyMapFactory().
        """
        BAMGenomeArray.__init__(self,bamfiles,mapping=mapping)

    def get_reads_and_hashed_counts(self,roi,roi_order=True):
        """Returns reads covering a region, and a dict of count vectors mapping
        reads to specific positions in the region, following rules specified by
        self.map_fn

        Parameters
        ----------
        roi : |GenomicSegment|
            Region of interest


        Returns
        -------
        list
            List of reads (as :class:`pysam.AlignedSegment`)
            covering region of interest

        dict<numpy.ndarray>
            Counts at each position of `roi`, keyed according to the map_fn
        """
        # fetch all reads
        if roi.chrom not in self.chroms():
            return [], {k: numpy.zeros(len(roi)) for k in self.map_fn.read_keys}
        reads = itertools.chain.from_iterable((X.fetch(reference=roi.chrom,
                                                       start=roi.start,
                                                       end=roi.end) for X in self.bamfiles))

        # filter by strand
        if roi.strand == "+":
            reads = itertools.ifilter(lambda x: not x.is_reverse, reads)
        elif roi.strand == "-":
            reads = itertools.ifilter(lambda x: x.is_reverse, reads)

        # Pass through additional filters (e.g. size filters, if they have
        # been added)
        for my_filter in self._filters.values():
            reads = filter(my_filter, reads)

        # retrieve selected parts of regions
        count_dict = self.map_fn(reads, roi)

        # normalize to reads per million of normalization flag is set
        if self._normalize:
            count_dict = {k: v*1.0e6/self.sum() for (k, v) in count_dict.iteritems()}

        if roi_order and roi.strand == "-":
            count_dict = {k: v[::-1] for (k, v) in count_dict.iteritems()}

        return reads, count_dict

    def get_reads_and_counts(self,roi,roi_order=True):
        """Return :term:`read alignments` covering a |GenomicSegment|, and a
        count vector mapping reads to each positions in the |GenomicSegment|,
        following the rule specified by :meth:`~BAMGenomeArray.set_mapping`.
        Reads are strand-matched to `roi` by default. To obtain unstranded reads,
        set the value of `roi.strand` to `'.'`

        Parameters
        ----------
        roi : |GenomicSegment|
            Region of interest


        Returns
        -------
        list
            List of reads (as :class:`pysam.AlignedSegment`)
            covering region of interest

        numpy.ndarray
            Counts at each position of `roi`, under whatever mapping rule
            was set by :py:meth:`~BAMGenomeArray.set_mapping`


        Raises
        ------
        ValueError
            if bamfiles not sorted or not indexed
        """
        reads, count_dict = self.get_reads_and_hashed_counts(roi, roi_order)
        return reads, sum(count_dict.values(), numpy.zeros(len(roi)))


def get_hashed_counts(segchain, hashedgnd, stranded=True):
    """Returns a dict of counts of IVCollection as a list of positions, in
       transcript coordinates, keyed according to hashedgnd

    Parameters
    ----------
    segchain : plastid.roitools.SegmentChain
        Segment chain indicating positions from which to fetch counts

    hashedgnd : HashedReadBAMGenomeNDArray
        GenomeArray from which to fetch hashed counts

    stranded : bool, optional
        If `True` and the SegmentChain is on the minus strand,
        count order will be reversed relative to genome so that the
        array positions march from the 5' to 3' end of the chain.
        (Default: `True`)


    Returns
    -------
    numpy.ndarray
        Array of counts from `gnd` covering `self`
    """
    dcnt = {k: [] for k in hashedgnd.map_fn.read_keys}
    for iv in segchain:
        _, dtmp = hashedgnd.get_reads_and_hashed_counts(iv, roi_order=False)
        for (k, v) in dcnt.iteritems():
            v.extend(dtmp[k])
    if stranded and segchain.strand == "-":
        for v in dcnt.itervalues():
            v.reverse()
    return dcnt
