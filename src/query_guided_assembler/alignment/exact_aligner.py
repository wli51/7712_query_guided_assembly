from typing import List, Dict
import bisect

def _build_suffix_array(
        query: str
        ) -> List[int]:
    """
    Helper function that constructs a suffix array for the given query sequence.

    :param query: Sequence to align reads to
    :type query: str
    :return: List of sorted suffix start positions
    :rtype: List[int]
    """

    # this may be memory intensive for large sequences
    suffixes = sorted((query[i:], i) for i in range(len(query)))
    suffix_array = [s[1] for s in suffixes]  # Store suffix starting positions

    return suffix_array

def _binary_search(
        query: str, 
        suffix_array: List[int], 
        read: str
        ) -> List[int]:
    """
    Helper function that uses binary search to efficiently find exact 
    matches of a read in the query sequence provided a pre-constructed suffix_array.

    :param query: Sequence to align reads to
    :type query: str
    :param suffix_array: Precomputed suffix array
    :type suffix_array: List[int]
    :param read: Sequence to be aligned to query
    :type read: str
    :return: List of match positions
    :rtype: List[int]
    """

    # find left-most occurence of read in query in sorted suffix array
    left = bisect.bisect_left(suffix_array, read, key=lambda i: query[i:i+len(read)])

    # find right-most occurence of read in query in sorted suffix array
    right = bisect.bisect_right(suffix_array, read, key=lambda i: query[i:i+len(read)])

    return suffix_array[left:right]  # Return all index position of exact matches

def exact_align_variable_reads(
        query: str, 
        reads: List[str]
        ) -> Dict[str, List]:
    """
    Aligns variable-length short reads exactly onto the long query using a suffix array.

    :param query: Sequence to align reads to
    :type query: str
    :param reads: A list of sequences to be aligned to query
    :type reads: List[str]
    :return: Dictionary mapping reads to a list of match positions with one key
    corresponding to each read. If no alignment is found, the list will be empty.
    :rtype: Dict[str, List]
    """
    suffix_array = _build_suffix_array(query)
    return {read: _binary_search(query, suffix_array, read) for read in reads}