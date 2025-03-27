from typing import List

def read_length_filter(
    reads: List[str], 
    min_length: int
) -> List[str]:
    """
    Filters reads based on a minimum length threshold.

    :param reads: List of DNA sequences (reads).
    :type reads: list[str]
    :param min_length: The minimum length a read must have to be retained.
    :type min_length: int
    :return: Filtered reads that meet the minimum length requirement.
    :rtype: list[str]
    """
    return [read for read in reads if len(read) >= min_length]
