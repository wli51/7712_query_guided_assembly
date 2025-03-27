from typing import List

def read_trimmer(
    reads: List[str], 
    trim_size: int = 0
    ) -> List[str]:
    """
    Trims a specified number of bases from both ends of each read.

    :param reads: List of DNA sequences (reads).
    :type reads: list[str]
    :param trim_size: Number of bases to trim from both ends (default 0).
    When set to 0, the function will return the original reads.
    :type trim_size: int
    :return: Trimmed reads.
    :rtype: list[str]
    """

    # When trim size is 0, return the original reads
    if trim_size == 0:
        return reads

    return [
        read[trim_size:-trim_size] if len(read) > 2 * trim_size else \
        # if the read is shorter than 2 * trim_size, return the original untrimmed read
        read for read in reads
        ]
