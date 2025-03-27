from typing import List
from collections import Counter

def _is_low_complexity(
    read: str,
    low_complexity_threshold: float = 0.8
    ) -> bool:
    """
    Helper function that checks if a read is low-complexity (dominated by a single nucleotide).

    :param read: Single short read sequence.
    :type read: str
    :param low_complexity_threshold: Fraction of most common base allowed before discarding (default 0.8).
    When set to 0, the function will always return False (disabling low complexity threshold).
    :type low_complexity_threshold: float, optional
    :return: True if the read is low-complexity, False otherwise.
    """

    if len(read) == 0:
        # empty strings are treated as low complexity as well
        return True

    if low_complexity_threshold == 0:
        # when the threshold is 0, the low complexity filter is disabled
        return False

    counts = Counter(read)
    _, freq = counts.most_common(1)[0]
    return freq / len(read) >= low_complexity_threshold

def _passes_kmer_freq_filter(
    read: str,
    kmer_counts: dict[str, int],
    min_kmer_freq: int,
    k: int
    ) -> bool:
    """
    Helper function that checks if the read has sufficient high-frequency k-mers.

    :param read: Single short read sequence.
    :type read: str
    :param kmer_counts: Dictionary where keys are k-mers and values are their frequencies.
    :type kmer_counts: dict[str, int]
    :param min_kmer_freq: Minimum k-mer frequency threshold to retain a read.
    :type min_kmer_freq: int
    :param k: K-mer size.
    :type k: int
    :return: True if the read passes the k-mer frequency filter, False otherwise.
    """
    return all(kmer_counts[read[i:i + k]] >= min_kmer_freq for i in range(len(read) - k + 1))

def _compute_kmer_frequencies(
    reads: List[str], 
    k: int = 31
    ) -> dict[str, int]:
    """
    Computes k-mer frequencies from a list of reads and a specified k.

    :param reads: List of DNA sequences (reads).
    :type reads: list[str]
    :param k: K-mer size (default 31).
    :type k: int
    :return: Dictionary where keys are k-mers and values are their frequencies.
    :rtype: dict[str, int]
    """
    kmer_counts = Counter()
    for read in reads:
        if len(read) < k:
            continue
        for i in range(len(read) - k + 1):
            kmer = read[i:i + k]
            kmer_counts[kmer] += 1
    return kmer_counts

def read_kmer_freq_filter(
    reads: List[str], 
    k: int, 
    min_kmer_freq: int, 
    low_complexity_threshold: float = 0.8
    ) -> List[str]:
    """
    Filters reads based on k-mer frequency and complexity.

    :param reads: List of DNA sequences (reads).
    :type reads: list[str]
    :param k: K-mer size.
    :type k: int
    :param min_kmer_freq: Minimum k-mer frequency threshold to retain a read.
    :type min_kmer_freq: int
    :param low_complexity_threshold: Fraction of most common base allowed before discarding (default 0.8).
    If set to 0, the low complexity filter is disabled.
    :type low_complexity_threshold: float, optional
    :return: Filtered reads that pass both the k-mer frequency and complexity checks.
    :rtype: list[str]
    """
    kmer_counts = _compute_kmer_frequencies(reads, k)    

    return [
        read for read in reads if \
            _passes_kmer_freq_filter(read, kmer_counts, min_kmer_freq, k) \
            and (not _is_low_complexity(read, low_complexity_threshold))
        ]
