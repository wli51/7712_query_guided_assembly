import pytest
from query_guided_assembler.preprocessing.read_kmer_freq_filter import read_kmer_freq_filter, _compute_kmer_frequencies, _passes_kmer_freq_filter

def test_read_kmer_freq_filter_basic():
    """
    Test filtering reads based on k-mer frequency threshold.
    """
    reads = ["ATCGATCG",
             "TCGATCGA",
             "CGATCGAT",
             "GATCGATC",
             "AAAA"]
    k = 4
    min_kmer_freq = 2
    low_complexity_threshold = 0
    
    expected = ["ATCGATCG", "TCGATCGA", "CGATCGAT", "GATCGATC"]

    assert read_kmer_freq_filter(reads, k, min_kmer_freq, low_complexity_threshold) == expected

def test_read_kmer_freq_filter_edge_case():
    """
    Test when a read containing kmers whose frequency is one below/exactly at the threshold
    """
    reads = ["ATCGATCG",
             "TCGATCGA",
             "CGATCGAT",
             "GATCGATC",
             "GGGGGG" # kmer frequency of GGGG is now 3, one below the threshold
             ] 
    k = 4
    min_kmer_freq = 4
    low_complexity_threshold = 0

    expected = reads[:-1] # last read should be filtered out
    assert read_kmer_freq_filter(reads, k, min_kmer_freq, low_complexity_threshold) == expected

    reads[-1] += "G" # adding another G to the last read making the kmer frequency 4
    expected = reads # now the last read should not be filtered out
    assert read_kmer_freq_filter(reads, k, min_kmer_freq, low_complexity_threshold) == expected


def test_read_kmer_freq_filter_low_complexity():
    """
    Test that low-complexity sequences are removed.
    """

    k = 4
    min_kmer_freq = 0
    low_complexity_threshold = 0.9

    reads = ["AAAAAAA", "CCCCCCCC", "GCTAGCTA"]
    expected = [reads[-1]] # only the last read should be retained
    assert read_kmer_freq_filter(reads, k, min_kmer_freq, low_complexity_threshold) == expected

def test_read_kmer_freq_filter_all_fail():
    """
    Test when all reads contain rare k-mers (relative to the threshold).
    """

    k = 4
    min_kmer_freq = 10 # artificially high threshold
    low_complexity_threshold = 0.0

    reads = ["ACGTACGT", "GTTGCCA", "AGAGAG"]
    assert read_kmer_freq_filter(reads, k, min_kmer_freq, low_complexity_threshold) == []

def test_read_kmer_freq_filter_empty_list():
    """
    Test handling of an empty input list.
    """

    k = 4
    min_kmer_freq = 2
    low_complexity_threshold = 0.8

    assert read_kmer_freq_filter([], k, min_kmer_freq, low_complexity_threshold) == []
