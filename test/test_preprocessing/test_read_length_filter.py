import pytest
from query_guided_assembler.preprocessing.read_length_filter import read_length_filter

def test_read_length_filter_basic():
    """
    Test filtering reads based on a minimum length threshold.
    """

    # Define a list of reads with the first two falling below the threshold
    # and the last one meeting the threshold
    reads = ["ATCG", "A", "GCTAGCTAGC"]
    assert read_length_filter(reads, min_length=5) == ["GCTAGCTAGC"]

def test_read_length_filter_empty_list():
    """
    Test handling of an empty input list.
    """
    assert read_length_filter([], min_length=5) == []

def test_read_length_filter_all_pass():
    """
    Test when all reads meet the length requirement.
    """
    reads = ["ATCGAT", "GCTAGCTAGC", "TTGACA"]
    assert read_length_filter(reads, min_length=3) == reads

def test_read_length_filter_all_fail():
    """
    Test when no reads meet the length requirement.
    """
    reads = ["A", "T", "G"]
    assert read_length_filter(reads, min_length=2) == []

def test_read_length_filter_edge_case():
    """
    Test when a read exactly meets the threshold.
    """
    reads = ["ATCGA", "TGCA", "AAAAA"]
    assert read_length_filter(reads, min_length=5) == ["ATCGA", "AAAAA"]
