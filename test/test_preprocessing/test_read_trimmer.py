import pytest
from query_guided_assembler.preprocessing.read_trimmer import read_trimmer

def test_read_trimmer_basic():
    """
    Test trimming bases from both ends of the reads.
    """
    reads = ["ATCGATCG", "GCTAGCTA"]
    assert read_trimmer(reads, trim_size=2) == ["CGAT", "TAGC"]

def test_read_trimmer_trim_all():
    """
    Test trimming more bases than the read length.
    Should return the original untrimmed reads.
    """
    reads = ["ATCG", "GC"]
    assert read_trimmer(reads, trim_size=3) == reads

def test_read_trimmer_trim_none():
    """
    Test when trim_size is 0 (should return original reads).
    """
    reads = ["ATCGATCG", "GCTAGCTA"]
    assert read_trimmer(reads, trim_size=0) == reads

def test_read_trimmer_edge_case():
    """
    Test trimming when the reads lengths are exactly 2 * trim_size.
    Should not trim and return the original reads.
    """
    reads = ["ATCGATCG", "GCGCGCGC"]
    assert read_trimmer(reads, trim_size=4) == reads

def test_read_trimmer_mixed():
    """
    Test trimming when the first read is trimmed but the second is not.
    """
    reads = ["ATCGATCGATCG", "GCGC"]
    assert read_trimmer(reads, trim_size=4) == ["ATCG", "GCGC"]

def test_read_trimmer_default():
    """
    Test trimming with the default trim_size (which is 0)
    Should return the original reads.
    """
    reads = ["TCGA", "GG"]
    assert read_trimmer(reads) == reads

def test_read_trimmer_empty_list():
    """
    Test handling of an empty input list.
    """
    assert read_trimmer([], trim_size=2) == []
