import pytest
from query_guided_assembler.alignment.exact_aligner import _build_suffix_array, _binary_search, exact_align_variable_reads

# Tests for _build_suffix_array

def test_build_suffix_array_simple():
    query = "banana"
    expected_suffix_array = [5, 3, 1, 0, 4, 2]  # sorted suffixes
    assert _build_suffix_array(query) == expected_suffix_array

def test_build_suffix_array_empty():
    assert _build_suffix_array("") == []

# Tests for _binary_search

def test_binary_search_single_match():
    query = "banana"
    suffix_array = _build_suffix_array(query)
    read = "ana"
    positions = _binary_search(query, suffix_array, read)
    assert sorted(positions) == [1, 3]

def test_binary_search_no_match():
    query = "banana"
    suffix_array = _build_suffix_array(query)
    read = "apple"
    assert _binary_search(query, suffix_array, read) == []

def test_binary_search_entire_match():
    query = "banana"
    suffix_array = _build_suffix_array(query)
    read = "banana"
    assert _binary_search(query, suffix_array, read) == [0]

def test_binary_search_partial_overlap():
    query = "aaaaa"
    suffix_array = _build_suffix_array(query)
    read = "aaa"
    positions = _binary_search(query, suffix_array, read)
    assert sorted(positions) == [0, 1, 2]

# Tests for exact_align_variable_reads

def test_exact_align_variable_reads_basic():
    query = "banana"
    reads = ["ana", "na", "ban"]
    result = exact_align_variable_reads(query, reads)
    assert set(result["ana"]) == set([1, 3])
    assert set(result["na"]) == set([2, 4])
    assert set(result["ban"]) == set([0])

def test_exact_align_variable_reads_no_match():
    query = "banana"
    reads = ["cat", "dog"]
    result = exact_align_variable_reads(query, reads)
    assert result == {"cat": [], "dog": []}

def test_exact_align_variable_reads_empty_query():
    result = exact_align_variable_reads("", ["a", "b"])
    assert result == {"a": [], "b": []}

def test_exact_align_variable_reads_empty_reads():
    result = exact_align_variable_reads("banana", [])
    assert result == {}

def test_exact_align_variable_reads_overlap():
    query = "aaaa"
    reads = ["aa"]
    result = exact_align_variable_reads(query, reads)
    assert set(result["aa"]) == set([0, 1, 2])