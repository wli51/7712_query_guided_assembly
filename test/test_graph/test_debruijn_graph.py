import pytest
from query_guided_assembler.graph import DeBruijnGraph

def test_k_must_be_integer():
    with pytest.raises(TypeError):
        dbg = DeBruijnGraph(k="not_an_int")

def test_sequences_must_be_list():
    dbg = DeBruijnGraph(k=3)
    with pytest.raises(TypeError):
        dbg.build_graph("not_a_list")
    with pytest.raises(TypeError):
        dbg.build_graph(123)
    with pytest.raises(TypeError):
        dbg.build_graph({})
    with pytest.raises(TypeError):
        dbg.build_graph([1, 2, 3])

def test_init_sets_k_and_empty_graph():
    dbg = DeBruijnGraph(k=3)
    assert dbg.k == 3
    assert dbg.graph == {}

def test_add_edge_creates_new_entry():
    dbg = DeBruijnGraph(k=3)
    assert dbg.graph == {}
    dbg._add_edge("AT", "TG")
    assert dbg.graph == {"AT": ["TG"]}

def test_add_edge_appends_to_existing():
    dbg = DeBruijnGraph(k=3)
    assert dbg.graph == {}
    dbg._add_edge("AT", "TG")
    dbg._add_edge("AT", "TC")
    assert "AT" in dbg.graph
    assert set(dbg.graph["AT"]) == set(["TG", "TC"])

def test_add_edge_repeatedly():
    dbg = DeBruijnGraph(k=3)
    assert dbg.graph == {}
    dbg._add_edge("AT", "TG")
    assert dbg.graph == {"AT": ["TG"]}
    dbg._add_edge("AT", "TG")
    assert dbg.graph == {"AT": ["TG", "TG"]}
    dbg._add_edge("AT", "TG")
    assert dbg.graph == {"AT": ["TG", "TG", "TG"]}

def test_build_graph_single_sequence():
    dbg = DeBruijnGraph(k=4)
    dbg.build_graph(["ATGCA"])
    # kmers of size 3 from "ATGCA" = ATG -> TGC -> GCA
    expected = {
        "ATG": ["TGC"],
        "TGC": ["GCA"]
    }
    assert dbg.graph == expected

def test_build_graph_multiple_sequences():
    dbg = DeBruijnGraph(k=4)
    # kmers of size 3 from "ATGCA" = ATG -> TGC -> GCA
    # kmers of size 3 from "TGCAA" = TGC -> GCA -> CAA
    dbg.build_graph(["ATGCA", "TGCAA"])
    expected = {
        "ATG": ["TGC"],
        "TGC": ["GCA", "GCA"],
        "GCA": ["CAA"]
    }
    assert dbg.graph == expected

def test_build_graph_variable_k():

    sequences = ["ATGCA", "TGCAA"]
    # kmers of size 4 from "ATGCA" = ATGC -> TGCA
    # kmers of size 4 from "TGCAA" = TGCA -> GCAA
    # kmers of size 3 from "ATGCA" = ATG -> TGC -> GCA
    # kmers of size 3 from "TGCAA" = TGC -> GCA -> CAA
    # kmers of size 2 from "ATGCA" = AT -> TG -> GC -> CA
    # kmers of size 2 from "TGCAA" = TG -> GC -> CA -> AA

    dbg = DeBruijnGraph(k=5)
    dbg.build_graph(sequences)
    expected = {
        "ATGC": ["TGCA"],
        "TGCA": ["GCAA"]
    }
    assert dbg.graph == expected

    dbg = DeBruijnGraph(k=4)
    dbg.build_graph(sequences)
    expected = {
        "ATG": ["TGC"],
        "TGC": ["GCA", "GCA"],
        "GCA": ["CAA"]
    }
    assert dbg.graph == expected

    dbg = DeBruijnGraph(k=3)
    dbg.build_graph(sequences)
    expected = {
        "AT": ["TG"],
        "TG": ["GC", "GC"],
        "GC": ["CA", "CA"],
        "CA": ["AA"]
    }
    assert dbg.graph == expected

def test_build_graph_empty_sequence_list():
    dbg = DeBruijnGraph(k=3)
    dbg.build_graph([])
    assert dbg.graph == {}

def test_build_graph_sequence_too_short():
    dbg = DeBruijnGraph(k=5)
    dbg.build_graph(["ATG"])  # length < k
    assert dbg.graph == {}