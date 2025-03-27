import pytest

from query_guided_assembler.graph import DeBruijnGraph
from query_guided_assembler.walk import EulerianWalk, StochasticEulerianWalk, StochasticGreedyEulerianWalk

@pytest.fixture
def empty_graph():
    return DeBruijnGraph(k=3)  # assumes constructor accepts dict adjacency

@pytest.fixture
def simple_graph():
    graph = DeBruijnGraph(k=4)
    graph.build_graph(["ATGCA"])
    return graph

@pytest.fixture
def start_node():
    return 'ATG'

@pytest.mark.parametrize("WalkClass", [
    EulerianWalk,
    StochasticEulerianWalk,
    StochasticGreedyEulerianWalk
])

def test_walk_on_empty_graph_returns_empty_list(WalkClass, empty_graph):
    
    walker = WalkClass(empty_graph)
    path = walker.walk()
    
    assert isinstance(path, list), f"{WalkClass.__name__}.walk() should return a list"
    assert len(path) == 0, f"{WalkClass.__name__}.walk() should return an empty path on empty graph"

@pytest.mark.parametrize("WalkClass", [
    EulerianWalk,
    StochasticEulerianWalk,
    StochasticGreedyEulerianWalk
])
def test_simple_walk(WalkClass, simple_graph, start_node):

    walker = WalkClass(simple_graph)
    path = walker.walk(start_node=start_node)

    assert len(path) == 3
    expected_path = ["ATG", "TGC", "GCA"]
    assert path[:3] == expected_path, f"Expected path to start with {expected_path}, but got {path[:3]}"