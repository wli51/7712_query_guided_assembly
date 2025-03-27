import copy
from typing import Union, List
from abc import ABC, abstractmethod
from collections import Counter
import logging

from ..graph import DeBruijnGraph

class AbstractWalk(ABC):
    def __init__(self, graph: DeBruijnGraph, verbose: bool = False, **kwargs):
        """
        :param graph: adjacency list of the form Dict[str, List[str]]
        """

        if isinstance(graph, DeBruijnGraph):
            pass
        else:
            raise TypeError("Input graph must be an instance of DeBruijnGraph.")
        
        self.original_graph = graph.graph

        # self.graph = copy.deepcopy(graph.graph)
        self.edge_weights = graph.edge_weights()
        self.graph = None
        self.visited_edges = ModuleNotFoundError

        logging.basicConfig(level=logging.INFO, format="%(message)s")
        self._logger = logging.getLogger(f"{__name__}.id_{id(self)}")

        if verbose:
            self._logger.setLevel(logging.DEBUG)
        else:
            self._logger.setLevel(logging.INFO)

    def _remove_edge(self, u, v):
        """
        Remove a single edge from the graph and update edge_weights.
        """
        if v in self.graph.get(u, []):
            self.graph[u].remove(v)
            self.edge_weights[(u, v)] -= 1

    def _get_neighbors(self, node):
        return self.graph.get(node, [])
    
    def _pre_walk_check(self, start_node: Union[str, None] = None, **kwargs) -> Union[str, None]:
        """
        Helper function to perform pre-walk checks to be performed before starting the walk.

        :param start_node: Node to start the walk from.
        :type start_node: Union[str, None]
        :param kwargs: Additional keyword arguments.
        :return: Start node for the walk or None if graph is empty.
        """

        # reset visited edges
        self.visited_edges = Counter()
        # reset graph if needed
        self.graph = {k: v.copy() for k, v in self.original_graph.items()} # shallow copy

        if not self.graph or len(self.graph) == 0:
            self._logger.info("Graph is empty. No walk performed.")
            return None
        
        if start_node is not None:
            if start_node not in self.graph:
                raise ValueError("Start node not found in the graph.")
            self._logger.info(f"Starting Eulerian walk from provided node: {start_node}")
        else:
            start_node = max(self.graph, key=lambda node: len(self.graph[node]))
            self._logger.info(f"No start node provided. Starting from node with highest out-degree: {start_node}")

        return start_node
    
    @abstractmethod
    def walk(self, start_node: Union[str, None] = None):
        """
        Abstract method to be implemented by all traversal strategies.
        Should return a list of nodes (path).
        """
        pass