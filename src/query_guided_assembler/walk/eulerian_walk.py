import copy
from typing import List, Optional
from collections import Counter
from random import Random

from ..graph import DeBruijnGraph
from .abstract_walk import AbstractWalk

class EulerianWalk(AbstractWalk):
    def __init__(self, graph: DeBruijnGraph, **kwargs):
        """
        Initialize the Eulerian walk class.
        Meant to be used with graph created by DeBruijnGraph.
        Performs an Eulerian walk through the graph, starting from a given node and 
            always choosing the edge with the least visited count.

        :param graph: adjacency list of the form Dict[str, List[str]]
        """
        super().__init__(graph, **kwargs)

    def walk(self, start_node: Optional[str] = None) -> List[str]:
        """
        Perform an Eulerian walk through the graph.
        """

        start_node = self._pre_walk_check(start_node)
        if not start_node:
            return []
        
        path = [start_node]
        current_node = start_node

        step= 0

        while self._get_neighbors(current_node):

            neighbors = self._get_neighbors(current_node)
            # choose the path with the least visited edge, when starting fresh, select the first neighbor
            next_node = min(neighbors, key=lambda n: self.visited_edges[(current_node, n)])

            self._logger.info(f"Step {step}: Moving from {current_node} -> {next_node}. "
                              f"(Visited {self.visited_edges[(current_node, next_node)]} times)"
                              )

            # update path and remove visited edge from graph
            self.visited_edges[(current_node, next_node)] += 1
            self._remove_edge(current_node, next_node)
            path.append(next_node)
            current_node = next_node

            step += 1

        self._logger.info(f"Eulerian walk completed. Path length: {len(path)}")

        return path

class StochasticEulerianWalk(AbstractWalk):
    def __init__(self, graph: DeBruijnGraph, **kwargs):
        """
        Initialize the stochastic eulerian walk class.
        Meant to be used with graph created by DeBruijnGraph.
        Performs a walk through the graph, starting from a given node and 
            chooses the edge stochastically weighted by out-degree.

        :param graph: adjacency list of the form Dict[str, List[str]]
        """
        super().__init__(graph, **kwargs)

    def walk(self, 
             start_node: Optional[str] = None, 
             random_seed: Optional[int] = None) -> List[str]:
        
        start_node = self._pre_walk_check(start_node)
        if not start_node:
            return []
        
        path = [start_node]
        current_node = start_node

        if random_seed is not None:            
            random = Random(random_seed)
        else:
            random = Random()

        step = 0

        while self._get_neighbors(current_node):

            neighbors = self._get_neighbors(current_node)
            weights = [self.edge_weights[(current_node, n)] for n in neighbors]

            # choose the path stochastically weighted by out-degree
            next_node = random.choices(neighbors, weights=weights, k=1)[0]

            self._logger.info(f"Step {step}: Moving from {current_node} -> {next_node}. "
                              f"(Visited {self.visited_edges[(current_node, next_node)]} times)"
                              )

            # update path and remove visited edge from graph
            self._remove_edge(current_node, next_node)
            path.append(next_node)
            current_node = next_node

            step += 1

        self._logger.info(f"Stochastic Eulerian walk completed. Path length: {len(path)}")

        return path

class StochasticGreedyEulerianWalk(AbstractWalk):
    def __init__(self, graph: DeBruijnGraph, **kwargs):
        """
        Initialize the stochastic eulerian walk class that uses a local greedy heuristic.
        Meant to be used with graph created by DeBruijnGraph.
        Performs a walk through the graph, starting from a given node and 
            chooses the edge stochastically weighted by out-degree and a greedy heuristic.

        :param graph: adjacency list of the form Dict[str, List[str]]
        """
        super().__init__(graph, **kwargs)

    def walk(self, 
             max_lookahead: int = 2, # large lookahead can cause very slow walk
             start_node: Optional[str] = None,
             random_seed: Optional[int] = None) -> List[str]:
        
        start_node = self._pre_walk_check(start_node)
        if not start_node:
            return []
        
        path = [start_node]
        current_node = start_node

        if random_seed is not None:            
            random = Random(random_seed)
        else:
            random = Random()

        step = 0

        while self._get_neighbors(current_node):

            neighbors = self._get_neighbors(current_node)

            # Score each neighbor by a local greedy heuristic
            neighbor_scores = {
                n: self._score_future_length(n, max_depth=max_lookahead)
                for n in neighbors
            }

            # Weight next step with both the local greedy score and out weights
            weights = [
                self.edge_weights[(current_node, n)] * (1 + neighbor_scores[n])
                for n in neighbors
            ]

            next_node = random.choices(neighbors, weights=weights, k=1)[0]

            self._logger.info(f"Step {step}: Moving from {current_node} -> {next_node}. "
                              f"(Visited {self.visited_edges[(current_node, next_node)]} times)"
                              )

            self._remove_edge(current_node, next_node)
            path.append(next_node)
            current_node = next_node

            step += 1

        self._logger.info(f"Stochastic Greedy Eulerian walk completed. Path length: {len(path)}")

        return path

    def _score_future_length(self, node: str, max_depth: int):
        """
        Recursively scores how far we can go from 'node' up to max_depth.
        """
        if max_depth == 0 or not self._get_neighbors(node):
            return 0

        return 1 + max(
            self._score_future_length(n, max_depth - 1)
            for n in self._get_neighbors(node)
        )