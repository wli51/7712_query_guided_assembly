from typing import List
from collections import Counter

from tqdm import tqdm

class DeBruijnGraph:
    """
    De Bruijin graph builder class that stores edges as
    adjacency list
    """
    def __init__(self, k: int):
        """
        Initialize the De-Bruijn graph.
        :param k: Length of k-mers used to construct the graph.
        :type k: int
        """

        if not isinstance(k, int):
            raise TypeError("k must be an integer.")
        if k < 1:
            raise ValueError("k must be greater than 0.")

        self._k = k  # k-mer size
        self._graph = {}  # Adjacency list representation of the graph
        self._edge_weights = None # Edge weights for the graph

    def _add_edge(self, 
                 kmer1: str, 
                 kmer2: str):
        """
        Add a directed edge from kmer1 to kmer2 in the graph.
        :param kmer1: Starting k-mer.
        :type kmer1: str
        :param kmer2: Ending k-mer.
        :type kmer2: str
        """
        if kmer1 not in self.graph:
            self.graph[kmer1] = []
        self.graph[kmer1].append(kmer2)

    def build_graph(self, 
                    sequences: List[str]
                    ):
        """
        Build the De-Bruijn graph from a list of sequences.
        :param sequences: List of DNA sequences to extract k-mers from.
        :type sequences: List[str]
        """

        self._edge_weights = None # Reset edge weights

        if isinstance(sequences, list):
            if not all(isinstance(sequence, str) for sequence in sequences):
                raise TypeError("Input sequences must be a list of strings.")
            else:
                pass
        else:
            raise TypeError("Input sequences must be a list of strings.")

        for sequence in tqdm(sequences, desc="Building De-Bruijn Graph"):
            for i in range(len(sequence) - self.k + 1):
                kmer1 = sequence[i:i + self.k - 1]
                kmer2 = sequence[i + 1:i + self.k]
                self._add_edge(kmer1, kmer2)

    def _compute_edge_weights(self):
        """
        Compute edge weights for the graph.
        """
        self._edge_weights = Counter(
            (u, v)
            for u in tqdm(self.graph, desc="Computing Graph Outdegree")
            for v in self.graph[u]
        )

    @property
    def graph(self):
        """
        Return the current De-Bruijn graph as an adjacency list.
        """
        return self._graph
    
    def edge_weights(self):
        """
        Return the edge weights for the graph.
        """
        if self._edge_weights is None:
            self._compute_edge_weights()
        return self._edge_weights
    
    @property
    def k(self):
        return self._k