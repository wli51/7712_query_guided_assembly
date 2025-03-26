from typing import List

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
        self._k = k  # k-mer size
        self._graph = {}  # Adjacency list representation of the graph

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
        for sequence in sequences:
            for i in range(len(sequence) - self.k + 1):
                kmer1 = sequence[i:i + self.k - 1]
                kmer2 = sequence[i + 1:i + self.k]
                self._add_edge(kmer1, kmer2)

    @property
    def graph(self):
        """
        Return the current De-Bruijn graph as an adjacency list.
        """
        return self._graph
    
    @property
    def k(self):
        return self._k