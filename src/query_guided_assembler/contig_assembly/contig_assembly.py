from pathlib import Path

from ..utils.read_fasta import read_fasta_to_list
from ..preprocessing.read_length_filter import read_length_filter
from ..preprocessing.read_trimmer import read_trimmer
from ..preprocessing.get_reverse_complement import get_reverse_complement, get_reverse
from ..graph import DeBruijnGraph
from ..walk import (
    assemble_path,
    EulerianWalk,
    StochasticEulerianWalk,
    StochasticGreedyEulerianWalk
)

def run_query_guided_assembly(
    read_fasta_path: Path,
    query_fasta_path: Path,
    k: int = 31,
    min_length: int = 75,
    trim_size: int = 3,
    preprocess: bool = True,
    walk_type: str = "eulerian",
    n_runs: int = 10
) -> dict:
    """
    Run query-guided assembly using unstranded reads, extending both directions.

    :param read_fasta_path: Path to the input FASTA file containing reads.
    :param query_fasta_path: Path to the input FASTA file containing the query sequence.
    :param k: Length of k-mers used to construct the De Bruijn graph.
    :param min_length: Minimum length of reads to keep after filtering.
    :param trim_size: Number of bases to trim from both ends of each read.
    :param preprocess: Whether to preprocess the reads (length filter and trimming).
    :param walk_type: Type of walk to use for assembly. 
        Options are "eulerian", "stochastic", "stochastic_greedy".
    :param n_runs: Number of runs for stochastic walks. 
        Only the longest extension will be included in the final contig.
        Does not apply to eulerian walk because it is deterministic.
    :return: Dictionary containing the forward extension, reverse extension, and full contig.
    :rtype: dict
    :raises ValueError: If the walk_type is invalid or if the start k-mer is not found in the graph.
    :raises TypeError: If the input sequences are not in the expected format.
    """

    walk_cls_map = {
        "eulerian": EulerianWalk,
        "stochastic": StochasticEulerianWalk,
        "stochastic_greedy": StochasticGreedyEulerianWalk
    }
    walk_is_stochastic = {
        "eulerian": False,
        "stochastic": True,
        "stochastic_greedy": True
    }

    if walk_type not in walk_cls_map:
        raise ValueError(f"Invalid walk_type '{walk_type}'. Must be one of: {list(walk_cls_map.keys())}")

    walk_cls = walk_cls_map[walk_type]

    # Load reads and query
    reads = read_fasta_to_list(read_fasta_path)
    query = read_fasta_to_list(query_fasta_path)[0]

    # Optional preprocessing
    if preprocess:
        reads = read_length_filter(reads, min_length=min_length)
        reads = read_trimmer(reads, trim_size=trim_size)

    # Prepare reverse complements for forward graph
    reads_rc = get_reverse_complement(reads)

    # --- Forward Extension ---
    fwd_reads = reads + reads_rc
    dbg_fwd = DeBruijnGraph(k=k)
    dbg_fwd.build_graph(fwd_reads)

    start_kmer = query[-(k - 1):]
    if start_kmer not in dbg_fwd.graph:
        raise ValueError(f"Forward start k-mer '{start_kmer}' not found in forward graph.")

    walker_fwd = walk_cls(graph=dbg_fwd, verbose=False)
    if walk_is_stochastic[walk_type]:
        best_path_fwd = []
        for _ in range(n_runs):
            _path_fwd = walker_fwd.walk(start_node=start_kmer)
            if len(_path_fwd) > len(best_path_fwd):
                best_path_fwd = _path_fwd
        path_fwd = best_path_fwd
    else:
        path_fwd = walker_fwd.walk(start_node=start_kmer)

    if len(path_fwd) == 0:
        fwd_extension = ""
    else:
        fwd_extension = assemble_path(path_fwd)
        # trim away the first k-1 bases
        fwd_extension = fwd_extension[k - 1:]

    # --- Reverse Extension ---
    rev_reads = get_reverse(fwd_reads)  # reverse both reads and rc(reads)
    dbg_rev = DeBruijnGraph(k=k)
    dbg_rev.build_graph(rev_reads)

    reverse_start_kmer = query[:k - 1][::-1]
    if reverse_start_kmer not in dbg_rev.graph:
        raise ValueError(f"Reverse start k-mer '{reverse_start_kmer}' not found in reverse graph.")

    walker_rev = walk_cls(graph=dbg_rev, verbose=False)
    if walk_is_stochastic[walk_type]:
        best_path_rev = []
        for _ in range(n_runs):
            _path_rev = walker_rev.walk(start_node=reverse_start_kmer)
            if len(_path_rev) > len(best_path_rev):
                best_path_rev = _path_rev
        path_rev = best_path_rev
    else:
        path_rev = walker_rev.walk(start_node=reverse_start_kmer)

    if len(path_rev) == 0:
        rev_extension = ""
    else:
        rev_extension = assemble_path(path_rev)[::-1]  # reverse the result
        # trim away the last k-1 bases
        rev_extension = rev_extension[:-(k - 1)]

    # Combine into final contig
    full_contig = rev_extension + query + fwd_extension

    return {
        'forward_extension': path_fwd,
        'reverse_extension': path_rev,
        'full_contig': full_contig
    }
