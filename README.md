# 7712_query_guided_assembly

## Overview

This project is a course assignment to implement a short read sequence assembly tool that reconstructs the largest contig containing a given query sequence. Intended uses of such query guided assembly tools include detecting contamination and identifying sequence contexts around a given query sequence in biological samples.

## Input Files
The program expects the following fasta files as inputs:
1. `QUERY.fasta`
   - Contains the **query sequence** to be searched for within the reads.
2. `READS.fasta.gz`
   - Contains **sequencing reads**. May be gzipped.

## Output Files
The program produces:
1. `ALLELES.fasta`
   - A FASTA file containing the **longest assembled contig** that contains the query sequence.
2. `ALLELES.aln`
   - A **tab-delimited text file** describing the alignment of sequencing reads to the reconstructed contig(s). See table below for format:

| sseqid | qseqid | sstart | send | qstart | qend |
|--------|--------|--------|------|--------|------|
| Read_001 | Contig_1 | 5 | 25 | 1 | 21 |
| Read_002 | Contig_1 | 50 | 70 | 46 | 66 |
| Read_003 | Contig_1 | 120 | 100 | 110 | 90 |

- **sseqid:** Sequencing read name (from `READS.fasta.gz`).
- **qseqid:** Contig name (from `ALLELES.fasta`).
- **sstart, send:** Start and end coordinates of the sequencing read's alignment to the contig.
- **qstart, qend:** Start and end coordinates of the contig’s alignment to the sequencing read.
- Reverse strand alignments are indicated by `send < sstart`.

## Repository Structure

```
.
├── README.md                    # This file
├── LICENSE
├── environment.yml              # Conda environment specification
├── pytest.ini                   # Pytest configuration
├── src/
│   └── query_guided_assembler/
│       ├── alignment/
│       │   └── exact_aligner.py         # Suffix array-based exact aligner
│       ├── graph/
│       │   └── debruijin.py             # De Bruijn graph construction logic
│       ├── preprocessing/
│       │   ├── read_kmer_freq_filter.py # K-mer-based and low-complexity filtering
│       │   ├── read_length_filter.py    # Read length filtering
│       │   └── read_trimmer.py          # End trimming for reads
│       ├── utils/
│       │   └── read_fasta.py            # FASTA/FASTQ I/O utilities
│       └── walk/
│           ├── abstract_walk.py         # Base class for traversal strategies
│           ├── eulerian_walk.py         # Eulerian walk strategies
│           └── utils.py                 # Path assembly and walk scoring helpers
├── test/
│   ├── test_aligner/
│   │   └── test_alignment.py            # Unit tests for the exact aligner
│   ├── test_graph/
│   │   └── test_debruijn_graph.py       # Tests for De Bruijn graph construction
│   └── test_preprocessing/
│       ├── test_read_kmer_freq_filter.py
│       ├── test_read_length_filter.py
│       └── test_read_trimmer.py
```

The `src/` directory contains all functional modules, while the `test/` directory contains unit tests organized by module. This structure follows standard Python packaging conventions and is compatible with `pytest`.

## Implementation Detals

### Preprocessing Module

The preprocessing module includes functions for filtering and trimming reads before assembly. These functions help reduce noise and improve assembly accuracy by removing low-quality or uninformative reads.

#### Functions

- `read_trimmer(reads: List[str], trim_size: int) -> List[str]`  
  Trims a specified number of bases from both ends of each read. If `trim_size` is 0, returns the original reads. Meant to remove ends of the read sequence where the read quality is potentially low.

- `read_length_filter(reads: List[str], min_length: int) -> List[str]`  
  Filters out reads shorter than a specified minimum length.

- `read_kmer_freq_filter(reads: List[str], k: int, min_kmer_freq: int, low_complexity_threshold: float = 0.8) -> List[str]`  
  Filters reads based on:
  - **K-mer frequency**: retains reads only if all their k-mers appear at least `min_kmer_freq` times across the dataset. Meant to exclude reads containing low frequency kmers relative to the pool due to them being potentially erroneous reads. 
  - **Low-complexity check**: optionally removes reads dominated by a single nucleotide (controlled by `low_complexity_threshold`) to reduce the chance of downstream assembly to be stuck in short repeat cycles.

#### Helper Functions

- `_compute_kmer_frequencies(reads: List[str], k: int) -> dict[str, int]`  
  Computes the frequency of each k-mer in the list of reads.

- `_passes_kmer_freq_filter(read: str, kmer_counts: dict[str, int], min_kmer_freq: int, k: int) -> bool`  
  Returns True if all k-mers in a read meet the frequency threshold.

- `_is_low_complexity(read: str, low_complexity_threshold: float = 0.8) -> bool`  
  Returns True if a read is composed primarily of a single nucleotide. This helps remove artifacts and repetitive noise.

These preprocessing steps can be chained together to clean raw sequencing reads before constructing the De Bruijn graph.

### The core of this project consists of two components: **De Bruijn graph construction** and **graph traversal strategies** that enable query-guided sequence assembly.

### De Bruijn Graph Construction

Sequencing reads are decomposed into overlapping **k-mers**, and a **De Bruijn graph** is constructed using these k-mers:

- **Nodes:** (k-1)-mers
- **Edges:** Directed edges between (k-1)-mers that overlap by k-2 bases, derived from each k-mer  
  Example: the k-mer `ATG` creates an edge from `AT` to `TG`.
- The graph is implemented as an adjacency list, and multiple occurrences of the same edge are tracked via **edge weights**.

The graph construction is handled by the `DeBruijnGraph` class, which supports building graphs from lists of sequences as well as adding individual edges if needed.

### Graph Traversal Strategies

Once the graph is built, a walk algorithm is applied to reconstruct the longest contig that contains the query sequence.

The project currently supports three traversal strategies, implemented as classes that inherit from a common `AbstractWalk` abstract class:

#### Key Features

- **Graph Validation and Copying**: Accepts a `DeBruijnGraph` instance and internally deep-copies its adjacency list and edge weights to ensure isolation from the original graph object.
- **Edge Management**: Provides helper methods to retrieve neighbors and remove edges during traversal.
- **Edge Weight Tracking**: Uses `Counter` objects to track how often edges are visited — useful for algorithms like the Eulerian walk.
- **Logging Support**: Includes verbose/debug logging for inspecting walk behavior.
- **Walk Entry Logic**: Includes a `_pre_walk_check()` helper to determine a valid starting node, either user-provided or chosen based on out-degree.

#### Eulerian Walk

- Approximates an Eulerian path, aiming to visit all edges in the graph.
- Uses a least-visited edge heuristic to avoid overusing paths and reduce premature termination.
- This approach is deterministic.

#### Greedy StochasticEulerianWalk Walk

- Approximates an Eulerian path, selecting the next edge **randomly**, but with stochasticity where traversal is weighted by edge outdegree.
- Encourages exploration of alternative paths by avoiding deterministic behavior, while still favoring more frequently observed transitions.

#### StochasticGreedyEulerianWalk

- Introduces stochasticity while favoring longer walks:
  - At each node, the next neighbor is selected randomly, but weighted by edge count (probabilistic component).
  - Incorporates a greedy lookahead heuristic that recursively checks potential future routes, to enourage neighbors that are expected to lead to longer paths. This significantly adds to the run complexity.

Both walk classes output a list of nodes representing a path, which is then converted into a full assembled sequence (a **contig**) using the `assemble_path()` function. This utility reconstructs the sequence by merging overlaps between consecutive (k-1)-mers in the path.

## Example Usage
```python
from query_guided_assembler.utils.read_fasta import read_fasta_to_list
from query_guided_assembler.graph import DeBruijnGraph
from query_guided_assembler.walk import EulerianWalk, StochasticEulerianWalk, StochasticGreedyEulerianWalk, assemble_path

# Read inputs
data_dir = pathlib.Path('./data')
reads = read_fasta_to_list(data_dir / 'READS.fasta')
query = read_fasta_to_list(data_dir / 'QUERY.fasta')[0]

# Define k
k = 31
start_kmer = query[-(k-1):]

# Build Graph
dbg = DeBruijnGraph(k=k)
dbg.build_graph(reads)

# Different walk strategies
walker = EulerianWalk(graph=dbg)
path = walker.walk(start_node=start_kmer)

stochastic_walker = StochasticEulerianWalk(graph=dbg)
path2 = stochastic_walker.walk(start_node=start_kmer)

stochastic_greedy_walker = StochasticGreedyEulerianWalk(graph=dbg)
path3 = stochastic_greedy_walker.walk(start_node=start_kmer, max_lookahead=3)
```

## Installation
TODO