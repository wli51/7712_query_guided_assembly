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
├── README.md
├── LICENSE
├── environment.yml
├── pytest.ini
├── setup.py
├── src/
│   └── query_guided_assembler/
│       ├── alignment/
│       │   └── exact_aligner.py
│       ├── graph/
│       │   └── debruijin.py
│       ├── preprocessing/
│       │   ├── read_kmer_freq_filter.py
│       │   ├── read_length_filter.py
│       │   └── read_trimmer.py
│       ├── utils/
│       │   └── read_fasta.py
│       └── walk/
│           ├── abstract_walk.py
│           ├── eulerian_walk.py
│           └── utils.py
└── test/
    ├── test_aligner/
    │   └── test_alignment.py
    ├── test_graph/
    │   └── test_debruijn_graph.py
    └── test_preprocessing/
        ├── test_read_kmer_freq_filter.py
        ├── test_read_length_filter.py
        └── test_read_trimmer.py
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

#### `EulerianWalk`

- Approximates an Eulerian path, aiming to visit all edges in the graph.
- Uses a least-visited edge heuristic to avoid overusing paths and reduce premature termination.
- This approach is deterministic.

#### `StochasticEulerianWalk`

- Approximates an Eulerian path, selecting the next edge **randomly**, but with stochasticity where traversal is weighted by edge outdegree.
- Encourages exploration of alternative paths by avoiding deterministic behavior, while still favoring more frequently observed transitions.

#### `StochasticGreedyEulerianWalk`

- Introduces stochasticity while favoring longer walks:
  - At each node, the next neighbor is selected randomly, but weighted by edge count (probabilistic component).
  - Incorporates a greedy lookahead heuristic that recursively checks potential future routes, to enourage neighbors that are expected to lead to longer paths. This significantly adds to the run complexity.

All walk classes output a list of nodes representing a path, which is then converted into a full assembled sequence (a **contig**) using the `assemble_path()` function. This utility reconstructs the sequence by merging overlaps between consecutive (k-1)-mers in the path.

### Alignment Module

After assembling the longest contig containing the query, the tool performs **exact alignments** of the original sequencing reads against the contig using a **suffix array–based alignment algorithm**. This alignment process:

- First attempts to align each read and its reverse complement as a **full-length exact match** to the contig.
- If no full match is found, it performs a **sliding window search** of each read and its reverse complement to identify **local exact matches** using a configurable minimum match length (e.g., 20 bp).
- All matching regions are recorded in a tab-delimited alignment file (`ALLELES.aln`), with strand-aware coordinate formatting: matches on the reverse strand are reported with `send < sstart`.

The alignment implementation is efficient and deterministic, relying entirely on binary search over a suffix array built from the contig. This allows rapid identification of both full-length and partial matches across many reads. The alignment suite is located in `alignment/exact_aligner.py` and `alignment/aligner.py`.

The wrapping contig-read alignment function first attempts full read alignment with the contig and falls back to attempting to align smaller substrings of the reads against the contig. This is feasible due to the alignment itself being relatively efficient.

## Walk Example Usage (Programmatic)
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

To install the package and expose the CLI tool `query-guided-assembly`, run:

This will:
- Make the `query_guided_assembler` package importable in your Python environment
- Enable usage of the CLI command `query-guided-assembly` from anywhere in your terminal

```bash
pip install -e .
```
---

## Command-Line Usage

Once installed, the tool can be run via the `query-guided-assembly` command.

### Example:

```bash
query-guided-assembly \
  --reads data/READS.fasta \
  --query data/QUERY.fasta \
  --k 31 \
  --min_length 75 \
  --trim 3 \
  --walk_type stochastic_greedy \
  --n_runs 10 \
  --min_align_len 20
```

### Required Arguments:
- `--reads`: path to input FASTA file of sequencing reads (can be gzipped)
- `--query`: path to input FASTA file of query sequence

### Optional Arguments:
- `--k`: k-mer size (default: 31)
- `--min_length`: minimum read length to keep (default: 75)
- `--trim`: number of bases to trim from both ends of each read (default: 3)
- `--walk_type`: assembly walk type (`eulerian`, `stochastic`, or `stochastic_greedy`)
- `--n_runs`: number of stochastic walk repetitions (default: 10)
- `--min_align_len`: minimum substring length used for local alignment (default: 20)