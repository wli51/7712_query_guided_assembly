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

## Implementation Detals
TODO

## Installation
TODO

## Usage
TODO