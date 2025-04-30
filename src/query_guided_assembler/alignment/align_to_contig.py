from typing import List
from pathlib import Path
import csv

from tqdm import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from ..alignment.exact_aligner import _build_suffix_array, _binary_search

def align_reads_to_contig(
    contig: str,
    reads: List[str],
    read_names: List[str],
    output_aln_path: Path,
    contig_id: str = "contig1",
    min_align_len: int = 20
):
    """
    Aligns reads and their reverse complements to a contig using exact matching,
    and outputs alignment and contig files.

    :param contig: Assembled contig sequence
    :param reads: List of sequencing reads
    :param read_names: Corresponding names of reads
    :param output_aln_path: Path to save tab-delimited alignment file
    :param contig_id: Identifier for the contig
    :param min_align_len: Minimum substring size for local alignment
    """
    assert len(reads) == len(read_names)

    # Prepare suffix array
    suffix_array = _build_suffix_array(contig)

    # Open alignment file
    with open(output_aln_path, "w", newline="") as aln_file:
        writer = csv.writer(aln_file, delimiter="\t")
        writer.writerow(["sseqid", "qseqid", "sstart", "send", "qstart", "qend"])

        for read, read_id in tqdm(zip(reads, read_names), desc="Aligning reads", total=len(reads)):
            matches = []

            # Try full forward match
            for qpos in _binary_search(contig, suffix_array, read):
                matches.append((read_id, contig_id, 1, len(read), qpos + 1, qpos + len(read)))

            # Try full reverse match
            read_rc = str(Seq(read).reverse_complement())
            for qpos in _binary_search(contig, suffix_array, read_rc):
                matches.append((read_id, contig_id, len(read), 1, qpos + 1, qpos + len(read)))

            # Local alignment windowed if no full match
            if not matches:
                for start in range(len(read) - min_align_len + 1):
                    subseq = read[start:start + min_align_len]
                    for qpos in _binary_search(contig, suffix_array, subseq):
                        matches.append((read_id, contig_id, start + 1, start + min_align_len, qpos + 1, qpos + min_align_len))
                    subseq_rc = read_rc[start:start + min_align_len]
                    for qpos in _binary_search(contig, suffix_array, subseq_rc):
                        matches.append((read_id, contig_id, start + min_align_len, start + 1, qpos + 1, qpos + min_align_len))

            for match in matches:
                writer.writerow(match)
