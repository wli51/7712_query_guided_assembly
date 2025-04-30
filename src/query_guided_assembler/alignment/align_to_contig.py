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
                forward_matches = []
                reverse_matches = []

                for start in range(len(read) - min_align_len + 1):
                    subseq = read[start:start + min_align_len]
                    for qpos in _binary_search(contig, suffix_array, subseq):
                        forward_matches.append((read_id, contig_id, start + 1, start + min_align_len, qpos + 1, qpos + min_align_len))
                    subseq_rc = read_rc[start:start + min_align_len]
                    for qpos in _binary_search(contig, suffix_array, subseq_rc):
                        reverse_matches.append((read_id, contig_id, start + min_align_len, start + 1, qpos + 1, qpos + min_align_len))

                # Merge forward matches
                merged_forward_match = None
                for i, _match in enumerate(forward_matches):
                    if not merged_forward_match:
                        merged_forward_match = _match
                        continue
                    _prev = forward_matches[i - 1]
                    if (
                        _prev[2] + 1 == _match[2] and
                        _prev[3] + 1 == _match[3] and
                        _prev[4] + 1 == _match[4] and
                        _prev[5] + 1 == _match[5]
                    ):
                        merged_forward_match = (
                            merged_forward_match[0],  # sseqid
                            merged_forward_match[1],  # qseqid
                            merged_forward_match[2],  # sstart
                            _match[3],                 # new send
                            merged_forward_match[4],  # qstart
                            _match[5]                  # new qend
                        )
                    else:
                        matches.append(merged_forward_match)
                        merged_forward_match = _match
                if merged_forward_match:
                    matches.append(merged_forward_match)

                # Merge reverse matches
                merged_reverse_match = None
                for i, _match in enumerate(reverse_matches):
                    if not merged_reverse_match:
                        merged_reverse_match = _match
                        continue
                    _prev = reverse_matches[i - 1]
                    if (
                        _prev[2] + 1 == _match[2] and
                        _prev[3] + 1 == _match[3] and
                        _prev[4] + 1 == _match[4] and
                        _prev[5] + 1 == _match[5]
                    ):
                        merged_reverse_match = (
                            merged_reverse_match[0],  # sseqid
                            merged_reverse_match[1],  # qseqid
                            _match[2],                # new sstart
                            merged_reverse_match[3],  # send
                            merged_reverse_match[4],  # qstart
                            _match[5]                  # new qend
                        )
                    else:
                        matches.append(merged_reverse_match)
                        merged_reverse_match = _match
                if merged_reverse_match:
                    matches.append(merged_reverse_match)

            for match in matches:
                writer.writerow(match)
