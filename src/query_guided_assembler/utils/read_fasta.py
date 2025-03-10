from typing import List

from Bio import SeqIO

def read_fasta_to_list(
        fasta_path: str
        ) -> List[str]:
    """
    Reads a FASTA file and returns a list of sequences as strings.

    :param fasta_path: Path to the FASTA file
    :return: List of sequences (strings)
    """
    sequences = []
    with open(fasta_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(str(record.seq))
    return sequences
