from typing import List

from Bio.Seq import Seq

def get_reverse_complement(reads: List[str]) -> List[str]:
    """
    Given a list of nucleotide sequences, return their reverse complements.

    :param reads: List of rna seq reads
    :return: List of reverse-complemented strings
    """
    return [str(Seq(seq).reverse_complement()) for seq in reads]

def get_reverse(reads: List[str]) -> List[str]:
    """
    Given a list of nucleotide sequences, return their reverse.

    :param reads: List of rna seq reads
    :return: List of reverse-complemented strings
    """
    return [seq[::-1] for seq in reads]
