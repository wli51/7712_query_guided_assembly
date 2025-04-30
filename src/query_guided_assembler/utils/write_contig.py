from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def write_contig(
    contig: str,
    output_fasta_path: str,
    contig_id: str = "contig1",
):
    record = SeqRecord(Seq(contig), id=contig_id, description="")
    with open(output_fasta_path, "w") as fasta_out:
        SeqIO.write(record, fasta_out, "fasta")