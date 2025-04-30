import argparse
from pathlib import Path

from ..contig_assembly import run_query_guided_assembly
from ..alignment.align_to_contig import align_reads_to_contig
from ..utils.read_fasta import read_fasta_to_dict
from ..utils.write_contig import write_contig

def main():
    parser = argparse.ArgumentParser(
        description="Run query-guided assembly and align reads to the resulting contig."
    )
    parser.add_argument("--reads", type=str, required=True, help="Path to the input FASTA file of reads.")
    parser.add_argument("--query", type=str, required=True, help="Path to the input FASTA file of the query.")
    parser.add_argument("--k", type=int, default=31, help="k-mer size.")
    parser.add_argument("--min_length", type=int, default=75, help="Minimum read length to retain.")
    parser.add_argument("--trim", type=int, default=3, help="Trim this many bases from both ends of reads.")
    parser.add_argument("--walk_type", choices=["eulerian", "stochastic", "stochastic_greedy"], default="eulerian",
                        help="Type of walk for assembly.")
    parser.add_argument("--n_runs", type=int, default=10, help="Number of runs for stochastic walks.")
    parser.add_argument("--min_align_len", type=int, default=20, help="Minimum alignment length for local matches.")

    args = parser.parse_args()

    # Set output paths in current working directory
    cwd = Path.cwd()
    fasta_out = cwd / "ALLELES.fasta"
    aln_out = cwd / "ALLELES.aln"
    log_out = cwd / "assembly_alignment.log"

    print(f"[INFO] Starting assembly with reads from {args.reads} and query from {args.query}")
    print(f"[INFO] Outputs will be saved to: {fasta_out}, {aln_out}, {log_out}")

    # Run assembly
    assembly_output = run_query_guided_assembly(
        read_fasta_path=Path(args.reads),
        query_fasta_path=Path(args.query),
        k=args.k,
        min_length=args.min_length,
        trim_size=args.trim,
        preprocess=True,
        walk_type=args.walk_type,
        n_runs=args.n_runs,
        log_path=log_out
    )

    contig = assembly_output["full_contig"]
    write_contig(contig, fasta_out)

    # Align reads to the contig
    read_dict = read_fasta_to_dict(Path(args.reads))
    reads = list(read_dict.values())
    read_names = list(read_dict.keys())

    align_reads_to_contig(
        contig=contig,
        reads=reads,
        read_names=read_names,
        output_aln_path=aln_out,
        contig_id="contig1",
        min_align_len=args.min_align_len
    )

    print(f"[INFO] Finished. Results written to {fasta_out}, {aln_out}, and log saved at {log_out}")

if __name__ == "__main__":
    main()
