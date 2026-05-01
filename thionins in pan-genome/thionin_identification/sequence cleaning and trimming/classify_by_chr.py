#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: classify_by_chr.py

Description:
    Classify protein sequences into chromosome groups based on sequence IDs
    and write them into separate FASTA files.

    This script is designed for barley thionin analysis, where chromosome
    information (e.g. 1H–7H) is embedded in the sequence ID.

Input:
    - Directory containing FASTA files (e.g. intersection results)

Output:
    - One FASTA file per chromosome (e.g. 1H_thionin.fasta)
    - Optional summary table (CSV)

Usage:
    python classify_by_chr.py \
        --input_dir input_fasta_dir \
        --output_dir output_dir

Requirements:
    - Python >= 3.8
    - Biopython
"""

import argparse
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Classify FASTA sequences by chromosome based on sequence IDs."
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing input FASTA files. USER SET THIS PATH."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to save chromosome-classified FASTA files. USER SET THIS PATH."
    )

    parser.add_argument(
        "--output_prefix",
        default="thionin",
        help="Prefix for output FASTA files (default: thionin)"
    )

    parser.add_argument(
        "--summary",
        action="store_true",
        help="Generate summary CSV file"
    )

    return parser.parse_args()


def extract_chromosome(seq_id):
    """
    Extract chromosome ID from sequence ID.

    Expected format example:
        HORVU.xxx.xxx.6HGxxxx

    Returns:
        "1H"–"7H" or "Un"
    """
    try:
        field = seq_id.split(".")[3]

        if len(field) >= 2 and field[0].isdigit() and field[1] == "H":
            return field[:2]

    except IndexError:
        pass

    return "Un"


def main():
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    chrom_seqs = defaultdict(list)
    total_sequences = 0

    fasta_files = [
        f for f in input_dir.iterdir()
        if f.suffix.lower() in [".fa", ".fasta"]
    ]

    if not fasta_files:
        raise ValueError(f"No FASTA files found in {input_dir}")

    # =========================
    # Parse sequences
    # =========================
    for fasta_file in fasta_files:
        print(f"Processing: {fasta_file.name}")

        for record in SeqIO.parse(fasta_file, "fasta"):
            chrom = extract_chromosome(record.id)
            chrom_seqs[chrom].append(record)
            total_sequences += 1

    print(f"\nTotal sequences processed: {total_sequences}")

    # =========================
    # Write FASTA by chromosome
    # =========================
    summary_records = []

    for chrom, records in sorted(chrom_seqs.items()):
        out_file = output_dir / f"{chrom}_{args.output_prefix}.fasta"

        with open(out_file, "w") as f:
            SeqIO.write(records, f, "fasta")

        count = len(records)

        summary_records.append({
            "chromosome": chrom,
            "sequence_count": count,
            "output_file": str(out_file)
        })

        print(f"{chrom}: {count} sequences → {out_file}")

    # =========================
    # Save summary
    # =========================
    if args.summary:
        summary_df = pd.DataFrame(summary_records)
        summary_file = output_dir / "chromosome_summary.csv"
        summary_df.to_csv(summary_file, index=False)

        print(f"\nSummary written to: {summary_file}")


if __name__ == "__main__":
    main()
