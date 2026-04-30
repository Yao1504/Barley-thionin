#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
filter_blast_hits.py

Description
-----------
Merge BLASTP tabular output files, filter significant hits by E-value,
keep the best hit for each target protein, and extract corresponding
protein sequences from a target FASTA file.

Required input files
--------------------
1. BLASTP output files in tabular format:
   BLAST outfmt 6 with the following columns:

   qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

2. Target protein FASTA file:
   The FASTA sequence IDs must match the BLAST subject IDs, namely sseqid.

Output files
------------
1. <prefix>_hits.csv
   Filtered and non-redundant BLAST hit table.

2. <prefix>_homologs_filtered.fasta
   Protein sequences corresponding to filtered target hits.

3. <prefix>_missing_ids.txt
   IDs present in BLAST hits but missing from the FASTA file.
   This file is generated only if missing IDs exist.

Example
-------
python filter_blast_hits.py \
    --blast_dir ./output \
    --target_fasta ./input/barley_proteins.fasta \
    --output_dir ./output \
    --prefix Hv_thionin \
    --evalue 1e-5 \
    --blast_files arabidopsis_vs_morex.out rice_vs_morex.out others_vs_morex.out
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Filter BLASTP hits and extract target homolog protein sequences."
        )
    )

    parser.add_argument(
        "--blast_dir",
        required=True,
        help="Directory containing BLAST output files. USER SET THIS PATH.",
    )

    parser.add_argument(
        "--target_fasta",
        required=True,
        help=(
            "Target protein FASTA file used as the BLAST database, "
            "for example barley proteome FASTA. USER SET THIS PATH."
        ),
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory for output files. USER SET THIS PATH.",
    )

    parser.add_argument(
        "--prefix",
        default="Hv_thionin",
        help=(
            "Prefix for output files. USER CAN MODIFY THIS VALUE. "
            "Default: Hv_thionin"
        ),
    )

    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help=(
            "E-value cutoff used to filter BLAST hits. "
            "USER CAN MODIFY THIS VALUE. Default: 1e-5"
        ),
    )

    parser.add_argument(
        "--blast_files",
        nargs="+",
        required=True,
        help=(
            "BLAST output file names inside --blast_dir. "
            "USER SET THESE FILE NAMES."
        ),
    )

    return parser.parse_args()


def check_file_exists(file_path: Path, description: str):
    if not file_path.exists():
        sys.exit(f"Error: {description} not found: {file_path}")


def load_blast_file(file_path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        file_path,
        sep="\t",
        header=None,
        names=BLAST_COLUMNS,
        comment="#",
    )

    if df.empty:
        return df

    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
    df["source_file"] = file_path.name

    return df


def load_all_blast_results(blast_dir: Path, blast_files: list[str]) -> pd.DataFrame:
    all_tables = []

    for filename in blast_files:
        file_path = blast_dir / filename

        if not file_path.exists():
            print(f"Warning: BLAST file not found and skipped: {file_path}")
            continue

        df = load_blast_file(file_path)

        if df.empty:
            print(f"Warning: BLAST file is empty and skipped: {file_path}")
            continue

        all_tables.append(df)

    if not all_tables:
        sys.exit("Error: No valid BLAST result files were loaded.")

    return pd.concat(all_tables, ignore_index=True)


def filter_hits(all_hits: pd.DataFrame, evalue_cutoff: float) -> pd.DataFrame:
    filtered = all_hits[all_hits["evalue"] < evalue_cutoff].copy()

    if filtered.empty:
        sys.exit(
            f"Error: No BLAST hits passed the E-value cutoff: {evalue_cutoff}"
        )

    filtered.sort_values(
        by=["sseqid", "evalue", "bitscore"],
        ascending=[True, True, False],
        inplace=True,
    )

    unique_hits = filtered.drop_duplicates(
        subset=["sseqid"],
        keep="first",
    )

    return unique_hits


def extract_fasta_sequences(
    target_fasta: Path,
    selected_ids: set[str],
    fasta_out: Path,
) -> tuple[int, list[str]]:
    seq_records = SeqIO.to_dict(SeqIO.parse(str(target_fasta), "fasta"))

    written_count = 0
    missing_ids = []

    with open(fasta_out, "w", encoding="utf-8") as out_f:
        for seq_id in sorted(selected_ids):
            if seq_id in seq_records:
                SeqIO.write(seq_records[seq_id], out_f, "fasta")
                written_count += 1
            else:
                missing_ids.append(seq_id)

    return written_count, missing_ids


def main():
    args = parse_args()

    blast_dir = Path(args.blast_dir)
    target_fasta = Path(args.target_fasta)
    output_dir = Path(args.output_dir)

    check_file_exists(blast_dir, "BLAST directory")
    check_file_exists(target_fasta, "Target FASTA file")

    output_dir.mkdir(parents=True, exist_ok=True)

    hits_csv = output_dir / f"{args.prefix}_hits.csv"
    fasta_out = output_dir / f"{args.prefix}_homologs_filtered.fasta"
    missing_ids_out = output_dir / f"{args.prefix}_missing_ids.txt"

    print("Loading BLAST result files...")
    all_hits = load_all_blast_results(
        blast_dir=blast_dir,
        blast_files=args.blast_files,
    )

    print("Filtering BLAST hits...")
    unique_hits = filter_hits(
        all_hits=all_hits,
        evalue_cutoff=args.evalue,
    )

    print("Writing filtered hit table...")
    unique_hits.to_csv(hits_csv, index=False)

    print("Extracting target protein sequences...")
    selected_ids = set(unique_hits["sseqid"].astype(str))

    written_count, missing_ids = extract_fasta_sequences(
        target_fasta=target_fasta,
        selected_ids=selected_ids,
        fasta_out=fasta_out,
    )

    if missing_ids:
        with open(missing_ids_out, "w", encoding="utf-8") as f:
            for seq_id in missing_ids:
                f.write(f"{seq_id}\n")

    print("\nSummary")
    print("-------")
    print(f"Total BLAST hits loaded: {len(all_hits)}")
    print(f"Unique target hits after filtering: {len(unique_hits)}")
    print(f"Protein sequences written: {written_count}")
    print(f"Filtered hit table: {hits_csv}")
    print(f"Filtered FASTA file: {fasta_out}")

    if missing_ids:
        print(f"Missing FASTA IDs: {len(missing_ids)}")
        print(f"Missing ID list: {missing_ids_out}")


if __name__ == "__main__":
    main()
