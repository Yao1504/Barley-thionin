#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: filter_blastp_hits.py

Description:
    Filter BLASTP results and extract candidate barley thionin homolog
    protein sequences.

    This script is designed to be used after running BLASTP searches
    against the barley proteome. It reads BLASTP tabular output files,
    filters hits by E-value, removes duplicated barley protein hits,
    and extracts the corresponding barley protein sequences.


Note:
This script is typically executed as part of the BLASTP pipeline via:
    run_thionin_pipeline.bat
However, it can also be run independently if BLASTP output
files are already available.

Required input:
    1. A folder containing BLASTP output files in outfmt 6 format:
       - arabidopsis_vs_morex.out
       - rice_vs_morex.out
       - others_vs_morex.out

    2. A barley protein FASTA file.
       The FASTA sequence IDs must match the BLAST subject IDs, namely sseqid.

BLAST outfmt 6 columns:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

Output:
    1. <prefix>_hits.csv
       Filtered non-redundant BLASTP hit table.

    2. <prefix>_homologs_filtered.fasta
       Barley protein sequences corresponding to filtered BLASTP hits.

    3. <prefix>_missing_ids.txt
       IDs found in BLASTP results but missing from the FASTA file.
       This file is generated only when missing IDs exist.

Usage:
    python filter_blastp_hits.py BLAST_DIR BARLEY_FASTA

Example:
    python filter_blastp_hits.py output Hv_Morex.pgsb.Jul2020.aa.fa
"""

import os
import sys
import pandas as pd
from Bio import SeqIO


# =========================
# User-adjustable settings
# =========================

EVALUE_CUTOFF = 1e-5          # USER CAN MODIFY: E-value threshold
OUTPUT_PREFIX = "Hv_thionin" # USER CAN MODIFY: output file prefix

BLAST_FILE_NAMES = [          # USER CAN MODIFY: expected BLASTP result files
    "arabidopsis_vs_morex.out",
    "rice_vs_morex.out",
    "others_vs_morex.out"
]


# =========================
# BLAST outfmt 6 columns
# =========================

COLNAMES = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore"
]


def check_arguments():
    if len(sys.argv) != 3:
        print("Usage:")
        print("  python filter_blastp_hits.py BLAST_DIR BARLEY_FASTA")
        print("")
        print("Example:")
        print("  python filter_blastp_hits.py output Hv_Morex.pgsb.Jul2020.aa.fa")
        sys.exit(1)

    blast_dir = sys.argv[1]
    barley_fasta = sys.argv[2]

    if not os.path.isdir(blast_dir):
        print(f"Error: BLAST directory not found: {blast_dir}")
        sys.exit(1)

    if not os.path.exists(barley_fasta):
        print(f"Error: barley FASTA file not found: {barley_fasta}")
        sys.exit(1)

    return blast_dir, barley_fasta


def load_blast_results(blast_dir):
    hits = []

    for file_name in BLAST_FILE_NAMES:
        file_path = os.path.join(blast_dir, file_name)

        if not os.path.exists(file_path):
            print(f"Warning: BLASTP result file not found and skipped: {file_path}")
            continue

        df = pd.read_csv(
            file_path,
            sep="\t",
            header=None,
            names=COLNAMES
        )

        if df.empty:
            print(f"Warning: BLASTP result file is empty and skipped: {file_path}")
            continue

        df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
        df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
        df["source_file"] = file_name

        hits.append(df)

    if not hits:
        print("Error: no valid BLASTP result files were found.")
        sys.exit(1)

    return pd.concat(hits, ignore_index=True)


def filter_hits(all_hits):
    filtered_hits = all_hits[all_hits["evalue"] < EVALUE_CUTOFF].copy()

    if filtered_hits.empty:
        print(f"Error: no BLASTP hits passed the E-value cutoff: {EVALUE_CUTOFF}")
        sys.exit(1)

    filtered_hits.sort_values(
        by=["sseqid", "evalue", "bitscore"],
        ascending=[True, True, False],
        inplace=True
    )

    unique_hits = filtered_hits.drop_duplicates(
        subset=["sseqid"],
        keep="first"
    )

    return unique_hits


def extract_sequences(barley_fasta, selected_ids, fasta_out):
    seq_records = SeqIO.to_dict(SeqIO.parse(barley_fasta, "fasta"))

    written_count = 0
    missing_ids = []

    with open(fasta_out, "w") as out_f:
        for sid in sorted(selected_ids):
            if sid in seq_records:
                SeqIO.write(seq_records[sid], out_f, "fasta")
                written_count += 1
            else:
                missing_ids.append(sid)

    return written_count, missing_ids


def main():
    blast_dir, barley_fasta = check_arguments()

    csv_out = os.path.join(blast_dir, f"{OUTPUT_PREFIX}_hits.csv")
    fasta_out = os.path.join(blast_dir, f"{OUTPUT_PREFIX}_homologs_filtered.fasta")
    missing_out = os.path.join(blast_dir, f"{OUTPUT_PREFIX}_missing_ids.txt")

    print("Loading BLASTP results...")
    all_hits = load_blast_results(blast_dir)

    print("Filtering BLASTP hits...")
    unique_hits = filter_hits(all_hits)

    print("Saving filtered hit table...")
    unique_hits.to_csv(csv_out, index=False)

    print("Extracting barley protein sequences...")
    selected_ids = set(unique_hits["sseqid"])
    written_count, missing_ids = extract_sequences(
        barley_fasta,
        selected_ids,
        fasta_out
    )

    if missing_ids:
        with open(missing_out, "w") as f:
            for sid in missing_ids:
                f.write(f"{sid}\n")

    print("")
    print("Summary")
    print("=======")
    print(f"Total BLASTP hits loaded: {len(all_hits)}")
    print(f"Total unique hits after filtering: {len(unique_hits)}")
    print(f"Protein sequences written: {written_count}")
    print(f"Hit details written to: {csv_out}")
    print(f"Filtered homolog sequences written to: {fasta_out}")

    if missing_ids:
        print(f"Warning: {len(missing_ids)} IDs were not found in FASTA.")
        print(f"Missing IDs written to: {missing_out}")


if __name__ == "__main__":
    main()
