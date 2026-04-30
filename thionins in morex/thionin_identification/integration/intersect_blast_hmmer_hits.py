#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: intersect_blastp_hmmer_hits.py

Description:
    Identify candidate thionin proteins detected by both BLASTP and HMMER.

    This script takes filtered BLASTP and HMMER results, identifies the
    intersection of protein IDs, keeps only the primary isoform ending with
    ".1", and extracts the corresponding protein sequences.

    This step corresponds to the manuscript statement:
    "Candidate thionin proteins identified by both BLASTP and HMMER were
    retained for further analysis. Redundant protein isoforms originating
    from alternative splicing of the same gene locus were removed."

Required input:
    1. HMMER filtered CSV file
       Example: PF00321_hits.csv

    2. HMMER filtered FASTA file
       Example: PF00321_hits_filtered.fasta

    3. BLASTP filtered CSV file
       Example: Hv_thionin_hits.csv

    4. BLASTP filtered FASTA file
       Example: Hv_thionin_homologs_filtered.fasta

Output:
    1. all_intersection.txt
       Protein IDs detected by both BLASTP and HMMER

    2. intersection_1.txt
       Intersected protein IDs ending with ".1"

    3. intersection_1.fasta
       FASTA sequences of intersected ".1" isoforms

Usage:
    python intersect_blastp_hmmer_hits.py \
        <hmm_csv> \
        <hmm_fasta> \
        <blast_csv> \
        <blast_fasta> \
        <output_dir>

Example:
    python intersect_blastp_hmmer_hits.py \
        output_hmm/PF00321_hits.csv \
        output_hmm/PF00321_hits_filtered.fasta \
        output_blast/Hv_thionin_hits.csv \
        output_blast/Hv_thionin_homologs_filtered.fasta \
        output_intersection
"""

import os
import sys
from Bio import SeqIO
import pandas as pd


# =========================
# User-adjustable settings
# =========================

HMM_ID_COLUMN = "Query_name"   # USER CAN MODIFY if HMMER CSV column name differs
BLAST_ID_COLUMN = "sseqid"     # USER CAN MODIFY if BLASTP CSV column name differs
PRIMARY_ISOFORM_SUFFIX = ".1"  # USER CAN MODIFY if another isoform rule is used


def check_arguments():
    if len(sys.argv) != 6:
        print("Usage:")
        print("  python intersect_blastp_hmmer_hits.py <hmm_csv> <hmm_fasta> <blast_csv> <blast_fasta> <output_dir>")
        print("")
        print("Example:")
        print("  python intersect_blastp_hmmer_hits.py output_hmm/PF00321_hits.csv output_hmm/PF00321_hits_filtered.fasta output_blast/Hv_thionin_hits.csv output_blast/Hv_thionin_homologs_filtered.fasta output_intersection")
        sys.exit(1)

    hmm_csv_path = sys.argv[1]
    hmm_fasta_path = sys.argv[2]
    blast_csv_path = sys.argv[3]
    blast_fasta_path = sys.argv[4]
    output_dir = sys.argv[5]

    input_files = [
        ("HMMER CSV file", hmm_csv_path),
        ("HMMER FASTA file", hmm_fasta_path),
        ("BLASTP CSV file", blast_csv_path),
        ("BLASTP FASTA file", blast_fasta_path),
    ]

    for file_label, file_path in input_files:
        if not os.path.exists(file_path):
            print(f"Error: {file_label} not found: {file_path}")
            sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    return hmm_csv_path, hmm_fasta_path, blast_csv_path, blast_fasta_path, output_dir


def load_hit_ids(csv_path, id_column, label):
    df = pd.read_csv(csv_path)

    if id_column not in df.columns:
        print(f"Error: column '{id_column}' not found in {label} CSV: {csv_path}")
        print(f"Available columns: {', '.join(df.columns)}")
        sys.exit(1)

    ids = set(
        df[id_column]
        .dropna()
        .astype(str)
        .str.strip()
    )

    return ids


def write_id_list(ids, output_path):
    with open(output_path, "w") as f:
        for protein_id in sorted(ids):
            f.write(protein_id + "\n")


def extract_fasta_by_ids(input_fasta, selected_ids, output_fasta):
    count = 0
    missing_ids = set(selected_ids)

    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            record_id = record.id.split()[0].strip()

            if record_id in selected_ids:
                SeqIO.write(record, out_f, "fasta")
                count += 1
                missing_ids.discard(record_id)

    return count, missing_ids


def main():
    hmm_csv_path, hmm_fasta_path, blast_csv_path, blast_fasta_path, output_dir = check_arguments()

    all_txt_path = os.path.join(output_dir, "all_intersection.txt")
    intersection_1_txt_path = os.path.join(output_dir, "intersection_1.txt")
    intersection_1_fasta_path = os.path.join(output_dir, "intersection_1.fasta")
    missing_ids_path = os.path.join(output_dir, "intersection_1_missing_ids.txt")

    # =========================
    # Load IDs from CSV files
    # =========================

    hmm_ids = load_hit_ids(
        csv_path=hmm_csv_path,
        id_column=HMM_ID_COLUMN,
        label="HMMER"
    )

    blast_ids = load_hit_ids(
        csv_path=blast_csv_path,
        id_column=BLAST_ID_COLUMN,
        label="BLASTP"
    )

    # =========================
    # Intersect BLASTP and HMMER IDs
    # =========================

    all_intersection_ids = hmm_ids.intersection(blast_ids)

    print(f"Total intersected protein IDs detected by both BLASTP and HMMER: {len(all_intersection_ids)}")

    write_id_list(all_intersection_ids, all_txt_path)

    # =========================
    # Keep only primary isoform ending with .1
    # =========================

    intersection_1_ids = {
        protein_id
        for protein_id in all_intersection_ids
        if protein_id.endswith(PRIMARY_ISOFORM_SUFFIX)
    }

    print(f"Intersected protein IDs ending with '{PRIMARY_ISOFORM_SUFFIX}': {len(intersection_1_ids)}")

    write_id_list(intersection_1_ids, intersection_1_txt_path)

    # =========================
    # Extract FASTA sequences
    # =========================

    count, missing_ids = extract_fasta_by_ids(
        input_fasta=hmm_fasta_path,
        selected_ids=intersection_1_ids,
        output_fasta=intersection_1_fasta_path
    )

    if missing_ids:
        write_id_list(missing_ids, missing_ids_path)

    # =========================
    # Summary
    # =========================

    print("")
    print("Summary")
    print("=======")
    print(f"HMMER CSV: {hmm_csv_path}")
    print(f"BLASTP CSV: {blast_csv_path}")
    print(f"All intersected IDs: {len(all_intersection_ids)}")
    print(f"Primary isoform IDs: {len(intersection_1_ids)}")
    print(f"FASTA sequences written: {count}")
    print(f"All intersection ID list: {all_txt_path}")
    print(f"Primary isoform ID list: {intersection_1_txt_path}")
    print(f"Primary isoform FASTA: {intersection_1_fasta_path}")

    if missing_ids:
        print(f"Warning: {len(missing_ids)} selected IDs were not found in FASTA.")
        print(f"Missing ID list: {missing_ids_path}")


if __name__ == "__main__":
    main()
