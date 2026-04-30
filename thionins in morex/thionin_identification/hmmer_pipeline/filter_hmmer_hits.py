#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: filter_hmmer_hits.py

Description:
    Filter HMMER hmmscan --tblout results and extract corresponding
    protein sequences from the target protein FASTA file.

    This script is typically called by:
        run_hmmer_pipeline.bat

    It can also be run independently if the HMMER tblout file is already
    available.

Required input:
    1. HMMER hmmscan table output file generated with --tblout
    2. Target protein FASTA file
    3. Output directory

Output:
    1. PF00321_hits.csv
       Filtered HMMER hit table

    2. PF00321_hits_filtered.fasta
       Protein sequences corresponding to filtered HMMER hits

    3. PF00321_missing_ids.txt
       IDs found in HMMER results but missing from the FASTA file.
       This file is generated only if missing IDs exist.

Usage:
    python filter_hmmer_hits.py <tbl_file> <protein_fasta> <output_dir>

Example:
    python filter_hmmer_hits.py output_hmm/PF00321_vs_morex.tbl input/barley_proteins.fasta output_hmm
"""

import sys
import os
import csv
from Bio import SeqIO


# =========================
# User-adjustable settings
# =========================

EVALUE_CUTOFF = 1e-5          # USER CAN MODIFY: E-value threshold
OUTPUT_PREFIX = "PF00321"    # USER CAN MODIFY: output file prefix


def check_arguments():
    if len(sys.argv) != 4:
        print("Usage:")
        print("  python filter_hmmer_hits.py <tbl_file> <protein_fasta> <output_dir>")
        print("")
        print("Example:")
        print("  python filter_hmmer_hits.py output_hmm/PF00321_vs_morex.tbl input/barley_proteins.fasta output_hmm")
        sys.exit(1)

    tbl_file = sys.argv[1]
    protein_file = sys.argv[2]
    output_dir = sys.argv[3]

    if not os.path.exists(tbl_file):
        print(f"Error: HMMER tblout file not found: {tbl_file}")
        sys.exit(1)

    if not os.path.exists(protein_file):
        print(f"Error: protein FASTA file not found: {protein_file}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    return tbl_file, protein_file, output_dir


def main():
    tbl_file, protein_file, output_dir = check_arguments()

    # =========================
    # Load protein sequences
    # =========================

    seq_dict = SeqIO.to_dict(SeqIO.parse(protein_file, "fasta"))

    # =========================
    # Define output files
    # =========================

    fasta_out = os.path.join(output_dir, f"{OUTPUT_PREFIX}_hits_filtered.fasta")
    csv_out = os.path.join(output_dir, f"{OUTPUT_PREFIX}_hits.csv")
    missing_out = os.path.join(output_dir, f"{OUTPUT_PREFIX}_missing_ids.txt")

    unique_hits = set()
    missing_ids = []

    # =========================
    # Parse and filter HMMER hits
    # =========================

    with open(tbl_file, "r") as f, \
         open(fasta_out, "w") as fasta_f, \
         open(csv_out, "w", newline="") as csv_f:

        writer = csv.writer(csv_f)

        writer.writerow([
            "HMM_name",
            "Target_protein",
            "E-value",
            "Score",
            "Bias",
            "Description"
        ])

        for line in f:
            line = line.strip()

            if line.startswith("#") or len(line) == 0:
                continue

            parts = line.split()

            if len(parts) < 7:
                continue

            hmm_name = parts[0]
            target_name = parts[2]
            evalue = float(parts[4])
            score = parts[5]
            bias = parts[6]
            description = " ".join(parts[18:]) if len(parts) > 18 else ""

            # Filter by E-value and keep one record per target protein
            if evalue <= EVALUE_CUTOFF and target_name not in unique_hits:
                unique_hits.add(target_name)

                writer.writerow([
                    hmm_name,
                    target_name,
                    evalue,
                    score,
                    bias,
                    description
                ])

                if target_name in seq_dict:
                    SeqIO.write(seq_dict[target_name], fasta_f, "fasta")
                else:
                    print(f"Warning: {target_name} not found in FASTA")
                    missing_ids.append(target_name)

    # =========================
    # Save missing IDs if any
    # =========================

    if missing_ids:
        with open(missing_out, "w") as f:
            for target_name in missing_ids:
                f.write(f"{target_name}\n")

    # =========================
    # Summary
    # =========================

    print("")
    print("Summary")
    print("=======")
    print(f"HMMER tblout file: {tbl_file}")
    print(f"Protein FASTA file: {protein_file}")
    print(f"E-value cutoff: {EVALUE_CUTOFF}")
    print(f"Unique HMMER hits after filtering: {len(unique_hits)}")
    print(f"Filtered hit table written to: {csv_out}")
    print(f"Filtered FASTA written to: {fasta_out}")

    if missing_ids:
        print(f"Warning: {len(missing_ids)} IDs were not found in FASTA.")
        print(f"Missing IDs written to: {missing_out}")


if __name__ == "__main__":
    main()
