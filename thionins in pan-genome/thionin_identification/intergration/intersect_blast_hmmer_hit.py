#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: intersect_blast_hmmer_hit.py

Description:
    Identify thionin candidate proteins detected by both BLASTP and HMMER
    across multiple pangenome protein datasets.

    For each genome/accession, this script:
        1. Reads BLASTP hit IDs
        2. Reads HMMER hit IDs
        3. Calculates the intersection
        4. Writes intersected IDs to TXT
        5. Extracts corresponding protein sequences from the proteome FASTA

Input:
    1. Directory containing BLASTP output files
    2. Directory containing HMMER output files
    3. Directory containing pangenome protein FASTA files

Output:
    1. <sample>_intersection.txt
    2. <sample>_intersection.fa
    3. pangenome_intersection_summary.csv

Usage:
    python run_thionin_pangenome_blastp_hmmer_intersection.py \
        --blastp_dir output_blastp \
        --hmmsearch_dir output_hmmsearch \
        --proteome_dir pangenome_protein \
        --output_dir output_intersection
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Intersect BLASTP and HMMER hits across pangenome protein datasets."
    )

    parser.add_argument(
        "--blastp_dir",
        required=True,
        help="Directory containing BLASTP output files. USER SET THIS PATH."
    )

    parser.add_argument(
        "--hmmsearch_dir",
        required=True,
        help="Directory containing HMMER hmmsearch output files. USER SET THIS PATH."
    )

    parser.add_argument(
        "--proteome_dir",
        required=True,
        help="Directory containing pangenome protein FASTA files. USER SET THIS PATH."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory. USER SET THIS PATH."
    )

    parser.add_argument(
        "--blast_suffix",
        default="_blastp.txt",
        help="Suffix of BLASTP result files. Default: _blastp.txt"
    )

    parser.add_argument(
        "--hmm_suffix",
        default=".protein_hmmsearch.txt",
        help="Suffix of HMMER result files. Default: .protein_hmmsearch.txt"
    )

    parser.add_argument(
        "--proteome_suffix",
        default=".protein.fa",
        help="Suffix of proteome FASTA files. Default: .protein.fa"
    )

    return parser.parse_args()


def check_directory(path, label):
    if not Path(path).exists():
        sys.exit(f"Error: {label} not found: {path}")


def get_sample_name(filename, blast_suffix):
    if filename.endswith(blast_suffix):
        return filename.replace(blast_suffix, "")
    return Path(filename).stem


def read_blastp_ids(blast_file):
    ids = set()

    with open(blast_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            cols = line.strip().split("\t")

            if len(cols) < 2:
                continue

            subject_id = cols[1].strip()
            ids.add(subject_id)

    return ids


def read_hmmer_ids(hmm_file):
    ids = set()

    with open(hmm_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            cols = line.strip().split()

            if len(cols) < 1:
                continue

            target_name = cols[0].strip()
            ids.add(target_name)

    return ids


def write_id_list(ids, output_txt):
    with open(output_txt, "w") as out_f:
        for seq_id in sorted(ids):
            out_f.write(seq_id + "\n")


def extract_fasta_by_ids(proteome_fasta, selected_ids, output_fasta):
    written_count = 0
    missing_ids = set(selected_ids)

    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(proteome_fasta, "fasta"):
            record_id = record.id.strip()

            if record_id in selected_ids:
                SeqIO.write(record, out_f, "fasta")
                written_count += 1
                missing_ids.discard(record_id)

    return written_count, missing_ids


def main():
    args = parse_args()

    blastp_dir = Path(args.blastp_dir)
    hmmsearch_dir = Path(args.hmmsearch_dir)
    proteome_dir = Path(args.proteome_dir)
    output_dir = Path(args.output_dir)

    check_directory(blastp_dir, "BLASTP directory")
    check_directory(hmmsearch_dir, "HMMER directory")
    check_directory(proteome_dir, "Proteome directory")

    output_dir.mkdir(parents=True, exist_ok=True)

    summary_records = []

    blast_files = sorted([
        f for f in os.listdir(blastp_dir)
        if f.endswith(args.blast_suffix)
    ])

    if not blast_files:
        sys.exit(f"Error: no BLASTP files ending with '{args.blast_suffix}' found in {blastp_dir}")

    for blast_file in blast_files:
        sample_name = get_sample_name(blast_file, args.blast_suffix)

        blast_path = blastp_dir / blast_file
        hmm_path = hmmsearch_dir / f"{sample_name}{args.hmm_suffix}"
        proteome_path = proteome_dir / f"{sample_name}{args.proteome_suffix}"

        if not hmm_path.exists():
            print(f"Warning: HMMER file missing for {sample_name}, skipping.")
            continue

        if not proteome_path.exists():
            print(f"Warning: proteome FASTA file missing for {sample_name}, skipping.")
            continue

        blast_ids = read_blastp_ids(blast_path)
        hmm_ids = read_hmmer_ids(hmm_path)

        intersection_ids = blast_ids.intersection(hmm_ids)

        intersection_txt = output_dir / f"{sample_name}_intersection.txt"
        intersection_fasta = output_dir / f"{sample_name}_intersection.fa"
        missing_ids_txt = output_dir / f"{sample_name}_missing_ids.txt"

        write_id_list(intersection_ids, intersection_txt)

        written_count, missing_ids = extract_fasta_by_ids(
            proteome_fasta=proteome_path,
            selected_ids=intersection_ids,
            output_fasta=intersection_fasta
        )

        if missing_ids:
            write_id_list(missing_ids, missing_ids_txt)

        summary_records.append({
            "sample": sample_name,
            "blastp_hit_count": len(blast_ids),
            "hmmer_hit_count": len(hmm_ids),
            "intersection_count": len(intersection_ids),
            "fasta_sequence_count": written_count,
            "missing_id_count": len(missing_ids),
            "intersection_txt": str(intersection_txt),
            "intersection_fasta": str(intersection_fasta)
        })

        print(f"{sample_name}: {len(intersection_ids)} intersected sequences; {written_count} written to FASTA")

    summary_df = pd.DataFrame(summary_records)
    summary_out = output_dir / "pangenome_intersection_summary.csv"
    summary_df.to_csv(summary_out, index=False)

    print("")
    print("Pipeline finished.")
    print(f"Summary table written to: {summary_out}")


if __name__ == "__main__":
    main()
