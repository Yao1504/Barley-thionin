#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: run_thionin_pangenome_blastp_pipeline.py

Description:
    Run BLASTP searches using barley Morex thionin protein sequences as queries against
    multiple pangenome protein FASTA files, filter candidate homologs, and
    extract matched protein sequences.

Input:
    1. Query protein FASTA file
       Example: thionin_proteins.fasta

    2. Directory containing pangenome protein FASTA files
       Supported extensions: .fa, .fasta, .faa

Output:
    1. Raw BLASTP tabular results for each genome
    2. Filtered BLASTP hit table for each genome
    3. Filtered protein FASTA file for each genome
    4. Combined summary table across all genomes

Usage:
    python blastp_pangenome_pipeline.py \
        --query thionin_proteins.fasta \
        --pangenome_dir pangenome_protein \
        --output_dir output_pangenome_blastp \
        --evalue 1e-5 \
        --max_mismatch 10 \
        --min_alignment 70

Requirements:
    - NCBI BLAST+
    - Python >= 3.8
    - Biopython
    - pandas
"""

import argparse
import os
import subprocess
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
        description="Run BLASTP against multiple pangenome protein FASTA files and extract filtered hits."
    )

    parser.add_argument(
        "--query",
        required=True,
        help="Query protein FASTA file, e.g. thionin protein sequences."
    )

    parser.add_argument(
        "--pangenome_dir",
        required=True,
        help="Directory containing pangenome protein FASTA files."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory."
    )

    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="BLASTP E-value cutoff. Default: 1e-5."
    )

    parser.add_argument(
        "--max_mismatch",
        type=int,
        default=10,
        help="Maximum allowed mismatch number. Default: 10."
    )

    parser.add_argument(
        "--min_alignment",
        type=int,
        default=70,
        help="Minimum alignment length. Default: 70."
    )

    parser.add_argument(
        "--max_target_seqs",
        type=int,
        default=5000,
        help="Maximum target sequences reported by BLASTP. Default: 5000."
    )

    parser.add_argument(
        "--force_rebuild_db",
        action="store_true",
        help="Rebuild BLAST databases even if database files already exist."
    )

    return parser.parse_args()


def check_executable(command):
    try:
        subprocess.run(
            [command, "-version"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False
        )
    except FileNotFoundError:
        sys.exit(
            f"Error: '{command}' was not found. "
            "Please install NCBI BLAST+ and add it to PATH."
        )


def find_fasta_files(folder):
    fasta_extensions = {".fa", ".fasta", ".faa"}
    fasta_files = [
        file for file in Path(folder).iterdir()
        if file.is_file() and file.suffix.lower() in fasta_extensions
    ]

    if not fasta_files:
        sys.exit(f"Error: no FASTA files found in {folder}")

    return sorted(fasta_files)


def blast_db_exists(db_prefix):
    required_suffixes = [".pin", ".phr", ".psq"]
    return all(Path(str(db_prefix) + suffix).exists() for suffix in required_suffixes)


def run_makeblastdb(fasta_file, db_prefix, force_rebuild=False):
    if blast_db_exists(db_prefix) and not force_rebuild:
        print(f"BLAST database already exists, skipping: {db_prefix}")
        return

    command = [
        "makeblastdb",
        "-in", str(fasta_file),
        "-dbtype", "prot",
        "-out", str(db_prefix)
    ]

    print(f"Building BLAST database: {fasta_file.name}")
    subprocess.run(command, check=True)


def run_blastp(query_file, db_prefix, output_file, evalue, max_target_seqs):
    command = [
        "blastp",
        "-query", str(query_file),
        "-db", str(db_prefix),
        "-out", str(output_file),
        "-outfmt", "6",
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs)
    ]

    print(f"Running BLASTP: {output_file.name}")
    subprocess.run(command, check=True)


def load_and_filter_blast(blast_file, evalue, max_mismatch, min_alignment):
    if os.path.getsize(blast_file) == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    df = pd.read_csv(
        blast_file,
        sep="\t",
        header=None,
        names=BLAST_COLUMNS
    )

    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
    df["length"] = pd.to_numeric(df["length"], errors="coerce")
    df["mismatch"] = pd.to_numeric(df["mismatch"], errors="coerce")
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")

    filtered = df[
        (df["evalue"] <= evalue) &
        (df["mismatch"] <= max_mismatch) &
        (df["length"] >= min_alignment)
    ].copy()

    filtered.sort_values(
        by=["sseqid", "evalue", "bitscore"],
        ascending=[True, True, False],
        inplace=True
    )

    filtered = filtered.drop_duplicates(
        subset=["sseqid"],
        keep="first"
    )

    return filtered


def extract_sequences(genome_fasta, selected_ids, output_fasta):
    seq_dict = SeqIO.to_dict(SeqIO.parse(str(genome_fasta), "fasta"))

    written = 0
    missing_ids = []

    with open(output_fasta, "w", encoding="utf-8") as out_f:
        for seq_id in sorted(selected_ids):
            if seq_id in seq_dict:
                SeqIO.write(seq_dict[seq_id], out_f, "fasta")
                written += 1
            else:
                missing_ids.append(seq_id)

    return written, missing_ids


def safe_stem(file_path):
    name = file_path.name
    for suffix in [".fasta", ".faa", ".fa"]:
        if name.endswith(suffix):
            return name.replace(suffix, "")
    return file_path.stem


def main():
    args = parse_args()

    query_file = Path(args.query)
    pangenome_dir = Path(args.pangenome_dir)
    output_dir = Path(args.output_dir)

    if not query_file.exists():
        sys.exit(f"Error: query FASTA file not found: {query_file}")

    if not pangenome_dir.exists():
        sys.exit(f"Error: pangenome directory not found: {pangenome_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    raw_dir = output_dir / "raw_blastp"
    filtered_dir = output_dir / "filtered_tables"
    fasta_dir = output_dir / "filtered_fastas"
    db_dir = output_dir / "blast_databases"

    raw_dir.mkdir(exist_ok=True)
    filtered_dir.mkdir(exist_ok=True)
    fasta_dir.mkdir(exist_ok=True)
    db_dir.mkdir(exist_ok=True)

    check_executable("makeblastdb")
    check_executable("blastp")

    genome_fastas = find_fasta_files(pangenome_dir)

    summary_records = []
    all_filtered_tables = []

    for genome_fasta in genome_fastas:
        genome_name = safe_stem(genome_fasta)

        print("")
        print("=" * 70)
        print(f"Processing genome: {genome_name}")
        print("=" * 70)

        db_prefix = db_dir / f"{genome_name}_blastdb"
        raw_blast_out = raw_dir / f"{genome_name}_blastp.tsv"
        filtered_csv = filtered_dir / f"{genome_name}_blastp_filtered.csv"
        filtered_fasta = fasta_dir / f"{genome_name}_blastp_filtered.fasta"
        missing_ids_file = filtered_dir / f"{genome_name}_missing_ids.txt"

        run_makeblastdb(
            fasta_file=genome_fasta,
            db_prefix=db_prefix,
            force_rebuild=args.force_rebuild_db
        )

        run_blastp(
            query_file=query_file,
            db_prefix=db_prefix,
            output_file=raw_blast_out,
            evalue=args.evalue,
            max_target_seqs=args.max_target_seqs
        )

        filtered_df = load_and_filter_blast(
            blast_file=raw_blast_out,
            evalue=args.evalue,
            max_mismatch=args.max_mismatch,
            min_alignment=args.min_alignment
        )

        filtered_df["genome"] = genome_name
        filtered_df.to_csv(filtered_csv, index=False)

        selected_ids = set(filtered_df["sseqid"].astype(str))

        written_count, missing_ids = extract_sequences(
            genome_fasta=genome_fasta,
            selected_ids=selected_ids,
            output_fasta=filtered_fasta
        )

        if missing_ids:
            with open(missing_ids_file, "w", encoding="utf-8") as f:
                for seq_id in missing_ids:
                    f.write(f"{seq_id}\n")

        summary_records.append({
            "genome": genome_name,
            "raw_blast_file": str(raw_blast_out),
            "filtered_hit_count": len(filtered_df),
            "fasta_sequence_count": written_count,
            "missing_id_count": len(missing_ids),
            "filtered_csv": str(filtered_csv),
            "filtered_fasta": str(filtered_fasta)
        })

        if not filtered_df.empty:
            all_filtered_tables.append(filtered_df)

        print(f"Filtered hits: {len(filtered_df)}")
        print(f"Sequences written: {written_count}")

    summary_df = pd.DataFrame(summary_records)
    summary_out = output_dir / "pangenome_blastp_summary.csv"
    summary_df.to_csv(summary_out, index=False)

    if all_filtered_tables:
        combined_df = pd.concat(all_filtered_tables, ignore_index=True)
    else:
        combined_df = pd.DataFrame(columns=BLAST_COLUMNS + ["genome"])

    combined_out = output_dir / "pangenome_blastp_all_filtered_hits.csv"
    combined_df.to_csv(combined_out, index=False)

    print("")
    print("=" * 70)
    print("Pangenome BLASTP pipeline finished.")
    print("=" * 70)
    print(f"Summary table: {summary_out}")
    print(f"Combined filtered hits: {combined_out}")
    print(f"Filtered FASTA folder: {fasta_dir}")


if __name__ == "__main__":
    main()
