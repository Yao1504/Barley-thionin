#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: run_thionin_pangenome_hmmer_pipeline.py

Description:
    Run HMMER hmmsearch to identify thionin proteins (PF00321)
    across multiple pangenome protein FASTA files.

Input:
    1. HMM profile (e.g. PF00321.hmm)
    2. Directory of protein FASTA files

Output:
    - Raw hmmsearch tblout files
    - Summary table of hits per genome

Usage:
    python run_thionin_hmmsearch_pangenome_pipeline.py \
        --hmm PF00321.hmm \
        --pangenome_dir input/pangenome_protein \
        --output_dir output_hmmsearch \
        --evalue 1e-5
"""

import argparse
import subprocess
import sys
from pathlib import Path
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run HMMER hmmsearch on multiple pangenome FASTA files."
    )

    parser.add_argument("--hmm", required=True, help="HMM profile file (e.g. PF00321.hmm)")
    parser.add_argument("--pangenome_dir", required=True, help="Directory of protein FASTA files")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value threshold")

    return parser.parse_args()


def check_hmmer():
    try:
        subprocess.run(["hmmsearch", "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        sys.exit("Error: hmmsearch not found. Install HMMER and add to PATH.")


def find_fasta_files(folder):
    return [
        f for f in Path(folder).iterdir()
        if f.suffix.lower() in [".fa", ".fasta", ".faa"]
    ]


def run_hmmsearch(hmm_file, fasta_file, output_file, evalue):
    cmd = [
        "hmmsearch",
        "-E", str(evalue),
        "--tblout", str(output_file),
        str(hmm_file),
        str(fasta_file)
    ]

    print(f"Running hmmsearch on {fasta_file.name}")
    subprocess.run(cmd, check=True)


def parse_tblout(tbl_file):
    hits = []
    with open(tbl_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            hits.append({
                "target": parts[0],
                "query": parts[2],
                "evalue": float(parts[4]),
                "score": float(parts[5])
            })
    return pd.DataFrame(hits)


def main():
    args = parse_args()

    hmm_file = Path(args.hmm)
    pangenome_dir = Path(args.pangenome_dir)
    output_dir = Path(args.output_dir)

    if not hmm_file.exists():
        sys.exit(f"Error: HMM file not found: {hmm_file}")

    if not pangenome_dir.exists():
        sys.exit(f"Error: directory not found: {pangenome_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    raw_dir = output_dir / "raw_tblout"
    raw_dir.mkdir(exist_ok=True)

    check_hmmer()

    fasta_files = find_fasta_files(pangenome_dir)

    summary = []

    for fasta in fasta_files:
        genome = fasta.stem
        out_tbl = raw_dir / f"{genome}_hmmsearch.tbl"

        run_hmmsearch(hmm_file, fasta, out_tbl, args.evalue)

        df = parse_tblout(out_tbl)
        hit_count = len(df)

        summary.append({
            "genome": genome,
            "hits": hit_count,
            "tblout": str(out_tbl)
        })

        print(f"{genome}: {hit_count} hits")

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_dir / "hmmsearch_summary.csv", index=False)

    print("\nPipeline finished.")
    print(f"Summary saved to: {output_dir / 'hmmsearch_summary.csv'}")


if __name__ == "__main__":
    main()
