#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: run_mafft_batch_alignment.py

Description:
    Perform batch multiple sequence alignment using MAFFT for all FASTA files
    in a given directory.

    Each input FASTA file is aligned independently using MAFFT (L-INS-i mode),
    and the aligned sequences are written to the output directory.

Input:
    - Directory containing FASTA files

Output:
    - Aligned FASTA files with suffix "_aligned.fasta"

Usage:
    python run_mafft_batch_alignment.py \
        --input_dir input_fasta_dir \
        --output_dir output_alignment_dir \
        --mafft mafft

Requirements:
    - MAFFT installed and accessible in PATH
    - Python >= 3.8
"""

import argparse
import subprocess
from pathlib import Path
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch MAFFT alignment for FASTA files."
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing input FASTA files. USER SET THIS PATH."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to save aligned FASTA files. USER SET THIS PATH."
    )

    parser.add_argument(
        "--mafft",
        default="mafft",
        help="Path to MAFFT executable (default: mafft in PATH)"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for MAFFT (default: 1)"
    )

    return parser.parse_args()


def check_mafft(mafft_cmd):
    try:
        subprocess.run(
            [mafft_cmd, "--help"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
    except FileNotFoundError:
        sys.exit(f"Error: MAFFT not found: {mafft_cmd}")


def main():
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    if not input_dir.exists():
        sys.exit(f"Error: input directory not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    check_mafft(args.mafft)

    fasta_files = [
        f for f in input_dir.iterdir()
        if f.suffix.lower() in [".fa", ".fasta", ".faa", ".fna"]
    ]

    if not fasta_files:
        sys.exit(f"No FASTA files found in {input_dir}")

    print(f"Found {len(fasta_files)} FASTA files.\n")

    for fasta in fasta_files:
        output_file = output_dir / f"{fasta.stem}_aligned.fasta"

        print(f"Running MAFFT for: {fasta.name}")

        cmd = [
            args.mafft,
            "--localpair",
            "--maxiterate", "1000",
            "--thread", str(args.threads),
            str(fasta)
        ]

        try:
            with open(output_file, "w") as out_f:
                subprocess.run(
                    cmd,
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )

            print(f"Saved: {output_file}\n")

        except subprocess.CalledProcessError as e:
            print(f"Error running MAFFT on {fasta.name}")
            print(e.stderr)

    print("All alignments completed successfully.")


if __name__ == "__main__":
    main()
