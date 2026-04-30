#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: calculate_pairwise_identity.py

Description:
    Calculate pairwise sequence identity (%) from an aligned FASTA file
    and export the result as a CSV matrix.

    This script can be used for protein or nucleotide/CDS alignments.

Input:
    - Aligned FASTA file

Output:
    - Pairwise identity matrix in CSV format

Usage:
    python calculate_pairwise_identity.py \
        --input aligned_sequences.fasta \
        --output identity_matrix.csv

Example:
    python calculate_pairwise_identity.py \
        --input protein_aligned.fasta \
        --output protein_identity_matrix.csv
"""

import argparse
import sys
import numpy as np
import pandas as pd
from Bio import AlignIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate pairwise sequence identity from an aligned FASTA file."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input aligned FASTA file. USER SET THIS PATH."
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output CSV file for the identity matrix. USER SET THIS PATH."
    )

    return parser.parse_args()


def calculate_identity(seq1, seq2):
    """
    Calculate pairwise identity between two aligned sequences.

    Gap positions are ignored when either sequence contains '-'.

    Formula:
        identity (%) = identical non-gap positions / comparable non-gap positions * 100
    """
    matches = 0
    comparable_sites = 0

    for a, b in zip(seq1, seq2):
        if a != "-" and b != "-":
            comparable_sites += 1
            if a == b:
                matches += 1

    if comparable_sites == 0:
        return 0.0

    return round(matches / comparable_sites * 100, 2)


def pairwise_identity_matrix(alignment_file, output_csv):
    """
    Calculate pairwise identity matrix and write it to CSV.
    """
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
    except Exception as error:
        sys.exit(f"Error: failed to read alignment file: {alignment_file}\n{error}")

    n = len(alignment)

    if n == 0:
        sys.exit("Error: no sequences found in the alignment file.")

    matrix = np.zeros((n, n))

    for i in range(n):
        seq1 = str(alignment[i].seq)

        for j in range(n):
            seq2 = str(alignment[j].seq)
            matrix[i, j] = calculate_identity(seq1, seq2)

    sequence_names = [record.id for record in alignment]

    df = pd.DataFrame(
        matrix,
        index=sequence_names,
        columns=sequence_names
    )

    df.to_csv(output_csv)

    print("Pairwise identity analysis completed.")
    print(f"Number of sequences: {n}")
    print(f"Identity matrix saved to: {output_csv}")


def main():
    args = parse_args()

    pairwise_identity_matrix(
        alignment_file=args.input,
        output_csv=args.output
    )


if __name__ == "__main__":
    main()
