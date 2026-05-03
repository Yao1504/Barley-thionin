#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: classify_transcripts_by_chr.py

Description:
    Classify final transcript sequences by chromosome based on sequence IDs.

    This script reads FASTA files ending with "_final_sequences.fasta",
    extracts chromosome information from sequence IDs, groups transcripts
    by chromosome, and writes chromosome-specific FASTA files.

    Unknown or unanchored sequences are assigned to "Un".

    In addition, the script generates combined FASTA files:
        - seq_6H+Un.fasta
        - seq_7H+Un.fasta

Input:
    Directory containing final transcript FASTA files:
        <genotype>_final_sequences.fasta

Output:
    Chromosome-grouped FASTA files:
        seqs_1H.fasta
        seqs_2H.fasta
        ...
        seqs_7H.fasta
        seq_6H+Un.fasta
        seq_7H+Un.fasta

Usage:
    python classify_transcripts_by_chr.py \
        --input_dir final_sequences \
        --output_dir chromosome_grouped_transcripts
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Classify transcript FASTA sequences by chromosome."
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing *_final_sequences.fasta files."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to save chromosome-classified FASTA files."
    )

    parser.add_argument(
        "--input_suffix",
        default="_final_sequences.fasta",
        help="Suffix of input FASTA files. Default: _final_sequences.fasta"
    )

    return parser.parse_args()


def extract_chromosome_info(seq_id):
    """
    Extract chromosome information from a sequence ID.

    Recognized patterns:
        chr1H, chr2H, ..., chr7H  -> 1H, 2H, ..., 7H
        contig / CAJHDD           -> Un

    Returns:
        Chromosome name such as 1H, 2H, ..., 7H, Un, or Unknown.
    """
    chr_match = re.search(r"chr(\dH)", seq_id, re.IGNORECASE)

    if chr_match:
        return chr_match.group(1)

    if re.search(r"contig", seq_id, re.IGNORECASE):
        return "Un"

    if re.search(r"CAJHDD", seq_id, re.IGNORECASE):
        return "Un"

    if re.search(r"_[Cc]ontig_", seq_id):
        return "Un"

    return "Unknown"


def analyze_chromosome_patterns(input_dir, input_suffix):
    """
    Analyze chromosome patterns from a subset of sequence IDs.
    This is used as a diagnostic step.
    """
    print("\nChromosome pattern diagnostics:")

    fasta_files = sorted(Path(input_dir).glob(f"*{input_suffix}"))

    patterns = defaultdict(list)

    for fasta_file in fasta_files[:3]:
        try:
            sequences = list(SeqIO.parse(fasta_file, "fasta"))

            for seq in sequences[:10]:
                chrom = extract_chromosome_info(seq.id)
                patterns[chrom].append(seq.id)

        except Exception as error:
            print(f"  Warning: failed to inspect {fasta_file.name}: {error}")

    if not patterns:
        print("  No sequence IDs were inspected.")
        return

    for chrom, ids in sorted(patterns.items()):
        print(f"  {chrom}: {len(ids)} example sequence(s)")

        if ids:
            print(f"    Example: {ids[0]}")


def process_transcript_sequences(input_dir, output_dir, input_suffix):
    """
    Read final transcript FASTA files, classify sequences by chromosome,
    and write chromosome-specific FASTA outputs.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    chromosome_seqs = defaultdict(list)

    fasta_files = sorted(input_dir.glob(f"*{input_suffix}"))

    if not fasta_files:
        print(f"No FASTA files ending with '{input_suffix}' were found in: {input_dir}")
        return chromosome_seqs

    print(f"\nDetected {len(fasta_files)} FASTA file(s).")

    total_sequences = 0

    for fasta_file in fasta_files:
        genotype = fasta_file.name.replace(input_suffix, "")
        print(f"Processing genotype: {genotype}")

        try:
            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            print(f"  Loaded {len(sequences)} sequence(s).")

            for seq in sequences:
                chromosome = extract_chromosome_info(seq.id)
                chromosome_seqs[chromosome].append(seq)
                total_sequences += 1

        except Exception as error:
            print(f"  Warning: failed to process {fasta_file.name}: {error}")

    print(f"\nTotal sequences processed: {total_sequences}")

    print("Chromosome distribution before merging Unknown:")
    for chrom, seqs in sorted(chromosome_seqs.items()):
        print(f"  {chrom}: {len(seqs)} sequence(s)")

    if "Unknown" in chromosome_seqs:
        print(
            f"\nMerging {len(chromosome_seqs['Unknown'])} Unknown sequence(s) into Un."
        )

        chromosome_seqs["Un"].extend(chromosome_seqs["Unknown"])
        del chromosome_seqs["Unknown"]

    known_chromosomes = ["1H", "2H", "3H", "4H", "5H", "6H", "7H"]

    for chrom in known_chromosomes:
        if chrom in chromosome_seqs and chromosome_seqs[chrom]:
            output_file = output_dir / f"seqs_{chrom}.fasta"
            SeqIO.write(chromosome_seqs[chrom], output_file, "fasta")
            print(
                f"Created {output_file.name}: "
                f"{len(chromosome_seqs[chrom])} sequence(s)"
            )

    un_sequences = chromosome_seqs.get("Un", [])

    seqs_6h_un = chromosome_seqs.get("6H", []) + un_sequences

    if seqs_6h_un:
        output_file = output_dir / "seq_6H+Un.fasta"
        SeqIO.write(seqs_6h_un, output_file, "fasta")
        print(f"Created {output_file.name}: {len(seqs_6h_un)} sequence(s)")

    seqs_7h_un = chromosome_seqs.get("7H", []) + un_sequences

    if seqs_7h_un:
        output_file = output_dir / "seq_7H+Un.fasta"
        SeqIO.write(seqs_7h_un, output_file, "fasta")
        print(f"Created {output_file.name}: {len(seqs_7h_un)} sequence(s)")

    print(f"\nAll output files were saved to: {output_dir}")

    return chromosome_seqs


def verify_results(output_dir):
    """
    Verify generated FASTA files by counting sequences.
    """
    output_dir = Path(output_dir)

    print("\nVerifying output files:")

    output_files = sorted(output_dir.glob("*.fasta"))

    if not output_files:
        print("  No output FASTA files found.")
        return

    for file_path in output_files:
        try:
            sequences = list(SeqIO.parse(file_path, "fasta"))

            print(f"  {file_path.name}: {len(sequences)} sequence(s)")

            if sequences:
                sample_ids = [seq.id for seq in sequences[:2]]
                print(f"    Example IDs: {', '.join(sample_ids)}")

        except Exception as error:
            print(f"  Warning: failed to read {file_path.name}: {error}")


def main():
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    print("=" * 70)
    print("Transcript chromosome classification pipeline")
    print("=" * 70)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    analyze_chromosome_patterns(
        input_dir=input_dir,
        input_suffix=args.input_suffix
    )

    process_transcript_sequences(
        input_dir=input_dir,
        output_dir=output_dir,
        input_suffix=args.input_suffix
    )

    verify_results(output_dir)

    print("\nPipeline completed successfully.")


if __name__ == "__main__":
    main()
