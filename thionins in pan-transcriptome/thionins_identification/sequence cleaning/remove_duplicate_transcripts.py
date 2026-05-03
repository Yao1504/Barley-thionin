#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: remove_duplicate_transcripts.py

Description:
    Remove redundant transcript sequences from FASTA files.

    For each input FASTA file, transcripts are grouped by base gene ID.
    Within the same gene, transcripts with identical sequence content are
    considered redundant, and only the first transcript is retained.

    Gap characters "-" are ignored when comparing sequence identity.

Input:
    - Directory containing FASTA files

Output:
    - Deduplicated FASTA files with the same filenames

Usage:
    python remove_duplicate_transcripts.py \
        --input_dir input_fasta \
        --output_dir output_deduplicated

Optional:
    python remove_duplicate_transcripts.py \
        --input_dir input_fasta \
        --output_dir output_deduplicated \
        --target_gene GeneID
"""

import argparse
import hashlib
import re
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove duplicate transcript sequences within the same gene."
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing input FASTA files."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to save deduplicated FASTA files."
    )

    parser.add_argument(
        "--input_suffix",
        default=".fasta",
        help="Input FASTA suffix. Default: .fasta"
    )

    parser.add_argument(
        "--target_gene",
        default=None,
        help="Optional gene ID for duplicate diagnostics."
    )

    return parser.parse_args()


def get_base_gene_id(seq_id):
    """
    Extract base gene ID from sequence ID.

    This removes:
        1. sequence position suffix such as /47-457
        2. isoform suffix such as .1, .2, .3

    Example:
        GeneA.1/47-457 -> GeneA
    """
    base_id = re.sub(r"/\d+-\d+$", "", seq_id)
    base_id = re.sub(r"\.\d+$", "", base_id)

    return base_id


def get_sequence_hash(sequence):
    """
    Calculate MD5 hash of a sequence after removing gap characters.
    """
    clean_sequence = str(sequence).replace("-", "")
    return hashlib.md5(clean_sequence.encode()).hexdigest()


def extract_isoform_number(seq_id):
    """
    Extract isoform number from sequence ID for sorting.
    """
    match = re.search(r"\.(\d+)", seq_id)
    return int(match.group(1)) if match else 9999


def remove_duplicate_transcripts(fasta_file, output_file):
    """
    Remove duplicate transcript sequences within the same gene.
    """
    gene_sequences = defaultdict(list)

    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    print(f"  Loaded {len(sequences)} sequence(s).")

    for seq in sequences:
        base_id = get_base_gene_id(seq.id)
        gene_sequences[base_id].append(seq)

    unique_sequences = []
    removed_count = 0

    for base_id, seq_list in gene_sequences.items():
        if len(seq_list) == 1:
            unique_sequences.append(seq_list[0])
            continue

        seen_sequences = {}

        seq_list_sorted = sorted(
            seq_list,
            key=lambda record: extract_isoform_number(record.id)
        )

        for seq in seq_list_sorted:
            seq_hash = get_sequence_hash(seq.seq)

            if seq_hash in seen_sequences:
                removed_count += 1
                print(
                    f"    Removed duplicate: {seq.id} "
                    f"(same sequence as {seen_sequences[seq_hash].id})"
                )
                continue

            seen_sequences[seq_hash] = seq
            unique_sequences.append(seq)

    SeqIO.write(unique_sequences, output_file, "fasta")

    return removed_count, len(unique_sequences)


def analyze_target_gene(input_dir, input_suffix, target_gene):
    """
    Analyze duplicated transcripts for one target gene.
    """
    if not target_gene:
        return

    print("\nTarget gene duplicate diagnostics:")
    print(f"Target gene: {target_gene}")

    fasta_files = sorted(Path(input_dir).glob(f"*{input_suffix}"))

    for fasta_file in fasta_files:
        print(f"\nInspecting file: {fasta_file.name}")

        try:
            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            gene_sequences = defaultdict(list)

            for seq in sequences:
                base_id = get_base_gene_id(seq.id)
                gene_sequences[base_id].append(seq)

            if target_gene not in gene_sequences:
                print("  Target gene not found.")
                continue

            print(
                f"  Found {len(gene_sequences[target_gene])} transcript(s) "
                f"for {target_gene}"
            )

            seen_hashes = {}
            duplicates = []

            for seq in gene_sequences[target_gene]:
                seq_hash = get_sequence_hash(seq.seq)
                seq_string = str(seq.seq)

                print(
                    f"    - {seq.id}: length={len(seq.seq)}, "
                    f"hash={seq_hash[:8]}..."
                )
                print(f"      sequence preview: {seq_string[:50]}...{seq_string[-50:]}")

                if seq_hash in seen_hashes:
                    duplicates.append((seq.id, seen_hashes[seq_hash]))
                else:
                    seen_hashes[seq_hash] = seq.id

            if duplicates:
                print("  Duplicate transcripts detected:")
                for duplicate_id, original_id in duplicates:
                    print(f"    - {duplicate_id} duplicates {original_id}")
            else:
                print("  No duplicate transcripts detected.")

        except Exception as error:
            print(f"  Warning: failed to inspect {fasta_file.name}: {error}")


def verify_results(output_dir, input_suffix):
    """
    Verify that no identical transcript sequences remain within each gene.
    """
    print("\nVerifying output files:")

    output_files = sorted(Path(output_dir).glob(f"*{input_suffix}"))

    if not output_files:
        print("  No output FASTA files found.")
        return

    for output_file in output_files:
        try:
            sequences = list(SeqIO.parse(output_file, "fasta"))
            gene_sequences = defaultdict(list)

            for seq in sequences:
                base_id = get_base_gene_id(seq.id)
                gene_sequences[base_id].append(seq)

            duplicate_found = False

            for gene_id, seqs in gene_sequences.items():
                if len(seqs) > 1:
                    seq_hashes = [get_sequence_hash(seq.seq) for seq in seqs]

                    if len(set(seq_hashes)) != len(seq_hashes):
                        duplicate_found = True
                        print(
                            f"  Warning: duplicated sequence still found "
                            f"for gene {gene_id} in {output_file.name}"
                        )

            if not duplicate_found:
                print(f"  {output_file.name}: verified ({len(sequences)} unique sequence(s))")
            else:
                print(f"  {output_file.name}: verification failed")

        except Exception as error:
            print(f"  Warning: failed to verify {output_file.name}: {error}")


def process_all_files(input_dir, output_dir, input_suffix):
    """
    Process all FASTA files in the input directory.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(input_dir.glob(f"*{input_suffix}"))

    if not fasta_files:
        print(f"No FASTA files ending with '{input_suffix}' found in: {input_dir}")
        return 0, 0

    print(f"Detected {len(fasta_files)} FASTA file(s).")

    total_removed = 0
    total_kept = 0

    for input_file in fasta_files:
        output_file = output_dir / input_file.name

        print("\n" + "=" * 60)
        print(f"Processing file: {input_file.name}")
        print("=" * 60)

        try:
            removed, kept = remove_duplicate_transcripts(
                fasta_file=input_file,
                output_file=output_file
            )

            total_removed += removed
            total_kept += kept

            print(f"  Removed duplicate sequence(s): {removed}")
            print(f"  Retained unique sequence(s): {kept}")
            print(f"  Output: {output_file}")

        except Exception as error:
            print(f"  Warning: failed to process {input_file.name}: {error}")

    return total_removed, total_kept


def main():
    args = parse_args()

    input_dir = Path(args.input_dir)

    print("=" * 70)
    print("Duplicate transcript removal pipeline")
    print("=" * 70)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    analyze_target_gene(
        input_dir=args.input_dir,
        input_suffix=args.input_suffix,
        target_gene=args.target_gene
    )

    total_removed, total_kept = process_all_files(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        input_suffix=args.input_suffix
    )

    verify_results(
        output_dir=args.output_dir,
        input_suffix=args.input_suffix
    )

    total = total_removed + total_kept

    print("\n" + "=" * 70)
    print("Pipeline completed")
    print("=" * 70)
    print(f"Total duplicate sequences removed: {total_removed}")
    print(f"Total unique sequences retained: {total_kept}")

    if total > 0:
        print(f"Deduplication rate: {total_removed / total * 100:.2f}%")


if __name__ == "__main__":
    main()
