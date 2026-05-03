#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: run_mafft_transcript_alignment.py

Description:
    Run batch multiple sequence alignment using MAFFT for transcript FASTA files.

    This script reads all FASTA files from an input directory, checks sequence
    number and sequence length, runs MAFFT alignment for each file, and verifies
    the aligned output.

Input:
    - Directory containing FASTA files

Output:
    - Aligned FASTA files with prefix "aligned_"

Usage:
    python run_mafft_transcript_alignment.py \
        --input_dir chromosome_grouped_transcripts \
        --output_dir aligned_transcripts \
        --mafft mafft \
        --threads -1
"""

import argparse
import subprocess
from pathlib import Path

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch MAFFT alignment for transcript FASTA files."
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing input FASTA files."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to save aligned FASTA files."
    )

    parser.add_argument(
        "--mafft",
        default="mafft",
        help="Path to MAFFT executable. Default: mafft"
    )

    parser.add_argument(
        "--threads",
        default="-1",
        help="Number of MAFFT threads. Default: -1, using all available cores."
    )

    parser.add_argument(
        "--algorithm",
        default="--auto",
        help="MAFFT algorithm option. Default: --auto"
    )

    parser.add_argument(
        "--input_suffix",
        default=".fasta",
        help="Input FASTA suffix. Default: .fasta"
    )

    return parser.parse_args()


def check_mafft_installation(mafft_path):
    """
    Check whether MAFFT is available.
    """
    try:
        result = subprocess.run(
            [mafft_path, "--version"],
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            version_text = result.stdout.strip() or result.stderr.strip()
            print(f"MAFFT detected: {version_text}")
            return True

        print("MAFFT was found but did not run correctly.")
        print(result.stderr)
        return False

    except FileNotFoundError:
        print(f"MAFFT executable not found: {mafft_path}")
        return False

    except Exception as error:
        print(f"Error while checking MAFFT: {error}")
        return False


def test_mafft_simple(mafft_path, output_dir):
    """
    Run a small MAFFT test before processing real input files.
    """
    print("Running a simple MAFFT test...")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    test_input = output_dir / "test_input.fasta"
    test_output = output_dir / "test_output.fasta"

    test_sequences = """>seq1
ACGTACGTACGT
>seq2
ACGTACGTACGT
>seq3
ACGTACGTACGA
"""

    try:
        with open(test_input, "w", encoding="utf-8") as handle:
            handle.write(test_sequences)

        command = [
            mafft_path,
            "--auto",
            str(test_input)
        ]

        with open(test_output, "w", encoding="utf-8") as out_f:
            result = subprocess.run(
                command,
                stdout=out_f,
                stderr=subprocess.PIPE,
                text=True
            )

        if result.returncode == 0 and test_output.exists():
            print("MAFFT test completed successfully.")
            test_input.unlink(missing_ok=True)
            test_output.unlink(missing_ok=True)
            return True

        print("MAFFT test failed.")
        print(result.stderr)
        return False

    except Exception as error:
        print(f"MAFFT test failed: {error}")
        return False


def analyze_sequences(input_file):
    """
    Print basic information about input FASTA sequences.
    """
    try:
        sequences = list(SeqIO.parse(input_file, "fasta"))
        print(f"  Number of sequences: {len(sequences)}")

        if sequences:
            lengths = [len(seq.seq) for seq in sequences]
            print(f"  Sequence length range: {min(lengths)} - {max(lengths)} bp")
            print(f"  Mean sequence length: {sum(lengths) / len(lengths):.1f} bp")
            print(f"  Example IDs: {', '.join([seq.id for seq in sequences[:3]])}")

        return len(sequences)

    except Exception as error:
        print(f"  Failed to analyze sequence file: {error}")
        return 0


def run_mafft_alignment(input_file, output_file, mafft_path, algorithm, threads):
    """
    Run MAFFT alignment for one FASTA file.
    """
    command = [
        mafft_path,
        algorithm,
        "--thread",
        str(threads),
        str(input_file)
    ]

    print(f"  Running command: {' '.join(command)}")

    try:
        with open(output_file, "w", encoding="utf-8") as out_f:
            subprocess.run(
                command,
                stdout=out_f,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )

        print(f"  Alignment completed: {input_file.name}")
        return True

    except subprocess.CalledProcessError as error:
        print(f"  MAFFT alignment failed for {input_file.name}")
        print(error.stderr)
        return False

    except Exception as error:
        print(f"  Unexpected error during MAFFT alignment: {error}")
        return False


def verify_alignment(output_file):
    """
    Verify whether aligned FASTA output is valid.
    """
    if not output_file.exists() or output_file.stat().st_size == 0:
        print("  Output file was not created or is empty.")
        return False

    try:
        aligned_sequences = list(SeqIO.parse(output_file, "fasta"))
        print(f"  Number of aligned sequences: {len(aligned_sequences)}")

        if aligned_sequences:
            aligned_lengths = [len(seq.seq) for seq in aligned_sequences]

            if len(set(aligned_lengths)) == 1:
                print(f"  Aligned sequence length: {aligned_lengths[0]} bp")
            else:
                print(
                    "  Warning: aligned sequence lengths are not identical: "
                    f"{min(aligned_lengths)} - {max(aligned_lengths)} bp"
                )

        print(f"  Output saved to: {output_file}")
        return True

    except Exception as error:
        print(f"  Failed to verify output file: {error}")
        return False


def process_all_files(args):
    """
    Process all FASTA files in the input directory.
    """
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    if not check_mafft_installation(args.mafft):
        raise RuntimeError("MAFFT is not available. Please install MAFFT or check the path.")

    if not test_mafft_simple(args.mafft, output_dir):
        raise RuntimeError("MAFFT test failed. Please check MAFFT installation.")

    fasta_files = sorted(input_dir.glob(f"*{args.input_suffix}"))

    if not fasta_files:
        print(f"No FASTA files ending with '{args.input_suffix}' found in: {input_dir}")
        return

    print(f"Detected {len(fasta_files)} FASTA file(s).")

    success_count = 0
    failed_files = []

    for input_file in fasta_files:
        output_file = output_dir / f"aligned_{input_file.name}"

        print("\n" + "=" * 60)
        print(f"Processing file: {input_file.name}")
        print("=" * 60)

        seq_count = analyze_sequences(input_file)

        if seq_count < 2:
            print("  Warning: fewer than two sequences. Skipping alignment.")
            failed_files.append((input_file.name, "fewer than two sequences"))
            continue

        success = run_mafft_alignment(
            input_file=input_file,
            output_file=output_file,
            mafft_path=args.mafft,
            algorithm=args.algorithm,
            threads=args.threads
        )

        if success and verify_alignment(output_file):
            success_count += 1
        else:
            failed_files.append((input_file.name, "MAFFT alignment or output verification failed"))

    print("\n" + "=" * 60)
    print("Batch MAFFT alignment completed")
    print("=" * 60)
    print(f"Successfully aligned files: {success_count}")

    if failed_files:
        print("Failed files:")
        for filename, reason in failed_files:
            print(f"  - {filename}: {reason}")


def main():
    args = parse_args()

    input_dir = Path(args.input_dir)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    print("=" * 60)
    print("MAFFT transcript alignment pipeline")
    print("=" * 60)
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"MAFFT executable: {args.mafft}")
    print(f"Algorithm: {args.algorithm}")
    print(f"Threads: {args.threads}")

    process_all_files(args)


if __name__ == "__main__":
    main()
