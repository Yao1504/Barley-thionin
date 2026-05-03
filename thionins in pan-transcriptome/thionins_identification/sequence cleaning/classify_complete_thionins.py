#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: classify_complete_thionins.py

Description:
    Classify thionin candidate sequences as complete or incomplete based on:
        1. strict ORF detection
        2. ORF length
        3. conserved cysteine motif
        4. predicted signal peptide-like hydrophobic N-terminus

Input:
    - Directory containing FASTA files

Output:
    - complete.txt
    - lack.txt
    - report.csv

Usage:
    python classify_complete_thionins.py \
        --input_dir input_fasta \
        --output_dir output_completeness
"""

import argparse
import csv
import re
from pathlib import Path


CYS_MOTIF = re.compile(r"C.{2,20}C.{2,20}C.{2,20}C")


CODON_TABLE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "G", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Classify thionin candidates as complete or incomplete."
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing input FASTA files."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to save complete.txt, lack.txt, and report.csv."
    )

    parser.add_argument(
        "--input_suffix",
        default=".fasta",
        help="Input FASTA suffix. Default: .fasta"
    )

    parser.add_argument(
        "--min_complete_orf",
        type=int,
        default=380,
        help="Minimum ORF length for complete thionin classification. Default: 380"
    )

    parser.add_argument(
        "--hydrophobic_window",
        type=int,
        default=18,
        help="N-terminal amino acid window used for signal peptide-like check. Default: 18"
    )

    parser.add_argument(
        "--hydrophobic_min_count",
        type=int,
        default=8,
        help="Minimum hydrophobic residues in the N-terminal window. Default: 8"
    )

    return parser.parse_args()


def read_fasta(path):
    seqs = {}
    header = None
    seq_lines = []

    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq_lines)

                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)

        if header:
            seqs[header] = "".join(seq_lines)

    return seqs


def translate_dna(seq):
    aas = []

    seq = seq.upper()

    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]

        if len(codon) == 3:
            aas.append(CODON_TABLE.get(codon, "X"))

    return "".join(aas)


def find_strict_orf(seq):
    """
    Find the first strict ORF from ATG to an in-frame stop codon.

    Returns:
        tuple: (orf_nt, orf_aa, start, end, length)
        None: if no strict ORF is detected
    """
    seq = seq.upper()
    stops = {"TAA", "TAG", "TGA"}

    for frame in range(3):
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i + 3]

            if codon == "ATG":
                j = i + 3

                while j < len(seq) - 2:
                    stop = seq[j:j + 3]

                    if stop in stops:
                        nuc = seq[i:j + 3]
                        aa = translate_dna(nuc)
                        return nuc, aa, i, j + 3, j + 3 - i

                    j += 3

                return None

    return None


def has_signal_peptide_like_region(
    aa_seq,
    hydrophobic_window=18,
    hydrophobic_min_count=8
):
    """
    Simple signal peptide-like check based on N-terminal hydrophobicity.

    This is not a replacement for SignalP.
    """
    if len(aa_seq) < 25:
        return False

    hydrophobic = set("AILVFM")
    count = sum(aa in hydrophobic for aa in aa_seq[:hydrophobic_window])

    return count >= hydrophobic_min_count


def classify_sequences(args):
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    complete_ids = []
    incomplete_ids = []
    report_rows = []

    fasta_files = sorted(input_dir.glob(f"*{args.input_suffix}"))

    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files ending with '{args.input_suffix}' found in {input_dir}"
        )

    print(f"Detected {len(fasta_files)} FASTA file(s).")

    for fasta_file in fasta_files:
        print(f"Processing: {fasta_file.name}")

        seqs = read_fasta(fasta_file)

        for seq_id, seq in seqs.items():
            orf_info = find_strict_orf(seq)

            if orf_info is None:
                incomplete_ids.append(seq_id)
                report_rows.append({
                    "id": seq_id,
                    "source_file": fasta_file.name,
                    "orf_length": 0,
                    "cys_motif": False,
                    "signal_peptide_like": False,
                    "is_complete": False,
                    "reason": "no_strict_orf"
                })
                continue

            _, aa_orf, start, end, orf_len = orf_info

            has_cys = bool(CYS_MOTIF.search(aa_orf))

            has_signal = has_signal_peptide_like_region(
                aa_seq=aa_orf,
                hydrophobic_window=args.hydrophobic_window,
                hydrophobic_min_count=args.hydrophobic_min_count
            )

            is_complete = (
                orf_len >= args.min_complete_orf
                and has_cys
                and has_signal
            )

            if is_complete:
                complete_ids.append(seq_id)
                reason = "complete"
            else:
                incomplete_ids.append(seq_id)

                failed_reasons = []
                if orf_len < args.min_complete_orf:
                    failed_reasons.append("short_orf")
                if not has_cys:
                    failed_reasons.append("missing_cys_motif")
                if not has_signal:
                    failed_reasons.append("missing_signal_peptide_like_region")

                reason = ";".join(failed_reasons)

            report_rows.append({
                "id": seq_id,
                "source_file": fasta_file.name,
                "orf_length": orf_len,
                "cys_motif": has_cys,
                "signal_peptide_like": has_signal,
                "is_complete": is_complete,
                "reason": reason
            })

    return complete_ids, incomplete_ids, report_rows


def write_outputs(output_dir, complete_ids, incomplete_ids, report_rows):
    output_dir = Path(output_dir)

    complete_file = output_dir / "complete.txt"
    lack_file = output_dir / "lack.txt"
    report_file = output_dir / "report.csv"

    with open(complete_file, "w", encoding="utf-8") as handle:
        for seq_id in complete_ids:
            handle.write(seq_id + "\n")

    with open(lack_file, "w", encoding="utf-8") as handle:
        for seq_id in incomplete_ids:
            handle.write(seq_id + "\n")

    with open(report_file, "w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "id",
            "source_file",
            "orf_length",
            "cys_motif",
            "signal_peptide_like",
            "is_complete",
            "reason"
        ]

        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(report_rows)

    print(f"Complete thionins written to: {complete_file}")
    print(f"Incomplete thionins written to: {lack_file}")
    print(f"Detailed report written to: {report_file}")


def main():
    args = parse_args()

    print("=" * 70)
    print("Thionin completeness classification pipeline")
    print("=" * 70)

    complete_ids, incomplete_ids, report_rows = classify_sequences(args)

    write_outputs(
        output_dir=args.output_dir,
        complete_ids=complete_ids,
        incomplete_ids=incomplete_ids,
        report_rows=report_rows
    )

    print("\nSummary")
    print("-" * 70)
    print(f"Complete sequences: {len(complete_ids)}")
    print(f"Incomplete sequences: {len(incomplete_ids)}")
    print(f"Total sequences: {len(report_rows)}")
    print("Pipeline completed successfully.")


if __name__ == "__main__":
    main()
