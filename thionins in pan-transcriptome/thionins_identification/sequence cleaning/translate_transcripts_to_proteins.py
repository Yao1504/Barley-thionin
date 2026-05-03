#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: translate_transcripts_to_proteins.py

Description:
    Translate transcript/CDS FASTA files into protein sequences.

    For each sequence:
        1. Remove gap characters
        2. Search for a strict ORF from ATG to an in-frame stop codon
        3. If a strict ORF is found:
            - translate the ORF
            - keep the original sequence ID
        4. If no strict ORF is found:
            - label the sequence as pseudogene
            - translate the full cleaned sequence
        5. Write protein FASTA, cleaned CDS FASTA, and ORF report

Input:
    - Directory containing FASTA files

Output:
    - Protein FASTA files
    - Cleaned CDS FASTA files
    - ORF report files

Usage:
    python translate_transcripts_to_proteins.py \
        --input_dir input_fasta \
        --protein_dir output_proteins \
        --cds_dir output_cds
"""

import argparse
import re
from pathlib import Path


STOP_CODONS = {"TAA", "TAG", "TGA"}

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
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Translate transcript/CDS FASTA files into protein sequences."
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing input FASTA files."
    )

    parser.add_argument(
        "--protein_dir",
        required=True,
        help="Directory to save translated protein FASTA files."
    )

    parser.add_argument(
        "--cds_dir",
        required=True,
        help="Directory to save cleaned CDS FASTA files."
    )

    parser.add_argument(
        "--input_suffix",
        default=".fasta",
        help="Input FASTA suffix. Default: .fasta"
    )

    return parser.parse_args()


def clean_header(header):
    """
    Remove sequence coordinate suffix from FASTA header.

    Example:
        GeneA.1/47-457 -> GeneA.1
    """
    return re.sub(r"/\d+-\d+", "", header)


def read_fasta(path):
    """
    Read FASTA file into a dictionary.
    """
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
                    seqs[clean_header(header)] = "".join(seq_lines)

                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line)

        if header:
            seqs[clean_header(header)] = "".join(seq_lines)

    return seqs


def find_strict_orf(seq):
    """
    Find the first strict ORF from ATG to an in-frame stop codon.

    Returns:
        tuple: (start, end, nucleotide_orf)
        None: if no strict ORF is detected
    """
    seq = seq.upper().replace("-", "")

    for frame in range(3):
        for i in range(frame, len(seq) - 2, 3):
            if seq[i:i + 3] == "ATG":
                j = i + 3

                while j < len(seq) - 2:
                    if seq[j:j + 3] in STOP_CODONS:
                        return i, j + 3, seq[i:j + 3]

                    j += 3

                return None

    return None


def translate_dna(seq):
    """
    Translate DNA sequence into amino acid sequence.
    Unknown codons are translated as X.
    """
    seq = seq.upper().replace("-", "")
    amino_acids = []

    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        amino_acids.append(CODON_TABLE.get(codon, "X"))

    return "".join(amino_acids)


def process_fasta_file(input_file, protein_dir, cds_dir):
    """
    Process one FASTA file and generate protein, CDS, and ORF report outputs.
    """
    protein_dir = Path(protein_dir)
    cds_dir = Path(cds_dir)

    output_protein = protein_dir / input_file.name.replace(".fasta", "_protein.fasta")
    output_report = protein_dir / input_file.name.replace(".fasta", "_orf_report.txt")
    output_cds = cds_dir / input_file.name.replace(".fasta", "_cds_clean.fasta")

    print(f"Processing: {input_file.name}")

    seqs = read_fasta(input_file)

    protein_output = []
    cds_output = []
    report_output = ["ID\tORF_info\tStatus"]

    ok_count = 0
    pseudogene_count = 0

    for header, seq in seqs.items():
        seq_clean = seq.upper().replace("-", "")
        orf_info = find_strict_orf(seq_clean)

        if orf_info is None:
            new_id = f"{header}_pseudogene"

            cds_output.append(f">{new_id}\n{seq_clean}\n")

            aa = translate_dna(seq_clean)
            aa_clean = aa.replace("X", "").replace("-", "")

            protein_output.append(f">{new_id}\n{aa_clean}\n")
            report_output.append(f"{header}\tNO_ORF\tPSEUDOGENE")

            pseudogene_count += 1

        else:
            start, end, nuc = orf_info

            cds_output.append(f">{header}\n{seq_clean}\n")

            aa = translate_dna(nuc)
            aa_clean = aa.replace("-", "")

            protein_output.append(f">{header}\n{aa_clean}\n")
            report_output.append(f"{header}\t{start}-{end}\tOK")

            ok_count += 1

    with open(output_protein, "w", encoding="utf-8") as handle:
        handle.writelines(protein_output)

    with open(output_report, "w", encoding="utf-8") as handle:
        handle.write("\n".join(report_output))

    with open(output_cds, "w", encoding="utf-8") as handle:
        handle.writelines(cds_output)

    print(f"  Protein output: {output_protein}")
    print(f"  CDS output: {output_cds}")
    print(f"  ORF report: {output_report}")
    print(f"  OK sequences: {ok_count}")
    print(f"  Pseudogene-like sequences: {pseudogene_count}")

    return ok_count, pseudogene_count


def main():
    args = parse_args()

    input_dir = Path(args.input_dir)
    protein_dir = Path(args.protein_dir)
    cds_dir = Path(args.cds_dir)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    protein_dir.mkdir(parents=True, exist_ok=True)
    cds_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(input_dir.glob(f"*{args.input_suffix}"))

    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files ending with '{args.input_suffix}' found in {input_dir}"
        )

    total_ok = 0
    total_pseudogene = 0

    print("=" * 70)
    print("Transcript-to-protein translation pipeline")
    print("=" * 70)

    for fasta_file in fasta_files:
        ok_count, pseudogene_count = process_fasta_file(
            input_file=fasta_file,
            protein_dir=protein_dir,
            cds_dir=cds_dir
        )

        total_ok += ok_count
        total_pseudogene += pseudogene_count

    print("\nSummary")
    print("-" * 70)
    print(f"Processed FASTA files: {len(fasta_files)}")
    print(f"Sequences with strict ORF: {total_ok}")
    print(f"Pseudogene-like sequences without strict ORF: {total_pseudogene}")
    print("Pipeline completed successfully.")


if __name__ == "__main__":
    main()
