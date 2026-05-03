#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: run_alignment_trimming_pipeline.py

Description:
    Build a codon alignment from multiple CDS FASTA files.

    This pipeline:
        1. Merges CDS FASTA files
        2. Translates CDS sequences into protein sequences
        3. Aligns protein sequences using MAFFT
        4. Converts protein alignment back to codon alignment using PAL2NAL
        5. Trims the codon alignment using trimAl

Input:
    - Directory containing CDS FASTA files

Output:
    - all_cds.fasta
    - all_protein.fasta
    - all_protein_aln.fasta
    - all_cds_codon_aln.fasta
    - all_cds_trimmed.fasta

Usage:
    python run_codon_alignment_pipeline.py \
        --input_dir input_cds \
        --output_dir output_tree \
        --mafft mafft \
        --trimal trimal \
        --perl perl \
        --pal2nal pal2nal.pl
"""

import argparse
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate trimmed codon alignment from CDS FASTA files."
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing CDS FASTA files."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to save output files."
    )

    parser.add_argument(
        "--mafft",
        default="mafft",
        help="Path to MAFFT executable. Default: mafft"
    )

    parser.add_argument(
        "--trimal",
        default="trimal",
        help="Path to trimAl executable. Default: trimal"
    )

    parser.add_argument(
        "--perl",
        default="perl",
        help="Path to Perl executable. Default: perl"
    )

    parser.add_argument(
        "--pal2nal",
        required=True,
        help="Path to pal2nal.pl script."
    )

    return parser.parse_args()


def run_command(command, output_file=None):
    if output_file:
        with open(output_file, "w") as out_f:
            subprocess.run(
                command,
                stdout=out_f,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
    else:
        subprocess.run(
            command,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )


def main():
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    if not input_dir.exists():
        sys.exit(f"Error: input directory not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    all_cds = output_dir / "all_cds.fasta"
    all_protein = output_dir / "all_protein.fasta"
    protein_aln = output_dir / "all_protein_aln.fasta"
    cds_codon_aln = output_dir / "all_cds_codon_aln.fasta"
    cds_trimmed = output_dir / "all_cds_trimmed.fasta"

    print("==== Step 1: Merging CDS files ====")

    records = []
    file_count = 0

    for file in sorted(input_dir.iterdir()):
        if file.suffix.lower() in [".fa", ".fasta"]:
            file_count += 1
            print(f"  · Reading file: {file.name}")

            for rec in SeqIO.parse(file, "fasta"):
                records.append(rec)

    print(f"Total files: {file_count}")
    print(f"Total sequences before ID correction: {len(records)}\n")

    if not records:
        sys.exit("Error: no CDS sequences found.")

    seen_ids = {}

    for rec in records:
        if rec.id in seen_ids:
            seen_ids[rec.id] += 1
            rec.id = f"{rec.id}_dup{seen_ids[rec.id]}"
            rec.name = rec.id
            rec.description = rec.id
        else:
            seen_ids[rec.id] = 1

    SeqIO.write(records, all_cds, "fasta")
    print(f"Merged CDS written to: {all_cds}\n")

    print("==== Step 2: Translating CDS to proteins ====")

    prot_records = []

    for rec in records:
        prot = rec.seq.translate(to_stop=False)
        prot_records.append(
            SeqRecord(
                prot,
                id=rec.id,
                name=rec.id,
                description=""
            )
        )

    SeqIO.write(prot_records, all_protein, "fasta")
    print(f"Protein sequences written to: {all_protein}\n")

    print("==== Step 3: Running MAFFT ====")

    cmd_mafft = [
        args.mafft,
        "--auto",
        str(all_protein)
    ]

    try:
        run_command(cmd_mafft, output_file=protein_aln)
    except subprocess.CalledProcessError as error:
        sys.exit(f"Error: MAFFT failed.\n{error.stderr}")

    print(f"Protein alignment written to: {protein_aln}\n")

    print("==== Step 4: Reordering CDS based on protein alignment IDs ====")

    aligned_proteins = list(SeqIO.parse(protein_aln, "fasta"))
    cds_dict = {
        rec.id: rec
        for rec in SeqIO.parse(all_cds, "fasta")
    }

    reordered_cds = []

    for p in aligned_proteins:
        if p.id in cds_dict:
            reordered_cds.append(cds_dict[p.id])
        else:
            print(f"Warning: no matching CDS found for: {p.id}")

    SeqIO.write(reordered_cds, all_cds, "fasta")
    print("CDS sequences reordered based on protein alignment IDs.\n")

    print("==== Step 5: Running PAL2NAL ====")

    cmd_pal2nal = [
        args.perl,
        args.pal2nal,
        str(protein_aln),
        str(all_cds),
        "-output",
        "fasta"
    ]

    try:
        run_command(cmd_pal2nal, output_file=cds_codon_aln)
    except subprocess.CalledProcessError as error:
        sys.exit(f"Error: PAL2NAL failed.\n{error.stderr}")

    print(f"Codon alignment written to: {cds_codon_aln}\n")

    print("==== Step 6: Running trimAl ====")

    cmd_trimal = [
        args.trimal,
        "-in",
        str(cds_codon_aln),
        "-out",
        str(cds_trimmed),
        "-automated1"
    ]

    try:
        run_command(cmd_trimal)
    except subprocess.CalledProcessError as error:
        sys.exit(f"Error: trimAl failed.\n{error.stderr}")

    print(f"Trimmed codon alignment written to: {cds_trimmed}\n")

    print("======== Pipeline completed successfully ========")
    print(f"Input CDS sequences: {len(records)}")
    print(f"Trimmed codon alignment: {cds_trimmed}")
    print("===============================================")


if __name__ == "__main__":
    main()
