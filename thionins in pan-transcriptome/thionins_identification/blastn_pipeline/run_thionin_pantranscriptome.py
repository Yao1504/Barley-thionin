#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: run_thionin_pantranscriptome_pipeline.py

Description:
    Identify thionin-related transcripts in pan-transcriptome datasets using
    BLASTN and coordinate validation.

Workflow:
    1. Build BLAST database from transcript FASTA
    2. Run BLASTN using candidate thionin CDS as query
    3. Filter BLASTN hits
    4. Remove isoform suffixes and deduplicate IDs
    5. Extract pan-transcriptome coordinates from GTF
    6. Extract pangenome coordinates from GFF
    7. Compare genomic locations
    8. Retain transcripts with overlapping coordinates and consistent strand
    9. Extract final transcript sequences

Expected input files per genotype:
    <genotype>_cds.fasta
    <genotype>_id.txt
    <genotype>_pg.gff
    <genotype>_pt.gtf
    <genotype>_transcript.fasta

Usage:
    python run_thionin_pantranscriptome_pipeline.py \
        --cds_dir input/cds \
        --id_dir input/pangenome_id \
        --gff_dir input/pangenome_gff \
        --gtf_dir input/pantranscriptome_gtf \
        --transcript_dir input/pantranscriptome_transcript \
        --output_dir output_pantranscriptome \
        --threads 4
"""

import argparse
import glob
import os
import re
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from Bio import SeqIO


BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send",
    "evalue", "bitscore"
]


OUTPUT_SUBDIRS = {
    "blast_results": "1_filter_results",
    "filtered_ids": "2_filter_id",
    "pantranscriptome_coords": "3_filter_pantranscriptome_location",
    "pangenome_coords": "4_filter_pangenome_location",
    "location_comparison": "5_compare_location",
    "final_ids": "6_compare_filter_ids",
    "final_sequences": "7_final_sequences",
    "temp": "temp"
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify thionin transcripts using BLASTN and coordinate validation."
    )

    parser.add_argument("--cds_dir", required=True, help="Directory containing <genotype>_cds.fasta files.")
    parser.add_argument("--id_dir", required=True, help="Directory containing <genotype>_id.txt files.")
    parser.add_argument("--gff_dir", required=True, help="Directory containing <genotype>_pg.gff files.")
    parser.add_argument("--gtf_dir", required=True, help="Directory containing <genotype>_pt.gtf files.")
    parser.add_argument("--transcript_dir", required=True, help="Directory containing <genotype>_transcript.fasta files.")
    parser.add_argument("--output_dir", required=True, help="Output directory.")

    parser.add_argument("--evalue", type=float, default=1e-5, help="BLASTN E-value cutoff. Default: 1e-5.")
    parser.add_argument("--min_length", type=int, default=200, help="Minimum BLAST alignment length. Default: 200.")
    parser.add_argument("--max_mismatch", type=int, default=10, help="Maximum BLAST mismatches. Default: 10.")
    parser.add_argument("--threads", type=int, default=4, help="Number of genotypes processed in parallel. Default: 4.")

    return parser.parse_args()


def make_output_dirs(output_dir):
    paths = {}
    output_dir = Path(output_dir)

    for key, subdir in OUTPUT_SUBDIRS.items():
        path = output_dir / subdir
        path.mkdir(parents=True, exist_ok=True)
        paths[key] = path

    return paths


def get_genotype_name(filename):
    name = Path(filename).name

    if name.endswith("_cds.fasta"):
        return name.replace("_cds.fasta", "")
    if name.endswith("_cds.fa"):
        return name.replace("_cds.fa", "")

    return Path(filename).stem


def remove_isoform_and_deduplicate(ids):
    base_ids = set()

    for value in ids:
        base_id = re.sub(r"\.\d+$", "", str(value))
        if base_id:
            base_ids.add(base_id)

    return sorted(base_ids)


def create_blast_db(transcript_file, db_prefix):
    cmd = [
        "makeblastdb",
        "-in", str(transcript_file),
        "-dbtype", "nucl",
        "-out", str(db_prefix)
    ]

    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run_blastn(cds_file, db_prefix, raw_output, evalue, min_length, max_mismatch):
    cmd = [
        "blastn",
        "-query", str(cds_file),
        "-db", str(db_prefix),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-evalue", str(evalue),
        "-out", str(raw_output)
    ]

    subprocess.run(cmd, check=True)

    if not raw_output.exists() or raw_output.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    df = pd.read_csv(raw_output, sep="\t", header=None, names=BLAST_COLUMNS)

    filtered_df = df[
        (df["length"] > min_length) &
        (df["mismatch"] < max_mismatch)
    ].copy()

    return filtered_df


def parse_gtf_attributes(attributes):
    attr_dict = {}

    for attr in attributes.split(";"):
        attr = attr.strip()
        if attr and " " in attr:
            key, value = attr.split(" ", 1)
            attr_dict[key] = value.strip('"')

    return attr_dict


def parse_gff_attributes(attributes):
    attr_dict = {}

    for attr in attributes.split(";"):
        if "=" in attr:
            key, value = attr.split("=", 1)
            attr_dict[key] = value

    return attr_dict


def extract_gtf_coordinates(gtf_file, gene_ids, output_file):
    coordinates = []

    with open(gtf_file, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = parts

            if feature != "exon":
                continue

            attr_dict = parse_gtf_attributes(attributes)

            gene_id = attr_dict.get("gene_id", "")
            transcript_id = attr_dict.get("transcript_id", "")

            base_gene_id = re.sub(r"\.\d+$", "", gene_id)
            base_transcript_id = re.sub(r"\.\d+$", "", transcript_id)

            target_id = None
            if base_gene_id in gene_ids:
                target_id = base_gene_id
            elif base_transcript_id in gene_ids:
                target_id = base_transcript_id

            if target_id:
                coordinates.append({
                    "gene_id": gene_id,
                    "transcript_id": transcript_id,
                    "chromosome": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "feature": feature,
                    "base_id": target_id
                })

    transcript_coords = {}

    for coord in coordinates:
        transcript_id = coord["transcript_id"]

        if transcript_id not in transcript_coords:
            transcript_coords[transcript_id] = {
                "gene_id": coord["gene_id"],
                "transcript_id": transcript_id,
                "chromosome": coord["chromosome"],
                "strand": coord["strand"],
                "start": coord["start"],
                "end": coord["end"],
                "base_id": coord["base_id"]
            }
        else:
            transcript_coords[transcript_id]["start"] = min(
                transcript_coords[transcript_id]["start"],
                coord["start"]
            )
            transcript_coords[transcript_id]["end"] = max(
                transcript_coords[transcript_id]["end"],
                coord["end"]
            )

    coord_df = pd.DataFrame(list(transcript_coords.values()))

    if not coord_df.empty:
        coord_df.to_csv(output_file, index=False)

    return coord_df


def extract_gff_coordinates(gff_file, gene_ids, output_file):
    coordinates = []

    with open(gff_file, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes = parts

            if feature not in ["gene", "mRNA"]:
                continue

            attr_dict = parse_gff_attributes(attributes)

            gene_id = attr_dict.get("ID", "")
            parent_id = attr_dict.get("Parent", "")

            base_gene_id = re.sub(r"\.\d+$", "", gene_id)
            base_parent_id = re.sub(r"\.\d+$", "", parent_id)

            target_id = None
            if base_gene_id in gene_ids:
                target_id = base_gene_id
            elif base_parent_id in gene_ids:
                target_id = base_parent_id

            if target_id:
                coordinates.append({
                    "gene_id": gene_id,
                    "parent_id": parent_id,
                    "chromosome": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "feature": feature,
                    "base_id": target_id
                })

    coord_df = pd.DataFrame(coordinates)

    if not coord_df.empty:
        coord_df.to_csv(output_file, index=False)

    return coord_df


def compare_locations(pantranscriptome_coords, pangenome_coords, output_file):
    comparisons = []

    for _, pt_row in pantranscriptome_coords.iterrows():
        for _, pg_row in pangenome_coords.iterrows():
            if pt_row["chromosome"] != pg_row["chromosome"]:
                continue

            if pt_row["strand"] != pg_row["strand"]:
                continue

            overlap = max(pt_row["start"], pg_row["start"]) <= min(pt_row["end"], pg_row["end"])

            overlap_length = 0
            if overlap:
                overlap_start = max(pt_row["start"], pg_row["start"])
                overlap_end = min(pt_row["end"], pg_row["end"])
                overlap_length = overlap_end - overlap_start + 1

            comparisons.append({
                "transcriptome_gene_id": pt_row.get("gene_id", ""),
                "transcriptome_transcript_id": pt_row.get("transcript_id", ""),
                "pangenome_gene_id": pg_row.get("gene_id", ""),
                "chromosome": pt_row["chromosome"],
                "strand": pt_row["strand"],
                "transcriptome_start": pt_row["start"],
                "transcriptome_end": pt_row["end"],
                "pangenome_start": pg_row["start"],
                "pangenome_end": pg_row["end"],
                "has_overlap": overlap,
                "overlap_length": overlap_length,
                "strand_consistent": True
            })

    comp_df = pd.DataFrame(comparisons)

    if not comp_df.empty:
        comp_df.to_csv(output_file, index=False)

    return comp_df


def extract_transcript_sequences(transcript_file, gene_ids, output_file):
    sequences = []

    for record in SeqIO.parse(transcript_file, "fasta"):
        base_id = re.sub(r"\.\d+$", "", record.id)

        if base_id in gene_ids:
            sequences.append(record)

    if sequences:
        SeqIO.write(sequences, output_file, "fasta")

    return len(sequences)


def process_genotype(genotype, args_dict):
    cds_dir = Path(args_dict["cds_dir"])
    id_dir = Path(args_dict["id_dir"])
    gff_dir = Path(args_dict["gff_dir"])
    gtf_dir = Path(args_dict["gtf_dir"])
    transcript_dir = Path(args_dict["transcript_dir"])
    output_dirs = {key: Path(value) for key, value in args_dict["output_dirs"].items()}

    evalue = args_dict["evalue"]
    min_length = args_dict["min_length"]
    max_mismatch = args_dict["max_mismatch"]

    print(f"[{genotype}] Processing started")

    cds_file = cds_dir / f"{genotype}_cds.fasta"
    id_file = id_dir / f"{genotype}_id.txt"
    pangenome_gff = gff_dir / f"{genotype}_pg.gff"
    pantranscriptome_gtf = gtf_dir / f"{genotype}_pt.gtf"
    transcript_file = transcript_dir / f"{genotype}_transcript.fasta"

    required_files = [cds_file, id_file, pangenome_gff, pantranscriptome_gtf, transcript_file]
    missing_files = [str(file) for file in required_files if not file.exists()]

    if missing_files:
        return {
            "genotype": genotype,
            "status": "failed",
            "reason": "missing input files",
            "missing_files": ";".join(missing_files)
        }

    blast_db = output_dirs["temp"] / f"{genotype}_blast_db"
    blast_output = output_dirs["blast_results"] / f"{genotype}_blast_results.csv"
    raw_blast_output = output_dirs["blast_results"] / f"{genotype}_blast_raw.tsv"

    try:
        create_blast_db(transcript_file, blast_db)

        blast_results = run_blastn(
            cds_file=cds_file,
            db_prefix=blast_db,
            raw_output=raw_blast_output,
            evalue=evalue,
            min_length=min_length,
            max_mismatch=max_mismatch
        )

        if blast_results.empty:
            return {
                "genotype": genotype,
                "status": "completed",
                "blast_hits": 0,
                "filtered_ids": 0,
                "final_ids": 0,
                "final_sequences": 0,
                "reason": "no BLAST hits"
            }

        blast_results.to_csv(blast_output, index=False)

        transcript_ids = blast_results["sseqid"].unique().tolist()
        filtered_ids = remove_isoform_and_deduplicate(transcript_ids)

        filtered_ids_output = output_dirs["filtered_ids"] / f"{genotype}_filtered_ids.txt"
        with open(filtered_ids_output, "w", encoding="utf-8") as handle:
            for gene_id in filtered_ids:
                handle.write(f"{gene_id}\n")

        pt_coords_output = output_dirs["pantranscriptome_coords"] / f"{genotype}_pantranscriptome_coords.csv"
        pt_coords = extract_gtf_coordinates(
            gtf_file=pantranscriptome_gtf,
            gene_ids=set(filtered_ids),
            output_file=pt_coords_output
        )

        with open(id_file, "r", encoding="utf-8") as handle:
            pangenome_ids = [line.strip() for line in handle if line.strip()]

        pangenome_base_ids = remove_isoform_and_deduplicate(pangenome_ids)

        pg_coords_output = output_dirs["pangenome_coords"] / f"{genotype}_pangenome_coords.csv"
        pg_coords = extract_gff_coordinates(
            gff_file=pangenome_gff,
            gene_ids=set(pangenome_base_ids),
            output_file=pg_coords_output
        )

        if pt_coords.empty or pg_coords.empty:
            return {
                "genotype": genotype,
                "status": "completed",
                "blast_hits": len(blast_results),
                "filtered_ids": len(filtered_ids),
                "final_ids": 0,
                "final_sequences": 0,
                "reason": "missing coordinate information"
            }

        comparison_output = output_dirs["location_comparison"] / f"{genotype}_location_comparison.csv"
        comparison_results = compare_locations(
            pantranscriptome_coords=pt_coords,
            pangenome_coords=pg_coords,
            output_file=comparison_output
        )

        if comparison_results.empty:
            return {
                "genotype": genotype,
                "status": "completed",
                "blast_hits": len(blast_results),
                "filtered_ids": len(filtered_ids),
                "final_ids": 0,
                "final_sequences": 0,
                "reason": "no overlapping coordinates"
            }

        filtered_comparisons = comparison_results[
            (comparison_results["has_overlap"] == True) &
            (comparison_results["strand_consistent"] == True)
        ]

        if filtered_comparisons.empty:
            return {
                "genotype": genotype,
                "status": "completed",
                "blast_hits": len(blast_results),
                "filtered_ids": len(filtered_ids),
                "final_ids": 0,
                "final_sequences": 0,
                "reason": "no overlap with consistent strand"
            }

        valid_transcript_ids = remove_isoform_and_deduplicate(
            filtered_comparisons["transcriptome_gene_id"].tolist()
        )

        final_ids_output = output_dirs["final_ids"] / f"{genotype}_final_filtered_ids.txt"
        with open(final_ids_output, "w", encoding="utf-8") as handle:
            for gene_id in valid_transcript_ids:
                handle.write(f"{gene_id}\n")

        final_sequences_output = output_dirs["final_sequences"] / f"{genotype}_final_sequences.fasta"
        sequence_count = extract_transcript_sequences(
            transcript_file=transcript_file,
            gene_ids=set(valid_transcript_ids),
            output_file=final_sequences_output
        )

        return {
            "genotype": genotype,
            "status": "completed",
            "blast_hits": len(blast_results),
            "filtered_ids": len(filtered_ids),
            "final_ids": len(valid_transcript_ids),
            "final_sequences": sequence_count,
            "reason": "success"
        }

    except Exception as error:
        return {
            "genotype": genotype,
            "status": "failed",
            "reason": str(error)
        }

    finally:
        for temp_file in glob.glob(str(blast_db) + ".*"):
            try:
                os.remove(temp_file)
            except OSError:
                pass

        if raw_blast_output.exists():
            try:
                raw_blast_output.unlink()
            except OSError:
                pass


def main():
    args = parse_args()
    start_time = time.time()

    required_dirs = [
        args.cds_dir,
        args.id_dir,
        args.gff_dir,
        args.gtf_dir,
        args.transcript_dir
    ]

    for directory in required_dirs:
        if not Path(directory).exists():
            raise FileNotFoundError(f"Input directory not found: {directory}")

    output_dirs = make_output_dirs(args.output_dir)

    genotypes = []
    for file in Path(args.cds_dir).iterdir():
        if file.name.endswith("_cds.fasta") or file.name.endswith("_cds.fa"):
            genotypes.append(get_genotype_name(file.name))

    genotypes = sorted(genotypes)

    if not genotypes:
        raise ValueError(f"No *_cds.fasta or *_cds.fa files found in {args.cds_dir}")

    print(f"Detected {len(genotypes)} genotypes.")
    print(", ".join(genotypes))

    args_dict = {
        "cds_dir": args.cds_dir,
        "id_dir": args.id_dir,
        "gff_dir": args.gff_dir,
        "gtf_dir": args.gtf_dir,
        "transcript_dir": args.transcript_dir,
        "output_dirs": {key: str(value) for key, value in output_dirs.items()},
        "evalue": args.evalue,
        "min_length": args.min_length,
        "max_mismatch": args.max_mismatch
    }

    max_workers = min(args.threads, len(genotypes))
    results = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_genotype, genotype, args_dict): genotype
            for genotype in genotypes
        }

        for future in as_completed(futures):
            genotype = futures[future]
            try:
                result = future.result()
                results.append(result)
                print(f"[{genotype}] {result.get('status')} - {result.get('reason', '')}")
            except Exception as error:
                results.append({
                    "genotype": genotype,
                    "status": "failed",
                    "reason": str(error)
                })
                print(f"[{genotype}] failed - {error}")

    summary_df = pd.DataFrame(results)
    summary_file = Path(args.output_dir) / "pantranscriptome_pipeline_summary.csv"
    summary_df.to_csv(summary_file, index=False)

    elapsed = time.time() - start_time

    print("\nPipeline completed.")
    print(f"Processed genotypes: {len(genotypes)}")
    print(f"Summary file: {summary_file}")
    print(f"Elapsed time: {elapsed:.2f} seconds")


if __name__ == "__main__":
    main()
