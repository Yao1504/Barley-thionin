#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: map_pantranscriptome_to_pangenome_clusters.py

Description:
    Map pan-transcriptome transcripts to pangenome CDS clusters using BLASTN,
    then generate cluster-level and subcluster-level summary files.

Workflow:
    1. Load pangenome CDS cluster definitions
    2. Merge pangenome CDS FASTA files
    3. Clean CDS and transcript FASTA files by removing gap characters
    4. Build BLAST database from cleaned pangenome CDS
    5. Run BLASTN using pan-transcriptome sequences as query
    6. Filter BLAST hits by identity and coverage
    7. Assign each transcript to the best CDS hit and corresponding CDS cluster
    8. Export transcript lists for each cluster
    9. Build sequence-identical subclusters within each species-cluster group
    10. Generate cluster and subcluster summary Excel files

Expected input files:
    cds_dir:
        *.fa / *.fasta / *.fna

    cluster_dir:
        cluster*.txt
        singleton.txt  optional

    transcript_fasta:
        pan-transcriptome transcript FASTA file

Usage:
    python map_pantranscriptome_to_pangenome_clusters.py \
        --cds_dir input/all_pangenome_cds \
        --cluster_dir input/cluster_pangenome \
        --transcript_fasta input/all_transcripts.fasta \
        --output_dir output \
        --threads 8
"""

import argparse
import glob
import os
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "qlen", "slen", "bitscore"
]


OUTPUT_SUBDIRS = {
    "blastdb": "blastdb",
    "cluster_transcripts": "cluster_transcripts",
    "subcluster_transcripts": "subcluster_transcripts",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Map pan-transcriptome transcripts to pangenome CDS clusters using BLASTN."
    )

    parser.add_argument(
        "--cds_dir",
        required=True,
        help="Directory containing pangenome CDS FASTA files."
    )

    parser.add_argument(
        "--cluster_dir",
        required=True,
        help="Directory containing cluster*.txt files and optional singleton.txt."
    )

    parser.add_argument(
        "--transcript_fasta",
        required=True,
        help="Pan-transcriptome transcript FASTA file."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory."
    )

    parser.add_argument(
        "--min_pident",
        type=float,
        default=90.0,
        help="Minimum BLAST percent identity. Default: 90.0"
    )

    parser.add_argument(
        "--min_qcov",
        type=float,
        default=0.7,
        help="Minimum query coverage for full-length transcripts. Default: 0.7"
    )

    parser.add_argument(
        "--min_scov",
        type=float,
        default=0.7,
        help="Minimum subject coverage for full-length transcripts. Default: 0.7"
    )

    parser.add_argument(
        "--partial_length_cutoff",
        type=int,
        default=300,
        help="Transcripts shorter than this value are treated as partial sequences. Default: 300"
    )

    parser.add_argument(
        "--max_subcluster_cols",
        type=int,
        default=5,
        help="Maximum number of subcluster columns in summary.xlsx. Default: 5"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of BLAST threads. Default: 8"
    )

    parser.add_argument(
        "--max_target_seqs",
        type=int,
        default=20,
        help="Maximum number of target sequences reported by BLASTN. Default: 20"
    )

    parser.add_argument(
        "--evalue",
        default="1e-20",
        help="BLASTN e-value threshold. Default: 1e-20"
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing intermediate files."
    )

    return parser.parse_args()


def load_pangenome_clusters(cluster_dir):
    cds_to_cluster = {}
    cluster_to_cds = defaultdict(set)

    for path in glob.glob(os.path.join(cluster_dir, "cluster*.txt")):
        cluster_name = os.path.splitext(os.path.basename(path))[0]

        with open(path, encoding="utf-8") as f:
            for line in f:
                cds_id = line.strip()

                if cds_id:
                    cds_to_cluster[cds_id] = cluster_name
                    cluster_to_cds[cluster_name].add(cds_id)

    singleton_path = os.path.join(cluster_dir, "singleton.txt")

    if os.path.exists(singleton_path):
        with open(singleton_path, encoding="utf-8") as f:
            for line in f:
                cds_id = line.strip()

                if cds_id:
                    cds_to_cluster[cds_id] = cds_id
                    cluster_to_cds[cds_id].add(cds_id)

    print(f"[cluster] Loaded {len(cluster_to_cds)} clusters, including singletons.")
    return cds_to_cluster, cluster_to_cds


def merge_fasta(input_dir, output_fasta, overwrite=False):
    if os.path.exists(output_fasta) and not overwrite:
        print(f"[merge_fasta] Skip existing file: {output_fasta}")
        return

    fasta_files = (
        glob.glob(os.path.join(input_dir, "*.fa"))
        + glob.glob(os.path.join(input_dir, "*.fasta"))
        + glob.glob(os.path.join(input_dir, "*.fna"))
    )

    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in: {input_dir}")

    with open(output_fasta, "w", encoding="utf-8") as out:
        for fasta_file in fasta_files:
            with open(fasta_file, encoding="utf-8") as fin:
                out.write(fin.read())

            out.write("\n")

    print(f"[merge_fasta] Written: {output_fasta}")


def clean_fasta(input_fasta, output_fasta, overwrite=False):
    if os.path.exists(output_fasta) and not overwrite:
        print(f"[clean_fasta] Skip existing file: {output_fasta}")
        return

    with open(input_fasta, encoding="utf-8") as fin, \
            open(output_fasta, "w", encoding="utf-8") as fout:

        for line in fin:
            if line.startswith(">"):
                fout.write(line)
            else:
                fout.write(line.replace("-", ""))

    print(f"[clean_fasta] Written: {output_fasta}")


def get_all_transcript_ids(fasta_path):
    transcript_ids = []

    with open(fasta_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                transcript_ids.append(line[1:].strip().split()[0])

    return transcript_ids


def read_fasta_dict(fasta_path):
    seqs = {}
    header = None
    seq_lines = []

    with open(fasta_path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()

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


def get_species_from_id(transcript_id):
    if "_chr" in transcript_id:
        return transcript_id.split("_chr", 1)[0]

    return transcript_id.split("_", 1)[0]


def make_blastdb(fasta, db_prefix, overwrite=False):
    if os.path.exists(db_prefix + ".nin") and not overwrite:
        print(f"[makeblastdb] Skip existing database: {db_prefix}")
        return

    os.makedirs(os.path.dirname(db_prefix), exist_ok=True)

    cmd = [
        "makeblastdb",
        "-in", fasta,
        "-dbtype", "nucl",
        "-out", db_prefix,
    ]

    print("[makeblastdb] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_blast(query, db_prefix, output_file, threads, max_target_seqs, evalue, overwrite=False):
    if os.path.exists(output_file) and not overwrite:
        print(f"[blastn] Skip existing file: {output_file}")
        return

    cmd = [
        "blastn",
        "-query", query,
        "-db", db_prefix,
        "-outfmt", "6 qseqid sseqid pident length qlen slen bitscore",
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-num_threads", str(threads),
        "-out", output_file,
    ]

    print("[blastn] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_blast_and_pick_best(
    blast_file,
    min_pident,
    min_qcov,
    min_scov,
    partial_length_cutoff
):
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        return pd.DataFrame(
            columns=["transcript_id", "best_cds_id", "pident", "qcov", "scov"]
        )

    df = pd.read_csv(
        blast_file,
        sep="\t",
        header=None,
        names=BLAST_COLUMNS
    )

    if df.empty:
        return pd.DataFrame(
            columns=["transcript_id", "best_cds_id", "pident", "qcov", "scov"]
        )

    df["qcov"] = df["length"] / df["qlen"]
    df["scov"] = df["length"] / df["slen"]
    df["is_partial"] = df["qlen"] < partial_length_cutoff

    cond_full = (
        (~df["is_partial"])
        & (df["pident"] >= min_pident)
        & (df["qcov"] >= min_qcov)
        & (df["scov"] >= min_scov)
    )

    cond_partial = (
        (df["is_partial"])
        & (df["pident"] >= min_pident)
        & (df["qcov"] >= 0.9)
    )

    filtered = df[cond_full | cond_partial].copy()

    if filtered.empty:
        return pd.DataFrame(
            columns=["transcript_id", "best_cds_id", "pident", "qcov", "scov"]
        )

    filtered = filtered.sort_values(
        ["qseqid", "bitscore"],
        ascending=[True, False]
    )

    best = filtered.groupby("qseqid", as_index=False).head(1)

    best = best.rename(
        columns={
            "qseqid": "transcript_id",
            "sseqid": "best_cds_id",
        }
    )

    return best[["transcript_id", "best_cds_id", "pident", "qcov", "scov"]]


def build_mapping_dataframe(best_hits, all_transcript_ids, cds_to_cluster):
    all_transcript_set = set(all_transcript_ids)

    if best_hits.empty:
        mapping_df = pd.DataFrame(
            {
                "transcript_id": all_transcript_ids,
                "best_cds_id": [pd.NA] * len(all_transcript_ids),
                "pident": [pd.NA] * len(all_transcript_ids),
                "qcov": [pd.NA] * len(all_transcript_ids),
                "scov": [pd.NA] * len(all_transcript_ids),
                "cluster": ["unassigned"] * len(all_transcript_ids),
            }
        )

    else:
        mapping_df = best_hits.copy()

        mapping_df["cluster"] = mapping_df["best_cds_id"].map(
            lambda cds_id: cds_to_cluster.get(cds_id, "unassigned")
        )

        missed = list(all_transcript_set - set(mapping_df["transcript_id"]))

        if missed:
            extra = pd.DataFrame(
                {
                    "transcript_id": missed,
                    "best_cds_id": [pd.NA] * len(missed),
                    "pident": [pd.NA] * len(missed),
                    "qcov": [pd.NA] * len(missed),
                    "scov": [pd.NA] * len(missed),
                    "cluster": ["unassigned"] * len(missed),
                }
            )

            mapping_df = pd.concat([mapping_df, extra], ignore_index=True)

    mapping_df["species"] = mapping_df["transcript_id"].map(get_species_from_id)

    return mapping_df


def write_cluster_transcript_files(mapping_df, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for cluster, transcript_ids in mapping_df.groupby("cluster")["transcript_id"]:
        output_file = os.path.join(output_dir, f"{cluster}.txt")

        with open(output_file, "w", encoding="utf-8") as f:
            f.write("\n".join(transcript_ids))

    print(f"[cluster_transcripts] Written to: {output_dir}")


def build_subclusters(mapping_df, transcript_seqs, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    df = mapping_df.copy()
    df["species"] = df["transcript_id"].map(get_species_from_id)

    valid = df[df["cluster"] != "unassigned"].copy()

    transcript_to_subcluster = {}
    subcluster_rows = []
    group_index = Counter()

    for (species, cluster), subdf in valid.groupby(["species", "cluster"]):
        transcript_ids = list(subdf["transcript_id"])

        if not transcript_ids:
            continue

        seq_to_ids = defaultdict(list)

        for transcript_id in transcript_ids:
            seq = transcript_seqs.get(transcript_id, "")
            seq_to_ids[seq].append(transcript_id)

        for seq, members in seq_to_ids.items():
            if len(members) < 2:
                continue

            group_index[(species, cluster)] += 1
            idx = group_index[(species, cluster)]
            subcluster_name = f"subcluster_{idx}"

            output_file = os.path.join(
                output_dir,
                f"{species}_{cluster}_{subcluster_name}.txt"
            )

            with open(output_file, "w", encoding="utf-8") as f:
                f.write("\n".join(members))

            for transcript_id in members:
                transcript_to_subcluster[transcript_id] = subcluster_name

            subcluster_rows.append(
                {
                    "species": species,
                    "cluster": cluster,
                    "subcluster": subcluster_name,
                    "n_members": len(members),
                }
            )

    df["subcluster"] = df["transcript_id"].map(
        lambda transcript_id: transcript_to_subcluster.get(
            transcript_id,
            "unassigned"
        )
    )

    return df, pd.DataFrame(subcluster_rows)


def build_summary_excel(mapping_with_subclusters, output_xlsx, max_subcluster_cols):
    df = mapping_with_subclusters.copy()
    df["species"] = df["species"].astype(str)

    pattern = re.compile(r"^cluster_(\d+)$")

    df["cluster_num"] = df["cluster"].map(
        lambda cluster: int(pattern.match(cluster).group(1))
        if isinstance(cluster, str) and pattern.match(cluster)
        else np.nan
    )

    df_valid = df[~df["cluster_num"].isna()].copy()
    df_valid["cluster_num"] = df_valid["cluster_num"].astype(int)

    species_list = sorted(df_valid["species"].unique())
    cluster_nums = sorted(df_valid["cluster_num"].unique())

    cluster_summary_rows = []

    for species in species_list:
        row = {"Species": species}
        sub = df_valid[df_valid["species"] == species]

        for cluster_num in cluster_nums:
            row[f"Cluster_{cluster_num}"] = sub[
                sub["cluster_num"] == cluster_num
            ]["transcript_id"].nunique()

        cluster_summary_rows.append(row)

    cluster_summary = pd.DataFrame(cluster_summary_rows)

    subcluster_summary_rows = []

    for species in species_list:
        for cluster_num in cluster_nums:
            label = f"{species}_cluster{cluster_num}"

            sub = df_valid[
                (df_valid["species"] == species)
                & (df_valid["cluster_num"] == cluster_num)
            ]

            if sub.empty:
                row = {
                    "Species_Cluster": label,
                    **{
                        f"subcluster_{idx}": 0
                        for idx in range(1, max_subcluster_cols + 1)
                    },
                }

                subcluster_summary_rows.append(row)
                continue

            counts = (
                sub[sub["subcluster"] != "unassigned"]
                .groupby("subcluster")["transcript_id"]
                .nunique()
                .sort_values(ascending=False)
            )

            sizes = list(counts)

            row = {"Species_Cluster": label}

            for idx in range(1, max_subcluster_cols + 1):
                row[f"subcluster_{idx}"] = (
                    sizes[idx - 1] if idx <= len(sizes) else 0
                )

            subcluster_summary_rows.append(row)

    subcluster_summary = pd.DataFrame(subcluster_summary_rows)

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        cluster_summary.to_excel(
            writer,
            sheet_name="cluster_summary",
            index=False
        )

        subcluster_summary.to_excel(
            writer,
            sheet_name="subcluster_summary",
            index=False
        )

    print(f"[summary] Written: {output_xlsx}")


def prepare_output_paths(output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    subdirs = {}

    for key, value in OUTPUT_SUBDIRS.items():
        path = output_dir / value
        path.mkdir(parents=True, exist_ok=True)
        subdirs[key] = path

    paths = {
        "output_dir": output_dir,
        "merged_cds": output_dir / "all_pangenome_cds_merged.fasta",
        "merged_cds_clean": output_dir / "all_pangenome_cds_clean.fasta",
        "transcripts_clean": output_dir / "all_transcripts_clean.fasta",
        "blast_db_prefix": subdirs["blastdb"] / "pangenome_cds_db",
        "blast_output": output_dir / "blast_transcripts_vs_cds.tsv",
        "mapping_tsv": output_dir / "transcript_to_cds_cluster.tsv",
        "mapping_xlsx": output_dir / "transcript_to_cds_cluster.xlsx",
        "subcluster_stats_xlsx": output_dir / "subcluster_stats.xlsx",
        "transcript_cluster_subcluster_xlsx": output_dir / "transcript_cluster_subcluster.xlsx",
        "summary_xlsx": output_dir / "summary.xlsx",
        "cluster_transcripts_dir": subdirs["cluster_transcripts"],
        "subcluster_transcripts_dir": subdirs["subcluster_transcripts"],
    }

    return paths


def main():
    args = parse_args()
    paths = prepare_output_paths(args.output_dir)

    print("\n=== Step 1: Load pangenome clusters ===")
    cds_to_cluster, _ = load_pangenome_clusters(args.cluster_dir)

    print("\n=== Step 2: Merge and clean FASTA files ===")
    merge_fasta(
        args.cds_dir,
        str(paths["merged_cds"]),
        overwrite=args.overwrite
    )

    clean_fasta(
        str(paths["merged_cds"]),
        str(paths["merged_cds_clean"]),
        overwrite=args.overwrite
    )

    clean_fasta(
        args.transcript_fasta,
        str(paths["transcripts_clean"]),
        overwrite=args.overwrite
    )

    print("\n=== Step 3: Build BLAST database and run BLASTN ===")
    make_blastdb(
        str(paths["merged_cds_clean"]),
        str(paths["blast_db_prefix"]),
        overwrite=args.overwrite
    )

    run_blast(
        query=str(paths["transcripts_clean"]),
        db_prefix=str(paths["blast_db_prefix"]),
        output_file=str(paths["blast_output"]),
        threads=args.threads,
        max_target_seqs=args.max_target_seqs,
        evalue=args.evalue,
        overwrite=args.overwrite
    )

    print("\n=== Step 4: Parse BLAST results and assign clusters ===")
    best_hits = parse_blast_and_pick_best(
        blast_file=str(paths["blast_output"]),
        min_pident=args.min_pident,
        min_qcov=args.min_qcov,
        min_scov=args.min_scov,
        partial_length_cutoff=args.partial_length_cutoff
    )

    all_transcript_ids = get_all_transcript_ids(str(paths["transcripts_clean"]))

    mapping_df = build_mapping_dataframe(
        best_hits=best_hits,
        all_transcript_ids=all_transcript_ids,
        cds_to_cluster=cds_to_cluster
    )

    mapping_df.to_csv(paths["mapping_tsv"], sep="\t", index=False)
    mapping_df.to_excel(paths["mapping_xlsx"], index=False)

    print(f"[mapping] Written: {paths['mapping_tsv']}")
    print(f"[mapping] Written: {paths['mapping_xlsx']}")

    write_cluster_transcript_files(
        mapping_df,
        str(paths["cluster_transcripts_dir"])
    )

    print("\n=== Step 5: Build subclusters ===")
    transcript_seqs = read_fasta_dict(str(paths["transcripts_clean"]))

    mapping_with_subclusters, subcluster_stats = build_subclusters(
        mapping_df=mapping_df,
        transcript_seqs=transcript_seqs,
        output_dir=str(paths["subcluster_transcripts_dir"])
    )

    subcluster_stats.to_excel(paths["subcluster_stats_xlsx"], index=False)
    mapping_with_subclusters.to_excel(
        paths["transcript_cluster_subcluster_xlsx"],
        index=False
    )

    print(f"[subcluster_stats] Written: {paths['subcluster_stats_xlsx']}")
    print(
        "[transcript_cluster_subcluster] Written: "
        f"{paths['transcript_cluster_subcluster_xlsx']}"
    )

    print("\n=== Step 6: Generate summary.xlsx ===")
    build_summary_excel(
        mapping_with_subclusters=mapping_with_subclusters,
        output_xlsx=str(paths["summary_xlsx"]),
        max_subcluster_cols=args.max_subcluster_cols
    )

    print("\n=== ALL DONE ===\n")


if __name__ == "__main__":
    main()
