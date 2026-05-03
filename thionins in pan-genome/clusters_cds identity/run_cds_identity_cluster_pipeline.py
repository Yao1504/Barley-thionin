#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: run_cds_identity_cluster_pipeline.py

Description:
    Generate codon alignment from CDS sequences, calculate pairwise identity,
    cluster sequences by identity threshold, and summarize cluster/subcluster
    distribution across genotypes.

Workflow:
    1. Merge CDS FASTA files
    2. Translate CDS to protein
    3. Align protein sequences using MAFFT
    4. Generate codon alignment using PAL2NAL
    5. Calculate pairwise identity matrix
    6. Cluster sequences by identity threshold
    7. Build 100% identity subclusters
    8. Summarize cluster and subcluster counts by genotype

Usage:
    python run_cds_identity_cluster_pipeline.py \
        --input_dir input_cds \
        --output_dir output \
        --perl perl \
        --pal2nal pal2nal.pl \
        --cluster_threshold 90

Requirements:
    - Python >= 3.8
    - Biopython
    - pandas
    - numpy
    - MAFFT
    - Perl
    - PAL2NAL
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="CDS identity clustering and subclustering pipeline."
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing CDS FASTA files."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory for output files."
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

    parser.add_argument(
        "--cluster_threshold",
        type=float,
        default=90,
        help="Identity threshold for clustering. Default: 90"
    )

    parser.add_argument(
        "--genotype_field",
        type=int,
        default=1,
        help=(
            "Zero-based field index used to extract genotype name from "
            "dot-separated sequence ID. Default: 1"
        )
    )

    return parser.parse_args()


def get_genotype(seqid, field_index):
    parts = seqid.split(".")
    if len(parts) > field_index:
        return parts[field_index]
    return "Unknown"


def run_command(command, output_file=None, shell=False):
    if output_file:
        with open(output_file, "w") as out:
            subprocess.run(
                command,
                stdout=out,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                shell=shell,
            )
    else:
        subprocess.run(
            command,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            shell=shell,
        )


def calculate_identity_matrix(alignment):
    seq_ids = [rec.id for rec in alignment]
    matrix = np.zeros((len(seq_ids), len(seq_ids)))

    for i, j in combinations(range(len(seq_ids)), 2):
        s1 = alignment[i].seq
        s2 = alignment[j].seq

        identity = sum(a == b for a, b in zip(s1, s2)) / len(s1) * 100
        matrix[i, j] = matrix[j, i] = identity

    np.fill_diagonal(matrix, 100)

    return pd.DataFrame(matrix, index=seq_ids, columns=seq_ids)


def build_clusters(seq_ids, identity_df, threshold):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(a, b):
        pa, pb = find(a), find(b)
        if pa != pb:
            parent[pb] = pa

    for i in seq_ids:
        for j in seq_ids:
            if identity_df.loc[i, j] >= threshold:
                union(i, j)

    clusters = {}

    for sid in seq_ids:
        root = find(sid)
        clusters.setdefault(root, []).append(sid)

    return clusters


def build_subclusters(member_list, identity_df):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(a, b):
        pa, pb = find(a), find(b)
        if pa != pb:
            parent[pb] = pa

    for i in range(len(member_list)):
        for j in range(i + 1, len(member_list)):
            if identity_df.loc[member_list[i], member_list[j]] == 100:
                union(member_list[i], member_list[j])

    subclusters = {}

    for sid in member_list:
        root = find(sid)
        subclusters.setdefault(root, []).append(sid)

    return [group for group in subclusters.values() if len(group) >= 2]


def main():
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    if not input_dir.exists():
        sys.exit(f"Error: input directory not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    all_cds = output_dir / "all_cds.fasta"
    all_protein = output_dir / "all_protein.fasta"
    protein_aln = output_dir / "protein_aln.fasta"
    codon_aln = output_dir / "codon_aln.fasta"
    identity_excel = output_dir / "identity.xlsx"

    cluster_dir = output_dir / f"clusters_{int(args.cluster_threshold)}"
    subcluster_dir = output_dir / f"subclusters_{int(args.cluster_threshold)}"

    cluster_dir.mkdir(exist_ok=True)
    subcluster_dir.mkdir(exist_ok=True)

    print("Step 1: Merging CDS files")

    records = []

    for fasta_file in sorted(input_dir.iterdir()):
        if fasta_file.suffix.lower() in [".fa", ".fasta", ".fna"]:
            records.extend(list(SeqIO.parse(fasta_file, "fasta")))

    if not records:
        sys.exit("Error: no CDS FASTA records found.")

    SeqIO.write(records, all_cds, "fasta")
    print(f"Merged CDS sequences: {len(records)}")
    print(f"Output: {all_cds}")

    print("\nStep 2: Translating CDS to protein")

    protein_records = []

    for rec in records:
        prot = rec.seq.translate(to_stop=False)
        protein_records.append(
            SeqRecord(
                prot,
                id=rec.id,
                description=""
            )
        )

    SeqIO.write(protein_records, all_protein, "fasta")
    print(f"Output: {all_protein}")

    print("\nStep 3: Running MAFFT")

    mafft = MafftCommandline(input=str(all_protein))
    mafft.set_parameter("--auto", True)

    stdout, stderr = mafft()

    with open(protein_aln, "w") as out:
        out.write(stdout)

    print(f"Output: {protein_aln}")

    print("\nStep 4: Running PAL2NAL for codon alignment")

    cmd_pal2nal = (
        f'"{args.perl}" "{args.pal2nal}" '
        f'"{protein_aln}" "{all_cds}" -output fasta'
    )

    try:
        run_command(cmd_pal2nal, output_file=codon_aln, shell=True)
    except subprocess.CalledProcessError as error:
        sys.exit(f"Error: PAL2NAL failed.\n{error.stderr}")

    print(f"Output: {codon_aln}")

    print("\nStep 5: Calculating pairwise identity matrix")

    alignment = AlignIO.read(codon_aln, "fasta")
    seq_ids = [rec.id for rec in alignment]

    identity_df = calculate_identity_matrix(alignment)
    identity_df.to_excel(identity_excel)

    print(f"Output: {identity_excel}")

    print(f"\nStep 6: Building clusters at {args.cluster_threshold}% identity")

    clusters = build_clusters(
        seq_ids=seq_ids,
        identity_df=identity_df,
        threshold=args.cluster_threshold
    )

    cluster_map = {}
    singletons = []
    cluster_id = 1

    for root, members in clusters.items():
        if len(members) >= 2:
            cluster_map[cluster_id] = members

            with open(cluster_dir / f"cluster_{cluster_id}.txt", "w") as out:
                out.write("\n".join(members))

            cluster_id += 1
        else:
            singletons.append(members[0])

    with open(cluster_dir / "singletons.txt", "w") as out:
        out.write("\n".join(singletons))

    num_clusters = len(cluster_map)

    print(f"Valid clusters with size >= 2: {num_clusters}")
    print(f"Singleton sequences: {len(singletons)}")

    print("\nStep 7: Grouping cluster members by genotype")

    genotype_cluster_members = {}

    for cid, members in cluster_map.items():
        for sid in members:
            genotype = get_genotype(
                seqid=sid,
                field_index=args.genotype_field
            )

            genotype_cluster_members.setdefault(genotype, {})
            genotype_cluster_members[genotype].setdefault(cid, [])
            genotype_cluster_members[genotype][cid].append(sid)

    print("\nStep 8: Building 100% identity subclusters")

    genotype_cluster_subclusters = {}
    max_subcluster_num = 0

    for genotype, cluster_dict in genotype_cluster_members.items():
        genotype_cluster_subclusters[genotype] = {}#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: run_cds_identity_cluster_pipeline.py

Description:
    Generate codon alignment from CDS sequences, calculate pairwise identity,
    cluster sequences by identity threshold, and summarize cluster/subcluster
    distribution across genotypes.

Workflow:
    1. Merge CDS FASTA files
    2. Translate CDS to protein
    3. Align protein sequences using MAFFT
    4. Generate codon alignment using PAL2NAL
    5. Calculate pairwise identity matrix
    6. Cluster sequences by identity threshold
    7. Build 100% identity subclusters
    8. Summarize cluster and subcluster counts by genotype

Usage:
    python run_cds_identity_cluster_pipeline.py \
        --input_dir input_cds \
        --output_dir output \
        --perl perl \
        --pal2nal pal2nal.pl \
        --cluster_threshold 90

Requirements:
    - Python >= 3.8
    - Biopython
    - pandas
    - numpy
    - MAFFT
    - Perl
    - PAL2NAL
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="CDS identity clustering and subclustering pipeline."
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing CDS FASTA files."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory for output files."
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

    parser.add_argument(
        "--cluster_threshold",
        type=float,
        default=90,
        help="Identity threshold for clustering. Default: 90"
    )

    parser.add_argument(
        "--genotype_field",
        type=int,
        default=1,
        help=(
            "Zero-based field index used to extract genotype name from "
            "dot-separated sequence ID. Default: 1"
        )
    )

    return parser.parse_args()


def get_genotype(seqid, field_index):
    parts = seqid.split(".")
    if len(parts) > field_index:
        return parts[field_index]
    return "Unknown"


def run_command(command, output_file=None, shell=False):
    if output_file:
        with open(output_file, "w") as out:
            subprocess.run(
                command,
                stdout=out,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                shell=shell,
            )
    else:
        subprocess.run(
            command,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            shell=shell,
        )


def calculate_identity_matrix(alignment):
    seq_ids = [rec.id for rec in alignment]
    matrix = np.zeros((len(seq_ids), len(seq_ids)))

    for i, j in combinations(range(len(seq_ids)), 2):
        s1 = alignment[i].seq
        s2 = alignment[j].seq

        identity = sum(a == b for a, b in zip(s1, s2)) / len(s1) * 100
        matrix[i, j] = matrix[j, i] = identity

    np.fill_diagonal(matrix, 100)

    return pd.DataFrame(matrix, index=seq_ids, columns=seq_ids)


def build_clusters(seq_ids, identity_df, threshold):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(a, b):
        pa, pb = find(a), find(b)
        if pa != pb:
            parent[pb] = pa

    for i in seq_ids:
        for j in seq_ids:
            if identity_df.loc[i, j] >= threshold:
                union(i, j)

    clusters = {}

    for sid in seq_ids:
        root = find(sid)
        clusters.setdefault(root, []).append(sid)

    return clusters


def build_subclusters(member_list, identity_df):
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(a, b):
        pa, pb = find(a), find(b)
        if pa != pb:
            parent[pb] = pa

    for i in range(len(member_list)):
        for j in range(i + 1, len(member_list)):
            if identity_df.loc[member_list[i], member_list[j]] == 100:
                union(member_list[i], member_list[j])

    subclusters = {}

    for sid in member_list:
        root = find(sid)
        subclusters.setdefault(root, []).append(sid)

    return [group for group in subclusters.values() if len(group) >= 2]


def main():
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    if not input_dir.exists():
        sys.exit(f"Error: input directory not found: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    all_cds = output_dir / "all_cds.fasta"
    all_protein = output_dir / "all_protein.fasta"
    protein_aln = output_dir / "protein_aln.fasta"
    codon_aln = output_dir / "codon_aln.fasta"
    identity_excel = output_dir / "identity.xlsx"

    cluster_dir = output_dir / f"clusters_{int(args.cluster_threshold)}"
    subcluster_dir = output_dir / f"subclusters_{int(args.cluster_threshold)}"

    cluster_dir.mkdir(exist_ok=True)
    subcluster_dir.mkdir(exist_ok=True)

    print("Step 1: Merging CDS files")

    records = []

    for fasta_file in sorted(input_dir.iterdir()):
        if fasta_file.suffix.lower() in [".fa", ".fasta", ".fna"]:
            records.extend(list(SeqIO.parse(fasta_file, "fasta")))

    if not records:
        sys.exit("Error: no CDS FASTA records found.")

    SeqIO.write(records, all_cds, "fasta")
    print(f"Merged CDS sequences: {len(records)}")
    print(f"Output: {all_cds}")

    print("\nStep 2: Translating CDS to protein")

    protein_records = []

    for rec in records:
        prot = rec.seq.translate(to_stop=False)
        protein_records.append(
            SeqRecord(
                prot,
                id=rec.id,
                description=""
            )
        )

    SeqIO.write(protein_records, all_protein, "fasta")
    print(f"Output: {all_protein}")

    print("\nStep 3: Running MAFFT")

    mafft = MafftCommandline(input=str(all_protein))
    mafft.set_parameter("--auto", True)

    stdout, stderr = mafft()

    with open(protein_aln, "w") as out:
        out.write(stdout)

    print(f"Output: {protein_aln}")

    print("\nStep 4: Running PAL2NAL for codon alignment")

    cmd_pal2nal = (
        f'"{args.perl}" "{args.pal2nal}" '
        f'"{protein_aln}" "{all_cds}" -output fasta'
    )

    try:
        run_command(cmd_pal2nal, output_file=codon_aln, shell=True)
    except subprocess.CalledProcessError as error:
        sys.exit(f"Error: PAL2NAL failed.\n{error.stderr}")

    print(f"Output: {codon_aln}")

    print("\nStep 5: Calculating pairwise identity matrix")

    alignment = AlignIO.read(codon_aln, "fasta")
    seq_ids = [rec.id for rec in alignment]

    identity_df = calculate_identity_matrix(alignment)
    identity_df.to_excel(identity_excel)

    print(f"Output: {identity_excel}")

    print(f"\nStep 6: Building clusters at {args.cluster_threshold}% identity")

    clusters = build_clusters(
        seq_ids=seq_ids,
        identity_df=identity_df,
        threshold=args.cluster_threshold
    )

    cluster_map = {}
    singletons = []
    cluster_id = 1

    for root, members in clusters.items():
        if len(members) >= 2:
            cluster_map[cluster_id] = members

            with open(cluster_dir / f"cluster_{cluster_id}.txt", "w") as out:
                out.write("\n".join(members))

            cluster_id += 1
        else:
            singletons.append(members[0])

    with open(cluster_dir / "singletons.txt", "w") as out:
        out.write("\n".join(singletons))

    num_clusters = len(cluster_map)

    print(f"Valid clusters with size >= 2: {num_clusters}")
    print(f"Singleton sequences: {len(singletons)}")

    print("\nStep 7: Grouping cluster members by genotype")

    genotype_cluster_members = {}

    for cid, members in cluster_map.items():
        for sid in members:
            genotype = get_genotype(
                seqid=sid,
                field_index=args.genotype_field
            )

            genotype_cluster_members.setdefault(genotype, {})
            genotype_cluster_members[genotype].setdefault(cid, [])
            genotype_cluster_members[genotype][cid].append(sid)

    print("\nStep 8: Building 100% identity subclusters")

    genotype_cluster_subclusters = {}
    max_subcluster_num = 0

    for genotype, cluster_dict in genotype_cluster_members.items():
        genotype_cluster_subclusters[genotype] = {}

        for cid, members in cluster_dict.items():
            subclusters = build_subclusters(members, identity_df)
            genotype_cluster_subclusters[genotype][cid] = subclusters
            max_subcluster_num = max(max_subcluster_num, len(subclusters))

            sub_dir = subcluster_dir / f"{genotype}_cluster{cid}"
            sub_dir.mkdir(exist_ok=True)

            for i, group in enumerate(subclusters, start=1):
                with open(sub_dir / f"subcluster_{i}.txt", "w") as out:
                    out.write("\n".join(group))

    print("Subcluster construction completed.")

    print("\nStep 9: Generating cluster summary table")

    genotype_list = sorted(genotype_cluster_members.keys())

    sheet1_rows = []

    for genotype in genotype_list:
        row = [genotype]

        for cid in range(1, num_clusters + 1):
            row.append(
                len(genotype_cluster_members.get(genotype, {}).get(cid, []))
            )

        sheet1_rows.append(row)

    sheet1_df = pd.DataFrame(
        sheet1_rows,
        columns=["Genotype"] + [
            f"Cluster_{i}" for i in range(1, num_clusters + 1)
        ]
    )

    print("Step 10: Generating subcluster summary table")

    sheet2_rows = []

    for genotype in genotype_list:
        for cid in range(1, num_clusters + 1):
            subclusters = genotype_cluster_subclusters.get(
                genotype, {}
            ).get(cid, [])

            counts = [len(group) for group in subclusters]

            while len(counts) < max_subcluster_num:
                counts.append(0)

            label = f"{genotype}_cluster{cid}"
            sheet2_rows.append([label] + counts)

    sheet2_df = pd.DataFrame(
        sheet2_rows,
        columns=["Genotype_Cluster"] + [
            f"subcluster_{i}"
            for i in range(1, max_subcluster_num + 1)
        ]
    )

    print("\nStep 11: Writing summary Excel file")

    out_excel = output_dir / "cluster_subcluster_summary.xlsx"

    with pd.ExcelWriter(out_excel) as writer:
        sheet1_df.to_excel(
            writer,
            sheet_name="Cluster_Count",
            index=False
        )
        sheet2_df.to_excel(
            writer,
            sheet_name="Subcluster_Count",
            index=False
        )

    print("\nPipeline completed successfully.")
    print(f"Identity matrix: {identity_excel}")
    print(f"Cluster directory: {cluster_dir}")
    print(f"Subcluster directory: {subcluster_dir}")
    print(f"Summary Excel: {out_excel}")


if __name__ == "__main__":
    main()

        for cid, members in cluster_dict.items():
            subclusters = build_subclusters(members, identity_df)
            genotype_cluster_subclusters[genotype][cid] = subclusters
            max_subcluster_num = max(max_subcluster_num, len(subclusters))

            sub_dir = subcluster_dir / f"{genotype}_cluster{cid}"
            sub_dir.mkdir(exist_ok=True)

            for i, group in enumerate(subclusters, start=1):
                with open(sub_dir / f"subcluster_{i}.txt", "w") as out:
                    out.write("\n".join(group))

    print("Subcluster construction completed.")

    print("\nStep 9: Generating cluster summary table")

    genotype_list = sorted(genotype_cluster_members.keys())

    sheet1_rows = []

    for genotype in genotype_list:
        row = [genotype]

        for cid in range(1, num_clusters + 1):
            row.append(
                len(genotype_cluster_members.get(genotype, {}).get(cid, []))
            )

        sheet1_rows.append(row)

    sheet1_df = pd.DataFrame(
        sheet1_rows,
        columns=["Genotype"] + [
            f"Cluster_{i}" for i in range(1, num_clusters + 1)
        ]
    )

    print("Step 10: Generating subcluster summary table")

    sheet2_rows = []

    for genotype in genotype_list:
        for cid in range(1, num_clusters + 1):
            subclusters = genotype_cluster_subclusters.get(
                genotype, {}
            ).get(cid, [])

            counts = [len(group) for group in subclusters]

            while len(counts) < max_subcluster_num:
                counts.append(0)

            label = f"{genotype}_cluster{cid}"
            sheet2_rows.append([label] + counts)

    sheet2_df = pd.DataFrame(
        sheet2_rows,
        columns=["Genotype_Cluster"] + [
            f"subcluster_{i}"
            for i in range(1, max_subcluster_num + 1)
        ]
    )

    print("\nStep 11: Writing summary Excel file")

    out_excel = output_dir / "cluster_subcluster_summary.xlsx"

    with pd.ExcelWriter(out_excel) as writer:
        sheet1_df.to_excel(
            writer,
            sheet_name="Cluster_Count",
            index=False
        )
        sheet2_df.to_excel(
            writer,
            sheet_name="Subcluster_Count",
            index=False
        )

    print("\nPipeline completed successfully.")
    print(f"Identity matrix: {identity_excel}")
    print(f"Cluster directory: {cluster_dir}")
    print(f"Subcluster directory: {subcluster_dir}")
    print(f"Summary Excel: {out_excel}")


if __name__ == "__main__":
    main()
