#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: calculate_cluster_distance.py

Description:
    Convert pairwise sequence identity into distance and calculate
    inter-cluster and intra-cluster distance matrices.

Input:
    1. Pairwise identity matrix in Excel format
    2. Directory containing cluster_*.txt files
    3. Optional singletons.txt file

Output:
    1. distance.xlsx
       Inter-cluster distance matrix

    2. distance_in_cluster.xlsx
       Intra-cluster distance matrices, one sheet per cluster

    3. analysis_results.pkl
       Pickled object containing clusters, singletons, identity matrix,
       sequence-level distance matrix, and cluster-level distance matrix

Usage:
    python calculate_cluster_distance.py \
        --identity_file identity.xlsx \
        --cluster_dir clusters_90 \
        --output_dir output_distance
"""

import argparse
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate inter-cluster and intra-cluster distances from identity matrix."
    )

    parser.add_argument(
        "--identity_file",
        required=True,
        help="Pairwise identity matrix in Excel format."
    )

    parser.add_argument(
        "--cluster_dir",
        required=True,
        help="Directory containing cluster_*.txt files and optional singletons.txt."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory."
    )

    return parser.parse_args()


def load_identity_matrix(identity_file):
    print("Loading identity matrix...")
    identity_matrix = pd.read_excel(identity_file, index_col=0)
    print(f"Loaded identity matrix with {identity_matrix.shape[0]} sequences.")
    return identity_matrix


def load_cluster_members(cluster_dir):
    print("Loading cluster members...")

    cluster_dir = Path(cluster_dir)
    clusters = {}

    cluster_files = sorted(cluster_dir.glob("cluster_*.txt"))
    print(f"Detected {len(cluster_files)} cluster files.")

    for file in cluster_files:
        cluster_name = file.stem

        with open(file, "r") as f:
            members = [line.strip() for line in f if line.strip()]

        clusters[cluster_name] = members
        print(f"Loaded {cluster_name}: {len(members)} members")

    singletons_file = cluster_dir / "singletons.txt"

    if singletons_file.exists():
        with open(singletons_file, "r") as f:
            single_seqs = [line.strip() for line in f if line.strip()]

        print(f"Loaded singletons: {len(single_seqs)} sequences")
    else:
        single_seqs = []
        print("No singletons.txt found. Singleton set is empty.")

    print(f"Total: {len(clusters)} clusters and {len(single_seqs)} singletons.")

    return clusters, single_seqs


def calculate_distance_from_identity(identity_matrix):
    print("Converting identity matrix to distance matrix...")
    distance_matrix = 1 - identity_matrix / 100
    print("Distance matrix calculated.")
    return distance_matrix


def calculate_cluster_distances(clusters, single_seqs, distance_matrix):
    print("Calculating inter-cluster distance matrix...")

    all_elements = list(clusters.keys()) + single_seqs
    n = len(all_elements)

    cluster_distance_matrix = pd.DataFrame(
        np.zeros((n, n)),
        index=all_elements,
        columns=all_elements,
    )

    for i, elem1 in enumerate(tqdm(all_elements, desc="Inter-cluster distance")):
        for j, elem2 in enumerate(all_elements):
            if i == j:
                cluster_distance_matrix.iloc[i, j] = 0
                continue

            distances = []

            if elem1 in clusters:
                members1 = clusters[elem1]
            else:
                members1 = [elem1]

            if elem2 in clusters:
                members2 = clusters[elem2]
            else:
                members2 = [elem2]

            for m1 in members1:
                for m2 in members2:
                    if m1 in distance_matrix.index and m2 in distance_matrix.columns:
                        distances.append(distance_matrix.loc[m1, m2])

            if distances:
                cluster_distance_matrix.iloc[i, j] = np.mean(distances)

    print("Inter-cluster distance matrix calculated.")
    return cluster_distance_matrix


def calculate_intra_cluster_distances(clusters, distance_matrix, output_file):
    print("Calculating intra-cluster distance matrices...")

    with pd.ExcelWriter(output_file) as writer:
        for cluster_name, members in tqdm(clusters.items(), desc="Intra-cluster distance"):
            valid_members = [
                member for member in members
                if member in distance_matrix.index
            ]

            if len(valid_members) > 1:
                dist_df = distance_matrix.loc[valid_members, valid_members]
                sheet_name = cluster_name[:31]
                dist_df.to_excel(writer, sheet_name=sheet_name)

    print(f"Intra-cluster distance matrices saved to: {output_file}")


def save_analysis_results(
    clusters,
    single_seqs,
    identity_matrix,
    distance_matrix,
    cluster_distance_matrix,
    output_dir,
):
    print("Saving analysis results...")

    output_dir = Path(output_dir)

    analysis_data = {
        "clusters": clusters,
        "single_seqs": single_seqs,
        "identity_matrix": identity_matrix,
        "distance_matrix": distance_matrix,
        "cluster_distance_matrix": cluster_distance_matrix,
    }

    with open(output_dir / "analysis_results.pkl", "wb") as f:
        pickle.dump(analysis_data, f)

    cluster_distance_matrix.to_excel(output_dir / "distance.xlsx")

    print(f"Inter-cluster distance matrix saved to: {output_dir / 'distance.xlsx'}")
    print(f"Full analysis object saved to: {output_dir / 'analysis_results.pkl'}")


def main():
    args = parse_args()

    identity_file = Path(args.identity_file)
    cluster_dir = Path(args.cluster_dir)
    output_dir = Path(args.output_dir)

    if not identity_file.exists():
        raise FileNotFoundError(f"Identity file not found: {identity_file}")

    if not cluster_dir.exists():
        raise FileNotFoundError(f"Cluster directory not found: {cluster_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    start_time = time.time()

    print("=" * 60)
    print("Cluster distance calculation started")
    print("=" * 60)

    identity_matrix = load_identity_matrix(identity_file)
    clusters, single_seqs = load_cluster_members(cluster_dir)
    distance_matrix = calculate_distance_from_identity(identity_matrix)

    cluster_distance_matrix = calculate_cluster_distances(
        clusters=clusters,
        single_seqs=single_seqs,
        distance_matrix=distance_matrix,
    )

    calculate_intra_cluster_distances(
        clusters=clusters,
        distance_matrix=distance_matrix,
        output_file=output_dir / "distance_in_cluster.xlsx",
    )

    save_analysis_results(
        clusters=clusters,
        single_seqs=single_seqs,
        identity_matrix=identity_matrix,
        distance_matrix=distance_matrix,
        cluster_distance_matrix=cluster_distance_matrix,
        output_dir=output_dir,
    )

    elapsed = time.time() - start_time

    print("=" * 60)
    print("Cluster distance calculation completed")
    print(f"Elapsed time: {elapsed:.2f} seconds")
    print("=" * 60)


if __name__ == "__main__":
    main()
