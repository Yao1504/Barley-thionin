#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: plot_cluster_network.py

Description:
    Visualize cluster-level evolutionary relationships as a network graph
    based on inter-cluster distance matrix.

    This script loads analysis_results.pkl generated from
    calculate_cluster_distance.py and produces publication-quality
    network figures.

Input:
    - analysis_results.pkl

Output:
    - network_final.png
    - network_final.tif

Usage:
    python plot_cluster_network.py \
        --input_dir input_analysis \
        --output_dir output_network \
        --distance_threshold 0.3
"""

import argparse
import pickle
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import warnings

warnings.filterwarnings("ignore")


# ===============================
# 参数解析
# ===============================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot cluster network based on distance matrix"
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing analysis_results.pkl"
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory for network plots"
    )

    parser.add_argument(
        "--distance_threshold",
        type=float,
        default=0.3,
        help="Edge cutoff threshold (default: 0.3)"
    )

    return parser.parse_args()


# ===============================
# Morandi color palette
# ===============================

morandi_colors = [
    "#A8D5BA", "#6EC4B6", "#B5D3E7", "#89BFD7",
    "#6E9AB6", "#97A7D5", "#CAB7D5", "#E7C1B5",
    "#D5B8A8", "#B69A85", "#C5C3BF", "#9E9E9E",
    "#D1E0C3", "#A6B9C7", "#84A7C1", "#6D7F99",
    "#7E9A95", "#93A8AA", "#A8D0D5", "#B2C3D1",
    "#D5D3A8", "#D5C4A7", "#C8A89A", "#B39C84",
    "#808080", "#989898", "#A8B8C5", "#7F8D99",
    "#6E7B80", "#BDC3CA", "#D5D8DF"
]


# ===============================
# 加载数据
# ===============================

def load_analysis_results(input_dir):
    pkl_file = Path(input_dir) / "analysis_results.pkl"

    if not pkl_file.exists():
        raise FileNotFoundError(f"analysis_results.pkl not found: {pkl_file}")

    with open(pkl_file, "rb") as f:
        data = pickle.load(f)

    return data["clusters"], data["single_seqs"], data["cluster_distance_matrix"]


# ===============================
# Scale bar
# ===============================

def calculate_scale_length(pos, G):
    actual = []
    evo = []

    for u, v, data in G.edges(data=True):
        actual.append(np.linalg.norm(pos[u] - pos[v]))
        evo.append(data["distance"])

    ratio = np.sum(actual) / np.sum(evo)
    return ratio * 0.1


def draw_scale_bar(ax, pos, G):

    L = calculate_scale_length(pos, G)
    start_x, start_y = list(pos.values())[0]

    ax.plot([start_x, start_x + L], [start_y, start_y],
            color="black", linewidth=2)

    ax.text(start_x + L/2, start_y,
            "Evolutionary distance = 0.1",
            ha="center", fontsize=10)


# ===============================
# 绘图
# ===============================

def plot_network(clusters, single_seqs, matrix, output_dir, threshold):

    G = nx.Graph()
    all_nodes = list(clusters.keys()) + single_seqs

    node_sizes = {}
    node_colors = {}

    min_display, max_display = 150, 6000

    # cluster nodes
    for i, c in enumerate(clusters):
        size = len(clusters[c])
        scaled = min_display + (size / 600) * (max_display - min_display)
        node_sizes[c] = scaled
        node_colors[c] = morandi_colors[i % len(morandi_colors)]

    # singleton nodes
    offset = len(clusters)
    for i, s in enumerate(single_seqs):
        node_sizes[s] = min_display
        node_colors[s] = morandi_colors[(offset + i) % len(morandi_colors)]

    for n in all_nodes:
        G.add_node(n)

    # edges
    for i, a in enumerate(all_nodes):
        for j, b in enumerate(all_nodes):
            if i < j:
                d = matrix.loc[a, b]
                if 0 < d < threshold:
                    G.add_edge(a, b, distance=d)

    pos = nx.spring_layout(G, seed=42, iterations=300)

    fig, ax = plt.subplots(figsize=(18, 12))

    # edges
    for u, v in G.edges():
        ax.plot([pos[u][0], pos[v][0]],
                [pos[u][1], pos[v][1]],
                color="black", alpha=0.5)

    # nodes
    for n in all_nodes:
        ax.scatter(pos[n][0], pos[n][1],
                   s=node_sizes[n],
                   color=node_colors[n],
                   edgecolors="black")

    ax.axis("off")

    draw_scale_bar(ax, pos, G)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    out_png = output_dir / "network.png"
    out_tif = output_dir / "network.tif"

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_tif, dpi=300, bbox_inches="tight")

    print(f"Saved: {out_png}")
    print(f"Saved: {out_tif}")


# ===============================
# 主程序
# ===============================

def main():
    args = parse_args()

    clusters, single_seqs, matrix = load_analysis_results(args.input_dir)

    plot_network(
        clusters,
        single_seqs,
        matrix,
        args.output_dir,
        args.distance_threshold
    )


if __name__ == "__main__":
    main()
