#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: plot_tree_identity_heatmap.py

Description:
    Generate an integrated visualization combining:
        1. phylogenetic tree
        2. pairwise sequence identity heatmap
        3. cluster annotation bar

Input:
    - Newick tree file
    - Pairwise identity matrix in Excel format
    - Directory containing cluster_*.txt files

Output:
    - tree_heatmap_cluster.png
    - tree_heatmap_cluster.pdf
    - tree_only.png / .pdf
    - heatmap_only.png / .pdf
    - cluster_only.png / .pdf
    - gene_id_order.txt
    - cluster_assignment.csv
    - reordered_identity_matrix.csv

Usage:
    python plot_tree_identity_heatmap.py \
        --tree tree.nwk \
        --identity identity.xlsx \
        --cluster_dir clusters_90 \
        --output_dir output_tree_heatmap
"""

import argparse
import os
import glob
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
from Bio import Phylo


CLUSTER_COLORS = [
    "#EFD3AC", "#C9B4C7", "#EED1CC", "#D4E3DD", "#D6D9B9",
    "#A1B0AD", "#F4EEAC", "#C9DCC4", "#E69191", "#A09952",
    "#A8D5BA", "#6EC4B6", "#B5D3E7", "#89BFD7", "#6E9AB6",
    "#97A7D5", "#CAB7D5", "#E7C1B5", "#D5B8A8", "#B69A85"
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot phylogenetic tree, identity heatmap, and cluster annotation."
    )

    parser.add_argument(
        "--tree",
        required=True,
        help="Input phylogenetic tree file in Newick format."
    )

    parser.add_argument(
        "--identity",
        required=True,
        help="Pairwise sequence identity matrix in Excel format."
    )

    parser.add_argument(
        "--cluster_dir",
        required=True,
        help="Directory containing cluster_*.txt files."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory."
    )

    parser.add_argument(
        "--tree_scale",
        type=float,
        default=2.0,
        help="Scale factor for tree branch lengths in the plot. Default: 2.0"
    )

    parser.add_argument(
        "--vmin",
        type=float,
        default=40,
        help="Minimum value for identity heatmap color scale. Default: 40"
    )

    parser.add_argument(
        "--vmax",
        type=float,
        default=100,
        help="Maximum value for identity heatmap color scale. Default: 100"
    )

    return parser.parse_args()


def read_cluster_files(cluster_dir):
    """
    Read cluster_*.txt files and return a gene-to-cluster dictionary.
    """
    cluster_dir = Path(cluster_dir)

    if not cluster_dir.exists():
        raise FileNotFoundError(f"Cluster directory not found: {cluster_dir}")

    clusters = {}

    cluster_files = sorted(glob.glob(str(cluster_dir / "cluster*.txt")))

    if not cluster_files:
        cluster_files = sorted(glob.glob(str(cluster_dir / "*.txt")))

    if not cluster_files:
        raise FileNotFoundError(f"No cluster txt files found in: {cluster_dir}")

    for idx, cluster_file in enumerate(cluster_files, start=1):
        cluster_name = Path(cluster_file).stem

        with open(cluster_file, "r", encoding="utf-8") as f:
            genes = [line.strip() for line in f if line.strip()]

        for gene in genes:
            clusters[gene] = cluster_name

        print(f"Loaded {cluster_name}: {len(genes)} genes")

    return clusters


def get_tree_leaf_order(tree_file):
    """
    Read Newick tree and return leaf order.
    """
    with open(tree_file, "r", encoding="utf-8") as f:
        tree_content = f.read()

    tree = Phylo.read(StringIO(tree_content), "newick")
    leaves = [leaf.name for leaf in tree.get_terminals()]

    return leaves, tree


def scale_tree_branches(tree, scale_factor):
    """
    Scale tree branch lengths for visualization.
    """
    for clade in tree.find_clades():
        if clade.branch_length:
            clade.branch_length *= scale_factor

    return tree


def reorder_matrix_by_tree(identity_df, tree_leaves, clusters):
    """
    Reorder identity matrix according to tree leaf order.
    """
    available_leaves = [
        leaf for leaf in tree_leaves
        if leaf in identity_df.index and leaf in identity_df.columns
    ]

    if not available_leaves:
        raise ValueError(
            "No tree leaf IDs were found in the identity matrix. "
            "Please check whether sequence IDs match."
        )

    reordered_df = identity_df.loc[available_leaves, available_leaves]
    cluster_info = [
        clusters.get(leaf, "No_Cluster")
        for leaf in available_leaves
    ]

    print(f"Tree leaves: {len(tree_leaves)}")
    print(f"Matched leaves in identity matrix: {len(available_leaves)}")
    print(f"Genes with cluster assignment: {sum(c != 'No_Cluster' for c in cluster_info)}")

    return reordered_df, cluster_info, available_leaves


def draw_tree_on_axis(tree, ax, leaves_order):
    """
    Draw phylogenetic tree on a matplotlib axis.
    """
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda x: "",
        show_confidence=False,
        branch_labels=lambda x: "",
    )

    ax.set_xlim(0, tree.total_branch_length() * 1.1)
    ax.set_ylim(-0.5, len(leaves_order) - 0.5)
    ax.invert_yaxis()

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)

    for line in ax.lines:
        line.set_linewidth(0.5)


def build_cluster_color_map(cluster_info):
    """
    Assign colors to clusters.
    """
    unique_clusters = sorted(
        cluster for cluster in set(cluster_info)
        if cluster != "No_Cluster"
    )

    return {
        cluster: CLUSTER_COLORS[i % len(CLUSTER_COLORS)]
        for i, cluster in enumerate(unique_clusters)
    }


def save_gene_order(leaves_order, output_dir):
    output_file = Path(output_dir) / "gene_id_order.txt"

    with open(output_file, "w", encoding="utf-8") as f:
        for gene_id in leaves_order:
            f.write(f"{gene_id}\n")

    print(f"Gene order saved to: {output_file}")


def save_separate_plots(
    tree,
    identity_matrix,
    cluster_info,
    leaves_order,
    cluster_color_map,
    output_dir,
    vmin,
    vmax,
):
    """
    Save tree, heatmap, and cluster annotation separately.
    """
    output_dir = Path(output_dir)
    n_genes = len(leaves_order)

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "identity_cmap",
        ["#F2F7FB", "#9CB6DD", "#4A5989"],
        N=256,
    )

    # Tree only
    fig_tree = plt.figure(figsize=(8, max(4, n_genes * 0.05)))
    ax_tree = fig_tree.add_subplot(111)
    draw_tree_on_axis(tree, ax_tree, leaves_order)
    ax_tree.set_title("Phylogenetic Tree", fontsize=14, pad=20)

    fig_tree.savefig(output_dir / "tree_only.png", dpi=300, bbox_inches="tight")
    fig_tree.savefig(output_dir / "tree_only.pdf", bbox_inches="tight")
    plt.close(fig_tree)

    # Heatmap only
    heatmap_size = min(30, max(6, n_genes * 0.03))
    fig_heatmap = plt.figure(figsize=(heatmap_size, heatmap_size))
    ax_heatmap = fig_heatmap.add_subplot(111)

    im = ax_heatmap.imshow(
        identity_matrix.values.astype(float),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        aspect="auto",
        interpolation="nearest",
    )

    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])
    ax_heatmap.set_title("Sequence Identity Matrix", fontsize=14, pad=20)

    cbar = plt.colorbar(im, ax=ax_heatmap, fraction=0.046, pad=0.04)
    cbar.set_label("Sequence Identity (%)", rotation=270, labelpad=20)

    fig_heatmap.savefig(output_dir / "heatmap_only.png", dpi=300, bbox_inches="tight")
    fig_heatmap.savefig(output_dir / "heatmap_only.pdf", bbox_inches="tight")
    plt.close(fig_heatmap)

    # Cluster only
    fig_cluster = plt.figure(figsize=(2, max(4, n_genes * 0.05)))
    ax_cluster = fig_cluster.add_subplot(111)

    for i, cluster in enumerate(cluster_info):
        color = cluster_color_map.get(cluster, "#F0F0F0")
        rect = Rectangle(
            (0, i),
            1,
            1,
            facecolor=color,
            edgecolor="white",
            linewidth=0.3,
        )
        ax_cluster.add_patch(rect)

    ax_cluster.set_xlim(0, 1)
    ax_cluster.set_ylim(0, len(cluster_info))
    ax_cluster.set_xticks([])
    ax_cluster.set_yticks([])
    ax_cluster.invert_yaxis()
    ax_cluster.set_title("Cluster Assignment", fontsize=14, pad=20)

    fig_cluster.savefig(output_dir / "cluster_only.png", dpi=300, bbox_inches="tight")
    fig_cluster.savefig(output_dir / "cluster_only.pdf", bbox_inches="tight")
    plt.close(fig_cluster)


def create_combined_plot(
    tree,
    identity_matrix,
    cluster_info,
    leaves_order,
    output_dir,
    vmin,
    vmax,
):
    """
    Create integrated tree + heatmap + cluster annotation plot.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_genes = len(leaves_order)
    gene_size = 0.1
    heatmap_size = n_genes * gene_size

    fig_width = heatmap_size * 1.5
    fig_height = heatmap_size

    max_size = 50
    if fig_width > max_size or fig_height > max_size:
        scale_factor = max_size / max(fig_width, fig_height)
        fig_width *= scale_factor
        fig_height *= scale_factor

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.8, 2, 0.1], wspace=0.02)

    ax_tree = fig.add_subplot(gs[0])
    ax_tree.set_title("Phylogenetic Tree", fontsize=12, pad=10)
    draw_tree_on_axis(tree, ax_tree, leaves_order)

    ax_heatmap = fig.add_subplot(gs[1])
    ax_heatmap.set_title("Sequence Identity Matrix", fontsize=12, pad=10)

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "identity_cmap",
        ["#F2F7FB", "#9CB6DD", "#4A5989"],
        N=256,
    )

    im = ax_heatmap.imshow(
        identity_matrix.values.astype(float),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        aspect="auto",
        interpolation="nearest",
    )

    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])

    ax_cluster = fig.add_subplot(gs[2])
    ax_cluster.set_title("Cluster", fontsize=10, pad=10)

    cluster_color_map = build_cluster_color_map(cluster_info)

    for i, cluster in enumerate(cluster_info):
        color = cluster_color_map.get(cluster, "#F0F0F0")
        rect = Rectangle(
            (0, i),
            1,
            1,
            facecolor=color,
            edgecolor="white",
            linewidth=0.1,
        )
        ax_cluster.add_patch(rect)

    ax_cluster.set_xlim(0, 1)
    ax_cluster.set_ylim(0, len(cluster_info))
    ax_cluster.set_xticks([])
    ax_cluster.set_yticks([])
    ax_cluster.invert_yaxis()

    cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])
    cbar = plt.colorbar(im, cax=cbar_ax)
    cbar.set_label("Sequence Identity (%)", rotation=270, labelpad=15, fontsize=10)

    if cluster_color_map:
        legend_elements = [
            Rectangle(
                (0, 0),
                1,
                1,
                facecolor=color,
                edgecolor="black",
                label=cluster,
            )
            for cluster, color in cluster_color_map.items()
        ]

        n_clusters = len(legend_elements)
        n_cols = min(5, n_clusters)
        n_rows = (n_clusters + n_cols - 1) // n_cols
        legend_height = 0.05 * n_rows
        bottom_margin = 0.05 + legend_height

        legend_ax = fig.add_axes([0.25, 0.01, 0.5, legend_height])
        legend_ax.axis("off")
        legend = legend_ax.legend(
            handles=legend_elements,
            loc="center",
            ncol=n_cols,
            fontsize=8,
            frameon=True,
            handlelength=1.5,
            handleheight=1.5,
        )
        legend.set_title("Cluster Legend", prop={"size": 10})
    else:
        bottom_margin = 0.05

    plt.subplots_adjust(
        left=0.05,
        right=0.93,
        top=0.95,
        bottom=bottom_margin,
        wspace=0.02,
    )

    fig.savefig(output_dir / "tree_heatmap_cluster.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / "tree_heatmap_cluster.pdf", bbox_inches="tight")
    plt.close(fig)

    save_separate_plots(
        tree=tree,
        identity_matrix=identity_matrix,
        cluster_info=cluster_info,
        leaves_order=leaves_order,
        cluster_color_map=cluster_color_map,
        output_dir=output_dir,
        vmin=vmin,
        vmax=vmax,
    )

    return cluster_color_map


def main():
    args = parse_args()

    tree_file = Path(args.tree)
    identity_file = Path(args.identity)
    cluster_dir = Path(args.cluster_dir)
    output_dir = Path(args.output_dir)

    if not tree_file.exists():
        raise FileNotFoundError(f"Tree file not found: {tree_file}")

    if not identity_file.exists():
        raise FileNotFoundError(f"Identity matrix not found: {identity_file}")

    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading cluster information...")
    clusters = read_cluster_files(cluster_dir)

    print("Loading phylogenetic tree...")
    tree_leaves, tree = get_tree_leaf_order(tree_file)
    tree = scale_tree_branches(tree, scale_factor=args.tree_scale)

    print("Loading identity matrix...")
    identity_df = pd.read_excel(identity_file, index_col=0)

    print("Reordering identity matrix by tree order...")
    reordered_matrix, cluster_info, available_leaves = reorder_matrix_by_tree(
        identity_df=identity_df,
        tree_leaves=tree_leaves,
        clusters=clusters,
    )

    print("Saving gene order and metadata...")
    save_gene_order(available_leaves, output_dir)

    cluster_info_df = pd.DataFrame({
        "Gene": available_leaves,
        "Cluster": cluster_info,
    })

    cluster_info_df.to_csv(output_dir / "cluster_assignment.csv", index=False)
    reordered_matrix.to_csv(output_dir / "reordered_identity_matrix.csv")

    print("Creating visualization...")
    create_combined_plot(
        tree=tree,
        identity_matrix=reordered_matrix,
        cluster_info=cluster_info,
        leaves_order=available_leaves,
        output_dir=output_dir,
        vmin=args.vmin,
        vmax=args.vmax,
    )

    print("Analysis completed.")
    print(f"Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
