# CDS Identity Clustering, Distance Analysis, and Visualization Pipeline

## Overview

This module provides an integrated workflow for:

* CDS identity analysis
* Sequence clustering and CNV quantification
* Evolutionary distance calculation
* Network and phylogenetic visualization

The pipeline is designed for systematic analysis of **thionin gene family variation across multiple barley genotypes**.

The workflow consists of three main components:

1. Identity-based clustering and CNV analysis
2. Cluster-level evolutionary distance calculation
3. Visualization of cluster relationships and phylogenetic structure

---

## Pipeline Workflow

CDS sequences (multiple genotypes)
↓
run_cds_identity_cluster_pipeline.py
↓
Identity matrix + clusters + subclusters
↓
calculate_cluster_distance.py
↓
Cluster distance matrix + analysis results
↓
├── plot_cluster_network.py
└── plot_tree_identity_heatmap.py

---

## Step 1: CDS Identity Clustering

### Script

run_cds_identity_cluster_pipeline.py

### Description

Coding sequences (CDSs) of candidate thionin genes are processed through the following steps:

1. Translation of CDS into protein sequences
2. Multiple sequence alignment using MAFFT v7
3. Codon alignment using PAL2NAL
4. Calculation of pairwise nucleotide identity matrix

---

### Clustering Strategy

Sequences with identity ≥ 90% are assigned to the same cluster.

Sequences that do not meet this threshold with any other sequence are defined as **singletons**.

---

### CNV Quantification

For each genotype:

* The number of sequences assigned to each cluster is counted
* This provides a quantitative measure of **copy number variation (CNV)**

---

### Subclustering

Within each genotype:

* Sequences with **100% identity** are grouped into subclusters
* The number of subclusters reflects **fine-scale duplication patterns**

---

### Output

* identity.xlsx
* clusters_90/
* subclusters_90/
* cluster_subcluster_summary.xlsx

---

## Step 2: Cluster Distance Calculation

### Script

calculate_cluster_distance.py

### Description

Pairwise evolutionary distances are derived from the identity matrix.

---

### Distance Definition

Pairwise distance between two sequences:

d(a,b) = 1 − identity(a,b) / 100

---

### Cluster Distance

For two clusters Ci and Cj:

D(Ci, Cj) = mean pairwise distance between all sequences in Ci and Cj

---

### Special Cases

* Cluster vs singleton:
  mean distance between singleton and all members of the cluster

* Singleton vs singleton:
  directly obtained from pairwise distance matrix

---

### Output

* distance.xlsx
* distance_in_cluster.xlsx
* analysis_results.pkl

---

## Step 3: Cluster Network Visualization

### Script

plot_cluster_network.py

### Description

Constructs a network representation of cluster relationships.

---

### Network Definition

Nodes:

* clusters (size proportional to number of sequences)
* singletons

Edges:

* defined by evolutionary distance threshold

---

### Visualization Features

* Force-directed layout (NetworkX)
* Node size scaling
* Cluster color mapping
* Evolutionary distance scale bar

---

### Output

* network_final_2panels.png
* network_final_2panels.tif

---

## Step 4: Phylogenetic Tree and Identity Heatmap

### Script

plot_tree_identity_heatmap.py

### Description

Generates a combined visualization of:

* Phylogenetic tree
* Sequence identity heatmap
* Cluster annotation

---

### Key Features

* Identity matrix reordered according to tree topology

* Integrated visualization layout:

  * Tree (left)
  * Heatmap (center)
  * Cluster annotation (right)

* Separate panels are also exported

---

### Output

* tree_heatmap_cluster.png
* tree_only.png
* heatmap_only.png
* cluster_only.png
* cluster_assignment.csv
* gene_id_order.txt

---

## Requirements

### Python packages

* biopython
* pandas
* numpy
* matplotlib
* networkx
* openpyxl

---

### External tools

* MAFFT v7
* PAL2NAL
* Perl
* IQ-TREE

---

## Notes

* Sequence IDs must be consistent across:

  * FASTA files
  * identity matrix
  * phylogenetic tree

* Cluster files must follow naming convention:

cluster_1.txt
cluster_2.txt

* Identity matrix must be symmetric and complete

---

## Summary

This pipeline enables:

* Identification of thionin gene clusters
* Quantification of CNV across genotypes
* Assessment of sequence divergence
* Visualization of evolutionary relationships

It provides a reproducible framework for integrating **sequence identity, clustering, and phylogenetic analysis** in gene family studies.
