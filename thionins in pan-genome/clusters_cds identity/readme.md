# CDS Identity Clustering and Distance Analysis Pipeline

## Overview

This module performs CDS-level identity analysis, clustering, copy number variation (CNV) quantification, and evolutionary distance calculation for candidate thionin genes across multiple barley genotypes.

The workflow consists of two main steps:

1. Identity-based clustering and CNV analysis  
2. Cluster-level evolutionary distance calculation  

---

## Pipeline Workflow

CDS sequences (20 genotypes)  
↓  
run_cds_identity_cluster_pipeline.py  
↓  
Identity matrix + clusters + subclusters  
↓  
calculate_cluster_distance.py  
↓  
Cluster distance matrix  

---

## Step 1: CDS Identity Clustering

### Script

run_cds_identity_cluster_pipeline.py

### Description

Coding sequences (CDSs) of candidate thionin genes from 20 barley genotypes are processed as follows:

1. Translated into protein sequences  
2. Aligned using MAFFT v7  
3. Converted into codon alignment using PAL2NAL  
4. Pairwise nucleotide identity matrix is calculated  

---

### Clustering Strategy

Sequences with identity ≥ 90% are assigned to the same cluster.

Sequences that do not meet this threshold are defined as singletons.

---

### CNV Quantification

For each genotype:

Number of sequences per cluster = gene copy number

This is used to assess copy number variation (CNV).

---

### Subcluster Definition

Within each genotype:

Sequences with 100% identity are grouped into subclusters.

This reflects:
- exact duplicates  
- highly conserved copies  

---

## Step 2: Cluster Distance Analysis

### Script

calculate_cluster_distance.py

### Description

Cluster-level divergence is quantified using pairwise sequence distances derived from identity.

---

### Distance Definition

d(a,b) = 1 − identity(a,b) / 100

---

### Inter-cluster Distance

For two clusters Ci and Cj:

D(Ci, Cj) = (1 / (|Ci| × |Cj|)) × Σ d(a,b)

where:
- a ∈ Ci  
- b ∈ Cj  

Interpretation:

Average evolutionary distance between two clusters.

---

### Cluster–Singleton Distance

Mean distance between the singleton sequence and all sequences in the cluster.

---

### Singleton–Singleton Distance

Directly obtained from the pairwise distance matrix.

---

## Output Files

### From run_cds_identity_cluster_pipeline.py

- identity.xlsx → Pairwise CDS identity matrix  
- clusters_90/ → Cluster definition files  
- singletons.txt → Singleton sequences  
- cluster_subcluster_summary.xlsx → CNV and subcluster summary  

---

### From calculate_cluster_distance.py

- distance.xlsx → Inter-cluster distance matrix  
- distance_in_cluster.xlsx → Intra-cluster distance  
- analysis_results.pkl → Serialized analysis data  

---

## Biological Interpretation

This pipeline enables:

1. Gene family classification  
   → Identification of thionin clusters  

2. Copy number variation (CNV)  
   → Cluster size reflects gene copy number per genotype  

3. Sequence redundancy  
   → Subclusters represent identical CDS copies  

4. Evolutionary divergence  
   → Cluster distances quantify sequence divergence  

---

## Requirements

### Software

- Python ≥ 3.8  
- MAFFT v7  
- PAL2NAL  
- Perl  

### Python packages

pip install biopython pandas numpy openpyxl tqdm  

---

## Notes

- Identity is calculated at the CDS level  
- Distance is defined as p-distance  
- Clustering threshold (90%) can be adjusted  
- Sequence IDs must be consistent across datasets  

