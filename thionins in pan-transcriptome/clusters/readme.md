# Pan-transcriptome to Pangenome Cluster Mapping

## Overview

This pipeline maps pan-transcriptome sequences to predefined gene clusters derived from a pangenome using BLASTN-based similarity search and downstream processing.

It is designed for identifying thionin gene family members in barley but can be adapted to other gene families or species.

---

## Method Description

To assign barley pan-transcriptome transcripts to thionin clusters previously defined in the pan-genome, BLASTN searches were performed against curated thionin CDS sequences derived from the pan-genome. Because a subset of transcripts is incomplete, candidate sequences were divided into two categories based on transcript length: full-length candidates (≥ 300 bp) and short or incomplete candidates (< 300 bp). For full-length candidates, BLAST matches were retained only if all of the following criteria were satisfied: percent nucleotide identity (pident) ≥ 90%, transcript query coverage (qcov) ≥ 0.7, and CDS subject coverage (scov) ≥ 0.7. For short candidates, which may represent bona fide but truncated transcript isoforms, more stringent identity and query coverage thresholds were applied (pident ≥ 90% and qcov ≥ 0.9), whereas no minimum CDS coverage requirement was imposed in order to avoid excluding valid partial homologs. For each transcript, only the BLAST hit with the highest bitscore that satisfied the above criteria was retained as the best-matching CDS. The cluster identity of the matched CDS was then transferred to the corresponding transcript, thereby assigning each transcript to a predefined thionin cluster, while transcripts lacking valid matches were classified as unassigned. The number of transcripts assigned to each cluster was subsequently quantified across all genotypes. To further characterize redundancy and sequence variation within clusters, transcripts were grouped within each species and cluster based on exact nucleotide sequence identity, and groups containing at least two identical sequences were defined as subclusters. This approach is equivalent to defining subclusters based on 100% pairwise sequence identity and provides a computationally efficient alternative to alignment-based methods such as MAFFT, while preserving the biological interpretation of identical transcript copies within and across genotypes.

---

## Features

- BLAST-based transcript-to-CDS mapping
- Length-aware filtering strategy
- Automatic cluster assignment
- Subcluster detection (100% identical sequences)
- Multi-level outputs (transcript, cluster, subcluster)
- Excel summary reports

---

## Input Requirements

### CDS directory (`--cds_dir`)

Contains pangenome CDS FASTA files:

*.fa  
*.fasta  
*.fna  

---

### Cluster directory (`--cluster_dir`)

cluster_1.txt  
cluster_2.txt  
...  
singleton.txt (optional)  

Each file contains CDS IDs belonging to a cluster.

---

### Transcript FASTA (`--transcript_fasta`)

>transcript_id  
ATGC...

---

## Installation

### Python dependencies

pip install -r requirements.txt

requirements.txt:

numpy  
pandas  
openpyxl  

---

### External dependency

Install BLAST+ and ensure the following commands are available:

makeblastdb  
blastn  

---

## Usage

python map_pantranscriptome_to_pangenome_clusters.py \
    --cds_dir input/all_pangenome_cds \
    --cluster_dir input/cluster_pangenome \
    --transcript_fasta input/all_transcripts.fasta \
    --output_dir output \
    --threads 8

---

## Output Structure

output/  
├── blast_transcripts_vs_cds.tsv  
├── transcript_to_cds_cluster.tsv  
├── transcript_to_cds_cluster.xlsx  
├── cluster_transcripts/  
├── subcluster_transcripts/  
├── subcluster_stats.xlsx  
├── transcript_cluster_subcluster.xlsx  
└── summary.xlsx  

---

## Key Outputs

### Transcript-to-cluster mapping

transcript_to_cds_cluster.xlsx

Includes:

- transcript_id  
- best_cds_id  
- pident / qcov / scov  
- cluster  
- species  

---

### Cluster summary

summary.xlsx (cluster_summary)

- Rows: species  
- Columns: clusters  
- Values: transcript counts  

---

### Subcluster summary

summary.xlsx (subcluster_summary)

- Rows: species × cluster  
- Columns: subcluster sizes (top N)  

---

### Subcluster definitions

subcluster_transcripts/

Each file contains transcript IDs belonging to a subcluster.

---

## Assumptions

- Transcript IDs encode species information  
  (e.g. Morex_chr1_xxx → species = Morex)  
- CDS clusters are predefined  
- Sequence identity is exact (100%)  

---

## Limitations

- BLASTN only (no splice-aware alignment)  
- Subclusters require exact identity (no near-match clustering)  
- Species parsing depends on ID format  
