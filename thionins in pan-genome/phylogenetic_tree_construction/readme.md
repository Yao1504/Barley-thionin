# Alignment, Curation, and Phylogenetic Analysis Pipeline

## Overview

This module performs sequence alignment, manual curation, codon-level alignment, and phylogenetic tree reconstruction for candidate thionin genes.

It integrates two main scripts:

* `run_alignment_trimming_pipeline.py`
* `run_iqtree_codon_tree_pipeline.py`

These steps ensure structural accuracy of candidate sequences and enable robust phylogenetic inference across barley genotypes.


## Workflow

```text
Candidate protein sequences
        ↓
MAFFT alignment
        ↓
Manual curation & trimming
        ↓
CDS refinement
        ↓
Codon alignment
        ↓
Model selection (ModelFinder)
        ↓
ML tree construction (IQ-TREE)
        ↓
Bootstrap support evaluation
```

## Step 1: Alignment and structural curation

### Script

`run_alignment_trimming_pipeline.py`

### Description

Candidate protein sequences were aligned and curated to ensure structural integrity and biological consistency.

### Corresponding method description

> To assess sequence integrity, candidate protein sequences were aligned using MAFFT v7 and manually curated by comparison with the characteristic thionin structural features observed on different chromosomes in the Morex V3 reference genome.
> Misaligned or irregular regions were trimmed to restore typical thionin-like architecture.
> For candidates located on unanchored scaffolds (chrUn), chromosomal origin was inferred based on sequence similarity and structural characteristics.
> The curated protein sequences were subsequently used as templates to refine their corresponding coding sequences (CDSs), ensuring consistency between CDSs and protein structures.

### Function

* Performs multiple sequence alignment using MAFFT
* Enables identification of structurally inconsistent regions
* Supports manual trimming and curation
* Provides curated protein sequences for downstream CDS refinement


## Step 2: Codon alignment and phylogenetic analysis

### Script

`run_iqtree_codon_tree_pipeline.py`

### Description

Trimmed CDS sequences were used to construct phylogenetic trees using a maximum-likelihood framework.

### Corresponding method description

> The trimmed CDSs of candidate thionin genes identified across the 20 barley genotypes were subjected to ML model testing in IQ-TREE v3.0.1.
> The best-fit nucleotide substitution model was selected using ModelFinder based on the BIC, with TVM+R3 identified as the optimal model.
> ML phylogenetic trees were subsequently reconstructed under the selected model, and branch support was evaluated using 1,000 ultrafast bootstrap replicates (UFBoot).

### Function

* Performs model selection using ModelFinder
* Extracts best-fit substitution model
* Constructs maximum-likelihood phylogenetic trees
* Performs bootstrap analysis (UFBoot or standard bootstrap)
* Supports checkpoint-based resume for long computations


## Input

| Input                 | Description                         |
| --------------------- | ----------------------------------- |
| Protein FASTA         | Candidate thionin protein sequences |
| CDS FASTA             | Corresponding coding sequences      |
| Trimmed CDS alignment | Curated and trimmed CDS sequences   |


## Output

### Alignment and curation

```text
aligned_protein.fasta
trimmed_protein.fasta
refined_cds.fasta
```

### Phylogenetic analysis

```text
all_cds_trimmed.fasta
thionin.treefile
thionin.iqtree
best_model_report.txt
```

## Key methods

### Multiple sequence alignment

Performed using **MAFFT**
High-quality alignments are critical for preserving structural features of thionin proteins.


### Codon alignment

Protein-guided codon alignment ensures:

* Correct reading frame preservation
* Accurate evolutionary inference
* Compatibility with downstream phylogenetic models


### Model selection

Performed using ModelFinder implemented in **IQ-TREE**:

* Model selection based on Bayesian Information Criterion (BIC)
* Automatic identification of best-fit substitution model


### Phylogenetic reconstruction

* Maximum-likelihood (ML) framework
* Bootstrap support (1,000 replicates)
* Supports both:

  * Standard bootstrap (`-b`)
  * Ultrafast bootstrap (`-bb`)


## Notes

* Manual curation is an essential step and not automated
* Chromosomal inference for chrUn sequences relies on biological interpretation
* Codon alignment is required for accurate evolutionary analysis
* Checkpointing enables recovery from interrupted runs


## Relationship to the thionin project

This module corresponds to the final analytical stage of the thionin pipeline:

* Structural validation of candidate proteins
* Refinement of CDS sequences
* Evolutionary analysis across barley genotypes

It produces high-confidence phylogenetic relationships suitable for publication.


## Requirements

* Python >= 3.8
* MAFFT
* IQ-TREE
* trimAl
* PAL2NAL (optional, if codon alignment is used)

Install dependencies:

```bash
pip install biopython pandas
```


## Summary

This pipeline ensures that:

* Candidate sequences are structurally valid
* CDS and protein sequences are consistent
* Phylogenetic inference is robust and reproducible

