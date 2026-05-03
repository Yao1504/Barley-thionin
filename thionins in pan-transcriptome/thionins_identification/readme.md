# Thionin Identification from Barley Pan-Transcriptome

## Overview

This pipeline identifies and curates thionin genes from barley pan-transcriptome assemblies across 20 genotypes.

It integrates sequence similarity search, genomic coordinate validation, transcript sequence refinement, redundancy removal, and ORF-based functional classification.

---

## Pipeline Workflow

Pan-transcriptome assemblies (20 genotypes)---run_thionin_pantranscriptome.py---Candidate thionin transcripts (BLAST-based)---classify_transcripts_by_chr.py---Chromosome-based grouping---
run_mafft_transcript_alignment.py---Sequence alignment (MAFFT)---remove_duplicate_transcripts.py---Non-redundant transcript set---translate_transcripts_to_proteins.py---Final CDS + protein sequences + pseudogene annotation

---

## Step 1: Thionin Candidate Identification

### Script

`run_thionin_pantranscriptome.py`

### Description

Transcriptome assemblies and corresponding expression datasets for 20 barley genotypes were retrieved from the barley pan-transcriptomic resource (Guo et al., 2025).

Curated thionin CDSs derived from the barley pan-genome were used as query sequences.

For each genotype:

* A BLASTn database was constructed from its pan-transcriptome assembly
* Sequence similarity search was performed using:

  * E-value ≤ 1e-5
  * Alignment length ≥ 200 bp
  * Mismatch ≤ 10

### Output

* Filtered BLAST results
* Candidate thionin transcript IDs
* Extracted transcript sequences

---

## Step 2: Genomic Coordinate Validation

### Script

`run_thionin_pantranscriptome.py` (internal steps)

### Description

To ensure accurate genomic assignment:

* Transcript coordinates were extracted from pan-transcriptome **GTF**
* Genomic coordinates were retrieved from pan-genome **GFF**
* Candidates were retained only if:

  * Same chromosome
  * Same strand
  * Overlapping genomic coordinates

### Output

* High-confidence thionin transcripts per genotype

---

## Step 3: Chromosome-Based Classification

### Script

`classify_transcripts_by_chr.py`

### Description

All validated transcript sequences were grouped based on chromosome information extracted from sequence IDs.

* Chromosomes: 1H–7H
* Unassigned sequences → `Un`
* Combined sets generated:

  * `6H + Un`
  * `7H + Un`

### Output

* Chromosome-specific FASTA files
* Combined chromosome datasets

---

## Step 4: Sequence Alignment and Cleaning

### Script

`run_mafft_transcript_alignment.py`

### Description

Candidate sequences were aligned using MAFFT v7:

* Automatic algorithm selection (`--auto`)
* Multi-threaded execution

### Output

* Multiple sequence alignments for each chromosome group

---

## Step 5: Redundancy Removal

### Script

`remove_duplicate_transcripts.py`

### Description

Redundant transcript isoforms were removed based on sequence identity:

* Group sequences by gene ID
* Compare sequences using hash (sequence identity)
* Remove identical isoforms
* Retain:

  * One representative per identical group
  * All distinct isoforms

### Output

* Non-redundant transcript FASTA files

---

## Step 6: ORF Detection and Translation

### Script

`translate_transcripts_to_proteins.py`

### Description

To validate coding potential:

1. ORF detection:

   * Scan all 3 reading frames
   * Start codon: ATG
   * Stop codons: TAA, TAG, TGA

2. Classification:

   * Valid ORF → functional gene
   * No ORF → putative pseudogene

3. Translation:

   * Functional ORF → translated
   * Pseudogene → full-length translation

### Output

* Protein sequences
* Cleaned CDS sequences
* ORF report (per sequence):

  * ORF coordinates
  * Functional status
  * Pseudogene annotation

---

## Key Features

* End-to-end transcriptome-based gene discovery pipeline
* Integrates transcriptome and genome annotation
* Handles isoforms and redundancy explicitly
* Supports pseudogene identification
* Modular and reproducible

---

## Requirements

* Python ≥ 3.8
* MAFFT v7
* BLAST+

Python packages:

* Biopython
* pandas

---

## Example Usage

```bash
# Step 1: BLAST-based identification
python run_thionin_pantranscriptome.py

# Step 2: Chromosome classification
python classify_transcripts_by_chr.py

# Step 3: Alignment
python run_mafft_transcript_alignment.py

# Step 4: Remove redundancy
python remove_duplicate_transcripts.py

# Step 5: Translation and ORF detection
python translate_transcripts_to_proteins.py
```

---

## Notes

* Manual trimming may be required after alignment to remove:

  * intronic remnants
  * incomplete exons
  * abnormal extensions

* Multiple valid isoforms are retained for downstream analysis

