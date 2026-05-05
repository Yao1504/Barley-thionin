# qRT-PCR Delta-Delta Ct Analysis

## Overview

This script analyses qRT-PCR Ct data using the Delta-Delta Ct method and generates publication-style grouped boxplots for each reference gene.

It was designed for barley gene expression analysis under untreated, clip-cage control, and aphid infestation conditions, but can be adapted to other plant species, treatments, or reference genes.

---

## Method Description

Raw Ct values are first separated into target-gene sample measurements and reference-gene measurements according to the `type` column. For each genotype, treatment, biological repeat, and reference gene, technical replicates of the reference gene are averaged to obtain the mean reference Ct value. The Delta Ct value is then calculated for each sample as the difference between the target-gene Ct and the corresponding mean reference-gene Ct.

Technical replicates within each biological repeat are subsequently averaged to obtain one Delta Ct value per genotype, treatment, reference gene, and biological repeat. Relative expression is calculated using the Delta-Delta Ct method. By default, the global calibrator is defined as `Akashinriki` under `Untreated control`. For each reference gene, the mean Delta Ct of this calibrator group is used as the baseline. Delta-Delta Ct is calculated as the difference between each sample Delta Ct and the calibrator mean Delta Ct.

Relative expression is reported as `2^-DeltaDeltaCt`, and log2 fold change is calculated as `-DeltaDeltaCt`. Statistical summaries, including mean, median, standard error of the mean, standard deviation, and replicate number, are generated for each genotype, treatment, and reference gene. Publication-style grouped boxplots are generated using log2 fold change values, with black points representing biological repeats.

---

## Features

- Delta Ct calculation from raw Ct values
- Delta-Delta Ct relative expression analysis
- Support for multiple reference genes
- Automatic cleaning of common genotype and treatment name variants
- Biological-repeat-level summary
- Excel output with multiple analysis sheets
- Publication-style grouped boxplots
- PNG and TIF figure export
- Italicized aphid species names in plot legends

---

## Input Requirements

### Input Excel file (`--input_file`)

The input file must contain the following columns:

| Column | Description |
|---|---|
| `genotype` | Plant genotype or accession name |
| `CT` | qRT-PCR Ct value |
| `reference_gene` | Reference gene used for normalization |
| `repeat` | Biological repeat ID |
| `treatment` | Experimental treatment |
| `type` | Either `sample` or `reference_gene` |

Example:

```text
genotype        CT      reference_gene  repeat  treatment                  type
Akashinriki     23.927  Actin           1       Clip-cage control           sample
Akashinriki     23.907  Actin           1       Clip-cage control           reference_gene
Akashinriki     23.297  Actin           1       R.padi infestation          sample
Akashinriki     21.737  Actin           1       R.padi infestation          reference_gene
```

## Installation

Install required Python packages using pip: pandas, numpy, matplotlib, and openpyxl. This can be done by running: pip install pandas numpy matplotlib openpyxl

---

## Usage

Run the script from the command line using: python qrt_pcr_delta_delta_ct_analysis.py --input_file input.xlsx --output_dir output_directory

---

## Output

The script generates an Excel file named qRT-PCR_analysis_results.xlsx and a figures directory containing grouped boxplots for each reference gene in both PNG and TIF formats.

The Excel file includes the following sheets: Raw_Data (cleaned input data), ave_reference_CT (mean reference Ct values), Sample_with_Delta_Ct (Delta Ct values), Global_Control_Delta_Ct (calibrator baseline), Per_Repeat_Values (Delta Ct, Delta-Delta Ct, and log2 fold change per biological repeat), and Summary_Stats (aggregated statistics).

---

## Calculation

Delta Ct is calculated as Ct_target minus Ct_reference. Delta-Delta Ct is calculated as Delta Ct_sample minus the mean Delta Ct of the calibrator group. Relative expression is calculated as 2^-DeltaDeltaCt. Log2 fold change is calculated as -DeltaDeltaCt.

---

## Figure Description

Each reference gene produces a grouped boxplot showing log2 fold change across genotypes and treatments. The x-axis represents genotype, and the y-axis represents log2 fold change. Boxplots show the distribution of biological replicates, and black dots indicate individual replicates. Figures use Arial font, bold axis labels, no grid lines, and italicised species names for R. padi and M. persicae.

---

## Notes

Statistical testing is not performed by default. For significance testing, it is recommended to use log2-transformed values (log2 fold change or Delta Ct-derived values) rather than raw 2^-DeltaDeltaCt values, which do not follow a normal distribution.
