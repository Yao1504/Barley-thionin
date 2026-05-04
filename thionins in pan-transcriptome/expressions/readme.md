# Expression Processing and Subcluster-Level Quantification

## Overview

This pipeline retrieves and summarizes transcript expression levels for candidate thionin transcripts identified from the barley pan-transcriptome. It extracts expression values using candidate transcript identifiers, calculates tissue-level mean expression values across biological replicates, merges expression values for transcripts belonging to the same subcluster, and generates genotype-level and combined expression tables for downstream analysis.

## Method Description

Transcript expression levels were retrieved based on the candidate thionin transcript identifiers obtained from the pan-transcriptome. The pan-transcriptome expression dataset contains expression profiles for five tissues: caryopsis, coleoptile, inflorescence, root, and shoot. In this pipeline, these tissues correspond to the expression column prefixes Ca, Co, In, Ro, and Sh, respectively. For each tissue, biological replicate columns are identified automatically using column names such as Ca1, Ca2, Ca3, Co1, Co2, Co3, In1, In2, In3, Ro1, Ro2, Ro3, and Sh1, Sh2, Sh3. The mean expression value across all available biological replicates is calculated for each tissue and genotype. Therefore, genotypes with only two available root replicates, such as HOR10350, HOR7552, and HOR8148, are handled automatically because the mean is calculated from all existing Ro replicate columns. As transcripts within the same subcluster share 100% nucleotide sequence identity, expression values cannot be unambiguously assigned to individual subcluster members. Therefore, expression values of all transcripts belonging to the same subcluster are summed to represent the total expression level of that subcluster, and the number of transcripts within each subcluster is recorded. Transcripts not assigned to any subcluster are retained as individual entries with number = 1. The resulting expression tables can be used to calculate total expression values for thionin clusters within each genotype.

## Workflow

1. Normalize transcript identifiers by removing "_pseudogene" suffixes and final isoform suffixes such as ".1" or ".2".
2. Generate an old_id to new_id mapping table.
3. Extract candidate thionin transcript expression values from genotype-specific expression matrices.
4. Calculate mean expression values for caryopsis, coleoptile, inflorescence, root, and shoot.
5. Load subcluster membership files.
6. Sum expression values for transcripts belonging to the same subcluster.
7. Record the number of transcripts represented by each subcluster.
8. Retain non-subcluster transcripts as individual records.
9. Generate a processing report.
10. Combine all genotype-level subcluster expression files into one final table.

## Input Files

### Candidate transcript ID file

A plain text file containing one candidate thionin transcript ID per line, for example:

Akashinriki_chr1HG00001.1  
HOR10350_chr6HG26769.2_pseudogene  
Morex_CAJHDD010000030.1G00264.1  

### Expression matrix directory

A directory containing genotype-specific expression matrices in CSV format, for example:

Genes_TPM_Akashinriki.csv  
Genes_TPM_B1K-04-12.csv  
Genes_TPM_HOR10350.csv  

The genotype name is inferred from the filename as the last underscore-separated field (e.g. Akashinriki from Genes_TPM_Akashinriki.csv).

Each expression file should contain a gene_id column. If not present, the first column is treated as gene_id.

Expression columns should follow naming patterns such as:

Ca1, Ca2, Ca3  
Co1, Co2, Co3  
In1, In2, In3  
Ro1, Ro2, Ro3  
Sh1, Sh2, Sh3  

If column names include genotype prefixes (e.g. Akashinriki.Ca1), the prefix is automatically removed.

### Subcluster directory

A directory containing subcluster membership files. Each file name follows the format:

<genotype>_cluster_<cluster_number>_subcluster_<subcluster_number>.txt

For example:

Akashinriki_cluster_1_subcluster_1.txt

Each file contains one transcript ID per line. These IDs must correspond to the original transcript IDs (old_id).

## Installation

Install required Python packages:

pip install -r requirements.txt

Example requirements.txt:

numpy  
pandas  
openpyxl  

## Usage

Run the pipeline using:

python run_expression_subcluster_pipeline.py \
    --ids_file input/8_all_id.txt \
    --expression_dir input/all_expression \
    --subcluster_dir input/subcluster_transcripts \
    --output_dir output_expression

## Output Structure

output_expression/  
1_old_new_id.txt  
2_expression_ext/  
3_expression_ext_ave/  
4_expression_ext_ave_subcluster/  
expression_pipeline_report.txt  
combined_expression_ext_ave_subcluster.xlsx  

## Output Description

1_old_new_id.txt  
Mapping table between original transcript IDs (old_id) and normalized IDs (new_id).

2_expression_ext/  
Contains extracted expression tables for each genotype.

3_expression_ext_ave/  
Contains average expression values for each tissue per genotype. Columns include gene_id, Ca_ave, Co_ave, In_ave, Ro_ave, Sh_ave, notes, and number.

4_expression_ext_ave_subcluster/  
Contains subcluster-merged expression tables. Expression values are summed within each subcluster, transcript IDs are merged in the notes column, and number indicates the number of transcripts represented.

expression_pipeline_report.txt  
Summary report listing total, processed, and unprocessed transcript IDs.

combined_expression_ext_ave_subcluster.xlsx  
Combined expression table across all genotypes, including a variety column.

## Tissue Abbreviations

Ca = Caryopsis  
Co = Coleoptile  
In = Inflorescence  
Ro = Root  
Sh = Shoot  

## Notes

Mean expression values are calculated using all available biological replicates. Missing replicates are allowed. Subcluster-level expression values are calculated as the sum of all member transcripts. Non-subcluster transcripts are retained individually. Total cluster-level expression can be obtained by summing subcluster-level values within each genotype.
