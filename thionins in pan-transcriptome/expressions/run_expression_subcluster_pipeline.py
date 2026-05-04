#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: run_expression_subcluster_pipeline.py

Description:
    Extract expression values for selected transcript IDs, calculate condition-level
    average expression values, merge expression records according to subcluster
    definitions, and generate a combined expression matrix.

Workflow:
    1. Normalize old transcript IDs and generate old_id -> new_id mapping
    2. Extract expression values from genotype-specific expression matrices
    3. Calculate average expression values for Ca, Co, In, Ro and Sh conditions
    4. Merge expression records according to subcluster membership
    5. Generate a processing report
    6. Combine all genotype-level subcluster expression files

Expected input files:
    ids_file:
        A text file containing one old transcript ID per line.

    expression_dir:
        Genotype-specific expression matrices in CSV format.
        Example:
            Genes_TPM_Akashinriki.csv
            Genes_TPM_B1K-04-12.csv

    subcluster_dir:
        Subcluster member files generated from the cluster-matching pipeline.
        Example:
            Akashinriki_cluster_1_subcluster_1.txt

Usage:
    python run_expression_subcluster_pipeline.py \
        --ids_file input/8_all_id.txt \
        --expression_dir input/all_expression \
        --subcluster_dir input/subcluster_transcripts \
        --output_dir output_expression
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


CONDITIONS = ["Ca", "Co", "In", "Ro", "Sh"]


OUTPUT_SUBDIRS = {
    "expression_ext": "2_expression_ext",
    "expression_ave": "3_expression_ext_ave",
    "expression_subcluster": "4_expression_ext_ave_subcluster",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract expression values, calculate condition-level averages, "
            "and merge expression records according to subcluster definitions."
        )
    )

    parser.add_argument(
        "--ids_file",
        required=True,
        help="Text file containing old transcript IDs, one ID per line."
    )

    parser.add_argument(
        "--expression_dir",
        required=True,
        help="Directory containing genotype-specific expression CSV files."
    )

    parser.add_argument(
        "--subcluster_dir",
        required=True,
        help="Directory containing subcluster member TXT files."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory."
    )

    return parser.parse_args()


def prepare_output_paths(output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    subdirs = {}

    for key, subdir_name in OUTPUT_SUBDIRS.items():
        path = output_dir / subdir_name
        path.mkdir(parents=True, exist_ok=True)
        subdirs[key] = path

    paths = {
        "base_out": output_dir,
        "id_map_file": output_dir / "1_old_new_id.txt",
        "expression_ext_dir": subdirs["expression_ext"],
        "expression_ave_dir": subdirs["expression_ave"],
        "expression_subcluster_dir": subdirs["expression_subcluster"],
        "report_file": output_dir / "expression_pipeline_report.txt",
        "combined_expression_file": output_dir / "combined_expression_ext_ave_subcluster.xlsx",
    }

    return paths


def normalize_id(old_id):
    """
    Normalize transcript ID by:
        1. Removing trailing '_pseudogene'
        2. Removing the final isoform suffix '.number'

    Examples:
        Akashinriki_chr1HG00001.1
            -> Akashinriki_chr1HG00001

        HOR10350_chr6HG26769.2_pseudogene
            -> HOR10350_chr6HG26769

        HOR8148_contig_corrected_v1_8420G38821.1_pseudogene
            -> HOR8148_contig_corrected_v1_8420G38821

        Morex_CAJHDD010000030.1G00264.1
            -> Morex_CAJHDD010000030.1G00264
    """
    if not isinstance(old_id, str):
        return old_id

    normalized = re.sub(r"_pseudogene$", "", old_id)
    normalized = re.sub(r"\.\d+$", "", normalized)

    return normalized


def get_variety_from_expr_filename(csv_path):
    """
    Infer genotype/variety name from expression matrix filename.

    Examples:
        Genes_TPM_Akashinriki.csv -> Akashinriki
        Genes_TPM_B1K-04-12.csv   -> B1K-04-12
    """
    stem = csv_path.stem
    parts = stem.split("_")

    return parts[-1]


def parse_subcluster_filename(path):
    """
    Parse subcluster filename.

    Example:
        Akashinriki_cluster_1_subcluster_1.txt
            -> ("Akashinriki", "1", "1")
    """
    name = path.stem
    match = re.match(r"^(.+)_cluster_(\d+)_subcluster_(\d+)$", name)

    if not match:
        return None, None, None

    return match.group(1), match.group(2), match.group(3)


def make_id_mapping(ids_file, id_map_file):
    """
    Read old transcript IDs and generate an old_id -> new_id mapping table.

    One new_id may correspond to multiple old_id values.
    All mappings are retained without deduplication.
    """
    print("Step 1: Generate old_id -> new_id mapping ...")

    rows = []
    all_old_ids = []

    with open(ids_file, encoding="utf-8") as f:
        for line in f:
            old_id = line.strip()

            if not old_id:
                continue

            new_id = normalize_id(old_id)
            rows.append((old_id, new_id))
            all_old_ids.append(old_id)

    id_map_df = pd.DataFrame(rows, columns=["old_id", "new_id"])
    id_map_df.to_csv(id_map_file, sep="\t", index=False)

    print(f"  Number of old IDs: {len(all_old_ids)}")
    print(f"  Mapping file written: {id_map_file}")

    return id_map_df, set(all_old_ids)


def process_expression_files(
    id_map_df,
    processed_old_ids,
    expression_dir,
    expression_ext_dir
):
    """
    Process genotype-specific expression matrices.

    For each CSV file:
        1. Identify the gene_id column
        2. Infer genotype/variety name from filename
        3. Remove '<variety>.' prefix from expression column names
        4. Merge expression matrix with id_map_df by:
               gene_id == new_id
        5. Add notes column containing old_id
        6. Export extended expression file
    """
    print("\nStep 2: Extract expression values by new_id ...")

    expression_dir = Path(expression_dir)
    expression_files = sorted(expression_dir.glob("*.csv"))

    if not expression_files:
        print(f"  Warning: no CSV expression files found in {expression_dir}")

    variety_to_extfile = {}

    for csv_path in expression_files:
        print(f"  Processing expression file: {csv_path.name}")

        variety = get_variety_from_expr_filename(csv_path)

        df = pd.read_csv(csv_path)

        if "gene_id" in df.columns:
            id_col = "gene_id"
        else:
            id_col = df.columns[0]
            df = df.rename(columns={id_col: "gene_id"})
            id_col = "gene_id"

        prefix = variety + "."
        new_columns = {}

        for col in df.columns:
            if col == id_col:
                continue

            if str(col).startswith(prefix):
                new_columns[col] = str(col)[len(prefix):]

        df = df.rename(columns=new_columns)

        merged = df.merge(
            id_map_df,
            left_on=id_col,
            right_on="new_id",
            how="inner"
        )

        if merged.empty:
            print("    Warning: no overlap with ID mapping. Skipped.")
            continue

        merged["notes"] = merged["old_id"]
        merged = merged.drop(columns=["new_id", "old_id"])

        processed_old_ids.update(
            merged["notes"].dropna().astype(str).tolist()
        )

        columns = list(merged.columns)
        columns.remove("notes")
        columns.append("notes")
        merged = merged[columns]

        out_path = expression_ext_dir / f"{variety}_expression_ext.xlsx"
        merged.to_excel(out_path, index=False)
        variety_to_extfile[variety] = out_path

        print(f"    Extracted rows: {merged.shape[0]}")
        print(f"    Output: {out_path}")

    return variety_to_extfile


def compute_condition_means(variety_to_extfile, expression_ave_dir):
    """
    Calculate average expression values for Ca, Co, In, Ro and Sh.

    Input columns are expected to include replicate-style names:
        Ca1, Ca2, ...
        Co1, Co2, ...
        In1, In2, ...
        Ro1, Ro2, ...
        Sh1, Sh2, ...

    Output columns:
        gene_id, Ca_ave, Co_ave, In_ave, Ro_ave, Sh_ave, notes, number
    """
    print("\nStep 3: Calculate condition-level average expression ...")

    variety_to_avefile = {}

    for variety, ext_path in variety_to_extfile.items():
        print(f"  Processing variety: {variety}")

        df = pd.read_excel(ext_path, engine="openpyxl")

        if "gene_id" not in df.columns:
            df = df.rename(columns={df.columns[0]: "gene_id"})

        if "notes" not in df.columns:
            df["notes"] = pd.NA

        for condition in CONDITIONS:
            condition_columns = [
                col for col in df.columns
                if re.fullmatch(rf"{condition}\d+", str(col))
            ]

            if condition_columns:
                df[f"{condition}_ave"] = df[condition_columns].mean(
                    axis=1,
                    skipna=True
                )
            else:
                df[f"{condition}_ave"] = np.nan

        df["number"] = 1

        output_columns = (
            ["gene_id"]
            + [f"{condition}_ave" for condition in CONDITIONS]
            + ["notes", "number"]
        )

        df_out = df[output_columns].copy()

        out_path = expression_ave_dir / f"{variety}_expression_ext_ave.xlsx"
        df_out.to_excel(out_path, index=False)
        variety_to_avefile[variety] = out_path

        print(f"    Rows: {df_out.shape[0]}")
        print(f"    Output: {out_path}")

    return variety_to_avefile


def load_subcluster_files(subcluster_dir):
    """
    Load subcluster definition files.

    Expected filename format:
        Species_cluster_1_subcluster_1.txt

    Each TXT file should contain one old transcript ID per line.
    """
    print("\nStep 4: Load subcluster definitions ...")

    subcluster_dir = Path(subcluster_dir)
    sub_info_by_species = defaultdict(list)
    global_sub_old_ids = set()

    txt_files = sorted(subcluster_dir.glob("*.txt"))

    if not txt_files:
        print(f"  Warning: no TXT subcluster files found in {subcluster_dir}")

    for path in txt_files:
        species, cluster_num, subcluster_num = parse_subcluster_filename(path)

        if species is None:
            print(f"  Warning: cannot parse filename {path.name}. Skipped.")
            continue

        with open(path, encoding="utf-8") as f:
            members = [line.strip() for line in f if line.strip()]

        sub_info_by_species[species].append(
            (cluster_num, subcluster_num, members)
        )
        global_sub_old_ids.update(members)

    print(f"  Number of species with subclusters: {len(sub_info_by_species)}")
    print(f"  Number of old IDs in subclusters: {len(global_sub_old_ids)}")

    return sub_info_by_species, global_sub_old_ids


def merge_by_subcluster(
    variety_to_avefile,
    sub_info_by_species,
    expression_subcluster_dir
):
    """
    Merge expression records according to subcluster membership.

    For each variety:
        - Read average expression file
        - Match subcluster members using old_id values stored in notes
        - Sum Ca_ave, Co_ave, In_ave, Ro_ave and Sh_ave across members
        - Merge notes using ';'
        - Set number to the number of merged old IDs
        - Name merged gene_id as:
              <Species>_C<cluster_num>_S<subcluster_num>
        - Keep non-subcluster rows unchanged
    """
    print("\nStep 4: Merge expression records by subcluster ...")

    for variety, ave_path in variety_to_avefile.items():
        print(f"  Processing variety: {variety}")

        df = pd.read_excel(ave_path, engine="openpyxl")

        if "gene_id" not in df.columns:
            df = df.rename(columns={df.columns[0]: "gene_id"})

        if "notes" not in df.columns:
            df["notes"] = pd.NA

        if "number" not in df.columns:
            df["number"] = 1

        sub_list = sub_info_by_species.get(variety, [])

        out_path = (
            expression_subcluster_dir
            / f"{variety}_expression_ext_ave_subcluster.xlsx"
        )

        if not sub_list:
            df.to_excel(out_path, index=False)
            print(f"    No subcluster definition. Original file copied to: {out_path}")
            continue

        oldid_to_rows = defaultdict(list)

        for idx, row in df.iterrows():
            notes_value = str(row["notes"]) if pd.notna(row["notes"]) else ""

            for old_id in notes_value.split(";"):
                old_id = old_id.strip()

                if old_id:
                    oldid_to_rows[old_id].append(idx)

        species_sub_old_ids = set()

        for _, _, members in sub_list:
            species_sub_old_ids.update(members)

        merged_rows = []
        used_indices = set()

        for cluster_num, subcluster_num, members in sub_list:
            row_indices = set()
            merged_old_ids = set()

            for old_id in members:
                rows_for_old_id = oldid_to_rows.get(old_id, [])

                if rows_for_old_id:
                    row_indices.update(rows_for_old_id)
                    merged_old_ids.add(old_id)

            if not row_indices:
                continue

            used_indices.update(row_indices)

            row_data = {}
            row_data["gene_id"] = f"{variety}_C{cluster_num}_S{subcluster_num}"

            for condition in CONDITIONS:
                col = f"{condition}_ave"

                if col in df.columns:
                    row_data[col] = df.loc[list(row_indices), col].sum(
                        skipna=True
                    )
                else:
                    row_data[col] = np.nan

            row_data["notes"] = ";".join(sorted(merged_old_ids))
            row_data["number"] = len(merged_old_ids)

            merged_rows.append(row_data)

        df_merged = pd.DataFrame(merged_rows)

        df["notes"] = df["notes"].astype(str)
        mask_unassigned = ~df["notes"].isin(species_sub_old_ids)
        df_unassigned = df[mask_unassigned].copy()
        df_unassigned["number"] = 1

        condition_columns = [f"{condition}_ave" for condition in CONDITIONS]
        all_columns = ["gene_id"] + condition_columns + ["notes", "number"]

        frames = []

        if not df_merged.empty:
            frames.append(df_merged[all_columns])

        if not df_unassigned.empty:
            frames.append(df_unassigned[all_columns])

        if frames:
            final_df = pd.concat(frames, ignore_index=True)
        else:
            final_df = pd.DataFrame(columns=all_columns)

        final_df.to_excel(out_path, index=False)

        print(f"    Rows after merging: {final_df.shape[0]}")
        print(f"    Output: {out_path}")


def write_report(all_old_ids, processed_old_ids, report_file):
    """
    Write processing report.
    """
    total = len(all_old_ids)
    processed = len(processed_old_ids)
    unprocessed_ids = sorted(all_old_ids - processed_old_ids)
    unprocessed = len(unprocessed_ids)

    lines = []
    lines.append(f"Total input sequence IDs: {total}")
    lines.append(f"Processed sequence IDs: {processed}")
    lines.append(f"Unprocessed sequence IDs: {unprocessed}")
    lines.append("")
    lines.append("Unprocessed sequence ID list:")

    if unprocessed_ids:
        lines.extend(unprocessed_ids)
    else:
        lines.append("(None)")

    text = "\n".join(lines)

    print("\n========== Processing report ==========")
    print(text)

    with open(report_file, "w", encoding="utf-8") as f:
        f.write(text)

    print(f"\nReport written: {report_file}")


def combine_all_subcluster_results(
    expression_subcluster_dir,
    combined_expression_file
):
    """
    Combine all genotype-level subcluster expression files into one Excel file.
    """
    print("\nStep 6: Combine all subcluster expression results ...")

    files = sorted(
        Path(expression_subcluster_dir).glob(
            "*_expression_ext_ave_subcluster.xlsx"
        )
    )

    if not files:
        print("  No subcluster expression files found. Skipped.")
        return None

    dfs = []

    for file_path in files:
        df = pd.read_excel(file_path, engine="openpyxl")
        variety = file_path.stem.split("_expression_ext_ave_subcluster")[0]
        df["variety"] = variety
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_excel(combined_expression_file, index=False)

    print(f"  Combined rows: {combined.shape[0]}")
    print(f"  Output: {combined_expression_file}")

    return combined_expression_file


def main():
    args = parse_args()
    paths = prepare_output_paths(args.output_dir)

    processed_old_ids = set()

    id_map_df, all_old_ids = make_id_mapping(
        ids_file=args.ids_file,
        id_map_file=paths["id_map_file"]
    )

    variety_to_extfile = process_expression_files(
        id_map_df=id_map_df,
        processed_old_ids=processed_old_ids,
        expression_dir=args.expression_dir,
        expression_ext_dir=paths["expression_ext_dir"]
    )

    variety_to_avefile = compute_condition_means(
        variety_to_extfile=variety_to_extfile,
        expression_ave_dir=paths["expression_ave_dir"]
    )

    sub_info_by_species, _ = load_subcluster_files(
        subcluster_dir=args.subcluster_dir
    )

    merge_by_subcluster(
        variety_to_avefile=variety_to_avefile,
        sub_info_by_species=sub_info_by_species,
        expression_subcluster_dir=paths["expression_subcluster_dir"]
    )

    write_report(
        all_old_ids=all_old_ids,
        processed_old_ids=processed_old_ids,
        report_file=paths["report_file"]
    )

    combine_all_subcluster_results(
        expression_subcluster_dir=paths["expression_subcluster_dir"],
        combined_expression_file=paths["combined_expression_file"]
    )

    print("\n=== All steps completed. ===\n")


if __name__ == "__main__":
    main()
