#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: qrt_pcr__analysis.py

Description:
    Analyse qRT-PCR data using the △△Ct method and generate
    publication-style grouped boxplots for each reference gene.

Input:
    1. qRT-PCR Ct data in Excel format

       Required columns:
           genotype
           CT
           reference_gene
           repeat
           treatment
           type

       The 'type' column should contain:
           sample
           reference_gene

Output:
    1. qRT-PCR_analysis_results.xlsx
       Excel workbook containing cleaned data, reference Ct values,
       Delta Ct, Delta-Delta Ct, relative expression, log2 fold change,
       and summary statistics.

    2. figures/
       Grouped boxplots for each reference gene in PNG and TIF format.

Calculation:
    Delta Ct = Ct_target - Ct_reference

    DeltaDelta Ct = Delta Ct_sample - mean Delta Ct_calibrator

    Relative expression = 2^-DeltaDeltaCt

    log2 fold change = -DeltaDeltaCt

Default calibrator:
    Akashinriki untreated control

Usage:
    python qrt_pcr_delta_delta_ct_analysis.py \
        --input_file input.xlsx \
        --output_dir output
"""

import argparse
import os
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")


# ==============================
# Global plotting parameters
# ==============================

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False
mpl.rcParams["mathtext.default"] = "regular"


# ==============================
# Default analysis parameters
# ==============================

TREATMENT_ORDER = [
    "Untreated control",
    "Clip-cage control",
    "R.padi infestation",
    "M.persicae infestation",
]

GENOTYPE_ORDER = [
    "Akashinriki",
    "HOR10350",
    "HOR21599",
    "Morex",
]

REFERENCE_GENES = [
    "Actin",
    "Ubiquitin",
]

TREATMENT_COLORS = {
    "Untreated control": "#D9D9D9",
    "Clip-cage control": "#CFE8F3",
    "R.padi infestation": "#F6D6AD",
    "M.persicae infestation": "#D9E8C8",
}

GLOBAL_CALIBRATOR_GENOTYPE = "Akashinriki"
GLOBAL_CALIBRATOR_TREATMENT = "Untreated control"


# ==============================
# Argument parser
# ==============================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Analyse qRT-PCR data using the Delta-Delta Ct method "
            "and generate grouped boxplots."
        )
    )

    parser.add_argument(
        "--input_file",
        required=True,
        help="qRT-PCR Ct data in Excel format."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory."
    )

    return parser.parse_args()


# ==============================
# Data loading and cleaning
# ==============================

def load_input_data(input_file):
    print("Loading input data...")
    df = pd.read_excel(input_file)
    print(f"Loaded {df.shape[0]} rows and {df.shape[1]} columns.")
    return df


def clean_input_data(df):
    print("Cleaning input data...")

    required_columns = [
        "genotype",
        "CT",
        "reference_gene",
        "repeat",
        "treatment",
        "type",
    ]

    missing_columns = [
        col for col in required_columns
        if col not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Input file is missing required columns: {missing_columns}"
        )

    for col in ["genotype", "reference_gene", "treatment", "type"]:
        df[col] = df[col].astype(str).str.strip()

    df["CT"] = pd.to_numeric(df["CT"], errors="coerce")

    df["genotype"] = df["genotype"].replace({
        "Akainriki": "Akashinriki",
        "Akaishinriki": "Akashinriki",
        "HOR1599": "HOR21599",
    })

    df["reference_gene"] = df["reference_gene"].replace({
        "actin": "Actin",
        "ubiquitin": "Ubiquitin",
    })

    df["treatment"] = df["treatment"].replace({
        "No treatment": "Untreated control",
        "Empty cage": "Clip-cage control",
        "Cage control": "Clip-cage control",
        "R.padi": "R.padi infestation",
        "M.persicae": "M.persicae infestation",
    })

    before_drop = df.shape[0]
    df = df.dropna(subset=["CT"]).copy()
    after_drop = df.shape[0]

    if before_drop != after_drop:
        print(f"Removed {before_drop - after_drop} rows with missing Ct values.")

    df = df[df["treatment"].isin(TREATMENT_ORDER)].copy()
    df = df[df["genotype"].isin(GENOTYPE_ORDER)].copy()
    df = df[df["reference_gene"].isin(REFERENCE_GENES)].copy()

    print("Detected genotypes:", list(df["genotype"].dropna().unique()))
    print("Detected treatments:", list(df["treatment"].dropna().unique()))
    print("Detected reference genes:", list(df["reference_gene"].dropna().unique()))
    print("Detected data types:", list(df["type"].dropna().unique()))

    return df


# ==============================
# qRT-PCR calculation
# ==============================

def calculate_qpcr_values(df):
    print("Calculating mean reference Ct values...")

    reference_df = df[df["type"] == "reference_gene"].copy()

    ave_ref_ct = (
        reference_df
        .groupby(["genotype", "reference_gene", "repeat", "treatment"])["CT"]
        .mean()
        .reset_index()
        .rename(columns={"CT": "ave_reference_CT"})
    )

    print("Calculating Delta Ct values...")

    sample_df = df[df["type"] == "sample"].copy()

    sample_df = pd.merge(
        sample_df,
        ave_ref_ct,
        on=["genotype", "reference_gene", "repeat", "treatment"],
        how="left",
    )

    sample_df["Delta_Ct"] = (
        sample_df["CT"] - sample_df["ave_reference_CT"]
    )

    print("Averaging technical replicates within each biological repeat...")

    per_repeat_values = (
        sample_df
        .groupby(["genotype", "reference_gene", "repeat", "treatment"])[
            ["CT", "ave_reference_CT", "Delta_Ct"]
        ]
        .mean()
        .reset_index()
        .rename(columns={
            "CT": "ave_sample_CT",
            "Delta_Ct": "ave_Delta_Ct",
        })
    )

    print("Calculating global calibrator Delta Ct values...")

    global_control = (
        per_repeat_values[
            (per_repeat_values["genotype"] == GLOBAL_CALIBRATOR_GENOTYPE) &
            (per_repeat_values["treatment"] == GLOBAL_CALIBRATOR_TREATMENT)
        ]
        .groupby("reference_gene")["ave_Delta_Ct"]
        .mean()
        .reset_index()
        .rename(columns={
            "ave_Delta_Ct": "global_control_Delta_Ct_ave"
        })
    )

    per_repeat_values = pd.merge(
        per_repeat_values,
        global_control,
        on="reference_gene",
        how="left",
    )

    print("Calculating Delta-Delta Ct, relative expression, and log2 fold change...")

    per_repeat_values["Delta_Delta_Ct"] = (
        per_repeat_values["ave_Delta_Ct"] -
        per_repeat_values["global_control_Delta_Ct_ave"]
    )

    per_repeat_values["2^-Delta_Delta_Ct"] = (
        2 ** (-per_repeat_values["Delta_Delta_Ct"])
    )

    per_repeat_values["log2FC"] = (
        -per_repeat_values["Delta_Delta_Ct"]
    )

    print("Calculating summary statistics...")

    summary_stats = (
        per_repeat_values
        .groupby(["genotype", "reference_gene", "treatment"])[
            ["2^-Delta_Delta_Ct", "log2FC"]
        ]
        .agg(["mean", "median", "sem", "std", "count"])
        .reset_index()
    )

    summary_stats.columns = [
        "genotype",
        "reference_gene",
        "treatment",
        "mean_2^-Delta_Delta_Ct",
        "median_2^-Delta_Delta_Ct",
        "sem_2^-Delta_Delta_Ct",
        "std_2^-Delta_Delta_Ct",
        "n_2^-Delta_Delta_Ct",
        "mean_log2FC",
        "median_log2FC",
        "sem_log2FC",
        "std_log2FC",
        "n_log2FC",
    ]

    return (
        ave_ref_ct,
        sample_df,
        global_control,
        per_repeat_values,
        summary_stats,
    )


# ==============================
# Plotting
# ==============================

def format_treatment_label(treatment):
    if treatment == "R.padi infestation":
        return r"$\it{R.\ padi}$ infestation"

    if treatment == "M.persicae infestation":
        return r"$\it{M.\ persicae}$ infestation"

    return treatment


def draw_grouped_boxplot(ax, plot_df):
    group_spacing = 5.5
    offsets = [-1.2, -0.4, 0.4, 1.2]
    box_width = 0.55

    positions = []
    data_for_boxes = []
    colors_for_boxes = []
    group_centers = []

    for i, genotype in enumerate(GENOTYPE_ORDER):
        group_center = i * group_spacing
        group_centers.append(group_center)

        for j, treatment in enumerate(TREATMENT_ORDER):
            pos = group_center + offsets[j]

            temp = plot_df[
                (plot_df["genotype"] == genotype) &
                (plot_df["treatment"] == treatment)
            ]

            values = temp["log2FC"].dropna().values

            if len(values) == 0:
                values = [np.nan]

            positions.append(pos)
            data_for_boxes.append(values)
            colors_for_boxes.append(TREATMENT_COLORS[treatment])

    bp = ax.boxplot(
        data_for_boxes,
        positions=positions,
        widths=box_width,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.8),
        boxprops=dict(edgecolor="black", linewidth=1.3),
        whiskerprops=dict(color="black", linewidth=1.3),
        capprops=dict(color="black", linewidth=1.3),
    )

    for patch, color in zip(bp["boxes"], colors_for_boxes):
        patch.set_facecolor(color)
        patch.set_alpha(0.9)

    for i, genotype in enumerate(GENOTYPE_ORDER):
        group_center = i * group_spacing

        for j, treatment in enumerate(TREATMENT_ORDER):
            pos = group_center + offsets[j]

            temp = plot_df[
                (plot_df["genotype"] == genotype) &
                (plot_df["treatment"] == treatment)
            ]

            yvals = temp["log2FC"].dropna().values

            if len(yvals) > 0:
                jitter = (
                    np.linspace(-0.06, 0.06, len(yvals))
                    if len(yvals) > 1
                    else [0]
                )

                xvals = [pos + jj for jj in jitter]

                ax.plot(
                    xvals,
                    yvals,
                    "o",
                    markersize=5.5,
                    markerfacecolor="black",
                    markeredgecolor="black",
                    alpha=0.85,
                )

    ax.set_xticks(group_centers)
    ax.set_xticklabels(GENOTYPE_ORDER, fontsize=12, rotation=0)

    ax.set_xlabel("Genotype", fontsize=14, fontweight="bold")
    ax.set_ylabel(r"log$_2$ fold change", fontsize=14, fontweight="bold")

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.5)
        spine.set_color("black")

    ax.tick_params(
        axis="both",
        which="both",
        width=1.2,
        length=5,
        labelsize=12,
    )

    ax.grid(False)

    legend_handles = [
        Patch(
            facecolor=TREATMENT_COLORS[t],
            edgecolor="black",
            label=format_treatment_label(t),
        )
        for t in TREATMENT_ORDER
    ]

    ax.legend(
        handles=legend_handles,
        title=None,
        frameon=False,
        fontsize=11,
    )


# ==============================
# Output
# ==============================

def save_excel_results(
    output_excel,
    cleaned_df,
    ave_ref_ct,
    sample_df,
    global_control,
    per_repeat_values,
    summary_stats,
):
    print("Saving Excel results...")

    with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
        cleaned_df.to_excel(
            writer,
            sheet_name="Raw_Data",
            index=False,
        )

        ave_ref_ct.to_excel(
            writer,
            sheet_name="ave_reference_CT",
            index=False,
        )

        sample_df.to_excel(
            writer,
            sheet_name="Sample_with_Delta_Ct",
            index=False,
        )

        global_control.to_excel(
            writer,
            sheet_name="Global_Control_Delta_Ct",
            index=False,
        )

        per_repeat_values.to_excel(
            writer,
            sheet_name="Per_Repeat_Values",
            index=False,
        )

        summary_stats.to_excel(
            writer,
            sheet_name="Summary_Stats",
            index=False,
        )

    print(f"Excel results saved to: {output_excel}")


def save_reference_gene_figures(per_repeat_values, figure_dir):
    print("Generating figures...")

    for ref_gene in REFERENCE_GENES:
        sub_df = per_repeat_values[
            per_repeat_values["reference_gene"] == ref_gene
        ].copy()

        if sub_df.empty:
            print(f"Warning: no data found for {ref_gene}. Skipping.")
            continue

        fig, ax = plt.subplots(figsize=(14, 7))

        draw_grouped_boxplot(ax, sub_df)

        safe_ref = ref_gene.replace(" ", "_")

        png_path = figure_dir / f"{safe_ref}_grouped_boxplot_log2FC.png"
        tif_path = figure_dir / f"{safe_ref}_grouped_boxplot_log2FC.tif"

        plt.tight_layout()
        plt.savefig(png_path, dpi=300, bbox_inches="tight")
        plt.savefig(tif_path, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"{ref_gene} figure saved:")
        print(png_path)
        print(tif_path)


# ==============================
# Main workflow
# ==============================

def main():
    start_time = None

    args = parse_args()

    input_file = Path(args.input_file)
    output_dir = Path(args.output_dir)
    figure_dir = output_dir / "figures"

    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    raw_df = load_input_data(input_file)
    cleaned_df = clean_input_data(raw_df)

    (
        ave_ref_ct,
        sample_df,
        global_control,
        per_repeat_values,
        summary_stats,
    ) = calculate_qpcr_values(cleaned_df)

    output_excel = output_dir / "qRT-PCR_analysis_results.xlsx"

    save_excel_results(
        output_excel=output_excel,
        cleaned_df=cleaned_df,
        ave_ref_ct=ave_ref_ct,
        sample_df=sample_df,
        global_control=global_control,
        per_repeat_values=per_repeat_values,
        summary_stats=summary_stats,
    )

    save_reference_gene_figures(
        per_repeat_values=per_repeat_values,
        figure_dir=figure_dir,
    )

    print("Analysis completed.")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
