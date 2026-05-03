#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: run_iqtree_codon_tree_pipeline.py

Description:
    Run IQ-TREE phylogenetic analysis with checkpoint-aware resume support.

    This pipeline:
        1. Runs ModelFinder to select the best-fit substitution model
        2. Extracts the best model from IQ-TREE output
        3. Builds a maximum-likelihood phylogenetic tree
        4. Runs bootstrap analysis
        5. Saves checkpoints for interrupted runs
        6. Generates a model selection report

Input:
    - Trimmed sequence alignment in FASTA format

Output:
    - ModelFinder output files
    - Best model report
    - Maximum-likelihood tree files
    - Checkpoint files

Usage:
    python run_iqtree_codon_tree_pipeline.py \
        --alignment all_cds_trimmed.fasta \
        --output_dir output_tree \
        --iqtree iqtree3 \
        --prefix thionin \
        --threads 4 \
        --bootstrap 1000 \
        --bootstrap_type standard
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run IQ-TREE model selection and ML tree construction with checkpoint support."
    )

    parser.add_argument(
        "--alignment",
        required=True,
        help="Input trimmed alignment file in FASTA format."
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory."
    )

    parser.add_argument(
        "--iqtree",
        default="iqtree3",
        help="Path to IQ-TREE executable. Default: iqtree3"
    )

    parser.add_argument(
        "--prefix",
        default="thionin",
        help="Output prefix for tree construction. Default: thionin"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of CPU threads. Default: 4"
    )

    parser.add_argument(
        "--bootstrap",
        type=int,
        default=1000,
        help="Number of bootstrap replicates. Default: 1000"
    )

    parser.add_argument(
        "--bootstrap_type",
        choices=["standard", "ultrafast"],
        default="standard",
        help="Bootstrap type: standard uses -b, ultrafast uses -bb. Default: standard"
    )

    parser.add_argument(
        "--fallback_model",
        default="GTR+G",
        help="Fallback model if no model can be extracted. Default: GTR+G"
    )

    return parser.parse_args()


def run_command(command):
    subprocess.run(command, check=True)


def ensure_file_exists(path, label):
    if not Path(path).exists():
        sys.exit(f"Error: {label} not found: {path}")


def backup_checkpoint(source_ckp, backup_ckp):
    if Path(source_ckp).exists():
        shutil.copy2(source_ckp, backup_ckp)
        print(f"Checkpoint backed up: {source_ckp} -> {backup_ckp}")
        return True
    return False


def restore_checkpoint(backup_ckp, source_ckp):
    if Path(backup_ckp).exists():
        shutil.copy2(backup_ckp, source_ckp)
        print(f"Checkpoint restored: {backup_ckp} -> {source_ckp}")
        return True
    return False


def cleanup_checkpoint(ckp_file):
    if Path(ckp_file).exists():
        Path(ckp_file).unlink()
        print(f"Checkpoint removed: {ckp_file}")


def check_model_completion(prefix):
    iqtree_file = Path(str(prefix) + ".iqtree")

    if not iqtree_file.exists():
        return False

    if iqtree_file.stat().st_size < 1000:
        return False

    content = iqtree_file.read_text(errors="ignore")

    completion_indicators = [
        "Best-fit model:",
        "Best model according to BIC:",
        "Bayesian Information Criterion",
        "Akaike information criterion",
        "Model selection completed",
        "Total CPU time used",
    ]

    return any(indicator in content for indicator in completion_indicators)


def check_tree_completion(prefix):
    iqtree_file = Path(str(prefix) + ".iqtree")
    treefile = Path(str(prefix) + ".treefile")

    if not treefile.exists():
        return False

    if treefile.stat().st_size < 100:
        return False

    if iqtree_file.exists():
        content = iqtree_file.read_text(errors="ignore")

        tree_indicators = [
            "MAXIMUM LIKELIHOOD TREE",
            "Consensus tree written to",
            "Bootstrap consensus tree written to",
            "Total CPU time used",
        ]

        return any(indicator in content for indicator in tree_indicators)

    return True


def extract_best_model(iqtree_file):
    iqtree_file = Path(iqtree_file)

    if not iqtree_file.exists():
        raise FileNotFoundError(f"ModelFinder output not found: {iqtree_file}")

    content = iqtree_file.read_text(errors="ignore")

    patterns = [
        r"Best-fit model:\s*(\S+)",
        r"Best model according to BIC:\s*(\S+)",
        r"Best model:\s*(\S+)",
    ]

    for pattern in patterns:
        match = re.search(pattern, content)
        if match:
            return match.group(1)

    return None


def generate_model_report(iqtree_file, best_model, report_path, alignment, modeltest_dir):
    iqtree_file = Path(iqtree_file)
    report_path = Path(report_path)

    if not iqtree_file.exists():
        report_path.write_text("Error: ModelFinder .iqtree file was not found.\n")
        return report_path

    content = iqtree_file.read_text(errors="ignore")
    lines = content.splitlines()

    alignment_info = []
    model_scores = []
    headers = ""

    for line in lines:
        if "Alignment has" in line or "sequences with" in line:
            alignment_info.append(line.strip())
        elif "Model" in line and "LogL" in line and ("AIC" in line or "BIC" in line):
            headers = line.strip()
        elif re.match(r"^\s*(\S+)\s+[-+]?\d*\.?\d+\s+\d+\s+[-+]?\d*\.?\d+", line):
            model_scores.append(line.strip())

    report = []
    report.append("=" * 70)
    report.append("IQ-TREE ModelFinder Report")
    report.append("=" * 70)
    report.append(f"Generated time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append(f"Input alignment: {alignment}")
    report.append(f"ModelFinder directory: {modeltest_dir}")
    report.append("")
    report.append("Alignment information")
    report.append("-" * 50)

    if alignment_info:
        for info in alignment_info[:5]:
            report.append(info)
    else:
        report.append("No alignment summary was detected.")

    report.append("")
    report.append("Best-fit model")
    report.append("-" * 50)
    report.append(f"Best model: {best_model}")
    report.append("")
    report.append("Top model comparison")
    report.append("-" * 50)

    if headers:
        report.append(headers)

    if model_scores:
        for model_line in model_scores[:5]:
            marker = "  <-- best" if best_model and best_model in model_line else ""
            report.append(f"{model_line}{marker}")
    else:
        report.append("No model comparison table was detected.")

    report.append("")
    report.append("Generated files")
    report.append("-" * 50)

    modeltest_dir = Path(modeltest_dir)
    if modeltest_dir.exists():
        for file in sorted(modeltest_dir.iterdir()):
            if file.is_file():
                size_kb = file.stat().st_size / 1024
                report.append(f"{file.name} ({size_kb:.1f} KB)")

    report.append("")
    report.append("=" * 70)
    report.append("End of report")
    report.append("=" * 70)

    report_path.write_text("\n".join(report), encoding="utf-8")
    return report_path


def main():
    args = parse_args()

    alignment = Path(args.alignment)
    output_dir = Path(args.output_dir)

    ensure_file_exists(alignment, "Alignment file")

    output_dir.mkdir(parents=True, exist_ok=True)

    modeltest_dir = output_dir / "modeltest"
    checkpoint_dir = output_dir / "checkpoint"

    modeltest_dir.mkdir(exist_ok=True)
    checkpoint_dir.mkdir(exist_ok=True)

    model_prefix = modeltest_dir / "model"
    tree_prefix = output_dir / args.prefix

    model_iqtree = Path(str(model_prefix) + ".iqtree")
    treefile = Path(str(tree_prefix) + ".treefile")

    model_report = modeltest_dir / "best_model_report.txt"

    model_ckp_source = Path(str(model_prefix) + ".ckp.gz")
    model_ckp_backup = checkpoint_dir / "model_test.ckp.gz"

    tree_ckp_source = Path(str(tree_prefix) + ".ckp.gz")
    tree_ckp_backup = checkpoint_dir / "tree_build.ckp.gz"

    print("=" * 70)
    print("IQ-TREE phylogenetic analysis pipeline")
    print("=" * 70)
    print(f"Alignment: {alignment}")
    print(f"Output directory: {output_dir}")
    print(f"Threads: {args.threads}")
    print(f"Bootstrap: {args.bootstrap}")
    print(f"Bootstrap type: {args.bootstrap_type}")
    print("=" * 70)

    tree_completed = check_tree_completion(tree_prefix)
    model_completed = check_model_completion(model_prefix)

    if tree_completed:
        print("Tree analysis appears complete. No further action required.")
        return

    if not model_completed:
        print("\nStep 1: Running ModelFinder")

        restore_checkpoint(model_ckp_backup, model_ckp_source)

        if model_ckp_source.exists():
            command_model = [
                args.iqtree,
                "-s", str(alignment),
                "-pre", str(model_prefix),
                "-nt", str(args.threads),
                "-keep-ident",
            ]
        else:
            command_model = [
                args.iqtree,
                "-s", str(alignment),
                "-m", "MFP",
                "-nt", str(args.threads),
                "-keep-ident",
                "-pre", str(model_prefix),
            ]

        try:
            run_command(command_model)
            backup_checkpoint(model_ckp_source, model_ckp_backup)
            cleanup_checkpoint(model_ckp_source)
            model_completed = True
        except subprocess.CalledProcessError as error:
            backup_checkpoint(model_ckp_source, model_ckp_backup)
            sys.exit(f"Error: ModelFinder failed.\n{error}")

    else:
        print("\nModelFinder output already exists. Skipping model selection.")

    if not model_iqtree.exists():
        sys.exit(f"Error: ModelFinder output file not found: {model_iqtree}")

    print("\nStep 2: Extracting best-fit model")

    best_model = extract_best_model(model_iqtree)

    if best_model is None:
        best_model = args.fallback_model
        print(f"Warning: could not extract best model. Using fallback model: {best_model}")
    else:
        print(f"Best-fit model: {best_model}")

    print("\nStep 3: Generating ModelFinder report")

    generate_model_report(
        iqtree_file=model_iqtree,
        best_model=best_model,
        report_path=model_report,
        alignment=alignment,
        modeltest_dir=modeltest_dir,
    )

    print(f"Model report written to: {model_report}")

    print("\nStep 4: Running ML tree analysis")

    restore_checkpoint(tree_ckp_backup, tree_ckp_source)

    if tree_ckp_source.exists():
        command_tree = [
            args.iqtree,
            "-s", str(alignment),
            "-pre", str(tree_prefix),
            "-nt", str(args.threads),
            "-keep-ident",
        ]
    else:
        if args.bootstrap_type == "standard":
            bootstrap_flag = "-b"
        else:
            bootstrap_flag = "-bb"

        command_tree = [
            args.iqtree,
            "-s", str(alignment),
            "-m", best_model,
            bootstrap_flag, str(args.bootstrap),
            "-nt", str(args.threads),
            "-keep-ident",
            "-pre", str(tree_prefix),
        ]

    try:
        run_command(command_tree)
        backup_checkpoint(tree_ckp_source, tree_ckp_backup)
        cleanup_checkpoint(tree_ckp_source)
    except subprocess.CalledProcessError as error:
        backup_checkpoint(tree_ckp_source, tree_ckp_backup)
        sys.exit(f"Error: IQ-TREE tree construction failed.\n{error}")

    print("\n" + "=" * 70)
    print("Pipeline completed successfully.")
    print("=" * 70)
    print(f"Best-fit model: {best_model}")
    print(f"Tree file: {treefile}")
    print(f"Model report: {model_report}")
    print(f"Checkpoint directory: {checkpoint_dir}")


if __name__ == "__main__":
    main()
