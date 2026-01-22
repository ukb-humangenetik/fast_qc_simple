#!/usr/bin/env python3
import argparse

import pandas as pd

# Thresholds for short-read WGS germline QC
THRESH = {
    "q30_min_percent": 85.0,
    "mean_depth_min_x": 30.0,
    "mean_read_length_min_bp": 70.0,  # Spec: strictly > 70 bp
    "target20x_min_fraction": 0.80,
}


def make_wgs_germline_overview(input_csv: str) -> pd.DataFrame:
    df = pd.read_csv(input_csv)

    # Robust parsing for numeric columns
    num_cols = [
        "q30_percent",
        "mean_read_length_bp",
        "mean_depth_x",
        "target_fraction_ge_20x",
        "reads_r1",
        "reads_r2",
    ]
    for col in num_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Numeric checks
    df["q30_ok"] = df.get("q30_percent", pd.NA) >= THRESH["q30_min_percent"]
    df["depth_ok"] = df.get("mean_depth_x", pd.NA) >= THRESH["mean_depth_min_x"]
    df["readlen_ok"] = (
        df.get("mean_read_length_bp", pd.NA) > THRESH["mean_read_length_min_bp"]
    )
    df["target20x_ok"] = (
        df.get("target_fraction_ge_20x", pd.NA) >= THRESH["target20x_min_fraction"]
    )

    # FASTQ/readcount checks
    # - if PASS columns exist: use them
    # - otherwise: only check reads_r1 == reads_r2 (if available)
    if "readcount_match_pass" in df.columns and "fastq_structure_pass" in df.columns:
        df["fastq_ok"] = (
            (df.get("reads_r1", pd.NA) == df.get("reads_r2", pd.NA))
            & (df["readcount_match_pass"].astype(str).str.upper() == "PASS")
            & (df["fastq_structure_pass"].astype(str).str.upper() == "PASS")
        )
    else:
        if "reads_r1" in df.columns and "reads_r2" in df.columns:
            df["fastq_ok"] = df["reads_r1"] == df["reads_r2"]
        else:
            df["fastq_ok"] = pd.NA

    # Overall status (only from checks that actually exist)
    check_cols = [
        c
        for c in ["q30_ok", "depth_ok", "readlen_ok", "target20x_ok", "fastq_ok"]
        if c in df.columns
    ]
    df["qc_wgs_germline_overall"] = df[check_cols].all(axis=1)

    # Output: only relevant columns
    cols = [
        "sample_id",
        "q30_percent",
        "q30_ok",
        "mean_read_length_bp",
        "readlen_ok",
        "reads_r1",
        "reads_r2",
        "fastq_ok",
        "mean_depth_x",
        "depth_ok",
        "target_fraction_ge_20x",
        "target20x_ok",
        "qc_wgs_germline_overall",
    ]
    cols = [c for c in cols if c in df.columns]
    return df[cols].copy()


def main():
    parser = argparse.ArgumentParser(
        description="Create a QC overview for short-read WGS germline from a QC CSV."
    )
    parser.add_argument("input_csv", help="Input CSV (e.g. file.csv)")
    parser.add_argument(
        "--out-tsv",
        default="qc_overview.tsv",
        help="Output TSV (default: qc_overview.tsv)",
    )
    parser.add_argument(
        "--out-xlsx",
        default=None,
        help="Optional: output Excel (e.g. qc_overview.xlsx)",
    )
    args = parser.parse_args()

    out = make_wgs_germline_overview(args.input_csv)

    # Save outputs
    out.to_csv(args.out_tsv, sep="\t", index=False)
    if args.out_xlsx:
        out.to_excel(args.out_xlsx, index=False)

    # Also print to console (human-readable)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_colwidth", 120)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
