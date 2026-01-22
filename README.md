# fast_qc_simple

Simple single-sample QC for **WGS / short-reads / germline** (SE context).

This repository provides a compact script `fast_qc_simple.py` that performs a minimal set of QC checks for a single sample:

- FASTQ R1/R2 checks: read counts, FASTQ structure, mean read length, optional per-base Q30% calculation
- CRAM + BED: mean coverage (mean depth) and fraction of targeted bases >= threshold (default 20x)
- Optional: DRAGEN/Mapper mapping metrics for Q30% extraction

It writes a single-row CSV report and prints a short summary to stdout. When run with `--verbose` it provides running feedback and writes a plain-text logfile named `<sample_id>.log`.

Author / version

- Author: Stephan Rosecker
- Script version: 1.1.2

---

## Highlights (what's new / important)

- `--verbose` mode:
  
  - Prints progress and status messages to the console (colored if `rich` is installed).
  - Creates a plain-text logfile named `<sample_id>.log` in the current working directory and writes the same timestamped messages there.
  - Logs include: FASTQ parsing progress, DRAGEN Q30 discovery, samtools index creation, mosdepth start/finish and timing.

- Optional colored output via `rich`:
  
  - If `rich` is available in your environment, console messages are colored and formatted.
  - The script falls back to plain timestamped stderr output if `rich` is not installed.

- FASTQ progress interval:
  
  - To reduce noise, progress messages while streaming FASTQs are emitted every 1,000,000 reads (1M). This is intentionally coarse for large WGS FASTQs.

---

## Requirements

### Conda / micromamba

Recommended channels: `conda-forge` + `bioconda`.

Example `environment.yml` (includes recommended `rich` pin):

```/dev/null/environment.yml#L1-17
name: fast_qc_simple
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.9
  - mosdepth
  - samtools
  - rich=13.*
  - pandas
```

Install / update:

```/dev/null/install.sh#L1-4
micromamba create -f environment.yml
micromamba activate fast_qc_simple
# or: conda env update -f environment.yml -n fast_qc_simple
```

### Tools in PATH

- `mosdepth` (required)
- `samtools` (recommended; used if a CRAM index is missing)

---

## Inputs

Required:

- `--r1` FASTQ R1 (`.fastq` or `.fastq.gz`)
- `--r2` FASTQ R2
- `--cram` aligned CRAM
- `--bed` BED with the WGS loci (e.g. a small representative BED of ~400 loci)
- `--ref` reference FASTA (needed for CRAM/mosdepth)

Optional:

- `--dragen-metrics` path to DRAGEN/Mapper metrics for Q30%:
  - a `*mapping_metrics.csv` file
  - or a `*.metrics.tar.gz` containing `mapping_metrics.csv`
  - or a directory that contains either of the above (searched recursively)

Notes:

- CRAM decoding typically requires the exact reference used when writing the CRAM. Ensure `--ref` matches the CRAM and create a FASTA index if necessary (`samtools faidx ref.fa`).
- If a CRAM index is missing the script will attempt to create it using `samtools index -c` (requires `samtools` in PATH).

---

## Usage

Minimal:

```/dev/null/usage_min.sh#L1-8
./fast_qc_simple.py \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --cram sample.cram \
  --bed wgs_400loci.bed \
  --ref hg38.fa
```

With DRAGEN Q30 (preferred when available):

```/dev/null/usage_dragen.sh#L1-10
./fast_qc_simple.py \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --cram sample.cram \
  --bed wgs_400loci.bed \
  --ref hg38.fa \
  --dragen-metrics /path/to/mapping_metrics.csv
```

Verbose mode (console + logfile, colored if `rich` present):

```/dev/null/usage_verbose.sh#L1-6
./fast_qc_simple.py \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --cram sample.cram \
  --bed wgs_400loci.bed \
  --ref hg38.fa \
  --sample-id SAMPLE123 \
  --verbose
```

Custom output path / threads / min coverage threshold:

```/dev/null/usage_custom.sh#L1-8
./fast_qc_simple.py ... \
  --dragen-metrics /path/to/run_or_metrics_dir \
  --threads 16 \
  --mincov 20 \
  --out SAMPLE123.csv \
  --sample-id SAMPLE123
```

Output naming guidance:

- Recommended convention: set `--out <sample_id>.csv` (e.g. `--out SAMPLE123.csv`) so the filename matches the `sample_id` in the report.
- If you omit `--sample-id`, the script derives it from the CRAM filename (without `.cram/.bam/.sam`), and the logfile name follows that derived ID (`<sample_id>.log`).

Count duplicates (if you want to count duplicates in mosdepth results):

```/dev/null/usage_dup.sh#L1-2
./fast_qc_simple.py ... --count-dup
```

---

## QC overview (multi-sample)

This repository also includes `make_qc_overview.py`, which converts one or more single-sample QC CSVs into a compact overview table and adds convenience boolean columns (e.g. `q30_ok`, `depth_ok`, `fastq_ok`) plus an overall flag `qc_wgs_germline_overall`.

### Modellvorhaben QC (thresholds and overview support)

Reference (official QC document):

- https://www.bfarm.de/SharedDocs/Downloads/DE/Forschung/modellvorhaben-genomsequenzierung/Qs-durch-GRZ.pdf?__blob=publicationFile

For the Modellvorhaben short-read WGS germline QC, this project uses the following default criteria (as implemented in `fast_qc_simple.py` and reused by `make_qc_overview.py`):

- Q30%: `>= 85`
- mean read length: `> 70 bp`
- mean depth: `>= 30x`
- target fraction with coverage `>= --mincov` (default 20x): `>= 0.80`

How `make_qc_overview.py` helps:

- It reads the single-sample QC CSV(s) and computes boolean helper columns:
  - `q30_ok`, `readlen_ok`, `depth_ok`, `target20x_ok`, `fastq_ok`
- It combines them into a single flag: `qc_wgs_germline_overall` so you can quickly filter PASS/FAIL across many samples.
- It writes a compact overview table (`.tsv`, optional `.xlsx`) suitable for reporting and downstream checks.

Typical usage:

```/dev/null/qc_overview_usage.sh#L1-4
./make_qc_overview.py report.csv --out-tsv qc_overview.tsv
# Optional:
# ./make_qc_overview.py report.csv --out-xlsx qc_overview.xlsx
```

Note: `make_qc_overview.py` depends on `pandas` (included in the example `environment.yml` above).

---

## Output (CSV)

A single CSV is written (one row per run / sample). The header:

```/dev/null/csv_header#L1-1
sample_id,fastq_r1,fastq_r2,cram,bed,dragen_metrics,q30_percent,q30_pass,mean_read_length_bp,read_length_pass,reads_r1,reads_r2,readcount_match_pass,fastq_structure_pass,mean_depth_x,depth_pass,target_fraction_ge_20x,target20x_pass,qc_overall
```

Field notes:

- `q30_percent` – Q30% (prefer DRAGEN if provided; otherwise computed from FASTQ)
- `q30_pass` – PASS if Q30% ≥ 85
- `mean_read_length_bp` – mean read length across R1 + R2 (PASS if > 70 bp)
- `mean_depth_x` – mean depth from mosdepth summary (`total mean`) (PASS if ≥ 30x)
- `target_fraction_ge_20x` – fraction of BED bases with depth ≥ `--mincov` (default 20); PASS if ≥ 0.80 (80%)
- `*_pass` fields contain `PASS` or `FAIL`
- `qc_overall` will be `PASS` only if all individual checks pass

---

## QC criteria (defaults in the script)

- Q30%: `>= 85`
- mean read length: `> 70 bp`
- mean depth: `>= 30x`
- target fraction ≥ `--mincov` (default 20x): `>= 0.80`

---

## Logfile behavior

- When run with `--verbose`, the script:
  - Writes console messages (colored if `rich` is available).
  - Creates a plain-text logfile named `<sample_id>.log` in the current working directory (sample ID is taken from `--sample-id` or the CRAM filename if not provided).
  - Logfile entries are timestamped plain-text copies of the console messages (suitable for archival / automated parsing).
- On start (when `--verbose`) the script will log the chosen progress interval (1,000,000 reads) and the logfile path so you can immediately verify settings.

---

## Repository hygiene / generated files

This project can produce runtime artifacts in the working directory (for example `report.csv`, `<sample_id>.log`, `qc_overview.tsv`, `qc_overview.xlsx`). A `.gitignore` is included to avoid committing these generated files.

---

## Troubleshooting

- If you see no colored output, install `rich` (we pin `rich=13.*` in the `environment.yml`) or update your environment.
- `mosdepth` not found:
  - Ensure conda/micromamba environment is active, or install mosdepth into your environment.
- CRAM/reference mismatch:
  - Ensure `--ref` matches the CRAM's reference and that the .fai is present (`samtools faidx ref.fa`).
  - If CRAM index (.crai) is missing the script will attempt to create it using `samtools index -c <cram>`.