#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fast_qc_simple.py

Single-sample QC for WGS / short-read / germline (SE context), simplified:
Inputs: FASTQ R1, FASTQ R2, CRAM, BED (400 WGS loci), reference FASTA.
Optional: DRAGEN mapping_metrics (CSV or *.metrics.tar.gz or directory) for Q30%.

Outputs: one-row CSV report with PASS/FAIL flags and qc_overall.

Features:
- Optional colored/structured console output via `rich` (if installed)
- `--verbose` enables progress messages and writes a plain-text logfile named
  <sample_id>.log in the current working directory
- FASTQ parsing progress, samtools/mosdepth status and timing are logged
- Author metadata included
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import sys
import tarfile
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, TextIO, Tuple

# Metadata
__author__ = "Stephan Rosecker"
__email__ = "stephan.rosecker@ukbonn.de"
__version__ = "1.1.2"

# QC thresholds
Q30_REQUIRED = 85.0
MEAN_READ_LENGTH_REQUIRED_BP = 70.0  # strictly ">" in spec
MEAN_DEPTH_REQUIRED_X = 30.0
TARGET_FRACTION_GE_20X_REQUIRED = 0.80

REPORT_HEADER = [
    "sample_id",
    "fastq_r1",
    "fastq_r2",
    "cram",
    "bed",
    "dragen_metrics",
    "q30_percent",
    "q30_pass",
    "mean_read_length_bp",
    "read_length_pass",
    "reads_r1",
    "reads_r2",
    "readcount_match_pass",
    "fastq_structure_pass",
    "mean_depth_x",
    "depth_pass",
    "target_fraction_ge_20x",
    "target20x_pass",
    "qc_overall",
]


def _open_maybe_gz(path: Path):
    """Open file in binary mode; transparently handle .gz."""
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return open(path, "rb")


@dataclass
class FastqStats:
    reads: int
    total_bases: int
    mean_read_length_bp: float
    structure_ok: bool
    bases_q30: Optional[int] = None  # only set if computed


def count_reads(
    path: Path, verbose: bool = False, progress_interval: int = 1_000_000
) -> int:
    """
    Fast read counting for FASTQ files by counting newline bytes.

    Strategy:
    - For plain (uncompressed) files: use mmap to count b'\n' in one fast operation.
      This is generally the fastest approach for large files.
    - For .gz files: stream-read in larger chunks (64 MB) and count b'\n' occurrences.
    - Falls back to chunked read for mmap failures (e.g., empty files, platform limits).
    - Returns lines // 4 (reads). If an error occurs, returns 0.
    - Logs a warning if the total line count is not divisible by 4 (FASTQ invariants).

    Progress:
    - If verbose, logs compressed-byte progress for .gz inputs (best-effort).
      Note: for .gz this reflects compressed bytes consumed, not uncompressed bytes.
    """
    gz_chunk = 64 * 1024 * 1024  # 64 MB for gz streaming
    lines = 0
    try:
        if path.suffix == ".gz":
            filesize = None
            try:
                filesize = path.stat().st_size
            except Exception:
                filesize = None

            bytes_read_compressed = 0
            next_log_at = 256 * 1024 * 1024  # 256 MiB
            start_t = time.time()

            # gzip: can't mmap easily; stream in large chunks
            with gzip.open(path, "rb") as f:
                for chunk in iter(lambda: f.read(gz_chunk), b""):
                    bytes_read_compressed += len(chunk)
                    lines += chunk.count(b"\n")

                    if verbose and bytes_read_compressed >= next_log_at:
                        elapsed = max(0.001, time.time() - start_t)
                        mib_read = bytes_read_compressed / (1024 * 1024)
                        mib_s = mib_read / elapsed
                        if filesize:
                            pct = (bytes_read_compressed / filesize) * 100.0
                            log(
                                f"count_reads: read {mib_read:,.0f} MiB of {filesize / (1024 * 1024):,.0f} MiB ({pct:.1f}%) — {mib_s:,.1f} MiB/s",
                                "INFO",
                            )
                        else:
                            log(
                                f"count_reads: read {mib_read:,.0f} MiB (compressed) — {mib_s:,.1f} MiB/s",
                                "INFO",
                            )
                        next_log_at += 256 * 1024 * 1024
        else:
            # plain file: prefer mmap for fastest counting
            with open(path, "rb") as f:
                try:
                    import mmap as _mmap

                    # Map entire file (size 0 -> whole file) for read-only access
                    mm = _mmap.mmap(f.fileno(), 0, access=_mmap.ACCESS_READ)
                    try:
                        # Some type checkers do not know about mmap.mmap.count(...); use bytes.find loop instead.
                        data = mm[:]  # bytes
                        lines = data.count(b"\n")
                    finally:
                        mm.close()
                except (ValueError, OSError):
                    # Fallback: stream in moderately large chunks if mmap is not possible
                    chunk = 32 * 1024 * 1024  # 32 MB
                    f.seek(0)
                    for c in iter(lambda: f.read(chunk), b""):
                        lines += c.count(b"\n")
    except Exception as e:
        # Best-effort: log a warning if logger is available, then return 0
        try:
            log(f"count_reads: failed to count lines for {path}: {e}", "WARN")
        except Exception:
            pass
        return 0

    if lines % 4 != 0:
        try:
            log(
                f"count_reads: line count {lines} not divisible by 4 for {path}; using floor(lines/4)",
                "WARN",
            )
        except Exception:
            pass

    return max(0, lines // 4)


def _init_logger() -> Tuple[Callable[[str, str], None], bool]:
    """
    Initialize a console logger function.

    Returns:
      (log_fn, using_rich)

    log_fn(message, level) -> None
    level: "DEBUG", "INFO", "SUCCESS", "WARN", "ERROR"
    """
    try:
        # optional rich integration
        from rich.console import Console
        from rich.style import Style
        from rich.text import Text
        from rich.traceback import install as rich_traceback_install

        console = Console()
        rich_traceback_install(show_locals=False)

        level_styles = {
            "DEBUG": Style(color="cyan"),
            "INFO": Style(color="white"),
            "SUCCESS": Style(color="green", bold=True),
            "WARN": Style(color="yellow", bold=True),
            "ERROR": Style(color="red", bold=True),
        }

        def rich_log(message: str, level: str = "INFO") -> None:
            ts = time.strftime("%Y-%m-%d %H:%M:%S")
            lvl = (level or "INFO").upper()
            style = level_styles.get(lvl, Style(color="white"))
            ts_text = Text(f"[{ts}] ", style=Style(dim=True))
            lvl_text = Text(f"{lvl}: ", style=style)
            msg_text = Text(message)
            # Print combined text
            console.print(ts_text + lvl_text + msg_text)

        return rich_log, True

    except Exception:
        # fallback simple logger to stderr
        def simple_log(message: str, level: str = "INFO") -> None:
            ts = time.strftime("%Y-%m-%d %H:%M:%S")
            lvl = (level or "INFO").upper()
            sys.stderr.write(f"[{ts}] {lvl}: {message}\n")
            sys.stderr.flush()

        return simple_log, False


# initial console logger (may be wrapped later to also write to file)
log, USING_RICH = _init_logger()


def fastq_stats(
    path: Path,
    compute_q30: bool,
    verbose: bool = False,
    sample_reads: int = 0,
    tail_reads: int = 50000,
    total_reads: int = 0,
    progress_interval: int = 1_000_000,
) -> FastqStats:
    """
    Stream-parse FASTQ for:
      - reads count (exact)
      - structure checks (exact)
      - mean read length: either exact (if sample_reads==0) or sampled from first `sample_reads` reads
      - tail check: basic stats from last `tail_reads` reads (mean length)

    Notes:
      - compute_q30 is ignored here (Q30 should be provided by DRAGEN externally).
      - sample_reads: number of reads to sample from file start to estimate mean read length (0 = use all reads)
      - tail_reads: number of reads to keep in a rolling buffer to inspect the tail of the FASTQ
    """
    import collections

    start_t = time.time()
    reads = 0
    structure_ok = True

    # sampling accumulators
    sample_sum_len = 0
    sample_count = 0

    # If not sampling, we'll accumulate full sum (may be expensive)
    full_sum_len = 0

    # tail buffer to check read-ends
    tail_buf = collections.deque(maxlen=tail_reads)

    # Progress interval (coarse)
    progress_interval = (
        int(progress_interval)
        if progress_interval and progress_interval > 0
        else 1_000_000
    )

    try:
        with _open_maybe_gz(path) as fh:
            while True:
                header = fh.readline()
                if not header:
                    break  # EOF
                seq = fh.readline()
                plus = fh.readline()
                qual = fh.readline()

                if not seq or not plus or not qual:
                    structure_ok = False
                    break

                # strip line endings
                seq = seq.rstrip(b"\r\n")
                qual = qual.rstrip(b"\r\n")
                header = header.rstrip(b"\r\n")
                plus = plus.rstrip(b"\r\n")

                # Basic FASTQ checks
                if not header.startswith(b"@") or not plus.startswith(b"+"):
                    structure_ok = False
                if len(seq) != len(qual):
                    structure_ok = False

                read_len = len(seq)
                reads += 1

                # sampling logic: accumulate only for first sample_reads if sampling enabled
                if sample_reads and sample_count < sample_reads:
                    sample_sum_len += read_len
                    sample_count += 1
                elif not sample_reads:
                    full_sum_len += read_len

                # keep tail info
                if tail_reads > 0:
                    tail_buf.append(read_len)

                # progress logging
                if verbose and (reads % progress_interval == 0):
                    elapsed_s = max(0.001, time.time() - start_t)
                    throughput = reads / elapsed_s  # reads/s
                    if total_reads and total_reads > 0:
                        pct = (reads / total_reads) * 100.0
                        remaining = max(0, total_reads - reads)
                        eta_s = (remaining / throughput) if throughput > 0 else None
                        if eta_s is not None:
                            eta_m = int(round(eta_s / 60.0))
                            log(
                                f"FASTQ {path.name}: processed {reads:,}/{total_reads:,} ({pct:.1f}%) reads — {throughput:,.0f} r/s — ETA {eta_m}m",
                                "INFO",
                            )
                        else:
                            log(
                                f"FASTQ {path.name}: processed {reads:,}/{total_reads:,} ({pct:.1f}%) reads — {throughput:,.0f} r/s",
                                "INFO",
                            )
                    else:
                        log(
                            f"FASTQ {path.name}: processed {reads:,} reads — {throughput:,.0f} r/s",
                            "INFO",
                        )

    except Exception as e:
        structure_ok = False
        log(f"Error while reading FASTQ {path}: {e}", "ERROR")

    # determine mean read length: sampled if sample_count>0, otherwise full
    if sample_count > 0:
        mean_len = (sample_sum_len / sample_count) if sample_count > 0 else 0.0
        total_bases_for_report = int(
            mean_len * reads
        )  # approximate total bases for downstream fields if needed
    else:
        mean_len = (full_sum_len / reads) if reads > 0 else 0.0
        total_bases_for_report = full_sum_len

    # tail statistics (simple check: mean tail length)
    tail_mean = (sum(tail_buf) / len(tail_buf)) if len(tail_buf) > 0 else None

    elapsed = time.time() - start_t
    if verbose:
        tail_msg = f" tail_mean={tail_mean:.2f}" if tail_mean is not None else ""
        log(
            f"Finished FASTQ {path.name}: reads={reads:,} sample_count={sample_count} mean_len={mean_len:.2f} elapsed={elapsed:.1f}s{tail_msg}",
            "INFO",
        )

    # We store the sampled/approximated values in the FastqStats fields:
    return FastqStats(
        reads=reads,
        total_bases=total_bases_for_report,
        mean_read_length_bp=mean_len,
        structure_ok=structure_ok,
        bases_q30=None,
    )


def _read_q30_from_mapping_metrics_text(text: str) -> Optional[float]:
    """
    Find a CSV row where row[2] contains "Q30 bases" and parse row[4] as percent.
    Compatible with older DRAGEN mapping metrics layout.
    """
    import io

    reader = csv.reader(io.StringIO(text))
    for row in reader:
        if len(row) >= 5 and "Q30 bases" in row[2]:
            try:
                return float(row[4].strip().replace("%", "").replace(",", "."))
            except Exception:
                return None
    return None


def extract_dragen_q30_percent(
    path: Optional[Path], verbose: bool = False
) -> Optional[float]:
    """
    Accepts:
      - mapping_metrics.csv
      - *.metrics.tar.gz containing *.mapping_metrics.csv
      - directory containing one of the above (searched recursively)

    Returns Q30% or None if not found/parsable.
    """
    if not path:
        return None

    if path.is_dir():
        if verbose:
            log(f"Searching for mapping_metrics in directory {path}", "INFO")
        csv_cands = sorted(
            [p for p in path.rglob("*mapping_metrics.csv") if p.is_file()]
        )
        for p in csv_cands:
            q = extract_dragen_q30_percent(p, verbose=verbose)
            if q is not None:
                return q
        tar_cands = sorted([p for p in path.rglob("*.metrics.tar.gz") if p.is_file()])
        for p in tar_cands:
            q = extract_dragen_q30_percent(p, verbose=verbose)
            if q is not None:
                return q
        return None

    # file cases
    if path.name.endswith(".metrics.tar.gz"):
        try:
            with tarfile.open(path, "r:gz") as tf:
                for m in tf.getmembers():
                    if m.name.endswith(".mapping_metrics.csv") or m.name.endswith(
                        "mapping_metrics.csv"
                    ):
                        f = tf.extractfile(m)
                        if f:
                            txt = f.read().decode("utf-8", errors="replace")
                            q = _read_q30_from_mapping_metrics_text(txt)
                            if q is not None and 0.0 <= q <= 100.0:
                                return q
        except Exception as e:
            if verbose:
                log(f"Error reading tarball {path}: {e}", "ERROR")
        return None

    # text/csv file
    try:
        txt = path.read_text(encoding="utf-8", errors="replace")
        q = _read_q30_from_mapping_metrics_text(txt)
        if q is not None and 0.0 <= q <= 100.0:
            return q
        return None
    except Exception as e:
        if verbose:
            log(f"Error reading mapping metrics file {path}: {e}", "ERROR")
        return None


def ensure_cram_index(cram: Path, verbose: bool = False) -> None:
    """
    Ensure CRAM index exists. Tries to create it using samtools if missing.
    Raises RuntimeError if samtools not found or indexing fails.
    """
    crai1 = cram.with_suffix(cram.suffix + ".crai")  # sample.cram.crai
    crai2 = cram.with_suffix(".crai")  # sample.cram (less common)

    if crai1.exists() or crai2.exists():
        if verbose:
            log(f"CRAM index found for {cram.name}", "INFO")
        return

    cmd = ["samtools", "index", "-c", str(cram)]
    if verbose:
        log(f"CRAM index missing; running: {' '.join(cmd)}", "INFO")
    try:
        proc = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
    except FileNotFoundError:
        raise RuntimeError("CRAM index missing and samtools not found in PATH.")
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or "samtools index failed")
    if verbose:
        log(f"samtools index completed for {cram.name}", "INFO")


def run_mosdepth(
    cram: Path,
    fasta: Path,
    bed: Path,
    threads: int,
    mincov: int,
    count_dup: bool,
    verbose: bool = False,
) -> Tuple[Optional[float], Optional[float]]:
    """
    Run mosdepth to compute mean depth and fraction of target bases >= mincov.

    Returns (target_fraction_ge_mincov, mean_depth_x) or (None, None) on failure.
    """
    start_t = time.time()
    filter_flags = 772 if count_dup else 1796  # keep compatibility with older script
    with tempfile.TemporaryDirectory() as td:
        prefix = Path(td) / "mosdepth"

        cmd = [
            "mosdepth",
            "--threads",
            str(threads),
            "--fasta",
            str(fasta),
            "--no-per-base",
            "--by",
            str(bed),
            "--thresholds",
            str(mincov),
            "-F",
            str(filter_flags),
            str(prefix),
            str(cram),
        ]
        if verbose:
            log(f"Running mosdepth: {' '.join(cmd)}", "INFO")
        try:
            proc = subprocess.run(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
        except FileNotFoundError:
            log("mosdepth not found in PATH", "ERROR")
            return None, None
        if proc.returncode != 0:
            log(f"mosdepth error: {proc.stderr.strip()}", "ERROR")
            return None, None

        thresholds = prefix.with_suffix(".thresholds.bed.gz")
        summary = prefix.with_suffix(".mosdepth.summary.txt")
        if not thresholds.exists() or not summary.exists():
            log("mosdepth output files missing", "ERROR")
            return None, None

        # compute target fraction >= mincov
        covered = 0
        total = 0
        try:
            with gzip.open(thresholds, "rt") as f:
                for line in f:
                    if not line or line.startswith("#"):
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 4:
                        continue
                    try:
                        covered += int(cols[-1])
                        total += int(cols[2]) - int(cols[1])
                    except Exception:
                        continue
            tf = (covered / total) if total > 0 else 0.0
        except Exception as e:
            log(f"Error reading mosdepth thresholds: {e}", "ERROR")
            tf = None

        # read mean depth from summary
        mean_depth = None
        try:
            with open(summary, "r", encoding="utf-8", errors="replace") as f:
                header = f.readline().strip().split("\t")
                if "mean" in header:
                    mean_idx = header.index("mean")
                    for line in f:
                        if line.startswith("total"):
                            parts = line.strip().split("\t")
                            mean_depth = float(parts[mean_idx])
                            break
        except Exception as e:
            log(f"Error reading mosdepth summary: {e}", "ERROR")
            mean_depth = None

        elapsed = time.time() - start_t
        if verbose:
            log(
                f"mosdepth finished in {elapsed:.1f}s (tf={tf if tf is not None else 'NA'} mean_depth={mean_depth if mean_depth is not None else 'NA'})",
                "INFO",
            )
        return tf, mean_depth


def bool_str(v: Optional[bool]) -> str:
    if v is None:
        return ""
    return "PASS" if v else "FAIL"


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Simplified single-sample QC for WGS short-read germline: FASTQ R1/R2 + CRAM + BED + reference; optional DRAGEN mapping_metrics for Q30%."
    )
    ap.add_argument(
        "--r1", required=True, type=Path, help="FASTQ R1 (.fastq/.fastq.gz)"
    )
    ap.add_argument(
        "--r2", required=True, type=Path, help="FASTQ R2 (.fastq/.fastq.gz)"
    )
    ap.add_argument("--cram", required=True, type=Path, help="Aligned CRAM")
    ap.add_argument(
        "--bed", required=True, type=Path, help="BED with WGS loci (e.g. 400 loci)"
    )
    ap.add_argument(
        "--ref", required=True, type=Path, help="Reference FASTA for CRAM/mosdepth"
    )
    ap.add_argument(
        "--dragen-metrics",
        type=Path,
        default=None,
        help="mapping_metrics.csv or *.metrics.tar.gz or directory containing them (optional)",
    )
    ap.add_argument(
        "--threads", type=int, default=16, help="Threads for mosdepth (default: 16)"
    )
    ap.add_argument(
        "--mincov",
        type=int,
        default=20,
        help="Minimum coverage threshold for targets (default: 20)",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("report.csv"),
        help="Output CSV path (default: report.csv)",
    )
    ap.add_argument(
        "--sample-id",
        default=None,
        help="Sample ID for report (default: derived from CRAM filename)",
    )
    ap.add_argument(
        "--count-dup",
        action="store_true",
        help="Count duplicates (default: duplicates excluded)",
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose progress output (default: False)",
    )
    ap.add_argument(
        "--progress-interval",
        type=int,
        default=1_000_000,
        help="Progress logging interval in reads for FASTQ parsing (default: 1,000,000)",
    )

    args = ap.parse_args()
    verbose = bool(args.verbose)

    # Determine sample_id early so we can open sample_id.log when verbose
    sample_id = args.sample_id or args.cram.name
    sample_id = re.sub(
        r"(\.cram|\.bam|\.sam)(\.gz)?$", "", sample_id, flags=re.IGNORECASE
    )

    # If verbose, open a logfile named <sample_id>.log in cwd and wrap logger
    log_fh: Optional[TextIO] = None
    log_path: Optional[Path] = None
    if verbose:
        try:
            log_path = Path.cwd() / f"{sample_id}.log"
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_fh = open(log_path, "a", encoding="utf-8")
            orig_log = log

            def _dual_log(message: str, level: str = "INFO") -> None:
                # write to console (orig_log) first
                try:
                    orig_log(message, level)
                except Exception:
                    pass
                # also append plain text to logfile
                try:
                    if log_fh is None:
                        return
                    ts = time.strftime("%Y-%m-%d %H:%M:%S")
                    lvl = (level or "INFO").upper()
                    log_fh.write(f"[{ts}] {lvl}: {message}\n")
                    log_fh.flush()
                except Exception:
                    pass

            # Install dual logger globally
            globals()["log"] = _dual_log
            log(f"Logging to {log_path}", "INFO")
            log(
                f"Verbose logging enabled; progress interval = {args.progress_interval:,} reads",
                "INFO",
            )
        except Exception as e:
            log(f"Could not open logfile {log_path}: {e}", "WARN")
            log_fh = None

    # Basic path checks
    for p in [args.r1, args.r2, args.cram, args.bed, args.ref]:
        if not p.exists():
            log(f"ERROR: file not found: {p}", "ERROR")
            try:
                if log_fh:
                    log_fh.close()
            except Exception:
                pass
            return 2

    if verbose:
        log(f"Starting QC for sample {sample_id}", "INFO")
        log(
            f"Inputs: r1={args.r1}, r2={args.r2}, cram={args.cram}, bed={args.bed}, ref={args.ref}",
            "INFO",
        )
        if args.dragen_metrics:
            log(f"DRAGEN metrics input: {args.dragen_metrics}", "INFO")
        if USING_RICH:
            log("Using rich for colored/structured output", "DEBUG")
        else:
            log("rich not available; using plain logging", "DEBUG")

    # DRAGEN Q30 first (optional)
    dragen_q30 = (
        extract_dragen_q30_percent(args.dragen_metrics, verbose=verbose)
        if args.dragen_metrics
        else None
    )
    if verbose:
        if dragen_q30 is not None:
            log(f"Using DRAGEN Q30% = {dragen_q30:.2f}", "INFO")
        else:
            log(
                "No valid DRAGEN Q30% found; will compute from FASTQ if requested",
                "INFO",
            )

    # FASTQ stats: use sampled mean read length (1% of reads) and tail check; Q30 from DRAGEN only
    try:
        # Count reads in R1 to determine sampling size (1%)
        total_reads_r1 = count_reads(
            args.r1, verbose=verbose, progress_interval=args.progress_interval
        )
        sample_reads = (
            max(1, int(total_reads_r1 * 0.01)) if total_reads_r1 > 0 else 100000
        )
        if verbose:
            log(
                f"Using sample_reads={sample_reads} (1% of {total_reads_r1} reads) for mean read length estimation",
                "INFO",
            )

        r1_stats = fastq_stats(
            args.r1,
            compute_q30=False,
            verbose=verbose,
            sample_reads=sample_reads,
            tail_reads=50000,
            total_reads=total_reads_r1,
            progress_interval=args.progress_interval,
        )
        r2_total_reads = 0
        try:
            r2_total_reads = count_reads(
                args.r2, verbose=verbose, progress_interval=args.progress_interval
            )
        except Exception:
            r2_total_reads = 0
        r2_stats = fastq_stats(
            args.r2,
            compute_q30=False,
            verbose=verbose,
            sample_reads=sample_reads,
            tail_reads=50000,
            total_reads=r2_total_reads,
            progress_interval=args.progress_interval,
        )
    except KeyboardInterrupt:
        log("Interrupted while reading FASTQs", "WARN")
        try:
            if log_fh:
                log_fh.close()
        except Exception:
            pass
        return 130

    reads_r1 = r1_stats.reads
    reads_r2 = r2_stats.reads
    readcount_match_pass = (reads_r1 == reads_r2) and (reads_r1 > 0)

    fastq_structure_pass = (
        r1_stats.structure_ok
        and r2_stats.structure_ok
        and (reads_r1 > 0)
        and (reads_r2 > 0)
    )

    total_reads = reads_r1 + reads_r2
    # Compute weighted mean read length from the per-file estimated means
    mean_read_length_bp = (
        (
            (
                (r1_stats.mean_read_length_bp * reads_r1)
                + (r2_stats.mean_read_length_bp * reads_r2)
            )
            / total_reads
        )
        if total_reads > 0
        else 0.0
    )
    read_length_pass = mean_read_length_bp > MEAN_READ_LENGTH_REQUIRED_BP

    # Q30 percent: use DRAGEN metrics only (do not compute from FASTQ)
    q30_percent = float(dragen_q30) if dragen_q30 is not None else None
    q30_pass = (q30_percent is not None) and (q30_percent >= Q30_REQUIRED)

    # Coverage from CRAM
    try:
        ensure_cram_index(args.cram, verbose=verbose)
    except Exception as e:
        log(f"ERROR: {e}", "ERROR")
        tf = None
        mean_depth_x = None
    else:
        tf, mean_depth_x = run_mosdepth(
            cram=args.cram,
            fasta=args.ref,
            bed=args.bed,
            threads=args.threads,
            mincov=args.mincov,
            count_dup=args.count_dup,
            verbose=verbose,
        )

    depth_pass = (mean_depth_x is not None) and (mean_depth_x >= MEAN_DEPTH_REQUIRED_X)
    target20x_pass = (tf is not None) and (tf >= TARGET_FRACTION_GE_20X_REQUIRED)

    qc_overall = all(
        [
            bool(fastq_structure_pass),
            bool(readcount_match_pass),
            bool(read_length_pass),
            bool(q30_pass),
            bool(depth_pass),
            bool(target20x_pass),
        ]
    )

    row = {
        "sample_id": sample_id,
        "fastq_r1": str(args.r1),
        "fastq_r2": str(args.r2),
        "cram": str(args.cram),
        "bed": str(args.bed),
        "dragen_metrics": str(args.dragen_metrics) if args.dragen_metrics else "",
        "q30_percent": f"{q30_percent:.2f}" if q30_percent is not None else "",
        "q30_pass": bool_str(q30_pass),
        "mean_read_length_bp": f"{mean_read_length_bp:.2f}",
        "read_length_pass": bool_str(read_length_pass),
        "reads_r1": str(reads_r1),
        "reads_r2": str(reads_r2),
        "readcount_match_pass": bool_str(readcount_match_pass),
        "fastq_structure_pass": bool_str(fastq_structure_pass),
        "mean_depth_x": f"{mean_depth_x:.2f}" if mean_depth_x is not None else "",
        "depth_pass": bool_str(depth_pass),
        "target_fraction_ge_20x": f"{tf:.6f}" if tf is not None else "",
        "target20x_pass": bool_str(target20x_pass),
        "qc_overall": bool_str(qc_overall),
    }

    # Write CSV (create parent dirs if needed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=REPORT_HEADER)
        w.writeheader()
        w.writerow(row)

    # Also print a short summary to stdout
    print(f"qc_overall={row['qc_overall']} sample_id={sample_id} out={args.out}")
    if verbose:
        log(f"Report written to {args.out}", "INFO")
        log(f"QC overall: {row['qc_overall']}", "INFO")
        try:
            if log_fh:
                log(f"Closing logfile {log_path}", "DEBUG")
                log_fh.close()
        except Exception:
            pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
