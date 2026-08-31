#!/usr/bin/env python3

from pathlib import Path
import subprocess
import shutil
import sys

ALIGNMENT_DIR = Path("results/H5/alignments")
OUTPUT_DIR = Path("results/H5/rate4site")

ALIGNMENTS = [
    "CMAH_aligned.fasta",
    "CMAS_aligned.fasta",
    "GNE_aligned.fasta",
    "NANP_aligned.fasta",
    "NANS_aligned.fasta",
    "NEU1_aligned.fasta",
    "NEU3_aligned.fasta",
    "SLC35A1_aligned.fasta",
    "ST3GAL1_aligned.fasta",
    "ST3GAL4_aligned.fasta",
    "ST6GAL1_aligned.fasta",
    "ST8SIA2_aligned.fasta",
    "ST8SIA4_aligned.fasta",
    "CMAH_RefSeq_X2_aligned.fasta",
]

rate4site = shutil.which("rate4site")
if rate4site is None:
    sys.exit("ERROR: rate4site executable not found")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

for filename in ALIGNMENTS:
    alignment = ALIGNMENT_DIR / filename

    if not alignment.exists():
        sys.exit(f"ERROR: missing alignment: {alignment}")

    stem = filename.replace("_aligned.fasta", "")

    score_file = OUTPUT_DIR / f"{stem}_rate4site.txt"
    tree_file = OUTPUT_DIR / f"{stem}_tree.nwk"
    log_file = OUTPUT_DIR / f"{stem}_rate4site.log"

    cmd = [
        rate4site,
        "-s", str(alignment),
        "-o", str(score_file),
        "-x", str(tree_file),
        "-Mj",
        "-ib",
    ]

    print(f"Running {stem} ...", flush=True)

    with log_file.open("w") as log:
        result = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )

    if result.returncode != 0:
        sys.exit(
            f"ERROR: Rate4Site failed for {stem}. "
            f"See {log_file}"
        )

    if not score_file.exists() or score_file.stat().st_size == 0:
        sys.exit(f"ERROR: missing/empty score file for {stem}")

    if not tree_file.exists() or tree_file.stat().st_size == 0:
        sys.exit(f"ERROR: missing/empty tree file for {stem}")

    print(f"PASS {stem}")

print()
print(f"Rate4Site executable: {rate4site}")
print(f"Alignment sets completed: {len(ALIGNMENTS)}")
print("H5F Rate4Site production run: PASS")
