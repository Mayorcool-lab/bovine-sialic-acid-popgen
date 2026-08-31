#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter

files = {
    "H2": Path("results/phase_H_62_variant_population_frequencies.tsv"),
    "H3": Path("results/phase_H_62_variant_breed_frequencies.tsv"),
    "H4": Path("results/H4/H4_complete_62_variant_evidence.tsv"),
    "H5": Path("results/H5/H5_62_variant_conservation_evidence.tsv"),
    "H6": Path("results/H6/H6.2/processed/H6.2_complete_62_candidate_functional_context.tsv"),
}

def read_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def key(row):
    return (
        str(row["CHROM"]),
        str(row["POS"]),
        row["REF"],
        row["ALT"],
    )

tables = {name: read_tsv(path) for name, path in files.items()}
keysets = {}

print("=== ROW COUNTS ===")
for name, rows in tables.items():
    print(name, len(rows))

print("\n=== UNIQUE VARIANT KEYS ===")
for name, rows in tables.items():
    keys = [key(r) for r in rows]
    keysets[name] = set(keys)
    print(
        name,
        "unique=", len(set(keys)),
        "duplicate_rows=", len(keys) - len(set(keys)),
    )

reference = keysets["H2"]

print("\n=== KEYSET AGREEMENT AGAINST H2 ===")
for name in ["H4", "H5", "H6"]:
    print(
        name,
        "missing_from_layer=", len(reference - keysets[name]),
        "extra_in_layer=", len(keysets[name] - reference),
    )

print("\n=== H3 BREED-ROW ACCOUNTING ===")
h3_counts = Counter(key(r) for r in tables["H3"])
count_distribution = Counter(h3_counts.values())

print("H3 unique variants:", len(h3_counts))
print("rows-per-variant distribution:", dict(sorted(count_distribution.items())))
print(
    "H3 variants not represented by exactly 16 breed rows:",
    sum(1 for n in h3_counts.values() if n != 16),
)

print("\n=== SYMBOL CROSS-CHECK ===")
symbol_maps = {
    name: {key(r): r["SYMBOL"] for r in tables[name]}
    for name in ["H2", "H4", "H5", "H6"]
}

mismatches = []

for k in sorted(reference):
    symbols = {
        name: symbol_maps[name].get(k)
        for name in ["H2", "H4", "H5", "H6"]
    }
    if len(set(symbols.values())) != 1:
        mismatches.append((k, symbols))

print("Symbol mismatches across H2/H4/H5/H6:", len(mismatches))

passed = (
    len(reference) == 62
    and all(keysets[x] == reference for x in ["H4", "H5", "H6"])
    and len(h3_counts) == 62
    and all(n == 16 for n in h3_counts.values())
    and len(mismatches) == 0
)

print("\n=== FINAL JOINABILITY ===")
print("H2-H6 JOIN-KEY AUDIT:", "PASS" if passed else "FAIL")

if not passed:
    raise SystemExit(1)
