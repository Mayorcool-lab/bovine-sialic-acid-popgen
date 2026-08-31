#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

INFILE = Path(
    "results/H_integration/H3_62_candidate_16_breed_companion.tsv"
)

OUTFILE = Path(
    "results/H_integration/H3_62_candidate_breed_summary.tsv"
)


def read_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fmt_af(value):
    return f"{float(value):.8g}"


rows = read_tsv(INFILE)

by_candidate = defaultdict(list)

for r in rows:
    by_candidate[r["candidate_id"]].append(r)


def summarize_cohort(rows, cohort):
    x = [r for r in rows if r["COHORT"] == cohort]

    source_unavailable = [
        r for r in x
        if r["STATUS"] == "NOT_DETECTED_IN_SOURCE"
    ]

    absent_qc = [
        r for r in x
        if r["STATUS"] == "ABSENT_AFTER_QC"
    ]

    present = [
        r for r in x
        if r["STATUS"] == "PRESENT"
    ]

    present_ok = [
        r for r in present
        if r["N_FLAG"] == "OK"
    ]

    observed_ok = [
        r for r in present_ok
        if int(r["ALT_COUNT"]) > 0
    ]

    very_small = [
        r for r in present
        if r["N_FLAG"] == "VERY_SMALL_N"
    ]

    observed_very_small = [
        r for r in very_small
        if int(r["ALT_COUNT"]) > 0
    ]

    observed_ok = sorted(
        observed_ok,
        key=lambda r: r["BREED"]
    )

    observed_very_small = sorted(
        observed_very_small,
        key=lambda r: r["BREED"]
    )

    if observed_ok:
        afs = [float(r["AF"]) for r in observed_ok]
        min_af = min(afs)
        max_af = max(afs)

        max_breeds = sorted(
            r["BREED"]
            for r in observed_ok
            if float(r["AF"]) == max_af
        )

        observed_breeds = ";".join(
            r["BREED"] for r in observed_ok
        )

        observed_afs = ";".join(
            f'{r["BREED"]}:{fmt_af(r["AF"])}'
            for r in observed_ok
        )

        max_breed_text = ";".join(max_breeds)

    else:
        min_af = ""
        max_af = ""
        max_breed_text = ""
        observed_breeds = ""
        observed_afs = ""

    if observed_very_small:
        very_small_text = ";".join(
            f'{r["BREED"]}:{fmt_af(r["AF"])}'
            for r in observed_very_small
        )
    else:
        very_small_text = ""

    return {
        f"{cohort}_breed_rows": len(x),
        f"{cohort}_source_unavailable": len(source_unavailable),
        f"{cohort}_absent_after_qc": len(absent_qc),
        f"{cohort}_present": len(present),
        f"{cohort}_present_OK": len(present_ok),
        f"{cohort}_observed_OK": len(observed_ok),
        f"{cohort}_observed_OK_breeds": observed_breeds,
        f"{cohort}_observed_OK_AF": observed_afs,
        f"{cohort}_observed_OK_AF_min": (
            fmt_af(min_af) if min_af != "" else ""
        ),
        f"{cohort}_observed_OK_AF_max": (
            fmt_af(max_af) if max_af != "" else ""
        ),
        f"{cohort}_observed_OK_AF_max_breeds": max_breed_text,
        f"{cohort}_very_small_n_present": len(very_small),
        f"{cohort}_observed_very_small_n": len(observed_very_small),
        f"{cohort}_observed_very_small_n_AF": very_small_text,
    }


outrows = []

for candidate_id, sub in sorted(by_candidate.items()):
    if len(sub) != 16:
        raise RuntimeError(
            f"{candidate_id}: expected 16 breed rows, observed {len(sub)}"
        )

    symbols = {r["SYMBOL"] for r in sub}
    if len(symbols) != 1:
        raise RuntimeError(
            f"{candidate_id}: inconsistent SYMBOL values: {symbols}"
        )

    out = {
        "candidate_id": candidate_id,
        "SYMBOL": next(iter(symbols)),
    }

    out.update(
        summarize_cohort(sub, "South_Africa")
    )

    out.update(
        summarize_cohort(sub, "Europe")
    )

    outrows.append(out)


if len(outrows) != 62:
    raise RuntimeError(
        f"Expected 62 candidate summaries, observed {len(outrows)}"
    )

if len({r["candidate_id"] for r in outrows}) != 62:
    raise RuntimeError("Candidate IDs are not unique")


fieldnames = list(outrows[0].keys())

OUTFILE.parent.mkdir(parents=True, exist_ok=True)

with OUTFILE.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=fieldnames,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(outrows)


print("Output:", OUTFILE)
print("Rows:", len(outrows))
print("Columns:", len(fieldnames))
print(
    "Unique candidates:",
    len({r["candidate_id"] for r in outrows})
)
print(
    "South Africa breed-row counts:",
    sorted({r["South_Africa_breed_rows"] for r in outrows})
)
print(
    "Europe breed-row counts:",
    sorted({r["Europe_breed_rows"] for r in outrows})
)

ranking_terms = [
    c for c in fieldnames
    if "score" in c.lower()
    or "rank" in c.lower()
    or "priority" in c.lower()
]

print("Ranking/score columns:", ranking_terms)
print("AF values recalculated from genotypes:", 0)
print("H3 CANDIDATE SUMMARY BUILD: PASS")
