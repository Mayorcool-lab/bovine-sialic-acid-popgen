#!/usr/bin/env python3

"""
H5: Retrieve and freeze orthologue evidence for conservation analysis.

Bovine reference:
    Ensembl release 110 / ARS-UCD1.2

Comparative orthology:
    Ensembl Compara REST release 116

Design:
    - Fixed 18-species comparative framework.
    - Automatic inclusion only for one-to-one orthologues.
    - One-to-many relationships are flagged, never chosen automatically.
    - Gene-specific missingness is retained explicitly.
    - Bos taurus reference sequence is handled separately from the
      release-110 bovine peptide FASTA.
"""

import csv
import json
import time
from pathlib import Path
from urllib.request import Request, urlopen

EXPECTED_COMPARA_RELEASE = 116

GENES = {
    "CMAH": "ENSBTAG00000003892",
    "CMAS": "ENSBTAG00000010136",
    "GNE": "ENSBTAG00000001142",
    "NANP": "ENSBTAG00000000686",
    "NANS": "ENSBTAG00000019023",
    "NEU1": "ENSBTAG00000005674",
    "NEU3": "ENSBTAG00000025931",
    "SLC35A1": "ENSBTAG00000014057",
    "ST3GAL1": "ENSBTAG00000001156",
    "ST3GAL4": "ENSBTAG00000016506",
    "ST6GAL1": "ENSBTAG00000010357",
    "ST8SIA2": "ENSBTAG00000048501",
    "ST8SIA4": "ENSBTAG00000021445",
}

TARGET_SPECIES = [
    "ovis_aries",
    "capra_hircus",
    "bison_bison_bison",
    "bos_grunniens",
    "camelus_dromedarius",
    "sus_scrofa",
    "equus_caballus",
    "canis_lupus_familiaris",
    "felis_catus",
    "homo_sapiens",
    "mus_musculus",
    "rattus_norvegicus",
    "loxodonta_africana",
    "vombatus_ursinus",
    "gallus_gallus",
    "anolis_carolinensis",
    "xenopus_tropicalis",
    "danio_rerio",
]

OUTDIR = Path("results/H5/orthologues")


def get_json(url, max_attempts=3):
    req = Request(url, headers={"Content-Type": "application/json"})

    for attempt in range(1, max_attempts + 1):
        try:
            with urlopen(req, timeout=30) as response:
                return json.load(response)
        except Exception:
            if attempt == max_attempts:
                raise

            wait_seconds = 2 ** (attempt - 1)
            print(
                f"REST attempt {attempt} failed; "
                f"retrying in {wait_seconds}s ...",
                flush=True,
            )
            time.sleep(wait_seconds)

    raise RuntimeError("Unreachable retry state")


def get_ensembl_release():
    data = get_json(
        "https://rest.ensembl.org/info/software?content-type=application/json"
    )
    return int(data["release"])


def retrieve_gene(symbol, gene_id):
    url = (
        f"https://rest.ensembl.org/homology/id/bos_taurus/{gene_id}"
        "?type=orthologues;sequence=protein"
    )

    data = get_json(url)
    homs = data["data"][0]["homologies"]

    by_species = {}
    for h in homs:
        sp = h["target"]["species"]
        if sp in TARGET_SPECIES:
            by_species.setdefault(sp, []).append(h)

    rows = []

    for sp in TARGET_SPECIES:
        matches = by_species.get(sp, [])

        if not matches:
            rows.append({
                "gene": symbol,
                "species": sp,
                "status": "NOT_RETURNED",
                "orthology_type": "",
                "source_protein_id": "",
                "target_gene_id": "",
                "target_protein_id": "",
                "source_protein_length": "",
                "target_protein_length": "",
                "length_ratio": "",
                "source_perc_id": "",
                "target_perc_id": "",
                "source_perc_pos": "",
                "target_perc_pos": "",
            })
            continue

        one2one = [h for h in matches if h["type"] == "ortholog_one2one"]

        if len(one2one) == 1:
            h = one2one[0]
            source = h["source"]
            target = h["target"]

            source_seq = source.get("align_seq", "").replace("-", "")
            target_seq = target.get("align_seq", "").replace("-", "")

            source_len = len(source_seq)
            target_len = len(target_seq)

            rows.append({
                "gene": symbol,
                "species": sp,
                "status": "ONE2ONE",
                "orthology_type": h["type"],
                "source_protein_id": source.get("protein_id", ""),
                "target_gene_id": target.get("id", ""),
                "target_protein_id": target.get("protein_id", ""),
                "source_protein_length": source_len,
                "target_protein_length": target_len,
                "length_ratio": (
                    f"{target_len / source_len:.4f}" if source_len else ""
                ),
                "source_perc_id": source.get("perc_id", ""),
                "target_perc_id": target.get("perc_id", ""),
                "source_perc_pos": source.get("perc_pos", ""),
                "target_perc_pos": target.get("perc_pos", ""),
                "_target_sequence": target_seq,
            })
        else:
            types = sorted({h["type"] for h in matches})
            rows.append({
                "gene": symbol,
                "species": sp,
                "status": "AMBIGUOUS",
                "orthology_type": ",".join(types),
                "source_protein_id": "",
                "target_gene_id": ",".join(
                    sorted(h["target"].get("id", "") for h in matches)
                ),
                "target_protein_id": ",".join(
                    sorted(h["target"].get("protein_id", "") for h in matches)
                ),
                "source_protein_length": "",
                "target_protein_length": "",
                "length_ratio": "",
                "source_perc_id": "",
                "target_perc_id": "",
                "source_perc_pos": "",
                "target_perc_pos": "",
            })

    return rows


def main():
    print("H5 orthologue retrieval configuration")
    print("------------------------------------")
    print(f"Expected Compara release: {EXPECTED_COMPARA_RELEASE}")
    print(f"Genes: {len(GENES)}")
    print(f"Comparative target species: {len(TARGET_SPECIES)}")
    print(f"Output directory: {OUTDIR}")

    observed_release = get_ensembl_release()

    print()
    print(f"Observed Ensembl REST release: {observed_release}")

    if observed_release != EXPECTED_COMPARA_RELEASE:
        raise RuntimeError(
            f"Ensembl release mismatch: expected "
            f"{EXPECTED_COMPARA_RELEASE}, observed {observed_release}"
        )

    print("Release guard: PASS")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "gene",
        "species",
        "status",
        "orthology_type",
        "source_protein_id",
        "target_gene_id",
        "target_protein_id",
        "source_protein_length",
        "target_protein_length",
        "length_ratio",
        "source_perc_id",
        "target_perc_id",
        "source_perc_pos",
        "target_perc_pos",
    ]

    all_rows = []

    print()
    print("Retrieving orthologues:")

    for i, (symbol, gene_id) in enumerate(GENES.items(), start=1):
        print(f"[{i}/{len(GENES)}] {symbol} ... ", end="", flush=True)

        rows = retrieve_gene(symbol, gene_id)
        all_rows.extend(rows)

        outfile = OUTDIR / f"{symbol}_orthologue_metadata.tsv"

        with outfile.open("w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=fieldnames,
                delimiter="\t",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(rows)

        fasta_out = OUTDIR / f"{symbol}_one2one_sequences.fasta"

        with fasta_out.open("w") as fh:
            for row in rows:
                if row["status"] != "ONE2ONE":
                    continue

                seq = row.get("_target_sequence", "")
                if not seq:
                    raise RuntimeError(
                        f"Missing target sequence for "
                        f"{symbol} {row['species']}"
                    )

                protein_id = row["target_protein_id"]
                header = (
                    f">{symbol}|{row['species']}|"
                    f"{protein_id}|Compara{EXPECTED_COMPARA_RELEASE}"
                )

                fh.write(header + "\n")
                for j in range(0, len(seq), 60):
                    fh.write(seq[j:j + 60] + "\n")

        counts = {}
        for row in rows:
            counts[row["status"]] = counts.get(row["status"], 0) + 1

        summary = ", ".join(
            f"{status}={counts[status]}" for status in sorted(counts)
        )
        print(f"done ({summary})")

    combined = OUTDIR / "all_13_genes_orthologue_metadata.tsv"

    with combined.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(all_rows)

    print()
    print(f"Genes processed: {len(GENES)}")
    print(f"Expected rows: {len(GENES) * len(TARGET_SPECIES)}")
    print(f"Rows written: {len(all_rows)}")
    print(f"Combined metadata: {combined}")


if __name__ == "__main__":
    main()
