#!/usr/bin/env python3
"""
Phase G5/G6: parse corrected VEP annotation (10,261-variant union set,
post CMAH boundary correction), selecting the canonical transcript per
variant, and applying the gene-symbol fallback (using
gene_symbol_to_id.tsv) whenever VEP leaves SYMBOL blank but Gene ID is
populated. Produces one row per unique variant.
"""

import csv
import re
from collections import Counter

id_to_symbol = {}
with open("data/reference/gene_symbol_to_id.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        id_to_symbol[row["ensembl_gene_id"]] = row["gene_symbol"]

severity = {
    "transcript_ablation": 1, "splice_acceptor_variant": 2, "splice_donor_variant": 3,
    "stop_gained": 4, "frameshift_variant": 5, "stop_lost": 6, "start_lost": 7,
    "transcript_amplification": 8, "inframe_insertion": 9, "inframe_deletion": 10,
    "missense_variant": 11, "protein_altering_variant": 12, "splice_region_variant": 13,
    "incomplete_terminal_codon_variant": 14, "start_retained_variant": 15,
    "stop_retained_variant": 16, "synonymous_variant": 17, "coding_sequence_variant": 18,
    "mature_miRNA_variant": 19, "5_prime_UTR_variant": 20, "3_prime_UTR_variant": 21,
    "non_coding_transcript_exon_variant": 22, "intron_variant": 23,
    "NMD_transcript_variant": 24, "non_coding_transcript_variant": 25,
    "upstream_gene_variant": 26, "downstream_gene_variant": 27, "TFBS_ablation": 28,
    "TFBS_amplification": 29, "TF_binding_site_variant": 30,
    "regulatory_region_ablation": 31, "regulatory_region_amplification": 32,
    "regulatory_region_variant": 33, "intergenic_variant": 34,
}

vcf = "results/vep/all_variants_annotated_corrected.vcf"
outfile = "results/vep/phase_g_canonical_annotations_corrected.tsv"

fields = None
rows = []

with open(vcf) as f:
    for line in f:
        if line.startswith("##INFO=<ID=CSQ"):
            m = re.search(r"Format: ([^\">]+)", line)
            fields = m.group(1).split("|")
            continue
        if line.startswith("#"):
            continue
        parts = line.rstrip().split("\t")
        chrom, pos, vid, ref, alt = parts[:5]
        info = parts[7]
        csq_value = None
        for item in info.split(";"):
            if item.startswith("CSQ="):
                csq_value = item[4:]
                break
        if not csq_value:
            continue
        annotations = []
        for annotation in csq_value.split(","):
            values = annotation.split("|")
            values += [""] * (len(fields) - len(values))
            ann = dict(zip(fields, values))
            annotations.append(ann)
        canonical = [a for a in annotations if a.get("CANONICAL") == "YES"]
        candidates = canonical if canonical else annotations

        def rank_annotation(a):
            consequences = a.get("Consequence", "").split("&")
            return min((severity.get(c, 999) for c in consequences if c), default=999)

        best = min(candidates, key=rank_annotation)

        symbol = best.get("SYMBOL", "")
        gene_id = best.get("Gene", "")
        if not symbol and gene_id in id_to_symbol:
            symbol = id_to_symbol[gene_id]

        rows.append({
            "CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt,
            "SYMBOL": symbol, "Gene": gene_id, "Feature": best.get("Feature", ""),
            "Consequence": best.get("Consequence", ""), "IMPACT": best.get("IMPACT", ""),
            "BIOTYPE": best.get("BIOTYPE", ""), "EXON": best.get("EXON", ""),
            "INTRON": best.get("INTRON", ""), "HGVSc": best.get("HGVSc", ""),
            "HGVSp": best.get("HGVSp", ""), "Protein_position": best.get("Protein_position", ""),
            "Amino_acids": best.get("Amino_acids", ""), "Codons": best.get("Codons", ""),
            "Existing_variation": best.get("Existing_variation", ""),
            "CANONICAL": best.get("CANONICAL", ""),
        })

columns = ["CHROM","POS","REF","ALT","SYMBOL","Gene","Feature","Consequence","IMPACT",
           "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","Protein_position","Amino_acids",
           "Codons","Existing_variation","CANONICAL"]

with open(outfile, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Variants parsed: {len(rows)}")
print(f"Unique keys: {len({(r['CHROM'],r['POS'],r['REF'],r['ALT']) for r in rows})}")
print()
print("=== Genes (after fix) ===")
for gene, n in Counter(r["SYMBOL"] for r in rows).most_common():
    print(f"{gene or '[STILL NO SYMBOL]'}\t{n}")
