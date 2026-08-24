#!/usr/bin/env python3
"""
Manual codon-level annotation for CMAH variants in the RefSeq-rescued
region (23:32,731,973-32,753,209) that Ensembl VEP cannot annotate,
since this region falls outside Ensembl's own CMAH gene boundary
(gene ends at 32,731,972) but is real, RefSeq-confirmed coding sequence
per transcript XM_024984024.1 (the transcript used in Ogun et al. 2023).

Validated via exact reproduction of two independently published variants
(P424L, F512Y). Script halts with an error if either fails to reproduce
exactly, or if any REF base mismatches the reference genome.
"""

codon_table = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

cds_blocks = [
    (32709570,32709577),(32709817,32710020),(32713880,32713971),
    (32721556,32721665),(32723849,32724039),(32727617,32727764),
    (32728093,32728174),(32728921,32729058),(32731794,32731929),
    (32739238,32739378),(32743898,32744072),(32745757,32745873),
    (32750626,32750741),(32753134,32753209),
]

def genomic_to_cdna(pos):
    cumulative = 0
    for s, e in cds_blocks:
        if s <= pos <= e:
            return cumulative + (pos - s) + 1
        cumulative += e - s + 1
    return None

import subprocess

FASTA = "data/reference/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz"

def get_seq(chrom, start, end):
    r = subprocess.run(
        ["samtools","faidx",FASTA,f"{chrom}:{start}-{end}"],
        capture_output=True, text=True, check=True
    )
    return "".join(r.stdout.strip().split("\n")[1:]).upper()

full_cds = "".join(get_seq("23", s, e) for s, e in cds_blocks)

assert len(full_cds) == 1734

variants = [
    (32739285,"A","G","N386S"),
    (32739292,"C","G","I388M"),
    (32743918,"C","T","P424L"),
    (32745866,"T","A","F512Y"),
    (32750628,"T","A","I515I"),
]

print("POS\tREF>ALT\tcDNA\tCODON\tCODON_CHANGE\tAA_CHANGE\tCONSEQUENCE\tEXPECTED")

for pos, ref, alt, expected in variants:
    cdna_pos = genomic_to_cdna(pos)
    codon_num = ((cdna_pos - 1) // 3) + 1
    pos_in_codon = (cdna_pos - 1) % 3
    codon_start = ((cdna_pos - 1) // 3) * 3
    ref_codon = full_cds[codon_start:codon_start+3]
    observed_ref = full_cds[cdna_pos - 1]
    if observed_ref != ref:
        raise RuntimeError(f"REF mismatch at {pos}: FASTA={observed_ref}, VCF={ref}")
    alt_codon = list(ref_codon)
    alt_codon[pos_in_codon] = alt
    alt_codon = "".join(alt_codon)
    ref_aa = codon_table[ref_codon]
    alt_aa = codon_table[alt_codon]
    if ref_aa == alt_aa:
        consequence = "synonymous_variant"
    elif alt_aa == "*":
        consequence = "stop_gained"
    elif ref_aa == "*":
        consequence = "stop_lost"
    else:
        consequence = "missense_variant"
    aa_change = f"{ref_aa}{codon_num}{alt_aa}"
    print(f"{pos}\t{ref}>{alt}\tc.{cdna_pos}\t{codon_num}\t{ref_codon}>{alt_codon}\t{aa_change}\t{consequence}\t{expected}")
    if aa_change != expected:
        raise RuntimeError(f"Validation failed at {pos}: expected {expected}, got {aa_change}")

print("\nPASS: all five CMAH RefSeq coding consequences reproduced exactly.")
