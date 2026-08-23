# sialic_acid_popgen — Full Analysis Workflow

## Scientific question

Determine how naturally occurring SNP variation across the bovine sialic-acid
metabolic network differs among European and South African cattle populations,
which variants are predicted to affect biological function, and whether pathway
loci show unusual population differentiation or signatures of selection relative
to the genome-wide background.

## Discovery cohorts

South African cohort: 120 cattle (Afrikaner 40, Bonsmara 20, Nguni 40, Tuli 20)
European cohort: 226 verified cattle (of 289 total; 63 Portuguese animals
excluded pending manual ID verification - see notes/STATUS.md)
Total: 346 cattle
No German breed represented (investigated and closed - see notes/STATUS.md
for DSN audit and broader search results)

## Candidate network (13 provisional genes)

GNE, NANS, NANP, CMAS, SLC35A1, CMAH, ST8SIA2, ST8SIA4, ST6GAL1, ST3GAL1,
ST3GAL4, NEU1, NEU3

## Core rule: SNPs vs indels

Primary scientific analysis = SNPs, especially high-quality biallelic SNPs.
Indels are retained in the normalized archive for completeness but are not part
of the main population-genomic analysis unless an obviously high-impact
functional indel emerges.

RAW VCF -> NORMALIZED MASTER -> SNP analysis dataset (MAIN STUDY)
                              -> indels retained (secondary/archive only)

## PHASE A — Source-data verification

A1. South African source: DONE. ARS-UCD1.2, unphased, 240 animals in source
VCF, 120 target animals identified via optibov_id/vcf_id, FreeBayes v1.0.2
(header) vs v1.3.1 (papers' text).

A2. European source: DONE. Identical header discrepancy confirmed (same
FreeBayes v1.0.2 / fileDate 20240226 as African file).

## PHASE B — Targeted 13-gene extraction

B1. South African: DONE. 120 animals, 13 genes, 10,765 raw regional records.
B2. European: DONE. 226 verified animals, 13 genes, 7,628 raw regional records.

PHASE B COMPLETE: 346 total animals, 18,393 total raw regional records
across both cohorts, before normalization/QC.

## PHASE C — Normalization

C1. bcftools norm against ARS-UCD1.2: DONE both cohorts. SA: 10,896 records.
EU: 7,701 records. Full variant-type breakdown reconciled (SNP/indel/MNP/other).
C2. Normalized master preserved (SNPs + indels).
C3. Biallelic SNP dataset derived: SA 9,684 sites, EU 6,842 sites. Reconciled
exactly against normalization summary.

## PHASE D — QC characterization

COMPLETE. DP (median ~8-9x both cohorts), GQ (median ~48-54), Ts/Tv (2.48 SA,
2.52 EU - confirms real biological data), per-animal missingness (very low,
<4% both), per-site missingness (quantified distribution, median 0%, only
0.11-0.19% of sites >10% missing), site QUAL (0 sites below 20, upstream
filter already effective), heterozygosity (pooled cohorts showed Wahlund
effect; breed-stratified: Europe fully normalized across all 10 breeds,
South Africa partially - Bonsmara/Nguni normalized, Afrikaner/Tuli remained
elevated, honestly unresolved). Full details in
methods/phase_D_qc_characterization_summary.md

## PHASE E — Apply QC

COMPLETE. DP>=4, GQ>=20 (genotype-level masking), site call rate >=90% after
masking. Final QC-passed SNPs: SA 7,552, EU 6,439. Independently verified via
machine-regenerated site retention table. Cohort retention differs (SA 78%,
EU 94% overall) - documented as unattributed technical dataset characteristic,
NOT explained by sample size alone (that claim was checked and found
insufficient). Full details in methods/phase_E_qc_application_summary.md

## PHASE F — Functional annotation (VEP) — NOT YET STARTED

Run VEP against ARS-UCD1.2 on the QC-passed SNPs. For each SNP: gene,
transcript, consequence, codon change, amino-acid change, position, impact
category. Apply cattle-applicable functional prediction tools from V1 spec
for missense variants.

## PHASE G — Population-level allele analysis — NOT YET STARTED

For every pathway SNP: allele frequency per breed, per continent, private/
near-private alleles, breed-enriched alleles.

## Preliminary manuscript track vs. collaborator handoff (overview)

Phases A-E: Data preparation + normalization + empirical QC
  -> Phase F/G: VEP functional annotation of the 13-gene dataset
  -> Functional variant landscape (consequences, missense, LoF/splice,
     allele-frequency patterns, etc.)
  -> HANDOFF TO AIYESA
  -> Aiyesa builds genome-wide dataset using same QC methodology
     (DP>=4, GQ>=20, 90% call rate - established Phase D/E)
  -> Separate genome-wide population-genomics track (PCA / population
     structure / FST / selection analyses as appropriate)
  -> Integrate targeted functional results + genome-wide population-
     genetic context
  -> PREPRINT / PRELIMINARY MANUSCRIPT

Handoff package to Aiyesa (after Phase F/G completes) includes:
- 13 QC-filtered SNP VCFs, both cohorts (South Africa 7,552 / Europe 6,439)
- VEP annotation results (consequence, transcript, AA change, impact)
- 409-animal metadata (breed, country, dataset)
- README stating exact QC methodology applied and explicit scope limits:
  this package supports within-gene diversity and targeted breed-level
  differentiation, but NOT PCA/ADMIXTURE/genome-background FST/selection
  scans, which require his separately-built genome-wide dataset.

Idowu/Mtileni engage once Aiyesa's population-structure/differentiation
results exist, for biological/breed interpretation, not before.

## PHASE I onward — genome-wide population genetics (Aiyesa's track)

PCA, ADMIXTURE, genome-wide FST, nucleotide diversity, Tajima's D, phasing,
iHS, XP-EHH. Requires processing full genomes chromosome-by-chromosome, not
the 13-gene panel. Interpretation principle: selection evidence is strongest
when multiple independent statistics converge AND the signal is extreme
relative to genomic background.
