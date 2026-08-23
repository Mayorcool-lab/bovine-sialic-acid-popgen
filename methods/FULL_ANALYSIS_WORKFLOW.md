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
No German breed represented (investigated and closed - see notes/STATUS.md)

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

A2. European source: DONE. Identical header discrepancy confirmed.

## PHASE B — Targeted 13-gene extraction

B1. South African: DONE. 120 animals, 13 genes, 10,765 raw regional records.
B2. European: DONE. 226 verified animals, 13 genes, 7,628 raw regional records.

PHASE B COMPLETE: 346 total animals, 18,393 total raw regional records.

## PHASE C — Normalization

DONE both cohorts. SA: 10,896 records. EU: 7,701 records. Biallelic SNP
dataset derived: SA 9,684 sites, EU 6,842 sites.

## PHASE D — QC characterization

COMPLETE. Full details in methods/phase_D_qc_characterization_summary.md.
DP median ~8-9x, GQ median ~48-54, Ts/Tv 2.48 SA / 2.52 EU, per-animal
missingness <4% both, per-site missingness median 0% (0.11-0.19% >10%
missing), site QUAL all >=20, heterozygosity breed-stratified (Europe fully
normalized, South Africa partially - Afrikaner/Tuli remained elevated,
unresolved honestly).

## PHASE E — Apply QC

COMPLETE. Full details in methods/phase_E_qc_application_summary.md.
DP>=4, GQ>=20, site call rate >=90% after masking. Final QC-passed SNPs:
SA 7,552, EU 6,439. Cohort retention differs (SA 78%, EU 94%) - documented
as unattributed technical characteristic, not sample-size alone.

## PHASE F — Functional annotation — NOT YET STARTED

Run VEP against ARS-UCD1.2 on QC-passed SNPs. Gene, transcript, consequence,
codon change, amino-acid change, position, impact category. Apply cattle-
applicable functional prediction tools from V1 spec for missense variants.

## PHASE G — Population-level allele analysis — NOT YET STARTED

Per-breed and per-continent allele frequency, private/near-private alleles,
breed-enriched alleles. Biological story starts appearing here, does not
yet prove selection.

## Preliminary manuscript track vs. collaborator handoff (overview)

Phases A-E: Data preparation + normalization + empirical QC
  -> Phase F/G: VEP functional annotation of the 13-gene dataset
  -> Functional variant landscape
  -> HANDOFF TO AIYESA
  -> Aiyesa builds genome-wide dataset using same QC methodology
  -> Separate genome-wide population-genomics track
  -> Integrate targeted functional results + genome-wide context
  -> PREPRINT / PRELIMINARY MANUSCRIPT

Handoff package to Aiyesa (after Phase F/G completes):
- 13 QC-filtered SNP VCFs, both cohorts
- VEP annotation results
- 409-animal metadata
- README with exact QC methodology and explicit scope limits (supports
  within-gene diversity/differentiation, NOT PCA/ADMIXTURE/genome-wide FST/
  selection scans, which need his separately-built genome-wide dataset)

Idowu/Mtileni engage once Aiyesa's population-structure results exist.

## PHASE I — Genome-wide population-genomics dataset (Aiyesa's track)

Note on scale: substantially larger than any step completed so far. PCA/
ADMIXTURE/genome-background FST/iHS/XP-EHH cannot use candidate genes alone.
Requires processing full African (~75GB) and European (~60GB) genomes,
chromosome-by-chromosome.

## PHASE J — PCA

Genome-wide, LD-pruned SNPs. Question: what is the actual genetic structure
of these cattle? Do breeds cluster as expected, is there admixture, do
breed labels agree with genomic ancestry? Comes before final population
comparisons.

## PHASE K — ADMIXTURE

LD-pruned genome-wide SNP set. Evaluate ancestry structure across plausible
K values. PCA + ADMIXTURE inform 17f.

## PHASE L — Freeze population comparisons (17f)

Only after PCA/ADMIXTURE. Likely hierarchy: within South Africa, within
Europe, Africa vs Europe (interpreted cautiously). Aiyesa/Idowu advise;
final decision remains yours.

## PHASE M — Genome-wide FST

Sliding windows, 50kb/20kb step (Kim et al. 2025 precedent). Determine
where the 13 pathway genes fall in the genome-wide empirical distribution
(percentile, not raw FST value alone).

## PHASE N — Nucleotide diversity (π)

Genome-wide diversity. Is diversity unusually reduced/elevated around
candidate genes? Interpret alongside population history.

## PHASE O — Tajima's D

Genome-wide sliding-window. Strongly negative = excess rare variants,
potential sweep/expansion. Strongly positive = balancing selection/
structure. Interpret against genome-wide distribution.

## PHASE P — Phasing

Both datasets confirmed unphased via VCF headers. Phase genome-wide/
chromosome-scale data (Beagle or agreed method) before haplotype analysis.

## PHASE Q — iHS

Phased chromosome-scale data. Unusually long haplotype within a population
suggests recent positive selection. Compute genome/chromosome-wide first,
then examine candidate loci.

## PHASE R — XP-EHH

Phased chromosome-scale data, between populations. Haplotype approaching
fixation in one population but not another. Compute genome/chromosome-wide
first, then interrogate candidate loci.

## PHASE S — Integrate selection statistics

Do not declare selection from one extreme statistic alone. Strongest
evidence: convergence of high FST + unusual Tajima's D + reduced π +
extreme iHS + extreme XP-EHH at the same locus.
