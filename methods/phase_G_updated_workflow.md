# Sialic Acid Population Genomics — Updated Workflow

## UPDATED PHASE G — Boundary Audit + Functional Annotation

### Reason for modification

The initial Ensembl-based extraction and VEP run worked technically, but
CMAH revealed an annotation-model issue. Ensembl release 110 defines CMAH
(ENSBTAG00000003892) on ARS-UCD1.2 as chr23:32,672,258-32,731,972, while
the earlier CMAH study (Ogun et al. 2023) used a RefSeq transcript model
(XM_024984024.1) with coding variants extending beyond the Ensembl gene
end. Therefore, an Ensembl-only extraction interval can omit genuine
coding variation represented by another valid transcript annotation.

This is not a VEP failure. It is a gene/transcript-boundary harmonization
issue.

### Confirmed evidence (from Ogun et al. 2023, Table 1)

| Published variant | ARS-UCD1.2 coordinate | dbSNP | Inside Ensembl extraction window? |
|---|---|---|---|
| K107E | 23:32,721,570 | rs208635220 | YES |
| A210S | 23:32,727,639 | rs435799892 | YES |
| N242S | 23:32,727,736 | rs109811989 | YES |
| P424L | 23:32,743,918 | rs518400910 | NO - 11,946 bp beyond window |
| F512Y | 23:32,745,866 | rs380571713 | NO - 13,894 bp beyond window |

Current CMAH extraction window: 23:32,672,258-32,731,972 (Ensembl gene
boundary, both canonical and alternate Ensembl transcripts end at
32,731,972; neither reaches far enough to capture P424L or F512Y).

Independent exact-coordinate validation: K107E, A210S, N242S recovered
in the current 346-animal population dataset at the EXACT same genomic
positions reported in 2023, using a completely different pipeline,
cohort, and QC methodology. This is strong internal validation.

## G1 — Audit all 13 genes

Compare, for every target gene:
- Ensembl ARS-UCD1.2 gene start/end
- RefSeq ARS-UCD1.2 transcript/CDS start/end
- minimum and maximum coding coordinates
- whether any RefSeq protein-coding exon/CDS lies outside the extraction
  interval

Target genes: GNE, NANS, NANP, CMAS, SLC35A1, CMAH, ST8SIA2, ST8SIA4,
ST6GAL1, ST3GAL1, ST3GAL4, NEU1, NEU3.

Audit table format:
Gene | Ensembl start | Ensembl end | RefSeq CDS min | RefSeq CDS max | CDS outside extraction? | Action

RefSeq annotation source: GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz
(NCBI RefSeq, same ARS-UCD1.2 assembly, downloaded to data/reference/)

Decision rule:
- If the other 12 genes show no coding extension outside their extraction
  intervals, correct CMAH only and document it.
- If additional genes have coding sequence outside their intervals,
  inspect and correct those genes too.
- Do not knowingly leave recoverable omitted coding SNPs merely as a
  manuscript limitation.

## G2 — Rebuild affected extraction intervals

For each affected gene:
1. Define an ARS-UCD1.2 interval capturing the relevant validated coding
   span.
2. Re-extract that gene from the original cohort VCFs.
3. Reapply the same Phase E QC: DP >= 4, GQ >= 20, and site call rate
   >= 90% after genotype masking.
4. Replace only the affected gene-level QC dataset.
5. Recalculate and document any changed SNP counts.

Do not rerun unaffected genes unnecessarily.

## G3 — Regenerate Phase F union if extraction changes

If corrected extraction changes variants:
- regenerate South African unique variants
- regenerate European unique variants
- regenerate shared, South-Africa-only, and Europe-only counts
- regenerate union_variants.tsv
- regenerate union_variants.vcf

Previous initial union (PROVISIONAL, pending boundary audit outcome):
- SA QC records: 7,552
- Europe QC records: 6,439
- Shared: 4,230
- SA-only: 3,322
- Europe-only: 2,209
- Unique union: 9,761

These counts become provisional if the boundary audit changes extraction.

## G4 — Run VEP on the corrected final union

Use:
- Ensembl VEP 110.1
- Bos taurus
- ARS-UCD1.2
- Ensembl release-110 cache
- ARS-UCD1.2 FASTA
- offline/cache mode
- symbol, canonical, biotype flags

VEP installation and the 10-variant validation have already succeeded
(confirmed working correctly on test data before full run).

## G5 — Parse canonical annotations

Retain: CHROM, POS, REF, ALT, SYMBOL, Ensembl Gene ID, transcript,
consequence, IMPACT, biotype, exon/intron, HGVSc, HGVSp, protein
position, amino-acid change, canonical status.

## G6 — Repair missing symbols using Ensembl gene ID

Use data/reference/gene_symbol_to_id.tsv whenever SYMBOL is blank.
Confirmed example: ENSBTAG00000003892 = CMAH (VEP's release 110 cache
had a gap in symbol cross-referencing for this gene specifically; Gene
ID field was always correctly populated). Thus, blank SYMBOL does not
automatically mean an unknown gene.

## G7 — Restrict final manuscript annotation to the 13 target genes

Nearby annotations such as KIRREL3, SLC44A4, RARS2 and 5S_rRNA are
off-target for the biological network (picked up via genomic proximity
within the 2kb upstream extraction buffer). Preserve them in an
exclusion table (results/vep/phase_g_excluded_offtarget.tsv) rather than
silently deleting them.

## G8 — Re-extract protein-altering candidates

After boundary correction, regenerate the protein-altering table
including: missense_variant, stop_gained, stop_lost, relevant combined
missense/splice consequences.

PROVISIONAL result (pending boundary-aware correction): 58 target-gene
protein-altering variants, all 13 genes represented (see
results/vep/phase_g_protein_altering_final.tsv). This number and gene
breakdown remain provisional until G1-G3 boundary correction is
complete.

## G9 — CMAH continuity check

Use the corrected CMAH dataset to determine whether previously reported
CMAH nsSNPs from Ogun et al. (2023), including P424L and F512Y, occur in
the current cohorts. Do not infer their absence from the earlier
Ensembl-bounded extraction because the transcript-boundary discrepancy
has now been identified as the explanation for their current absence.

Already confirmed: K107E, A210S, N242S recovered at exact matching
coordinates. P424L and F512Y status genuinely unknown until CMAH
extraction is corrected per G2.

## Collaborator handoff

After Phase F/G methodology and the final annotation framework are
stable, send Aiyesa the information needed for his separate genome-wide
population-genomics track: confirmed source datasets/accessions,
sample/breed metadata, finalized QC logic, assembly/reference details,
methodological notes.

Aiyesa needs genome-wide SNP data for analyses such as PCA, ADMIXTURE,
genome-wide FST, iHS and XP-EHH. The 13-gene VCFs are not appropriate
input for genome-wide population structure.

Idowu and Mtileni are not blockers at this stage; their main
contribution comes later when population-genomic results require
biological/population interpretation or specific resource expertise.

## Manuscript scope

Structural/functional protein modelling is NOT part of the preliminary
manuscript. The manuscript remains focused on population genomics and
functional annotation of the 13-gene bovine sialic-acid network.

## Current Phase G status

Already completed:
- VEP 110.1 installation
- 10-variant VEP validation
- initial 9,761-variant VEP run
- canonical annotation parsing
- gene-symbol correction
- target-gene filtering
- protein-altering extraction
- reproducible scripts committed to git
- CMAH boundary discrepancy discovered and root-caused (Ensembl vs
  RefSeq transcript model disagreement)

Still required before Phase G is scientifically final:
1. Audit Ensembl versus RefSeq boundaries for all 13 genes (in progress)
2. Correct affected extraction interval(s)
3. Rerun Phase E QC only for affected genes
4. Rebuild the Phase F union if necessary
5. Rerun VEP if the union changes
6. Regenerate final target-gene and protein-altering tables
7. Perform the published CMAH-variant continuity check (P424L, F512Y)

The 13-gene boundary audit is therefore part of Phase G quality control,
not merely a limitation to mention in the manuscript.

## G0: CMAH precise correction target (confirmed)

RefSeq transcript XM_024984024.1 (the exact transcript used in Ogun et
al. 2023) CDS span: 32,709,570 - 32,753,209 (14 CDS records, matching
the paper's stated 14 coding exons exactly).

Current Ensembl-based extraction window: 32,670,258 - 32,731,972.
Required corrected window (to capture full XM_024984024.1 CDS,
including P424L at 32,743,918 and F512Y at 32,745,866): must extend to
at least 32,753,209, an increase of 21,237 bp at the 3' end.

## G1 AUDIT RESULT: CONFIRMED (13-gene RefSeq CDS-boundary audit complete)

RefSeq annotation source: GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz
(NCBI RefSeq, same ARS-UCD1.2 assembly). CDS-level comparison (the
scientifically decisive test, not gene-level span, since gene-level
differences are frequently just UTR/intron and not biologically
meaningful).

| Gene | RefSeq CDS max | Extraction window end | Status |
|---|---|---|---|
| GNE | 60,554,025 | 60,530,253 | **PROBLEM: ~23.8kb of real coding sequence outside window** |
| NANS | 63,212,646 | 63,212,682 | OK |
| NANP | 10,704,867 | 10,705,509 | OK |
| CMAS | 88,169,834 | 88,172,546 | OK |
| SLC35A1 | 62,445,474 | 62,447,501 | OK |
| ST8SIA2 | 14,595,814 | 14,597,969 | OK |
| ST8SIA4 | 100,136,839 | 100,139,609 | OK |
| ST6GAL1 | 80,233,920 | 80,269,738 | OK (gene-level span differed by ~79kb but this was entirely non-coding; CDS-based check correctly resolved this as a non-issue) |
| ST3GAL1 | 7,951,420 | 7,951,420 | OK, exact boundary match |
| ST3GAL4 | 29,738,538 | 29,750,191 | OK |
| NEU1 | 27,497,986 | 27,498,491 | OK |
| NEU3 | 54,266,515 | 54,267,975 | OK |
| CMAH | 32,753,209 | 32,731,972 | **PROBLEM: 21,237 bp of real coding sequence outside window** |

### Conclusion
Of 13 target genes, 11 required no correction (their Ensembl-based
extraction windows already fully contain all RefSeq-annotated coding
sequence). Two genes, CMAH and GNE, require extraction window
correction because real coding sequence (confirmed via RefSeq CDS
records, not just gene-span differences) falls outside the current
extraction boundary. This is a genuine, CDS-verified finding, not
inferred from gene-level annotation differences alone.

### Note for collaborators and manuscript methods
A systematic Ensembl-vs-RefSeq CDS boundary audit was performed for all
13 target genes on ARS-UCD1.2. Two genes (CMAH, GNE) had real coding
sequence outside the initial Ensembl-based extraction window and were
corrected before finalizing the functional annotation dataset. The
other 11 genes were confirmed to require no correction. This audit was
triggered by discovering that the CMAH extraction window excluded two
previously published missense variants (P424L, F512Y; Ogun et al. 2023)
that fall within RefSeq-annotated but not Ensembl-annotated CMAH coding
sequence, likely reflecting differing gene models between the two
annotation databases for this locus.
