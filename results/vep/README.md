# results/vep/ — File Guide

## AUTHORITATIVE, CURRENT (post-CMAH-boundary-correction)
- `phase_g_protein_altering_FINAL.tsv` — the definitive 62-variant
  protein-altering candidate table. USE THIS ONE.
- `phase_g_canonical_annotations_corrected.tsv` — full 10,261-variant
  canonical annotation, corrected CMAH window
- `phase_g_final_target_genes_corrected.tsv` — 13-gene table, VEP-assigned
- `cmah_rescued_region_override.tsv` — CMAH rescued-region variants
  (393), including 5 manually-annotated coding variants
- `cmah_refseq_rescued_coding_annotation.txt` — positive-control output
  from the RefSeq rescue validation (reproduces P424L, F512Y exactly)
- `phase_g_excluded_offtarget_corrected.tsv` — off-target/blank-symbol
  variants excluded from the 13-gene table (includes CMAH rescue region
  before override applied)
- `all_variants_annotated_corrected.vcf` — raw VEP output, corrected
  10,261-variant union

## SUPERSEDED (pre-CMAH-boundary-correction; kept for provenance only)
- `phase_g_protein_altering_final.tsv` (lowercase) — old 58-variant
  table, CMAH boundary not yet corrected
- `phase_g_canonical_annotations.tsv` / `_fixed.tsv` — old annotation
  passes, superseded by `_corrected.tsv`
- `phase_g_final_target_genes.tsv` — old 13-gene table
- `phase_g_excluded_offtarget.tsv` — old exclusion table
- `all_variants_annotated.vcf` (no `_corrected` suffix, if present) —
  original VEP run, before boundary fix

Do not use the superseded files for any new analysis, manuscript text,
or downstream handoff. See methods/phase_G_updated_workflow.md for the
full history of why the correction happened.
