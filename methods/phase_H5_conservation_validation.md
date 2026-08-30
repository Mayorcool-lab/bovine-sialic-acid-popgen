# Phase H5 — Evolutionary Conservation Analysis: Methods and Validation Record

## H5 objective

Phase H5 evaluates evolutionary conservation at protein positions affected by the 62 Phase G protein-altering candidate variants.

Conservation is treated as an independent evidence layer for later candidate prioritisation. It does not by itself demonstrate functional impairment, pathogenicity, phenotypic effect, or selection.

No structural modelling is performed in Phase H5.

## H5.1 — Candidate protein-coordinate inventory

Authoritative candidate input:

`results/vep/phase_g_protein_altering_FINAL.tsv`

Results:

- Total protein-altering candidates: 62
- Candidates with Ensembl Feature/transcript ID: 58
- Candidates without Ensembl Feature/transcript ID: 4
- Candidates with protein position: 62
- Candidates without protein position: 0

The four candidates without an Ensembl Feature are the manually rescued CMAH variants:

- CMAH 23:32739285 A>G — N386S
- CMAH 23:32739292 C>G — I388M
- CMAH 23:32743918 C>T — P424L
- CMAH 23:32745866 T>A — F512Y

These variants retain manually established RefSeq protein positions and must not be interpreted using the shorter Ensembl CMAH protein coordinate system.

## H5.2 — Local alignment-software inventory

The `vep_env` environment was checked before conservation analysis.

Available:

- MAFFT
- MUSCLE
- Clustal Omega

Not available at this stage:

- Rate4Site

MAFFT is the preferred multiple-sequence alignment tool for the planned reproducible conservation workflow. No conservation-scoring software was installed or run at this stage.

## H5.3 — Phase G transcript architecture

Across the 62 candidates, 13 unique Ensembl transcript IDs were represented.

Twelve genes used a single Ensembl transcript context for all candidates. CMAH contained one Ensembl transcript context plus the separate manual RefSeq rescue context.

The 13 Ensembl transcript contexts were:

- CMAH — ENSBTAT00000081903
- CMAS — ENSBTAT00000013373
- GNE — ENSBTAT00000001520
- NANP — ENSBTAT00000000918
- NANS — ENSBTAT00000025318
- NEU1 — ENSBTAT00000067415
- NEU3 — ENSBTAT00000036738
- SLC35A1 — ENSBTAT00000018680
- ST3GAL1 — ENSBTAT00000001538
- ST3GAL4 — ENSBTAT00000075324
- ST6GAL1 — ENSBTAT00000083430
- ST8SIA2 — ENSBTAT00000077891
- ST8SIA4 — ENSBTAT00000028590

## H5.4 — Ensembl transcript-to-protein mapping

Transcript-to-protein mapping was derived from:

`data/reference/Bos_taurus.ARS-UCD1.2.110.gtf.gz`

All 13 Phase G Ensembl transcript contexts mapped to exactly one protein ID.

Results:

- Unique Phase G Ensembl transcripts: 13
- Exactly one protein ID: 13
- Missing protein ID: 0
- Multiple protein IDs: 0

Mappings:

- CMAH — ENSBTAT00000081903 → ENSBTAP00000073225
- CMAS — ENSBTAT00000013373 → ENSBTAP00000013373
- GNE — ENSBTAT00000001520 → ENSBTAP00000001520
- NANP — ENSBTAT00000000918 → ENSBTAP00000000918
- NANS — ENSBTAT00000025318 → ENSBTAP00000025318
- NEU1 — ENSBTAT00000067415 → ENSBTAP00000073421
- NEU3 — ENSBTAT00000036738 → ENSBTAP00000036593
- SLC35A1 — ENSBTAT00000018680 → ENSBTAP00000018680
- ST3GAL1 — ENSBTAT00000001538 → ENSBTAP00000001538
- ST3GAL4 — ENSBTAT00000075324 → ENSBTAP00000073357
- ST6GAL1 — ENSBTAT00000083430 → ENSBTAP00000063809
- ST8SIA2 — ENSBTAT00000077891 → ENSBTAP00000057846
- ST8SIA4 — ENSBTAT00000028590 → ENSBTAP00000028590

## H5.5 — Release-matched bovine protein reference

To avoid mixing annotation releases, the Bos taurus Ensembl release 110 peptide reference corresponding to ARS-UCD1.2 was obtained:

`data/reference/Bos_taurus.ARS-UCD1.2.pep.all.fa.gz`

Source:

Ensembl release 110, Bos taurus ARS-UCD1.2 `pep.all` FASTA.

All 13 expected protein IDs were present exactly once.

Results:

- Exactly once: 13
- Missing: 0
- Multiple: 0

Reference protein lengths:

- CMAH — 383 aa
- CMAS — 434 aa
- GNE — 722 aa
- NANP — 248 aa
- NANS — 359 aa
- NEU1 — 416 aa
- NEU3 — 428 aa
- SLC35A1 — 337 aa
- ST3GAL1 — 339 aa
- ST3GAL4 — 395 aa
- ST6GAL1 — 405 aa
- ST8SIA2 — 375 aa
- ST8SIA4 — 359 aa

## H5.6 — Candidate reference-residue coordinate validation

Before orthologue retrieval or conservation scoring, Phase G protein positions were validated directly against the release-matched Ensembl 110 protein sequences.

The four manually rescued CMAH variants were excluded from this Ensembl-coordinate check because their positions derive from a separate RefSeq protein model.

For the 58 Ensembl-annotated candidates:

- Ordinary amino-acid positions checked: 57
- Reference amino-acid matches: 57/57
- Reference amino-acid mismatches: 0
- Stop-lost positions checked separately: 1
- Coordinate/reference errors: 0

The single stop-lost candidate was:

- ST3GAL4 29:29738209 G>T
- Consequence: stop_lost
- Reference protein length: 395 aa
- VEP stop position: 396
- Validation: PASS because the affected termination codon occurs immediately after the 395-aa reference protein.

Therefore, the apparent out-of-range position is expected for a stop-lost variant and is not a coordinate mismatch.

### CMAH annotation-model safeguard

The Ensembl release 110 CMAH protein ENSBTAP00000073225 is 383 aa long.

The four manually rescued CMAH variants occur at RefSeq-derived positions 386, 388, 424, and 512. These positions therefore cannot be mapped onto the 383-aa Ensembl protein and will remain on a separate RefSeq-coordinate conservation track.

This independently confirms the annotation-model difference previously encountered during Phase G/H4 and prevents inappropriate mixing of Ensembl and RefSeq CMAH protein coordinates.

## H5 pre-conservation validation status

PASS.

The Ensembl-coordinate candidates have an exact and validated chain:

genomic allele → Phase G transcript → Ensembl 110 protein → protein residue.

Orthologue retrieval and conservation scoring must preserve this mapping.

The four manually rescued CMAH variants require separate RefSeq-coordinate handling.

