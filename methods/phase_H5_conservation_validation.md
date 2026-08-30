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


## H5 orthologue-source and sampling validation

### Ensembl comparative-genomics release

The bovine variant and protein-coordinate framework remains anchored to Ensembl release 110 (ARS-UCD1.2), consistent with Phase G/H4 annotation and the release-matched bovine peptide FASTA.

The current Ensembl REST server reported:

- Ensembl REST release: 116

A test of the historical May 2023 REST endpoint returned HTTP 503 and was therefore not used for production comparative-genomics retrieval.

Accordingly, Ensembl release 110 remains the authoritative bovine coordinate/protein reference, whereas comparative orthology will be treated as a separately versioned evidence layer obtained from Ensembl Compara release 116. Exact orthologue sequences and metadata used for conservation analysis must be frozen locally to prevent future database updates from silently changing the analysis.

### ST6GAL1 orthologue retrieval pilot

ST6GAL1 (ENSBTAG00000010357) was used as the initial retrieval test.

The Ensembl REST homology endpoint returned:

- Total orthologue records: 254
- ortholog_one2one: 197
- ortholog_one2many: 57
- Unique one-to-one species labels: 197

Expected mammalian orthologues including sheep, goat, pig, horse, dog, human, and mouse were recovered with plausible protein lengths.

However, the returned set contained numerous Sus scrofa breed/assembly-specific labels. Therefore, Ensembl species labels cannot be treated automatically as independent evolutionary samples, and all returned orthologues must not be used without taxonomic filtering.

### Cross-gene one-to-one orthologue availability

One-to-one orthologue availability was evaluated across all 13 pathway genes.

Per-gene numbers of one-to-one species labels were:

- CMAH: 146
- CMAS: 197
- GNE: 205
- NANP: 204
- NANS: 185
- NEU1: 179
- NEU3: 154
- SLC35A1: 199
- ST3GAL1: 167
- ST3GAL4: 195
- ST6GAL1: 197
- ST8SIA2: 206
- ST8SIA4: 171

A total of 73 one-to-one species labels were shared across all 13 genes. This intersection was nevertheless taxonomically unbalanced and included multiple Sus scrofa assemblies/breeds.

Therefore, simple intersection of all available one-to-one orthologues will not be used as the production conservation panel.

### Balanced-panel feasibility test

A preliminary phylogenetically distributed species set was examined to determine whether clean one-to-one orthologues are broadly available across the pathway.

Complete 13/13 one-to-one coverage was observed for:

- Ovis aries
- Bison bison bison
- Bos grunniens
- Camelus dromedarius
- Sus scrofa
- Felis catus
- Mus musculus
- Rattus norvegicus
- Loxodonta africana
- Vombatus ursinus
- Xenopus tropicalis

Near-complete coverage included:

- Capra hircus: 12/13; CMAH not classified one-to-one
- Homo sapiens: 12/13; CMAH not classified one-to-one
- Canis lupus familiaris: 12/13; NANP not classified one-to-one
- Anolis carolinensis: 12/13; NEU3 not classified one-to-one
- Equus caballus: 11/13; CMAS and NEU1 not classified one-to-one
- Gallus gallus: 11/13; CMAH and NEU1 not classified one-to-one
- Danio rerio: 9/13; CMAS, NANS, NEU3, and ST3GAL1 not classified one-to-one

These results support use of a predefined phylogenetically balanced species framework with gene-specific missingness rather than requiring an identical orthologue set for every gene.

### Goat CMAH safeguard

Failure of Capra hircus CMAH to enter the strict one-to-one set does not indicate absence of CMAH in goat.

Ensembl Compara release 116 returned two Capra hircus CMAH-related models, both classified as ortholog_one2many relative to bovine CMAH:

- ENSCHIG00000016421: 333 aa; 71.4715% identity; 76.5766% positive similarity
- ENSCHIG00000011804: 237 aa; 95.7806% identity; 98.3122% positive similarity

Therefore, absence from a strict ortholog_one2one filter must not be interpreted as biological gene absence. One-to-many cases require explicit evaluation before inclusion in conservation alignments.

### Orthologue-design decision at this checkpoint

No conservation scores have yet been calculated.

Production conservation analysis will require:

1. a predefined phylogenetically balanced species-selection rule;
2. explicit handling of one-to-many and genuinely missing orthologues;
3. avoidance of redundant breed/assembly representations;
4. freezing of the exact Ensembl Compara release 116 sequences and metadata used;
5. preservation of the validated Ensembl 110 bovine residue-coordinate system; and
6. separate RefSeq-coordinate handling for the four manually rescued CMAH variants.

### Production taxonomic framework

A fixed 18-species phylogenetically distributed target framework was selected before examination of candidate-site conservation scores:

- Ovis aries
- Capra hircus
- Bison bison bison
- Bos grunniens
- Camelus dromedarius
- Sus scrofa
- Equus caballus
- Canis lupus familiaris
- Felis catus
- Homo sapiens
- Mus musculus
- Rattus norvegicus
- Loxodonta africana
- Vombatus ursinus
- Gallus gallus
- Anolis carolinensis
- Xenopus tropicalis
- Danio rerio

The framework intentionally spans close bovids/ruminants, other livestock and Laurasiatherian mammals, primates and rodents, broader mammalian lineages, and increasingly distant vertebrates.

The target framework is fixed independently of candidate-site conservation results. Gene-specific missingness is permitted and will be reported rather than forcing an orthologue assignment.

One-to-one Ensembl Compara orthologues will be used automatically where available. One-to-many relationships will not be selected automatically and will instead be flagged for explicit gene-specific evaluation. Redundant breed- or assembly-specific representations will not be added when the designated reference species is available.

Bos taurus itself will serve as the validated reference sequence and is not counted among the 18 comparative target species.

### Pilot sequence-completeness safeguard

Production retrieval was piloted with ST6GAL1 across the fixed 18-species
comparative framework. All 18 target species were returned by Ensembl
Compara release 116 as ortholog_one2one relationships.

Sequence-length inspection identified one important exception: the
Rattus norvegicus ST6GAL1 orthologue (ENSRNOG00000001823;
ENSRNOP00000041089) was represented by a 216-aa protein, compared with
the validated 405-aa bovine ST6GAL1 reference. The returned rat sequence
showed 72.22% target identity and 82.87% target positive similarity over
the represented sequence, supporting orthology but indicating that the
specific protein model is substantially shorter than the bovine
reference.

Therefore, orthology status and sequence completeness will be treated
as separate QC dimensions. An ortholog_one2one relationship alone will
not be sufficient for automatic inclusion in residue-level conservation
analysis. Protein-length/completeness checks will be applied
gene-specifically, and substantially shortened sequences will be flagged
for review rather than silently accepted or discarded.

This safeguard is particularly important for genes with complex
annotation histories, including CMAH.

### Ensembl release/transcript compatibility safeguard

The bovine protein used by Ensembl Compara release 116 was compared with
the exact Ensembl release-110 protein previously validated as the
candidate-coordinate reference for each of the 13 genes.

Only CMAS and NANS matched both protein identifier and protein length
exactly between the release-110 candidate reference and the Compara-116
source protein. Several genes showed substantial source-protein changes,
including CMAH (383 aa in the validated release-110 reference versus
587 aa in Compara 116), NANP (248 versus 306 aa), ST3GAL4 (395 versus
425 aa), GNE (722 versus 752 aa), and NEU3 (428 versus 472 aa).
ST6GAL1, ST8SIA2 and ST8SIA4 had different protein identifiers but the
same protein length between the two reference layers.

Consequently, Ensembl Compara 116 will be used to identify and retrieve
comparative orthologues, but its bovine source protein and precomputed
alignment will not define candidate residue coordinates. Candidate-site
mapping will remain anchored to the exact Ensembl release-110 bovine
protein validated against Phase G.

For each gene, selected comparative protein sequences will therefore be
combined with the exact release-110 bovine reference protein and aligned
locally using MAFFT. Conservation at candidate sites will subsequently
be mapped through the release-110 bovine sequence.

Length ratios calculated against the Compara-116 bovine source protein
are retained as retrieval metadata but will not, by themselves, be used
as sequence-exclusion criteria because source isoform differences
between Ensembl releases can create misleading apparent truncation.

### Production orthologue retrieval and sequence QC

The fixed 18-species framework was retrieved for all 13 genes using
Ensembl Compara release 116, giving 234 expected gene-by-species
relationships. Classification of the returned relationships produced:

- 222 ortholog_one2one relationships;
- 8 ambiguous relationships containing one-to-many orthology; and
- 4 target gene/species combinations not returned by the Compara query.

All 222 one-to-one target protein sequences were frozen locally from the
ungapped target `align_seq` returned by Ensembl Compara. Cross-checking
the per-gene FASTA files against the combined retrieval metadata
confirmed 222 expected and 222 observed sequences, with no missing or
unexpected records and no protein-identifier or sequence-length
mismatches.

The eight one-to-many cases were examined individually. No relationship
was resolved solely by choosing the sequence with the highest identity,
similarity, or protein length. These cases therefore remain AMBIGUOUS
unless independent annotation evidence later provides a unique
biological basis for selecting an orthologue. Similarly, NOT_RETURNED
denotes absence from this specific Compara retrieval and is not
interpreted as evidence of biological gene absence.

Sequence-completeness screening was repeated using the exact Ensembl
release-110 bovine proteins as the denominator rather than the
Compara-116 bovine source isoforms. Of the 222 one-to-one target
sequences, 220 had a target/reference length ratio >=0.80. Only two
sequences were substantially shorter than their release-110 bovine
references:

- Rattus norvegicus ST6GAL1: 216 aa versus 405 aa (ratio 0.5333);
- Vombatus ursinus NEU1: 262 aa versus 416 aa (ratio 0.6298).

Gene-model inspection showed that the rat ST6GAL1 gene
(ENSRNOG00000001823) contains two protein-coding translations of 214
and 216 aa, with the 216-aa Compara protein ENSRNOP00000041089 being
canonical. The wombat NEU1 gene (ENSVURG00010002968) contains a single
protein-coding transcript producing the 262-aa Compara protein
ENSVURP00010003709. Thus, neither short Compara sequence can be replaced
by a clearly longer alternative Ensembl translation from the same gene
model.

The relevant bovine candidate positions are 28, 153, 165 and 173 for
ST6GAL1, and 8, 13, 254, 275, 294, 380, 381, 382 and 400 for NEU1.
Consequently, overall protein length will not be used as an arbitrary
sequence-exclusion threshold. The two shortened models will remain
explicitly flagged, and their usability will be evaluated from local
alignment coverage and alignment quality at each candidate residue.

No candidate-site conservation scores have been calculated at this
checkpoint.
