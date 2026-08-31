# Phase H evidence integration design

## Purpose

Integrate the validated evidence for the same 62 protein-altering variants across:

- H2: South African versus European population-frequency evidence
- H3: 16-breed frequency evidence
- H4: VEP consequence and bovine SIFT evidence
- H5: evolutionary-conservation evidence
- H6: validated protein/domain/functional-context evidence

This stage synthesises existing evidence only. It does not generate new structural modelling, docking, molecular dynamics, or genome-wide population-genetic results.

## Canonical variant identity

Primary join key:

`CHROM + POS + REF + ALT`

`SYMBOL` is retained as an independent cross-check.

The join-key audit established:

- H2: 62 unique variants
- H3: 62 unique variants represented by 16 breed rows each
- H4: 62 unique variants
- H5: 62 unique variants
- H6: 62 unique variants
- no cross-layer allele-key mismatches
- no SYMBOL mismatches

Therefore the genomic allele key is authoritative for integration.

## Integration principle

No naive arithmetic composite score will be used.

Evidence dimensions remain separate so that a reviewer can see exactly why a candidate is or is not prioritised.

Missing evidence is never converted into negative evidence.

Overlapping annotations from related databases are not counted as independent biological votes.

## Evidence dimensions

### 1. Variant identity and molecular consequence

Retain:

- SYMBOL
- CHROM
- POS
- REF
- ALT
- HGVSp
- Consequence
- IMPACT

Stop-gained and stop-lost variants remain separate from missense-style interpretation.

### 2. H2 population-frequency evidence

Retain raw validated fields:

- SA_status
- SA_AF
- SA_N_genotyped
- SA_het
- SA_homalt
- EU_status
- EU_AF
- EU_N_genotyped
- EU_het
- EU_homalt

Rules:

- PRESENT permits frequency interpretation.
- ABSENT_AFTER_QC is retained explicitly and is not equivalent to a source-level absence.
- NOT_IN_SOURCE is retained explicitly and must never be converted to AF=0.
- No SA-versus-Europe difference will be calculated when either population lacks interpretable frequency data.

### 3. H3 breed-frequency evidence

The authoritative breed evidence remains in long format with 16 rows per variant.

Retain:

- COHORT
- BREED
- STATUS
- N_CALLABLE
- REF_COUNT
- ALT_COUNT
- AF
- N_FLAG

Rules:

- PRESENT permits breed-frequency interpretation.
- ABSENT_AFTER_QC remains distinct from source-level absence.
- NOT_DETECTED_IN_SOURCE remains explicit and is not treated as AF=0.
- VERY_SMALL_N is retained as a sample-size warning.
- Breed-level evidence will not be reduced to a single score.
- Any later breed-summary fields must remain traceable to the 16 underlying rows.

### 4. H4 functional-impact evidence

Retain:

- Consequence
- IMPACT
- SIFT_status
- SIFT_prediction
- SIFT_score

Observed states include:

- SIFT_AVAILABLE
- SIFT_UNAVAILABLE
- NOT_APPLICABLE_HIGH_IMPACT

Rules:

- SIFT_UNAVAILABLE is missing predictive evidence, not tolerated evidence.
- NOT_APPLICABLE_HIGH_IMPACT is interpreted according to the variant consequence rather than as missing missense evidence.
- deleterious_low_confidence remains distinct from deleterious.
- tolerated_low_confidence remains distinct from tolerated.
- SIFT prediction is evidence, not proof of biological effect.

### 5. H5 evolutionary-conservation evidence

Retain:

- coordinate_system
- reference_transcript
- reference_protein
- protein_position
- reference_aa
- comparative_total
- comparative_non_gap
- coverage_status
- rate4site_status
- rate4site_score
- rate4site_qq_interval
- rate4site_std

Rules:

- ASSESSABLE permits Rate4Site interpretation.
- NOT_ASSESSABLE is missing comparative evidence, not evidence of poor conservation.
- STOP_LOST_SPECIAL_CASE remains separate.
- Lower Rate4Site score means slower evolutionary rate / stronger conservation.
- No arbitrary high/moderate/low conservation threshold will be introduced.

### 6. H6 protein/domain/functional-context evidence

Retain the validated mapping and transfer-confidence fields:

- H5_coordinate_system
- H5_reference_protein
- H5_position
- H5_reference_aa
- external_database
- external_accession
- external_position
- external_aa
- mapping_status
- feature_transfer_status
- H6_2_branch
- context_level

Retain residue/domain-context evidence where available:

- uniprot_overlap
- uniprot_feature_types
- interpro_overlap
- interpro_entry_types
- interpro_source_databases
- interpro_entry_accessions
- ST3GAL4 reconciliation fields
- CMAH RefSeq X2 sequence-based fields
- protein_context
- context_evidence
- residue_specific_interpretation

Rules:

- RESIDUE_LEVEL_DIRECT_OR_MAPPED, RESIDUE_LEVEL_MULTI_REFERENCE_RECONCILED and RESIDUE_LEVEL_SEQUENCE_BASED remain distinguishable evidence types.
- PROTEIN_LEVEL_ONLY is lower residue-specific resolution, not absence of functional relevance.
- NO_DIRECT_UNIPROT_FEATURE_TRANSFER is a mapping/transfer constraint, not negative biological evidence.
- Missing feature overlap is not interpreted as functional neutrality.
- Domain/family overlap supplies context but does not prove that a variant changes function.
- Multiple InterPro/member-database hits describing overlapping hierarchy are not counted as independent evidence votes.

### 7. CMAH RefSeq X2 special rule

For:

- N386
- I388
- P424
- F512

the exact XP_024839792.1 sequence-based InterProScan audit established:

- 10 raw protein-level matches
- 7 signature libraries
- exactly one candidate-residue overlap per variant
- PANTHER PTHR46522 is the sole residue-overlapping signature
- no non-PANTHER narrower-domain overlap
- interpretation: BROAD_FAMILY_OVERLAP_ONLY

The complete protein-level match set remains available separately and must not be mistaken for ten residue-level overlaps.

### 8. Prioritisation principles

Candidate prioritisation will be based on transparent convergence among independent evidence dimensions, not score summation.

Factors that may support prioritisation include:

- biologically relevant molecular consequence
- damaging-direction functional prediction where assessable
- evolutionary constraint where assessable
- validated residue/domain/protein functional context
- informative population or breed segregation patterns
- suitability for downstream genotype-phenotype testing

Factors must remain individually visible.

A candidate will not be penalised merely because one evidence layer is structurally unavailable.

Extremely rare, fixed, source-absent, or poorly callable variants may have high mechanistic interest but lower practical informativeness for genotype-phenotype association; this will be considered explicitly rather than through an arbitrary AF threshold.

### 9. Genome-wide evidence

Genome-wide population structure, differentiation, and selection evidence generated independently by Aiyesa will be added later as a separate evidence dimension.

The pathway-side H2-H6 evidence matrix will therefore be complete and interpretable without waiting for genome-wide results.

## Planned authoritative outputs

1. One 62-row integrated candidate evidence matrix.
2. One retained 992-row breed-frequency companion table or traceable breed-summary representation.
3. QC proving one-to-one H2/H4/H5/H6 joins and 16 H3 breed rows per candidate.
4. A transparent candidate-prioritisation table generated only after the all-62 matrix has been inspected.
5. Approximately 5-10 final candidates for downstream genotype-informed early-lactation milk validation.

No candidate ranking is performed at the design stage.

## HGVSp cross-layer audit and integration rule

A cross-layer audit of Phase G, H2, H4 and H5 found:

- 62 shared genomic alleles
- 0 Consequence mismatches
- 58 variants with non-identical stored HGVSp values
- in those 58 cases, Phase G, H2 and H5 stored HGVSp as blank while H4 stored an amino-acid-change shorthand such as `R/C`, `T/S`, `E/*` or `*/L`

This pattern reflects annotation-field representation rather than disagreement in genomic variant identity or consequence.

Therefore:

- HGVSp is not part of the canonical integration key.
- HGVSp text is not required to be identical across phases.
- `CHROM + POS + REF + ALT` remains the authoritative candidate identity.
- SYMBOL and Consequence remain mandatory cross-layer consistency checks.
- For the integrated evidence matrix, the H4 amino-acid-change field is retained as the most informative validated protein-change representation when present.
- Missing HGVSp in Phase G, H2 or H5 is not interpreted as missing biological evidence because H5 independently retains validated protein position and reference amino acid.
- No upstream frozen file is modified to harmonise HGVSp retrospectively.

HGVSp representation differences therefore do not reopen H2, H4, H5 or H6.

## H3 candidate-level descriptive integration

The authoritative H3 evidence remains the lossless 992-row breed-frequency companion:

`results/H_integration/H3_62_candidate_16_breed_companion.tsv`

Each of the 62 candidates has exactly 16 breed rows:
- 4 South African breeds
- 12 European breeds

Breed-level evidence will not be converted into a composite score.

### Interpretation hierarchy

For each candidate and cohort, breed evidence is distinguished as:

1. source unavailable:
   `STATUS = NOT_DETECTED_IN_SOURCE`

2. source/site represented but removed after QC:
   `STATUS = ABSENT_AFTER_QC`

3. site available for breed-level interpretation:
   `STATUS = PRESENT`

4. ALT allele observed:
   `STATUS = PRESENT` and `ALT_COUNT > 0`

These states are not interchangeable.

`NOT_DETECTED_IN_SOURCE` is never interpreted as AF = 0.

`ABSENT_AFTER_QC` is never interpreted as AF = 0.

`PRESENT` with `ALT_COUNT = 0` represents an interrogated breed in which no ALT allele was observed among callable chromosomes.

### Sample-size handling

`N_FLAG = VERY_SMALL_N` is retained explicitly.

Mertolenga and Mirandesa contain systematically very-small-N observations in the present dataset. Their frequencies remain reportable and traceable but will not be allowed to drive adequately sampled breed-pattern summaries.

No new arbitrary sample-size threshold is introduced during integration. The validated H3 `N_FLAG` classification is retained.

### Candidate-level descriptive fields

For South Africa and Europe separately, candidate summaries may include:

- number of breed rows
- number with source unavailable
- number absent after QC
- number with STATUS = PRESENT
- number of PRESENT breeds with N_FLAG = OK
- number of adequately sampled breeds with ALT_COUNT > 0
- names of adequately sampled breeds with ALT_COUNT > 0
- AF values for those observed adequately sampled breeds
- minimum and maximum AF among those observed adequately sampled breeds
- breed or breeds carrying the maximum observed adequately sampled AF
- number of PRESENT breeds flagged VERY_SMALL_N
- names and AF values of observed VERY_SMALL_N breeds, retained separately

AF ranges are descriptive only and are not converted into evidence scores.

### Population-pattern interpretation boundary

Breed summaries may support descriptions such as:

- allele observed across multiple South African breeds
- allele observed across multiple European breeds
- allele observed in both population groups
- allele observed only in adequately sampled South African breeds within the available data
- allele observed only in adequately sampled European breeds within the available data

Such descriptions refer only to the analysed breed/source dataset.

They do not establish:
- continental fixation
- population specificity outside the sampled breeds
- selection
- adaptation
- functional causality
- association with milk sialic-acid phenotype

Aiyesa's independent genome-wide population-genetic analyses will later provide the appropriate broader population-structure, differentiation and selection context.

### Integration principle

The 992-row H3 companion remains authoritative for breed-level frequencies.

Any 62-row H3 summary is a descriptive convenience layer and must remain traceable to the underlying breed rows.

No candidate ranking is performed during H3 summarisation.

## Final candidate prioritisation framework

Candidate prioritisation is based on transparent convergence of independent evidence dimensions rather than a summed numerical score.

No arbitrary global cutoff is introduced for SIFT, Rate4Site, allele frequency or number of breeds.

### Evidence dimensions

Each candidate is evaluated across the following dimensions.

#### 1. Molecular consequence

High-impact stop-gained and stop-lost variants are considered separately from missense variants.

For missense variants, H4 SIFT is retained as one functional-prediction dimension.

- `deleterious` supports a functional-impact hypothesis
- `deleterious_low_confidence` is retained with its lower confidence
- `tolerated` does not prove functional neutrality
- `SIFT_UNAVAILABLE` is missing prediction evidence and is not interpreted as tolerated

#### 2. Evolutionary conservation

H5 Rate4Site is treated as a continuous phylogeny-aware conservation measure.

Lower Rate4Site scores indicate slower evolutionary change and stronger relative conservation.

No high/moderate/low conservation bins are created.

`NOT_ASSESSABLE` represents missing comparative evidence, not lack of conservation.

The stop-lost variant remains a separate H5 special case.

#### 3. Known protein/domain context

H6 provides the resolution of known curated functional/domain context.

Residue-level context can strengthen mechanistic plausibility but does not prove functional alteration.

`PROTEIN_LEVEL_ONLY` indicates lower residue-specific annotation resolution and is not interpreted as absence of biological relevance.

CMAH RefSeq-X2 sequence-based context remains distinct from direct UniProt residue transfer.

Overlapping InterPro/member-database annotations are not counted as independent evidence votes.

#### 4. Population and breed segregation

H2 population frequencies and H3 breed-level observations provide population-distribution and practical genotyping context.

Evidence is interpreted using validated source/QC states.

- `NOT_IN_SOURCE` and `NOT_DETECTED_IN_SOURCE` are never treated as AF = 0
- `ABSENT_AFTER_QC` is not treated as biological absence
- `PRESENT` with `ALT_COUNT = 0` indicates an interrogated cohort/breed with no ALT observed among callable chromosomes
- `VERY_SMALL_N` observations remain separate from adequately sampled breed evidence

Observation across multiple adequately sampled breeds can support practical genotype-phenotype follow-up.

Apparent population restriction is described only within the analysed dataset and is not interpreted as fixation, selection or adaptation.

#### 5. Aim 3 genotype-phenotype feasibility

A candidate must ultimately be useful for the proposed targeted-genotyping and early-lactation milk validation design.

Therefore, practical segregation is considered separately from molecular interest.

A variant may be highly mechanistically interesting but difficult to test in genotype-phenotype association if the ALT allele is not observed or is extremely restricted in accessible populations.

No allele-frequency threshold is imposed at this stage.

### Evidence convergence

Candidates are prioritised when multiple independent evidence dimensions support follow-up.

Examples of stronger convergence include:

- predicted functional impact plus evolutionary conservation
- predicted functional impact plus residue-level functional/domain context
- evolutionary conservation plus residue-level functional/domain context
- molecular evidence plus robust segregation in adequately sampled breeds
- high-impact molecular consequence plus practical segregation

Convergence does not require every evidence layer to be available.

Missing evidence is kept distinct from contradictory evidence.

### Contradictory evidence

Discordant evidence is retained explicitly rather than averaged away.

Examples include:

- deleterious SIFT with comparatively fast-evolving Rate4Site position
- strong conservation with tolerated SIFT
- strong molecular evidence but no ALT observed in adequately sampled breeds
- broad population segregation but weak residue-specific functional evidence

Such variants may remain biologically interesting but require a transparent rationale for inclusion or exclusion.

### Prioritisation outputs

The final prioritisation table will retain all 62 candidates and record:

- molecular-consequence evidence
- SIFT evidence
- conservation evidence
- H6 context resolution
- South African breed observation
- European breed observation
- practical genotype-phenotype informativeness
- evidence-convergence rationale
- evidence limitations or contradictions
- prioritisation decision

No arithmetic composite score will be generated.

The final shortlist is expected to contain approximately 5-10 candidates, selected for complementary biological and population-genetic evidence rather than numerical rank alone.

Aiyesa's genome-wide population-structure, differentiation and selection results will later be added as a separate independent evidence dimension and will not overwrite the pathway-side evidence.

### CMAH RefSeq X2 post-closure propagation repair

A post-closure audit identified a narrow representation defect in the authoritative H6.2 table for the four CMAH RefSeq X2 candidates (N386S, I388M, P424L and F512Y). The validated InterProScan 5.78-109.0 evidence was already preserved in dedicated CMAH X2 fields, but had not been propagated into the generic `protein_context`, `context_evidence` and `residue_specific_interpretation` fields. This was not a sequence-mapping, InterProScan, candidate-identity or integration-join error.

A guarded repair populated only those three generic fields for the four affected candidates using the validated CMAH X2 audit evidence: PANTHER PTHR46522 / IPR027033, interval 5-576, with `BROAD_FAMILY_OVERLAP_ONLY` interpretation and no narrower non-PANTHER signature overlapping the candidate residues. Column-aware comparison against the previously committed H6.2 table confirmed exactly 4 affected candidates and 12 changed cells, with no other H6.2 fields altered.

Post-repair validation confirmed correct propagation through the H2-H6 evidence matrix, complete H2-H6-H3 matrix and candidate-prioritisation review. Rebuilt factual evidence flags were unchanged, demonstrating that the repair restored descriptive H6 context without altering the underlying population, breed, SIFT, conservation or evidence-convergence landscape.
