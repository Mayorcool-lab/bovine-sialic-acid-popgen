# H9: Candidate-level pairwise breed differentiation

## Proposed manuscript Methods text

Pairwise genetic differentiation was evaluated among 16 cattle breeds represented by 346 animals at 62 protein-altering candidate variants. Breed membership was established using a verified sample manifest matched to the sample identifiers in the European and South African quality-controlled VCF datasets. Candidate genotypes were extracted from the corresponding gene-specific VCF files, and all 120 unique breed pairs were evaluated using the Weir–Cockerham FST estimator implemented in VCFtools. At least two called diploid genotypes were required per breed for each candidate comparison. Finite estimates were cross-checked against an independent implementation of the estimator, with a maximum permitted numerical difference of 1e-5. Negative FST estimates were retained. Comparisons with undefined estimates, insufficient genotype calls, or candidate variants absent from the relevant cohort VCF were explicitly classified without assigning artificial zero values. Breed pairs involving fewer than five animals in either breed were flagged for cautious interpretation. The candidate-level estimates and associated genotype statistics are provided in Supplementary Table S4.

## Reproducibility

- Analysis: scripts/H9_candidate_breed_fst.py
- Results: results/H_integration/H9_candidate_breed_fst/candidate_breed_pairwise_fst_62.tsv
- Manifest: results/H_integration/H9_candidate_breed_fst/breed_sample_manifest_346.tsv
- S4 generator: scripts/figures/generate_supplementary_table_s4.py
- S4 table: results/figures/source_data/Supplementary_Table_S4_candidate_breed_fst.tsv

## Validation summary

- Animals: 346
- Breeds: 16
- Candidate variants: 62
- Breed pairs: 120
- Records: 7440
- Finite FST estimates: 2028
- Undefined reference-fixed comparisons: 2463
- Absent breed A cohort VCF records: 2168
- Absent breed B cohort VCF records: 766
- Insufficient called genotypes: 15
- Small-breed comparison records: 1798
- Maximum estimator difference: 4.98e-7

## Manuscript integration

This text is saved for later integration into the existing manuscript Methods section. Check for duplication and confirm Supplementary Table S4 numbering before insertion. Candidate-level FST must not be described as genome-wide FST.
