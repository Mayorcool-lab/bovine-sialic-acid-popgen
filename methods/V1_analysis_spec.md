# V1 Analysis Specification — sialic_acid_popgen

Provisional, revisable. Advisory input from Galuska (gene panel) and Aiyesa
(windows/methods) may still adjust this; nothing here is final.

## Gene panel (V1, written biological rationale)

| Biological component | Genes |
|---|---|
| Sialic-acid biosynthesis | GNE, NANS, NANP |
| CMP-sialic-acid activation/transport | CMAS, SLC35A1 |
| Neu5Gc production | CMAH |
| Sialylation | ST6GAL1, ST3GAL1, ST3GAL4, ST8SIA2, ST8SIA4 |
| Desialylation | NEU1, NEU3 |

Hypothesis-driven pathway panel covering every major step of sialic acid
metabolism. NOT an arbitrary gene list. Coordinates (release 110, ARS-UCD1.2)
in data/reference/candidate_gene_coordinates_provisional.tsv, all 13 of 13
resolved and verified.

## V1 genomic analysis scopes

1. Functional variant discovery
   - Current provisional pathway genes (13, above)
   - Gene body + 2 kb strand-aware upstream sequence
   - VEP + seven-tool functional consensus (verify cattle-applicability of
     each tool before treating results as final)

2. Population structure
   - Genome-wide harmonised, LD-pruned SNPs
   - PCA and ADMIXTURE, run BEFORE defining final population contrasts (17f)

3. FST / nucleotide diversity / Tajima's D
   - Genome-wide harmonised SNP dataset, NOT restricted to gene ±1Mb
   - Provisional sliding-window spec: 50 kb window, 20 kb step
   - Verified citations: Frontiers in Genetics (2022), Dengchuan cattle,
     https://doi.org/10.3389/fgene.2022.833475; and BMC Genomics (2026),
     Guanling cattle, https://doi.org/10.1186/s12864-025-12480-2
   - Candidate pathway genes interrogated against the genome-wide empirical
     distribution, not interpreted in isolation

4. iHS / XP-EHH
   - REQUIRES phased, chromosome-scale (not gene-window) data — mathematical
     requirement of the statistic itself, not a design choice
   - selscan-based computation, genome-wide normalization
   - Candidate pathway loci evaluated AFTER genome/chromosome-wide computation
   - selscan default max EHH extension = 1 Mb — internal computational
     parameter, NOT justification for restricting input data
   - PRACTICAL IMPLICATION: Aiyesa needs full-chromosome phased data for
     whichever chromosomes the 13 genes sit on, not just small extracts

5. Interpretation principle
   - Selection evidence strongest when multiple independent statistics
     converge AND the signal is extreme relative to genomic background

## References (verified directly)

1. Kim, H.J., Moghaddar, N., Clark, S., van der Werf, J.H.J., de Las
   Heras-Saldana, S. (2025). Genome-wide scan for signatures of selection in
   Hanwoo and Angus cattle using whole-genome sequence data. PLOS ONE, 20(5),
   e0324034. https://doi.org/10.1371/journal.pone.0324034
   USE FOR: post-hoc candidate-gene annotation only. Do NOT cite as
   justification for restricting FST/iHS/XP-EHH computation itself.

2. Liu, D., Chen, Z., Zhao, W., Guo, L., Sun, H., Zhu, K., Liu, G., Shen, X.,
   Zhao, X., Wang, Q., Ma, P., Pan, Y. (2021). Genome-wide selection
   signatures detection in Shanghai Holstein cattle population identified
   genes related to adaption, health and reproduction traits. BMC Genomics,
   22(1), 747. https://doi.org/10.1186/s12864-021-08042-x
   FINDING: ST3GAL4 independently identified via |iHS| and ROH (NOT XP-EHH).
   CORRECT: "Previous cattle selection-signature analysis independently
   identified ST3GAL4 within a candidate region under selection (Liu et al.,
   2021)."

3. Frontiers in Genetics (2022). Whole-Genome Sequencing of Endangered
   Dengchuan Cattle Reveals Its Genomic Diversity and Selection Signatures.
   https://doi.org/10.3389/fgene.2022.833475
   USE FOR: 50kb/20kb sliding window precedent for nucleotide diversity.

4. BMC Genomics (2026). Genetic admixture and adaptive signatures of
   Guanling cattle revealed by whole-genome sequence.
   https://doi.org/10.1186/s12864-025-12480-2
   USE FOR: 50kb/20kb sliding window precedent for FST/XP-EHH.

## Still open
- [ ] Confirm which chromosomes need full-phased extraction for iHS/XP-EHH
