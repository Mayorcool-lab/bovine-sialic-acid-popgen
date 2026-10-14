# Population Genomics of the Bovine Sialic Acid Network

Analysis code and supporting outputs for 13 selected genes in
346 cattle: 226 European animals from 12 breeds and 120 South African
animals from 4 breeds.

The study identified 62 protein-altering candidates and prioritized 15
using functional prediction, evolutionary conservation and protein
context, with population representation as supporting evidence.

## Data sources

- Ginja et al. (2025), Scientific Data 12, 1897.
  DOI: 10.1038/s41597-025-06188-x
  European variants: EVA PRJEB98152.
  Raw reads: ENA PRJEB90816.
- Dlamini et al. (2026), Scientific Data 13, 906.
  DOI: 10.1038/s41597-026-07458-y
  African variants: EVA PRJEB110435 and PRJEB94226.
  Raw reads: ENA PRJEB90914.

## Repository contents

- scripts/: analysis scripts.
- scripts/H_integration/: integrated evidence and prioritization.
- scripts/figures/: figure and supplementary-table preparation.
- data/metadata/: dataset and sample metadata.
- data/processed/: processed variants and intermediate files.
- data/external/aiyesa_popgen/: collaborator inputs and provenance.
- results/: analysis outputs and figures.
- envs/: environment definitions.

Historical files may describe earlier sample sets or intermediate
results. They should not be interpreted as final manuscript results.

## Environments and execution

envs/environment.yml specifies core analysis dependencies.
envs/sialic_plot_environment.yml records the plotting environment.
Neither file alone specifies every tool used across all analysis stages.

To recreate the plotting environment:

    conda env create -f envs/sialic_plot_environment.yml
    conda activate sialic_plot

Run scripts from the repository root. Scripts represent individual
stages, not a single automated workflow; inspect their input paths
before execution.

## Collaborator inputs

Leke V. Aiyesa supplied the initial population-genomic analyses.
See data/external/aiyesa_popgen/README.md for file roles and provenance.

Verify these inputs with:

    (cd data/external/aiyesa_popgen && sha256sum -c SHA256SUMS)

The original pop_gen_script.R is retained unchanged as an exploratory
analysis record, not a standalone final reproduction workflow.
Final Figure 8 uses validated H10 diversity estimates.

## Interpretation

Diversity is reported per ascertained pathway SNP, not per genomic base.
Differentiation and Tajima's D are interpreted descriptively.
Functional prioritization does not establish experimental effects.

## Contributions and citation

Oluwamayowa J. Ogun developed the pathway analyses, integrated evidence,
validated and revised analyses, and prepared manuscript outputs.
Leke V. Aiyesa performed the initial genome-wide population-genomic analyses,
including PCA, differentiation, diversity and Tajima's D, and supplied the
associated scripts and outputs. Oluwamayowa J. Ogun integrated the retained
population-genomic results with the pathway-level findings.
Full author contributions are described in the manuscript.

Citation metadata and a Zenodo DOI will accompany the archived release.
