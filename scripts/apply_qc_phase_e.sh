#!/bin/bash
set -euo pipefail

genes="CMAH GNE NANS NANP CMAS SLC35A1 ST8SIA2 ST8SIA4 ST6GAL1 ST3GAL1 ST3GAL4 NEU1 NEU3"
mkdir -p data/processed/qc_filtered/south_africa
mkdir -p data/processed/qc_filtered/europe

for cohort in south_africa europe; do
    if [ "$cohort" = "south_africa" ]; then
        indir="data/processed/snps/south_africa"
    else
        indir="data/processed/snps/europe"
    fi
    outdir="data/processed/qc_filtered/${cohort}"

    echo ">>> $cohort"
    for gene in $genes; do
        infile="${indir}/${gene}_biallelic_snps.vcf.gz"
        outfile="${outdir}/${gene}_qc.vcf.gz"

        # Step 1: set genotypes below DP4/GQ20 to missing (not remove the site,
        # just mark the individual genotype call as unreliable)
        bcftools filter -S . -e 'FMT/DP<4 | FMT/GQ<20' -O z -o "${outdir}/${gene}_gtmasked.vcf.gz" "$infile"
        bcftools index "${outdir}/${gene}_gtmasked.vcf.gz"

        # Step 2: NOW remove sites where, after masking, call rate drops below 90%
        n_samples=$(bcftools query -l "${outdir}/${gene}_gtmasked.vcf.gz" | wc -l)
        bcftools view -i "F_MISSING < 0.10" -O z -o "$outfile" "${outdir}/${gene}_gtmasked.vcf.gz"
        bcftools index "$outfile"

        before=$(bcftools view -H "$infile" | wc -l)
        after=$(bcftools view -H "$outfile" | wc -l)
        echo "  $gene: $before -> $after sites after QC"

        rm "${outdir}/${gene}_gtmasked.vcf.gz" "${outdir}/${gene}_gtmasked.vcf.gz.csi"
    done
done
