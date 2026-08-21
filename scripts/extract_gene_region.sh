#!/bin/bash 
set -euo pipefail
# Usage: ./extract_gene_region.sh <gene_symbol> <vcf_url> <output_dir>
# Reads coordinates from data/reference/candidate_gene_coordinates_provisional.tsv
# Adds strand-aware 2kb upstream buffer for functional-annotation extraction

GENE=$1
VCF_URL=$2
OUTDIR=$3
UPSTREAM_BP=2000

COORDS=$(grep -P "^${GENE}\t" data/reference/candidate_gene_coordinates_provisional.tsv)

if [ -z "$COORDS" ]; then
    echo "ERROR: Gene $GENE not found in coordinates file"
    exit 1
fi

CHR=$(echo "$COORDS" | cut -f3)
START=$(echo "$COORDS" | cut -f4)
END=$(echo "$COORDS" | cut -f5)
STRAND=$(echo "$COORDS" | cut -f6)

# Strand-aware: upstream means BEFORE the gene in its own reading direction
if [ "$STRAND" == "+" ]; then
    EXT_START=$((START - UPSTREAM_BP))
    EXT_END=$END
else
    EXT_START=$START
    EXT_END=$((END + UPSTREAM_BP))
fi

if [ "$EXT_START" -lt 1 ]; then
    EXT_START=1
fi

echo "Gene: $GENE | Strand: $STRAND | Gene body: ${CHR}:${START}-${END} | Extraction window (with 2kb upstream): ${CHR}:${EXT_START}-${EXT_END}"

mkdir -p "$OUTDIR"
bcftools view -r "${CHR}:${EXT_START}-${EXT_END}" "$VCF_URL" -O z -o "${OUTDIR}/${GENE}.vcf.gz"
bcftools index "${OUTDIR}/${GENE}.vcf.gz"

echo "Done: ${OUTDIR}/${GENE}.vcf.gz"
bcftools stats "${OUTDIR}/${GENE}.vcf.gz" | grep "number of records:"
