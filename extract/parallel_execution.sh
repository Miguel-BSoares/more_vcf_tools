#!/bin/bash

set -euo pipefail

usage() {
    cat << EOF
Usage: $0 -v <input_vcf> -c <coord_file> -o <output_prefix> -l <left_flank> -r <right_flank> -t <threads> [-s <chunk_size>] [-h]

Options:
  -v   Input VCF file (bgzipped and indexed)
  -c   Coordinates file (chrom pos)
  -o   Output prefix
  -l   Left flank coordinate (e.g. -10000)
  -r   Right flank coordinate (e.g. 10000)
  -t   Number of parallel threads to run
  -s   Chunk size: number of variants per chunk (default: 100000)
  -h   Show this help message and exit

Example:
  $0 -v your_large.vcf.gz -c snps.txt -o variants -l -10000 -r 10000 -t 8 -s 50000
EOF
    exit 0
}

# Defaults
CHUNK_SIZE=100000

while getopts ":v:c:o:l:r:t:s:h" opt; do
  case $opt in
    v) VCF="$OPTARG" ;;
    c) COORD="$OPTARG" ;;
    o) OUTPREFIX="$OPTARG" ;;
    l) LEFT_FLANK="$OPTARG" ;;
    r) RIGHT_FLANK="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    s) CHUNK_SIZE="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Validate required args
if [ -z "${VCF-}" ] || [ -z "${COORD-}" ] || [ -z "${OUTPREFIX-}" ] || [ -z "${LEFT_FLANK-}" ] || [ -z "${RIGHT_FLANK-}" ] || [ -z "${THREADS-}" ]; then
    echo "Error: Missing required arguments." >&2
    usage
fi

CHUNKS_DIR=chunks
mkdir -p "$CHUNKS_DIR"

echo "=== Step 1: Splitting VCF into chunks of $CHUNK_SIZE variants ==="
bcftools split -c "$CHUNK_SIZE" -Oz -o "$CHUNKS_DIR/${OUTPREFIX}_chunk_%03d.vcf.gz" "$VCF"

echo "Indexing chunked VCF files..."
for f in "$CHUNKS_DIR/${OUTPREFIX}_chunk_"*.vcf.gz; do
    tabix -p vcf "$f"
done

echo "=== Step 2: Running parallel extraction on chunks ==="
ls "$CHUNKS_DIR/${OUTPREFIX}_chunk_"*.vcf.gz | parallel -j "$THREADS" \
  python extract_variants_cyvcf2.py {} "$COORD" {}_output "$LEFT_FLANK" "$RIGHT_FLANK"

echo "=== Step 3: Merging outputs ==="

bcftools concat -Oz -o "${OUTPREFIX}_all_extracted.vcf.gz" "$CHUNKS_DIR/"*_output_extracted.vcf
tabix -p vcf "${OUTPREFIX}_all_extracted.vcf.gz"

head -1 "$CHUNKS_DIR/${OUTPREFIX}_chunk_000_output_allele_distribution.tsv" > "${OUTPREFIX}_all_allele_distribution.tsv"
tail -q -n +2 "$CHUNKS_DIR/"*_output_allele_distribution.tsv >> "${OUTPREFIX}_all_allele_distribution.tsv"

cat "$CHUNKS_DIR/"*_output_no_snps.txt > "${OUTPREFIX}_all_no_snps.txt"

echo "Pipeline complete!"
