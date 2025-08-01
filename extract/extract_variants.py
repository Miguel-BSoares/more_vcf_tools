#!/usr/bin/env python3

import sys
import pysam

def parse_coordinates(file_path):
    coords = []
    with open(file_path) as f:
        for line in f:
            if line.strip():
                chrom, pos = line.strip().split()
                coords.append((chrom, int(pos)))
    return coords

def extract_variants(vcf_path, coords, out_prefix, flank_left, flank_right):
    vcf = pysam.VariantFile(vcf_path)

    extracted_vcf_path = f"{out_prefix}_extracted.vcf"
    allele_dist_path = f"{out_prefix}_allele_distribution.tsv"
    no_snp_path = f"{out_prefix}_no_snps.txt"

    extracted_vcf = pysam.VariantFile(extracted_vcf_path, 'w', header=vcf.header)
    allele_file = open(allele_dist_path, 'w')
    no_snp_file = open(no_snp_path, 'w')

    allele_file.write("CHROM\tPOS\tREF\tALT\tAF\tSample_GTs\n")

    for chrom, pos in coords:
        found = False
        try:
            for record in vcf.fetch(chrom, pos + flank_left, pos + flank_right + 1):
                found = True
                extracted_vcf.write(record)

                afs = record.info.get('AF')
                if afs is not None:
                    af_val = afs[0] if isinstance(afs, (list, tuple)) else afs
                    try:
                        allele_freq = round(float(af_val), 3)
                    except Exception:
                        allele_freq = 'NA'
                else:
                    allele_freq = 'NA'

                gts = []
                for sample in record.samples:
                    gt = record.samples[sample].get('GT')
                    if gt:
                        gts.append("/".join(str(x) for x in gt))
                    else:
                        gts.append("./.")

                allele_file.write(
                    f"{record.chrom}\t{record.pos}\t{record.ref}\t{','.join(str(a) for a in record.alts)}"
                    f"\t{allele_freq}\t{','.join(gts)}\n"
                )

            if not found:
                no_snp_file.write(f"{chrom}\t{pos}\n")
        except ValueError as e:
            no_snp_file.write(f"{chrom}\t{pos}\t# ERROR: {str(e)}\n")

    allele_file.close()
    no_snp_file.close()
    extracted_vcf.close()
    vcf.close()

def main():
    if len(sys.argv) != 6:
        print("Usage: python extract_variants.py <vcf_file> <coordinate_file> <output_prefix> <left_flank> <right_flank>")
        sys.exit(1)

    vcf_path = sys.argv[1]
    coord_path = sys.argv[2]
    out_prefix = sys.argv[3]
    flank_left = int(sys.argv[4])
    flank_right = int(sys.argv[5])

    coords = parse_coordinates(coord_path)
    extract_variants(vcf_path, coords, out_prefix, flank_left, flank_right)

if __name__ == "__main__":
    main()
