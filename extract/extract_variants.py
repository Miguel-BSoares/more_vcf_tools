#!/usr/bin/env python3
import sys
import os
from cyvcf2 import VCF, Writer

def parse_coords(coord_file, left_flank, right_flank):
    """
    Parse coordinate file and return list of (chrom, start, end, coord_str) tuples.
    coord_str is for reporting (e.g. "chr1:12345")
    """
    regions = []
    with open(coord_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            start = max(1, pos + int(left_flank))  # ensure start >= 1
            end = pos + int(right_flank)
            if end < start:
                start, end = end, start
            coord_str = f"{chrom}:{pos}"
            regions.append((chrom, start, end, coord_str))
    return regions

def get_genotype_code(gt):
    """
    gt: tuple/list of two alleles, e.g. (0,0), (0,1), (1,1)
    Returns:
        0 for homozygous reference (0/0)
        1 for homozygous alternative (1/1)
        2 for heterozygous (0/1 or 1/0)
        -1 if missing or other
    """
    if gt is None or len(gt) != 2:
        return -1
    a1, a2 = gt
    if a1 is None or a2 is None:
        return -1
    if a1 == 0 and a2 == 0:
        return 0
    if a1 == 1 and a2 == 1:
        return 1
    if (a1 == 0 and a2 == 1) or (a1 == 1 and a2 == 0):
        return 2
    return -1

def main():
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} <vcf_file> <coord_file> <output_prefix> <left_flank> <right_flank>", file=sys.stderr)
        sys.exit(1)

    vcf_file = sys.argv[1]
    coord_file = sys.argv[2]
    output_prefix = sys.argv[3]
    left_flank = int(sys.argv[4])
    right_flank = int(sys.argv[5])

    regions = parse_coords(coord_file, left_flank, right_flank)
    if not regions:
        print("No valid regions found in coordinate file.", file=sys.stderr)
        sys.exit(1)

    vcf = VCF(vcf_file)
    # Prepare output VCF writer, using the header from input VCF
    output_vcf_path = f"{output_prefix}_extracted.vcf"
    w = Writer(output_vcf_path, vcf)

    allele_dist_path = f"{output_prefix}_allele_distribution.tsv"
    no_snps_path = f"{output_prefix}_no_snps.txt"

    # Prepare allele distribution header
    # We'll output: CHROM, POS, ID (if available), REF, ALT, AlleleFreqRef, AlleleFreqAlt, GenotypeCodes (sample-wise)
    samples = vcf.samples
    with open(allele_dist_path, 'w') as ad_out, open(no_snps_path, 'w') as no_snp_out:
        header_cols = ["CHROM", "POS", "ID", "REF", "ALT", "AlleleFreqRef", "AlleleFreqAlt"] + [f"{s}_GTcode" for s in samples]
        ad_out.write("\t".join(header_cols) + "\n")

        for chrom, start, end, coord_str in regions:
            # Check if chrom exists in VCF contigs, skip if not
            if chrom not in vcf.seqnames:
                no_snp_out.write(f"{coord_str}\tchromosome_not_found\n")
                continue

            found_variants = False
            # cyvcf2 returns 0-based start, so region is start-1 to end (1-based inclusive)
            for variant in vcf(f"{chrom}:{start}-{end}"):
                found_variants = True
                # Safely get variant info
                chrom_v = variant.CHROM if hasattr(variant, "CHROM") else "."
                pos_v = variant.POS if hasattr(variant, "POS") else "."
                id_v = variant.ID if (hasattr(variant, "ID") and variant.ID is not None) else "."
                ref_v = variant.REF if hasattr(variant, "REF") else "."
                alt_v = ",".join(variant.ALT) if hasattr(variant, "ALT") else "."

                # Allele counts and freq
                # Allele counts are total counts of alleles observed in samples:
                # Use variant.format('AD') if available or calculate from genotypes
                # cyvcf2 provides variant.num_alleles etc.

                # We can use variant.format("AF") if present; otherwise compute from genotypes
                afs = variant.INFO.get("AF")
                if afs is not None and len(afs) > 0:
                    af = afs[0]
                    allele_freq_alt = af
                    allele_freq_ref = 1 - af
                else:
                    # calculate allele freq manually from genotypes
                    allele_counts = [0, 0]  # ref, alt
                    for gt in variant.genotypes:  # gt: [allele1, allele2, phased]
                        a1, a2, phased = gt
                        if a1 is not None and a1 >= 0:
                            if a1 == 0:
                                allele_counts[0] += 1
                            elif a1 == 1:
                                allele_counts[1] += 1
                        if a2 is not None and a2 >= 0:
                            if a2 == 0:
                                allele_counts[0] += 1
                            elif a2 == 1:
                                allele_counts[1] += 1
                    total_alleles = sum(allele_counts)
                    if total_alleles > 0:
                        allele_freq_ref = allele_counts[0] / total_alleles
                        allele_freq_alt = allele_counts[1] / total_alleles
                    else:
                        allele_freq_ref = 0.0
                        allele_freq_alt = 0.0

                # Genotype codes per sample
                gt_codes = []
                for gt in variant.genotypes:
                    code = get_genotype_code(gt[:2])  # ignore phased flag
                    gt_codes.append(str(code))

                # Write variant to output VCF
                w.write_record(variant)

                # Write allele distribution line
                out_fields = [chrom_v, str(pos_v), id_v, ref_v, alt_v,
                              f"{allele_freq_ref:.5f}", f"{allele_freq_alt:.5f}"] + gt_codes
                ad_out.write("\t".join(out_fields) + "\n")

            if not found_variants:
                no_snp_out.write(f"{coord_str}\n")

    w.close()
    vcf.close()

if __name__ == "__main__":
    main()
