#!/usr/bin/env python3

import argparse
import gzip
from collections import defaultdict
from intervaltree import IntervalTree
from cyvcf2 import VCF
from pyfaidx import Fasta
from Bio.Seq import Seq

def parse_gff(gff_file):
    features = defaultdict(lambda: defaultdict(IntervalTree))
    with gzip.open(gff_file, 'rt') if gff_file.endswith('.gz') else open(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, _, feature_type, start, end, _, strand, _, attributes = parts
            start, end = int(start) - 1, int(end)  # GFF is 1-based inclusive
            attr_dict = {k: v.strip('"') for k, v in [x.strip().split('=') if '=' in x else ('ID', x) for x in attributes.split(';') if x]}
            gene = attr_dict.get('gene_name') or attr_dict.get('Name') or attr_dict.get('gene') or attr_dict.get('ID')
            features[feature_type][chrom].addi(start, end, {'strand': strand, 'gene': gene})
    return features

def load_methylation_bed(bed_paths):
    methylation_trees = defaultdict(IntervalTree)
    for bed_path in bed_paths:
        with open(bed_path) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                level = parts[3] if len(parts) > 3 else 'Yes'
                methylation_trees[chrom].addi(start, end, level)
    return methylation_trees

def get_region(pos, chrom, features):
    for region_type in ['CDS', 'five_prime_UTR', 'three_prime_UTR', 'UTR', 'intron']:
        for interval in features.get(region_type, {}).get(chrom, []):
            if interval.begin <= pos < interval.end:
                return region_type, interval.data
    if features.get('gene', {}).get(chrom):
        for interval in features['gene'][chrom]:
            if interval.begin <= pos < interval.end:
                return 'Intron', interval.data
    return 'Intergenic', {'gene': None, 'strand': '+'}

def is_synonymous(ref, alt, pos, chrom, strand, fasta):
    codon_start = (pos // 3) * 3
    codon_seq = fasta[chrom][codon_start:codon_start + 3].seq
    if strand == '-':
        codon_seq = str(Seq(codon_seq).reverse_complement())
        alt_seq = list(codon_seq)
        alt_seq[pos % 3] = str(Seq(alt).complement())
    else:
        alt_seq = list(codon_seq)
        alt_seq[pos % 3] = alt
    ref_aa = str(Seq(codon_seq).translate())
    alt_aa = str(Seq(''.join(alt_seq)).translate())
    return ref_aa == alt_aa

def main():
    parser = argparse.ArgumentParser(description="Annotate SNPs using GFF and FASTA, with optional methylation BED support.")
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--gff', required=True, help='GFF3 annotation file')
    parser.add_argument('--fasta', required=True, help='Reference genome FASTA')
    parser.add_argument('--out', default='annotated_variants.tsv', help='Output file')
    parser.add_argument('--log', default='annotation_errors.log', help='Log file for skipped variants')
    parser.add_argument('--methylation', nargs='*', help='Optional methylation BED file(s)')
    args = parser.parse_args()

    fasta = Fasta(args.fasta)
    features = parse_gff(args.gff)
    methylation_data = load_methylation_bed(args.methylation) if args.methylation else {}

    vcf = VCF(args.vcf)
    with open(args.out, 'w') as out_f, open(args.log, 'w') as log_f:
        out_f.write('\t'.join(['Chrom', 'Pos', 'Ref', 'Alt', 'RegionType', 'Synonymous', 'Gene', 'Methylation']) + '\n')
        for variant in vcf:
            try:
                if len(variant.REF) != 1 or len(variant.ALT[0]) != 1:
                    continue  # Skip non-SNPs

                chrom = variant.CHROM
                pos = variant.POS - 1
                ref = variant.REF
                alt = variant.ALT[0]

                region, data = get_region(pos, chrom, features)
                gene = data.get('gene', '') if data else ''
                strand = data.get('strand', '+') if data else '+'

                synonymous = ''
                if region == 'CDS':
                    try:
                        synonymous = 'Yes' if is_synonymous(ref, alt, pos, chrom, strand, fasta) else 'No'
                    except Exception:
                        synonymous = ''

                # Methylation
                meth_status = 'No'
                if chrom in methylation_data:
                    hits = methylation_data[chrom][pos:pos+1]
                    if hits:
                        meth_status = ';'.join(sorted({str(h.data) for h in hits}))

                out_f.write(f'{chrom}\t{variant.POS}\t{ref}\t{alt}\t{region}\t{synonymous}\t{gene}\t{meth_status}\n')

            except Exception as e:
                log_f.write(f'{variant.CHROM}:{variant.POS} - {str(e)}\n')

if __name__ == '__main__':
    main()
