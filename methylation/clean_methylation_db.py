#!/usr/bin/env python3

import argparse
import re

def process_bed_like(infile, outfile, min_level, max_level, drop_missing):
    for line in infile:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 4:
            if drop_missing:
                continue
            else:
                parts.append("NA")
        try:
            level = float(parts[3])
            if level > 1.0:
                level /= 100.0
            if min_level <= level <= max_level:
                outfile.write(f'{parts[0]}\t{parts[1]}\t{parts[2]}\t{round(level, 3)}\n')
        except ValueError:
            if not drop_missing:
                outfile.write(f'{parts[0]}\t{parts[1]}\t{parts[2]}\tNA\n')

def process_wig(infile, outfile, min_level, max_level):
    chrom = None
    span = 1
    is_fixed = False
    pos = None

    for line in infile:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        if line.startswith('fixedStep'):
            is_fixed = True
            chrom = re.search(r'chrom=(\S+)', line).group(1)
            pos = int(re.search(r'start=(\d+)', line).group(1)) - 1
            span_match = re.search(r'span=(\d+)', line)
            span = int(span_match.group(1)) if span_match else 1
            continue

        if line.startswith('variableStep'):
            is_fixed = False
            chrom = re.search(r'chrom=(\S+)', line).group(1)
            span_match = re.search(r'span=(\d+)', line)
            span = int(span_match.group(1)) if span_match else 1
            continue

        try:
            if is_fixed:
                level = float(line)
                end = pos + span
                if level > 1.0:
                    level /= 100.0
                if min_level <= level <= max_level:
                    outfile.write(f'{chrom}\t{pos}\t{end}\t{round(level, 3)}\n')
                pos += span
            else:
                wig_parts = line.split()
                if len(wig_parts) != 2:
                    continue
                pos = int(wig_parts[0]) - 1
                level = float(wig_parts[1])
                if level > 1.0:
                    level /= 100.0
                if min_level <= level <= max_level:
                    outfile.write(f'{chrom}\t{pos}\t{pos + span}\t{round(level, 3)}\n')
        except ValueError:
            continue

def is_wig(file_handle):
    for _ in range(20):
        line = file_handle.readline()
        if line.startswith("fixedStep") or line.startswith("variableStep"):
            file_handle.seek(0)
            return True
    file_handle.seek(0)
    return False

def main():
    parser = argparse.ArgumentParser(description="Clean and convert methylation WIG or BED-like files into BED format.")
    parser.add_argument('--input', required=True, help='Input WIG or BED-like methylation file')
    parser.add_argument('--output', required=True, help='Output cleaned BED file')
    parser.add_argument('--min_level', type=float, default=0.0, help='Minimum methylation level to include (default: 0.0)')
    parser.add_argument('--max_level', type=float, default=1.0, help='Maximum methylation level to include (default: 1.0)')
    parser.add_argument('--drop_missing', action='store_true', help='Drop rows with missing or malformed methylation values')
    args = parser.parse_args()

    with open(args.input) as infile, open(args.output, 'w') as outfile:
        if is_wig(infile):
            process_wig(infile, outfile, args.min_level, args.max_level)
        else:
            process_bed_like(infile, outfile, args.min_level, args.max_level, args.drop_missing)

if __name__ == '__main__':
    main()
