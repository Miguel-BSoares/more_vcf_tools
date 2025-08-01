import pandas as pd
import numpy as np

# Function to compute allele distribution

def compute_allele_distribution(vcf_file):
    # Read the VCF file
    df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
    allele_distribution = []

    # Iterate through each variant in the VCF
    for index, row in df.iterrows():
        chrom = row[0]
        pos = row[1]
        id_ = row[2]
        ref = row[3]
        alt = row[4]
        genotypes = row[9:]  # Assuming genotypes are in the 10th column onward

        # Initialize counts
        count_0_0 = 0
        count_0_1 = 0
        count_1_1 = 0
        count_dot_dot = 0

        # Count genotypes
        for genotype in genotypes:
            if genotype == '0/0':
                count_0_0 += 1
            elif genotype in ['0/1', '1/0']:
                count_0_1 += 1
            elif genotype == '1/1':
                count_1_1 += 1
            elif genotype == './.':
                count_dot_dot += 1

        # Append results
        allele_distribution.append([
            chrom, pos, id_, ref, alt,
            count_0_0, count_0_1, count_1_1, count_dot_dot
        ])

    # Create a DataFrame
    output_df = pd.DataFrame(allele_distribution, columns=[
        'CHROM', 'POS', 'ID', 'REF', 'ALT',
        'GT_0/0', 'GT_0/1', 'GT_1/1', 'GT_./.'
    ])

    return output_df
