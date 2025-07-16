# Variant Extraction Pipeline

## Overview

This pipeline extracts variants from a large VCF file within user-defined coordinate regions plus optional flanking windows. It outputs:

- Extracted variants in VCF format
- Allele frequency distributions and genotype coding per sample
- Report of coordinate regions with no variants found

It uses [cyvcf2](https://brentp.github.io/cyvcf2/) for fast VCF parsing and supports parallel execution to handle large datasets efficiently.

---

## Installation & Dependencies

- Python 3.6+
- [cyvcf2](https://pypi.org/project/cyvcf2/)

```bash
pip install cyvcf2

Assume these files exist in your working directory:

extract_variants.py 

variants.vcf.gz

variants.vcf.gz.tbi

snps.txt 



## extract variants 

chmod +x extract_variants_cyvcf2.py
python extract_variants_cyvcf2.py your_chunk.vcf.gz snps.txt variants_output -10000 10000

## paralellize execution in large vcf files

chmod +x run_full_pipeline.sh
./run_full_pipeline.sh -v your_large.vcf.gz -c snps.txt -o variants -l -10000 -r 10000 -t 8 -s 50000
