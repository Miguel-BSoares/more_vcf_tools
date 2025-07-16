# 🧬 Variant Extraction Pipeline

## 📂 Overview

This pipeline extracts variants from a large VCF file within user-defined coordinate regions plus optional flanking windows. It outputs:

- Extracted variants in VCF format
- Allele frequency distributions and genotype coding per sample
- Report of coordinate regions with no variants found

It uses [cyvcf2](https://brentp.github.io/cyvcf2/) for fast VCF parsing and supports parallel execution to handle large datasets efficiently.

---

# 📂 Installation & Dependencies

- Python 3.6+
- [cyvcf2](https://pypi.org/project/cyvcf2/)

pip install cyvcf2

# ✅ Required files:

extract_variants.py 

variants.vcf.gz

variants.vcf.gz.tbi

snps.txt 

## ▶️ Example Usage

## 🔧 Arguments

| Flag    | Description |
|---------|-------------|
| `-c`    | Input file with coordinates of SNP of interest |
| `-o`    | output file |
| `-l`    | positions to the left of target SNP (only for `parallel_execution.sh`) |
| `-r`    | positions to the right of target SNP (only for `parallel_execution.sh`) |
| `-t`    | parallelization (only for `parallel_execution.sh`)|
| `-s`    | variants to process by vcf file (only for `parallel_execution.sh`) |

```bash
## extract variants 

chmod +x extract_variants.py
python extract_variants.py file.vcf.gz snps.txt variants_output -10000 10000

## paralellize execution in large vcf files

chmod +x parallel_execution.sh
./parallel_execution.sh -v your_large.vcf.gz -c snps.txt -o variants -l -10000 -r 10000 -t 8 -s 50000
