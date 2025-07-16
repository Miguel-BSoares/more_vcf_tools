# 🧬 Variant Annotation Script with Methylation Support

This script annotates SNPs from a VCF file using a GFF3 genome annotation and reference FASTA. It supports optional annotation from methylation BED files.

---

## ✅ Features

- Identifies whether SNPs fall in:
  - Coding regions (CDS)
  - UTRs (5' or 3')
  - Introns
  - Intergenic regions
- Classifies coding SNPs as:
  - **Synonymous** (does not change amino acid)
  - **Non-synonymous**
- Optional: Annotates whether variants overlap **known methylated regions**

---

## 🔧 Arguments

| Flag           | Description |
|----------------|-------------|
| `--vcf`        | Input VCF file |
| `--gff`        | GFF3 genome annotation file |
| `--fasta`      | Reference genome FASTA file |
| `--out`        | Output file (default: `annotated_variants.tsv`) |
| `--log`        | Error log file (default: `annotation_errors.log`) |
| `--methylation`| One or more optional BED files containing methylated regions |

---

## 📂 Methylation File Format

- Standard BED format:  

- Columns:
- `chrom`, `start`, `end`, `[optional methylation level or label]`
- If no value is given, "Yes" is used by default.

---

## ▶️ Example Usage

```bash
python annotate_variants.py \
--vcf mydata.vcf \
--gff annotation.gff3.gz \
--fasta genome.fa \
--methylation tissue1.bed tissue2.bed \
--out output.tsv \
--log errors.log
