# 🧬 clean_methylation_bed.py — Methylation Data Cleaner & Converter

This script processes raw DNA methylation data files in **BED-like** or **WIG** format, filters by methylation level, and outputs a clean **BED file** ready for integration with genome annotation pipelines (e.g., for SNP annotation).

---

## ✅ Features

- Supports **WIG (fixedStep / variableStep)** and **BED-like** formats  
- Auto-scales methylation values (0–100 → 0–1 if needed)  
- Filters by user-defined methylation level thresholds  
- Optionally skips rows with missing or invalid data  
- Produces BED format: `chr`, `start`, `end`, `methylation_level`

---

## 📂 Input File Requirements

### 1. BED-like tab-delimited (default mode)

chr1 10468 10469 0.92
chr2 20000 20001 87.5 # will be auto-scaled to 0.875


### 2. WIG formats (auto-detected)
#### fixedStep

fixedStep chrom=chr1 start=1000 step=1 span=1
0.6
0.7


#### variableStep
variableStep chrom=chr1 span=1

1000 0.65
1002 0.9


---

## 🚀 Usage

## in bash 

python clean_methylation_bed.py \
  --input path/to/methylation_data.wig \
  --output output/cleaned_methylation.bed \
  --min_level 0.3 \
  --max_level 1.0 \
  --drop_missing


Filter and clean a BED-like file:

python clean_methylation_bed.py \
  --input human_chr1_methylation.bed \
  --output chr1_cleaned.bed \
  --min_level 0.2 \
  --max_level 0.9 \
  --drop_missing

Convert a WIG file to BED with only high methylation:
  
python clean_methylation_bed.py \
  --input brain_methylation_chr5.wig \
  --output brain_chr5_filtered.bed \
  --min_level 0.8

