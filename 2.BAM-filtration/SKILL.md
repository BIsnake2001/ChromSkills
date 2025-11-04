---
name: BAM-filtration
description: Perform the essential alignment-level filtering steps for ChIP-seq/ATAC-seq BAMs (mitochondrial, blacklist, duplicates via Picard, non-primary, unmapped, multi-mapped). The skill auto-detects valid BAMs (sorted + RG), infers genome build from BAM headers, validates required assets (blacklist), and outputs clean, indexed BAMs ready for peak calling and downstream analysis.
---

# BAM Filtration for ChIP-seq / ATAC-seq

## Overview
- Discover input BAMs in the current directory (or those matching a target token), and **only select BAMs that are already coordinate-sorted and contain read group (RG) information**.
- Auto-detect the **genome build** (e.g., hg38, hg19, mm10) from BAM headers; then verify required blacklist files exist in **Resources** provided in this skill. If missing/ambiguous, **prompt the user** to provide them.
- **Make sure the consistence** of the chromosome name between BAM file and blacklist.
- **Remove reads mapping to mitochondrial DNA** (robust to `chrM`/`MT` naming via BED).
- **Remove duplicates using Picard MarkDuplicates** and save duplication metrics.
- **Remove unmapped, non-primary/supplementary, and multi-mapped reads** (via SAM flags + MAPQ threshold).
- **Remove reads overlapping ENCODE blacklist regions** for the inferred genome.
- Re-sort (if needed), **index**, and emit final `*.filtered.bam` (plus `.bai`) as the deliverable for peak calling and downstream analyses.
- **Organizes outputs into**: 
  - **temp/** → intermediate temporary files 
  - **current directory** → essential files for downstream steps (final filtered BAMs and their indexes)
---

## Decision Tree

### 1) Select inputs & target scoping
- **Input discovery**:  
  - If no target is given, select all `*.bam` in the current directory.  
  - If a **target token** is provided (e.g., a TF name like `CTCF`, or a cell line like `K562`), select BAMs whose **filename contains** the token (case-insensitive).  
- **Idempotency**: If `sample.filtered.bam` exists and is newer than its source, it may be skipped (configurable).

### 2) Validate sorting & read group (RG) info
- **Requirement**: Each candidate BAM must be **coordinate-sorted** and have at least one **@RG** entry.  
- **If NOT sorted** → sort + index:
  ```bash
  samtools sort -o sample.sorted.bam sample.bam
  samtools index sample.sorted.bam
  ```
- **If missing RG** → add a minimal RG (recommend users set real metadata):
  ```bash
  picard AddOrReplaceReadGroups     I=sample.sorted.bam O=sample.rg.bam     RGID=grp1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=sample1     VALIDATION_STRINGENCY=SILENT
  samtools index sample.rg.bam
  ```
- **Skill behavior**:  
  - If BAM is already sorted + RG present → keep as-is.  
  - If only one condition is unmet → repair **only that** condition (no redundant work).  
  - Downstream steps reference the most recent valid file (e.g., `sample.rg.bam` if RG was added).

### 3) Auto-detect genome & verify assets
- **Genome inference** from BAM header:
  - Inspect `@SQ` lines: `SN` naming (`chr1` vs `1`, `chrM` vs `MT`), `AS` (assembly), `UR`/`M5` hints, presence of decoy/scaffold patterns that are build-specific (e.g., `KI/GL` for hg38).
  - Infer among common builds (hg38/hg19/mm10/mm9). If ambiguous, prompt user to specify.
- **Assets required** (in current directory or `assets/`):
  - **ENCODE blacklist** for the inferred genome (e.g., `hg38.blacklist.bed`).
- **If assets missing or genome unknown** → **STOP** with a clear prompt:  
  > “Genome could not be determined or required assets are missing. Please provide `<genome>-blacklist.bed(.gz)` and `mito.bed` (or confirm mitochondrial name) for the intended genome (e.g., hg38/hg19/mm10).”

### 4) Remove mitochondrial reads
```bash
samtools idxstats sample.bam | cut -f1 | grep -v -E '^(chrM|MT)$' | xargs samtools view -b sample.bam > sample.noMito.bam
samtools index sample.noMito.bam
```

### 5) Remove duplicates (Picard)
```bash
picard MarkDuplicates   I=sample.noMito.bam O=sample.dedup.bam   M=sample.dedup.metrics.txt REMOVE_DUPLICATES=true   ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT
samtools index sample.dedup.bam
```

### 6) Remove unmapped, non-primary/supplementary, and multi-mapped reads
```bash
samtools view -h -F 260 -F 2048 -q 30 -b sample.dedup.bam > sample.primary.q30.bam
samtools sort -o sample.primary.q30.sorted.bam sample.primary.q30.bam
samtools index sample.primary.q30.sorted.bam
```

### 7) Remove reads overlapping ENCODE blacklist

Before performing the remove step, make sure the chromosome name in blacklist file is consistent to those in the BAM file (with or without "chr").

```bash
bedtools intersect -v -abam sample.primary.q30.sorted.bam -b <genome>.blacklist.bed > sample.filtered.bam
samtools index sample.filtered.bam
```

---

## Outputs
- `temp/sample.dedup.metrics.txt` — Picard duplicate metrics  
- `temp/sample.primary.q30.sorted.bam` & `.bai` — primary, high-MAPQ reads (intermediate)  
- `sample.filtered.bam` & `.bai` — **final deliverable BAM** for downstream analysis 

## Resources

Use BED files in `assets` for blacklist files

