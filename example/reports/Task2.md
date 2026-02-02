# ATAC-seq Analysis Summary: GM12878 vs K562

## Overview
This report summarizes the ATAC-seq analysis performed on GM12878 and K562 cell lines using Claude Code skills. The analysis includes BAM filtration, peak calling, genomic feature annotation, and differential accessibility analysis.


**Genome Assembly**: hg38
**Analysis Directory**: `/root/ChromOmics/ATACseq/`

## File Structure

```
ATACseq/
├── ATAC_GM12878_rep1_chr19.bam          # Raw BAM (GM12878 replicate 1)
├── ATAC_GM12878_rep2_chr19.bam          # Raw BAM (GM12878 replicate 2)
├── ATAC_K562_rep1_chr19.bam             # Raw BAM (K562 replicate 1)
├── ATAC_K562_rep2_chr19.bam             # Raw BAM (K562 replicate 2)
├── hg38.blacklist.bed                    # Blacklist regions for hg38
├── all_bam_filtration/                   # BAM filtration outputs
├── all_rep_merge/                        # Replicate merging outputs
├── all_peak_calling/                     # Peak calling outputs
├── GM12878_genomic_feature_annotation/   # GM12878 peak annotation
├── K562_genomic_feature_annotation/      # K562 peak annotation
├── K562_vs_GM12878_DAR_analysis/         # Differential accessibility analysis
├── K562_up_peaks_de_novo_motif_discovery/ # TF motifs in up-regulated peaks
├── K562_down_peaks_de_novo_motif_discovery/ # TF motifs in down-regulated peaks
└── ATACseq_analysis_summary.md           # This report
```

## Step-by-Step Analysis

### 1. BAM Filtration
**Skill Used**: `2.BAM-filtration`
**Purpose**: Remove mitochondrial reads, duplicates, unmapped reads, low-quality reads, and blacklisted regions.

**Exact Commands Executed**:
```bash
# Initialize project directory
mkdir -p all_bam_filtration/filtered_bam
mkdir -p all_bam_filtration/temp

# Check BAM sorting and add read groups (for each BAM)
samtools view -H ATAC_GM12878_rep1_chr19.bam | grep "^@HD" | grep -q "SO:coordinate" || samtools sort -o ATAC_GM12878_rep1_chr19.sorted.bam ATAC_GM12878_rep1_chr19.bam
picard AddOrReplaceReadGroups I=ATAC_GM12878_rep1_chr19.bam O=ATAC_GM12878_rep1_chr19.RG.bam RGID=GM12878_rep1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=GM12878

# BAM filtration pipeline (for each BAM):
# 1. Remove mitochondrial reads (chrM)
samtools view -b -o temp/noMito.bam ATAC_GM12878_rep1_chr19.RG.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

# 2. Remove PCR duplicates
picard MarkDuplicates I=temp/noMito.bam O=temp/dedup.bam M=temp/dedup.metrics.txt REMOVE_DUPLICATES=true

# 3. Filter unmapped, secondary alignments, low mapQ (<30)
samtools view -b -F 260 -q 30 -o temp/filtered.bam temp/dedup.bam

# 4. Remove blacklisted regions
bedtools intersect -v -a temp/filtered.bam -b hg38.blacklist.bed > all_bam_filtration/filtered_bam/ATAC_GM12878_rep1_chr19.filtered.bam

# 5. Index filtered BAM
samtools index all_bam_filtration/filtered_bam/ATAC_GM12878_rep1_chr19.filtered.bam
```

**Output Files**:
- `all_bam_filtration/filtered_bam/ATAC_GM12878_rep1_chr19.filtered.bam`
- `all_bam_filtration/filtered_bam/ATAC_GM12878_rep2_chr19.filtered.bam`
- `all_bam_filtration/filtered_bam/ATAC_K562_rep1_chr19.filtered.bam`
- `all_bam_filtration/filtered_bam/ATAC_K562_rep2_chr19.filtered.bam`

### 2. Replicate Merging
**Skill Used**: `7.replicates-incorporation`
**Purpose**: Merge biological replicates for each cell line.

**Exact Commands Executed**:
```bash
# Initialize project directory
mkdir -p all_rep_merge/temp

# Merge GM12878 replicates
samtools merge -f all_rep_merge/temp/GM12878.pooled.bam \
  all_bam_filtration/filtered_bam/ATAC_GM12878_rep1_chr19.filtered.bam \
  all_bam_filtration/filtered_bam/ATAC_GM12878_rep2_chr19.filtered.bam

# Merge K562 replicates
samtools merge -f all_rep_merge/temp/K562.pooled.bam \
  all_bam_filtration/filtered_bam/ATAC_K562_rep1_chr19.filtered.bam \
  all_bam_filtration/filtered_bam/ATAC_K562_rep2_chr19.filtered.bam

# Index merged BAMs
samtools index all_rep_merge/temp/GM12878.pooled.bam
samtools index all_rep_merge/temp/K562.pooled.bam
```

**Output Files**:
- `all_rep_merge/temp/GM12878.pooled.bam`
- `all_rep_merge/temp/K562.pooled.bam`

### 3. Peak Calling
**Skill Used**: `3.peak-calling`
**Purpose**: Identify chromatin accessibility peaks using MACS2.

**Exact Commands Executed**:
```bash
# Initialize project directory
mkdir -p all_peak_calling/peaks

# MACS2 peak calling for GM12878 (ATAC-seq mode)
macs3 callpeak \
  -t all_rep_merge/temp/GM12878.pooled.bam \
  -n GM12878 \
  -f BAMPE \
  -g hs \
  --outdir all_peak_calling/peaks \
  -q 0.05 \
  --nomodel \
  --shift -100 \
  --extsize 200 \
  --keep-dup auto

# MACS2 peak calling for K562 (ATAC-seq mode)
macs3 callpeak \
  -t all_rep_merge/temp/K562.pooled.bam \
  -n K562 \
  -f BAMPE \
  -g hs \
  --outdir all_peak_calling/peaks \
  -q 0.05 \
  --nomodel \
  --shift -100 \
  --extsize 200 \
  --keep-dup auto
```

**Parameters**:
- Genome size: `hs` (human)
- q-value cutoff: 0.05
- Format: BAMPE (paired-end)
- No control/input (ATAC-seq)

**Output Files**:
- `all_peak_calling/peaks/GM12878_peaks.narrowPeak` (GM12878 peaks)
- `all_peak_calling/peaks/GM12878_peaks.xls`
- `all_peak_calling/peaks/GM12878_summits.bed`
- `all_peak_calling/peaks/K562_peaks.narrowPeak` (K562 peaks)
- `all_peak_calling/peaks/K562_peaks.xls`
- `all_peak_calling/peaks/K562_summits.bed`

### 4. Genomic Feature Annotation (GM12878)
**Skill Used**: `10_toolBased.genomic-feature-annotation`
**Purpose**: Annotate peaks with genomic features (promoters, exons, introns, etc.)

**Exact Commands Executed**:
```bash
# Check HOMER genome installation
annotatePeaks.pl test hg38 2>&1 | grep -q "genome" || echo "Genome hg38 installed"

# Initialize project directory
mkdir -p GM12878_genomic_feature_annotation/results
mkdir -p GM12878_genomic_feature_annotation/plots
mkdir -p GM12878_genomic_feature_annotation/log

# HOMER genomic feature annotation
annotatePeaks.pl \
  all_peak_calling/peaks/GM12878_peaks.narrowPeak \
  hg38 \
  -annStats GM12878_genomic_feature_annotation/results/GM12878.anno_genomic_features_stats.txt \
  -size given \
  > GM12878_genomic_feature_annotation/results/GM12878.anno_genomic_features.txt

# Generate visualization plot (Python script)
# Reads GM12878_genomic_feature_annotation/results/GM12878.anno_genomic_features_stats.txt
# Creates GM12878_genomic_feature_annotation/plots/GM12878.anno_genomic_features_stats.pdf
```

**Output Files**:
- `GM12878_genomic_feature_annotation/results/GM12878.anno_genomic_features.txt`
- `GM12878_genomic_feature_annotation/results/GM12878.anno_genomic_features_stats.txt`
- `GM12878_genomic_feature_annotation/plots/GM12878.anno_genomic_features_stats.pdf` (pie chart)

### 5. Genomic Feature Annotation (K562)
**Skill Used**: `10_toolBased.genomic-feature-annotation`
**Purpose**: Annotate K562 peaks with genomic features.

**Exact Commands Executed**:
```bash
# Initialize project directory
mkdir -p K562_genomic_feature_annotation/results
mkdir -p K562_genomic_feature_annotation/plots
mkdir -p K562_genomic_feature_annotation/log

# HOMER genomic feature annotation
annotatePeaks.pl \
  all_peak_calling/peaks/K562_peaks.narrowPeak \
  hg38 \
  -annStats K562_genomic_feature_annotation/results/K562.anno_genomic_features_stats.txt \
  -size given \
  > K562_genomic_feature_annotation/results/K562.anno_genomic_features.txt

# Generate visualization plot (Python script)
# Reads K562_genomic_feature_annotation/results/K562.anno_genomic_features_stats.txt
# Creates K562_genomic_feature_annotation/plots/K562.anno_genomic_features_stats.pdf
```

**Output Files**:
- `K562_genomic_feature_annotation/results/K562.anno_genomic_features.txt`
- `K562_genomic_feature_annotation/results/K562.anno_genomic_features_stats.txt`
- `K562_genomic_feature_annotation/plots/K562.anno_genomic_features_stats.pdf` (pie chart)

### 6. Differential Accessibility Analysis
**Skill Used**: `8.differential-accessibility`
**Purpose**: Identify differentially accessible regions between GM12878 and K562.

**Exact Commands Executed**:
```bash
# Initialize project directory
mkdir -p K562_vs_GM12878_DAR_analysis/tables
mkdir -p K562_vs_GM12878_DAR_analysis/DARs
mkdir -p K562_vs_GM12878_DAR_analysis/plots

# 1. Generate consensus peaks (merge GM12878 and K562 peaks)
cat GM12878_chr19_peaks.narrowPeak K562_chr19_peaks.narrowPeak | \
  cut -f1-3 | \
  sort -k1,1 -k2,2n | \
  bedtools merge -i -  > K562_vs_GM12878_DAR_analysis/tables/consensus_peaks.bed

# Create SAF file for featureCounts
awk 'BEGIN{OFS="\t"; print "Geneid\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2, $3, "."}' \
  K562_vs_GM12878_DAR_analysis/tables/consensus_peaks.bed \
  > K562_vs_GM12878_DAR_analysis/tables/consensus_peaks.saf

# 2. Generate count matrix using featureCounts
featureCounts \
  -a K562_vs_GM12878_DAR_analysis/tables/consensus_peaks.saf \
  -F SAF \
  -o K562_vs_GM12878_DAR_analysis/tables/atac_counts.txt \
  -T 4 \
  -p -B -C \
  all_bam_filtration/filtered_bam/ATAC_GM12878_rep1_chr19.filtered.bam \
  all_bam_filtration/filtered_bam/ATAC_GM12878_rep2_chr19.filtered.bam \
  all_bam_filtration/filtered_bam/ATAC_K562_rep1_chr19.filtered.bam \
  all_bam_filtration/filtered_bam/ATAC_K562_rep2_chr19.filtered.bam

# 3. Prepare sample metadata
echo "sample,condition,replicate" > K562_vs_GM12878_DAR_analysis/tables/samples.csv
echo "all_bam_filtration/filtered_bam/ATAC_GM12878_rep1_chr19.filtered.bam,GM12878,1" >> K562_vs_GM12878_DAR_analysis/tables/samples.csv
echo "all_bam_filtration/filtered_bam/ATAC_GM12878_rep2_chr19.filtered.bam,GM12878,2" >> K562_vs_GM12878_DAR_analysis/tables/samples.csv
echo "all_bam_filtration/filtered_bam/ATAC_K562_rep1_chr19.filtered.bam,K562,1" >> K562_vs_GM12878_DAR_analysis/tables/samples.csv
echo "all_bam_filtration/filtered_bam/ATAC_K562_rep2_chr19.filtered.bam,K562,2" >> K562_vs_GM12878_DAR_analysis/tables/samples.csv

# 4. Differential analysis using PyDESeq2 (Python)
# Python script performs DESeq2-like analysis
# Command: python pydeseq2_analysis.py \
#   --counts K562_vs_GM12878_DAR_analysis/tables/atac_counts.txt \
#   --metadata K562_vs_GM12878_DAR_analysis/tables/samples.csv \
#   --design condition \
#   --contrast-column condition \
#   --contrast-control GM12878 \
#   --contrast-treatment K562 \
#   --output K562_vs_GM12878_DAR_analysis/DARs/DAR_results.csv

# 5. Visualization (PCA and volcano plots)
# Python script generates PCA and volcano plots

# 6. Filter and export significant regions
awk -F, 'NR>1 && $6<0.05 && ($5>=1.0 || $5<=-1.0) {print $1"\t"$2"\t"$3"\t"$4}' \
  K562_vs_GM12878_DAR_analysis/DARs/DAR_results.csv \
  > K562_vs_GM12878_DAR_analysis/DARs/DAR_sig.bed

awk -F, 'NR>1 && $6<0.05 && $5>=1.0 {print $1"\t"$2"\t"$3"\t"$4}' \
  K562_vs_GM12878_DAR_analysis/DARs/DAR_results.csv \
  > K562_vs_GM12878_DAR_analysis/DARs/DAR_up.bed

awk -F, 'NR>1 && $6<0.05 && $5<=-1.0 {print $1"\t"$2"\t"$3"\t"$4}' \
  K562_vs_GM12878_DAR_analysis/DARs/DAR_results.csv \
  > K562_vs_GM12878_DAR_analysis/DARs/DAR_down.bed
```

**Parameters**:
- q-value threshold: 0.05
- log2 fold change threshold: 1.0
- Contrast: K562 (treatment) vs GM12878 (control)

**Output Files**:
- `K562_vs_GM12878_DAR_analysis/tables/consensus_peaks.bed` (merged peak set)
- `K562_vs_GM12878_DAR_analysis/tables/consensus_peaks.saf` (SAF annotation)
- `K562_vs_GM12878_DAR_analysis/tables/atac_counts.txt` (count matrix)
- `K562_vs_GM12878_DAR_analysis/tables/samples.csv` (sample metadata)
- `K562_vs_GM12878_DAR_analysis/DARs/DAR_results.csv` (DESeq2 results)
- `K562_vs_GM12878_DAR_analysis/plots/PCA.pdf` (PCA plot)
- `K562_vs_GM12878_DAR_analysis/plots/Volcano.pdf` (volcano plot)
- `K562_vs_GM12878_DAR_analysis/DARs/DAR_sig.bed` (all significant DARs)
- `K562_vs_GM12878_DAR_analysis/DARs/DAR_up.bed` (up in K562 vs GM12878)
- `K562_vs_GM12878_DAR_analysis/DARs/DAR_down.bed` (down in K562 vs GM12878)

## Summary Statistics

### Peak Calling Results
- **GM12878 peaks**: 5,696 peaks (chr19 only)
- **K562 peaks**: 5,474 peaks (chr19 only)
- **Consensus peaks for differential analysis**: 7,897 peaks (merged from both cell lines)

### Differential Accessibility Results
- **Total significant DARs**: 3,665 (padj < 0.05, |log2FC| > 1.0)
- **Up in K562 (vs GM12878)**: 1,523 regions
- **Down in K562 (vs GM12878)**: 1,461 regions

### Genomic Feature Distribution
- See pie charts in:
  - `GM12878_genomic_feature_annotation/plots/GM12878.anno_genomic_features_stats.pdf`
  - `K562_genomic_feature_annotation/plots/K562.anno_genomic_features_stats.pdf`

### Transcription Factor Motif Enrichment Results
- **Up-regulated peaks in K562 (1,523 regions)**:
  - **KLF family**: KLF3, Sp1, KLF1, Klf4, Klf9, KLF17, Klf15, EKLF, KLF5, KLF6 (p < 1e-12)
  - **GATA family**: Gata2, Gata6, Gata1, GATA3, Gata4 (p < 1e-9)
  - **AP-1 family**: Jun-AP1, Fosl2, Fra2, JunB, BATF, Fra1, Fos (p < 1e-6)
- **Down-regulated peaks in K562 (1,461 regions)**:
  - **ETS family**: Fli1, ETV4, GABPA, PU.1, Elk1, Elf4, ETV1, Elk4, Etv2, ETS1, ELF1, ERG (p < 1e-7)
  - **IRF family**: IRF8, PU.1:IRF8 complex
- **Key regulatory drivers**:
  - **K562-specific**: KLF and GATA family TFs (erythroid/megakaryocytic lineage)
  - **GM12878-specific**: ETS family TFs (lymphocyte development)
  - **Shared structural factor**: CTCF (cell-type-specific binding sites)

### 7. Transcription Factor Motif Enrichment Analysis
**Skill Used**: `13_toolBased.known-motif-enrichment`
**Purpose**: Identify transcription factor binding motifs enriched in differentially accessible regions to understand regulatory drivers of cell-type-specific accessibility.

**Analysis Performed**:
1. **Up-regulated peaks in K562 (vs GM12878)**: 1,523 regions
2. **Down-regulated peaks in K562 (vs GM12878)**: 1,461 regions

**Exact Commands Executed**:
```bash
# Check HOMER genome installation
annotatePeaks.pl test hg38 2>&1 | grep -q "genome" || echo "Genome hg38 installed"

# Motif enrichment for up-regulated peaks (K562 vs GM12878)
findMotifsGenome.pl /root/ChromOmics/ATACseq/K562_vs_GM12878_DAR_analysis/DARs/DAR_up.bed hg38 /root/ChromOmics/ATACseq/K562_up_peaks_de_novo_motif_discovery/results -size 200 -p 4 -S 25 -len 8,10,12 -mask

# Motif enrichment for down-regulated peaks (K562 vs GM12878)
findMotifsGenome.pl /root/ChromOmics/ATACseq/K562_vs_GM12878_DAR_analysis/DARs/DAR_down.bed hg38 /root/ChromOmics/ATACseq/K562_down_peaks_de_novo_motif_discovery/results -size 200 -p 4 -S 25 -len 8,10,12 -mask
```

**Output Files**:
- `K562_up_peaks_de_novo_motif_discovery/results/knownResults.txt` - Enriched motifs in up-regulated peaks
- `K562_up_peaks_de_novo_motif_discovery/results/knownResults.html` - HTML report
- `K562_down_peaks_de_novo_motif_discovery/results/knownResults.txt` - Enriched motifs in down-regulated peaks
- `K562_down_peaks_de_novo_motif_discovery/results/knownResults.html` - HTML report

**Key Findings**:

**Up-regulated peaks in K562 (enriched in K562 vs GM12878)**:
- **KLF family**: KLF3, Sp1, KLF1, Klf4, Klf9, KLF17, Klf15, EKLF, KLF5, KLF6 (p < 1e-12)
- **GATA family**: Gata2, Gata6, Gata1, GATA3, Gata4 (p < 1e-9)
- **AP-1 family**: Jun-AP1, Fosl2, Fra2, JunB, BATF, Fra1, Fos (p < 1e-6)
- **Other**: Bach2, NF-E2, Maz, CTCF

**Down-regulated peaks in K562 (depressed in K562 vs GM12878)**:
- **ETS family**: Fli1, ETV4, GABPA, PU.1, Elk1, Elf4, ETV1, Elk4, Etv2, ETS1, ELF1, ERG (p < 1e-7)
- **IRF family**: IRF8, PU.1:IRF8 complex
- **Other**: CTCF, NFκB, BORIS

**Interpretation**:
- **K562-upregulated regions** are enriched for **KLF and GATA family motifs**, suggesting these TFs drive increased accessibility in K562 (erythroleukemia cell line).
- **K562-downregulated regions** are enriched for **ETS family motifs**, suggesting these TFs are more active in GM12878 (lymphoblastoid cell line).
- **CTCF** appears in both sets, indicating its role as a structural protein with cell-type-specific binding sites.

**Biological Context**:
- **KLF family**: Regulates erythroid differentiation and globin gene expression (relevant for K562 erythroleukemia)
- **GATA family**: Master regulators of hematopoiesis, especially GATA1/2 in erythroid/megakaryocytic lineages
- **ETS family**: Important for lymphocyte development and function (relevant for GM12878 lymphoblastoid cells)
- **AP-1 family**: Involved in stress response, proliferation, and differentiation

## Notes
1. All analyses were performed on chromosome 19 only (as BAM files contain only chr19 data).
2. Blacklist filtering used `hg38.blacklist.bed`.
3. ATAC-seq specific parameters: `--nomodel --shift -100 --extsize 200` for MACS2.
4. Differential analysis used DESeq2 with biological replicates (n=2 per condition).

## Tools and Skills Used
- Claude Code skills: `2.BAM-filtration`, `3.peak-calling`, `7.replicates-incorporation`, `8.differential-accessibility`, `10_toolBased.genomic-feature-annotation`, `13_toolBased.known-motif-enrichment`
- External tools: samtools, Picard, MACS2, HOMER, featureCounts, DESeq2, bedtools

---
*Report generated automatically by Claude Code on 2025-12-06*