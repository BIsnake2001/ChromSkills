---
name: differential-accessibility
description: Perform differential chromatin accessibility analysis using DESeq2 from BAM or count matrix, including peak merging, count generation, statistical testing, and visualization.
---

# ATAC-seq Differential Accessibility Analysis with DESeq2

## Overview

This skill performs differential accessibility analysis between conditions using **DESeq2**.  
Main steps include:

1. Merge peaks across replicates or samples to build a consensus peak set.  
2. Generate read count matrix over peaks using featureCounts or bedtools.  
3. Prepare sample metadata file describing conditions and replicates.  
4. Perform differential analysis using DESeq2.  
5. Visualize and interpret results (PCA, volcano plot).  

---

## Workflow Decision Tree

### 1. Input Type

- **If starting from BAM files** → Generate consensus peaks and count matrix.  
- **If starting from existing count matrix** → Go directly to DESeq2 analysis.  
- **If multiple conditions or batches** → Include batch/condition in design formula.  

---

## Step-by-Step Workflow

### Step 1: Generate Consensus Peaks

Combine peaks from replicates to define a shared feature space.

```bash
# Merge replicate peaks
bedtools merge -i <rep1_peaks.bed> <rep2_peaks.bed> > consensus_peaks.bed
```

Output: `consensus_peaks.bed`

---

### Step 2: Generate Count Matrix

Count reads overlapping each consensus peak.

```bash
# Option 1: featureCounts (recommended)
featureCounts -T 8 -p -a consensus_peaks.bed     -o atac_counts.txt sample1.bam sample2.bam sample3.bam

# Option 2: bedtools multicov
bedtools multicov -bams sample1.bam sample2.bam sample3.bam     -bed consensus_peaks.bed > atac_counts.txt
```

Output: `atac_counts.txt`

---

### Step 3: Prepare Metadata

Prepare `samples.csv` describing condition and replicate information.

```csv
sample,condition,replicate
sample1,K562,1
sample2,K562,2
sample3,GM12878,1
sample4,GM12878,2
```

---

### Step 4: Differential Accessibility with DESeq2

Run DESeq2 in R for normalization, dispersion estimation, and statistical testing.

```r
library(DESeq2)

counts <- read.table("atac_counts.txt", header=TRUE, row.names=1)
colData <- read.csv("samples.csv", row.names=1)
counts_matrix <- counts[, -c(1:6)]   # remove peak annotation columns if any

dds <- DESeqDataSetFromMatrix(countData=counts_matrix,
                              colData=colData,
                              design=~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","K562","GM12878"))
write.csv(as.data.frame(res), "differential_accessibility_results.csv")
```

Output: `differential_accessibility_results.csv`

---

### Step 5: Visualization and QC

#### PCA Plot

```r
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")
```

#### Volcano Plot

```r
library(EnhancedVolcano)
EnhancedVolcano(res, lab=rownames(res),
                x='log2FoldChange', y='pvalue',
                pCutoff=0.05, FCcutoff=1,
                title='K562 vs GM12878')
```

---

## Advanced Usage

- **Batch effects**: `design = ~ batch + condition`
- **Multi-group comparison**: `contrast=c("condition","A","B")`
- **Time series**: `DESeq(dds, test="LRT", reduced=~1)`
- **Filter low counts**: `dds[rowSums(counts(dds)) >= 20, ]`

---

## Output Files

| File | Description |
|------|--------------|
| `consensus_peaks.bed` | Unified peak set |
| `atac_counts.txt` | Count matrix of reads per peak |
| `samples.csv` | Sample metadata |
| `differential_accessibility_results.csv` | DESeq2 results (log2FC, p-values) |
| `PCA.png`, `Volcano.png` | QC and visualization outputs |

---

## Required Tools

- **featureCounts / bedtools** — read counting  
- **R packages** — `DESeq2`, `ggplot2`, `EnhancedVolcano`, `pheatmap`  

---

## Notes & Troubleshooting

| Issue | Solution |
|-------|-----------|
| Very low counts | Increase threshold (`rowSums >= 20`) |
| Batch effect | Add batch term to design |
| Non-converging model | Use `fitType="local"` or `betaPrior=FALSE` |
| Mismatched sample names | Ensure count column names match metadata rows |
