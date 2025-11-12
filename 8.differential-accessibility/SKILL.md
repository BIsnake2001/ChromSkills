---
name: differential-accessibility
description: The differential-accessibility pipeline is designed to identify genomic regions exhibiting significant differences in chromatin accessibility between experimental conditions (e.g., treatment vs. control). It integrates DESeq2 for robust statistical testing and supports workflows starting from either BAM files (aligned ATAC-seq/DNase-seq reads) or a precomputed count matrix. Use it when you aim to identify changes in chromatin accessibility across biological conditions, treatments, or cell types. It is particularly suited for ATAC-seq, DNase-seq, or other open-chromatin profiling datasets.
---

# ATAC-seq Differential Accessibility Analysis with DESeq2

## Overview

This skill performs differential accessibility analysis between conditions.  
Main steps include:

- Refer to the **Inputs & Outputs** section to check inputs and build the output architecture.
- **Always wait the user feedback** if required files are not available in the current working directory by asking "${files} not available, provide required files or skip and proceed ?" 
- Merge peaks across replicates or samples to build a consensus peak set.  
- Generate read count matrix over peaks using featureCounts or bedtools.  
- Prepare sample metadata file describing conditions and replicates.  
- Perform differential analysis using DESeq2.  
- Visualize and interpret results (PCA, volcano plot). 
- Output significantly up and down accessible regions.

---

## When to use this sill
Use the differential-accessibility pipeline when you aim to identify changes in chromatin accessibility across biological conditions, treatments, or cell types. It is particularly suited for ATAC-seq, DNase-seq, or other open-chromatin profiling datasets.

Recommended scenarios include:

- Comparing treated vs. control samples to identify regulatory regions responsive to a drug, signaling molecule, or environmental change.
- Investigating cell differentiation or developmental trajectories to reveal dynamic chromatin remodeling.
- Analyzing disease vs. normal tissues to pinpoint dysregulated enhancer or promoter accessibility.
- Integrating with RNA-seq or ChIP-seq data to connect chromatin accessibility with transcriptional or epigenetic regulation.

The pipeline performs best with datasets containing biological replicates (≥2 per condition) and moderate to high sequencing depth (~20–50 million reads per sample).

---

## Inputs & Outputs

### Inputs (choose one)

- If starting from BAM files and BED peak files → Generate consensus peaks and count matrix.  
- If starting from existing count matrix → Go directly to DESeq2 analysis.  
- If multiple conditions or batches → Include batch/condition in design 

### Outputs

```bash
DAR_analysis/
    tables/
      all_peaks.bed
      consensus_peaks.bed # Unified peak set
      atac_counts.txt # Count matrix of reads per peak
      samples.csv # Sample metadata
    DARs/
      DAR_results.csv # DESeq2 results (log2FC, p-values)
      DAR_sig.bed # Significantly diffential accessible regions
      DAR_up.bed
      DAR_down.bed  
    plots/ # visualization outputs
      PCA.pdf
      Volcano.pdf
    logs/ # analysis logs 
    temp/ # other temp files
```
---

## Decision Tree

### Step 1: Generate Consensus Peaks

Combine peaks from replicates to define a shared feature space.

```bash
# Merge replicate peaks
cat <rep1_peaks.bed>  <rep2_peaks.bed> | sort -k1,1 -k2,,2n input.bed > all_peaks.bed
bedtools merge -i <rep1_peaks.bed> <rep2_peaks.bed> > consensus_peaks.bed
```

Output: `consensus_peaks.bed`

---

### Step 2: Generate Count Matrix

Count reads overlapping each consensus peak.

```bash
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"} {print "peak_"NR, $1, $2+1, $3, "."}' consensus_peaks.bed > consensus_peaks.saf

# the -p -B -C parameter for featureCounts is for pair-end reads
featureCounts -a consensus_peaks.saf -F SAF -o atac_counts.txt -T 8  -p -B -C sample1.bam sample2.bam sample3.bam
```

Output: `atac_counts.txt`

---

### Step 3: Prepare Metadata

Prepare `samples.csv` describing condition and replicate information.

```csv
sample,condition,replicate
sample1,c1,1
sample2,c1,2
sample3,c2,1
sample4,c2,2
```

---

### Step 4: Differential Accessibility with DESeq2

Run DESeq2 in R for normalization, dispersion estimation, and statistical testing.

```r
library(DESeq2)

counts <- read.table("atac_counts.txt", header=TRUE, row.names=1)
colData <- read.csv("samples.csv", row.names=1)
counts_matrix <- counts[, -c(1:6)]   # Keep only the count columns

dds <- DESeqDataSetFromMatrix(countData=counts_matrix,
                              colData=colData,
                              design=~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","c1","c2"))

region_info <- region_info[rownames(dds), ]
res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)
res_with_coords <- cbind(region_info, res_df)

write.csv(as.data.frame(res), "DAR_results.csv") # output DAR results with region loci

```

Output: `DAR_results.csv`

---

### Step 5: Visualization and QC

1) PCA Plot

```r
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")
```

2) Volcano Plot

```r
library(EnhancedVolcano)
EnhancedVolcano(res, lab=rownames(res),
                x='log2FoldChange', y='pvalue',
                pCutoff=0.05, FCcutoff=1,
                title='"c1 vs c2')
```

---

### Step 6: Output significantly up and down accessible regions

```r
sig <- subset(res_with_coords, !is.na(padj) & padj < 0.05)

to_bed4 <- function(df) {
  data.frame(chr = df$Chr,
             start = as.integer(df$Start - 1),  # SAF -> BED
             end = as.integer(df$End),
             name = df$GeneID,
             check.names = FALSE)
}

write.table(to_bed4(sig),  "DAR_sig.bed",  sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(to_bed4(subset(sig, log2FoldChange > 0)), "DAR_up.bed",  sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(to_bed4(subset(sig, log2FoldChange < 0)), "DAR_down.bed",sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
```
Output: `DAR_sig.bed` `DAR_up.bed` `DAR_down.bed`

## Advanced Usage

- **Batch effects**: `design = ~ batch + condition`
- **Multi-group comparison**: `contrast=c("condition","A","B")`
- **Time series**: `DESeq(dds, test="LRT", reduced=~1)`
- **Filter low counts**: `dds[rowSums(counts(dds)) >= 20, ]`

---

## Notes & Troubleshooting

| Issue | Solution |
|-------|-----------|
| Very low counts | Increase threshold (`rowSums >= 20`) |
| Batch effect | Add batch term to design |
| Non-converging model | Use `fitType="local"` or `betaPrior=FALSE` |
| Mismatched sample names | Ensure count column names match metadata rows |
