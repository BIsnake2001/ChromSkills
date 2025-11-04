---
name: TF-differential-binding
description: This skill performs differential transcription factor (TF) binding analysis from ChIP-seq data using the DiffBind package in R. It provides workflows for sample loading, contrast definition, normalization, statistical testing, and identification of differentially bound regions between experimental conditions.
---

# DiffBind TF Differential Binding Analysis

## Overview

This skill enables comprehensive differential TF binding analysis using **DiffBind** in R. DiffBind integrates read counting, normalization, and statistical modeling (via DESeq2 or edgeR) to identify differentially bound peaks between conditions.

To perform DiffBind differential binding analysis:

1. **Prepare input data**: Provide a sample sheet with ChIP-seq peak files and corresponding BAM files for each sample.
2. **Load data into DiffBind**: Construct a `DBA` object from the sample sheet.
3. **Count reads**: Compute read counts over consensus peak regions.
4. **Define contrasts**: Specify experimental conditions (e.g., treatment vs. control or cell_type_A vs. cell_type_B).
5. **Perform differential analysis**: Run statistical tests to identify differentially bound regions.
6. **Visualize and interpret**: Generate correlation heatmaps, PCA plots, and volcano plots; extract significant binding events.

---

## Quick Start

### 1. Prepare Input Data

Create a CSV sample sheet (`samplesheet.csv`) with the following columns:

| SampleID | Condition | Replicate | bamReads | Peaks | PeakCaller |
|-----------|------------|------------|-----------|--------|-------------|
| TF_A_1    | Control    | 1          | Control1.bam | Control1_peaks.narrowPeak | macs |
| TF_A_2    | Control    | 2          | Control2.bam | Control2_peaks.narrowPeak | macs |
| TF_B_1    | Treated    | 1          | Treated1.bam | Treated1_peaks.narrowPeak | macs |
| TF_B_2    | Treated    | 2          | Treated2.bam | Treated2_peaks.narrowPeak | macs |

### 2. Load Data and Build the DiffBind Object

```r
library(DiffBind)

# Load the sample sheet
samples <- read.csv("samplesheet.csv")

# Create the DBA object
dbObj <- dba(sampleSheet=samples)
```

**Key parameters:**
- `sampleSheet`: CSV file with BAM and peak information
- `PeakCaller`: must match the software used (e.g., "macs", "macs2", "bed")
- Supports both narrowPeak and broadPeak formats

---

## Decision Tree

### Step 1: Read Counting and Consensus Peak Generation

Count reads overlapping consensus peaks across samples:

```r
# Generate a consensus peakset
dbObj <- dba.count(dbObj, summits=250)
```

**Notes:**
- `summits`: re-centers peaks ±250 bp around summits for consistency.
- The resulting matrix contains normalized counts for all samples.

---

### Step 2: Contrast Definition

Define conditions for comparison:

```r
# Define experimental contrasts (e.g., Treated vs Control)
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION, minMembers=2)
```

**Alternatives:**
- For multifactor experiments: use `DBA_TISSUE`, `DBA_TREATMENT`, or custom metadata.
- Check contrasts:
  ```r
  dba.show(dbObj, bContrasts=TRUE)
  ```

---

### Step 3: Differential Binding Analysis

Perform differential analysis using DESeq2 or edgeR:

```r
# Perform analysis
dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2)

# Extract results
diffResults <- dba.report(dbObj, th=0.05, fold=1)
```

**Parameters:**
- `method`: choose `DBA_DESEQ2` (default) or `DBA_EDGER`
- `th`: FDR threshold (default 0.05)
- `fold`: minimum log2 fold change
- `bUsePval=TRUE`: use p-values instead of FDR cutoff

---

### Step 4: Visualization and Quality Control

#### Correlation Heatmap

```r
dba.plotHeatmap(dbObj, correlations=TRUE, scale="row")
```

#### PCA Plot

```r
dba.plotPCA(dbObj, attributes=DBA_CONDITION, label=DBA_ID)
```

#### MA Plot / Volcano Plot

```r
# MA plot
plot(diffResults$Fold, -log10(diffResults$FDR),
     pch=16, col=ifelse(diffResults$FDR < 0.05, "red", "grey"),
     main="Differential Binding MA Plot")

# Volcano plot
with(diffResults, plot(Fold, -log10(FDR),
     col=ifelse(FDR < 0.05 & abs(Fold) > 1, "red", "grey"),
     pch=16, main="Volcano Plot"))
```

---

### Step 5: Result Extraction

Export significant differential peaks:

```r
# Get all significant binding sites
sigSites <- dba.report(dbObj, method=DBA_DESEQ2, th=0.05)

# Export as BED
library(rtracklayer)
export(sigSites, "differential_binding_peaks.bed")
```

---

## Interpretation and Biological Insights

### Significance Criteria

- **FDR < 0.05** → statistically significant  
- **|log2FC| > 1** → biologically meaningful difference  
- **Consistent replicates** → at least two replicates per condition recommended

### Typical Biological Interpretations

- **Increased binding** in treated condition → potential activation or recruitment of TFs
- **Decreased binding** → loss of TF affinity or chromatin closing
- Combine with RNA-seq to correlate with target gene expression.

---

## Example End-to-End Workflow

```r
library(DiffBind)

# Load sample sheet
dbObj <- dba(sampleSheet="samplesheet.csv")

# Count reads and generate consensus peaks
dbObj <- dba.count(dbObj, summits=250)

# Define contrast
dbObj <- dba.contrast(dbObj, categories=DBA_CONDITION, minMembers=2)

# Run DESeq2-based analysis
dbObj <- dba.analyze(dbObj, method=DBA_DESEQ2)

# Extract differential peaks
diffSites <- dba.report(dbObj, th=0.05, fold=1)

# Visualization
dba.plotPCA(dbObj)
dba.plotHeatmap(dbObj)
```

---

## Troubleshooting

| Problem | Possible Cause | Solution |
|----------|----------------|-----------|
| No differential peaks found | Insufficient replicates or low coverage | Increase sequencing depth or lower FDR threshold |
| Errors in sample sheet | Column names incorrect or missing | Use standard DiffBind column format |
| Inconsistent genome build | Mixed genome assemblies | Ensure all BAM and peak files use the same genome reference |
| Over-normalization | Strong batch effects | Include batch term in design or run `dba.contrast(..., block=...)` |

---

## Advanced Topics

### Multi-factor Designs

Include multiple covariates:

```r
dbObj <- dba.contrast(dbObj, categories=c(DBA_TREATMENT, DBA_TISSUE))
```

### Custom Contrasts

Manually define contrasts:

```r
dbObj <- dba.contrast(dbObj, contrast=c("Treated", "Control"), name1="Treated", name2="Control")
```

### Alternative Statistical Engines

```r
# Use edgeR instead of DESeq2
dbObj <- dba.analyze(dbObj, method=DBA_EDGER)
```

### Export Results for Downstream Analysis

```r
write.table(as.data.frame(diffSites), "DiffBind_results.tsv", sep="\t", quote=FALSE)
```
