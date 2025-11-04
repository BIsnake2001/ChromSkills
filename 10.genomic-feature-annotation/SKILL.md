---
name: genomic-feature-annotation
description: This skill is used to perform genomic feature annotation and visualization for any file containing genomic region information using Homer (Hypergeometric Optimization of Motif EnRichment). It annotates regions such as promoters, exons, introns, intergenic regions, and TSS proximity, and generates visual summaries of feature distributions.
---

# Genomic Feature Annotation and Visualization with Homer

## Overview

1. Prepare genomic region files in BED or other supported formats. Ensure that the input genomic regions are provided in a valid BED format (chromosome, start, end). If the file does not meet this format, extract the required columns to create a valid BED file.regions file.
2. Identify and specify the correct genome assembly for annotation.  
3. Annotate the genomic regions using Homer's `annotatePeaks.pl`.  
4. Generate annotation statistics and feature distribution summaries.  
5. Visualize annotation results (e.g., pie charts, barplots).  

---

## Decision Tree

### 1. Input File Preparation

Accepts any genomic region file containing chromosome, start, and end coordinates:  
- **Formats supported:** BED, narrowPeak, broadPeak, or Homer-style peak files.  
- **Command example:**
  ```bash
  head -n 5 regions.bed
  ```
  Ensure files have at least three columns: `chr`, `start`, `end`. 
  If more than three columns, make sure the 6th column is the strand information.

---

### 2. Genome Assembly Selection

Select a genome that matches the coordinates of your regions:  
- Human: `hg38`, `hg19`  
- Mouse: `mm10`, `mm9`  
- Others can be installed via `configureHomer.pl`.  

**Command example:**  
```bash
# List available genomes
find $HOMER_DATA/genomes -maxdepth 1 -type d
```

If genome not installed:  
```bash
configureHomer.pl -install hg38
```

---

### 3. Genomic Feature Annotation

Use Homer's `annotatePeaks.pl` to assign each region to genomic features.  

**Command example:**  
```bash
annotatePeaks.pl regions.bed hg38 -annStats annotation_stats.txt > annotated_regions.txt
```

**Key parameters:**  
- `-annStats`: Generate annotation statistics file.  
- `-size given`: Keep original region sizes.  
- `-hist 10`: Create histogram of distances to TSS (10 bp bins).  
- `-CpG`: Include CpG information.  

---

### 4. Visualization of Annotation Results

After annotation, visualize the feature composition using R or Python.  
A typical visualization could be a **pie chart** or **barplot** showing feature proportions.

**Example command (R):**
```R
data <- read.table("annotation_stats.txt", header=TRUE, sep="\t")
pie(data$Number.of.Peaks, labels=data$Annotation, main="Genomic Feature Distribution")
```

**Example command (Python):**
```python
import pandas as pd
import matplotlib.pyplot as plt
data = pd.read_csv("annotation_stats.txt", sep="\t")
plt.pie(data['# of peaks'], labels=data['Annotation'], autopct='%1.1f%%')
plt.title('Genomic Feature Annotation Distribution')
plt.show()
```

---

### 5. Interpretation of Results

Typical annotation categories:
- **Promoter**: -1 kb to +100 bp from TSS  
- **5' UTR**, **Exon**, **Intron**, **3' UTR**, **Intergenic**, **TTS**  

Quality indicators:
- **Annotation rate**: % of peaks successfully annotated.  
- **Promoter fraction**: Often high in TF ChIP-seq.  
- **Intergenic fraction**: Reflects enhancer-rich or noncoding regions.  

---

### 6. Example Use Cases

**(1) Transcription Factor Binding Sites (TF ChIP-seq):**  
```bash
annotatePeaks.pl TF_binding.bed hg38 -size 1000 -hist 50 > TF_annotation.txt
```

**(2) Histone Modification ChIP-seq:**  
```bash
annotatePeaks.pl H3K27ac_peaks.bed hg38 -annStats H3K27ac_stats.txt > H3K27ac_annotation.txt
```

**(3) ATAC-seq Accessible Regions:**  
```bash
annotatePeaks.pl ATAC_peaks.bed hg38 -size 5000 -hist 100 > ATAC_annotation.txt
```

---

### 7. Troubleshooting

| Issue | Description | Solution |
|--------|--------------|-----------|
| Genome not found | Homer cannot find genome reference | Run `configureHomer.pl -install <genome>` |
| Format error | BED file malformed | Use `awk '{if(NF<3) print NR}' regions.bed` |
| Empty output | Wrong genome assembly | Check chromosome naming conventions (e.g., chr1 vs 1) |

---

## Output Files

| File | Description |
|------|--------------|
| `annotated_regions.txt` | Annotated genomic regions |
| `annotation_stats.txt` | Summary of peak annotation counts |
| `annotation_plot.png` (optional) | Visualization of feature distribution |

---

## Best Practices

- Use high-confidence regions (e.g., IDR-filtered peaks).  
- Ensure genome naming convention matches input files.  
- Use visualization to assess annotation patterns across datasets.  
- Save annotation parameters and plots for reproducibility.  

---

## Example Workflow Summary

```bash
# Step 1: Annotate genomic regions
annotatePeaks.pl regions.bed hg38 -annStats annotation_stats.txt > annotated_regions.txt

# Step 2: Visualize annotation
Rscript plot_annotation.R annotation_stats.txt annotation_plot.png
```
