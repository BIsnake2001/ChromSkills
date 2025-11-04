# Methodology Details

## Overview

This document provides detailed methodology for the cell type transcription factor driver analysis pipeline.

## Data Processing Steps

### 1. BED File Processing

**Input Validation:**
- Verify BED file format (chromosome, start, end columns)
- Check coordinate validity (start < end, positive coordinates)
- Validate chromosome naming convention matches genome assembly

**Peak Processing:**
- Remove duplicate peaks
- Filter peaks by size (typically 100-1000 bp)
- Merge overlapping peaks within 100bp

### 2. Differential Accessibility Analysis

**Methods:**
- Use BEDTools intersect to find overlapping peaks
- Calculate peak overlap statistics
- Identify gained peaks (present in B but not A)
- Identify lost peaks (present in A but not B)
- Calculate fold change: FC = (count_B + pseudocount) / (count_A + pseudocount)

**Statistical Testing:**
- Fisher's exact test for peak overlap significance
- Multiple testing correction using Benjamini-Hochberg
- Significance threshold: FDR < 0.05

### 3. Enhancer Activity Analysis

**H3K27ac Integration:**
- Intersect ATAC-seq peaks with H3K27ac peaks
- Define active enhancers as regions with both ATAC and H3K27ac signal
- Calculate enhancer activity score: log2(H3K27ac_signal + 1)

**Cell-type Specific Enhancers:**
- Identify enhancers specific to each cell type
- Calculate specificity score using Shannon entropy

### 4. Motif Enrichment Analysis

**Sequence Extraction:**
- Extract genomic sequences using BEDTools getfasta
- Use 500bp flanking regions for motif scanning
- Generate background sequences from random genomic regions

**Motif Scanning:**
- Use FIMO (Find Individual Motif Occurrences) from MEME Suite
- Parameters: p-value threshold 1e-4, output format TSV
- Scan both foreground (differential regions) and background sequences

**Enrichment Calculation:**
- Count motif occurrences in foreground vs. background
- Calculate enrichment score: E = (foreground_count / foreground_size) / (background_count / background_size)
- Statistical significance: Fisher's exact test
- Multiple testing correction: Benjamini-Hochberg FDR

### 5. Transcription Factor Prioritization

**Composite Scoring:**
- Motif enrichment score (log2 fold enrichment)
- Accessibility change score (log2 fold change)
- Enhancer activity score (H3K27ac signal)
- Expression change score (if RNA-seq available)

**Weighted Scoring Formula:**
```
Composite_Score = w1 * motif_enrichment + w2 * accessibility_change +
                  w3 * enhancer_activity + w4 * expression_change
```

**Default Weights:**
- motif_enrichment: 0.4
- accessibility_change: 0.3
- enhancer_activity: 0.2
- expression_change: 0.1

## Quality Control Metrics

### Data Quality
- Peak calling quality scores
- Sequence read depth
- Peak size distribution
- Chromosome coverage

### Analysis Quality
- Background sequence GC content matching
- Motif scanning sensitivity
- Statistical power estimation
- Reproducibility between replicates

## Parameter Settings

### BEDTools Parameters
- intersect: -f 0.5 -r (50% reciprocal overlap)
- subtract: -A (remove entire feature if overlap)
- merge: -d 100 (merge within 100bp)

### FIMO Parameters
- p-value threshold: 1e-4
- output format: TSV
- max-stored-scores: 100000

### Statistical Parameters
- FDR threshold: 0.05
- Fold change threshold: 1.5
- Pseudocount: 1

## References

1. Quinlan, A.R. and Hall, I.M., 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), pp.841-842.
2. Grant, C.E., Bailey, T.L. and Noble, W.S., 2011. FIMO: scanning for occurrences of a given motif. Bioinformatics, 27(7), pp.1017-1018.
3. Bailey, T.L., Johnson, J., Grant, C.E. and Noble, W.S., 2015. The MEME suite. Nucleic acids research, 43(W1), pp.W39-W49.