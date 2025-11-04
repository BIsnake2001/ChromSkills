---
name: celltype-tf-driver
description: This skill should be used when users need to identify key driver transcription factors between two cell types using chromatin accessibility (ATAC-seq), histone modification (H3K27ac), and optional RNA-seq data. It processes BED files, genome FASTA files, and MEME motif databases to find transcription factors that drive cell type transitions.
---

# Celltype Transcription Factor Driver Analysis

## Overview

This skill implements a bioinformatics pipeline to identify key driver transcription factors (TFs) that regulate transitions between two cell types. The analysis integrates chromatin accessibility data (ATAC-seq), histone modification data (H3K27ac), transcription factor motif databases, and optionally RNA-seq data for differential expression analysis.

## Workflow Decision Tree

To determine the appropriate analysis path:

1. **Check available data files**
   - Required: ATAC_peak_{celltype_A}.bed, ATAC_peak_{celltype_B}.bed, H3K27ac_peak_{celltype_A}.bed, H3K27ac_peak_{celltype_B}.bed, {genome}.fa, {species}.meme
   - Optional: RNA-seq files for DEG analysis

2. **Choose analysis mode**
   - Basic: Chromatin accessibility + motif analysis only
   - Enhanced: Include differential expression analysis if RNA-seq data available

## Input Requirements

### Required Files
- `ATAC_peak_{celltype_A}.bed` - ATAC-seq peaks for cell type A
- `ATAC_peak_{celltype_B}.bed` - ATAC-seq peaks for cell type B
- `H3K27ac_peak_{celltype_A}.bed` - H3K27ac peaks for cell type A
- `H3K27ac_peak_{celltype_B}.bed` - H3K27ac peaks for cell type B
- `{genome}.fa` - Reference genome FASTA file
- `{species}.meme` - MEME format motif database

### Optional Files
- RNA-seq count matrices or BAM files for differential expression analysis

### Parameters
- Cell type A name (e.g., "fibroblast", "stem_cell")
- Cell type B name (e.g., "cardiomyocyte", "neuron")
- Genome assembly (e.g., "hg38", "mm10")
- Species (e.g., "human", "mouse")
- Output directory path

## Core Analysis Steps

### Step 1: Data Validation and Preparation

Validate input files and prepare data for analysis:

- Check that all required BED files exist and are properly formatted
- Verify genome FASTA and MEME database files
- Create output directory structure
- Generate summary statistics for each input file

### Step 2: Identify Gained Regions

Use `bedtools subtract` and `bedtools intersect` to identify *gained enhancers*

- `Gained_ATAC = ATAC_B - ATAC_A`
- `Gained_H3K27ac = H3K27ac_B - H3K27ac_A`
- `Gained_Enhancers = intersect(Gained_ATAC, Gained_H3K27ac)`

### Step 3: Transcription Factor Motif Analysis

Scan differential regions for TF binding motifs using HOMER

- Use HOMER findMotifsGenome.pl for motif enrichment analysis with proper background generation
- HOMER automatically handles GC-content normalization and background matching
- Calculate enrichment ratios and statistical significance with Benjamini-Hochberg correction
- Parse HOMER knownResults.txt to extract TF names, p-values, q-values, and enrichment metrics
- Fallback to FIMO analysis if HOMER fails for robust motif scanning
- Extract TF names from motif identifiers using pattern matching (e.g., "MA0004.1(AP2)" → "AP2")

### Step 4: Integration with Expression Data (Optional)

If RNA-seq data is available:

- Perform differential expression analysis
- Identify transcription factors with both motif enrichment and expression changes
- Prioritize TFs that show coordinated regulation

### Step 5: Key Driver Identification

Integrate all evidence to identify key driver TFs:

- Combine motif enrichment scores, accessibility changes, and expression data
- Calculate composite scores for each transcription factor
- Apply statistical thresholds and filtering
- Generate ranked list of candidate driver TFs

## Output Generation

### Primary Outputs
- `key_driver_tfs.txt` - List of identified key driver transcription factors with scores
- `motif_enrichment_results.csv` - Detailed motif enrichment statistics
- `differential_accessibility_results.bed` - Genomic coordinates of differential regions
- `enhancer_analysis_results.csv` - Enhancer activity analysis results

### Visualization Files
- `tf_enrichment_plot.png` - Bar plot of top enriched transcription factors
- `volcano_plot.png` - Volcano plot of accessibility changes vs. significance
- `heatmap.png` - Heatmap of TF activity across conditions

### Summary Reports
- `analysis_summary.md` - Comprehensive analysis summary
- `methodology_details.md` - Detailed description of methods used
- `quality_control_report.md` - QC metrics and data quality assessment

## Tools and Methods

### Required Bioinformatics Tools
- BEDTools for genomic interval operations and background generation
- HOMER (findMotifsGenome.pl) for motif enrichment with GC correction
- FIMO for fallback motif scanning when HOMER fails
- SAMtools/BCFtools for sequence handling and genome size extraction
- R/Bioconductor for statistical analysis and visualization

### Statistical Methods
- Fisher's exact test for motif enrichment
- Benjamini-Hochberg correction for multiple testing
- Fold change calculations for accessibility
- Enrichment ratio calculations (target % / background %)
- Log2 enrichment scores for ranking
- Composite scoring for driver identification

## Quality Control

### Data Quality Checks
- Verify BED file format and coordinate validity
- Check genome sequence completeness
- Validate MEME database format
- Assess peak calling quality metrics

### Analysis Quality Metrics
- Background region selection validation
- Motif scanning sensitivity assessment
- Statistical power estimation
- Reproducibility measures

## Troubleshooting

### Common Issues
- Missing or malformed BED files
- Genome sequence mismatches
- MEME database compatibility
- Memory limitations for large datasets

### Solutions
- Use BEDTools validate for file checking
- Verify genome assembly version matches
- Convert motif databases if needed
- Implement chunked processing for large files

## References

For detailed methodology and implementation details, refer to:
- `references/methodology.md` - Comprehensive analysis methods
- `references/tools_config.md` - Tool configurations and parameters
- `references/interpretation_guide.md` - Result interpretation guidelines