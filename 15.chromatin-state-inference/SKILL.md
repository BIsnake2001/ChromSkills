---
name: cchromatin-state-inference
description: This skill should be used when users need to infer chromatin states from histone modification ChIP-seq data using chromHMM. It provides workflows for chromatin state segmentation, model training, state annotation, and comparative analysis across multiple samples or conditions.
---

# ChromHMM Chromatin State Inference

## Overview

This skill enables comprehensive chromatin state analysis using chromHMM for histone modification ChIP-seq data. ChromHMM uses a multivariate Hidden Markov Model to segment the genome into discrete chromatin states based on combinatorial patterns of histone modifications.

## Quick Start

To perform chromatin state inference:

1. **Prepare input data**: Aligned BAM or BED files for histone modifications (H3K4me3, H3K27ac, H3K27me3, H3K36me3, etc.)
2. **Determine genome assembly**: hg38, mm10, hg19, mm9, etc.
3. **Set analysis parameters**: Number of states, bin size, genome regions
4. **Run chromHMM workflow**: Binarization → Learning → Segmentation → Annotation
5. **Interpret results**: State characterization, genomic feature enrichment, visualization

## Core Workflows

### 1. Data Preparation and Binarization

Convert BAM or BED files to chromHMM binarized format:

```bash
ChromHMM BinarizeBam \
  -b 200 \
  CHROMSIZES_FILE \
  INPUT_DIR \
  CELLMARKFILE \
  OUTPUT_DIR
```

```bash
ChromHMM BinarizeBed \
  -b 200 \
  CHROMSIZES_FILE \
  INPUT_DIR \
  CELLMARKFILE \
  OUTPUT_DIR
```

**Parameters:**
- `-b 200`: Bin size (200bp recommended)
- `CHROMSIZES_FILE`: Genome chromosome sizes file
- `INPUT_DIR`: Directory containing BAM files
- `CELLMARKFILE`: Cell mark file defining histone modifications
- `OUTPUT_DIR`: Output directory for binarized data

### 2. Model Learning

Train chromatin state model:

```bash
ChromHMM LearnModel \
  -p 4 \
  BINARIZED_DIR \
  OUTPUT_MODEL_DIR \
  15 \
  ASSEMBLY
```

**Parameters:**
- `-p 4`: Number of processors
- `BINARIZED_DIR`: Directory with binarized data
- `OUTPUT_MODEL_DIR`: Output directory for model
- `15`: Number of chromatin states
- `ASSEMBLY`: Genome assembly (hg38, mm10, etc.)

### 3. State Segmentation

Generate genome segmentation:

```bash
ChromHMM MakeSegmentation \
  MODEL_FILE \
  BINARIZED_DIR \
  OUTPUT_SEGMENTATION_DIR
```

### 4. State Annotation and Enrichment

Annotate states with genomic features:

```bash
ChromHMM OverlapEnrichment \
  SEGMENTATION_FILE \
  ANNOTATION_FILE \
  OUTPUT_ENRICHMENT_FILE
```

## File Format Handling

### Input Formats
- **BAM files**: Aligned reads for histone modifications
- **Cell mark file**: Tab-delimited file defining marks and cell types
- **Chromosome sizes**: Two-column file with chromosome names and sizes

### Output Formats
- **Emission parameters**: State-specific histone mark probabilities
- **Transition parameters**: State transition probabilities
- **Segmentation BED**: Genome-wide chromatin state assignments
- **Enrichment files**: Genomic feature overlap statistics

## Parameter Optimization

### Number of States
- **8 states**: Basic chromatin states
- **15 states**: Standard comprehensive states
- **25 states**: High-resolution states
- **Optimization**: Use Bayesian Information Criterion (BIC)

### Bin Size
- **200bp**: Standard resolution
- **100bp**: High resolution (requires more memory)
- **500bp**: Low resolution (faster computation)

## State Interpretation

### Common Chromatin States
1. **Active Promoter**: H3K4me3, H3K27ac
2. **Weak Promoter**: H3K4me3
3. **Poised Promoter**: H3K4me3, H3K27me3
4. **Strong Enhancer**: H3K27ac, H3K4me1
5. **Weak Enhancer**: H3K4me1
6. **Insulator**: CTCF
7. **Transcribed**: H3K36me3
8. **Repressed**: H3K27me3
9. **Heterochromatin**: Low signal across marks

## Quality Control

### Model Assessment
- Check convergence of learning algorithm
- Verify emission parameters are biologically meaningful
- Assess state transition probabilities
- Validate with known genomic annotations

### State Validation
- Compare with ENCODE chromatin states
- Validate with gene expression data
- Check enrichment at known functional elements
- Assess reproducibility across replicates

## Multi-sample Analysis

### Comparative Analysis
Compare chromatin states across conditions:

```bash
ChromHMM CompareModels \
  MODEL1_DIR \
  MODEL2_DIR \
  OUTPUT_COMPARISON_DIR
```

### Differential States
Identify condition-specific chromatin states:

```bash
ChromHMM DifferentialEnrichment \
  SEGMENTATION1 \
  SEGMENTATION2 \
  ANNOTATION_FILE \
  OUTPUT_DIFF_FILE
```

## Visualization

### State Emission Heatmaps
Visualize histone mark patterns for each state:

```bash
ChromHMM MakeHeatmap \
  SEGMENTATION_FILE \
  REFERENCE_POINTS \
  OUTPUT_HEATMAP_DIR
```

### Genome Browser Tracks
Generate tracks for visualization:

```bash
ChromHMM MakeBrowserFiles \
  SEGMENTATION_FILE \
  ASSEMBLY \
  OUTPUT_TRACKS_DIR
```

## Troubleshooting

### Common Issues
1. **Memory errors**: Reduce bin size or number of states
2. **Convergence problems**: Increase iterations or adjust learning rate
3. **Uninterpretable states**: Check input data quality and mark combinations
4. **Missing chromosomes**: Verify chromosome naming consistency

### Error Handling
- Ensure Java memory allocation is sufficient
- Verify all input files are properly formatted
- Check chromosome names match genome assembly
- Validate mark file contains correct paths to BAM files
