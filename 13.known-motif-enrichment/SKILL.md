---
name: known-motif-enrichment
description: This skill should be used when users need to perform known motif enrichment analysis on ChIP-seq, ATAC-seq, or other genomic peak files using HOMER (Hypergeometric Optimization of Motif EnRichment). It identifies enrichment of known transcription factor binding motifs from established databases in genomic regions.
---

# HOMER Known Motif Enrichment

## Overview

This skill enables comprehensive known motif enrichment analysis using HOMER tools for genomic peak files. It identifies enrichment of known transcription factor binding motifs from established databases in genomic regions.

## Quick Start

To perform known motif enrichment analysis:

1. **Identify input files**: BED, narrowPeak, or broadPeak files containing genomic regions. Detect and standardize chromosome names (chr1 ↔ 1, chrM ↔ MT) before analysis to match the HOMER genome format and avoid missing sequences.
2. **Determine genome assembly**: hg38, mm10, hg19, mm9, etc.
3. **Choose analysis type**: Standard known motif analysis or custom motif checking
4. **Run HOMER known motif enrichment command**

## Core Capabilities

### Known Motif Enrichment

Use `findMotifsGenome.pl` with known motif databases:

```bash
findMotifsGenome.pl <peakfile> <genome> <output_dir> -known
```

**Essential Options:**
- `-known`: Include known motif analysis
- `-nomotif`: Skip de novo motif finding, only do known motifs
- `-size <size>`: Region size for motif finding (default: 200)
- `-mask`: Mask repeat regions
- `-p <threads>`: Number of processors to use

**Custom Motif Options:**
- `-mknown <motif_file>`: Use custom motif file
- `-mcheck <motif_file>`: Check enrichment for specific motifs

**Example:**
```bash
findMotifsGenome.pl peaks.bed hg38 motif_output -known -nomotif -size 200 -p 8
```

### Custom Motif Enrichment

For enrichment analysis with custom motif databases:

```bash
findMotifsGenome.pl <peakfile> <genome> <output_dir> -mknown <motif_file> -nomotif
```

**Example:**
```bash
findMotifsGenome.pl peaks.bed hg38 motif_output -mknown custom_motifs.motif -nomotif -size 200
```

### Specific Motif Checking

For checking enrichment of specific motifs:

```bash
findMotifsGenome.pl <peakfile> <genome> <output_dir> -mcheck <motif_file> -nomotif
```

**Example:**
```bash
findMotifsGenome.pl peaks.bed hg38 motif_output -mcheck specific_motifs.motif -nomotif -size 200
```

### Comparative Analysis

For comparing known motif enrichment between two sets of peaks:

```bash
findMotifsGenome.pl <target_peaks> <genome> <output_dir> -bg <background_peaks> -known -nomotif
```

**Example:**
```bash
findMotifsGenome.pl treatment_peaks.bed hg38 motif_output -bg control_peaks.bed -known -nomotif -size 200
```

## File Format Handling

### Supported Input Formats
- **BED files**: Standard genomic interval format
- **narrowPeak**: ENCODE narrow peak format
- **broadPeak**: ENCODE broad peak format
- **HOMER peak files**: Output from HOMER peak calling

### Converting Between Formats

Convert MACS2 output to HOMER format:
```bash
pos2bed.pl macs2_peaks.xls > peaks.bed
```

Convert narrowPeak to BED:
```bash
cut -f1-6 narrowPeak_file.narrowPeak > peaks.bed
```

## Genome Assembly Support

HOMER supports multiple genome assemblies. Common assemblies include:
- **Human**: hg38, hg19, hg18
- **Mouse**: mm10, mm9
- **Other**: dm6 (fly), ce10 (worm), rn6 (rat)

To check available genomes:
```bash
findMotifsGenome.pl -list
```

## Quality Control and Best Practices

### Pre-processing Steps
1. **Filter peaks**: Remove low-quality or artifact peaks
2. **Size selection**: Use appropriate region size (-size parameter)
3. **Background selection**: Choose appropriate background for enrichment analysis
4. **Repeat masking**: Use `-mask` for cleaner motif analysis

### Parameter Optimization
- **Region size**: Typically 200-500bp for transcription factors
- **Threads**: Use available CPU cores for faster processing
- **Background**: Use appropriate control peaks for comparative analysis

## Output Interpretation

### Key Output Files
- `knownResults.html`: Known motif enrichment results
- `knownResults.txt`: Tab-delimited known motif results
- `seq.autonorm.tsv`: Sequence composition statistics
- `motifFindingParameters.txt`: Parameters used for analysis

### Important Metrics
- **p-value**: Statistical significance of motif enrichment
- **% of targets**: Percentage of input sequences containing motif
- **% of background**: Percentage of background sequences containing motif
- **Log P-value**: -log10(p-value) for visualization
- **Fold enrichment**: Ratio of target vs background motif occurrence

## Known Motif Databases

### Built-in Databases
HOMER includes several built-in motif databases:
- **Vertebrate motifs**: Common transcription factors
- **Mouse motifs**: Mouse-specific transcription factors
- **Human motifs**: Human-specific transcription factors
- **Fly motifs**: Drosophila transcription factors
- **Worm motifs**: C. elegans transcription factors

### Custom Databases
Users can create custom motif databases in HOMER format:
```
>MOTIF_NAME
A [counts] C [counts] G [counts] T [counts]
```

## Troubleshooting

### Common Issues
1. **No known motifs found**: Check if genome assembly is supported
2. **Memory errors**: Reduce region size
3. **Slow performance**: Use `-p` option for parallel processing
4. **Genome not found**: Verify genome assembly name and installation

### Error Handling
- Ensure HOMER is properly installed and configured
- Check that genome data is downloaded and accessible
- Verify input file formats and chromosome naming
- Ensure sufficient disk space for output files

## Advanced Workflows

### Combined De Novo and Known Motif Analysis
```bash
findMotifsGenome.pl peaks.bed hg38 output -known -size 200 -p 8
```

### Custom Motif Database Analysis
```bash
findMotifsGenome.pl peaks.bed hg38 output -mknown custom_db.motif -nomotif -size 200
```

### Specific Transcription Factor Analysis
```bash
findMotifsGenome.pl peaks.bed hg38 output -mcheck tf_motifs.motif -nomotif -size 200
```

### Comparative Analysis with Control Peaks
```bash
findMotifsGenome.pl treatment.bed hg38 output -bg control.bed -known -nomotif -size 200
```