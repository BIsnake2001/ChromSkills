---
name: De-novo-motif-discovery
description: This skill should be used when users need to perform de novo motif discovery on ChIP-seq, ATAC-seq, or other genomic peak files using HOMER (Hypergeometric Optimization of Motif EnRichment). It discovers novel transcription factor binding motifs from genomic regions without prior knowledge of motif patterns.
---

# HOMER De Novo Motif Discovery

## Overview

This skill enables comprehensive de novo motif discovery using HOMER tools for genomic peak files. It discovers novel transcription factor binding motifs from genomic regions without requiring prior knowledge of motif patterns.

## Quick Start

To perform de novo motif discovery:

1. **Identify input files**: BED, narrowPeak, or broadPeak files containing genomic regions. Detect and standardize chromosome names (chr1 ↔ 1, chrM ↔ MT) before analysis to match the HOMER genome format and avoid missing sequences.
2. **Determine genome assembly**: hg38, mm10, hg19, mm9, etc.
3. **Set analysis parameters**: Region size, number of motifs, motif lengths
4. **Run HOMER de novo motif discovery command**

## Core Capabilities

### De Novo Motif Discovery

Use `findMotifsGenome.pl` for discovering novel motifs in peak regions:

```bash
findMotifsGenome.pl <peakfile> <genome> <output_dir> [options]
```

**Essential Options:**
- `-size <size>`: Region size for motif finding (default: 200)
- `-mask`: Mask repeat regions
- `-p <threads>`: Number of processors to use
- `-S <num>`: Number of motifs to find (default: 25)
- `-len <length>`: Motif lengths to search (e.g., 8,10,12)

**Advanced Options:**
- `-cpg`: Enrich for CpG islands
- `-chopify`: Chop sequences into smaller fragments
- `-norevopp`: Don't search reverse complement
- `-rna`: For RNA motif finding
- `-bits`: Set information content threshold

**Example:**
```bash
findMotifsGenome.pl peaks.bed hg38 motif_output -size 200 -mask -p 8 -S 25 -len 8,10,12 -noknown
```

### Comparative Analysis

For comparing motif discovery between two sets of peaks:

```bash
findMotifsGenome.pl <target_peaks> <genome> <output_dir> -bg <background_peaks>  -noknown
```

**Background options:**
- `-bg <background_file>`: Use specific background peaks
- `-nlen <length>`: Use random genomic regions as background
- `-gc`: Match GC content in background regions

**Example:**
```bash
findMotifsGenome.pl treatment_peaks.bed hg38 motif_output -bg control_peaks.bed -size 200 -p 8  -noknown
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
4. **Repeat masking**: Use `-mask` for cleaner motif discovery

### Parameter Optimization
- **Region size**: Typically 200-500bp for transcription factors
- **Motif length**: 8-12bp for most transcription factors
- **Number of motifs**: 10-25 for initial discovery
- **Threads**: Use available CPU cores for faster processing

## Output Interpretation

### Key Output Files
- `homerResults.html`: De novo motif discovery results
- `motif<number>.motif`: Individual motif files
- `seq.autonorm.tsv`: Sequence composition statistics
- `motifFindingParameters.txt`: Parameters used for analysis

### Important Metrics
- **p-value**: Statistical significance of motif enrichment
- **% of targets**: Percentage of input sequences containing motif
- **% of background**: Percentage of background sequences containing motif
- **Log P-value**: -log10(p-value) for visualization
- **Information content**: Measure of motif specificity

## Troubleshooting

### Common Issues
1. **Memory errors**: Reduce region size or number of motifs
2. **Slow performance**: Use `-p` option for parallel processing
3. **No motifs found**: Check input file format and region size
4. **Genome not found**: Verify genome assembly name and installation

### Error Handling
- Ensure HOMER is properly installed and configured
- Check that genome data is downloaded and accessible
- Verify input file formats and chromosome naming
- Ensure sufficient disk space for output files

## Advanced Workflows

### Multi-length Motif Discovery
```bash
findMotifsGenome.pl peaks.bed hg38 output -len 6,8,10,12 -size 200 -p 8 -noknown
```

### High-throughput Discovery
```bash
findMotifsGenome.pl peaks.bed hg38 output -S 50 -p 16 -size 300 -noknown
```

### Comparative Discovery with GC-matched Background
```bash
findMotifsGenome.pl treatment.bed hg38 output -bg control.bed -gc -size 200 -p 8 -noknown
```