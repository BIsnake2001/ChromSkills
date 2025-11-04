---
name: motif-scanning
description: This skill should be used when users need to perform motif scanning on ChIP-seq, ATAC-seq, or other genomic peak files using HOMER (Hypergeometric Optimization of Motif EnRichment). It scans genomic regions for specific transcription factor binding motifs using position-specific scoring matrices and identifies exact motif locations.
---

# HOMER Motif Scanning

## Overview

This skill enables comprehensive motif scanning using HOMER tools for genomic peak files. It scans genomic regions for specific transcription factor binding motifs using position-specific scoring matrices and identifies exact motif locations.

## Quick Start

To perform motif scanning:

1. **Identify input files**: BED, narrowPeak, or broadPeak files containing genomic regions. Detect and standardize chromosome names (chr1 ↔ 1, chrM ↔ MT) before analysis to match the HOMER genome format and avoid missing sequences.
2. **Determine genome assembly**: hg38, mm10, hg19, mm9, etc.
3. **Prepare motif files**: Position-specific scoring matrices (PSSM) in HOMER format
4. **Set scanning parameters**: Region size, score thresholds, output format
5. **Run HOMER motif scanning command**

## Core Capabilities

### Motif Scanning with annotatePeaks.pl

Use `annotatePeaks.pl` to scan for specific motifs in peak regions:

```bash
annotatePeaks.pl <peakfile> <genome> -m <motif_file> -size <size> [options] > output.txt
```

**Essential Options:**
- `-m <motif_file>`: Motif file to scan for
- `-size <size>`: Region size around peak centers
- `-nmotifs`: Number of motifs to report per peak
- `-mbed`: Output in BED format with motif locations
- `-mscore`: Include motif scores in output

**Advanced Options:**
- `-cpu <threads>`: Number of processors for parallel processing
- `-bedGraph`: Output in bedGraph format
- `-hist <bins>`: Include histone modification data

**Example:**
```bash
annotatePeaks.pl peaks.bed hg38 -m known.motifs -size 200 -nmotifs 3 > motif_hits.txt
```

### Genome-wide Motif Scanning

Use `scanMotifGenomeWide.pl` for genome-wide motif scanning:

```bash
scanMotifGenomeWide.pl <motif_file> <genome> [options] > output.bed
```

**Options:**
- `-bed`: Output in BED format
- `-p <threads>`: Number of processors to use
- `-mask`: Mask repeat regions
- `-keepAll`: Keep all matches (not just best per position)

**Example:**
```bash
scanMotifGenomeWide.pl motif.motif hg38 -bed > genome_scan.bed
```

### BED Format Output

Generate BED files with exact motif locations:

```bash
annotatePeaks.pl peaks.bed hg38 -m motif.motif -size 200 -mbed > motif_locations.bed
```

## File Format Handling

### Supported Input Formats
- **BED files**: Standard genomic interval format
- **narrowPeak**: ENCODE narrow peak format
- **broadPeak**: ENCODE broad peak format
- **HOMER peak files**: Output from HOMER peak calling

### Motif File Formats
- **HOMER motif format**: Position-specific scoring matrices
- **MEME motif format**: MEME suite motif format
- **TRANSFAC format**: TRANSFAC database format

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
3. **Motif quality**: Use high-quality position-specific scoring matrices
4. **Score thresholds**: Set appropriate motif score cutoffs

### Parameter Optimization
- **Region size**: Typically 200-500bp for transcription factors
- **Number of motifs**: Report top 1-5 motifs per peak
- **Score thresholds**: Use default or optimize based on motif quality
- **Threads**: Use available CPU cores for faster processing

## Output Interpretation

### Key Output Files
- **Tab-delimited output**: Peak annotations with motif hits
- **BED format**: Exact motif locations in genomic coordinates
- **bedGraph format**: Continuous motif score tracks

### Important Metrics
- **Motif score**: Position-specific scoring matrix match score
- **Position**: Exact genomic location of motif match
- **Strand**: DNA strand where motif was found
- **Sequence**: Actual DNA sequence at motif location

## Motif File Preparation

### HOMER Motif Format
Create motif files in HOMER position-specific scoring matrix format:

```
>MOTIF_NAME
A [counts] C [counts] G [counts] T [counts]
A [counts] C [counts] G [counts] T [counts]
...
```

**Example:**
```
>CTCF
0.173 0.348 0.348 0.130
0.161 0.339 0.375 0.125
0.179 0.330 0.366 0.125
0.155 0.330 0.384 0.131
0.118 0.360 0.360 0.162
0.118 0.360 0.360 0.162
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
0.125 0.348 0.384 0.143
```

### Converting Motif Formats
Convert from MEME to HOMER format:
```bash
meme2meme meme_output/meme.txt -homer > homer_motifs.motif
```

## Troubleshooting

### Common Issues
1. **No motif hits found**: Check motif file format and region size
2. **Memory errors**: Reduce region size or use fewer threads
3. **Slow performance**: Use `-cpu` option for parallel processing
4. **Genome not found**: Verify genome assembly name and installation

### Error Handling
- Ensure HOMER is properly installed and configured
- Check that genome data is downloaded and accessible
- Verify input file formats and chromosome naming
- Ensure motif files are in correct format
- Check sufficient disk space for output files

## Advanced Workflows

### Multi-motif Scanning
```bash
annotatePeaks.pl peaks.bed hg38 -m multiple_motifs.motif -size 200 -nmotifs 5 > all_hits.txt
```

### BED Output for Visualization
```bash
annotatePeaks.pl peaks.bed hg38 -m motif.motif -size 200 -mbed > motif_locations.bed
```

### Genome-wide Scanning
```bash
scanMotifGenomeWide.pl motif.motif hg38 -bed -p 8 > genome_scan.bed
```

### Score-based Filtering
```bash
annotatePeaks.pl peaks.bed hg38 -m motif.motif -size 200 -mscore | awk '$NF > 8.0' > high_score_hits.txt
```

### Strand-specific Analysis
```bash
annotatePeaks.pl peaks.bed hg38 -m motif.motif -size 200 | grep "+" > forward_strand_hits.txt
```