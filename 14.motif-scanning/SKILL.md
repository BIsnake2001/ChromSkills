---
name: motif-scanning
description: This skill identifies the locations of known transcription factor (TF) binding motifs within genomic regions such as ChIP-seq or ATAC-seq peaks. It utilizes HOMER to search for specific sequence motifs defined by position-specific scoring matrices (PSSMs) from known motif databases. Use this skill when you need to detect the presence and precise genomic coordinates of known TF binding motifs within experimentally defined regions such as ChIP-seq or ATAC-seq peaks.
---

# Motif Scanning

## Overview

This skill enables comprehensive motif scanning using HOMER tools for genomic peak files. It scans genomic regions for specific transcription factor binding motifs using position-specific scoring matrices and identifies exact motif locations.To perform motif scanning:

- Always refer to the **Inputs & Outputs** section to check inputs and build the output architecture.
- Genome assembly: Always returned from user feedback (hg38, mm10, hg19, mm9, etc), never determined by yourself.
- Check chromosome names: Standardize chromosome names to format with "chr" (1 -> chr1, MT -> chrM).
- Prepare motif files: Position-specific scoring matrices (PSSM) in HOMER format, saved in ${HOMER_data}/knownTFs/motifs/${tf}.motif, and "tf" should be in lower case.
4. Set scanning parameters: Region size, score thresholds, output format
5. Run HOMER motif scanning command

---

## When to use this skill

- Scanning ChIP-seq or ATAC-seq peaks for known motifs to validate TF binding specificity.
- Testing whether co-factor motifs (e.g., TAL1, KLF1, SPI1) co-occur within TF-bound or accessible regions to infer cooperative binding.
- Evaluating motif distribution patterns relative to genomic landmarks such as transcription start sites (TSS) or enhancers.
- Generating motif-annotated BED files for visualization in genome browsers or subsequent feature analysis.

---

## Inputs & Outputs

### Inputs
(1) Peak formats supported
- **BED files**: Standard genomic interval format
- **narrowPeak**: ENCODE narrow peak format
- **broadPeak**: ENCODE broad peak format
- **HOMER peak files**: Output from HOMER peak calling
(2) Motif formats supported
- **HOMER motif format**: Position-specific scoring matrices
- **MEME motif format**: MEME suite motif format
- **TRANSFAC format**: TRANSFAC database format

Convert from MEME to HOMER format:

```bash
meme2meme meme_output/meme.txt -homer > homer_motifs.motif
```

### Outputs
```bash
motif_scan/
    results/
      ${motif}.anno # Peak annotations with motif hits
      ${motif}.bed # Exact motif locations in genomic coordinates
      ${motif}.bedgraph # Continuous motif score tracks
    logs/ # analysis logs 
        motif_scan.log
    temp/ # other temp files
```
---

## Decision Tree

### Standardize chromosome names 

From `1` format to `chr1` format
From `MT` format to `chrM` format

If the chromosome name not startswith "chr", then run:

```bash
awk '{print "chr"$1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' <peakfile> > <new_peakfile> 
```

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
scanMotifGenomeWide.pl ${HOMER_data}/knownTFs/motifs/ctcf.motif hg38 -bed > genome_scan.bed
```

### BED Format Output

Generate BED files with exact motif locations:

```bash
annotatePeaks.pl peaks.bed hg38 -m ${HOMER_data}/knownTFs/motifs/ctcf.motif -size 200 -mbed > motif_locations.bed
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

### Important Metrics
- **Motif score**: Position-specific scoring matrix match score
- **Position**: Exact genomic location of motif match
- **Strand**: DNA strand where motif was found
- **Sequence**: Actual DNA sequence at motif location


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
annotatePeaks.pl peaks.bed hg38 -m ${HOMER_data}/knownTFs/known.motifs -size 200 -nmotifs 5 > all_hits.txt
```

### BED Output for Visualization
```bash
annotatePeaks.pl peaks.bed hg38 -m ${HOMER_data}/knownTFs/motifs/ctcf.motif -size 200 -mbed > motif_locations.bed
```

### Genome-wide Scanning
```bash
scanMotifGenomeWide.pl ${HOMER_data}/knownTFs/motifs/ctcf.motif hg38 -bed -p 8 > genome_scan.bed
```

### Score-based Filtering
```bash
annotatePeaks.pl peaks.bed hg38 -m ${HOMER_data}/knownTFs/motifs/ctcf.motif -size 200 -mscore | awk '$NF > 8.0' > high_score_hits.txt
```

### Strand-specific Analysis
```bash
annotatePeaks.pl peaks.bed hg38 -m ${HOMER_data}/knownTFs/motifs/ctcf.motif -size 200 | grep "+" > forward_strand_hits.txt
```