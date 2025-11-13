---
name: De-novo-motif-discovery
description: This skill identifies novel transcription factor binding motifs directly from genomic regions of interest such as ChIP-seq peaks, ATAC-seq accessible sites, or differentially acessible regions. It employs HOMER (Hypergeometric Optimization of Motif Enrichment) to detect both known and previously uncharacterized sequence motifs enriched within the supplied genomic intervals. Use the skill when you need to uncover sequence motifs enriched or want to know which TFs might regulate the target regions.
---

# HOMER De Novo Motif Discovery

## Overview

This skill enables comprehensive de novo motif discovery using HOMER tools for genomic peak files. It discovers novel transcription factor binding motifs from genomic regions without requiring prior knowledge of motif patterns. To perform de novo motif discovery:

- Always refer to the **Inputs & Outputs** section to check inputs and build the output architecture.
- Genome assembly: Always returned from user feedback (hg38, mm10, hg19, mm9, etc), never determined by yourself.
- Check chromosome names: Standardize chromosome names to format with "chr" (1 -> chr1, MT -> chrM).
- Set analysis parameters: Region size, number of motifs, motif lengths
- Run HOMER de novo motif discovery command

---

## When to use this skill
Use this skill when you need to uncover sequence motifs enriched in a set of genomic regions, such as peaks from ChIP-seq or ATAC-seq, without prior assumptions about which transcription factors are involved. Use this skill after generating a peak or region set that represents a meaningful biological comparison (e.g., differentially accessible or bound sites). The resulting motifs provide mechanistic insights into which TFs or sequence features might drive the observed chromatin or binding differences. Typical use cases include:

- Identifying unknown or cell-type–specific TF binding motifs in differential binding regions.
- Performing motif enrichment analysis on ATAC-seq peaks to infer potential transcriptional regulators of accessible chromatin regions.
- Exploring novel sequence patterns associated with enhancers, promoters, or other regulatory elements identified from functional genomics data.
- Building mechanistic hypotheses by linking enriched motifs with TFs, expression data, or chromatin state information.

---

## Inputs & Outputs

### Inputs
- **BED files**: Standard genomic interval format
- **narrowPeak**: ENCODE narrow peak format
- **broadPeak**: ENCODE broad peak format
- **HOMER peak files**: Output from HOMER peak calling

### Outputs
```bash
denovo_motif/
    results/
        homerResults.html # De novo motif discovery results
        motif<number>.motif # Individual motif files
        seq.autonorm.tsv # Sequence composition statistics
        motifFindingParameters.txt # Parameters used for analysis
    logs/ # analysis logs 
        motif.log
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
findMotifsGenome.pl <target_peaks> <genome> <output_dir> -size 200 -mask -p 8 -S 25 -len 8,10,12 -noknown > motif.log 2>&1
```

### Comparative Analysis

For comparing motif discovery between two sets of peaks:

```bash
findMotifsGenome.pl <target_peaks> <genome> <output_dir> -bg <background_peaks> -noknown > motif.log 2>&1
```

**Background options:**
- `-bg <background_file>`: Use specific background peaks
- `-nlen <length>`: Use random genomic regions as background
- `-gc`: Match GC content in background regions

**Example:**
```bash
findMotifsGenome.pl <target_peaks> <genome> <output_dir> -bg <control_peaks> -size 200 -p 8 -noknown > motif.log 2>&1
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
findMotifsGenome.pl <target_peaks> <genome> <output_dir> -len 6,8,10,12 -size 200 -p 8 -noknown > motif.log 2>&1
```

### High-throughput Discovery
```bash
findMotifsGenome.pl <target_peaks> <genome> <output_dir> -S 50 -p 16 -size 300 -noknown > motif.log 2>&1
```

### Comparative Discovery with GC-matched Background
```bash
findMotifsGenome.pl <target_peaks> <genome> <output_dir> -bg <control_peaks> -gc -size 200 -p 8 -noknown > motif.log 2>&1
```