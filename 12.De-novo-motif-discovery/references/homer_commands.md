# HOMER De Novo Motif Discovery Command Reference

## Core Commands

### findMotifsGenome.pl
Main command for de novo motif discovery on genomic regions.

**Basic Syntax:**
```bash
findMotifsGenome.pl <peakfile> <genome> <output_dir> [options]
```

**Essential Options for De Novo Discovery:**
- `-size <size>`: Region size around peak centers (default: 200)
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

**Background Options:**
- `-bg <background_file>`: Use specific background peaks
- `-nlen <length>`: Use random genomic regions as background
- `-gc`: Match GC content in background regions

## File Conversion Commands

### pos2bed.pl
Convert MACS2 peak files to BED format:
```bash
pos2bed.pl macs2_peaks.xls > peaks.bed
```

### peak2bed.pl
Convert HOMER peak files to BED format:
```bash
peak2bed.pl homer_peaks.txt > peaks.bed
```

## Utility Commands

### findMotifs.pl
Motif analysis for non-genomic sequences:
```bash
findMotifs.pl <sequence_file> fasta <output_dir> [options]
```

## Common Workflow Examples

### Basic De Novo Motif Discovery
```bash
findMotifsGenome.pl peaks.bed hg38 output_dir -size 200 -mask -p 8
```

### Comprehensive Motif Discovery
```bash
findMotifsGenome.pl peaks.bed hg38 output_dir -size 300 -S 30 -len 8,10,12 -p 16
```

### Comparative Analysis with Background
```bash
findMotifsGenome.pl treatment.bed hg38 output_dir -bg control.bed -size 200 -p 8
```

### GC-matched Background Analysis
```bash
findMotifsGenome.pl peaks.bed hg38 output_dir -bg background.bed -gc -size 200 -p 8
```

## Output Files Explained

### findMotifsGenome.pl Output
- `homerResults.html`: De novo motif discovery results
- `motif<number>.motif`: Individual motif files
- `seq.autonorm.tsv`: Sequence composition statistics
- `motifFindingParameters.txt`: Parameters used for analysis

## Parameter Optimization Guidelines

### Region Size (-size)
- **Transcription factors**: 200-500 bp
- **Histone marks**: 1000-2000 bp
- **ATAC-seq**: 200-500 bp
- **DNase-seq**: 200-500 bp

### Motif Length (-len)
- **Most TFs**: 8-12 bp
- **Longer motifs**: 12-20 bp for complex factors
- **Multiple lengths**: 8,10,12 for comprehensive search

### Number of Motifs (-S)
- **Initial discovery**: 10-25 motifs
- **Comprehensive**: 25-50 motifs
- **Memory considerations**: Reduce if encountering memory errors

### Threads (-p)
- **Standard**: 4-8 threads
- **High-performance**: 16-32 threads
- **Available cores**: Use `nproc` to check available cores

## Troubleshooting Common Issues

### Memory Errors
- Reduce region size (-size)
- Reduce number of motifs (-S)
- Use fewer threads (-p)
- Increase system memory

### Slow Performance
- Use more threads (-p)
- Reduce region size (-size)
- Use faster storage (SSD)
- Check system load

### No Motifs Found
- Check input file format
- Increase region size (-size)
- Try different motif lengths (-len)
- Verify genome assembly

### Genome Not Found
- Run `findMotifsGenome.pl -list` to check available genomes
- Verify genome installation
- Check genome naming conventions
- Download missing genome data