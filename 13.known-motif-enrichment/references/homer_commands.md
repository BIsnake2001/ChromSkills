# HOMER Known Motif Enrichment Command Reference

## Core Commands

### findMotifsGenome.pl
Main command for known motif enrichment analysis on genomic regions.

**Basic Syntax:**
```bash
findMotifsGenome.pl <peakfile> <genome> <output_dir> [options]
```

**Essential Options for Known Motif Analysis:**
- `-known`: Include known motif analysis
- `-nomotif`: Skip de novo motif finding (only known motifs)
- `-size <size>`: Region size around peak centers (default: 200)
- `-mask`: Mask repeat regions
- `-p <threads>`: Number of processors to use

**Custom Motif Options:**
- `-mknown <motif_file>`: Use custom motif file
- `-mcheck <motif_file>`: Check enrichment for specific motifs
- `-noknown`: Skip known motif analysis (opposite of -known)

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

## Common Workflow Examples

### Basic Known Motif Enrichment
```bash
findMotifsGenome.pl peaks.bed hg38 output_dir -known -nomotif -size 200 -p 8
```

### Custom Motif Database Analysis
```bash
findMotifsGenome.pl peaks.bed hg38 output_dir -mknown custom_motifs.motif -nomotif -size 200
```

### Specific Motif Checking
```bash
findMotifsGenome.pl peaks.bed hg38 output_dir -mcheck specific_tfs.motif -nomotif -size 200
```

### Comparative Analysis with Background
```bash
findMotifsGenome.pl treatment.bed hg38 output_dir -bg control.bed -known -nomotif -size 200 -p 8
```

### Combined De Novo and Known Analysis
```bash
findMotifsGenome.pl peaks.bed hg38 output_dir -known -size 200 -p 8
```

## Output Files Explained

### findMotifsGenome.pl Output
- `knownResults.html`: Known motif enrichment results (HTML format)
- `knownResults.txt`: Known motif enrichment results (tab-delimited)
- `seq.autonorm.tsv`: Sequence composition statistics
- `motifFindingParameters.txt`: Parameters used for analysis

### Known Results File Format
Tab-delimited columns in `knownResults.txt`:
- **Motif Name**: Name of the transcription factor motif
- **Consensus**: Consensus sequence
- **P-value**: Statistical significance
- **Log P-value**: -log10(p-value)
- **q-value (Benjamini)**: Multiple testing corrected p-value
- **% of Targets**: Percentage of input sequences with motif
- **% of Background**: Percentage of background sequences with motif
- **Fold Enrichment**: Target/background ratio

## Known Motif Databases

### Built-in Databases
HOMER includes several motif databases:
- **Vertebrate motifs**: Common vertebrate transcription factors
- **Mouse motifs**: Mouse-specific transcription factors
- **Human motifs**: Human-specific transcription factors
- **Fly motifs**: Drosophila transcription factors
- **Worm motifs**: C. elegans transcription factors

### Custom Motif Format
Create custom motif files in HOMER format:
```
>MOTIF_NAME
A [counts] C [counts] G [counts] T [counts]
```

Example:
```
>MY_TF
0.1 0.2 0.6 0.1
0.7 0.1 0.1 0.1
0.1 0.1 0.1 0.7
```

## Parameter Optimization Guidelines

### Region Size (-size)
- **Transcription factors**: 200-500 bp
- **Histone marks**: 1000-2000 bp
- **ATAC-seq**: 200-500 bp
- **DNase-seq**: 200-500 bp

### Threads (-p)
- **Standard**: 4-8 threads
- **High-performance**: 16-32 threads
- **Available cores**: Use `nproc` to check available cores

### Background Selection
- **Control peaks**: Use matched control samples
- **Random regions**: Use `-nlen` for random genomic background
- **GC matching**: Use `-gc` for GC-content matched background

## Interpretation Guidelines

### Significant Enrichment
- **p-value < 0.05**: Statistically significant
- **q-value < 0.1**: Multiple testing corrected significance
- **Fold enrichment > 2**: Biologically meaningful enrichment
- **% targets > 10%**: Substantial motif presence

### False Positives
- High GC-content motifs may show false enrichment
- Repeat-associated motifs may be artifacts
- Check background rates carefully

## Troubleshooting Common Issues

### No Known Motifs Found
- Verify genome assembly is supported
- Check if known motif databases are installed
- Try different region sizes
- Use `-list` to check available genomes

### Memory Errors
- Reduce region size (-size)
- Use fewer threads (-p)
- Increase system memory

### Slow Performance
- Use more threads (-p)
- Reduce region size (-size)
- Use faster storage (SSD)
- Check system load

### Genome Not Found
- Run `findMotifsGenome.pl -list` to check available genomes
- Verify genome installation
- Check genome naming conventions
- Download missing genome data

### Custom Motif Issues
- Verify motif file format
- Check motif file permissions
- Ensure motif file is accessible
- Validate motif position weight matrix format