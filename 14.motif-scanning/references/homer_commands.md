# HOMER Motif Scanning Command Reference

## Core Commands

### annotatePeaks.pl
Main command for motif scanning and peak annotation.

**Basic Syntax:**
```bash
annotatePeaks.pl <peakfile> <genome> [options] > output.txt
```

**Essential Motif Scanning Options:**
- `-m <motif_file>`: Scan for specific motifs
- `-size <size>`: Region size around peak centers
- `-nmotifs`: Number of motifs to report per peak
- `-mbed`: Output in BED format with motif locations
- `-mscore`: Include motif scores in output

**Advanced Options:**
- `-cpu <threads>`: Number of processors for parallel processing
- `-bedGraph`: Output in bedGraph format
- `-hist <bins>`: Include histone modification data
- `-d <tag_directory>`: Include tag density information
- `-dfile <file>`: Include custom data tracks

### scanMotifGenomeWide.pl
Command for genome-wide motif scanning.

**Basic Syntax:**
```bash
scanMotifGenomeWide.pl <motif_file> <genome> [options] > output.bed
```

**Options:**
- `-bed`: Output in BED format
- `-p <threads>`: Number of processors to use
- `-mask`: Mask repeat regions
- `-keepAll`: Keep all matches (not just best per position)
- `-score <threshold>`: Set score threshold
- `-match <sequence>`: Match specific sequence

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

### meme2meme
Convert MEME motif format to HOMER format:
```bash
meme2meme meme_output/meme.txt -homer > homer_motifs.motif
```

## Common Workflow Examples

### Basic Motif Scanning
```bash
annotatePeaks.pl peaks.bed hg38 -m known.motifs -size 200 -nmotifs 3 > motif_hits.txt
```

### BED Format Output
```bash
annotatePeaks.pl peaks.bed hg38 -m ${tf}.motif -size 200 -mbed > motif_locations.bed
```

### Genome-wide Scanning
```bash
scanMotifGenomeWide.pl ${tf}.motif hg38 -bed -p 8 > genome_scan.bed
```

### Score-based Filtering
```bash
annotatePeaks.pl peaks.bed hg38 -m ${tf}.motif -size 200 -mscore | awk '$NF > 8.0' > high_score_hits.txt
```

### Multi-motif Scanning
```bash
annotatePeaks.pl peaks.bed hg38 -m multiple_motifs.motif -size 200 -nmotifs 5 > all_hits.txt
```

### Strand-specific Analysis
```bash
annotatePeaks.pl peaks.bed hg38 -m ${tf}.motif -size 200 | grep "+" > forward_strand_hits.txt
```

## Output Files Explained

### annotatePeaks.pl Output
Tab-delimited columns in standard output:
- **PeakID**: Unique identifier for each peak
- **chr**: Chromosome
- **start**: Start position
- **end**: End position
- **strand**: DNA strand
- **Peak Score**: Peak score from input file
- **Focus Ratio/Region Size**: Additional peak metrics
- **Annotation**: Genomic feature annotation
- **Detailed Annotation**: Detailed genomic context
- **Distance to TSS**: Distance to nearest transcription start site
- **Nearest PromoterID**: ID of nearest promoter
- **Entrez ID**: Entrez gene ID
- **Nearest Unigene**: Nearest Unigene cluster
- **Gene Name**: Gene symbol
- **Gene Alias**: Gene aliases
- **Gene Description**: Gene description
- **Gene Type**: Gene biotype
- **Motif Hits**: Motif match information (when using -m)

### BED Format Output (-mbed)
Standard BED format columns plus:
- **name**: Motif name and score
- **score**: Motif match score
- **strand**: DNA strand of motif match
- **thickStart/End**: Motif match coordinates
- **itemRgb**: Color coding
- **blockCount**: Number of blocks
- **blockSizes**: Block sizes
- **blockStarts**: Block start positions

### bedGraph Format Output (-bedGraph)
Continuous score tracks:
- **chr**: Chromosome
- **start**: Start position
- **end**: End position
- **score**: Motif score

## Motif File Formats

### HOMER Motif Format
Motif file for all TFs is ${HOMER_data}/knownTFs/all.motifs and motif file for known TFs are in ${HOMER_data}/knownTFs/known.motifs. Separate motif files are under ${HOMER_data}/knownTFs/motifs.
Position-specific scoring matrix format:
```
>MOTIF_NAME
A [counts] C [counts] G [counts] T [counts]
A [counts] C [counts] G [counts] T [counts]
...
```

Example CTCF motif:
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

### Supported Formats
- **HOMER format**: Native HOMER position-specific scoring matrices
- **MEME format**: MEME suite motif format
- **TRANSFAC format**: TRANSFAC database format
- **JASPAR format**: JASPAR database format

## Parameter Optimization Guidelines

### Region Size (-size)
- **Transcription factors**: 200-500 bp
- **Histone marks**: 1000-2000 bp
- **ATAC-seq**: 200-500 bp
- **DNase-seq**: 200-500 bp

### Number of Motifs (-nmotifs)
- **Standard**: 1-3 motifs per peak
- **Comprehensive**: 5-10 motifs per peak
- **Memory considerations**: Higher numbers increase memory usage

### Score Thresholds
- **High confidence**: > 8.0 bits
- **Medium confidence**: 6.0-8.0 bits
- **Low confidence**: < 6.0 bits
- **Motif-specific**: Adjust based on motif information content

### Threads (-cpu / -p)
- **Standard**: 4-8 threads
- **High-performance**: 16-32 threads
- **Available cores**: Use `nproc` to check available cores

## Interpretation Guidelines

### Motif Score Interpretation
- **> 10 bits**: Very high confidence match
- **8-10 bits**: High confidence match
- **6-8 bits**: Medium confidence match
- **< 6 bits**: Low confidence match

### Position Analysis
- **Peak centers**: Motifs near peak centers are more likely functional
- **Strand bias**: Some motifs show strand-specific enrichment
- **Cluster analysis**: Multiple motifs in same peak may indicate cooperativity

## Troubleshooting Common Issues

### No Motif Hits Found
- Check motif file format
- Verify region size is appropriate
- Check motif score thresholds
- Verify genome assembly

### Memory Errors
- Reduce region size (-size)
- Use fewer threads (-cpu)
- Reduce number of motifs reported (-nmotifs)
- Increase system memory

### Slow Performance
- Use more threads (-cpu)
- Reduce region size (-size)
- Use faster storage (SSD)
- Check system load

### Genome Not Found
- Run `findMotifsGenome.pl -list` to check available genomes
- Verify genome installation
- Check genome naming conventions
- Download missing genome data

### Motif File Issues
- Verify motif file format
- Check motif file permissions
- Ensure motif file is accessible
- Validate position weight matrix format

### Output Format Issues
- Check file permissions for output
- Ensure sufficient disk space
- Verify output redirection syntax
- Check for special characters in file names
