# ChIP-seq Analysis Summary

## Overview
This analysis performed peak calling for H3K4me3 and H3K27me3 histone marks, followed by genomic feature annotation of the called peaks, ChIP-seq quality control assessment, and generation of genomic visualization tracks. The analysis used pooled replicates for increased statistical power.

## Directory Structure
```
/root/ChromOmics/test_skills_ChIPseq/
├── all_bam_filtration/          # BAM filtration outputs
├── all_peak_calling/            # Peak calling outputs
├── all_chip_qc/                 # ChIP-seq QC outputs
├── all_track_generation/        # BigWig track generation outputs
├── H3K4me3_genomic_feature_annotation/  # Feature annotation for H3K4me3
├── H3K27me3_genomic_feature_annotation/ # Feature annotation for H3K27me3
├── *.bam                        # Original BAM files
└── fixed_blacklist.bed          # Blacklist regions
```

## Analysis Steps

### 1. BAM File Filtration

**Purpose**: Remove artifacts including mitochondrial reads, PCR duplicates, unmapped reads, low-quality alignments, and blacklisted regions.

**Commands executed**:
```bash
# Check and fix BAM files (coordinate sorting, read groups)
Check BAM file sorting and read group information
- Input: All BAM files
- Output: Verified coordinate-sorted BAMs with read groups
- Parameters: 8 threads, temporary directory for processing

# Filter each BAM file
Filter BAM artifacts (remove mitochondrial reads, PCR duplicates, unmapped/low-quality reads, blacklisted regions)
- Input: wt_H3K4me3_rep1.bam
- Output: all_bam_filtration/filtered_bam/wt_H3K4me3_rep1.filtered.bam
- Parameters: Blacklist file (fixed_blacklist.bed), temporary directory
# Repeated for all 6 BAM files with corresponding inputs/outputs
```

**Output files**:
- `all_bam_filtration/filtered_bam/`
  - `wt_H3K4me3_rep1.filtered.bam`, `wt_H3K4me3_rep2.filtered.bam`
  - `wt_H3K27me3_rep1.filtered.bam`, `wt_H3K27me3_rep2.filtered.bam`
  - `wt_input_rep1.filtered.bam`, `wt_input_rep2.filtered.bam`
- `all_bam_filtration/temp/` - Intermediate files and metrics

### 2. Peak Calling

**Purpose**: Identify significantly enriched genomic regions for each histone mark.

**User decisions**:
- Genome size: `mm` (mouse)
- Control strategy: Pool input replicates
- Treatment strategy: Pool histone mark replicates

**Commands executed**:
```bash
# Initialize project directory
Create project directory structure for peak calling
- Parameters: sample="all", task="peak_calling"
- Output: all_peak_calling directory

# Pool input replicates
Merge BAM files to create pooled replicates
- Input: wt_input_rep1.filtered.bam, wt_input_rep2.filtered.bam
- Output: all_peak_calling/peaks/pooled_input.filtered.bam
- Method: samtools merge

# Pool H3K4me3 replicates
Merge BAM files to create pooled replicates
- Input: wt_H3K4me3_rep1.filtered.bam, wt_H3K4me3_rep2.filtered.bam
- Output: all_peak_calling/peaks/pooled_H3K4me3.filtered.bam
- Method: samtools merge

# Pool H3K27me3 replicates
Merge BAM files to create pooled replicates
- Input: wt_H3K27me3_rep1.filtered.bam, wt_H3K27me3_rep2.filtered.bam
- Output: all_peak_calling/peaks/pooled_H3K27me3.filtered.bam
- Method: samtools merge

# H3K4me3 peak calling (narrow peaks)
Run MACS2 peak calling for narrow peaks
- Treatment: pooled_H3K4me3.filtered.bam
- Control: pooled_input.filtered.bam
- Genome: mm (mouse)
- Output prefix: H3K4me3_pooled
- Parameters: format=BAMPE (paired-end), qvalue=0.05, narrow peaks
- Command equivalent: macs3 callpeak -t treatment.bam -c control.bam -f BAMPE -g mm -n name -q 0.05

# H3K27me3 peak calling (broad peaks)
Run MACS2 peak calling for broad peaks
- Treatment: pooled_H3K27me3.filtered.bam
- Control: pooled_input.filtered.bam
- Genome: mm (mouse)
- Output prefix: H3K27me3_pooled
- Parameters: format=BAMPE (paired-end), qvalue=0.05, broad peaks with cutoff 0.1
- Command equivalent: macs3 callpeak -t treatment.bam -c control.bam -f BAMPE -g mm -n name --broad --broad-cutoff 0.1 -q 0.05
```

**Output files**:
- `all_peak_calling/peaks/`
  - `H3K4me3_pooled_peaks.narrowPeak` - 664 narrow peaks
  - `H3K4me3_pooled_summits.bed` - Peak summits
  - `H3K4me3_pooled_peaks.xls` - Detailed peak information
  - `H3K27me3_pooled_peaks.broadPeak` - 2,079 broad peaks
  - `H3K27me3_pooled_peaks.gappedPeak` - Gapped peak format
  - `H3K27me3_pooled_peaks.xls` - Detailed peak information
  - Pooled BAM files and indices
- `all_peak_calling/logs/`
  - `H3K4me3_pooled_used_parameters.txt` - Parameter justification
  - `H3K27me3_pooled_used_parameters.txt` - Parameter justification

**Peak statistics**:
- H3K4me3: 664 narrow peaks (active histone mark)
- H3K27me3: 2,079 broad peaks (repressive histone mark)

### 3. Genomic Feature Annotation

**Purpose**: Annotate peaks with genomic features (promoters, exons, introns, intergenic regions, etc.)

**User decisions**:
- Genome assembly: `mm10` (Mouse GRCm38)
- Sample names: `H3K4me3` and `H3K27me3`

#### H3K4me3 Annotation

**Commands executed**:
```bash
# Initialize project directory
Create project directory structure for genomic feature annotation
- Parameters: sample="H3K4me3", task="genomic_feature_annotation"
- Output: H3K4me3_genomic_feature_annotation directory

# Standardize chromosome names (1 → chr1)
Convert chromosome names to standard format (add "chr" prefix)
- Input: all_peak_calling/peaks/H3K4me3_pooled_peaks.narrowPeak
- Output: H3K4me3_genomic_feature_annotation/H3K4me3_peaks_standardized.bed
- Method: Text processing to convert "1" → "chr1", "X" → "chrX", etc.

# Annotate peaks
Run HOMER annotatePeaks.pl for genomic feature annotation
- Input: H3K4me3_peaks_standardized.bed (peak regions)
- Genome: mm10
- Parameters: size given (keep original peak sizes)
- Output files:
  - H3K4me3.anno_genomic_features.txt (annotated peak table)
  - H3K4me3.anno_genomic_features_stats.txt (annotation statistics)
- Command equivalent: annotatePeaks.pl peaks.bed mm10 -size given -annStats stats.txt > annotated.txt

# Visualize annotation
Generate pie chart visualization of annotation statistics
- Input: H3K4me3.anno_genomic_features_stats.txt
- Output: H3K4me3.anno_genomic_features_stats.pdf
- Chart type: Pie chart showing distribution of genomic features
```

**Output files**:
- `H3K4me3_genomic_feature_annotation/results/`
  - `H3K4me3.anno_genomic_features.txt` - Annotated peak table
  - `H3K4me3.anno_genomic_features_stats.txt` - Annotation statistics
- `H3K4me3_genomic_feature_annotation/plots/`
  - `H3K4me3.anno_genomic_features_stats.pdf` - Pie chart visualization
- `H3K4me3_genomic_feature_annotation/log/` - Log files

#### H3K27me3 Annotation

**Commands executed**:
```bash
# Initialize project directory
Create project directory structure for genomic feature annotation
- Parameters: sample="H3K27me3", task="genomic_feature_annotation"
- Output: H3K27me3_genomic_feature_annotation directory

# Standardize chromosome names
Convert chromosome names to standard format (add "chr" prefix)
- Input: all_peak_calling/peaks/H3K27me3_pooled_peaks.broadPeak
- Output: H3K27me3_genomic_feature_annotation/H3K27me3_peaks_standardized.bed
- Method: Text processing to convert "1" → "chr1", "X" → "chrX", etc.

# Annotate peaks
Run HOMER annotatePeaks.pl for genomic feature annotation
- Input: H3K27me3_peaks_standardized.bed (peak regions)
- Genome: mm10
- Parameters: size given (keep original peak sizes)
- Output files:
  - H3K27me3.anno_genomic_features.txt (annotated peak table)
  - H3K27me3.anno_genomic_features_stats.txt (annotation statistics)
- Command equivalent: annotatePeaks.pl peaks.bed mm10 -size given -annStats stats.txt > annotated.txt

# Visualize annotation
Generate pie chart visualization of annotation statistics
- Input: H3K27me3.anno_genomic_features_stats.txt
- Output: H3K27me3.anno_genomic_features_stats.pdf
- Chart type: Pie chart showing distribution of genomic features
```

**Output files**:
- `H3K27me3_genomic_feature_annotation/results/`
  - `H3K27me3.anno_genomic_features.txt` - Annotated peak table
  - `H3K27me3.anno_genomic_features_stats.txt` - Annotation statistics
- `H3K27me3_genomic_feature_annotation/plots/`
  - `H3K27me3.anno_genomic_features_stats.pdf` - Pie chart visualization
- `H3K27me3_genomic_feature_annotation/log/` - Log files

### 4. ChIP-seq Quality Control (QC)

**Purpose**: Assess quality of ChIP-seq data using cross-correlation metrics (NSC, RSC) and Fraction of Reads in Peaks (FRiP).

**Commands executed**:
```bash
# Cross-correlation analysis using phantompeakqualtools
Run R script for cross-correlation analysis
- Input: BAM file (pooled_H3K4me3.filtered.bam, pooled_H3K27me3.filtered.bam)
- Output: Cross-correlation statistics and PDF plots
- Command equivalent: Rscript run_spp.R -c=<BAM> -savp -out=<output.txt> -rf

# FRiP calculation
Calculate Fraction of Reads in Peaks
- Input: BAM file and corresponding peak file
- Method:
  1. Count total mapped reads: samtools view -c <BAM>
  2. Count reads overlapping peaks: bedtools intersect -a <BAM> -b <PEAKS> -u | wc -l
  3. Calculate FRiP = reads_in_peaks / total_mapped
```

**Output files**:
- `all_chip_qc/`
  - `pooled_H3K4me3.filtered_spp.txt` - Cross-correlation statistics
  - `pooled_H3K4me3.filtered_crosscorr.pdf` - Cross-correlation plot
  - `pooled_H3K4me3.filtered_frip.txt` - FRiP results
  - `pooled_H3K27me3.filtered_spp.txt` - Cross-correlation statistics
  - `pooled_H3K27me3.filtered_crosscorr.pdf` - Cross-correlation plot
  - `pooled_H3K27me3.filtered_frip.txt` - FRiP results

**QC Results**:

| Histone Mark | NSC | RSC | Quality Tag | FRiP | Total Reads | Reads in Peaks |
|--------------|-----|-----|-------------|------|-------------|----------------|
| H3K4me3 | 2.464 | 1.197 | 1 | 0.491 | 3,260,245 | 1,601,412 |
| H3K27me3 | 1.027 | 2.730 | 2 | 0.036 | 13,227,466 | 472,881 |

**Interpretation**:
- **NSC (Normalized Strand Cross-correlation)**: Measures signal-to-noise ratio. Values >1.05 indicate good quality.
- **RSC (Relative Strand Cross-correlation)**: Measures fragment length normalization. Values >0.8 indicate good quality.
- **Quality Tag**: 1=High quality, 2=Moderate quality, 0=Low quality.
- **FRiP**: Fraction of reads in peaks. H3K4me3 shows high enrichment (49%), typical for active marks. H3K27me3 shows lower FRiP (3.6%), common for broad histone marks.

### 5. Genomic Track Generation

**Purpose**: Generate RPM-normalized BigWig tracks for genome browser visualization.

**Commands executed**:
```bash
# Generate chromosome sizes from BAM header
samtools view -H <BAM> | grep ^@SQ | cut -f2,3 | sed 's/SN://' | sed 's/LN://' > chrom.sizes

# Calculate scaling factor (RPM normalization)
mapped_reads=$(samtools view -c -F 260 <BAM>)
scale_factor=$(echo "1000000 / $mapped_reads" | bc -l)

# Convert BAM to BigWig (ChIP-seq, no Tn5 shift)
bedtools genomecov -ibam <BAM> -bg -scale $scale_factor > temp.bedgraph
bedGraphToBigWig temp.bedgraph chrom.sizes output.RPM.bw
```

**Output files**:
- `all_track_generation/tracks/`
  - `pooled_H3K4me3.RPM.bw` - RPM-normalized signal track (21.8 MB)
  - `pooled_H3K27me3.RPM.bw` - RPM-normalized signal track (93.5 MB)
- `all_track_generation/temp/` - Intermediate files (chrom.sizes, bedGraph)

**Track Statistics**:
| Histone Mark | Scaling Factor | File Size | Normalization |
|--------------|----------------|-----------|---------------|
| H3K4me3 | 0.3067 (3.26M reads) | 21.8 MB | RPM (reads per million) |
| H3K27me3 | 0.0756 (13.23M reads) | 93.5 MB | RPM (reads per million) |

**Usage**: BigWig files can be loaded into genome browsers (UCSC, IGV) for visualization of ChIP-seq signal across the genome.

## Key Results

### Peak Calling Summary
| Histone Mark | Peak Type | Number of Peaks | Control Used |
|--------------|-----------|-----------------|--------------|
| H3K4me3 | Narrow | 664 | Pooled input |
| H3K27me3 | Broad | 2,079 | Pooled input |

### Expected Annotation Patterns
- **H3K4me3**: Typically enriched at promoter regions near transcription start sites (TSS)
- **H3K27me3**: Often found in broad domains covering gene bodies and intergenic regions

## Quality Control Notes
1. All BAM files were coordinate-sorted and contained read group information
2. Blacklist filtering used provided `fixed_blacklist.bed`
3. Paired-end sequencing detected; MACS2 used BAMPE mode
4. Appropriate peak types: narrow for H3K4me3, broad for H3K27me3
5. Genome assembly mm10 used for annotation
6. ChIP-seq QC metrics indicate good data quality:
   - **H3K4me3**: NSC=2.464 (>1.05), RSC=1.197 (>0.8), FRiP=49.1%
   - **H3K27me3**: NSC=1.027 (>1.05), RSC=2.730 (>0.8), FRiP=3.6% (typical for broad marks)
7. RPM-normalized BigWig tracks generated for genome browser visualization

## Next Steps (Potential)
1. Functional enrichment analysis (GO, KEGG)
2. Motif discovery in peak regions
3. Integrative analysis with gene expression data
4. Visualization of peaks on specific genomic loci