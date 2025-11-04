# Alignment QC Metrics and Thresholds for ChIP-seq/ATAC-seq

## Overview

This document provides comprehensive guidelines for assessing alignment quality of ChIP-seq and ATAC-seq BAM files. Use these metrics and thresholds to determine whether BAM files are suitable for downstream peak calling analysis.

## Core QC Metrics

### 1. Total Reads
**Description**: Total number of reads in the BAM file
**Interpretation**: Indicates sequencing depth

**Thresholds**:
- **ChIP-seq**:
  - PASS: ≥10 million reads (for TFs), ≥20 million reads (for histone marks)
  - WARN: 5-10 million reads (TFs), 10-20 million reads (histones)
  - FAIL: <5 million reads (TFs), <10 million reads (histones)
- **ATAC-seq**:
  - PASS: ≥50 million reads
  - WARN: 25-50 million reads
  - FAIL: <25 million reads

### 2. Mapped Reads Percentage
**Description**: Percentage of reads successfully aligned to the reference genome
**Interpretation**: Overall alignment efficiency

**Thresholds**:
- **PASS**: ≥90%
- **WARN**: 80-90%
- **FAIL**: <80%

### 3. Properly Paired Reads (Paired-end Only)
**Description**: Percentage of read pairs that are properly oriented and mapped
**Interpretation**: Library quality and mapping accuracy

**Thresholds**:
- **PASS**: ≥85%
- **WARN**: 70-85%
- **FAIL**: <70%

### 4. Duplicate Rate
**Description**: Percentage of PCR duplicates
**Interpretation**: Library complexity and potential over-amplification

**Thresholds**:
- **PASS**: ≤20%
- **WARN**: 20-40%
- **FAIL**: >40%

### 5. Mitochondrial Reads Percentage
**Description**: Percentage of reads mapping to mitochondrial DNA
**Interpretation**: Nuclear enrichment efficiency

**Thresholds**:
- **PASS**: ≤5%
- **WARN**: 5-15%
- **FAIL**: >15%

### 6. Insert Size Metrics
**Description**: Average insert size and distribution
**Interpretation**: Library fragment size distribution

**Thresholds**:
- **ChIP-seq**:
  - PASS: 100-300 bp (typical range)
  - WARN: Outside typical range but consistent
  - FAIL: Highly variable or inappropriate for target
- **ATAC-seq**:
  - PASS: Clear nucleosome-free (<100 bp) and mononucleosome (~200 bp) peaks
  - WARN: Weak periodicity
  - FAIL: No clear fragment size distribution

## ChIP-seq Specific Considerations

### Transcription Factor (TF) ChIP-seq
- **Read depth**: 10-20 million mapped reads typically sufficient
- **Duplicate rate**: More tolerant of higher duplicates (up to 30%) if library complexity is good
- **Fragment size**: Should match expected TF binding site size

### Histone Mark ChIP-seq
- **Read depth**: 20-40 million mapped reads recommended
- **Broad marks (H3K27me3, H3K36me3)**: Require higher depth (40+ million)
- **Sharp marks (H3K4me3, H3K27ac)**: 20-30 million typically sufficient

## ATAC-seq Specific Considerations

### Library Complexity
- **Non-redundant fraction**: Should be >0.5 for good complexity
- **TSS enrichment**: Should show clear enrichment at transcription start sites
- **Fragment size distribution**: Should show clear nucleosome-free and nucleosome-bound peaks

### Quality Metrics
- **FRiP score** (Fraction of reads in peaks):
  - PASS: ≥0.2
  - WARN: 0.1-0.2
  - FAIL: <0.1
- **TSS enrichment**:
  - PASS: ≥8
  - WARN: 4-8
  - FAIL: <4

## Peak-calling Readiness Assessment

### PASS Criteria
All of the following must be met:
1. Mapped reads ≥90%
2. Duplicate rate ≤20%
3. Mitochondrial reads ≤5%
4. Read depth meets minimum requirements for experiment type
5. Proper pairing ≥85% (if paired-end)

### WARN Criteria
One or more of the following:
1. Mapped reads 80-90%
2. Duplicate rate 20-40%
3. Mitochondrial reads 5-15%
4. Read depth borderline for experiment type

### FAIL Criteria
One or more of the following:
1. Mapped reads <80%
2. Duplicate rate >40%
3. Mitochondrial reads >15%
4. Insufficient read depth for experiment type
5. Severe technical artifacts

## Troubleshooting Common Issues

### Low Mapping Rate
- Check reference genome compatibility
- Verify sequencing quality (FastQC)
- Consider adapter contamination

### High Duplicate Rate
- May indicate low library complexity
- Consider downsampling if necessary
- Evaluate whether duplicates should be removed

### High Mitochondrial Reads
- Indicates poor nuclear enrichment
- For ATAC-seq, may require more stringent nuclei isolation
- Consider computational filtering if enrichment was suboptimal

### Abnormal Insert Size Distribution
- Check library preparation protocol
- Verify size selection steps
- Consider whether distribution matches experimental expectations

## MultiQC Integration

MultiQC automatically collects and visualizes these metrics:
- **General Statistics**: Summary table of key metrics
- **Samtools Flagstat**: Mapping statistics
- **Samtools Stats**: Detailed alignment metrics
- **Insert Size**: Fragment size distributions
- **Coverage**: Read depth across genome

Use the MultiQC report to quickly identify samples that require closer inspection and to compare quality across multiple samples.

## Final Recommendations

### For PASS Samples
- Proceed directly to peak calling
- No additional processing needed

### For WARN Samples
- Consider duplicate removal if duplicate rate is high
- Evaluate whether additional sequencing is needed
- May proceed to peak calling with caution

### For FAIL Samples
- Do not proceed to peak calling
- Investigate and address underlying issues
- Consider re-sequencing if problems cannot be resolved computationally