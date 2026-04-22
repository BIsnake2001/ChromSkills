# H3K4me3 and H3K27me3 ChIP-seq Analysis Report

## Analysis Summary

**Date**: April 16, 2026  
**Project Directory**: `/root/ChromOmics/check_cost/ChIPseq2`

Successfully identified H3K4me3 and H3K27me3 peaks from input BAM files and generated genome-wide signal tracks for IGV visualization.

## Input Data

### Biological Replicates
- **H3K4me3**: 2 replicates (wt_H3K4me3_rep1.bam, wt_H3K4me3_rep2.bam)
- **H3K27me3**: 2 replicates (wt_H3K27me3_rep1.bam, wt_H3K27me3_rep2.bam)
- **Input Control**: 2 replicates (wt_input_rep1.bam, wt_input_rep2.bam)

### Genome Assembly
- **Assembly**: Mouse (mm)
- **Total genome size**: ~2.73 Gb
- **Chromosome naming**: Numerical (1, 2, 3, ... X, Y, MT)

## Analysis Workflow

### Step 1: Initialize Analysis and Pool Replicates
**MCP Tool**: `mcp__project-init-tools__project_init`
- **Parameters**:
  - `sample`: H3K4me3, H3K27me3, input
  - `task`: rep_merge
- **Output Directories**:
  - `/root/ChromOmics/check_cost/ChIPseq2/H3K4me3_rep_merge`
  - `/root/ChromOmics/check_cost/ChIPseq2/H3K27me3_rep_merge`
  - `/root/ChromOmics/check_cost/ChIPseq2/input_rep_merge`

**MCP Tool**: `mcp__bw-tools__pool_bams`
- **Parameters for H3K4me3**:
  - `bam_files`: ["wt_H3K4me3_rep1.filtered.bam", "wt_H3K4me3_rep2.filtered.bam"]
  - `output_bam`: `H3K4me3_rep_merge/temp/H3K4me3.pooled.bam`
- **Parameters for H3K27me3**:
  - `bam_files`: ["wt_H3K27me3_rep1.filtered.bam", "wt_H3K27me3_rep2.filtered.bam"]
  - `output_bam`: `H3K27me3_rep_merge/temp/H3K27me3.pooled.bam`
- **Parameters for Input Control**:
  - `bam_files`: ["wt_input_rep1.filtered.bam", "wt_input_rep2.filtered.bam"]
  - `output_bam`: `input_rep_merge/temp/input.pooled.bam`

### Step 2: Peak Calling with MACS2

#### H3K4me3 (Narrow Peaks)
**MCP Tool**: `mcp__project-init-tools__project_init`
- **Parameters**: `sample`: H3K4me3, `task`: peak_calling
- **Output Directory**: `/root/ChromOmics/check_cost/ChIPseq2/H3K4me3_peak_calling`

**MCP Tool**: `mcp__macs2-tools__run_macs2`
- **Parameters**:
  - `treatment_file`: `H3K4me3_rep_merge/temp/H3K4me3.pooled.bam`
  - `control_file`: `input_rep_merge/temp/input.pooled.bam`
  - `genome_size`: mm (mouse)
  - `name`: H3K4me3
  - `out_dir`: `H3K4me3_peak_calling/peaks`
  - `broad`: false (narrow peaks)
  - `qvalue`: 0.05 (5% FDR)
  - `format`: BAMPE (paired-end)
- **Output Files**:
  - `H3K4me3_peaks.narrowPeak` - 664 peaks
  - `H3K4me3_peaks.xls` - Detailed peak information
  - `H3K4me3_summits.bed` - Peak summits

#### H3K27me3 (Broad Peaks)
**MCP Tool**: `mcp__project-init-tools__project_init`
- **Parameters**: `sample`: H3K27me3, `task`: peak_calling
- **Output Directory**: `/root/ChromOmics/check_cost/ChIPseq2/H3K27me3_peak_calling`

**MCP Tool**: `mcp__macs2-tools__run_macs2`
- **Parameters**:
  - `treatment_file`: `H3K27me3_rep_merge/temp/H3K27me3.pooled.bam`
  - `control_file`: `input_rep_merge/temp/input.pooled.bam`
  - `genome_size`: mm (mouse)
  - `name`: H3K27me3
  - `out_dir`: `H3K27me3_peak_calling/peaks`
  - `broad`: true (broad histone mark)
  - `broad_cutoff`: 0.1
  - `qvalue`: 0.05 (5% FDR)
  - `format`: BAMPE (paired-end)
- **Output Files**:
  - `H3K27me3_peaks.broadPeak` - 2,079 broad domains
  - `H3K27me3_peaks.gappedPeak` - Gapped peak format
  - `H3K27me3_peaks.xls` - Detailed peak information

### Step 3: Genome-wide Signal Track Generation

**MCP Tool**: `mcp__project-init-tools__project_init`
- **Parameters**: `sample`: all, `task`: track_generation
- **Output Directory**: `/root/ChromOmics/check_cost/ChIPseq2/all_track_generation`

**MCP Tool**: `mcp__bw-tools__generate_chrom_sizes`
- **Parameters**:
  - `bam_file`: `H3K4me3_rep_merge/temp/H3K4me3.pooled.bam`
  - `output_path`: `all_track_generation/temp/all.chrom.sizes`

**MCP Tool**: `mcp__bw-tools__calculate_scaling_factor`
- **Parameters for H3K4me3**: `bam_file`: `H3K4me3_rep_merge/temp/H3K4me3.pooled.bam`
  - **Scale factor**: 0.3067 (RPM normalization)
- **Parameters for H3K27me3**: `bam_file`: `H3K27me3_rep_merge/temp/H3K27me3.pooled.bam`
  - **Scale factor**: 0.0756 (RPM normalization)
- **Parameters for Input**: `bam_file`: `input_rep_merge/temp/input.pooled.bam`
  - **Scale factor**: 0.3174 (RPM normalization)

**MCP Tool**: `mcp__bw-tools__bam_to_bigwig`
- **Parameters for H3K4me3**:
  - `bam_file`: `H3K4me3_rep_merge/temp/H3K4me3.pooled.bam`
  - `chrom_sizes`: `all_track_generation/temp/all.chrom.sizes`
  - `output_bw`: `all_track_generation/tracks/H3K4me3.RPM.bw`
  - `scale_factor`: 0.3067
  - `shift_tn5`: false (ChIP-seq, no Tn5 shift)
- **Parameters for H3K27me3**:
  - `bam_file`: `H3K27me3_rep_merge/temp/H3K27me3.pooled.bam`
  - `chrom_sizes`: `all_track_generation/temp/all.chrom.sizes`
  - `output_bw`: `all_track_generation/tracks/H3K27me3.RPM.bw`
  - `scale_factor`: 0.0756
  - `shift_tn5`: false (ChIP-seq, no Tn5 shift)
- **Parameters for Input Control**:
  - `bam_file`: `input_rep_merge/temp/input.pooled.bam`
  - `chrom_sizes`: `all_track_generation/temp/all.chrom.sizes`
  - `output_bw`: `all_track_generation/tracks/input.RPM.bw`
  - `scale_factor`: 0.3174
  - `shift_tn5`: false (ChIP-seq, no Tn5 shift)

## Results Summary

### Peak Statistics
| Mark | Peak Type | Number of Peaks | Output File |
|------|-----------|----------------|-------------|
| H3K4me3 | Narrow | 664 | `H3K4me3_peak_calling/peaks/H3K4me3_peaks.narrowPeak` |
| H3K27me3 | Broad | 2,079 | `H3K27me3_peak_calling/peaks/H3K27me3_peaks.broadPeak` |

### Signal Tracks (BigWig)
| Sample | File Size | Normalization | File Path |
|--------|-----------|---------------|-----------|
| H3K4me3 | 20.8 MB | RPM (Reads Per Million) | `all_track_generation/tracks/H3K4me3.RPM.bw` |
| H3K27me3 | 89.3 MB | RPM (Reads Per Million) | `all_track_generation/tracks/H3K27me3.RPM.bw` |
| Input Control | 29.4 MB | RPM (Reads Per Million) | `all_track_generation/tracks/input.RPM.bw` |

### Quality Control Notes
- All BAM files were filtered and cleaned prior to analysis
- Paired-end sequencing data detected for all samples
- RPM normalization applied to all signal tracks
- Appropriate peak calling parameters used for each histone mark:
  - H3K4me3: Narrow peaks with summit calling
  - H3K27me3: Broad peaks with domain detection

## Output Directory Structure

```
/root/ChromOmics/check_cost/ChIPseq2/
├── all_bam_filtration/           # Filtered BAM files
├── H3K4me3_rep_merge/           # H3K4me3 replicate pooling
├── H3K27me3_rep_merge/          # H3K27me3 replicate pooling
├── input_rep_merge/             # Input control replicate pooling
├── H3K4me3_peak_calling/        # H3K4me3 peak calling results
│   └── peaks/
│       ├── H3K4me3_peaks.narrowPeak
│       ├── H3K4me3_peaks.xls
│       └── H3K4me3_summits.bed
├── H3K27me3_peak_calling/       # H3K27me3 peak calling results
│   └── peaks/
│       ├── H3K27me3_peaks.broadPeak
│       ├── H3K27me3_peaks.gappedPeak
│       └── H3K27me3_peaks.xls
├── all_track_generation/        # Signal track generation
│   ├── tracks/
│   │   ├── H3K4me3.RPM.bw
│   │   ├── H3K27me3.RPM.bw
│   │   └── input.RPM.bw
│   └── temp/
└── ChIPseq_analysis_report.md   # This report
```

## Usage Instructions for IGV

1. **Load Genome**: Use mouse genome (mm10) with appropriate chromosome naming (numerical)
2. **Load Signal Tracks**:
   - `all_track_generation/tracks/H3K4me3.RPM.bw`
   - `all_track_generation/tracks/H3K27me3.RPM.bw`
   - `all_track_generation/tracks/input.RPM.bw`
3. **Load Peak Annotations**:
   - `H3K4me3_peak_calling/peaks/H3K4me3_peaks.narrowPeak`
   - `H3K27me3_peak_calling/peaks/H3K27me3_peaks.broadPeak`

## Analysis Parameters Justification

1. **Genome Selection**: Mouse (mm) based on chromosome naming and size (~2.73 Gb)
2. **Peak Type Selection**:
   - H3K4me3: Active histone mark → narrow peaks
   - H3K27me3: Repressive histone mark → broad domains
3. **Statistical Threshold**: q-value = 0.05 (5% FDR) for both marks
4. **Normalization**: RPM (Reads Per Million) for consistent signal scaling
5. **Control Usage**: Pooled input control used for both peak calling analyses

## Notes
- No pseudo-replicate generation needed (only 2 biological replicates each)
- No IDR analysis performed (sufficient with 2 replicates and pooled analysis)
- Signal tracks are RPM-normalized for cross-sample comparison
- All intermediate files preserved for reproducibility

---
**Analysis completed successfully on April 16, 2026**