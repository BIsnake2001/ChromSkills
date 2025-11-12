---
name: ATACseq-QC
description: Perform comprehensive ATAC-seq quality control using ataqv, samtools, and IDR, automatically detecting genome type and generating an integrated PDF report with key metrics and visualizations.
---

# ATAC-seq Quality Control

## Overview

This skill performs complete ATAC-seq data quality control from BAM and peak files, automatically inferring genome type and generating visualized reports.

**Main Steps:**
1. Detect genome assembly (hg38, hg19, mm10, etc.) from BAM header.  
2. Collect basic alignment-level statistics using `samtools`.  
3. Assess fragment size distribution and library complexity.  
4. Compute TSS enrichment and FRiP using `ataqv`.  
5. Evaluate replicate reproducibility using `idr` and correlation metrics.  
6. Summarize all QC metrics and visualizations into a final PDF report in `qc_results/`.

---

## Decision Tree

### 1. Genome Detection and Setup
Automatically determine genome from BAM header (e.g., `SN:chr1`, `SN:1`).
```bash
genome=$(samtools view -H sample.bam | grep "@SQ" | head -1 | grep -oE 'SN:[^[:space:]]+' | head -1 | cut -d':' -f2)
# Map pattern to known genomes
# human (hg38/hg19) or mouse (mm10/mm9)
```
TSS annotation path is automatically assigned:
```
$(dirname $(which homer))/../data/genomes/${genome}/${genome}.tss
```

### 2. Alignment-Level QC (samtools)
Evaluate mapping quality, duplication, and mitochondrial reads.
```bash
samtools flagstat sample.bam > qc_results/sample.flagstat.txt
samtools idxstats sample.bam > qc_results/sample.idxstats.txt
```
**Metrics Extracted:**
- Total reads, mapped reads, properly paired (%)
- Mitochondrial fraction (MT%)
- Duplicate reads (if marked)

---

### 3. Fragment Size Distribution and Library Complexity
Assess fragment periodicity and complexity:
```bash
samtools view sample.bam | awk '{print $9}' | awk '{if($1>0 && $1<1000) print $1}' > qc_results/insert_sizes.txt
picard CollectInsertSizeMetrics I=sample.bam O=qc_results/insert_metrics.txt H=qc_results/insert_hist.pdf M=0.5
```
**Metrics:**
- Fragment length periodicity (mono/di/tri-nucleosome)
- NRF (unique reads / total reads)
- PBC (distinct read positions / unique reads)

---

### 4. TSS Enrichment and FRiP (ataqv)
Quantify transcription start site (TSS) enrichment and peak-associated signal.
```bash
ataqv --tss-file $(dirname $(which homer))/../data/genomes/${genome}/${genome}.tss       --metrics-file qc_results/sample_metrics.json       --peak-file peaks/sample.narrowPeak       ${genome} sample.bam
mkarv qc_results/sample_report.html qc_results/sample_metrics.json
```
**Metrics:**
- TSS enrichment score  
- FRiP (Fraction of Reads in Peaks)  
- Nucleosome-free vs mono-/di-nucleosome ratios  

---

### 5. Reproducibility Across Replicates (IDR + Correlation)
Perform reproducibility assessment between replicates.

**Example:**
```bash
idr --samples rep1_peaks.narrowPeak rep2_peaks.narrowPeak     --input-file-type narrowPeak     --rank p.value --output-file qc_results/rep1_rep2.idr.txt     --plot
```
**Outputs:**
- IDR value (≤0.05 indicates high reproducibility)
- Pearson/Spearman correlation between replicates

Correlation example:
```bash
multiBigwigSummary bins -b rep1.bw rep2.bw -o qc_results/summary.npz
plotCorrelation -in qc_results/summary.npz --corMethod pearson   --whatToPlot scatterplot -o qc_results/rep_corr.pdf
```

---

### 6. Generate Final QC Report (PDF)
All QC data (TSS enrichment, MT%, NRF, PBC, insert size, FRiP, IDR, correlation) are summarized into a single `qc_results/ATACseq_QC_Report.pdf`.

Example script call:
```bash
python scripts/generate_qc_report.py --input qc_results/                              --output qc_results/ATACseq_QC_Report.pdf                              --metrics TSS_enrichment MT_fraction NRF PBC FRiP IDR Correlation
```
Report contents:
- **Table**: Key QC metrics with thresholds  
- **Plots**: Insert size histogram, TSS profile, correlation heatmap, IDR plot  
- **Interpretation**: Automated text summary on data quality

---

## Quality Thresholds (Recommended)

| Metric | Excellent | Acceptable | Poor |
|---------|------------|-------------|------|
| TSS enrichment | >10 | 5–10 | <5 |
| Mitochondrial (%) | <20 | 20–50 | >50 |
| NRF | >0.8 | 0.5–0.8 | <0.5 |
| PBC | >0.9 | 0.5–0.9 | <0.5 |
| FRiP | >0.3 | 0.1–0.3 | <0.1 |
| IDR | <0.05 | 0.05–0.2 | >0.2 |
| Correlation | >0.9 | 0.7–0.9 | <0.7 |

---

## Output Structure

```
qc_results/
 ├── sample.flagstat.txt
 ├── sample.idxstats.txt
 ├── insert_hist.pdf
 ├── sample_metrics.json
 ├── sample_report.html
 ├── rep1_rep2.idr.txt
 ├── rep_corr.pdf
 └── ATACseq_QC_Report.pdf
```

---

## Notes

- The skill automatically detects genome type and uses matching TSS files.  
- If annotation is missing, it prompts the user to install the genome with `configureHomer.pl --install <genome>`.  
- Optional blacklist filtering can be added if available.  
- Intermediate files (temp results) can be safely removed after PDF report generation.