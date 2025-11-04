---
name: ChIPseq-QC
description: This skill performs comprehensive ChIP-seq quality control using phantompeakqualtools, samtools, and IDR. It generates QC metrics (NSC, RSC, NRF, PBC, FRiP, Insert size, Correlation, and IDR), visualizes key results, and outputs a complete PDF summary in `qc_results/`.
---

# Comprehensive ChIP-seq QC Pipeline

## Overview

This skill performs a full ChIP-seq quality control analysis from aligned BAM files.  
The pipeline includes the following main steps:

1. **Compute basic alignment statistics** using `samtools flagstat` and `samtools idxstats`.  
2. **Estimate library complexity** and calculate **NRF/PBC** using `phantompeakqualtools` (from `spp`).  
3. **Perform cross-correlation analysis** to calculate **NSC** and **RSC** using `phantompeakqualtools`.  
4. **Inspect fragment length distribution** (Insert size) using `picard CollectInsertSizeMetrics` (or samtools).  
5. **Compute FRiP (Fraction of Reads in Peaks)** using peak files and aligned BAMs.  
6. **Evaluate replicate reproducibility** using **IDR** between replicates.  
7. **Assess signal correlation** between replicates (Spearman/Pearson).  
8. **Generate a final PDF report** summarizing all metrics, plots, and interpretations.  

All intermediate and result files are saved in the directory `qc_results/`.

---

## Decision Tree

### Step 1. Alignment Summary with samtools
Assess mapping efficiency and basic alignment quality.

**Command example:**
```bash
samtools flagstat sample.bam > qc_results/sample.flagstat.txt
samtools idxstats sample.bam > qc_results/sample.idxstats.txt
```
**Key metrics:** total reads, mapped reads, properly paired reads, duplication rate.

---

### Step 2. Library Complexity (NRF, PBC)
Compute **Non-Redundant Fraction (NRF)** and **PCR Bottlenecking Coefficients (PBC)** using `phantompeakqualtools`.

**Command example:**
```bash
run_spp.R -c=sample.bam -savp -out=qc_results/sample_crosscorr.txt
```
**Definitions:**  
- **NRF = unique reads / total reads**  
- **PBC = distinct read pairs / total read pairs**  

**Recommended thresholds:** NRF > 0.8; PBC > 0.8

---

### Step 3. Cross-Correlation Metrics (NSC, RSC)
Calculate **NSC** and **RSC** via phantompeakqualtools.

**Command example:**
```bash
run_spp.R -c=sample.bam -p=8 -savp -out=qc_results/sample_spp.txt
```
**Thresholds:** NSC > 1.05 (acceptable), > 1.1 (good); RSC > 0.8 (acceptable), > 1.0 (good)  
Outputs include `.pdf` plot and `.txt` metrics in `qc_results/`.

---

### Step 4. Insert Size Distribution
Inspect fragment length distribution using Picard (or `samtools stats`).

**Command example (Picard):**
```bash
picard CollectInsertSizeMetrics       I=sample.bam       O=qc_results/sample.insert_metrics.txt       H=qc_results/sample.insert_histogram.pdf       M=0.5
```
Interpret the modal insert size (TFs often ~150–250 bp; histone marks broader).

---

### Step 5. FRiP (Fraction of Reads in Peaks)
Calculate the fraction of reads falling within peak regions.

**Command example:**
```bash
bedtools intersect -u -a sample.bam -b peaks.bed | samtools view -c - > qc_results/sample.inpeak.count
samtools view -c sample.bam > qc_results/sample.total.count
awk '{printf \"%.6f\n\", $1/$2}' qc_results/sample.inpeak.count qc_results/sample.total.count > qc_results/sample.frip.txt
```
**Suggested thresholds:** TFs: FRiP > 0.01 (acceptable), > 0.05 (good); Histone marks: FRiP > 0.1 (good).

---

### Step 6. Replicate Reproducibility (IDR)
Use **IDR** to evaluate reproducibility between replicates.

**Command example:**
```bash
idr --samples rep1_peaks.narrowPeak rep2_peaks.narrowPeak     --input-file-type narrowPeak     --rank p.value     --plot     --output-file qc_results/replicates.idr.txt     --log-output-file qc_results/replicates.idr.log
```
**Interpretation:** IDR < 0.05 → strong; 0.05–0.1 → moderate; > 0.1 → weak.

---

### Step 7. Replicate Correlation
Compute correlation between replicates using deepTools or R.

**Example (deepTools):**
```bash
multiBamSummary bins --bamfiles rep1.bam rep2.bam -o qc_results/multibam.npz
plotCorrelation -in qc_results/multibam.npz --corMethod spearman     --whatToPlot heatmap --skipZeros     -o qc_results/replicates_correlation.pdf     --outFileCorMatrix qc_results/replicates_correlation.tsv
```
**Rule of thumb:** correlation coefficient > 0.8 indicates good reproducibility.

---

### Step 8. Final PDF Report
Generate a summarized PDF report combining all metrics and plots.

The report includes:
- A summary table for each sample (NSC, RSC, NRF, PBC, FRiP, IDR, Correlation)
- QC plots (cross-correlation, insert size, correlation heatmap)
- Interpretations and pass/fail indicators

**Example:**
```bash
Rscript scripts/compile_chipseq_qc_report.R qc_results/ > qc_results/compile.log 2>&1
# Output: qc_results/ChIPseq_QC_Report.pdf
```

---

## Output Structure

```
qc_results/
├── *.flagstat.txt
├── *.idxstats.txt
├── *spp.txt, *crosscorr.txt, *crosscorr.pdf
├── *.insert_metrics.txt, *.insert_histogram.pdf
├── *.frip.txt
├── replicates.idr.txt, replicates.idr.log, replicates_correlation.pdf
└── ChIPseq_QC_Report.pdf
```

---

## Quality Standards Summary

| Metric | Description | Tool | Recommended Threshold |
|---------|--------------|------|------------------------|
| NSC | Normalized Strand Cross-correlation | phantompeakqualtools | >1.05 (ok), >1.10 (good) |
| RSC | Relative Strand Cross-correlation | phantompeakqualtools | >0.8 (ok), >1.0 (good) |
| NRF | Non-Redundant Fraction | spp | >0.8 |
| PBC | PCR Bottleneck Coefficient | spp | >0.8 |
| Insert Size | Fragment length distribution | Picard/samtools | TFs modal ~150–250 bp |
| FRiP | Fraction of reads in peaks | bedtools | TFs >0.05 good; Histone >0.1 good |
| IDR | Reproducibility (replicate peaks) | idr | <0.05 strong |
| Correlation | Replicate signal correlation | deepTools/R | r > 0.8 good |

---

## Interpretation Guidelines

- **Excellent**: High NSC/RSC, high FRiP, low duplication, consistent replicates (high r, low IDR).  
- **Moderate**: Acceptable NSC/RSC with lower FRiP or moderate IDR.  
- **Poor**: Low RSC/NSC, high duplication, poor reproducibility → consider additional filtering or re-sequencing.

---

## Notes and Best Practices

- Provide input/control when available.  
- Ensure consistent genome builds across samples.  
- Organize all outputs under `qc_results/`.  
- The included report script is minimal and can be customized to parse lab-specific outputs.

---

## References

- phantompeakqualtools: https://github.com/kundajelab/phantompeakqualtools  
- samtools: http://www.htslib.org/  
- idr: https://github.com/nboley/idr  
- picard: https://broadinstitute.github.io/picard/  
- deepTools: https://deeptools.readthedocs.io/