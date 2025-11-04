---
name: alignment-QC
description: This skill performs alignment-level quality control for ChIP-seq or ATAC-seq BAM files using samtools, Picard, and MultiQC. It checks (and only if needed, fixes) BAM sorting and read group information, collects QC metrics, and organizes outputs so that QC artifacts go to qc_results/, required downstream files stay in the current directory, and transient files go to temp/.
---

# Alignment Quality Control for ChIP-seq/ATAC-seq

## Overview

Perform comprehensive **preliminary alignment-level quality control** for ChIP-seq and ATAC-seq BAM files using **samtools**, **Picard**, and **MultiQC**.  
This skill:
- **Does not re-sort or add read groups if they already exist** (fully conditional).
- **Run preliminary QC metrics**
- **Generate MultiQC report**  
- **Organizes outputs into**:
  - **qc_results/** → all QC outputs and reports  
  - **temp/** → intermediate temporary files created only when fixes are required  
  - **current directory** → essential files for downstream steps (final working BAMs and their indexes)

---

## Workflow Decision Tree

Follow the workflow below for each `*.bam` in the current directory.

1. **Check sorting and read group information (conditional fixes only)**  
   Create output folders once:
   ```bash
   mkdir -p qc_results temp
   ```

   Set the input and working BAM for the sample (replace `sample.bam` in a loop as needed):
   ```bash
   in_bam="sample.bam"
   work_bam="$in_bam"
   sample_basename="$(basename "$in_bam" .bam)"
   ```

   - **Sorting check** (no re-sorting if already coordinate-sorted):
     ```bash
     if samtools view -H "$work_bam" | grep -q '^@HD.*SO:coordinate'; then
       echo "[OK] $work_bam already coordinate-sorted"
     else
       echo "[FIX] Sorting $work_bam"
       samtools sort -@ 8 -o "temp/${sample_basename}.sorted.bam" "$work_bam"
       # keep the sorted BAM in current directory for downstream steps
       mv "temp/${sample_basename}.sorted.bam" "${sample_basename}.sorted.bam"
       work_bam="${sample_basename}.sorted.bam"
     fi
     ```

   - **Read group (@RG) check** (no modification if @RG exists):
     ```bash
     if samtools view -H "$work_bam" | grep -q '^@RG'; then
       echo "[OK] Read group present in $work_bam"
     else
       echo "[FIX] Adding read group to $work_bam"
       picard AddOrReplaceReadGroups          I="$work_bam"          O="temp/${sample_basename}.RG.bam"          RGID="${sample_basename}" RGLB="lib1" RGPL="illumina" RGPU="unit1" RGSM="${sample_basename}"
       # keep the RG-annotated BAM in current directory for downstream steps
       mv "temp/${sample_basename}.RG.bam" "${sample_basename}.RG.bam"
       work_bam="${sample_basename}.RG.bam"
     fi
     ```

2. **Index BAM files (always for the final working BAM in current directory)**  
   The index is needed downstream; keep it next to the final working BAM:
   ```bash
   samtools index "$work_bam"
   ```

3. **Run preliminary QC metrics**  
   All QC outputs go to **qc_results/**. No BAM content is changed here.
   ```bash
   # Basic alignment summaries
   samtools flagstat "$work_bam" > "qc_results/${sample_basename}.flagstat.txt"
   samtools stats    "$work_bam" > "qc_results/${sample_basename}.stats.txt"

   picard CollectInsertSizeMetrics      I="$work_bam" O="qc_results/${sample_basename}.insertsize_metrics.txt"      H="qc_results/${sample_basename}.insertsize_histogram.pdf" M=0.5

   picard MarkDuplicates      I="$work_bam" O="temp/${sample_basename}.markdup.bam"      M="qc_results/${sample_basename}.dup_metrics.txt" REMOVE_DUPLICATES=false
   ```

4. **Generate MultiQC report**  
   Aggregate all QC outputs from `qc_results/` into a single HTML report saved under `qc_results/`:
   ```bash
   multiqc qc_results/ --filename qc_results/alignment_qc_report.html
   ```

5. **Evaluate results**  
   Assess mapping rate, duplication rate, insert size, and mitochondrial content using thresholds in `references/qc_metrics.md`.

6. **Generate summary report**  
   Summarize each BAM’s PASS/WARN/FAIL against thresholds and write recommendations to:
   ```
   qc_results/qc_summary.txt
   ```

---

## Directory Structure

After running this skill, the layout will be:

```
project_directory/
├── sample.bam                        # original input (unchanged)
├── sample.sorted.bam                 # only present if sorting was required
├── sample.sorted.bam.bai             # index for the current working BAM
├── sample.RG.bam                     # only present if RG was added
├── sample.RG.bam.bai                 # index for the current working BAM
├── qc_results/
│   ├── sample.flagstat.txt
│   ├── sample.stats.txt
│   ├── sample.alignment_metrics.txt
│   ├── sample.insertsize_metrics.txt
│   ├── sample.dup_metrics.txt
│   ├── sample.insertsize_histogram.pdf
│   └── alignment_qc_report.html
└── temp/
    ├── sample.markdup.bam            # intermediate for duplication estimation only
    └── intermediate_files_if_any...
```

> ⚠️ **Conditional behavior:** If the input BAM is already coordinate-sorted and contains @RG, no new sorted/RG BAMs are created; `work_bam` remains the original input, only its index (`.bai`) is produced (if missing). Temporary files under `temp/` are created **only** when a fix or a QC intermediate is required.

---

## Quick Start

To perform alignment QC on all BAM files in the current directory:

1. Ensure **samtools**, **Picard**, and **MultiQC** are installed and available in PATH.

2. Run the basic QC workflow. The skill will:
   - Verify BAM sorting and read group information if necessary
   - Compute QC metrics
   - Save results into target directory as required

3. Run MultiQC to generate the integrated report:
  ```bash
  # Aggregate report → qc_results/alignment_qc_report.html
  multiqc qc_results/ --filename qc_results/alignment_qc_report.html
  ```

---

## Quality Assessment

### Key QC Metrics

- **Total reads** – overall sequencing depth  
- **Mapped reads (%)** – alignment efficiency  
- **Properly paired (%)** – valid pair fraction (paired-end)  
- **Duplicate rate (%)** – PCR duplication estimate  
- **Mitochondrial reads (%)** – mitochondrial contamination  
- **Insert size distribution** – fragment length profile  

All metrics are derived from `samtools`/`Picard` and summarized by MultiQC.

---

### Quality Thresholds

| Category | Criteria | Interpretation |
|-----------|-----------|----------------|
| **Pass** | All metrics within recommended thresholds | Suitable for downstream analysis |
| **Warn** | One or more borderline metrics | Likely acceptable; review recommended |
| **Fail** | Critical metrics outside acceptable ranges | Re-sequencing or reprocessing suggested |

---

## MultiQC Integration

All QC outputs in `qc_results/` are aggregated into a **single HTML report**:
```bash
multiqc qc_results/ --filename qc_results/alignment_qc_report.html
```

---

## Report Generation

After MultiQC completes, generate a sample-wise summary (PASS/WARN/FAIL) per thresholds in `references/qc_metrics.md` and save it as:
```
qc_results/qc_summary.txt
```

---

## Resources

Use `references/qc_metrics.md` for:
- Metric definitions and recommended thresholds  
- Troubleshooting guidance  
- Readiness criteria for peak calling  
- Pointers to ENCODE/nf-core QC standards
