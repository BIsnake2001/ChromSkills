---
name: peak-calling
description: Automatically perform peak calling for ChIP-seq or ATAC-seq data using MACS2, with intelligent genome and parameter detection.
---

# MACS2 Peak Calling Auto Skill

## Overview

This skill automatically performs core **peak calling** with **MACS2** for **ChIP-seq** and **ATAC-seq** data, based on the BAM files in the current directory.  
It includes automatic experiment recognition, genome detection, and parameter selection.

### Steps
1. **Identify BAM files** and classify treatment/control samples. Always use filtered BAM file (`filtered.bam`) if availabel  
2. **Detect sequencing type** (single-end or paired-end) using SAM/BAM flags.  
3. **Determine genome assembly** automatically from BAM header.  
4. **Detect experiment type** (TF, histone mark, or ATAC-seq).  
5. **Automatically decide** whether to call narrow or broad peaks, whether control is needed, and which parameters to use.  
6. **Generate a parameter log file** (`used_parameters.txt`) with justification for each chosen option.  
7. **Perform MACS2 peak calling** accordingly and save peaks in a dedicated folder.

---

## Decision Tree

### Step 1. Identify and Classify BAM Files

**Command Example**
```bash
find . -name "*.bam" | sort
```

- **Treatment BAMs**: filenames contain TF names or histone marks (e.g., `CTCF`, `H3K27me3`, `ATAC`).
- **Control BAMs**: filenames contain “input”, “control”, or “IgG”.
- Prefer **deduplicated** BAMs (`dedup_*.bam`, `rmdup_*.bam`) if available.

---

### Step 2. Detect Sequencing Type (Single-End or Paired-End)

**Command Example**
```bash
samtools flagstat sample.bam | egrep "properly paired|singletons"
```

- If `properly paired > 0` → Paired-end (`-f BAMPE`)
- If `singletons ≈ total` → Single-end (`-f BAM`)

---

### Step 3. Determine Genome Assembly Automatically

**Command Example**
```bash
samtools view -H sample.bam | grep "@SQ" | head
```

Typical recognition rules:
| Pattern in header | Genome detected | MACS2 parameter |
|--------------------|------------------|-----------------|
| `chr1`, `chrX` | Human (hg19/hg38) | `-g hs` |
| `1`, `X` | Mouse (mm9/mm10) | `-g mm` |
| `chrIII` | Yeast | `-g 1.2e7` |

If unclear, ask user:
> “Unable to determine genome assembly from BAM headers. Please specify (hg19/hg38/mm10/etc.):”

---

### Step 4. Detect Experiment Type and Choose Peak Mode

| Detected Pattern | Experiment Type | Peak Type | Parameter Key Options |
|------------------|-----------------|------------|------------------------|
| TF name (CTCF, GATA1, MYC, TP53…) | TF ChIP-seq | Narrow | `--call-summits -q 0.01` |
| Active histone marks (H3K4me3, H3K27ac, H3K9ac) | Histone (sharp) | Narrow | `--call-summits -q 0.05` |
| Broad histone marks (H3K27me3, H3K9me3, H3K36me3) | Histone (broad) | Broad | `--broad --broad-cutoff 0.1 -q 0.05` |
| H3K4me1 | Intermediate | Narrow | `--call-summits -q 0.05` (optional `--broad`) |
| ATAC | ATAC-seq | Narrow | `--nomodel --shift -100 --extsize 200 -q 0.05` |

If experiment type cannot be inferred, ask:
> “Could not determine experiment type. Is it TF / histone modification / ATAC-seq?”

---

### Step 5. Determine Whether Control File is Needed

- **TF or histone ChIP-seq**: control BAM required if available.  
  Example:
  ```bash
  macs2 callpeak -t TF_treatment.bam -c TF_input.bam ...
  ```
- **ATAC-seq**: no control used.

If no control found for ChIP-seq:
> “No control file detected. Continue without control using relaxed threshold?”

---

### Step 6. Execute MACS2 with Auto Parameters

Example conditional execution (pseudo-script):

```bash
if [ $EXPERIMENT_TYPE == "TF" ]; then
  macs2 callpeak -t treatment.bam -c control.bam -f $FORMAT -g $GENOME -n TF_peaks   --outdir peaks/ --call-summits -q 0.01

elif [ $EXPERIMENT_TYPE == "BROAD" ]; then
  macs2 callpeak -t treatment.bam -c control.bam -f $FORMAT -g $GENOME -n broad_peaks   --outdir peaks/ --broad --broad-cutoff 0.1 -q 0.05

elif [ $EXPERIMENT_TYPE == "ATAC" ]; then
  macs2 callpeak -t atac.bam -f $FORMAT -g $GENOME -n ATAC_peaks   --outdir peaks/ --nomodel --shift -100 --extsize 200 -q 0.05
fi
```

---

### Step 7. Generate Parameter Log File

After auto-selection, the skill writes a log file:
```
peaks/used_params_peak_calling.txt
```

**Example content:**
```
Genome detected: hg38
Experiment type: H3K27me3 (broad histone)
Sequencing type: paired-end
Control used: input_control.bam
MACS2 mode: --broad --broad-cutoff 0.1 -q 0.05
Reasoning:
- Broad mark (H3K27me3) requires domain-level detection
- Control detected and applied
- Genome identified as human; using -g hs
- Paired-end library; use -f BAMPE
```

If any parameter is unknown, the script prompts the user for clarification before execution.

---

## Example Workflows

**CTCF ChIP-seq with Control (TF narrow peaks)**
```bash
macs2 callpeak -t CTCF_treatment.bam -c input.bam -f BAMPE -g hs -q 0.01 --call-summits -n CTCF
```

**H3K27me3 (broad peaks)**
```bash
macs2 callpeak -t H3K27me3.bam -c input.bam -f BAMPE -g hs --broad --broad-cutoff 0.1 -n H3K27me3
```

**ATAC-seq**
```bash
macs2 callpeak -t ATAC.bam -f BAMPE -g hs --nomodel --shift -100 --extsize 200 -q 0.05 -n ATAC
```

---

## Best Practices

1. Always use deduplicated BAMs (`rmdup_` or `dedup_`).
2. Provide input controls when available.
3. Use `--call-summits` for TF and ATAC-seq.
4. Keep broad domains for H3K27me3/H3K9me3 using `--broad`.
5. Document all parameters in `used_parameters.txt`.
6. Visually check peak quality in IGV or UCSC Browser.

---

## Output Files

| File | Description |
|------|--------------|
| `*_peaks.narrowPeak` / `*_peaks.broadPeak` | Peak coordinates |
| `*_summits.bed` | Peak summits (narrow mode only) |
| `*_peaks.xls` | MACS2 statistics |
| `used_parameters.txt` | Summary of chosen genome, type, mode, and reasoning |
| `peaks/` | Output directory for all peak files |

---

## Resources
Use `references/used_params_peak_calling` as a temptate for generating log file