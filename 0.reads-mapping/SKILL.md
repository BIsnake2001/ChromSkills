---
name: reads-mapping
description: Align ChIP-seq or ATAC-seq FASTQ files to a reference genome using Bowtie2, with strict input validation, library layout detection, output organization and logging. Use it when raw sequencing reads must be converted into sorted/indexed BAM files before downstream QC, peak calling, or footprinting.
---

# ChIP-seq / ATAC-seq Sequence Alignment (Bowtie2)

## Overview

This skill performs core **sequence alignment** for **ChIP-seq** and **ATAC-seq** data starting from FASTQ files using **Bowtie2**. It is designed for autonomous execution with explicit user confirmation for biologically important parameters that must not be guessed.

Main steps include:
- Refer to the **Inputs & Outputs** section and create the output architecture in Step 0.
- **Always ask the user** for the reference genome/index to use. **Never infer genome build by filename alone.**
- **Always ask the user** for assay type if it cannot be confidently determined from file names or project context.
- Detect whether input data are **paired-end** or **single-end**.
- Group FASTQ files into samples using file naming conventions.
- Validate that the Bowtie2 index exists before running.
- Run alignment with **Bowtie2**.
- Convert SAM to BAM, then **sort**, **index**, and generate **flagstat** and **idxstats** reports.
- For ATAC-seq, preserve metadata and log that downstream duplicate handling and Tn5 shifting occur in later skills unless explicitly requested here.
- Write a per-sample parameter log file documenting all chosen options and their justification.

---

## Inputs & Outputs

### Inputs

Accepted FASTQ naming patterns include:

```bash
${sample}.fastq.gz
${sample}.fq.gz
${sample}_R1.fastq.gz
${sample}_R2.fastq.gz
${sample}_1.fastq.gz
${sample}_2.fastq.gz
```

Reference resources must be supplied by the user as one of the following:

```bash
/path/to/reference.fa
/path/to/bowtie2_index_prefix
```

### Outputs

```bash
all_alignment/
  aligned_bam/
    ${sample}.sorted.bam
    ${sample}.sorted.bam.bai
  logs/
    ${sample}_alignment.log
    ${sample}_used_parameters.txt
  temp/
```

All outputs must be placed under `${proj_dir}` returned in Step 0.

---

## Required User Inputs

The agent must ask for the following when missing:

- **Reference genome build** (for example hg38, mm10, dm6)
- **Reference asset path**:
  - FASTA path for Bowtie2 index construction, or
  - existing Bowtie2 index prefix
- **Assay type**: ChIP-seq or ATAC-seq, if not obvious from file names or prior context
- **Number of threads** if the user has a compute preference; otherwise use a reproducible default
- Whether to keep intermediate SAM files; default is **no**

The agent must **not** guess:
- genome build
- reference annotation/version
- assay type when ambiguous
- single-end vs paired-end if file pairing is incomplete or inconsistent

---

## Decision Logic

### Assay Type

Use file and sample names only for **tentative classification**:

- Names containing `ATAC`, `OmniATAC`, `scATAC` → likely **ATAC-seq**
- Names containing TF or histone mark identifiers such as `CTCF`, `MYC`, `H3K27ac`, `H3K4me3`, `H3K27me3` → likely **ChIP-seq**

If naming is ambiguous, **ask the user**.

### Library Layout

Use FASTQ grouping rules:

- If files appear as `${sample}_R1` and `${sample}_R2` or `${sample}_1` and `${sample}_2` → paired-end
- If only one FASTQ exists for a sample → single-end
- If an R1 file exists without its mate → stop and ask the user to resolve missing mates

### Threads

- Default to `8` threads unless the user specifies otherwise

---

## Step-by-Step Workflow

### Step 0: Initialize Project

Create a task directory for alignment outputs.

Suggested call:

- `mcp__project-init-tools__project_init`

with:

- `sample`: all
- `task`: alignment
- `genome`: provided by user

The tool will return `${proj_dir}`. Use it for all output placement.

If a project-init MCP tool is not available in the runtime, create this directory structure manually:

```bash
all_alignment/
  aligned_bam/
  qc/
  logs/
  temp/
```

Set `${proj_dir}` to `all_alignment`.

---

### Step 1: Detect and Group FASTQ Files

Call:

- `mcp__bowtie2-tools__detect_fastq_samples`

with:

- `input_dir`: directory containing FASTQ files

The tool will:
- find supported FASTQ files
- group them into samples
- detect single-end vs paired-end layout
- report any missing mate files
- return a machine-readable sample summary

Rules:
- Prefer compressed FASTQ (`*.fastq.gz`, `*.fq.gz`)
- Use consistent sample grouping
- Fail if file naming is inconsistent or ambiguous

---

### Step 2: Validate or Build Bowtie2 Reference Assets

If the user supplied a FASTA, call:

- `mcp__bowtie2-tools__build_bowtie2_index`

with:

- `reference_fasta`: user-provided FASTA path
- `index_prefix`: desired Bowtie2 index prefix

If the user supplied an existing Bowtie2 prefix, validate it before alignment by calling:

- `mcp__bowtie2-tools__validate_bowtie2_index`

with:

- `index_prefix`: user-provided index prefix

Stop and ask the user to correct the path if validation fails.

---

### Step 3: Align Each Sample

For each detected sample, call:

- `mcp__bowtie2-tools__run_bowtie2_alignment`

with:

- `sample_name`: sample identifier
- `fastq_r1`: path to R1 FASTQ or single-end FASTQ
- `fastq_r2`: path to R2 FASTQ for paired-end data, otherwise omit
- `assay_type`: `chipseq` or `atacseq`
- `index_prefix`: validated Bowtie2 index prefix
- `out_dir`: `${proj_dir}/aligned_bam`
- `log_dir`: `${proj_dir}/logs`
- `threads`: user-specified or default `8`
- `keep_sam`: `false` by default

Tool behavior:
- run Bowtie2
- write an alignment log
- convert to BAM
- sort BAM
- index BAM
- remove intermediate SAM unless `keep_sam=true`

Expected output:

```bash
${proj_dir}/aligned_bam/${sample}.sorted.bam
${proj_dir}/aligned_bam/${sample}.sorted.bam.bai
${proj_dir}/logs/${sample}_alignment.log
```

Notes:
- For ATAC-seq, this skill performs alignment only. Duplicate marking/removal, mitochondrial filtering, proper-pair filtering, and Tn5 shifting belong to downstream preprocessing or peak-calling skills unless explicitly requested in another skill.
- Do not remove duplicates here unless the user specifically requests an alternate alignment workflow.

---

### Step 4: Write Parameter Log

For each sample, the agent must write:

```bash
${proj_dir}/logs/${sample}_used_parameters.txt
```

Example content:

```text
Sample: ATAC_rep1
Assay type: ATAC-seq
Library layout: paired-end
Aligner: bowtie2
Reference genome build: hg38
Reference index: /refs/hg38/bowtie2/hg38
Threads: 8
Intermediate SAM kept: no

Reasoning:
- Sample name contains ATAC, so assay classified as ATAC-seq
- Paired FASTQ mates were detected automatically
- User provided hg38 Bowtie2 index
- Alignment-only workflow selected; duplicate handling and Tn5 shifting deferred to downstream preprocessing/peak-calling
```

---

## Failure Handling

Stop execution and ask the user for correction if any of the following occurs:

- no FASTQ files found
- inconsistent file naming prevents sample grouping
- paired-end mate missing
- reference FASTA or index prefix missing
- `bowtie2` or `bowtie2-build` executable not found in PATH
- `samtools` not found in PATH
- alignment command returns non-zero exit status
- sorted BAM or BAM index is not created

Do not continue to downstream QC if alignment fails.

---

## Exact MCP Tool Calls Required

1. `mcp__bowtie2-tools__detect_fastq_samples`
2. `mcp__bowtie2-tools__build_bowtie2_index` (only when user provides FASTA)
3. `mcp__bowtie2-tools__validate_bowtie2_index`
4. `mcp__bowtie2-tools__run_bowtie2_alignment`

---

## When the Agent Must Ask the User

The agent must ask before execution when any of the following are missing or ambiguous:

- genome build
- reference FASTA or Bowtie2 index prefix
- assay type when not inferable from context
- whether incomplete FASTQ pairs should be excluded or fixed
- thread count if a project-specific compute policy exists

The agent must not invent these values.
