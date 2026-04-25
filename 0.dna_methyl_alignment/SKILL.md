---
name: dna-methylation-alignment-bismark
description: Align bisulfite sequencing DNA methylation reads using Bismark only, with explicit validation of reference preparation, library layout detection, output organization, logging, and alignment QC. Use it for WGBS, RRBS, or other bisulfite-converted DNA methylation sequencing data when raw FASTQ files must be aligned before methylation extraction and downstream analysis.
---

# DNA Methylation Sequence Alignment with Bismark

## Overview

This skill performs **bisulfite-aware sequence alignment** for **DNA methylation sequencing** using **Bismark only**. It is designed for autonomous execution from FASTQ input through aligned BAM generation and basic QC, while preventing unsafe assumptions about genome build, library layout, or assay design.

Main steps include:
- Refer to the **Inputs & Outputs** section and create the output architecture in Step 0.
- **Always ask the user** for the reference genome directory to use for Bismark. **Never infer genome build from filenames alone.**
- **Always ask the user** whether the assay is **WGBS**, **RRBS**, or another bisulfite-based methylation assay if that affects trimming or downstream interpretation and is not already known.
- Detect whether input data are **paired-end** or **single-end**.
- Group FASTQ files into samples using naming conventions.
- Validate that the Bismark genome folder has been prepared.
- Run **Bismark** alignment.
- Sort and index BAM output with **samtools**.
- Generate basic alignment QC reports and write a parameter log file for every sample.
- Keep the workflow alignment-focused; methylation extraction belongs in a downstream skill.

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
/path/to/reference_genome_folder
/path/to/reference.fa
```

If the user supplies a FASTA, prepare a Bismark genome folder before alignment.

### Outputs

```bash
all_methylation_alignment/
  aligned_bam/
    ${sample}.sorted.bam
    ${sample}.sorted.bam.bai
  qc/
    ${sample}.flagstat.txt
    ${sample}.idxstats.txt
  logs/
    ${sample}_alignment.log
    ${sample}_used_parameters.txt
  temp/
```

All outputs must be placed under `${proj_dir}` returned in Step 0.

---

## Required User Inputs

The agent must ask for the following when missing:

- **Reference genome build**
- **Bismark genome folder** or a **reference FASTA** to prepare one
- **Assay type** when biologically relevant and not already known:
  - WGBS
  - RRBS
  - other bisulfite-based methylation assay
- **Number of threads** if the user has a compute preference; otherwise use a reproducible default
- Whether to keep intermediate BAM files; default is **no**

The agent must **not** guess:
- genome build
- assay type when ambiguous
- single-end vs paired-end when file pairing is incomplete or inconsistent
- whether a genome folder prepared for another build/version is acceptable

---

## Decision Logic

### Assay Type

Use file and sample names only for **tentative classification**:

- Names containing `WGBS` → likely whole-genome bisulfite sequencing
- Names containing `RRBS` → likely reduced representation bisulfite sequencing
- Names containing only generic terms such as `methylation`, `BSseq`, `bisulfite` → ambiguous

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

Create a task directory for methylation alignment outputs.

Suggested call:

- `mcp__project-init-tools__project_init`

with:

- `sample`: all
- `task`: methylation_alignment
- `genome`: provided by user

The tool will return `${proj_dir}`. Use it for all output placement.

If a project-init MCP tool is not available in the runtime, create this directory structure manually:

```bash
all_methylation_alignment/
  aligned_bam/
  logs/
  temp/
```

Set `${proj_dir}` to `all_methylation_alignment`.

---

### Step 1: Detect and Group FASTQ Files

Call:

- `mcp__bismark-tools__detect_fastq_samples`

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

### Step 2: Prepare or Validate the Bismark Genome Folder

If the user supplied a FASTA and wants the agent to prepare a Bismark genome folder, call:

- `mcp__bismark-tools__prepare_bismark_genome`

with:

- `reference_fasta`: user-provided FASTA path
- `genome_folder`: destination directory for Bismark genome preparation

If the user supplied an existing Bismark genome folder, validate it before alignment by calling:

- `mcp__bismark-tools__validate_bismark_genome`

with:

- `genome_folder`: user-provided genome folder

Stop and ask the user to correct the path if validation fails.

---

### Step 3: Align Each Sample with Bismark

For each detected sample, call:

- `mcp__bismark-tools__run_bismark_alignment`

with:

- `sample_name`: sample identifier
- `fastq_r1`: path to R1 FASTQ or single-end FASTQ
- `fastq_r2`: path to R2 FASTQ for paired-end data, otherwise omit
- `genome_folder`: validated Bismark genome folder
- `out_dir`: `${proj_dir}/aligned_bam`
- `log_dir`: `${proj_dir}/logs`
- `temp_dir`: `${proj_dir}/temp`
- `threads`: user-specified or default `8`
- `keep_intermediate_bam`: `false` by default

Tool behavior:
- run Bismark only
- write an alignment log
- locate the Bismark BAM output
- sort BAM
- index BAM
- remove intermediate unsorted BAM unless `keep_intermediate_bam=true`

Expected output:

```bash
${proj_dir}/aligned_bam/${sample}.sorted.bam
${proj_dir}/aligned_bam/${sample}.sorted.bam.bai
${proj_dir}/logs/${sample}_alignment.log
```

This skill is alignment-only. Methylation extraction, deduplication policy, and cytosine report generation belong to downstream skills unless explicitly requested elsewhere.

---



### Step 4: Write Parameter Log

For each sample, the agent must write:

```bash
${proj_dir}/logs/${sample}_used_parameters.txt
```

Example content:

```text
Sample: WGBS_rep1
Assay type: WGBS
Library layout: paired-end
Aligner: Bismark
Reference genome build: hg38
Bismark genome folder: /refs/hg38/bismark_genome
Threads: 8
Intermediate BAM kept: no

Reasoning:
- Sample name contains WGBS, so assay classified as whole-genome bisulfite sequencing
- Paired FASTQ mates were detected automatically
- User provided hg38 Bismark genome folder
- Alignment-only workflow selected; methylation extraction deferred to downstream analysis
```

---

## Failure Handling

Stop execution and ask the user for correction if any of the following occurs:

- no FASTQ files found
- inconsistent file naming prevents sample grouping
- paired-end mate missing
- Bismark genome folder missing or not prepared
- requested executable not found in PATH
- Bismark returns non-zero exit status
- sorted BAM or BAM index is not created

Do not continue to downstream QC if alignment fails.

---

## Exact MCP Tool Calls Required

1. `mcp__bismark-tools__detect_fastq_samples`
2. `mcp__bismark-tools__prepare_bismark_genome` (only when user provides FASTA)
3. `mcp__bismark-tools__validate_bismark_genome`
4. `mcp__bismark-tools__run_bismark_alignment`

---

## When the Agent Must Ask the User

The agent must ask before execution when any of the following are missing or ambiguous:

- genome build
- Bismark genome folder or reference FASTA
- assay type when not inferable from context
- whether incomplete FASTQ pairs should be excluded or fixed
- thread count if a project-specific compute policy exists

The agent must not invent these values.
