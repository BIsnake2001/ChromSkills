---
name: nested-TAD-detection
description: This skill detects hierarchical (nested) TAD structures from Hi-C contact maps (in .cool or mcool format) using OnTAD, starting from multi-resolution .mcool files. It extracts a user-specified chromosome and resolution, converts the data to a dense matrix, runs OnTAD, and organizes TAD calls and logs for downstream 3D genome analysis.
---

# Nested TAD Detection from .mcool Using OnTAD

## Overview

This skill performs nested TAD (hierarchical TAD/subTAD) detection from Hi-C data using **OnTAD**, starting from a .mcool or ..cool file.

Main steps include:

- Refer to the **Inputs & Outputs** section to verify required files and output structure.
- Inspect the `.mcool` file to list available resolutions and alway remember to confirm the chromosome name and analysis resolution with the user.
- Extract a **balanced or raw dense Hi-C matrix** for a selected chromosome and resolution from the `.mcool` file.
- Ensure matrix quality (symmetry, no all-zero rows/columns, reasonable contact decay).
- Run **OnTAD** to call TADs and nested TAD structures.
- Parse and standardize OnTAD output into BED-like tables and hierarchical annotation files.

---

## When to use this skill

Use this skill when you want to **identify TADs and nested sub-TADs** from high- or mid-resolution Hi-C data, especially when your contact maps are stored as **Cooler multi-resolution files (.mcool)** and you need **chromosome- and resolution-specific** OnTAD calls.

Typical biological questions / use-cases:

- Comparing **TAD hierarchy** between cell types (e.g., GM12878 vs K562) or conditions (control vs treated).
- Investigating whether **inner subTADs** are enriched for active regulatory elements, specific histone marks, or gene expression.
- Studying **boundary usage**, **boundary sharing**, or **hierarchical TAD levels** around key loci (e.g., HOX clusters, oncogenes).
- Integrating nested TAD structure with **ChIP-seq**, **ATAC-seq**, or **WGBS** to understand spatial regulatory architecture.

Data quality & replication assumptions:

- Hi-C experiments should have **sufficient depth** for the target resolution (e.g., ≥5–10 kb typically requires deep sequencing).
- Preferably, use **biological replicates** per condition and call TADs on either:
  - individual replicate matrices (and then merge/consensus), or  
  - replicate-merged matrices (if justified).
- The `.mcool` should be **properly normalized or at least QC’d** (ICE/balanced weights available if using `--balanced`).

---

## Inputs & Outputs

### Inputs

Required core inputs:

- **Hi-C matrix file**
  - Multi-resolution Cooler file (`.mcool`), e.g.:
    - `sample.mcool`
- **User-supplied parameters (must come from user feedback)**
  - Chromosome name: e.g., `chr1`, `chr2`, `chrX`
  - Resolution of interest: e.g., `10000`, `25000`, `40000` (in bp)
  - Chromosome length: e.g., 133275309
- **Software and environment**
  - `cooler` command-line utilities (or Python `cooler` module) available in `$PATH` or environment.
  - `OnTAD` installed and callable (e.g., `OnTAD` or `OnTAD_linux` command).
  - `python` (3.x) and basic scientific Python stack if using Python-based extraction.

Optional entry point:

- **Precomputed dense Hi-C matrix (OnTAD-ready)**
  - Plain text dense square matrix (no header), e.g. `chr1_10kb_dense.matrix`.  
  - If this is already present, you may **skip the `.mcool` → dense matrix conversion** and jump directly to OnTAD.

Operational rules for missing inputs:

- If `.mcool` file is missing:  
  `"sample.mcool not available, provide required files or skip and proceed ?"`
- If chromosome list not specified:  
  Ask the user explicitly rather than assuming default.
- If OnTAD executable is not found:  
  Ask user to install/locate OnTAD before proceeding.

---

### Outputs

Default output directory structure:

```bash
nested_TAD_detection/
    scripts/ # scripts generated in this section
    config/
        params.yaml
    matrices/
        ${chromosome}_${res}_dense.txt
        ${chromosome}_${res}_dense.log
    Nested_TADs/
        ${chromosome}_${res}.tad
        ${chromosome}_${res}.bed
```
---

## Decision Tree


### Step 1: Check required inputs, list resolutions and check chromosome length
```bash
MCool=sample.mcool
cooler ls "$MCool" # List groups/resolutions inside .mcool
RES=25000  # chosen by user
COOL="${MCool}::/resolutions/${RES}"

cooler dump -t chroms sample.mcool::/resolutions/${RES} # check chromosome length
```

### Step 2: Check environment, inputs, params

```bash

which OnTAD || echo "OnTAD executable not found."
mkdir -p nested_TAD_OnTAD/{config,matrices,OnTAD_calls}
```

And record parameters in YAML.

---

### Step 3: Extract dense matrix from `.mcool`

- **example code**

```python
import cooler, numpy as np
c = cooler.Cooler(f"sample.mcool::/resolutions/{resolution}")
mat = c.matrix(balance=True).fetch(chrom)
mat = np.nan_to_num(mat, nan=0.0)
np.savetxt("${chromosome}_${res}_dense.matrix", mat,
           fmt="%.6f", delimiter="\t")
```

---

### Step 4: Run OnTAD

```bash
OnTAD ${chromosome}_${res}_dense.matrix -penalty 0.1 -maxsz 200 -o ${chromosome}_${res}_OnTAD.tad -bedout <chr> <chrlength> <resolution> > OnTAD_${chromosome}_${res}.log 2>&1
```

## Advanced Usage

detailed usage:
```bash
OnTAD <Hi-C matrix> [-penalty <float>] [-maxsz <int>] [-minsz <int>] [-ldiff <float>] [-lsize <int>] [-bedout <chrnum int> <chrlength int> <int>] [-log2] [-o output_file] [-hic_norm <NONE/VC/VC_SQRT/KR>]
```

