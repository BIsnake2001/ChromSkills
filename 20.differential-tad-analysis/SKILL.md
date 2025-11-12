---
name: differential-tad-analysis
description: This skill performs differential topologically associating domain (TAD) analysis using HiCExplorer's hicDifferentialTAD tool. It compares Hi-C contact matrices between two conditions based on existing TAD definitions to identify significantly altered chromatin domains.
---

# Differential TAD Analysis with HiCExplorer

## Overview

This skill identifies differentially interacting TADs between two experimental conditions using HiCExplorer.  
It assumes that TADs have already been called for the target condition.

Steps:
1. Normalize Hi-C matrices between conditions using `hicNormalize`. Modify the chromosome name in the .mcool file if not started with "chr".
2. Prepare TAD domains (BED file) obtained from `hicFindTADs` of the target sample. Make sure the consistence of the chromosame names between .mcool files and BED files. Modify the chromosome name in the BED file if not consistent with the .mcool file.
3. Perform differential TAD analysis using `hicDifferentialTAD`.
4. Visualize and interpret significant TAD changes.

Example commands to harmonize chromosome naming

```python
import cooler
clr = cooler.Cooler('input.mcool::/resolutions/100000')

rename_dict = {chrom:f"chr{chrom}" for chrom in clr.chromnames if not chrom.startswith('chr')}
cooler.rename_chroms(clr, rename_dict)
```


## When to Use This Skill

Use this skill when:
- You have already called TADs for one condition.
- You want to detect TADs that show significant interaction changes between two Hi-C matrices.
- You are comparing chromatin architecture between experimental conditions (e.g., treated vs. control, different cell types).

---

## Decision Tree


### Step 1. Normalize Hi-C matrices

To ensure both matrices have comparable sequencing depth and coverage, perform normalization before analysis.

Example:
```bash
hicNormalize --matrices target.mcool::/resolutions/25000 control.mcool::/resolutions/25000 --normalize smallest -o target_norm.cool control_norm.cool
```

**Notes:**
- Use the same resolution for both matrices.
- `--normalize smallest` scales both matrices to the same total number of contacts.

---

### Step 2. Prepare TAD domains

You must provide a BED file of TAD domains obtained from a prior `hicFindTADs` analysis.  
Only the TAD file from the **target** condition is required.

Example format:
```
chr1    1000000    2000000    TAD_1    0    .
chr1    2000000    3200000    TAD_2    0    .
```

---

### Step 3. Run differential TAD analysis

Use `hicDifferentialTAD` to detect TADs with statistically different intra- and inter-TAD interactions between the normalized matrices.

Example command:
```bash
hicDifferentialTAD   -tm target_norm.cool   -cm control_norm.cool   -td target_tads.bed   -o diffTAD_results   -p 0.05   -m all   -mr one   -t 4
```

**Parameter Notes:**
- `-tm` : Target Hi-C matrix.
- `-cm` : Control Hi-C matrix.
- `-td` : TAD domains from the target sample.
- `-p`  : Significance cutoff (default 0.05).
- `-m`  : Region types used for testing (`intra-TAD`, `left-inter-TAD`, `right-inter-TAD`, or `all`).
- `-mr` : Reject rule: `one` (default) or `all` for stricter criteria.
- `-t`  : Number of threads for chromosome-level parallelization.

---

### Step 4. Visualization

Visualize the Hi-C matrices and the identified differential TADs to better interpret the results.

#### 4.1 Visualize Hi-C contact maps
Use `hicPlotMatrix` to display the normalized Hi-C contact matrices and overlay TAD domains.

Example:
```bash
hicPlotMatrix   --matrix target_norm.cool::/resolutions/25000   --region chr1:1,000,000-5,000,000   --perChr   --log1p   --dpi 300   --outFileName target_contactmap_chr1.png   --title "Target condition contact map"
```

Do the same for the control matrix:
```bash
hicPlotMatrix   --matrix control_norm.cool::/resolutions/25000   --region chr1:1,000,000-5,000,000   --perChr   --log1p   --dpi 300   --outFileName control_contactmap_chr1.png   --title "Control condition contact map"
```

#### 4.2 Overlay differential TADs
Plot the differential TADs identified by `hicDifferentialTAD` using `hicPlotTADs`.

Example:
```bash
hicPlotTADs --tracks target_tads.bed diffTAD_results_accepted.diff_tad --region chr1:1,000,000-5,000,000 --outFileName diffTAD_visualization_chr1.png   --dpi 300   --title "Differential TADs (accepted regions)"
```

**Notes:**
- The `diffTAD_results_accepted.diff_tad` file highlights TADs with significant differential interactions.
- Compare accepted and rejected TADs to assess the degree of chromatin reorganization.

---

## Output Summary

| File | Description |
|------|--------------|
| `*_accepted.diff_tad` | TADs showing significant differential interactions |
| `*_rejected.diff_tad` | Non-significant TADs |
| `log.txt` | Analysis log with per-chromosome statistics |
| `*_contactmap_chrN.png` | Normalized Hi-C contact maps |
| `*_diffTAD_visualization_chrN.png` | Visualization of differential TADs |

---

## Quality Control

- Ensure both Hi-C matrices are balanced and have identical resolutions.
- Confirm that the TAD BED file is derived from the target condition at the same resolution.
- Verify consistent normalization methods across all samples.
- Inspect the resulting visualizations to confirm biological plausibility.
