---
name: hic-loop-calling
description: This skill performs chromatin loop detection from Hi-C .mcool files using cooltools.
---

# Hi-C Loop Calling

## Overview

This skill provides a minimal and efficient workflow for detecting chromatin loops from Hi-C data stored in .mcool format and preparing results for visualization in IGV. The key steps involved include:
- Refer to the **Inputs & Outputs** section to verify required files and output structure.
- **Always prompt user** for genome assembly used.
- **Always prompt user** for resolution used to call cloops.
- **Rename chromosomes** in the .mcool or .cool file to satisfy the chromosome format with "chr".
- Generate chromosome-arm view files for compartment calling **after changing the chromosome name**.
- **Extract contact matrices** from .mcool files at the desired resolution.
- **Detect chromatin loops** using `cooltools`.

---

## When to Use This Skill

Use this skill when:

- You need to identify chromatin loops from Hi-C data in .mcool format.

--

## Inputs & Outputs

### Inputs

- **File format:** .mcool (Hi-C data file).
- **Resolution:** Choose the desired resolution for loop calling (e.g., 5 kb, 10 kb, etc.).

### Outputs

```bash
loop_calling/
    loops/
        ${sample}_loops_${res}.bedpe  # Detected chromatin loops in BEDPE format.
    temp/
        view_${genome}.tsv
        ${sample}.expected.cis.${res}.tsv 
```
---

## Decision Tree

### Step 1: Extract Contact Matrix at a Given Resolution

Select the desired resolution (e.g., 5 kb) and extract it using `cooler`. First, check available resolutions in the .mcool file:

```bash
cooler ls <input.mcool>
```

### Step 2: Modify the chromosome name in .mcool file

**Python Example:**
```python
import cooler
clr = cooler.Cooler(f'input.mcool::/resolutions/{resoluton}')
rename_dict = {chrom:f"chr{chrom}" for chrom in clr.chromnames if not chrom.startswith('chr')}
cooler.rename_chroms(clr, rename_dict)
```

### Step 3: Create Chromosome-Arm View

Generate a view file defining chromosome arms (based on centromere positions). Change chromosome name in the .mcool or .cool file if not consistent with those in the .fa file using `cooler.rename_chroms`.

**Python Example:**
```python
import bioframe, cooler

hg38_chromsizes = bioframe.fetch_chromsizes('hg38') # change the genome name according to user feedback
hg38_cens = bioframe.fetch_centromeres('hg38')
view_hg38 = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

clr = cooler.Cooler(f'input.mcool::/resolutions/{res}')
view_hg38 = view_hg38[view_hg38.chrom.isin(clr.chromnames)]
view_hg38.to_csv('view_hg38.tsv', sep='	', header=False, index=False)
```

### Step 4: Detect Chromatin Loops with `cooltool`

```bash
cooltools expected-cis --nproc 6 -o <sample>.expected.cis.<res>.tsv --view "view_hg38.tsv" <sample>.mcool::resolutions/<res>
cooltools dots --nproc 6 -o <sample>_loops_<res>.bedpe --view view_hg38.tsv <sample>.mcool::resolutions/<res> <sample>.expected.cis.<res>.tsv
# Remove header (first line) from the bedpe file, if present
sed -i '1{/^chrom1\tstart1\tend1\tchrom2\tstart2\tend2.*/d}' <sample>_loops_<res>.bedpe
```

**Notes:**
- The output file `<sample>_loops_<res>.bedpe` contains chromatin loop anchors and their scores.

---

## Best Practices

- Choose a resolution that is appropriate for your sequencing depth (e.g., 5–100 kb).
