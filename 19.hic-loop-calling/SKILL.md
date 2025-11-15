---
name: hic-loop-calling
description: This skill performs chromatin loop detection from Hi-C .mcool files using HiCExplorer's hicDetectLoops and generates IGV-ready interaction tracks for visualization.
---

# Hi-C Loop Calling

## Overview

This skill provides a minimal and efficient workflow for detecting chromatin loops from Hi-C data stored in .mcool format and preparing results for visualization in IGV. The key steps involved include:

- Refer to the **Inputs & Outputs** section to verify required files and output structure.
- **Data Preparation**: Ensure .mcool files are formatted correctly and resolutions are verified.
- **Always prompt user** for resolution used to call loops.
- **Extract contact matrices** from .mcool files at the desired resolution.
- **Detect chromatin loops** using `hicDetectLoops` from HiCExplorer.

---

## When to Use This Skill

Use this skill when:

- You need to identify chromatin loops from Hi-C data in .mcool format.

--

## Inputs & Outputs

### Inputs

- **File format:** .mcool (Hi-C data file).
- **Resolution:** Choose the desired resolution for loop calling (e.g., 5 kb, 10 kb, etc.).
- **Target region:** Genome region for loop detection, if applicable.

### Outputs

```bash
loop_calling/
    ${sample}_loops_${res}.bedpe  # Detected chromatin loops in BEDPE format.
```
---

## Decision Tree

### Step 1: Extract Contact Matrix at a Given Resolution

Select the desired resolution (e.g., 5 kb) and extract it using `cooler`. First, check available resolutions in the .mcool file:

```bash
cooler ls <input.mcool>
```

This ensures the selected resolution is available for loop detection.

### Step 2: Detect Chromatin Loops with `hicDetectLoops`

Use HiCExplorer's `hicDetectLoops` to identify statistically significant loops.

```bash
hicDetectLoops -m input.mcool::/resolutions/<res> -o <sample>_loops_<res>.bedpe -p 0.05
```

**Notes:**
- The output file `<sample>_loops_<res>.bedpe` contains chromatin loop anchors and their scores.

---

## Best Practices

- Choose a resolution that is appropriate for your sequencing depth (e.g., 5–100 kb).
