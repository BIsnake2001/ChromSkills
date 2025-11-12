---
name: hic-loop-calling
description: This skill performs chromatin loop detection from Hi-C .mcool files using HiCExplorer's hicDetectLoops and generates IGV-ready interaction tracks for visualization.
---

# Hi-C Loop Calling and Visualization

## Overview

This skill provides a minimal and efficient workflow for detecting chromatin loops from Hi-C data stored in .mcool format and preparing results for visualization in IGV.  
It includes the following key steps:

1. Extract contact matrices from .mcool files at the desired resolution.  
2. Detect chromatin loops using `hicDetectLoops` from HiCExplorer.  
---

## When to Use This Skill

Use this skill when:
- You need to identify chromatin loops from Hi-C data in `.mcool` format.  
- You want to visualize loop structures directly in IGV.  
- You do not require compartment or expected contact calculations.

---

## Decision Tree

### Step 1: Extract Contact Matrix at a Given Resolution

Select the desired resolution (e.g., 5 kb) and extract it using `cooler`:

```bash
h5ls <input.mcool>resolutions
```

This ensures that the selected resolution is available for loop detection.

---

### Step 2: Detect Chromatin Loops with hicDetectLoops

Use HiCExplorer to identify statistically significant loops:

```bash
hicDetectLoops -m input.mcool::/resolutions/5000 -o loops_5kb.bedpe
```

**Notes:**
- The output `loops_5kb.bedpe` contains chromatin loop anchors and scores.  
- You can adjust the distance parameters based on sequencing depth and desired scale.

---

### Step 3: Visualize Loop Overlays on Contact Maps

Generate loop-annotated heatmaps to verify loop quality:

```bash
hicPlotMatrix --matrix input.mcool::/resolutions/5000 --loops loops_5kb.bedpe -o loops_heatmap.png
```

This produces a contact heatmap with loop arcs overlaid.

---

## Outputs

| File | Description |
|------|--------------|
| `loops_5kb.bedpe` | Detected chromatin loops in BEDPE format. |

---

## Best Practices

- Use ICE or KR normalized matrices when available in `.mcool` format.  
- Choose a resolution appropriate for your sequencing depth (e.g., 5–100 kb).  
- Verify that chromosome names in `.mcool` match your genome reference.  
- Adjust `--min-dist` and `--max-dist` for the genomic scale of interest.  
- For IGV visualization, ensure you open the `.interact` file with the correct genome assembly loaded.